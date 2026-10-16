"""Region segmentation by normal-coherent region growing + primitive fits."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import networkx as nx
import numpy as np
import trimesh

from .fit import (Cone, Cylinder, Plane, Sphere, fit_cylinder, fit_plane, fit_sphere,
                  refine_cylinder, refine_sphere, fit_cone, refine_cone, fit_cone_about_axis)


@dataclass
class Region:
    faces: np.ndarray
    surface: Plane | Cylinder | Sphere | Cone | None  # None = faceted fallback
    area: float


def _neighbors(m: trimesh.Trimesh):
    nb = [[] for _ in range(len(m.faces))]
    for a, b in m.face_adjacency:
        nb[a].append(b)
        nb[b].append(a)
    return nb


def segment(m: trimesh.Trimesh, tol: float, angle_deg: float = 20.0,
            min_area_frac: float = 0.002, cylinders: bool = True,
            spheres: bool = True, cones: bool = True) -> list[Region]:
    nb = _neighbors(m)
    # segmentation runs on a denoised copy; fits use original geometry
    sm = m.copy()
    trimesh.smoothing.filter_taubin(sm, iterations=10)
    fn, area, cen = sm.face_normals, m.area_faces, m.triangles_center
    tri = m.triangles
    label = -np.ones(len(m.faces), int)
    cos_t = np.cos(np.radians(angle_deg))
    dist_tol = max(3.0 * tol, 1e-9)
    regions: list[Region] = []

    # --- planes, seeded from largest faces
    for seed in np.argsort(-area):
        if label[seed] >= 0:
            continue
        rid = len(regions)
        members = [seed]
        label[seed] = rid
        plane = Plane(fn[seed], float(fn[seed] @ cen[seed]))
        q = deque([seed])
        grown = 0
        while q:
            f = q.popleft()
            for g in nb[f]:
                if label[g] >= 0:
                    continue
                if fn[g] @ plane.normal < cos_t:
                    continue
                if plane.distance(cen[g][None])[0] > dist_tol:
                    continue
                label[g] = rid
                members.append(g)
                q.append(g)
                grown += 1
                if grown % 50 == 0:
                    plane = _refit(tri, area, members, plane)
        members = np.array(members)
        plane = _refit(tri, area, members, plane)
        regions.append(Region(members, plane, float(area[members].sum())))

    total = area.sum()
    keep, leftover = [], []
    for r in regions:
        vp = m.vertices[np.unique(m.faces[r.faces])]
        ok = (r.area >= min_area_frac * total and _minor_extent(m, r) > 0.05 * m.scale
              and np.percentile(r.surface.distance(vp), 90) < 2.5 * tol)
        (keep if ok else leftover).append(r)

    keep, leftover, strip_cyls = _claim_strip_cylinders(m, sm, keep, leftover, tol)
    keep, more = _strip_chains(m, sm, keep, tol)
    leftover += more

    # --- leftovers: give them to planes, then fit cylinders or spheres
    out = _merge_planes(m, keep, tol) + strip_cyls
    if leftover:
        left = np.concatenate([r.faces for r in leftover])
        out = _absorb(m, nb, out + [Region(left, None, float(area[left].sum()))], tol)
        out = [r for r in out if r.surface is not None or len(r.faces)]
        rest = [r for r in out if r.surface is None]
        out = [r for r in out if r.surface is not None]
        if rest:
            left = np.concatenate([r.faces for r in rest])
            out += _cylinder_regions(m, sm, nb, left, tol, cylinders, spheres, cones)
    for _ in range(3):
        out = _merge_planes(m, out, tol)
        out = _merge_cylinders(m, out, tol)
        if spheres:
            out = _promote_spheres(m, sm, out, tol)
            out = _merge_spheres(m, out, tol)
        if cones:
            out = _promote_cones(m, sm, out, tol)
            out = _merge_cones(m, out, tol)
        out = _absorb(m, nb, out, tol)
        out = _dissolve_small(m, nb, out, tol)
        out = [r for r in out if len(r.faces)]
    out = _absorb(m, nb, out, tol, force_frac=0.005)
    if cones:
        out = _reclaim_coaxial_cones(m, sm, nb, out, tol)
    return [r for r in out if len(r.faces)]


def _strip_chains(m, sm, regs, tol, frac=0.05):
    """Chains of small adjacent planes that together form a cylinder
    (a faceted / scanned cylinder segmented as strips)."""
    lab = -np.ones(len(m.faces), int)
    for i, r in enumerate(regs):
        lab[r.faces] = i
    small = {i for i, r in enumerate(regs) if r.area < frac * m.area}
    pairs = set()
    for a, b in lab[m.face_adjacency]:
        if a != b and a in small and b in small:
            pairs.add((min(a, b), max(a, b)))
    g = nx.Graph()
    g.add_nodes_from(small)
    g.add_edges_from(pairs)
    moved = set()
    for comp in nx.connected_components(g):
        if len(comp) < 3:
            continue
        comp = set(comp)
        for _ in range(5):
            if len(comp) < 3:
                break
            faces = np.concatenate([regs[i].faces for i in comp])
            cyl = fit_cylinder(m.triangles_center[faces], sm.face_normals[faces])
            if cyl is None or cyl.radius > 0.5 * m.scale:
                break
            fit = {i for i in comp if np.percentile(
                cyl.distance(m.vertices[np.unique(m.faces[regs[i].faces])]), 90) < 3 * tol}
            if fit == comp:
                moved |= comp
                break
            comp = fit
    return ([r for i, r in enumerate(regs) if i not in moved],
            [r for i, r in enumerate(regs) if i in moved])


def _claim_strip_cylinders(m, sm, keep, leftover, tol, frac=0.05):
    """Turn chains of tangent planar strips into a cylinder region.

    A noisy cylinder is often region-grown as strips. Two adjacent strips
    already determine a radius; isolated strips on that same cylinder are
    taken with them so a later absorb cannot grow those planes across the
    wall. Cone-like strips (normals tilted toward the axis) are left alone.
    """
    if len(keep) < 2:
        return keep, leftover, []
    lab = -np.ones(len(m.faces), int)
    for i, r in enumerate(keep):
        if len(r.faces):
            lab[r.faces] = i
    small = {i for i, r in enumerate(keep)
             if 0 < r.area < frac * m.area and isinstance(r.surface, Plane)}
    if len(small) < 2:
        return keep, leftover, []
    g = nx.Graph()
    g.add_nodes_from(small)
    for a, b in lab[m.face_adjacency]:
        a, b = int(a), int(b)
        if a != b and a in small and b in small:
            g.add_edge(a, b)
    cands = []
    for comp in nx.connected_components(g):
        comp = set(comp)
        if len(comp) < 2:
            continue
        cyl = None
        for _ in range(5):
            if len(comp) < 2:
                break
            faces = np.concatenate([keep[i].faces for i in comp])
            cyl = fit_cylinder(m.triangles_center[faces], sm.face_normals[faces])
            if cyl is None or not 0 < cyl.radius <= 0.5 * m.scale:
                comp = set()
                break
            fit = {i for i in comp if np.percentile(
                cyl.distance(m.vertices[np.unique(m.faces[keep[i].faces])]), 90) < 3 * tol}
            if fit == comp:
                break
            comp = fit
            cyl = None
        if cyl is None or len(comp) < 2:
            continue
        dots = [abs(float(keep[i].surface.normal @ cyl.axis)) for i in comp]
        if float(np.median(dots)) > 0.35:
            continue
        cands.append((cyl, comp))
    if not cands:
        return keep, leftover, []
    cands.sort(key=lambda c: -sum(keep[i].area for i in c[1]))
    groups = []
    for cyl, comp in cands:
        for grp in groups:
            c0 = grp[0][0]
            coaxial = (abs(float(cyl.axis @ c0.axis)) > 0.98
                       and abs(cyl.radius - c0.radius) < 0.15 * max(c0.radius, 1e-9)
                       and np.linalg.norm(np.cross(cyl.point - c0.point, c0.axis)) < 0.15 * c0.radius)
            if coaxial:
                grp.append((cyl, comp))
                break
        else:
            groups.append([(cyl, comp)])
    drop = set()
    extra = []
    for grp in groups:
        comp = set()
        for _, c in grp:
            comp |= c
        if comp & drop:
            comp -= drop
        if len(comp) < 2:
            continue
        faces = np.concatenate([keep[i].faces for i in comp])
        cyl = fit_cylinder(m.triangles_center[faces], sm.face_normals[faces])
        if cyl is None or not 0 < cyl.radius <= 0.5 * m.scale:
            continue
        for i in list(small):
            if i in comp or i in drop:
                continue
            if abs(float(keep[i].surface.normal @ cyl.axis)) > 0.35:
                continue
            pts = m.vertices[np.unique(m.faces[keep[i].faces])]
            if np.percentile(cyl.distance(pts), 90) < 3 * tol:
                comp.add(i)
        faces = np.concatenate([keep[i].faces for i in comp])
        pts = m.vertices[np.unique(m.faces[faces])]
        cyl = refine_cylinder(cyl, pts)
        if not np.isfinite(cyl.radius) or not 0 < cyl.radius <= 0.5 * m.scale:
            continue
        if np.percentile(cyl.distance(pts), 90) > max(4 * tol, 0.02 * cyl.radius):
            continue
        # a spherical cap also breaks into strips; don't steal it when the
        # sphere explains those faces at least as well as the cylinder
        cens = m.triangles_center[faces]
        sph = fit_sphere(cens, sm.face_normals[faces])
        if sph is not None and np.percentile(sph.distance(cens), 90) <= np.percentile(cyl.distance(cens), 90):
            continue
        drop |= comp
        extra.append(Region(faces, cyl, float(m.area_faces[faces].sum())))
    if not drop:
        return keep, leftover, []
    drop_faces = set(np.concatenate([keep[i].faces for i in drop]).tolist())
    new_keep = [r for i, r in enumerate(keep) if i not in drop]
    new_left = []
    for r in leftover:
        faces = r.faces[~np.isin(r.faces, list(drop_faces))] if len(r.faces) else r.faces
        if len(faces):
            new_left.append(Region(faces, r.surface, float(m.area_faces[faces].sum())))
    return new_keep, new_left, extra


def _angular_gap(axis, origin, pts):
    """Largest gap in the angle of ``pts`` around ``axis`` (radians)."""
    axis = axis / max(np.linalg.norm(axis), 1e-12)
    u = np.cross(axis, [1.0, 0.0, 0.0])
    if np.linalg.norm(u) < 0.1:
        u = np.cross(axis, [0.0, 1.0, 0.0])
    u = u / np.linalg.norm(u)
    v = np.cross(axis, u)
    rel = pts - origin
    ang = np.arctan2(rel @ v, rel @ u)
    ang = np.sort(ang)
    if len(ang) < 3:
        return 2 * np.pi
    gaps = np.diff(ang)
    wrap = ang[0] + 2 * np.pi - ang[-1]
    return float(max(gaps.max(), wrap))


def _reclaim_coaxial_cones(m, sm, nb, regs, tol):
    """Reassemble a narrow cone that noise left as strips beside a cylinder.

    Countersinks and chamfers are coaxial with the neighbouring cylinder and
    perpendicular to the cap. Under noise the band splits into a few planar
    strips and leaks into the cylinder and the cap. Faces that fit a cone
    about that axis, and that the cylinder or cap do not explain better, are
    gathered back into one conical region.
    """
    if not any(isinstance(r.surface, Cylinder) for r in regs):
        return regs
    label = -np.ones(len(m.faces), int)
    for i, r in enumerate(regs):
        if len(r.faces):
            label[r.faces] = i
    raw_n = m.face_normals
    sm_n = sm.face_normals
    area = m.area_faces
    total = float(area.sum()) or 1.0
    claimed = np.zeros(len(m.faces), bool)
    found = []
    cyl_ids = [i for i, r in enumerate(regs) if isinstance(r.surface, Cylinder) and len(r.faces)]
    for ci in cyl_ids:
        cyl = regs[ci].surface
        axis0 = cyl.axis / max(np.linalg.norm(cyl.axis), 1e-12)
        neigh = set()
        for f in regs[ci].faces:
            for g in nb[f]:
                j = int(label[g])
                if j >= 0 and j != ci:
                    neigh.add(j)
        axis = axis0.copy()
        cap = None
        best_align = 0.0
        for j in neigh:
            if not isinstance(regs[j].surface, Plane) or regs[j].area < 0.01 * total:
                continue
            n = regs[j].surface.normal
            align = abs(float(n @ axis0))
            if align > best_align and align > 0.9:
                best_align = align
                cap = n
        # a short hole's axis wanders under noise; the cap plane normal does not
        if cap is not None and best_align < 0.999:
            axis = cap * (1.0 if float(cap @ axis0) >= 0 else -1.0)
            axis = axis / np.linalg.norm(axis)
        if _coaxial_cone_present(m, regs, neigh, axis, cyl, tol):
            continue
        # strips beside the cylinder are required; stray noisy wall faces are not a cone
        strips = []
        for j in neigh:
            r = regs[j]
            if not 0 < r.area <= 0.08 * total:
                continue
            if isinstance(r.surface, Cylinder):
                continue
            if r.surface is not None and not isinstance(r.surface, (Plane, Cone, Sphere)):
                continue
            dots = np.abs(raw_n[r.faces] @ axis)
            med = float(np.median(dots))
            # above ~0.85 the region is the cap (a leaked triangle or two pulls
            # the median down from 1). Real chamfer/countersink strips sit lower.
            if 0.30 <= med <= 0.82:
                strips.extend(int(f) for f in r.faces)
        if len(strips) < 8:
            continue
        seed = np.unique(np.asarray(strips, int))
        seed = seed[~claimed[seed]]
        # a strip often keeps going down the hole. Faces that still sit on the
        # cylinder AND look like a wall are not the cone; a chamfer's inner
        # edge sits on the cylinder too, but its normals stay cone-like.
        if len(seed):
            cd = cyl.distance(m.triangles_center[seed])
            wallish = (cd <= max(2.0 * tol, 0.30)) & (np.abs(raw_n[seed] @ axis) < 0.40)
            seed = seed[~wallish]
        if len(seed) < 8 or float(area[seed].sum()) < 0.0015 * total:
            continue
        dots = np.abs(raw_n[seed] @ axis)
        band = (dots > 0.25) & (dots < 0.95)
        if int(band.sum()) >= 8:
            med = float(np.median(dots[band]))
            inl = seed[np.abs(dots - med) < 0.2]
            if len(inl) >= 8:
                seed_fit = inl
            else:
                seed_fit = seed[band]
        else:
            seed_fit = seed
        if len(seed_fit) < 8:
            continue
        # a tilted cylinder fit puts cyl.point off the hole centre, and the
        # mean of the wall centroids shifts when a sector was labelled as strips.
        # A circle in the plane perpendicular to the axis does not.
        wall = regs[ci].faces
        if len(wall):
            near = cyl.distance(m.triangles_center[wall]) < max(3.0 * tol, 0.40)
            near &= np.abs(raw_n[wall] @ axis) < 0.45
            if int(near.sum()) >= 8:
                wall = wall[near]
        origin = _circle_center(m.triangles_center[wall], axis)
        if origin is None:
            origin = m.triangles_center[wall].mean(0)
        pts = m.vertices[np.unique(m.faces[seed_fit])]
        cens = m.triangles_center[seed_fit]
        cone = fit_cone_about_axis(origin, axis, cens, raw_n[seed_fit])
        if cone is None:
            cone = fit_cone_about_axis(origin, axis, pts, raw_n[seed_fit])
        if cone is None:
            continue
        cone = _lock_cone_axis(cone, axis, cens, raw_n[seed_fit])
        if cone is None:
            continue
        err = float(np.percentile(cone.distance(cens), 90))
        h = (cens - cone.apex) @ cone.axis
        rmax = float(cone.radius_at(np.max(h))) if np.max(h) > 0 else 0.0
        # strips are contaminated, so this gate only rejects a hopeless fit
        if err > max(14.0 * tol, 0.2 * max(rmax, 1.0)) or rmax > 0.5 * m.scale:
            continue
        chosen = _grow_cone_faces(m, nb, regs, label, ci, cyl, cone, seed_fit, claimed, tol, total, raw_n, sm_n)
        chosen = _absorb_cone_slivers(m, nb, regs, label, chosen, cone, tol, total)
        # growth can leak down the cylinder wall. Centroids stay near the seed
        # band; a vertex that reaches the far cap is a wall triangle.
        if len(chosen) and len(seed_fit):
            href = (m.triangles_center[seed_fit] - cone.apex) @ cone.axis
            lo, hi = np.percentile(href, [10, 90])
            pad = max(3.0 * tol, 0.35 * max(float(hi - lo), tol))
            hc = (m.triangles_center[chosen] - cone.apex) @ cone.axis
            chosen = chosen[(hc >= lo - pad) & (hc <= hi + pad)]
        if len(chosen) < 8 or float(area[chosen].sum()) < 0.002 * total:
            continue
        pts = m.vertices[np.unique(m.faces[chosen])]
        cone = fit_cone_about_axis(origin, axis, m.triangles_center[chosen], raw_n[chosen])
        if cone is None:
            cone = fit_cone_about_axis(origin, axis, pts, raw_n[chosen])
        if cone is None:
            continue
        cone = _lock_cone_axis(cone, axis, pts, raw_n[chosen])
        if cone is None:
            continue
        # boundary vertices are shared with the cap and the cylinder, so the
        # vertex residual sits above the noise; centroids track the band
        err_v = float(np.percentile(cone.distance(pts), 90))
        err_c = float(np.percentile(cone.distance(m.triangles_center[chosen]), 90))
        h = (pts - cone.apex) @ cone.axis
        if float(np.percentile(h, 5)) <= 0:
            continue
        rmax = float(cone.radius_at(float(np.max(h))))
        if (err_v > max(5.0 * tol, 0.04 * max(rmax, 1e-6))
                or err_c > max(3.5 * tol, 0.025 * max(rmax, 1e-6))
                or rmax > 0.5 * m.scale):
            continue
        if abs(float(cone.axis @ axis)) < 0.995:
            continue
        cens = m.triangles_center[chosen]
        if _angular_gap(cone.axis, cone.apex, cens) > np.radians(75):
            continue
        nagree = max(float(np.mean(np.abs(np.sum(n * cone.normal(cens), axis=1))))
                     for n in (raw_n[chosen], sm_n[chosen]))
        if nagree < 0.72:
            continue
        cyl_try = fit_cylinder(cens, sm_n[chosen])
        if cyl_try is not None and np.percentile(cyl_try.distance(pts), 90) < 0.75 * err_c:
            continue
        claimed[chosen] = True
        found.append((chosen, cone))
    if not found:
        return regs
    out = []
    for r in regs:
        if not len(r.faces):
            continue
        keep = r.faces[~claimed[r.faces]]
        if len(keep) == 0:
            continue
        if len(keep) == len(r.faces):
            out.append(r)
            continue
        surf = r.surface
        pts = m.vertices[np.unique(m.faces[keep])]
        if isinstance(surf, Cylinder) and len(keep) >= 8:
            surf = refine_cylinder(surf, pts)
        elif isinstance(surf, Plane) and len(keep) >= 3:
            surf = _refit(m.triangles, area, keep, surf)
        elif isinstance(surf, Cone) and len(keep) >= 8:
            surf = refine_cone(surf, pts)
        out.append(Region(keep, surf, float(area[keep].sum())))
    for faces, cone in found:
        out.append(Region(faces, cone, float(area[faces].sum())))
    out = _merge_cones(m, [r for r in out if len(r.faces)], tol)
    out = _merge_tiny_into_cones(m, nb, out, total, tol)
    out = _peel_cone_wall(m, nb, out, tol)
    out = _return_wall_to_cylinders(m, nb, out, total, tol)
    out = _absorb_coplanar_slivers(m, out, tol)
    return _realign_cones_to_cylinders(m, nb, out, tol)


def _merge_tiny_into_cones(m, nb, regs, total, tol, frac=0.02):
    """Pull leftover seam faces back onto the one cone they touch.

    Cutting the band out of a strip leaves corners between the cone and the
    cylinder. Those faces are taken one at a time, and only when they lie on
    the cone within an axial pad of the band. The rest of the strip is the
    wall and is left for the cylinder.
    """
    regs = [r for r in regs if len(r.faces)]
    area = m.area_faces
    for _ in range(6):
        label = -np.ones(len(m.faces), int)
        for i, r in enumerate(regs):
            if len(r.faces):
                label[r.faces] = i
        cones = [i for i, r in enumerate(regs) if isinstance(r.surface, Cone) and len(r.faces)]
        if not cones:
            break
        spans = {}
        for i in cones:
            href = (m.triangles_center[regs[i].faces] - regs[i].surface.apex) @ regs[i].surface.axis
            spans[i] = (float(np.percentile(href, 5)), float(np.percentile(href, 95)))
        small = set()
        for i, r in enumerate(regs):
            if i in spans or not len(r.faces) or r.area > frac * total:
                continue
            if r.surface is not None and not isinstance(r.surface, (Plane, Cone)):
                continue
            small.update(int(f) for f in r.faces)
        claim = {i: [] for i in cones}
        for f in small:
            touch = {int(label[g]) for g in nb[f] if int(label[g]) in spans}
            if len(touch) != 1:
                continue
            j = next(iter(touch))
            cone = regs[j].surface
            if abs(float(m.face_normals[f] @ cone.axis)) > 0.82:
                continue
            hc = float((m.triangles_center[f] - cone.apex) @ cone.axis)
            lo, hi = spans[j]
            pad = max(4.0 * tol, 1.20)
            if hc < lo - pad or hc > hi + pad:
                continue
            if float(cone.distance(m.triangles_center[f][None])[0]) > max(6.0 * tol, 0.75):
                continue
            claim[j].append(f)
        if not any(claim.values()):
            break
        taken = {f for fs in claim.values() for f in fs}
        new = []
        for i, r in enumerate(regs):
            if i in claim and claim[i]:
                faces = np.concatenate([r.faces, np.asarray(claim[i], int)])
                new.append(Region(faces, r.surface, float(area[faces].sum())))
                continue
            if not len(r.faces):
                continue
            keep = np.asarray([int(f) for f in r.faces if int(f) not in taken], int)
            if len(keep) == 0:
                continue
            if len(keep) == len(r.faces):
                new.append(r)
            else:
                new.append(Region(keep, r.surface, float(area[keep].sum())))
        regs = new
    return [r for r in regs if len(r.faces)]


def _return_wall_to_cylinders(m, nb, regs, total, tol, frac=0.02):
    """Give a strip's leftover wall faces back to the cylinder they touch.

    A narrow cone band and the hole below it arrive as one planar strip. After
    the cone is cut out, the tail still bridges the cylinder and the far cap,
    which leaves the cylinder boundary non-manifold.
    """
    regs = [r for r in regs if len(r.faces)]
    area = m.area_faces
    label = -np.ones(len(m.faces), int)
    for i, r in enumerate(regs):
        if len(r.faces):
            label[r.faces] = i
    claim = {}
    for i, r in enumerate(regs):
        if not isinstance(r.surface, Cylinder) or not len(r.faces):
            continue
        cyl = r.surface
        axis = cyl.axis / max(float(np.linalg.norm(cyl.axis)), 1e-12)
        got = []
        seen = set()
        for f in r.faces:
            for g in nb[f]:
                j = int(label[g])
                if j < 0 or j == i or j in seen:
                    continue
                o = regs[j]
                if not len(o.faces) or o.area > frac * total:
                    continue
                if o.surface is not None and not isinstance(o.surface, (Plane, Cone)):
                    continue
                seen.add(j)
                cd = cyl.distance(m.triangles_center[o.faces])
                dots = np.abs(m.face_normals[o.faces] @ axis)
                take = (cd <= max(3.0 * tol, 0.50)) & (dots < 0.45)
                got.extend(int(f) for f in o.faces[take])
        if got:
            claim[i] = got
    if not claim:
        return regs
    taken = {f for fs in claim.values() for f in fs}
    new = []
    for i, r in enumerate(regs):
        if i in claim:
            faces = np.unique(np.concatenate([r.faces, np.asarray(claim[i], int)]))
            pts = m.vertices[np.unique(m.faces[faces])]
            surf = refine_cylinder(r.surface, pts) if len(faces) >= 8 else r.surface
            new.append(Region(faces, surf, float(area[faces].sum())))
            continue
        keep = np.asarray([int(f) for f in r.faces if int(f) not in taken], int)
        if len(keep) == 0:
            continue
        if len(keep) == len(r.faces):
            new.append(r)
        else:
            new.append(Region(keep, r.surface, float(area[keep].sum())))
    return [r for r in new if len(r.faces)]


def _peel_cone_wall(m, nb, regs, tol):
    """Move cylinder-wall faces that growth pulled into a cone back out.

    The band is the cone faces whose normals match the half-angle. Faces
    outside that axial span that still sit on the neighbouring cylinder are
    the hole wall; leaving them in the cone pinches the cylinder boundary.
    """
    regs = [r for r in regs if len(r.faces)]
    label = -np.ones(len(m.faces), int)
    for i, r in enumerate(regs):
        if len(r.faces):
            label[r.faces] = i
    area = m.area_faces
    out = list(regs)
    for i, r in enumerate(regs):
        if not isinstance(r.surface, Cone) or len(r.faces) < 8:
            continue
        cyls = set()
        for f in r.faces:
            for g in nb[f]:
                j = int(label[g])
                if j >= 0 and isinstance(regs[j].surface, Cylinder):
                    cyls.add(j)
        if len(cyls) != 1:
            continue
        cj = next(iter(cyls))
        # a manifold cylinder already has its two loops; peeling faces out of
        # the cone would open a hole in a solid that already closes
        from .brep import _loops
        if _loops(m, label, cj) is not None:
            continue
        cyl = regs[cj].surface
        cone = r.surface
        dots = np.abs(m.face_normals[r.faces] @ cone.axis)
        band = r.faces[(dots > 0.40) & (dots < 0.88)]
        if len(band) < 8:
            continue
        href = (m.triangles_center[band] - cone.apex) @ cone.axis
        lo, hi = np.percentile(href, [10, 90])
        pad = max(2.5 * tol, 0.80)
        hc = (m.triangles_center[r.faces] - cone.apex) @ cone.axis
        cd = cyl.distance(m.triangles_center[r.faces])
        outside = (hc < lo - pad) | (hc > hi + pad)
        wall = outside & (cd <= max(3.0 * tol, 0.45)) & (dots < 0.50)
        if int(wall.sum()) < 1 or int((~wall).sum()) < 8:
            continue
        keep = r.faces[~wall]
        moved = r.faces[wall]
        out[i] = Region(keep, cone, float(area[keep].sum()))
        faces = np.concatenate([out[cj].faces, moved])
        out[cj] = Region(faces, cyl, float(area[faces].sum()))
    return [r for r in out if len(r.faces)]


def _absorb_coplanar_slivers(m, regs, tol, frac=0.008):
    """Fold a tiny plane into the larger plane its faces actually lie on.

    A handful of cap triangles left beside a cone get a nonsense fitted
    normal, so the ordinary plane merge misses them, and they split the
    cone's loop. Membership is decided from the face normals and the
    distance to the neighbour, not from that fit.
    """
    regs = [r for r in regs if len(r.faces)]
    area = m.area_faces
    total = float(area.sum()) or 1.0
    changed = True
    while changed:
        changed = False
        for i, r in enumerate(regs):
            if not isinstance(r.surface, Plane) or not len(r.faces) or r.area > frac * total:
                continue
            med = np.median(m.face_normals[r.faces], axis=0)
            nrm = float(np.linalg.norm(med))
            med = med / nrm if nrm > 1e-9 else med
            cens = m.triangles_center[r.faces]
            host = None
            host_dist = None
            for j, o in enumerate(regs):
                if j == i or not isinstance(o.surface, Plane) or o.area < max(r.area * 3, 1e-6):
                    continue
                dist = float(np.percentile(o.surface.distance(cens), 90))
                if dist > max(3.0 * tol, 0.15):
                    continue
                # noisy triangles on the cap have a scrambled fitted normal but
                # still lie on it; a real cone strip does not
                aligned = nrm > 1e-9 and float(med @ o.surface.normal) > 0.85
                if not aligned and dist > max(tol, 0.25):
                    continue
                if host is None or dist < host_dist:
                    host, host_dist = j, dist
            if host is None:
                continue
            faces = np.concatenate([regs[host].faces, r.faces])
            regs[host] = Region(faces, regs[host].surface, float(area[faces].sum()))
            regs[i] = Region(np.array([], int), r.surface, 0.0)
            changed = True
            break
        regs = [r for r in regs if len(r.faces)]
    return regs


def _circle_center(pts, axis):
    """Center of a circle in the plane perpendicular to ``axis``.

    The returned point lies on the axis line (not necessarily near ``pts``).
    """
    pts = np.asarray(pts, float)
    axis = np.asarray(axis, float)
    axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
    if len(pts) < 6:
        return None
    u = np.cross(axis, [1.0, 0.0, 0.0])
    if np.linalg.norm(u) < 0.1:
        u = np.cross(axis, [0.0, 1.0, 0.0])
    u = u / np.linalg.norm(u)
    v = np.cross(axis, u)
    x, y = pts @ u, pts @ v
    A = np.c_[2 * x, 2 * y, np.ones(len(pts))]
    sol, *_ = np.linalg.lstsq(A, x * x + y * y, rcond=None)
    cx, cy, k = sol
    if not np.isfinite(k) or k + cx * cx + cy * cy <= 0:
        return None
    return cx * u + cy * v


def _realign_cones_to_cylinders(m, nb, regs, tol):
    """Put a reclaimed cone and its cylinder on one axis.

    The exact seam is a circle only when the two axes coincide and the cone
    meets the cap plane squarely. A short hole's own axis wanders under noise;
    the cap normal does not, so both surfaces are locked to it when that still
    fits. Otherwise the cone follows the cylinder.
    """
    label = -np.ones(len(m.faces), int)
    for i, r in enumerate(regs):
        if len(r.faces):
            label[r.faces] = i
    total = float(m.area_faces.sum()) or 1.0
    out = list(regs)
    for i, r in enumerate(regs):
        if not isinstance(r.surface, Cone) or not len(r.faces):
            continue
        cyls = set()
        planes = set()
        for f in r.faces:
            for g in nb[f]:
                j = int(label[g])
                if j < 0:
                    continue
                if isinstance(regs[j].surface, Cylinder):
                    cyls.add(j)
                elif isinstance(regs[j].surface, Plane):
                    planes.add(j)
        if len(cyls) != 1:
            continue
        cj = next(iter(cyls))
        cyl = regs[cj].surface
        for f in regs[cj].faces:
            for g in nb[f]:
                j = int(label[g])
                if j >= 0 and isinstance(regs[j].surface, Plane):
                    planes.add(j)
        axis0 = cyl.axis / max(np.linalg.norm(cyl.axis), 1e-12)
        if float(r.surface.axis @ axis0) < 0:
            axis0 = -axis0
        cap = None
        best = 0.0
        for j in planes:
            if regs[j].area < 0.01 * total:
                continue
            n = regs[j].surface.normal
            align = abs(float(n @ axis0))
            if align > best and align > 0.95:
                best = align
                cap = n
        axis = axis0
        if cap is not None:
            axis = cap * (1.0 if float(cap @ axis0) >= 0 else -1.0)
            axis = axis / np.linalg.norm(axis)
        wall = regs[cj].faces
        origin = _circle_center(m.triangles_center[wall], axis)
        if origin is None:
            origin = cyl.point
        pts = m.vertices[np.unique(m.faces[r.faces])]
        cpts = m.vertices[np.unique(m.faces[wall])]
        locked = fit_cone_about_axis(origin, axis, pts, m.face_normals[r.faces])
        if locked is None or abs(float(locked.axis @ axis)) < 0.995:
            locked = None
        rad = None
        if len(cpts):
            radial = cpts - origin
            radial = radial - np.outer(radial @ axis, axis)
            rad = float(np.median(np.linalg.norm(radial, axis=1)))
        old_c = float(np.percentile(r.surface.distance(pts), 90))
        old_y = float(np.percentile(cyl.distance(cpts), 90)) if len(cpts) else 0.0
        use_cap = (locked is not None and cap is not None and rad is not None and rad > 0
                   and float(np.percentile(locked.distance(pts), 90)) <= max(old_c * 1.5, old_c + 0.05))
        if use_cap:
            trial = Cylinder(axis, origin, rad)
            new_y = float(np.percentile(trial.distance(cpts), 90))
            if new_y <= max(4.0 * tol, 0.12 * rad, old_y + 0.05):
                out[i] = Region(r.faces, locked, r.area)
                out[cj] = Region(regs[cj].faces, trial, regs[cj].area)
                continue
        # the cylinder axis was already good enough; only the cone has to follow it
        axis = axis0
        origin = cyl.point
        locked = fit_cone_about_axis(origin, axis, pts, m.face_normals[r.faces])
        if locked is None or abs(float(locked.axis @ axis)) < 0.995:
            continue
        new_c = float(np.percentile(locked.distance(pts), 90))
        if new_c <= max(old_c * 1.5, old_c + 0.05):
            out[i] = Region(r.faces, locked, r.area)
    return out


def _coaxial_cone_present(m, regs, neigh, axis, cyl, tol):
    """True when this cylinder already meets a cone good enough to keep."""
    for j in neigh:
        c = regs[j].surface
        if not isinstance(c, Cone) or abs(float(c.axis @ axis)) < 0.995:
            continue
        pts = m.vertices[np.unique(m.faces[regs[j].faces])]
        if len(pts) < 6:
            continue
        h = (pts - c.apex) @ c.axis
        if float(np.max(h)) <= 0:
            continue
        rmax = float(c.radius_at(float(np.max(h))))
        err = float(np.percentile(c.distance(pts), 90))
        if err <= max(4.0 * tol, 0.02 * max(rmax, cyl.radius)):
            return True
    return False


def _lock_cone_axis(cone, axis, pts, normals):
    """Keep the cone coaxial. A free refine may tilt a narrow band off-axis
    and the exact frustum builder then rejects the face."""
    axis = np.asarray(axis, float)
    axis = axis / max(np.linalg.norm(axis), 1e-12)
    if float(cone.axis @ axis) < 0:
        axis = -axis
    if abs(float(cone.axis @ axis)) >= 0.9995:
        return cone
    locked = fit_cone_about_axis(cone.apex, axis, pts, normals)
    if locked is None:
        return None
    e_free = float(np.percentile(cone.distance(pts), 90))
    e_lock = float(np.percentile(locked.distance(pts), 90))
    if e_lock > max(e_free * 1.35, e_free + 1e-6):
        return None
    return locked


def _absorb_cone_slivers(m, nb, regs, label, chosen, cone, tol, total, frac=0.015):
    """Close gaps in a narrow band without swallowing the cap.

    Growth leaves corners of a strip behind. A corner that still lies near the
    cone is part of the band; a face that lies on a large plane is the cap and
    must stay there. A face walled in on every side is a one-triangle hole.
    """
    if len(chosen) == 0:
        return chosen
    chosen_set = set(int(f) for f in chosen)
    limit = max(6.0 * tol, 0.75)

    def on_cap(f):
        tri = m.vertices[m.faces[f]]
        n = m.face_normals[f]
        for r in regs:
            if not isinstance(r.surface, Plane) or r.area <= 0.02 * total:
                continue
            if (float(r.surface.distance(tri).max()) < max(2.5 * tol, 0.25)
                    and abs(float(n @ r.surface.normal)) > 0.85):
                return True
        return False

    changed = True
    while changed:
        changed = False
        neigh = {}
        for f in chosen_set:
            for g in nb[f]:
                if g in chosen_set:
                    continue
                j = int(label[g])
                if j >= 0:
                    neigh.setdefault(j, set()).add(g)
        extra = []
        for j, faces in neigh.items():
            r = regs[j]
            small = 0 < r.area <= frac * total and (r.surface is None or isinstance(r.surface, (Plane, Cone)))
            for f in faces:
                walled = sum(1 for g in nb[f] if g in chosen_set) >= 3
                if walled:
                    # a one-triangle hole; leaving it splits the cap loop
                    extra.append(f)
                    continue
                if not small or on_cap(f):
                    continue
                tri = m.vertices[m.faces[f]]
                if float(cone.distance(tri).max()) <= limit:
                    extra.append(f)
        if extra:
            chosen_set.update(extra)
            changed = True
    return np.asarray(sorted(chosen_set), int)


def _grow_cone_faces(m, nb, regs, label, ci, cyl, cone, seed, claimed, tol, total, raw_n, sm_n):
    """Faces of the conical band: normal agrees, vertices lie on the cone,
    and neither the cylinder wall nor a large cap explains them better."""
    sin_a = abs(float(np.sin(cone.half_angle)))
    r_ref = max(float(cyl.radius), 1e-6)
    tol_d = max(4.0 * tol, 0.03 * r_ref)
    chosen = []
    chosen_set = set()
    q = deque(int(f) for f in seed)
    seen = set(q)

    def ok(f, dist_limit, norm_limit):
        if claimed[f]:
            return False
        dnorm = min(abs(abs(float(raw_n[f] @ cone.axis)) - sin_a),
                    abs(abs(float(sm_n[f] @ cone.axis)) - sin_a))
        if dnorm > norm_limit:
            return False
        tri = m.vertices[m.faces[f]]
        if float(cone.distance(tri).max()) > dist_limit:
            return False
        j = int(label[f])
        if j >= 0 and isinstance(regs[j].surface, Plane):
            on_plane = float(regs[j].surface.distance(tri).max()) < max(3.0 * tol, 0.40)
            # a cap is perpendicular to the axis whatever its area; the side
            # faces of the solid are the large planes
            cap_like = abs(float(regs[j].surface.normal @ cone.axis)) > 0.90
            large = regs[j].area > 0.05 * total and abs(float(raw_n[f] @ regs[j].surface.normal)) > 0.9
            if on_plane and (cap_like or large):
                return False
        if j == ci or (j >= 0 and isinstance(regs[j].surface, Cylinder)):
            cd = float(cyl.distance(tri).max())
            if cd + max(tol, 0.05) < float(cone.distance(tri).max()):
                return False
            if cd < max(2.5 * tol, 0.3) and abs(float(raw_n[f] @ cone.axis)) < sin_a * 0.55:
                return False
        return True

    while q:
        f = q.popleft()
        if not ok(f, tol_d, 0.30):
            continue
        chosen.append(f)
        chosen_set.add(f)
        for g in nb[f]:
            if g not in seen:
                seen.add(g)
                q.append(g)
    # pinholes: a face pinched between cone faces, slightly noisier than the gate
    for _ in range(3):
        extra = []
        fringe = set()
        for f in chosen:
            for g in nb[f]:
                if g not in chosen_set and not claimed[g]:
                    fringe.add(g)
        for g in fringe:
            n_touch = sum(1 for h in nb[g] if h in chosen_set)
            if n_touch >= 2 and ok(g, tol_d * 1.35, 0.42):
                extra.append(g)
        if not extra:
            break
        for g in extra:
            chosen.append(g)
            chosen_set.add(g)
    return np.asarray(chosen, int)


def _minor_extent(m, r):
    pts = m.triangles_center[r.faces]
    if len(pts) < 3:
        return 0.0
    n = r.surface.normal
    u = np.cross(n, [1.0, 0, 0])
    if np.linalg.norm(u) < 0.1:
        u = np.cross(n, [0, 1.0, 0])
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    q = np.c_[pts @ u, pts @ v]
    q -= q.mean(0)
    s = np.linalg.svd(q, compute_uv=False)
    return float(4 * s[-1] / np.sqrt(len(q)))  # ~ full width for uniform strip


def _merge_planes(m, regs, tol, ang_deg=3.0):
    """Merge plane regions lying on the same plane (split by noise/edges)."""
    planes = [r for r in regs if isinstance(r.surface, Plane)]
    rest = [r for r in regs if not isinstance(r.surface, Plane)]
    planes.sort(key=lambda r: -r.area)
    cos_t = np.cos(np.radians(ang_deg))
    groups = []
    for r in planes:
        for g in groups:
            p = g[0].surface
            pts = m.triangles_center[r.faces]
            if r.surface.normal @ p.normal > cos_t and np.percentile(p.distance(pts), 90) < 3 * tol:
                g.append(r)
                break
        else:
            groups.append([r])
    out = []
    for g in groups:
        faces = np.concatenate([r.faces for r in g])
        surf = g[0].surface if len(g) == 1 else _refit(m.triangles, m.area_faces, faces, g[0].surface)
        out.append(Region(faces, surf, float(m.area_faces[faces].sum())))
    return out + rest


def _merge_spheres(m, regs, tol):
    sph = [r for r in regs if isinstance(r.surface, Sphere)]
    rest = [r for r in regs if not isinstance(r.surface, Sphere)]
    sph.sort(key=lambda r: -r.area)
    groups = []
    for r in sph:
        for g in groups:
            s = g[0].surface
            same = (abs(r.surface.radius - s.radius) < 0.15 * s.radius
                    and np.linalg.norm(r.surface.center - s.center) < 0.15 * s.radius)
            if same or np.percentile(s.distance(m.triangles_center[r.faces]), 90) < 3 * tol:
                g.append(r)
                break
        else:
            groups.append([r])
    out = []
    for g in groups:
        faces = np.concatenate([r.faces for r in g])
        surf = refine_sphere(g[0].surface, m.vertices[np.unique(m.faces[faces])])
        out.append(Region(faces, surf, float(m.area_faces[faces].sum())))
    return rest + out


def _accept_sphere(m, sm, faces, tol) -> Sphere | None:
    """Fit a sphere to faces and accept it only when it beats a plane or cylinder."""
    if len(faces) < 8 or m.area_faces[faces].sum() < 0.002 * m.area:
        return None
    cen = m.triangles_center[faces]
    nrm = sm.face_normals[faces]
    sph = fit_sphere(cen, nrm)
    if sph is None:
        return None
    pts = m.vertices[np.unique(m.faces[faces])]
    sph = refine_sphere(sph, pts)
    if not np.isfinite(sph.radius) or sph.radius <= 0:
        return None
    ext = float(np.linalg.norm(pts - pts.mean(0), axis=1).max())
    if ext < 1e-9 or sph.radius > max(8.0 * ext, 2.0 * m.scale):
        return None
    rho = min(ext, sph.radius * 0.999)
    sag = sph.radius - float(np.sqrt(max(sph.radius ** 2 - rho ** 2, 0.0)))
    if sag < max(3.0 * tol, 0.015 * sph.radius):
        return None
    err = float(np.percentile(sph.distance(pts), 95))
    if err > max(5.0 * tol, 0.04 * sph.radius):
        return None
    radial = cen - sph.center
    radial /= np.maximum(np.linalg.norm(radial, axis=1, keepdims=True), 1e-12)
    if float(np.mean(np.abs(np.sum(nrm * radial, axis=1)))) < 0.8:
        return None
    cyl = fit_cylinder(cen, nrm)
    if cyl is not None and 0 < cyl.radius < 0.5 * m.scale:
        cerr = float(np.percentile(cyl.distance(pts), 95))
        if cerr < max(4.0 * tol, 0.02 * cyl.radius) and cerr < err * 0.75:
            return None
    return sph


def _promote_spheres(m, sm, regs, tol, frac=0.05, accept=None):
    """Merge a connected cluster of small regions when they lie on one sphere.

    Clean caps arrive as a single unfitted component (handled in
    _cylinder_regions). Noisy / faceted caps are split into small planes
    (and sometimes a spurious cylinder); those clusters are promoted here.
    """
    lab = -np.ones(len(m.faces), int)
    for i, r in enumerate(regs):
        if len(r.faces):
            lab[r.faces] = i
    small = {i for i, r in enumerate(regs) if 0 < r.area < frac * m.area}
    if len(small) < 1:
        return regs
    pairs = set()
    for a, b in lab[m.face_adjacency]:
        a, b = int(a), int(b)
        if a != b and a in small and b in small:
            pairs.add((min(a, b), max(a, b)))
    g = nx.Graph()
    g.add_nodes_from(small)
    g.add_edges_from(pairs)
    drop = set()
    extra = []
    for comp in nx.connected_components(g):
        comp = set(comp)
        kinds = [regs[i].surface for i in comp]
        if len(comp) < 2 and not (len(comp) == 1 and kinds[0] is None):
            continue
        if len(comp) == 1 and isinstance(kinds[0], (Sphere, Cone)):
            continue
        faces = np.concatenate([regs[i].faces for i in comp])
        sph = (accept or _accept_sphere)(m, sm, faces, tol)
        if sph is None:
            continue
        drop |= comp
        extra.append(Region(faces, sph, float(m.area_faces[faces].sum())))
    if not drop:
        return regs
    return [r for i, r in enumerate(regs) if i not in drop] + extra


def _merge_cylinders(m, regs, tol):
    cyl = [r for r in regs if isinstance(r.surface, Cylinder)]
    rest = [r for r in regs if not isinstance(r.surface, Cylinder)]
    cyl.sort(key=lambda r: -r.area)
    groups = []
    for r in cyl:
        for g in groups:
            c = g[0].surface
            same_axis = (abs(r.surface.axis @ c.axis) > 0.98
                         and abs(r.surface.radius - c.radius) < 0.15 * c.radius
                         and np.linalg.norm(np.cross(r.surface.point - c.point, c.axis)) < 0.15 * c.radius)
            if same_axis or np.percentile(c.distance(m.triangles_center[r.faces]), 90) < 3 * tol:
                g.append(r)
                break
        else:
            groups.append([r])
    out = []
    for g in groups:
        faces = np.concatenate([r.faces for r in g])
        surf = refine_cylinder(g[0].surface, m.vertices[np.unique(m.faces[faces])])
        out.append(Region(faces, surf, float(m.area_faces[faces].sum())))
    return rest + out


def _dissolve_small(m, nb, regs, tol, frac=0.05):
    """Small fitted regions that neighbours explain anyway get absorbed."""
    for i in np.argsort([r.area for r in regs]):
        r = regs[i]
        if r.surface is None or r.area >= frac * m.area or not len(r.faces):
            continue
        trial = list(regs)
        trial[i] = Region(r.faces, None, r.area)
        trial = _absorb(m, nb, trial, tol)
        if 0 < len(trial[i].faces) < 0.15 * len(r.faces) or r.area < 0.002 * m.area:
            trial = _absorb(m, nb, trial, tol, max_d=np.inf)
        if len(trial[i].faces) == 0:
            regs = trial
    return regs


def _absorb(m, nb, regs, tol, rounds=50, force_frac=0.0, max_d=None):
    """Give faceted faces to an adjacent fitted region when they fit it."""
    label = -np.ones(len(m.faces), int)
    for i, r in enumerate(regs):
        label[r.faces] = i
    fitted = np.array([r.surface is not None for r in regs])
    tri = m.triangles
    if force_frac > 0:
        # only force-absorb faceted connected bits that are tiny
        small = np.zeros(len(m.faces), bool)
        for comp in trimesh.graph.connected_components(
                m.face_adjacency[(~fitted[label[m.face_adjacency]]).all(1)], nodes=np.where(~fitted[label])[0]):
            if m.area_faces[comp].sum() < force_frac * m.area:
                small[comp] = True
    else:
        small = None
    for _ in range(rounds):
        changed = False
        cand = np.where(~fitted[label])[0]
        if small is not None:
            cand = cand[small[cand]]
        for f in cand:
            best, bd = None, np.inf
            for g in nb[f]:
                j = label[g]
                if fitted[j]:
                    d = regs[j].surface.distance(tri[f]).max()
                    if d < bd:
                        best, bd = j, d
            if best is not None and (bd < (5 * tol if max_d is None else max_d) or force_frac > 0):
                label[f] = best
                changed = True
        if not changed:
            break
    return [Region(np.where(label == i)[0], r.surface, float(m.area_faces[label == i].sum()))
            for i, r in enumerate(regs)]


def _refit(tri, area, members, plane):
    pts = tri[members].reshape(-1, 3)
    p = fit_plane(pts, np.repeat(area[members], 3))
    if p.normal @ plane.normal < 0:
        p = Plane(-p.normal, -p.d)
    return p


def _cylinder_regions(m, sm, nb, faces, tol, cylinders, spheres=True, cones=True, crease=None):
    fset = set(faces.tolist())
    seen = set()
    regs = []
    for s in faces:
        if s in seen:
            continue
        comp, q = [], deque([s])
        seen.add(s)
        while q:
            f = q.popleft()
            comp.append(f)
            for g in nb[f]:
                if g in fset and g not in seen and (crease is None or (min(f, g), max(f, g)) not in crease):
                    seen.add(g)
                    q.append(g)
        comp = np.array(comp)
        a = float(m.area_faces[comp].sum())
        surf = None
        if cylinders and len(comp) >= 8:
            cyl = fit_cylinder(m.triangles_center[comp], sm.face_normals[comp])
            if cyl is not None:
                pts = m.vertices[np.unique(m.faces[comp])]
                if (cyl.radius < 0.5 * m.scale and a >= 0.002 * m.area
                        and np.percentile(cyl.distance(pts), 95) < max(4 * tol, 0.02 * cyl.radius)):
                    surf = cyl
        if spheres:
            sph = _accept_sphere(m, sm, comp, tol)
            if sph is not None and (surf is None or np.percentile(
                    sph.distance(m.vertices[np.unique(m.faces[comp])]), 95)
                    < 0.75 * np.percentile(surf.distance(m.vertices[np.unique(m.faces[comp])]), 95)):
                surf = sph
        if cones:
            cone = _accept_cone(m, sm, comp, tol)
            if cone is not None:
                if surf is None:
                    surf = cone
                else:
                    # a narrow chamfer sits about as close to a large sphere as to
                    # its cone; only cylinders must win by a clear margin
                    pts = m.vertices[np.unique(m.faces[comp])]
                    factor = 0.95 if isinstance(surf, Sphere) else 0.75
                    if np.percentile(cone.distance(pts), 95) < factor * np.percentile(surf.distance(pts), 95):
                        surf = cone
        if surf is None and len(comp) >= 3:
            pl = fit_plane(m.triangles[comp].reshape(-1, 3))
            if np.percentile(pl.distance(m.vertices[np.unique(m.faces[comp])]), 90) < 2.5 * tol:
                if pl.normal @ m.face_normals[comp].mean(0) < 0:
                    pl = Plane(-pl.normal, -pl.d)
                surf = pl
        if surf is None and cones and crease is None and len(comp) >= 16:
            # coaxial cylinder + cone (countersink, chamfered boss) arrive as one
            # unfitted component: split it along sharp creases and fit the pieces
            cr = _creases(m, sm, comp)
            sub = [] if cr is None else _cylinder_regions(m, sm, nb, comp, tol, cylinders,
                                                         spheres, cones, crease=cr)
            if len(sub) > 1 and any(r.surface is not None for r in sub):
                regs += sub
                continue
        regs.append(Region(comp, surf, a))
    return regs


def _accept_cone(m, sm, faces, tol) -> Cone | None:
    """Fit a cone (chamfer, countersink, tapered boss); accept it only when it
    fits well and a cylinder would not."""
    if len(faces) < 8 or m.area_faces[faces].sum() < 0.002 * m.area:
        return None
    cen = m.triangles_center[faces]
    pts = m.vertices[np.unique(m.faces[faces])]
    # smoothed normals are robust to noise but blur narrow bands (chamfers);
    # try both and keep the better refined fit
    best = None
    for nrm in (sm.face_normals[faces], m.face_normals[faces]):
        c = fit_cone(cen, nrm)
        if c is None:
            continue
        c = refine_cone(c, pts)
        e = float(np.percentile(c.distance(pts), 90))
        if best is None or e < best[0]:
            best = (e, c, nrm)
    if best is None:
        return None
    _, cone, nrm = best
    h = (pts - cone.apex) @ cone.axis
    if h.min() <= 0:
        return None
    rmax = float(cone.radius_at(h.max()))
    if rmax > 0.5 * m.scale or np.linalg.norm(cone.apex - pts.mean(0)) > 2 * m.scale:
        return None
    # p90: boundary vertices of noisy strips belong half to the neighbours
    err = float(np.percentile(cone.distance(pts), 90))
    if err > max(4.0 * tol, 0.02 * rmax):
        return None
    # normals must agree with the cone's (either orientation: holes flip them)
    if max(float(np.mean(np.abs(np.sum(n * cone.normal(cen), axis=1))))
           for n in (sm.face_normals[faces], m.face_normals[faces])) < 0.8:
        return None
    cyl = fit_cylinder(cen, sm.face_normals[faces])
    if cyl is not None and np.percentile(cyl.distance(pts), 90) < 0.75 * err:
        return None
    return cone


def _merge_cones(m, regs, tol):
    cs = [r for r in regs if isinstance(r.surface, Cone)]
    rest = [r for r in regs if not isinstance(r.surface, Cone)]
    cs.sort(key=lambda r: -r.area)
    groups = []
    for r in cs:
        for g in groups:
            if np.percentile(g[0].surface.distance(m.triangles_center[r.faces]), 90) < 3 * tol:
                g.append(r)
                break
        else:
            groups.append([r])
    out = []
    for g in groups:
        faces = np.concatenate([r.faces for r in g])
        surf = g[0].surface if len(g) == 1 else refine_cone(g[0].surface, m.vertices[np.unique(m.faces[faces])])
        out.append(Region(faces, surf, float(m.area_faces[faces].sum())))
    return rest + out



def _creases(m, sm, comp):
    """Split a coaxial revolved component by normal-to-axis angle.

    For a surface of revolution the axis is the eigenvector of the normal
    covariance whose eigenvalue stands apart; a cylinder has n.a ~ 0, a cone
    n.a ~ const, so a 1-D 2-means on n.a separates them. Returns the set of
    adjacency pairs crossing the split (or None)."""
    n = sm.face_normals[comp]
    w, V = np.linalg.eigh(np.cov(n.T))
    k = int(np.argmax([abs(w[0] - w[1]) + abs(w[0] - w[2]), 0.0, abs(w[2] - w[1]) + abs(w[2] - w[0])]))
    t = n @ V[:, k]
    c = np.percentile(t, [10, 90])
    for _ in range(30):
        lab = np.abs(t - c[0]) > np.abs(t - c[1])
        if lab.all() or not lab.any():
            return None
        c = np.array([t[~lab].mean(), t[lab].mean()])
    if abs(c[1] - c[0]) < 0.2:
        return None
    side = np.full(len(m.faces), -1)
    side[comp] = lab
    fa = m.face_adjacency
    cross = fa[(side[fa[:, 0]] >= 0) & (side[fa[:, 1]] >= 0) & (side[fa[:, 0]] != side[fa[:, 1]])]
    return {(int(min(a, b)), int(max(a, b))) for a, b in cross}


def _promote_cones(m, sm, regs, tol, frac=0.05):
    """Merge clusters of small planar/unfitted strips that tile a cone.

    Noisy chamfers and drafted walls segment into thin plane strips. A cluster
    often also touches the cap plane next to it, so besides the whole cluster
    each leave-one-out subset is tried and the largest accepted cone wins."""
    lab = -np.ones(len(m.faces), int)
    for i, r in enumerate(regs):
        if len(r.faces):
            lab[r.faces] = i
    small = {i for i, r in enumerate(regs)
             if 0 < r.area < frac * m.area and (r.surface is None or isinstance(r.surface, Plane))}
    if len(small) < 2:
        return regs
    g = nx.Graph()
    g.add_nodes_from(small)
    for a, b in lab[m.face_adjacency]:
        a, b = int(a), int(b)
        if a != b and a in small and b in small:
            g.add_edge(a, b)
    drop, extra = set(), []
    for comp in nx.connected_components(g):
        comp = sorted(comp)
        if len(comp) < 2:
            continue
        subsets = [comp]
        if 3 <= len(comp) <= 16:
            subsets += [[j for j in comp if j != i] for i in comp]
        best = None
        for sub in subsets:
            if best is not None and len(sub) < len(best[0]):
                break
            faces = np.concatenate([regs[i].faces for i in sub])
            if not nx.is_connected(g.subgraph(sub)):
                continue
            cone = _accept_cone(m, sm, faces, tol)
            if cone is None:
                continue
            area = float(m.area_faces[faces].sum())
            if best is None or area > best[2]:
                best = (sub, cone, area, faces)
        if best is None:
            continue
        drop |= set(best[0])
        extra.append(Region(best[3], best[1], best[2]))
    if not drop:
        return regs
    return [r for i, r in enumerate(regs) if i not in drop] + extra
