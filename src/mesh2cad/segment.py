"""Region segmentation by normal-coherent region growing + primitive fits."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import networkx as nx
import numpy as np
import trimesh

from .fit import (Cylinder, Plane, Sphere, fit_cylinder, fit_plane, fit_sphere,
                  refine_cylinder, refine_sphere)


@dataclass
class Region:
    faces: np.ndarray
    surface: Plane | Cylinder | Sphere | None  # None = faceted fallback
    area: float


def _neighbors(m: trimesh.Trimesh):
    nb = [[] for _ in range(len(m.faces))]
    for a, b in m.face_adjacency:
        nb[a].append(b)
        nb[b].append(a)
    return nb


def segment(m: trimesh.Trimesh, tol: float, angle_deg: float = 20.0,
            min_area_frac: float = 0.002, cylinders: bool = True,
            spheres: bool = True) -> list[Region]:
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

    keep, more = _strip_chains(m, sm, keep, tol)
    leftover += more

    # --- leftovers: give them to planes, then fit cylinders or spheres
    out = _merge_planes(m, keep, tol)
    if leftover:
        left = np.concatenate([r.faces for r in leftover])
        out = _absorb(m, nb, out + [Region(left, None, float(area[left].sum()))], tol)
        out = [r for r in out if r.surface is not None or len(r.faces)]
        rest = [r for r in out if r.surface is None]
        out = [r for r in out if r.surface is not None]
        if rest:
            left = np.concatenate([r.faces for r in rest])
            out += _cylinder_regions(m, sm, nb, left, tol, cylinders, spheres)
    for _ in range(3):
        out = _merge_planes(m, out, tol)
        out = _merge_cylinders(m, out, tol)
        if spheres:
            out = _promote_spheres(m, sm, out, tol)
            out = _merge_spheres(m, out, tol)
        out = _absorb(m, nb, out, tol)
        out = _dissolve_small(m, nb, out, tol)
        out = [r for r in out if len(r.faces)]
    out = _absorb(m, nb, out, tol, force_frac=0.005)
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


def _promote_spheres(m, sm, regs, tol, frac=0.05):
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
        if len(comp) == 1 and isinstance(kinds[0], Sphere):
            continue
        faces = np.concatenate([regs[i].faces for i in comp])
        sph = _accept_sphere(m, sm, faces, tol)
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


def _cylinder_regions(m, sm, nb, faces, tol, cylinders, spheres=True):
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
                if g in fset and g not in seen:
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
        if surf is None and len(comp) >= 3:
            pl = fit_plane(m.triangles[comp].reshape(-1, 3))
            if np.percentile(pl.distance(m.vertices[np.unique(m.faces[comp])]), 90) < 2.5 * tol:
                if pl.normal @ m.face_normals[comp].mean(0) < 0:
                    pl = Plane(-pl.normal, -pl.d)
                surf = pl
        regs.append(Region(comp, surf, a))
    return regs
