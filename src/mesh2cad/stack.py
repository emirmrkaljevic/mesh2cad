"""Union of partial-depth extrusions.

The intersection of full-length silhouettes is the wrong solid when a profile
keeps its shape for only part of the length: a base, then a taller block on
part of it, a roof that is a tilted plane, a blind pocket, a slot, or a hole
that opens on an intermediate face. This module builds that union. It is only
tried after the intersection model has already missed the mesh.
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
import trimesh

from .extrude import _plane_coord, _shoelace, _unit
from .fit import Cone, Cylinder, Plane
from .parametric import _nice, _principal_axis

# A joint between two extrusions is overlapped by this much so the union is one
# solid. Outer ends are not grown, so the part does not get longer.
_GROW = 0.02
# How far a 45 degree edge bevel is pushed back into the body before its plane
# clips it. Face contact alone does not fuse.
_BEVEL_IN = 0.15


def _plane(region):
    n = np.asarray(region.surface.normal, float)
    length = float(np.linalg.norm(n))
    if length < 1e-12:
        return None
    return n / length, float(region.surface.d) / length


def _rings(mesh, region):
    """Boundary rings of one fitted face, as coordinate arrays. Largest first."""
    counts = defaultdict(int)
    for tri in mesh.faces[region.faces]:
        a, b, c = (int(tri[0]), int(tri[1]), int(tri[2]))
        for u, v in ((a, b), (b, c), (c, a)):
            counts[(u, v) if u < v else (v, u)] += 1
    nxt = defaultdict(list)
    for tri in mesh.faces[region.faces]:
        a, b, c = (int(tri[0]), int(tri[1]), int(tri[2]))
        for u, v in ((a, b), (b, c), (c, a)):
            if counts[(u, v) if u < v else (v, u)] == 1:
                nxt[u].append(v)
    used = set()
    rings = []
    for start in list(nxt):
        if start in used:
            continue
        ring = [start]
        used.add(start)
        prev, cur = None, start
        guard = 0
        while guard < 100000:
            guard += 1
            nxts = [p for p in nxt[cur] if p != prev and p not in used]
            if not nxts:
                if start in nxt[cur]:
                    ring.append(start)
                break
            prev, cur = cur, nxts[0]
            ring.append(cur)
            used.add(cur)
        if len(ring) >= 4:
            rings.append(mesh.vertices[ring])
    rings.sort(key=lambda p: -_ring_area(p))
    return rings


def _ring_area(pts):
    pts = np.asarray(pts, float)
    if len(pts) >= 2 and float(np.linalg.norm(pts[0] - pts[-1])) < 1e-6:
        pts = pts[:-1]
    if len(pts) < 3:
        return 0.0
    centered = pts - pts.mean(0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    uv = centered @ vt[:2].T
    return _loop_area(uv)


def _loop_area(xy):
    if len(xy) < 3:
        return 0.0
    pts = [(float(p[0]), float(p[1])) for p in xy]
    return abs(_shoelace(pts))


def _section_loops(mesh, z):
    sec = mesh.section(plane_origin=[0, 0, float(z)], plane_normal=[0, 0, 1])
    if sec is None:
        return []
    loops = []
    for path in sec.discrete:
        xy = np.asarray(path[:, :2], float)
        if len(xy) < 4 or float(np.linalg.norm(xy[0] - xy[-1])) > 1e-3:
            continue
        loops.append((abs(_shoelace([(float(p[0]), float(p[1])) for p in xy[:-1]])), xy[:-1]))
    loops.sort(key=lambda t: -t[0])
    return loops


def _polyline_distance(a, b):
    """Max distance from points of each polyline to the other."""
    def one(src, dst):
        if len(dst) < 2 or len(src) == 0:
            return 1e9
        worst = 0.0
        for p in src:
            d = np.linalg.norm(dst - p, axis=1)
            # also the segments, not just the vertices
            best = float(d.min())
            for i in range(len(dst)):
                q = dst[(i + 1) % len(dst)] - dst[i]
                L = float(np.dot(q, q))
                if L < 1e-12:
                    continue
                t = float(np.dot(p - dst[i], q) / L)
                t = min(1.0, max(0.0, t))
                best = min(best, float(np.linalg.norm(p - (dst[i] + t * q))))
            worst = max(worst, best)
        return worst
    return max(one(a, b), one(b, a))


def _simplify(pts, lim):
    seq = [(float(p[0]), float(p[1])) for p in pts]
    if seq and abs(seq[0][0] - seq[-1][0]) < 1e-9 and abs(seq[0][1] - seq[-1][1]) < 1e-9:
        seq = seq[:-1]
    changed = True
    while changed and len(seq) > 3:
        changed = False
        n = len(seq)
        for i in range(n):
            a = np.asarray(seq[(i - 1) % n], float)
            b = np.asarray(seq[i], float)
            c = np.asarray(seq[(i + 1) % n], float)
            ac = c - a
            norm = float(np.linalg.norm(ac))
            if norm < 1e-9:
                dist = float(np.linalg.norm(b - a))
            else:
                dist = abs(float(ac[0] * (a[1] - b[1]) - ac[1] * (a[0] - b[0]))) / norm
            if dist <= lim or float(np.linalg.norm(b - a)) < lim:
                del seq[i]
                changed = True
                break
    if _shoelace(seq) < 0:
        seq = list(reversed(seq))
    return seq


def _shelves(mesh, regions):
    """Horizontal caps large enough to be a real change of cross-section."""
    found = []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 40.0:
            continue
        got = _plane(region)
        if got is None or abs(got[0][2]) < 0.97:
            continue
        z = _plane_coord(region.surface, 2)
        found.append((float(z), float(region.area)))
    if not found:
        return []
    found.sort()
    bands = [[found[0]]]
    for z, area in found[1:]:
        if abs(z - bands[-1][-1][0]) <= 0.4:
            bands[-1].append((z, area))
        else:
            bands.append([(z, area)])
    out = []
    for band in bands:
        area = sum(a for _, a in band)
        z = sum(z * a for z, a in band) / area
        out.append((float(z), float(area)))
    return out


def _stable_outline(mesh, z0, z1, tol):
    """Outer loop when the cross-section does not move between two shelves."""
    if z1 - z0 < max(0.4, 4.0 * tol):
        return None
    heights = [z0 + 0.25 * (z1 - z0), 0.5 * (z0 + z1), z1 - 0.25 * (z1 - z0)]
    # stay off a chamfer or a fillet right on the shelf
    heights = [min(z1 - 0.3, max(z0 + 0.3, h)) for h in heights]
    loops = []
    for z in heights:
        got = _section_loops(mesh, z)
        if not got:
            return None
        loops.append(got[0][1])
    lim = max(4.0 * tol, 0.2)
    for other in loops[1:]:
        if _polyline_distance(loops[0], other) > lim:
            return None
        if abs(_loop_area(loops[0]) - _loop_area(other)) > 0.04 * max(_loop_area(loops[0]), 1.0):
            return None
    # A clean tessellation keeps a round in the outline (a 0.02 mm chord).
    # Noise is allowed to move a vertex by about the tolerance.
    from .segment import _analytic_mesh
    lim = 0.02 if _analytic_mesh(mesh) else max(tol, 0.05)
    return _simplify(loops[1], lim)


def _bbox(pts):
    p = np.asarray(pts, float)
    return p.min(0), p.max(0)


def _face_points(mesh, region):
    return mesh.vertices[np.unique(mesh.faces[region.faces])]


def _bottom_chamfer(mesh, regions, outline, z0):
    """Leg of a 45 degree chamfer on the base bottom, and which edges carry it.

    A sharp edge of the extruded outline sits on the mesh. A chamfered one
    leaves the corner off the surface by about ``leg / sqrt(2)``.
    """
    legs = []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 0.05:
            continue
        got = _plane(region)
        if got is None:
            continue
        n, _d = got
        if abs(abs(n[2]) - 0.707) > 0.08:
            continue
        verts = _face_points(mesh, region)
        if float(verts[:, 2].min()) > z0 + 0.3 or float(verts[:, 2].max()) > z0 + 1.2:
            continue
        legs.append(float(verts[:, 2].max() - verts[:, 2].min()))
    if not legs:
        return 0.0, []
    leg = float(np.median(legs))
    if leg < 0.05:
        return 0.0, []
    seq = outline
    mids = []
    keep = []
    for i, a in enumerate(seq):
        b = seq[(i + 1) % len(seq)]
        if math.hypot(b[0] - a[0], b[1] - a[1]) < 0.02:
            continue
        mids.append(((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, z0))
        keep.append(i)
    if not mids:
        return 0.0, []
    gap = np.abs(trimesh.proximity.closest_point(mesh, np.asarray(mids, float))[1])
    thresh = max(0.08, 0.35 * leg)
    chosen = [i for i, g in zip(keep, gap) if g > thresh]
    if not chosen:
        return 0.0, []
    return leg, chosen


def _upper_faces(mesh, regions, z0):
    roofs, bevels, vertical = [], [], []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 8.0:
            continue
        got = _plane(region)
        if got is None:
            continue
        n, d = got
        verts = _face_points(mesh, region)
        zmin, zmax = float(verts[:, 2].min()), float(verts[:, 2].max())
        if zmax < z0 + 0.4:
            continue
        aligned = int(np.argmax(np.abs(n)))
        if abs(n[aligned]) > 0.98 and aligned != 2:
            if zmin >= z0 - 0.3:
                vertical.append((region, n, d, verts))
            continue
        extent = verts.max(0) - verts.min(0)
        thin = float(np.min(extent)) < 3.0
        comps = sorted((abs(float(n[0])), abs(float(n[1])), abs(float(n[2]))), reverse=True)
        bevelish = abs(comps[0] - 0.707) < 0.08 and abs(comps[1] - 0.707) < 0.08 and comps[2] < 0.15
        if thin and bevelish:
            bevels.append((region, n, d, verts))
        elif (not thin) and abs(n[2]) > 0.25 and zmin >= z0 - 0.5:
            roofs.append((region, n, d, verts))
    return roofs, bevels, vertical


def _tower_box(mesh, regions, z0, z1):
    """Axis-aligned block standing on the shelf, from the top cap plus its walls."""
    caps = []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 40.0:
            continue
        got = _plane(region)
        if got is None or abs(got[0][2]) < 0.97:
            continue
        z = _plane_coord(region.surface, 2)
        if abs(z - z1) > 0.4:
            continue
        caps.append(region)
    if not caps:
        return None
    cap = max(caps, key=lambda r: r.area)
    c0, c1 = _bbox(_face_points(mesh, cap))
    margin = 6.0
    box0 = np.array([c0[0] - margin, c0[1] - margin, z0])
    box1 = np.array([c1[0] + margin, c1[1] + margin, z1])
    acc0, acc1 = c0.copy(), c1.copy()
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 15.0:
            continue
        got = _plane(region)
        if got is None:
            continue
        n, _d = got
        if abs(n[2]) > 0.2 or max(abs(n[0]), abs(n[1])) < 0.98:
            continue
        verts = _face_points(mesh, region)
        if float(verts[:, 2].min()) < z0 - 0.3:
            continue
        # ignore vertices that run off along a neighbouring roof; they are not
        # part of this block
        near = verts[np.all((verts[:, :2] >= c0[:2] - margin) & (verts[:, :2] <= c1[:2] + margin), axis=1)]
        if len(near) < 3:
            continue
        center = near.mean(0)
        if not (box0[0] <= center[0] <= box1[0] and box0[1] <= center[1] <= box1[1]):
            continue
        acc0 = np.minimum(acc0, near.min(0))
        acc1 = np.maximum(acc1, near.max(0))
    if acc1[0] - acc0[0] < 1.0 or acc1[1] - acc0[1] < 1.0:
        return None
    return (float(acc0[0]), float(acc1[0]), float(acc0[1]), float(acc1[1]))


def _clip_dict(n, d):
    return {"x": float(n[0]), "y": float(n[1]), "z": float(n[2]), "d": float(d)}


def _bevel_prism(n, d, verts, miters):
    """Triangle prism of one edge bevel, clipped by its plane and its miters."""
    axes = [i for i in range(3) if abs(n[i]) > 0.3]
    if len(axes) != 2:
        return None
    long_axis = [i for i in range(3) if i not in axes][0]
    u, v = axes
    uv = verts[:, [u, v]]
    u0, v0 = uv.min(0)
    u1, v1 = uv.max(0)
    corners = [(u0, v0), (u0, v1), (u1, v0), (u1, v1)]

    def side(p):
        q = np.zeros(3)
        q[u], q[v] = p
        return float(n @ q)

    inward = min(corners, key=side)
    others = [p for p in corners if p != inward and p != (u0 if inward[0] == u1 else u1, v0 if inward[1] == v1 else v1)]
    # the two corners adjacent to the inward one, not the opposite corner
    adj = []
    for p in corners:
        if p == inward:
            continue
        if p[0] == inward[0] or p[1] == inward[1]:
            adj.append(p)
    if len(adj) != 2:
        return None
    # push the right angle into the body
    du = -1.0 if n[u] > 0 else 1.0
    dv = -1.0 if n[v] > 0 else 1.0
    corner = (inward[0] + du * _BEVEL_IN, inward[1] + dv * _BEVEL_IN)
    outline = [corner, adj[0], adj[1]]
    # profile stored for the long axis: X stores (y, z), Y stores (x, z), Z stores (x, y)
    order = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[long_axis]
    prof = []
    for p in outline:
        xy = {u: p[0], v: p[1]}
        prof.append((float(xy[order[0]]), float(xy[order[1]])))
    if _shoelace(prof) < 0:
        prof = list(reversed(prof))
    lo = float(verts[:, long_axis].min())
    hi = float(verts[:, long_axis].max())
    clips = [_clip_dict(n, d)]
    center = verts.mean(0)
    for mn, md in miters:
        if abs(float(mn @ n)) > 0.98:
            continue
        on = np.abs(verts @ mn - md)
        if int(np.sum(on < 0.1)) < 2:
            continue
        if float(mn @ center) > md + 0.2:
            continue
        clips.append(_clip_dict(mn, md))
    return {
        "axis": "XYZ"[long_axis],
        "axis_i": long_axis,
        "origin": lo,
        "length": hi - lo,
        "outline": prof,
        "cutouts": [],
        "clips": clips,
        "grow_lo": 0.0,
        "grow_hi": 0.0,
    }


def _roof_prism(mesh, region, n, d, verts, z_shelf):
    rings = _rings(mesh, region)
    if not rings:
        return None
    outline = _simplify(rings[0][:, :2], 0.05)
    if len(outline) < 3:
        return None
    z0 = float(verts[:, 2].min())
    z1 = float(verts[:, 2].max())
    if abs(z0 - z_shelf) < 0.6:
        z0 = float(z_shelf)
    if z1 - z0 < 0.4:
        return None
    return {
        "axis": "Z",
        "axis_i": 2,
        "origin": z0,
        "length": z1 - z0,
        "outline": outline,
        "cutouts": [],
        "clips": [_clip_dict(n, d)],
        "grow_lo": _GROW,
        "grow_hi": 0.0,
    }


def _slot_cuts(mesh, regions):
    walls = []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 4.0:
            continue
        got = _plane(region)
        if got is None:
            continue
        n, d = got
        ax = int(np.argmax(np.abs(n)))
        if abs(n[ax]) < 0.98:
            continue
        verts = _face_points(mesh, region)
        coord = float(d / n[ax])
        walls.append((ax, coord, verts.min(0), verts.max(0), region.area))
    cuts = []
    used = set()
    for i, (ax, coord, b0, b1, _area) in enumerate(walls):
        best = None
        for j, (jax, jcoord, c0, c1, _ja) in enumerate(walls):
            if i >= j or ax != jax:
                continue
            gap = abs(coord - jcoord)
            if gap < 0.25 or gap > 1.2:
                continue
            # a slot is two faces of one rib, not a small floor against a large cap
            if max(_area, _ja) > 2.5 * max(min(_area, _ja), 1e-6):
                continue
            lo = np.maximum(b0, c0)
            hi = np.minimum(b1, c1)
            other = [k for k in range(3) if k != ax]
            if any(hi[k] - lo[k] < 1.5 for k in other):
                continue
            if best is None or gap < best[0]:
                best = (gap, j, lo, hi, coord, jcoord)
        if best is None:
            continue
        _gap, j, lo, hi, d, jd = best
        key = (i, j) if i < j else (j, i)
        if key in used:
            continue
        used.add(key)
        mid = 0.5 * (lo + hi)
        mid[ax] = 0.5 * (d + jd)
        if bool(mesh.contains([mid])[0]):
            continue
        # profile is the two axes that are not the thin one, in outline order
        order = {0: (1, 2), 1: (0, 2), 2: (0, 1)}[ax]
        p0, p1 = order
        outline = [
            (float(lo[p0]), float(lo[p1])),
            (float(hi[p0]), float(lo[p1])),
            (float(hi[p0]), float(hi[p1])),
            (float(lo[p0]), float(hi[p1])),
        ]
        origin = float(min(d, jd))
        cuts.append({
            "axis": "XYZ"[ax],
            "axis_i": ax,
            "origin": origin,
            "length": float(abs(d - jd)),
            "outline": outline,
            "cutouts": [],
        })
    return cuts


def _as_circle(pts, tol):
    p = np.asarray(pts, float)
    if len(p) < 6:
        return None
    c = p.mean(0)
    rho = np.linalg.norm(p - c, axis=1)
    radius = float(np.median(rho))
    if radius < 0.3 or float(np.max(np.abs(rho - radius))) > max(0.25, 4.0 * tol):
        return None
    return radius


def _pocket_mouth(mesh, regions, floor, ring, axis, coord):
    """Opening coordinate from the side walls that share the floor's outer ring.

    A rib standing inside the pocket is a closer parallel plane than the real
    mouth. The walls around the floor end at the mouth, and a wall of the
    surrounding body is much longer, so the short walls name the opening.
    """
    ends = []
    step = max(1, len(ring) // 12)
    for other in regions:
        if other is floor or not isinstance(other.surface, Plane) or other.area < 1.0:
            continue
        got = _plane(other)
        if got is None or abs(got[0][axis]) > 0.35:
            continue
        pts = _face_points(mesh, other)
        if len(pts) == 0:
            continue
        sample = pts if len(pts) <= 40 else pts[:: max(1, len(pts) // 40)]
        near = False
        for p in ring[::step]:
            if float(np.linalg.norm(sample - p, axis=1).min()) <= 0.12:
                near = True
                break
        if not near:
            continue
        lo = float(pts[:, axis].min())
        hi = float(pts[:, axis].max())
        if hi - lo < 0.3 or hi - lo > 8.0:
            continue
        if coord < lo - 0.15 or coord > hi + 0.15:
            continue
        ends.append(hi if abs(lo - coord) <= abs(hi - coord) else lo)
    if len(ends) < 3:
        return None
    return float(np.median(ends))


def _pocket_cuts(mesh, regions, tol):
    """Blind pockets from the floor face out to the mouth.

    The floor is a modest cap with a hole or an island in it. Circular holes
    are left to the hole features; any other inner ring is an island that stays.
    """
    cuts = []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 8.0 or region.area > 80.0:
            continue
        got = _plane(region)
        if got is None:
            continue
        n, _d = got
        ax = int(np.argmax(np.abs(n)))
        if abs(n[ax]) < 0.97:
            continue
        rings = _rings(mesh, region)
        if len(rings) < 2:
            continue
        uv = [i for i in range(3) if i != ax]
        outer = rings[0][:, uv]
        if _as_circle(outer, tol) is not None or _loop_area(outer) < 8.0:
            continue
        coord = _plane_coord(region.surface, ax)
        mouth = _pocket_mouth(mesh, regions, region, rings[0], ax, coord)
        if mouth is None:
            continue
        depth = abs(float(mouth) - float(coord))
        if depth < 0.3 or depth > 12.0:
            continue
        islands = []
        for ring in rings[1:]:
            pts = ring[:, uv]
            if _as_circle(pts, tol) is not None or _loop_area(pts) < 0.5:
                continue
            simplified = _simplify(pts, max(tol, 0.04))
            if len(simplified) >= 3:
                islands.append(simplified)
        outline = _simplify(outer, max(tol, 0.04))
        if len(outline) < 3:
            continue
        # break the skin at the mouth and stop on the floor
        if mouth < coord:
            origin, length = float(mouth) - 0.02, depth + 0.02
        else:
            origin, length = float(coord), depth + 0.02
        cuts.append({
            "axis": "XYZ"[ax],
            "axis_i": ax,
            "origin": origin,
            "length": float(length),
            "outline": outline,
            "cutouts": islands,
        })
    return cuts


def _measured_span(mesh, regions, hole):
    ax = "XYZ".index(hole["axis"])
    point = np.zeros(3)
    for i, coord in enumerate("xyz"):
        if coord in hole:
            point[i] = float(hole[coord])
    along = []
    for region in regions:
        surface = region.surface
        if isinstance(surface, Cylinder):
            if _principal_axis(surface.axis) != ax:
                continue
            if _perp(ax, surface.point, point) > 0.8:
                continue
        elif isinstance(surface, Cone):
            if _principal_axis(surface.axis) != ax:
                continue
            if _perp(ax, surface.apex, point) > 0.8:
                continue
        else:
            continue
        verts = _face_points(mesh, region)
        along.append(verts[:, ax])
    if not along:
        return None
    coords = np.concatenate(along)
    lo, hi = float(coords.min()), float(coords.max())
    # snap outward onto a horizontal (or axial) cap the hole actually meets
    caps = []
    for region in regions:
        if not isinstance(region.surface, Plane):
            continue
        got = _plane(region)
        if got is None or abs(got[0][ax]) < 0.97 or region.area < 20.0:
            continue
        caps.append(_plane_coord(region.surface, ax))
    part_lo = float(mesh.vertices[:, ax].min())
    part_hi = float(mesh.vertices[:, ax].max())
    caps.extend([part_lo, part_hi])
    for cap in caps:
        if lo - 0.55 <= cap <= lo:
            lo = float(cap)
        if hi <= cap <= hi + 0.55:
            hi = float(cap)
    if hi - lo < 0.3:
        return None
    return lo, hi


def _perp(axis, a, b):
    aa = np.asarray(a, float) - np.asarray(b, float)
    aa[axis] = 0.0
    return float(np.linalg.norm(aa))


def _round_loop(loop, tol):
    return [(float(_nice(p[0], tol)), float(_nice(p[1], tol))) for p in loop]


def _round_cut(spec, tol):
    spec = dict(spec)
    spec["origin"] = float(_nice(spec["origin"], tol))
    spec["length"] = float(_nice(max(spec["length"], 0.0), tol))
    spec["outline"] = _round_loop(spec["outline"], tol)
    spec["cutouts"] = [_round_loop(c, tol) for c in spec.get("cutouts") or []]
    return spec


def _overlap_stacked(positives):
    """Overlap Z extrusions that meet, so the fuse is one solid.

    The named length stays the true length. The builder grows the joint only.
    """
    ordered = sorted((p for p in positives if p["axis"] == "Z"), key=lambda p: p["origin"])
    for lower, upper in zip(ordered, ordered[1:]):
        if abs((lower["origin"] + lower["length"]) - upper["origin"]) > 0.25:
            continue
        lower["grow_hi"] = max(float(lower.get("grow_hi") or 0.0), _GROW)
        upper["grow_lo"] = max(float(upper.get("grow_lo") or 0.0), _GROW)


def extract_stack(mesh, regions, holes, tol):
    """A union of partial extrusions, or None when this vocabulary does not apply."""
    shelves = _shelves(mesh, regions)
    if len(shelves) < 2:
        return None
    positives = []
    for (z0, _a0), (z1, _a1) in zip(shelves, shelves[1:]):
        outline = _stable_outline(mesh, z0, z1, tol)
        if outline is None or len(outline) < 3:
            continue
        positives.append({
            "axis": "Z",
            "axis_i": 2,
            "origin": float(z0),
            "length": float(z1 - z0),
            "outline": outline,
            "cutouts": [],
            "clips": [],
            "grow_lo": 0.0,
            "grow_hi": 0.0,
        })
    if not positives:
        return None
    # the lowest band carries a bottom chamfer when the mesh has one
    base = min(positives, key=lambda p: p["origin"])
    leg, edges = _bottom_chamfer(mesh, regions, base["outline"], base["origin"])
    if edges:
        base["chamfer"] = float(leg)
        base["chamfer_at"] = edges

    top = max(p["origin"] + p["length"] for p in positives)
    # a band that was not a straight extrusion may still be roofs, a block, bevels
    upper = None
    for z0, z1 in zip([s[0] for s in shelves], [s[0] for s in shelves[1:]]):
        if any(abs(p["origin"] - z0) < 0.2 and abs(p["origin"] + p["length"] - z1) < 0.2 for p in positives):
            continue
        if z1 <= top + 0.2 and z0 >= min(p["origin"] for p in positives) - 0.2:
            upper = (z0, z1)
            break
    # also the range above the tallest stable band, up to the highest shelf
    if upper is None and shelves[-1][0] > top + 0.4:
        upper = (top, shelves[-1][0])
    if upper is not None:
        z0, z1 = upper
        roofs, bevels, _vertical = _upper_faces(mesh, regions, z0)
        box = _tower_box(mesh, regions, z0, z1)
        tower = None
        if box is not None:
            x0, x1, y0, y1 = box
            clips = []
            for region, n, d, verts in roofs:
                xy0, xy1 = verts[:, :2].min(0), verts[:, :2].max(0)
                inside = xy0[0] >= x0 - 0.5 and xy1[0] <= x1 + 0.5 and xy0[1] >= y0 - 0.5 and xy1[1] <= y1 + 0.5
                if inside:
                    clips.append(_clip_dict(n, d))
            tower = {
                "axis": "Z",
                "axis_i": 2,
                "origin": float(z0),
                "length": float(z1 - z0),
                "outline": [(x0, y0), (x1, y0), (x1, y1), (x0, y1)],
                "cutouts": [],
                "clips": clips,
                "grow_lo": _GROW,
                "grow_hi": 0.0,
            }
            positives.append(tower)
        miters = [(n, d) for _r, n, d, _v in roofs + bevels]
        # vertical 45 degree corners are miters too, not extra prisms
        for region in regions:
            if not isinstance(region.surface, Plane) or region.area < 8.0:
                continue
            got = _plane(region)
            if got is None:
                continue
            n, d = got
            if abs(n[2]) > 0.2:
                continue
            if abs(abs(n[0]) - 0.707) < 0.08 and abs(abs(n[1]) - 0.707) < 0.08:
                miters.append((n, d))
        for region, n, d, verts in bevels:
            prism = _bevel_prism(n, d, verts, miters)
            if prism is not None:
                positives.append(prism)
        for region, n, d, verts in roofs:
            if tower is not None and any(c["d"] == d and abs(c["x"] - n[0]) < 1e-9 for c in tower["clips"]):
                continue
            prism = _roof_prism(mesh, region, n, d, verts, z0)
            if prism is not None:
                positives.append(prism)
        # grow the shelf under this upper work so the joint fuses
        for p in positives:
            if p["axis"] == "Z" and abs(p["origin"] + p["length"] - z0) < 0.2:
                p["grow_hi"] = _GROW

    _overlap_stacked(positives)
    cuts = _slot_cuts(mesh, regions) + _pocket_cuts(mesh, regions, tol)
    out_holes = []
    for hole in holes:
        h = dict(hole)
        span = _measured_span(mesh, regions, h)
        if span is not None:
            h["span0"], h["span1"] = span
        out_holes.append(h)

    # round to the fit tolerance, but keep the arc. 0.02 mm is enough to name
    # a coordinate without moving a faceted round.
    nice = max(min(tol, 0.02), 1e-4)
    rounded = []
    for p in positives:
        q = _round_cut(p, nice)
        q["axis_i"] = p["axis_i"]
        q["clips"] = []
        for c in p.get("clips") or []:
            q["clips"].append({
                "x": float(_nice(c["x"], 1e-4)),
                "y": float(_nice(c["y"], 1e-4)),
                "z": float(_nice(c["z"], 1e-4)),
                "d": float(_nice(c["d"], 1e-3)),
            })
        if p.get("chamfer_at"):
            q["chamfer"] = float(_nice(p["chamfer"], nice))
            q["chamfer_at"] = list(p["chamfer_at"])
        q["grow_lo"] = float(p.get("grow_lo") or 0.0)
        q["grow_hi"] = float(p.get("grow_hi") or 0.0)
        rounded.append(q)
    cuts = [_round_cut(c, nice) for c in cuts]
    for h in out_holes:
        if "span0" in h:
            h["span0"] = float(_nice(h["span0"], nice))
            h["span1"] = float(_nice(h["span1"], nice))
    return {"positive": rounded, "cuts": cuts, "holes": out_holes, "combine": "union"}
