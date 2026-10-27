"""Fallback outline and chamfer cuts for a combined extrusion.

Used only after the silhouette intersection has already missed the mesh.
Two things the silhouettes miss:

* A hole that breaks out through a slot becomes part of the outer loop, so
  the fused cap follows the cylinder instead of the walls around it. The
  outline is rebuilt from a raster of the mesh: edges snap to vertical wall
  planes, and a run that rides a hole cylinder is closed by the chord of
  that hole. A flare that sits outside the cylinder is left on the loop.
* A bevel that does not run the full length of the part is not in any
  silhouette. It becomes a triangular cut extrusion along the edge, spanning
  the wide part of the face. A 45 degree cone that was not already used as a
  hole chamfer is cut as its own cone.

Neither change is applied to a part the silhouette boolean already accepts.
"""
from __future__ import annotations

import math

import numpy as np

from .csg import (
    _PROFILE_AXES, _cone_ends, _points_inside, _raster_outline, _shoelace,
)
from .fit import Cone, Cylinder, Plane
from .parametric import _principal_axis


def _unit(v):
    v = np.asarray(v, float)
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-12 else v


def _profile_of(point, axis):
    u, v = _PROFILE_AXES[axis]
    p = np.asarray(point, float)
    return (float(p @ u), float(p @ v))


def _walls(regions, axis, min_area):
    """Vertical wall lines in the profile plane: ``(nu, nv, d, key)``."""
    u, v = _PROFILE_AXES[axis]
    found = []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < min_area:
            continue
        normal = _unit(region.surface.normal)
        if abs(float(normal[axis])) > 0.08:
            continue
        nu, nv = float(normal @ u), float(normal @ v)
        scale = math.hypot(nu, nv)
        if scale < 0.5:
            continue
        nu, nv = nu / scale, nv / scale
        d = float(region.surface.d) / max(float(np.linalg.norm(region.surface.normal)), 1e-12)
        d = d / scale
        key = (round(nu, 3), round(nv, 3), round(d, 2))
        if any(w[3] == key for w in found):
            continue
        found.append((nu, nv, d, key))
    return found


def _hole_circles(holes, axis):
    """``(centre_u, centre_v, radius)`` for holes along ``axis``, including counterbores."""
    circles = []
    for hole in holes:
        if hole.get("axis") != "XYZ"[axis]:
            continue
        point = [float(hole.get(c, 0.0)) for c in "xyz"]
        cu, cv = _profile_of(point, axis)
        radii = [float(hole["d"]) / 2.0]
        if "cbore_d" in hole:
            radii.append(float(hole["cbore_d"]) / 2.0)
        for radius in radii:
            circles.append((cu, cv, radius))
    return circles


def _near_circle(point, circles, gap):
    for cu, cv, radius in circles:
        if abs(math.hypot(point[0] - cu, point[1] - cv) - radius) <= gap:
            return (cu, cv, radius)
    return None


def _line_circle(wall, centre, radius):
    nu, nv, dd = wall[0], wall[1], wall[2]
    dist = dd - (nu * centre[0] + nv * centre[1])
    half = radius * radius - dist * dist
    if half < 0:
        return []
    root = math.sqrt(half)
    fx, fy = centre[0] + nu * dist, centre[1] + nv * dist
    tx, ty = -nv, nu
    return [(fx + tx * root, fy + ty * root), (fx - tx * root, fy - ty * root)]


def _intersect(w1, w2):
    det = w1[0] * w2[1] - w1[1] * w2[0]
    if abs(det) < 1e-8:
        return None
    return (
        (w2[1] * w1[2] - w1[1] * w2[2]) / det,
        (w1[0] * w2[2] - w2[0] * w1[2]) / det,
    )


def _same_wall(a, b):
    return a is not None and b is not None and a[3] == b[3]


def _on_line(point, wall, gap):
    return abs(wall[0] * point[0] + wall[1] * point[1] - wall[2]) <= gap


def _coalesce_wall_runs(runs, pts, n):
    """Join pieces of one wall that a raster staircase split apart.

    A marching-squares step is a short unassigned edge whose vertices still
    lie on the wall. It is not a notch. A gap that leaves the wall is kept,
    so a slot or a flare stays on the loop.
    """
    if not runs:
        return runs

    def gap_stays(wall, start, end):
        for t in range(start, end + 1):
            if not _on_line(pts[t % n], wall, 0.55):
                return False
        return True

    merged = []
    i = 0
    count = len(runs)
    while i < count:
        wall, a, b = runs[i]
        if wall is None:
            merged.append(runs[i])
            i += 1
            continue
        end_b = b
        j = i + 1
        while j < count:
            wall2, a2, b2 = runs[j]
            if _same_wall(wall2, wall):
                end_b = b2
                j += 1
                continue
            if wall2 is not None:
                break
            k = j
            while k < count and runs[k][0] is None:
                k += 1
            if k >= count or not _same_wall(runs[k][0], wall):
                break
            if not gap_stays(wall, a2, runs[k - 1][2]):
                break
            j = k
        merged.append((wall, a, end_b))
        i = j
    if len(merged) < 2 or merged[-1][0] is None:
        return merged
    wall = merged[-1][0]
    front = 0
    while front < len(merged) - 1 and merged[front][0] is None:
        front += 1
    if front == 0 and _same_wall(merged[0][0], wall):
        head = merged[0]
        return [(wall, merged[-1][1], head[2])] + merged[1:-1]
    if front and _same_wall(merged[front][0], wall) and gap_stays(wall, merged[0][1], merged[front - 1][2]):
        return [(wall, merged[-1][1], merged[front][2])] + merged[front + 1:]
    return merged


def wall_outline(mesh, regions, holes, axis, tol):
    """Outer profile snapped to wall planes, with a hole breakout closed.

    ``None`` when the raster has no loop or the snapped loop is far from it.
    """
    pitch = max(0.1, 2.0 * tol)
    loops = _raster_outline(mesh, axis, pitch, max(0.05, tol), fill=True)
    if not loops:
        return None
    raw = [(float(p[0]), float(p[1])) for p in loops[0]]
    if len(raw) < 4:
        return None
    raw_area = abs(_shoelace(raw))
    circles = _hole_circles(holes, axis)
    walls = _walls(regions, axis, min_area=max(1.2, 4.0 * tol * tol))
    if len(walls) < 3:
        return None

    def on_wall(point, gap):
        return _near_circle(point, circles, gap) is not None

    n = len(raw)
    flags = [on_wall(p, 0.1) for p in raw]
    pts = [raw[i] for i in range(n) if not (flags[i] and flags[i - 1] and flags[(i + 1) % n])]
    if len(pts) < 4:
        return None

    sign = 1.0 if _shoelace(pts) > 0 else -1.0

    def assign(a, b):
        dx, dy = b[0] - a[0], b[1] - a[1]
        length = math.hypot(dx, dy)
        if length < 0.3:
            return None
        # a chord across a hole is the breakout, not a wall
        if on_wall(a, 0.35) and on_wall(b, 0.35):
            return None
        ox, oy = sign * dy / length, sign * -dx / length
        best = None
        for wall in walls:
            para = abs(wall[0] * dx + wall[1] * dy) / length
            if para > 0.4:
                continue
            da = abs(wall[0] * a[0] + wall[1] * a[1] - wall[2])
            db = abs(wall[0] * b[0] + wall[1] * b[1] - wall[2])
            if max(da, db) > 0.45:
                continue
            if wall[0] * ox + wall[1] * oy < 0.35:
                continue
            score = max(da, db) + 0.5 * para
            if best is None or score < best[0]:
                best = (score, wall)
        return None if best is None else best[1]

    n = len(pts)
    edges = [assign(pts[i], pts[(i + 1) % n]) for i in range(n)]

    def same(a, b):
        return a is not None and b is not None and a[3] == b[3]

    k = next((i for i in range(n) if not same(edges[i - 1], edges[i])), 0)
    pts = pts[k:] + pts[:k]
    edges = edges[k:] + edges[:k]
    runs = []
    i = 0
    while i < n:
        if edges[i] is None:
            runs.append((None, i, i + 1))
            i += 1
            continue
        j = i + 1
        while j < n and same(edges[j], edges[i]):
            j += 1
        runs.append((edges[i], i, j))
        i = j

    def near_line(point, wall):
        return abs(wall[0] * point[0] + wall[1] * point[1] - wall[2]) <= 0.4

    cleaned = []
    for idx, (wall, a, b) in enumerate(runs):
        if wall is None or b - a > 1:
            cleaned.append((wall, a, b))
            continue
        length = math.hypot(pts[b % n][0] - pts[a][0], pts[b % n][1] - pts[a][1])
        ends = (pts[a], pts[b % n])
        prev_w = runs[idx - 1][0]
        next_w = runs[(idx + 1) % len(runs)][0]
        if length < 0.9 and prev_w is not None and all(near_line(p, prev_w) for p in ends):
            continue
        if length < 0.9 and next_w is not None and all(near_line(p, next_w) for p in ends):
            continue
        cleaned.append((wall, a, b))
    runs = _coalesce_wall_runs(cleaned, pts, n)

    verts = []

    def push(point):
        if not verts or math.hypot(verts[-1][0] - point[0], verts[-1][1] - point[1]) > 0.08:
            verts.append(point)

    i = 0
    while i < len(runs):
        wall, a, b = runs[i]
        if wall is None:
            j = i
            while j < len(runs) and runs[j][0] is None:
                j += 1
            prev_w = runs[i - 1][0]
            next_w = runs[j][0] if j < len(runs) else None
            start_p = pts[a]
            end_i = runs[j - 1][2]
            end_p = pts[end_i % n]
            dist = math.hypot(end_p[0] - start_p[0], end_p[1] - start_p[1])
            if prev_w is not None and next_w is not None and dist < 3.5 and not same(prev_w, next_w):
                hit_a = _near_circle(start_p, circles, 0.12)
                hit_b = _near_circle(end_p, circles, 0.12)
                corner = _intersect(prev_w, next_w)
                if corner is None and hit_a is not None and hit_b is not None:
                    for side, anchor, hit in ((prev_w, start_p, hit_a), (next_w, end_p, hit_b)):
                        choices = _line_circle(side, (hit[0], hit[1]), hit[2])
                        if not choices:
                            push(anchor)
                        else:
                            push(min(choices, key=lambda t: math.hypot(t[0] - anchor[0], t[1] - anchor[1])))
                    i = j
                    continue
                if corner is not None and math.hypot(corner[0] - start_p[0], corner[1] - start_p[1]) < 3.0:
                    push(corner)
                    i = j
                    continue
            for t in range(a, end_i + 1):
                push(pts[t % n])
            i = j
            continue
        nxt = runs[(i + 1) % len(runs)]
        if nxt[0] is not None and not same(wall, nxt[0]):
            shared = pts[b % n]
            corner = _intersect(wall, nxt[0])
            if corner is not None and math.hypot(corner[0] - shared[0], corner[1] - shared[1]) < 1.8:
                push(corner)
            else:
                push(shared)
        i += 1

    if len(verts) > 2 and math.hypot(verts[0][0] - verts[-1][0], verts[0][1] - verts[-1][1]) < 0.15:
        verts.pop()
    if len(verts) < 3:
        return None
    area = abs(_shoelace(verts))
    if raw_area <= 1e-6 or area < 0.7 * raw_area or area > 1.45 * raw_area:
        return None
    return verts


def _radial_facet(centre, normal, regions, holes):
    """True when the face is a chord of a hole or a fitted cylinder.

    A bevel next to a hole sits near the same radius, but its normal is not
    the radius, so it is kept.
    """
    point = np.asarray(centre, float)
    direction = _unit(normal)

    def radial(axis, origin, radius, limit):
        axis_u = _unit(axis)
        rel = point - np.asarray(origin, float)
        rel = rel - axis_u * float(rel @ axis_u)
        rho = float(np.linalg.norm(rel))
        if rho < 1e-6 or abs(rho - radius) > limit:
            return False
        return abs(float(direction @ (rel / rho))) >= 0.92

    for region in regions:
        if isinstance(region.surface, Cylinder) and radial(
            region.surface.axis, region.surface.point, float(region.surface.radius), 0.2
        ):
            return True
    for hole in holes:
        hax = "XYZ".index(hole["axis"])
        axis = np.zeros(3)
        axis[hax] = 1.0
        origin = np.zeros(3)
        for i, coord in enumerate("xyz"):
            if coord in hole and i != hax:
                origin[i] = float(hole[coord])
        radii = [float(hole["d"]) / 2.0]
        if "cbore_d" in hole:
            radii.append(float(hole["cbore_d"]) / 2.0)
        if any(radial(axis, origin, radius, 0.25) for radius in radii):
            return True
    return False


def _wide_span(verts, edge, u, v):
    """Edge interval where the face is a real strip, not a single tip vertex.

    A chamfer that runs into a corner ends on one point. Extruding the full
    triangle out to that point cuts material past the face.
    """
    coords = verts[:, edge]
    order = np.argsort(coords)
    groups = []
    for idx in order:
        t = float(coords[idx])
        if not groups or t - groups[-1][1] > 0.12:
            groups.append([t, t, [int(idx)]])
        else:
            groups[-1][1] = t
            groups[-1][2].append(int(idx))
    prof = [(float(verts[i, u]), float(verts[i, v])) for i in range(len(verts))]
    full = 0.0
    for a in prof:
        for b in prof:
            full = max(full, math.hypot(a[0] - b[0], a[1] - b[1]))
    if full < 0.15:
        return None
    wide = []
    for start, end, idxs in groups:
        width = 0.0
        for i in idxs:
            for j in idxs:
                width = max(width, math.hypot(prof[i][0] - prof[j][0], prof[i][1] - prof[j][1]))
        if width >= 0.45 * full:
            wide.append((start, end))
    lo, hi = float(coords.min()), float(coords.max())
    if not wide:
        return lo, hi
    start, end = wide[0][0], wide[-1][1]
    # A dense cluster in the middle of a long straight bevel is where the
    # face meets another edge, not the only place the bevel exists.
    if end - start < 0.55 * (hi - lo):
        return lo, hi
    return start, end


def chamfer_cuts(mesh, regions, holes):
    """Planar bevels as triangular prisms along the edge they follow.

    The section is the right triangle of the face's bounding box, so the two
    walls it joins have to be axis-aligned. A face that is a chord of a hole
    cylinder is left to the hole cutter. A bevel that tapers to a point is
    extruded only where the strip is wide: the tip is one vertex, and the
    full triangle does not exist there.
    """
    found = []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 0.4:
            continue
        normal = _unit(region.surface.normal)
        edge = int(np.argmin(np.abs(normal)))
        if abs(normal[edge]) > 0.18:
            continue
        others = [i for i in range(3) if i != edge]
        # 45 degrees is 0.707. A bevel between about 27 and 63 degrees still
        # has a right-triangle section when the two walls it joins are
        # axis-aligned. Shallower than that is a draft, not this cut.
        if not all(abs(normal[i]) >= 0.45 for i in others):
            continue
        weights = mesh.area_faces[region.faces]
        centres = mesh.triangles_center[region.faces]
        centre = (centres * weights[:, None]).sum(0) / max(float(weights.sum()), 1e-9)
        if _radial_facet(centre, normal, regions, holes):
            continue
        verts = mesh.vertices[np.unique(mesh.faces[region.faces])]
        u, v = others
        spanned = _wide_span(verts, edge, u, v)
        if spanned is None:
            continue
        origin, end = spanned
        span = end - origin
        if span < 0.4:
            continue
        du = float(verts[:, u].max() - verts[:, u].min())
        dv = float(verts[:, v].max() - verts[:, v].min())
        if min(du, dv) < 0.15 or max(du, dv) > 8.0 or max(du, dv) / min(du, dv) > 1.8:
            continue
        trial = centre - normal * 0.3
        inside = _points_inside(mesh, trial.reshape(1, 3))
        if len(inside) and not bool(inside[0]):
            normal = -normal
        u0, u1 = float(verts[:, u].min()), float(verts[:, u].max())
        v0, v1 = float(verts[:, v].min()), float(verts[:, v].max())
        cu = u1 if normal[u] > 0 else u0
        cv = v1 if normal[v] > 0 else v0
        tri = [(cu, cv), (u0 if cu == u1 else u1, cv), (cu, v0 if cv == v1 else v1)]
        mid = np.zeros(3)
        mid[edge] = 0.5 * (origin + end)
        mid[u] = sum(p[0] for p in tri) / 3.0
        mid[v] = sum(p[1] for p in tri) / 3.0
        inn = _points_inside(mesh, mid.reshape(1, 3))
        if len(inn) and bool(inn[0]):
            continue
        prof = _PROFILE_AXES[edge]
        pax = [int(np.argmax(np.abs(prof[0]))), int(np.argmax(np.abs(prof[1])))]
        outline = []
        for pu, pv in tri:
            coord = {u: pu, v: pv}
            outline.append((coord[pax[0]], coord[pax[1]]))
        if abs(_shoelace(outline)) < 1e-4:
            continue
        found.append({
            "axis": "XYZ"[edge],
            "axis_i": edge,
            "origin": float(origin),
            "length": float(span),
            "outline": outline,
            "cutouts": [],
        })
    return found


def cone_cuts(mesh, regions):
    """45 degree cones as finite cone cutters. A hole chamfer may be repeated."""
    found = []
    for region in regions:
        if not isinstance(region.surface, Cone) or region.area < 0.3:
            continue
        if abs(float(region.surface.half_angle) - math.pi / 4.0) > 0.08:
            continue
        axis = _principal_axis(region.surface.axis)
        if axis is None:
            continue
        lo, hi, r0, r1 = _cone_ends(mesh, region, axis)
        if hi - lo < 0.08 or min(r0, r1) < 0.05:
            continue
        origin = np.asarray(region.surface.apex, float).copy()
        origin[axis] = lo
        found.append({
            "axis": "XYZ"[axis],
            "origin": [float(origin[0]), float(origin[1]), float(origin[2])],
            "r0": float(r0),
            "r1": float(r1),
            "length": float(hi - lo),
        })
    return found
