"""Straight extrusions of a 2D outline, for parametric CadQuery and OpenSCAD.

A part qualifies when it is a prism along X, Y, or Z: two end caps, side
walls parallel to that axis, optional polygonal through-openings in the
outline, and optional round through-holes. A through-hole along the
extrusion may be counterbored. On a Z extrusion a vertical through-hole may
also be countersunk, opening on the top or the bottom face (``hole_N_csk_d``,
``hole_N_csk_angle``), including a 45 degree mouth, which is a 90 degree
countersink. A Z extrusion may also carry straight round bosses on its top
face (``boss_N_d``, ``boss_N_h``, ``boss_N_x``, ``boss_N_y``). When every
vertical corner of a Z extrusion, convex and concave, shares one radius, that
radius is ``outline_fillet_r`` and ``outline`` stays the sharp polygon the
walls meet at. The outline vertices and the extrusion length become named
parameters. The box/pad vocabulary stays in charge for parts it already
explains; this module is the fallback.
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
import trimesh

from .fit import Cone, Cylinder, Plane, Sphere, fit_principal_cylinder
from .parametric import Segment, Unsupported, _nice, _plain_or_counterbore, _principal_axis, _snap
from .segment import _partial_span

# Right-handed profile coordinates, in the order stored in ``outline``.
# X extrusion: (y, z). Y extrusion: (x, z). Z extrusion: (x, y).
_PROFILE_AXES = {
    0: (np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
    1: (np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])),
    2: (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])),
}
_AXIS_NAME = "XYZ"


class _Reject(Exception):
    """This extrusion axis does not fit. ``rank`` is how far recognition got."""

    def __init__(self, message: str, rank: int):
        super().__init__(message)
        self.rank = rank


def _unit(v):
    a = np.asarray(v, float)
    n = float(np.linalg.norm(a))
    if n < 1e-12:
        return a
    return a / n


def _attitude(normal, axis: int, ang_deg: float) -> str:
    c = abs(float(_unit(normal)[axis]))
    if c >= math.cos(math.radians(ang_deg)):
        return "cap"
    if c <= math.sin(math.radians(ang_deg)):
        return "wall"
    return "tilted"


def _plane_coord(pl: Plane, axis: int) -> float:
    n = _unit(pl.normal)
    return float(pl.d / n[axis])


def _eps(tol: float, floor: float = 0.15) -> float:
    return max(6.0 * tol, floor)


def _min_area(m) -> float:
    return max(0.05, 1e-4 * float(m.area))


def _loops(mesh, face_ids) -> list[list[int]]:
    """Boundary loops of a set of triangles, as vertex-index rings."""
    counts: dict[tuple[int, int], int] = defaultdict(int)
    faces = mesh.faces
    tris = []
    for fi in face_ids:
        tri = faces[int(fi)]
        abc = (int(tri[0]), int(tri[1]), int(tri[2]))
        tris.append(abc)
        for a, b in ((abc[0], abc[1]), (abc[1], abc[2]), (abc[2], abc[0])):
            counts[(a, b) if a < b else (b, a)] += 1
    nxt: dict[int, list[int]] = defaultdict(list)
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            if counts[(u, v) if u < v else (v, u)] == 1:
                nxt[u].append(v)
    pending = {u: set(vs) for u, vs in nxt.items()}
    loops = []
    for start in list(pending):
        while pending.get(start):
            loop = [start]
            cur = start
            guard = 0
            while pending.get(cur) and guard < 100000:
                guard += 1
                nb = pending[cur].pop()
                if not pending[cur]:
                    del pending[cur]
                if nb == start:
                    break
                loop.append(nb)
                cur = nb
            if len(loop) >= 3 and loop[0] != loop[-1]:
                loops.append(loop)
    return loops


def _project(pts, U, V):
    p = np.asarray(pts, float)
    return np.column_stack([p @ U, p @ V])


def _shoelace(pts) -> float:
    s = 0.0
    for i, (x1, y1) in enumerate(pts):
        x2, y2 = pts[(i + 1) % len(pts)]
        s += x1 * y2 - x2 * y1
    return 0.5 * s


def _ccw(pts):
    pts = [(float(p[0]), float(p[1])) for p in pts]
    if _shoelace(pts) < 0:
        pts = list(reversed(pts))
    k = min(range(len(pts)), key=lambda i: (pts[i][0], pts[i][1]))
    return pts[k:] + pts[:k]


def _drop_near(pts, tol):
    if len(pts) < 2:
        return pts
    out = [pts[0]]
    lim = max(2.0 * tol, 0.05)
    for p in pts[1:]:
        if math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) >= lim:
            out.append(p)
    if len(out) > 1 and math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) < lim:
        out.pop()
    return out


def _drop_collinear(pts, tol):
    pts = _drop_near(list(pts), tol)
    if len(pts) <= 3:
        return pts
    lim = max(4.0 * tol, 0.12)
    changed = True
    while changed and len(pts) > 3:
        changed = False
        n = len(pts)
        for i in range(n):
            a = np.asarray(pts[(i - 1) % n], float)
            b = np.asarray(pts[i], float)
            c = np.asarray(pts[(i + 1) % n], float)
            ac = c - a
            norm = float(np.linalg.norm(ac))
            if norm < 1e-9:
                dist = float(np.linalg.norm(b - a))
            else:
                dist = abs(float(ac[0] * (a[1] - b[1]) - ac[1] * (a[0] - b[0]))) / norm
            if dist <= lim:
                del pts[i]
                changed = True
                break
    return pts


def _wall_line(pl: Plane, U, V):
    n = _unit(pl.normal)
    a, b = float(n @ U), float(n @ V)
    norm = math.hypot(a, b)
    if norm < 1e-8:
        return None
    return (a / norm, b / norm, float(pl.d) / norm)


def _intersect(l1, l2):
    a1, b1, c1 = l1
    a2, b2, c2 = l2
    det = a1 * b2 - a2 * b1
    if abs(det) < 0.05:  # lines already unit-normal; ~3°
        return None
    return ((c1 * b2 - c2 * b1) / det, (a1 * c2 - a2 * c1) / det)


def _merge_walls(planes, tol):
    groups: list[tuple[Plane, list]] = []
    ang = math.cos(math.radians(2.0))
    gap = _eps(tol, 0.08)
    for pl in planes:
        n = _unit(pl.normal)
        d = float(pl.d)
        # distance is n·x; flip so the first component pattern is stable
        if (n[0] < -1e-9) or (abs(n[0]) <= 1e-9 and n[1] < -1e-9) or (
                abs(n[0]) <= 1e-9 and abs(n[1]) <= 1e-9 and n[2] < 0):
            n, d = -n, -d
            pl = Plane(n, d)
        for g, members in groups:
            ng = _unit(g.normal)
            dot = float(n @ ng)
            if dot >= ang and abs(d - float(g.d)) <= gap:
                members.append(pl)
                break
            if dot <= -ang and abs(d + float(g.d)) <= gap:
                members.append(pl)
                break
        else:
            groups.append((pl, [pl]))
    return [g for g, _ in groups]


def _loop_polygon(mesh, loop, walls, U, V, tol, round_corners=False):
    """Corners of one boundary loop from the wall planes it runs along.

    Returns (points, "polygon") or (points2d, "curve") when the loop does not
    lie on walls (a round hole). None when the loop is too broken to use.
    ``round_corners`` keeps the wall intersection when a vertical round puts
    the mesh vertex at the tangent instead of at that sharp corner.
    """
    if len(loop) < 3:
        return None
    pts = mesh.vertices[np.asarray(loop, int)]
    xy = _project(pts, U, V)
    lines = []
    for w in walls:
        ln = _wall_line(w, U, V)
        if ln is not None:
            lines.append(ln)
    if not lines:
        return xy, "curve"
    # label each boundary edge by the wall it lies on
    labels = []
    lim = _eps(tol, 0.2)
    for i in range(len(loop)):
        mid = 0.5 * (xy[i] + xy[(i + 1) % len(xy)])
        best, second = 1e9, 1e9
        best_i = None
        for j, (a, b, c) in enumerate(lines):
            dist = abs(a * mid[0] + b * mid[1] - c)
            if dist < best:
                second = best
                best, best_i = dist, j
            elif dist < second:
                second = dist
        if best_i is None or best > lim or (second < 1e8 and second - best < 0.5 * best and best > 0.05):
            labels.append(None)
        else:
            labels.append(best_i)
    n_wall = sum(lab is not None for lab in labels)
    if n_wall < 0.6 * len(labels):
        return xy, "curve"
    # compress consecutive labels; absorb a one-edge flicker between two runs
    runs = []
    for i, lab in enumerate(labels):
        if lab is None:
            continue
        if runs and runs[-1][0] == lab:
            runs[-1][1].append(i)
        else:
            runs.append([lab, [i]])
    if len(runs) >= 2 and runs[0][0] == runs[-1][0]:
        runs[0][1] = runs[-1][1] + runs[0][1]
        runs.pop()
    # drop tiny runs (a single mislabeled edge at a corner)
    min_len = max(2.0 * tol, 0.3)
    cleaned = []
    for lab, idxs in runs:
        length = 0.0
        for i in idxs:
            length += float(np.linalg.norm(xy[(i + 1) % len(xy)] - xy[i]))
        if length < min_len and cleaned:
            continue
        if cleaned and cleaned[-1][0] == lab:
            cleaned[-1][1].extend(idxs)
        else:
            cleaned.append([lab, idxs])
    if len(cleaned) >= 2 and cleaned[0][0] == cleaned[-1][0]:
        cleaned[0][1] = cleaned[-1][1] + cleaned[0][1]
        cleaned.pop()
    if len(cleaned) < 3:
        return xy, "curve"
    corners = []
    for i, (lab, idxs) in enumerate(cleaned):
        lab2 = cleaned[(i + 1) % len(cleaned)][0]
        hit = _intersect(lines[lab], lines[lab2])
        if hit is None:
            return None
        # vertex where the next wall's edges begin; the corner sits on it
        nxt = cleaned[(i + 1) % len(cleaned)][1][0]
        vtx = xy[nxt]
        if math.hypot(hit[0] - float(vtx[0]), hit[1] - float(vtx[1])) > max(12.0 * tol, 0.8):
            if not round_corners:
                return None
        corners.append(hit)
    corners = _drop_collinear(corners, tol)
    if len(corners) < 3:
        return None
    return corners, "polygon"


def _pip(pt, poly) -> bool:
    x, y = pt
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            xcross = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < xcross:
                inside = not inside
    return inside


def _inside_profile(pt, outline, cutouts) -> bool:
    if not _pip(pt, outline):
        return False
    return not any(_pip(pt, c) for c in cutouts)


def _scan_intervals(poly, fixed_index: int, value: float):
    """Intervals of the free coordinate where the scanline is inside ``poly``.

    ``fixed_index`` 0 holds u fixed (intervals are v); 1 holds v fixed
    (intervals are u).
    """
    vary = 1 - fixed_index
    hits = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        a, b = p[fixed_index], q[fixed_index]
        if abs(a - b) < 1e-12:
            continue
        if (a <= value < b) or (b <= value < a):
            t = (value - a) / (b - a)
            hits.append(p[vary] + t * (q[vary] - p[vary]))
    hits.sort()
    return [(hits[i], hits[i + 1]) for i in range(0, len(hits) - 1, 2)]


def _material_intervals(outline, cutouts, fixed_index: int, value: float):
    """Free-coordinate intervals of material at a constant profile coordinate."""
    outer = _scan_intervals(outline, fixed_index, value)
    if not cutouts:
        return outer
    blocked = []
    for c in cutouts:
        blocked.extend(_scan_intervals(c, fixed_index, value))
    if not blocked:
        return outer
    # subtract blocked from outer
    events = []
    for a, b in outer:
        events.append((a, 1))
        events.append((b, -1))
    for a, b in blocked:
        events.append((a, -1))
        events.append((b, 1))
    events.sort()
    depth = 0
    start = None
    out = []
    for pos, delta in events:
        prev = depth
        depth += delta
        if prev <= 0 and depth > 0:
            start = pos
        elif prev > 0 and depth <= 0 and start is not None:
            if pos - start > 1e-6:
                out.append((start, pos))
            start = None
    return out


def _align_error(a, b):
    if len(a) != len(b) or not a:
        return None
    best = None
    for k in range(len(b)):
        rot = b[k:] + b[:k]
        err = max(math.hypot(p[0] - q[0], p[1] - q[1]) for p, q in zip(a, rot))
        if best is None or err < best[0]:
            best = (err, [(0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])) for p, q in zip(a, rot)])
    return best


def _nice_loop(pts, tol):
    return _ccw([( _nice(p[0], tol), _nice(p[1], tol)) for p in pts])


def _polygons_from_cap(mesh, face_ids, walls, U, V, tol, round_corners=False):
    loops = _loops(mesh, face_ids)
    if not loops:
        return None
    built = []
    for loop in loops:
        got = _loop_polygon(mesh, loop, walls, U, V, tol, round_corners)
        if got is None:
            return None
        pts, kind = got
        area = abs(_shoelace(pts if kind == "polygon" else _ccw(
            [(float(p[0]), float(p[1])) for p in pts[::max(1, len(pts) // 24)]])))
        # area of a curve loop: shoelace of the raw projected ring
        if kind == "curve":
            raw = [(float(p[0]), float(p[1])) for p in pts]
            area = abs(_shoelace(raw))
            built.append((area, kind, raw))
        else:
            built.append((abs(_shoelace(pts)), kind, pts))
    built.sort(key=lambda t: -t[0])
    return built


def _same_loops_ordered(low, high, tol):
    """Average matching polygonal loops. Curve loops are checked loosely."""
    if low is None or high is None or len(low) != len(high):
        return None
    lim = _eps(tol, 0.35)
    averaged = []
    for (a1, k1, p1), (a2, k2, p2) in zip(low, high):
        if k1 != k2:
            return None
        if k1 == "curve":
            # both rings should sit on the same circle; compare centroids and area
            c1 = np.mean(np.asarray(p1, float), axis=0)
            c2 = np.mean(np.asarray(p2, float), axis=0)
            if float(np.linalg.norm(c1 - c2)) > lim or abs(a1 - a2) > max(lim * max(a1, 1.0) ** 0.5, 1.0):
                # area scales with radius squared; allow a loose area gap and a tight center
                if float(np.linalg.norm(c1 - c2)) > lim:
                    return None
            averaged.append((k1, p1))
            continue
        p1, p2 = _ccw(p1), _ccw(p2)
        hit = _align_error(p1, p2)
        if hit is None or hit[0] > lim:
            return None
        averaged.append((k1, hit[1]))
    return averaged


def _same_loops(low, high, tol):
    """Match the two end caps.

    Area order is enough when both ends carry the same circles. A countersink
    or a counterbore makes the mouth wider than the through-hole, so the two
    ends can list those circles in a different order. Pair by centre then,
    and leave a match that already succeeded on the area order alone.
    """
    matched = _same_loops_ordered(low, high, tol)
    if matched is not None or not low or not high or len(low) != len(high):
        return matched
    lim = _eps(tol, 0.35)
    used = [False] * len(high)
    ordered = []
    for _area, kind, pts in low:
        best_i = None
        best_err = None
        for i, (_area2, kind2, pts2) in enumerate(high):
            if used[i] or kind2 != kind:
                continue
            if kind == "curve":
                c1 = np.mean(np.asarray(pts, float), axis=0)
                c2 = np.mean(np.asarray(pts2, float), axis=0)
                err = float(np.linalg.norm(c1 - c2))
            else:
                hit = _align_error(_ccw(pts), _ccw(pts2))
                if hit is None:
                    continue
                err = hit[0]
            if best_err is None or err < best_err:
                best_err = err
                best_i = i
        if best_i is None or best_err > lim:
            return None
        used[best_i] = True
        ordered.append(high[best_i])
    return _same_loops_ordered(low, ordered, tol)


def _on_cylinder(pts, cyl, axis: int, tol) -> bool:
    if _principal_axis(cyl.axis) != axis:
        return False
    rel = np.asarray(pts, float) - np.asarray(cyl.point, float)
    rel[:, axis] = 0.0
    rho = np.linalg.norm(rel, axis=1)
    return float(np.max(np.abs(rho - float(cyl.radius)))) <= _eps(tol, 0.3)


def _outward(mesh, region, cyl) -> bool:
    ax = _unit(cyl.axis)
    c = np.asarray(cyl.point, float)
    fc = mesh.triangles_center[region.faces]
    radial = fc - c
    # drop the axial component
    radial = radial - np.outer(radial @ ax, ax)
    nrm = np.linalg.norm(radial, axis=1, keepdims=True)
    radial = radial / np.maximum(nrm, 1e-12)
    return float(np.mean(np.sum(mesh.face_normals[region.faces] * radial, axis=1))) > 0


def _cone_radius_at(cone, z: float, axis: int) -> float:
    """Cone radius at world coordinate ``z`` along ``axis``.

    ``cone.axis`` is a unit vector from the apex into the cone, so the radius
    is ``h * tan(half_angle)`` with ``h`` measured from the apex.
    """
    apex = np.asarray(cone.apex, float)
    direction = np.asarray(cone.axis, float)
    comp = float(direction[axis])
    if abs(comp) < 1e-9:
        return float("nan")
    h = (float(z) - float(apex[axis])) / comp
    return h * math.tan(float(cone.half_angle))


def _on_cone(pts, cone, axis: int, tol) -> bool:
    """True when ``pts`` sit on ``cone`` and the cone runs along ``axis``."""
    if _principal_axis(cone.axis) != axis:
        return False
    pts = np.asarray(pts, float)
    apex = np.asarray(cone.apex, float)
    direction = np.asarray(cone.axis, float)
    comp = float(direction[axis])
    if abs(comp) < 1e-9:
        return False
    h = (pts[:, axis] - float(apex[axis])) / comp
    if float(np.min(h)) < -_eps(tol, 0.3):
        return False
    expected = h * math.tan(float(cone.half_angle))
    axial = apex + np.outer(h, direction)
    rho = np.linalg.norm(pts - axial, axis=1)
    return float(np.max(np.abs(rho - expected))) <= _eps(tol, 0.3)


def _cone_outward(mesh, region, cone) -> bool:
    ax = _unit(cone.axis)
    c = np.asarray(cone.apex, float)
    fc = mesh.triangles_center[region.faces]
    radial = fc - c
    radial = radial - np.outer(radial @ ax, ax)
    nrm = np.linalg.norm(radial, axis=1, keepdims=True)
    radial = radial / np.maximum(nrm, 1e-12)
    return float(np.mean(np.sum(mesh.face_normals[region.faces] * radial, axis=1))) > 0


def _full_turn(mesh, surface, faces) -> bool:
    """False for a partial cylinder (a fillet). None span is kept, as elsewhere."""
    span = _partial_span(mesh, surface, faces)
    return span is None or span >= np.radians(250.0)


def _turn(prev, corner, nxt):
    """Signed turn at ``corner`` of a CCW outline, and the direction to its fillet centre.

    Positive turn is a convex left turn. The unit vector points from the sharp
    corner toward the centre: into the material for a convex corner, into the
    notch for a concave one.
    """
    din = (corner[0] - prev[0], corner[1] - prev[1])
    dout = (nxt[0] - corner[0], nxt[1] - corner[1])
    lin = math.hypot(din[0], din[1])
    lout = math.hypot(dout[0], dout[1])
    if lin < 1e-9 or lout < 1e-9:
        return None
    din = (din[0] / lin, din[1] / lin)
    dout = (dout[0] / lout, dout[1] / lout)
    cross = din[0] * dout[1] - din[1] * dout[0]
    dot = din[0] * dout[0] + din[1] * dout[1]
    theta = math.atan2(cross, max(-1.0, min(1.0, dot)))
    bis = (-din[1] - dout[1], din[0] + dout[0])
    bl = math.hypot(bis[0], bis[1])
    if bl < 1e-8:
        return None
    sign = 1.0 if theta > 0.0 else -1.0
    return theta, sign * bis[0] / bl, sign * bis[1] / bl


def _fillet_center_at(prev, corner, nxt, radius):
    """``((x, y), turn)`` of the round at one corner, or None."""
    got = _turn(prev, corner, nxt)
    if got is None:
        return None
    theta, ux, uy = got
    half = abs(theta) / 2.0
    # A nearly straight corner, or a hairpin, is not this fillet.
    if half < math.radians(8.0) or abs(abs(theta) - math.pi) < math.radians(4.0):
        return None
    sine = math.sin(half)
    if sine < 1e-6:
        return None
    dist = float(radius) / sine
    return (corner[0] + dist * ux, corner[1] + dist * uy), theta


def _outline_pts(outline):
    pts = [(float(p[0]), float(p[1])) for p in outline]
    if len(pts) < 3:
        return None
    return pts


def outline_fillet_centers(outline, radius: float):
    """Centre of the vertical round at each corner, or None if one corner cannot take it."""
    pts = _outline_pts(outline)
    if pts is None or float(radius) <= 0.0:
        return None
    centers = []
    n = len(pts)
    for i, corner in enumerate(pts):
        got = _fillet_center_at(pts[i - 1], corner, pts[(i + 1) % n], radius)
        if got is None:
            return None
        centers.append(got[0])
    return centers


def outline_fillet_limit(outline) -> float | None:
    """Largest radius that still leaves a straight bit on every edge."""
    pts = _outline_pts(outline)
    if pts is None:
        return None
    tans = []
    n = len(pts)
    for i, corner in enumerate(pts):
        got = _turn(pts[i - 1], corner, pts[(i + 1) % n])
        if got is None:
            return None
        half = abs(got[0]) / 2.0
        if half < math.radians(8.0) or abs(abs(got[0]) - math.pi) < math.radians(4.0):
            return None
        tans.append(math.tan(half))
    limit = None
    for i in range(n):
        nxt = pts[(i + 1) % n]
        length = math.hypot(nxt[0] - pts[i][0], nxt[1] - pts[i][1])
        denom = tans[i] + tans[(i + 1) % n]
        if denom < 1e-9:
            return None
        cap = length / denom
        if limit is None or cap < limit:
            limit = cap
    return limit


def outline_fillet_conflict(outline, radius: float):
    """None when ``radius`` fits.

    ``"corner"`` when a corner cannot take a fillet. A one-tuple is the length
    of an edge whose two rounds meet.
    """
    pts = _outline_pts(outline)
    if pts is None:
        return "corner"
    tans = []
    n = len(pts)
    for i, corner in enumerate(pts):
        got = _turn(pts[i - 1], corner, pts[(i + 1) % n])
        if got is None:
            return "corner"
        half = abs(got[0]) / 2.0
        if half < math.radians(8.0) or abs(abs(got[0]) - math.pi) < math.radians(4.0):
            return "corner"
        tans.append(math.tan(half))
    for i in range(n):
        nxt = pts[(i + 1) % n]
        length = math.hypot(nxt[0] - pts[i][0], nxt[1] - pts[i][1])
        if float(radius) * (tans[i] + tans[(i + 1) % n]) >= length - 1e-9:
            return (length,)
    return None


def _piled_end(verts, sign: int, tol: float) -> float:
    """Median Z of the vertices piled at a free end.

    The single furthest vertex sits past the cap once the mesh is noisy, and
    the cap triangles are often absorbed into the cylinder region.
    """
    along = np.asarray(verts, float)[:, 2]
    extreme = float(along.min() if sign < 0 else along.max())
    band = max(6.0 * tol, 0.4)
    near = along[along <= extreme + band] if sign < 0 else along[along >= extreme - band]
    if near.size < 8:
        return extreme
    return float(np.median(near))


def _absorbed_cap(verts, z_end: float, point, radius: float, tol: float) -> bool:
    """True when vertices pile up on a flat disk at ``z_end`` inside the cylinder."""
    band = max(6.0 * tol, 0.45)
    near = np.asarray(verts, float)
    near = near[np.abs(near[:, 2] - float(z_end)) <= band]
    if len(near) < 12:
        return False
    if float(np.std(near[:, 2])) > max(4.0 * tol, 0.2):
        return False
    rel = near[:, :2] - np.asarray(point, float)[:2]
    rho = np.linalg.norm(rel, axis=1)
    return float(np.percentile(rho, 90)) <= float(radius) + _eps(tol, 0.5)


def _disk_cap(mesh, region, point, radius, tol) -> bool:
    """A Z cap whose vertices stay inside one cylinder (a boss end, or an annulus)."""
    if not isinstance(region.surface, Plane):
        return False
    if _attitude(region.surface.normal, 2, 8.0) != "cap":
        return False
    verts = mesh.vertices[np.unique(mesh.faces[region.faces])]
    if len(verts) < 3:
        return False
    rel = verts - np.asarray(point, float)
    rel[:, 2] = 0.0
    rho = np.linalg.norm(rel, axis=1)
    return float(np.percentile(rho, 95)) <= float(radius) + _eps(tol, 0.45)


def _circle_inside(center, radius, outline, cutouts, slack) -> bool:
    """True when the circle lies in the outline and clear of every cutout.

    ``slack`` shrinks the tested radius, so a fit that drifted by about that
    much still counts as inside.
    """
    rad = float(radius) - float(slack)
    cx, cy = float(center[0]), float(center[1])
    if not _inside_profile((cx, cy), outline, cutouts):
        return False
    if rad <= 0.05:
        return True
    for k in range(48):
        ang = 2.0 * math.pi * k / 48.0
        if not _inside_profile((cx + rad * math.cos(ang), cy + rad * math.sin(ang)),
                               outline, cutouts):
            return False
    return True


def _outward_groups(mesh, cylinders, tol, axis):
    """Coaxial outward cylinders of one radius, along ``axis``."""
    groups = []
    gtol = _eps(tol, 0.35)
    for region in cylinders:
        surface = region.surface
        if not isinstance(surface, Cylinder) or region.area < _min_area(mesh):
            continue
        if _principal_axis(surface.axis) != axis or not _outward(mesh, region, surface):
            continue
        if not _full_turn(mesh, surface, region.faces):
            continue
        point = np.asarray(surface.point, float)
        for group in groups:
            if (_perp_dist(axis, group["point"], point) <= gtol
                    and abs(group["radius"] - float(surface.radius)) <= gtol):
                group["regs"].append(region)
                group["point"] = group["point"] + point
                group["n"] += 1
                break
        else:
            groups.append({
                "point": point.copy(), "radius": float(surface.radius),
                "regs": [region], "n": 1,
            })
    for group in groups:
        group["point"] = group["point"] / group["n"]
    return groups


def _z_boss_plan(mesh, caps, cylinders, tol):
    """Straight top bosses on a Z extrusion, and the end caps that belong to them.

    Returns ``(previews, drop_caps)``. A side boss, a bottom boss, or an
    outward cylinder that is not a straight boss on one top face raises
    :class:`_Reject`. No outward cylinder returns empty lists. A partial
    cylinder is left for the rounded-edge check.
    """
    for region in cylinders:
        surface = region.surface
        if not isinstance(surface, Cylinder) or region.area < _min_area(mesh):
            continue
        ax = _principal_axis(surface.axis)
        if ax is None or ax == 2 or not _outward(mesh, region, surface):
            continue
        if not _full_turn(mesh, surface, region.faces):
            continue
        raise _Reject("a boss on a side face of an extruded outline is not supported yet", 6)
    groups = _outward_groups(mesh, cylinders, tol, 2)
    if not groups:
        return [], []
    gap = _eps(tol, 0.35)
    previews = []
    drop = []
    for group in groups:
        verts = np.vstack([
            mesh.vertices[np.unique(mesh.faces[region.faces])] for region in group["regs"]
        ])
        z0 = float(verts[:, 2].min())
        z1 = float(verts[:, 2].max())
        z_lo = _piled_end(verts, -1, tol)
        z_hi = _piled_end(verts, 1, tol)
        if z_hi - z_lo <= gap:
            raise _Reject("a boss on an extruded outline is not supported yet", 6)
        point = group["point"]
        radius = group["radius"]
        near_lo, near_hi = [], []
        for cap in caps:
            coord = _plane_coord(cap.surface, 2)
            if abs(coord - z_lo) <= gap or abs(coord - z0) <= gap:
                near_lo.append(cap)
            if abs(coord - z_hi) <= gap or abs(coord - z1) <= gap:
                near_hi.append(cap)
        hi_disk = [cap for cap in near_hi if _disk_cap(mesh, cap, point, radius, tol)]
        lo_disk = [cap for cap in near_lo if _disk_cap(mesh, cap, point, radius, tol)]
        lo_ids = {id(cap) for cap in lo_disk}
        hi_ids = {id(cap) for cap in hi_disk}
        lo_broad = [cap for cap in near_lo if id(cap) not in lo_ids]
        hi_broad = [cap for cap in near_hi if id(cap) not in hi_ids]
        flat_top = bool(hi_disk) or _absorbed_cap(verts, z_hi, point, radius, tol)
        # A broad face under the cylinder is the plate top. An extra disk there
        # is the underside of a boss that hangs off that face, not a second part.
        # The end disk is often merged into the cylinder, so the flat top can
        # also be just the vertices piled at z_hi.
        if lo_broad and flat_top:
            previews.append({
                "point": np.asarray(point, float),
                "radius": float(radius),
                "z0": z_lo,
                "z1": z1,
            })
            drop.extend(hi_disk)
            continue
        if lo_disk and hi_broad and not hi_disk:
            raise _Reject(
                "a boss on the bottom face of an extruded outline is not supported yet", 6)
        raise _Reject("a boss on an extruded outline is not supported yet", 6)
    base = previews[0]["z0"]
    if any(abs(item["z0"] - base) > gap for item in previews):
        raise _Reject("a boss does not start on the top face of the extruded outline", 6)
    return previews, drop


def _within_bosses(points, previews, gap) -> bool:
    """True when every point sits inside some preview boss."""
    if len(points) == 0:
        return True
    for point in np.asarray(points, float):
        covered = False
        for boss in previews:
            if point[2] < boss["z0"] - gap or point[2] > boss["z1"] + gap:
                continue
            if math.hypot(point[0] - boss["point"][0], point[1] - boss["point"][1]) <= boss["radius"] + gap:
                covered = True
                break
        if not covered:
            return False
    return True


def _bosses_outside(previews, outline, cutouts, tol) -> bool:
    slack = max(4.0 * tol, 0.2)
    for boss in previews:
        center = (float(boss["point"][0]), float(boss["point"][1]))
        if not _circle_inside(center, boss["radius"], outline, cutouts, slack):
            return True
    return False


def _loops_of(built):
    """Polygon entries and the raw point lists of curve entries."""
    if built is None:
        return None, None
    polygons = [item for item in built if item[1] == "polygon"]
    curves = [item[2] for item in built if item[1] == "curve"]
    return polygons, curves


def _merge_equal_radii(segs, gap):
    """Join coaxial cylinders of one radius that already meet."""
    merged = []
    for seg in segs:
        if (merged and merged[-1].kind == "cylinder" and seg.kind == "cylinder"
                and abs(merged[-1].zb - seg.za) <= gap
                and abs(merged[-1].ra - seg.ra) <= gap
                and abs(merged[-1].rb - seg.rb) <= gap):
            merged[-1].zb = seg.zb
            continue
        merged.append(seg)
    return merged


def _boss_host(xy, radius, bosses, gap):
    """``(index, kind)`` for a hole circle against the boss disks.

    ``kind`` is ``inside`` (fully in exactly one boss), ``overlap`` (meets a
    boss but is not inside it), ``many``, or ``clear``.
    """
    x, y = float(xy[0]), float(xy[1])
    inside = []
    overlap = False
    for index, boss in enumerate(bosses):
        dist = math.hypot(x - boss["x"], y - boss["y"])
        limit = float(boss["d"]) / 2.0
        if dist + float(radius) <= limit + gap:
            inside.append(index)
        elif dist < limit + float(radius) - gap:
            overlap = True
    if len(inside) > 1:
        return None, "many"
    if len(inside) == 1 and not overlap:
        return inside[0], "inside"
    if overlap or inside:
        return None, "overlap"
    return None, "clear"


def _perp_dist(axis: int, p, q) -> float:
    d = np.asarray(p, float) - np.asarray(q, float)
    d[axis] = 0.0
    return float(np.linalg.norm(d))


def _hole_segments(m, regs, ax, tol):
    """Axial pieces of one coaxial hole. The caller snaps the ends onto planes."""
    segs = []
    for r in regs:
        s = r.surface
        verts = m.vertices[np.unique(m.faces[r.faces])]
        along = verts[:, ax]
        segs.append(Segment("cylinder", float(along.min()), float(along.max()),
                            float(s.radius), float(s.radius)))
    segs.sort(key=lambda s: s.za)
    gap = _eps(tol)
    for a, b in zip(segs, segs[1:]):
        if abs(a.zb - b.za) <= gap:
            mid = 0.5 * (a.zb + b.za)
            a.zb = b.za = mid
    return segs


def _top_bosses(mesh, outward_regs, tol, axis, origin, length, outline, cutouts, plane_coords):
    """Straight cylindrical bosses on the top face of a Z extrusion.

    Other outward cylinders are refused here, with a message that says which.
    """
    if not outward_regs:
        return []
    if axis != 2:
        for region in outward_regs:
            if _principal_axis(region.surface.axis) != axis:
                raise Unsupported(
                    "a boss on a side face of an extruded outline is not supported yet")
        raise Unsupported("a boss on an extruded outline is not supported yet")
    for region in outward_regs:
        if _principal_axis(region.surface.axis) != 2:
            raise Unsupported("a boss on a side face of an extruded outline is not supported yet")
    groups = _outward_groups(mesh, outward_regs, tol, 2)
    gap = _eps(tol)
    top = float(origin + length)
    slack = max(4.0 * tol, 0.2)
    bosses = []
    for group in groups:
        verts = np.vstack([
            mesh.vertices[np.unique(mesh.faces[region.faces])] for region in group["regs"]
        ])
        # A noisy cap plane can sit further from the true end than the vertices
        # piled on it. Snap only when the plane is already that close.
        snap = max(2.0 * tol, 0.08)
        z0 = _snap(_piled_end(verts, -1, tol), plane_coords[2], snap)
        z1 = _snap(_piled_end(verts, 1, tol), plane_coords[2], snap)
        if abs(z0 - top) > gap:
            if z1 <= origin + gap and z0 < origin - gap:
                raise Unsupported(
                    "a boss on the bottom face of an extruded outline is not supported yet")
            raise Unsupported("a boss does not start on the top face of the extruded outline")
        height = z1 - top
        if height <= gap:
            raise Unsupported("a boss on an extruded outline is not supported yet")
        center = (float(group["point"][0]), float(group["point"][1]))
        if not _circle_inside(center, group["radius"], outline, cutouts, slack):
            raise Unsupported("a boss is partly outside the extruded outline")
        bosses.append({
            "d": 2.0 * float(group["radius"]),
            "h": float(height),
            "axis": "Z",
            "from": 1,
            "x": center[0],
            "y": center[1],
        })
    bosses.sort(key=lambda boss: (boss["x"], boss["y"], boss["d"]))
    return bosses


def _hole_through_boss(segs, point, bosses, origin, length, gap):
    """1-based boss index when ``segs`` is a hole that comes out that boss top.

    None when the hole does not run past the extrusion. A hole that meets a
    boss and does not come out the top is refused.
    """
    if not bosses or len(segs) != 1 or getattr(segs[0], "kind", "cylinder") != "cylinder":
        return None
    seg = segs[0]
    top = float(origin + length)
    if not (abs(float(seg.za) - float(origin)) <= gap and float(seg.zb) > top + gap):
        return None
    host, kind = _boss_host((float(point[0]), float(point[1])), float(seg.ra), bosses, gap)
    if kind == "clear":
        return None
    if kind != "inside":
        raise Unsupported("a hole meets a boss but does not pass through it")
    if abs(float(seg.zb) - (top + float(bosses[host]["h"]))) <= gap:
        return host + 1
    raise Unsupported("a hole that stops inside a boss is not supported yet")


def _apply_cbore(feat, stepped):
    if "cbore_d" in stepped:
        feat["cbore_d"] = stepped["cbore_d"]
        feat["cbore_depth"] = stepped["cbore_depth"]
        feat["from"] = stepped["from"]
    return feat


def _apply_step(feat, stepped):
    """Copy a counterbore or a countersink onto ``feat``."""
    feat = _apply_cbore(feat, stepped)
    if "csk_d" in stepped:
        feat["csk_d"] = stepped["csk_d"]
        feat["csk_angle"] = stepped["csk_angle"]
        feat["from"] = stepped["from"]
    return feat


def _csk_angle(seg) -> float | None:
    cone = getattr(seg, "cone", None)
    if cone is not None:
        return 2.0 * math.degrees(float(cone.half_angle))
    height = float(seg.zb - seg.za)
    if height <= 1e-6:
        return None
    return 2.0 * math.degrees(math.atan2(abs(float(seg.rb) - float(seg.ra)), height))


def _countersunk_through(segs, lo: float, hi: float, gap: float) -> dict | None:
    """A through-hole with a conical mouth on one end, spanning ``lo``..``hi``.

    The cone opens toward the face (radius grows that way) and meets the
    cylinder. The included angle is ``hole_N_csk_angle``, the same name a box
    uses. A 45 degree mouth is the 90 degree case of that angle.
    """
    if len(segs) != 2:
        return None
    if abs(float(segs[0].za) - lo) > gap or abs(float(segs[-1].zb) - hi) > gap:
        return None
    if abs(float(segs[0].zb) - float(segs[1].za)) > gap:
        return None
    first, second = segs
    if first.kind == "cylinder" and second.kind == "cone":
        # mouth on the high face
        if float(second.rb) <= float(second.ra) + gap:
            return None
        if abs(float(second.ra) - float(first.ra)) > gap:
            return None
        angle = _csk_angle(second)
        if angle is None:
            return None
        return {"d": 2.0 * float(first.ra), "csk_d": 2.0 * float(second.rb),
                "csk_angle": angle, "from": 1}
    if first.kind == "cone" and second.kind == "cylinder":
        if float(first.ra) <= float(first.rb) + gap:
            return None
        if abs(float(first.rb) - float(second.ra)) > gap:
            return None
        angle = _csk_angle(first)
        if angle is None:
            return None
        return {"d": 2.0 * float(second.ra), "csk_d": 2.0 * float(first.ra),
                "csk_angle": angle, "from": -1}
    return None


def _cone_segment(m, region, ax: int):
    """One cone as a segment. Radii are refreshed after the ends snap."""
    cone = region.surface
    verts = m.vertices[np.unique(m.faces[region.faces])]
    along = verts[:, ax]
    za, zb = float(along.min()), float(along.max())
    seg = Segment("cone", za, zb, _cone_radius_at(cone, za, ax), _cone_radius_at(cone, zb, ax))
    seg.cone = cone
    return seg


def _refresh_cone_radii(segs, ax: int) -> None:
    for seg in segs:
        cone = getattr(seg, "cone", None)
        if cone is None:
            continue
        seg.ra = _cone_radius_at(cone, seg.za, ax)
        seg.rb = _cone_radius_at(cone, seg.zb, ax)


def _cylinder_features(m, regions, tol, axis, origin, length, outline, cutouts, cones=(),
                       skip=()):
    """Through-holes, counterbored through-holes, and straight top bosses on Z.

    A Z extrusion may also carry a countersunk through-hole along Z. A hole
    that continues through one of the top bosses is the same through-hole,
    with the cutter extended by that boss height. A hole that meets a boss and
    stops short is refused. ``skip`` is the vertical corner rounds already
    named as ``outline_fillet_r``. Any other partial cylinder stays refused.
    """
    U, V = _PROFILE_AXES[axis]
    skip_ids = set(skip)
    cyl_regs = []
    outward_regs = []
    for r in regions:
        if id(r) in skip_ids:
            continue
        s = r.surface
        if not isinstance(s, Cylinder) or not len(r.faces):
            continue
        if r.area < _min_area(m):
            continue
        ax = _principal_axis(s.axis)
        if ax is None:
            raise Unsupported("curved surface whose axis is not aligned with X, Y or Z")
        span = _partial_span(m, s, r.faces)
        if span is not None and span < np.radians(250.0):
            raise Unsupported("a rounded edge is not part of an editable outline")
        if _outward(m, r, s):
            outward_regs.append(r)
            continue
        cyl_regs.append(r)
    # snap targets: planes perpendicular to each axis
    plane_coords = {0: [], 1: [], 2: []}
    for r in regions:
        s = r.surface
        if isinstance(s, Plane) and r.area >= _min_area(m):
            for ax in range(3):
                if _attitude(s.normal, ax, 4.0) == "cap":
                    plane_coords[ax].append(_plane_coord(s, ax))
    bosses = _top_bosses(m, outward_regs, tol, axis, origin, length, outline, cutouts, plane_coords)
    groups = []
    gtol = _eps(tol, 0.35)
    for r in cyl_regs:
        s = r.surface
        ax = _principal_axis(s.axis)
        pt = np.asarray(s.point, float)
        for g in groups:
            if g["axis"] == ax and _perp_dist(ax, g["point"], pt) <= gtol:
                g["regs"].append(r)
                g["point"] = g["point"] + pt
                g["n"] += 1
                break
        else:
            groups.append({"axis": ax, "point": pt.copy(), "regs": [r], "n": 1, "cones": []})
    for region in cones:
        cone = region.surface
        cax = _principal_axis(cone.axis)
        tip = np.asarray(cone.apex, float)
        host = None
        for g in groups:
            if g["axis"] == cax == 2 and _perp_dist(cax, g["point"] / g["n"], tip) <= gtol:
                host = g
                break
        if host is None:
            raise Unsupported(
                "a countersink, chamfer or draft on an extruded outline is not supported yet")
        host["cones"].append(region)
    holes = []
    gap = _eps(tol)
    for g in groups:
        ax = g["axis"]
        pt = g["point"] / g["n"]
        segs = _hole_segments(m, g["regs"], ax, tol)
        cone_regs = g.get("cones") or ()
        if cone_regs:
            # A countersink is a vertical cone on a Z extrusion. Anything else
            # (a draft, a chamfer on an X or Y hole) keeps the old refusal.
            if ax != 2 or axis != 2:
                raise Unsupported(
                    "a countersink, chamfer or draft on an extruded outline is not supported yet")
            for region in cone_regs:
                segs.append(_cone_segment(m, region, ax))
            segs.sort(key=lambda s: (s.za, 0 if s.kind == "cylinder" else 1))
        for s in segs:
            if s.kind == "cone" and ax == axis:
                # Snap the mouth onto an end face. The junction is not a cap,
                # and another hole's counterbore floor must not pull it there.
                ends = (float(origin), float(origin + length))
                s.za = _snap(s.za, ends, gap)
                s.zb = _snap(s.zb, ends, gap)
            else:
                s.za = _snap(s.za, plane_coords[ax], gap)
                s.zb = _snap(s.zb, plane_coords[ax], gap)
        if cone_regs:
            _refresh_cone_radii(segs, ax)
        for a, b in zip(segs, segs[1:]):
            if abs(a.zb - b.za) <= gap:
                mid = 0.5 * (a.zb + b.za)
                a.zb = b.za = mid
        if cone_regs:
            _refresh_cone_radii(segs, ax)
        segs = _merge_equal_radii(segs, gap)
        name = _AXIS_NAME[ax]
        if ax == axis:
            uv = (float(pt @ U), float(pt @ V))
            if not _inside_profile(uv, outline, cutouts):
                raise Unsupported("hole is outside the extruded outline")
            through = _hole_through_boss(segs, pt, bosses, origin, length, gap)
            stepped = _plain_or_counterbore(segs, origin, origin + length, gap)
            if stepped is None and axis == 2:
                stepped = _countersunk_through(segs, origin, origin + length, gap)
            if through is not None:
                stepped = {"d": 2.0 * float(segs[0].ra)}
            elif stepped is not None and bosses and _boss_host(
                    (float(pt[0]), float(pt[1])), float(segs[0].ra), bosses, gap)[1] != "clear":
                raise Unsupported("a hole meets a boss but does not pass through it")
            if stepped is None:
                if len(segs) == 1 and segs[0].kind == "cylinder":
                    raise Unsupported("hole along the extrusion does not pass through both ends")
                if any(s.kind == "cone" for s in segs):
                    raise Unsupported(
                        "a countersink, chamfer or draft on an extruded outline is not supported yet")
                raise Unsupported("stepped hole is not a counterbored through-hole")
            feat = {"d": stepped["d"], "axis": name, "span": "extrude"}
            for i, coord in enumerate("xyz"):
                if i != ax:
                    feat[coord] = float(pt[i])
            if through is not None:
                feat["through_boss"] = through
            holes.append(_apply_step(feat, stepped))
            continue
        profile_axes = [i for i in range(3) if i != axis]
        fixed_world = profile_axes[0] if profile_axes[1] == ax else profile_axes[1]
        fixed_index = profile_axes.index(fixed_world)
        value = float(pt[fixed_world])
        intervals = _material_intervals(outline, cutouts, fixed_index, value)
        if not intervals:
            raise Unsupported("hole is outside the extruded outline")
        lim = _eps(tol, 0.35)
        za = min(s.za for s in segs)
        zb = max(s.zb for s in segs)
        if not all(za - lim <= a and b <= zb + lim for a, b in intervals):
            raise Unsupported("hole does not pass through the extruded profile")
        if not (origin - lim <= float(pt[axis]) <= origin + length + lim):
            raise Unsupported("hole is outside the extruded profile")
        # one contiguous run of material: the counterbore opens on one side of it
        lo = min(a for a, _b in intervals)
        hi = max(b for _a, b in intervals)
        stepped = _plain_or_counterbore(segs, lo, hi, gap)
        if stepped is None:
            if len(segs) == 1:
                raise Unsupported("hole does not pass through the extruded profile")
            raise Unsupported("stepped hole is not a counterbored through-hole")
        feat = {"d": stepped["d"], "axis": name, "span": ("outline", 1 - fixed_index)}
        for i, coord in enumerate("xyz"):
            if i != ax:
                feat[coord] = float(pt[i])
        holes.append(_apply_cbore(feat, stepped))
    holes.sort(key=lambda h: (h["axis"], h.get("x", 0.0), h.get("y", 0.0), h.get("z", 0.0)))
    return holes, bosses


def profile_extent(hole: dict, lang: str) -> tuple[str, str]:
    """``(low, high)`` of the material a hole passes through, as script expressions."""
    span = hole["span"]
    if span == "extrude":
        hi = "(extrude_0 + extrude_len)"
        boss = hole.get("through_boss")
        if boss:
            hi = f"(extrude_0 + extrude_len + boss_{boss}_h)"
        return "extrude_0", hi
    idx = span[1]
    if lang == "cadquery":
        lo = f"min(p[{idx}] for p in outline)"
        hi = f"max(p[{idx}] for p in outline)"
    else:
        lo = f"min([for (p = outline) p[{idx}]])"
        hi = f"max([for (p = outline) p[{idx}]])"
    return lo, hi


def profile_cutter(hole: dict, lang: str) -> tuple[str, str]:
    """``(start, length)`` of a through-hole cutter along its axis, past the ends by 1."""
    lo, hi = profile_extent(hole, lang)
    return f"({lo} - 1)", f"({hi} - ({lo}) + 2)"


def profile_cbore_cutter(hole: dict, name: str, lang: str) -> tuple[str, str]:
    """``(start, length)`` of the counterbore, overlapping the face it opens on by 1."""
    lo, hi = profile_extent(hole, lang)
    if hole.get("from", 1) < 0:
        return f"({lo} - 1)", f"{name}_cbore_depth + 1"
    return f"({hi} - {name}_cbore_depth)", f"{name}_cbore_depth + 1"


def profile_csk_exprs(hole: dict, name: str, lang: str) -> tuple[str, str, str, str, str] | None:
    """Countersink cone as ``(height, start, length, r1, r2)`` expressions.

    The cone runs 1 mm past the face it opens on, and the radius follows the
    countersink slope over that extra millimetre, the same construction the
    box vocabulary uses. ``lang`` is ``cadquery`` or ``openscad``. FreeCAD
    uses the CadQuery expressions. The height is assigned to ``{name}_csk_h``
    before the cutter, so the start of a top countersink can name it.
    """
    if "csk_d" not in hole:
        return None
    lo, hi = profile_extent(hole, lang)
    if lang == "openscad":
        slope = f"tan({name}_csk_angle / 2)"
    else:
        slope = f"math.tan(math.radians({name}_csk_angle / 2))"
    height = f"({name}_csk_d - {name}_d) / 2 / {slope}"
    if hole.get("from", 1) < 0:
        start = f"({lo} - 1)"
        r1 = f"{name}_csk_d / 2 + {slope}"
        r2 = f"{name}_d / 2"
    else:
        start = f"(({hi}) - {name}_csk_h)"
        r1 = f"{name}_d / 2"
        r2 = f"{name}_csk_d / 2 + {slope}"
    length = f"({name}_csk_h + 1)"
    return height, start, length, r1, r2


def profile_parameters(feats: dict, tol: float) -> list[tuple[str, float, str]]:
    """Float parameters of an extruded profile, in emitter order."""
    pr = feats["profile"]
    ax = pr["axis"]
    P = [
        ("extrude_len", _nice(pr["length"], tol), f"length along {ax}"),
        ("extrude_0", _nice(pr["origin"], tol), f"start along {ax}"),
    ]
    if pr.get("fillet_r"):
        P.append(("outline_fillet_r", _nice(pr["fillet_r"], tol),
                  "radius of every vertical corner"))
    for i, d in enumerate(feats.get("bosses") or [], 1):
        n = f"boss_{i}"
        P.append((f"{n}_d", _nice(d["d"], tol), "diameter at the top face"))
        P.append((f"{n}_h", _nice(d["h"], tol), "height"))
        P.append((f"{n}_x", _nice(d["x"], tol), ""))
        P.append((f"{n}_y", _nice(d["y"], tol), ""))
    for i, d in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        P.append((f"{n}_d", _nice(d["d"], tol), f"diameter (axis {d['axis']})"))
        for coord in "xyz":
            if coord != d["axis"].lower() and coord in d:
                P.append((f"{n}_{coord}", _nice(d[coord], tol), ""))
        if "csk_d" in d:
            end = f"{'+' if d.get('from', 1) > 0 else '-'}{d['axis']}"
            P.append((f"{n}_csk_d", _nice(d["csk_d"], tol), f"countersink diameter at the {end} end"))
            # Half a degree, the same rounding the box vocabulary uses for this angle.
            P.append((f"{n}_csk_angle", _nice(d["csk_angle"], 0.5), "countersink angle"))
        if "cbore_d" in d:
            end = f"{'+' if d.get('from', 1) > 0 else '-'}{d['axis']}"
            P.append((f"{n}_cbore_d", _nice(d["cbore_d"], tol), "counterbore diameter"))
            P.append((f"{n}_cbore_depth", _nice(d["cbore_depth"], tol),
                      f"counterbore depth from the {end} end"))
    return P


def profile_point_lists(feats: dict) -> list[tuple[str, list, str]]:
    pr = feats["profile"]
    plane, coords = {"X": ("YZ", "y, z"), "Y": ("XZ", "x, z"), "Z": ("XY", "x, y")}[pr["axis"]]
    rows = [("outline", pr["outline"], f"outline in the {plane} plane, ({coords})")]
    for i, loop in enumerate(pr["cutouts"], 1):
        rows.append((f"cutout_{i}", loop, "through opening in the outline"))
    return rows


def _cap_faces(m, cluster, axis):
    ids = []
    for pl, _coord, region in cluster:
        ids.extend(int(i) for i in region.faces)
    return ids


def _vertical_fillets(m, cylinders, tol, origin, length):
    """Partial Z cylinders that can be the vertical corners of one extrusion.

    A nearly full turn is a hole. A fragment on that same axis is part of the
    hole. Pieces of one round are merged, and the round has to run the height
    of the extrusion.
    """
    full = []
    parts = []
    for region in cylinders:
        surface = region.surface
        if not isinstance(surface, Cylinder) or region.area < _min_area(m):
            continue
        if _principal_axis(surface.axis) != 2:
            continue
        span = _partial_span(m, surface, region.faces)
        point = np.asarray(surface.point, float)
        if span is None or span >= np.radians(250.0):
            full.append((float(point[0]), float(point[1]), float(surface.radius)))
            continue
        if span < np.radians(8.0):
            continue
        verts = m.vertices[np.unique(m.faces[region.faces])]
        parts.append({
            "radius": float(surface.radius),
            "x": float(point[0]),
            "y": float(point[1]),
            "span": float(span),
            "area": float(region.area),
            "z0": float(verts[:, 2].min()),
            "z1": float(verts[:, 2].max()),
            "ids": [id(region)],
        })
    gap = _eps(tol, 0.35)
    kept = []
    for part in parts:
        if any(math.hypot(part["x"] - fx, part["y"] - fy) <= gap and abs(part["radius"] - fr) <= gap
               for fx, fy, fr in full):
            continue
        kept.append(part)
    kept.sort(key=lambda item: -item["area"])
    merged = []
    near = max(4.0 * tol, 0.4)
    for part in kept:
        for group in merged:
            if (math.hypot(group["x"] - part["x"], group["y"] - part["y"]) <= near
                    and abs(group["radius"] - part["radius"]) <= near
                    and group["span"] + part["span"] <= np.radians(210.0)):
                area = group["area"] + part["area"]
                group["x"] = (group["x"] * group["area"] + part["x"] * part["area"]) / area
                group["y"] = (group["y"] * group["area"] + part["y"] * part["area"]) / area
                group["radius"] = (group["radius"] * group["area"] + part["radius"] * part["area"]) / area
                group["span"] += part["span"]
                group["area"] = area
                group["z0"] = min(group["z0"], part["z0"])
                group["z1"] = max(group["z1"], part["z1"])
                group["ids"].extend(part["ids"])
                break
        else:
            merged.append(dict(part))
    zgap = max(6.0 * tol, 0.55)
    out = []
    for group in merged:
        if group["span"] < np.radians(12.0) or group["span"] > np.radians(200.0):
            continue
        if group["z0"] > float(origin) + zgap or group["z1"] < float(origin + length) - zgap:
            continue
        out.append(group)
    return out


def _line_side(line, x, y) -> float:
    a, b, c = line
    return a * x + b * y - c


def _corner_of_fillet(fillet, lines, tol):
    """Sharp corner where the two walls tangent to ``fillet`` meet.

    The walls are shifted so each is exactly the fillet radius from the
    centre. A noisy plane would otherwise miss the corner by a millimetre.
    """
    fx, fy, fr = float(fillet["x"]), float(fillet["y"]), float(fillet["radius"])
    tol_d = max(4.0 * tol, 0.45)
    scored = []
    for line in lines:
        signed = _line_side(line, fx, fy)
        err = abs(abs(signed) - fr)
        if err <= tol_d:
            scored.append((err, line, signed))
    if len(scored) < 2:
        return None
    scored.sort(key=lambda item: item[0])
    pair = None
    for i in range(len(scored)):
        for j in range(i + 1, len(scored)):
            a1, b1, _c1 = scored[i][1]
            a2, b2, _c2 = scored[j][1]
            if abs(a1 * b2 - a2 * b1) >= 0.2:
                pair = (scored[i], scored[j])
                break
        if pair is not None:
            break
    if pair is None:
        return None
    corrected = []
    for _err, line, signed in pair:
        a, b, _c = line
        sign = 1.0 if signed >= 0.0 else -1.0
        corrected.append((a, b, a * fx + b * fy - sign * fr))
    return _intersect(corrected[0], corrected[1])


def _ordered_fillet_corners(mesh, face_ids, fillets, lines, U, V, tol):
    """Sharp corners in boundary order, one per fillet on the outer loop."""
    loops = _loops(mesh, face_ids)
    if not loops or not fillets:
        return None
    best = None
    for loop in loops:
        xy = _project(mesh.vertices[np.asarray(loop, int)], U, V)
        area = abs(_shoelace([(float(p[0]), float(p[1])) for p in xy]))
        if best is None or area > best[0]:
            best = (area, xy)
    if best is None:
        return None
    band = max(4.0 * tol, 0.45)
    seq = []
    for p in best[1]:
        hit, hit_err = None, None
        for i, fil in enumerate(fillets):
            err = abs(math.hypot(float(p[0]) - fil["x"], float(p[1]) - fil["y"]) - fil["radius"])
            if err <= band and (hit_err is None or err < hit_err):
                hit, hit_err = i, err
        if hit is None:
            continue
        if not seq or seq[-1] != hit:
            seq.append(hit)
    if len(seq) >= 2 and seq[0] == seq[-1]:
        seq.pop()
    if len(seq) != len(fillets) or len(set(seq)) != len(fillets):
        return None
    corners = []
    for index in seq:
        hit = _corner_of_fillet(fillets[index], lines, tol)
        if hit is None:
            return None
        corners.append(hit)
    return corners


def _average_outline(low, high, tol):
    """Mean of two outlines, trying the reverse winding of ``high``."""
    if low is None:
        return high
    if high is None:
        return low
    lim = _eps(tol, 0.4)
    best = None
    for candidate in (high, list(reversed(high))):
        hit = _align_error(low, candidate)
        if hit is None:
            continue
        if best is None or hit[0] < best[0]:
            best = hit
    if best is None or best[0] > lim:
        return None
    return best[1]


def _rebuild_fillet_outline(m, low, high, axis, wall_planes, fillets, U, V, tol):
    """Sharp outline from the fillet centres when the wall intersections drift."""
    lines = []
    for wall in wall_planes:
        line = _wall_line(wall, U, V)
        if line is not None:
            lines.append(line)
    if len(lines) < 3:
        return None
    low_pts = _ordered_fillet_corners(m, _cap_faces(m, low, axis), fillets, lines, U, V, tol)
    high_pts = _ordered_fillet_corners(m, _cap_faces(m, high, axis), fillets, lines, U, V, tol)
    averaged = _average_outline(low_pts, high_pts, tol)
    if not averaged or len(averaged) < 3:
        return None
    outline = _nice_loop(_drop_collinear(averaged, tol), tol)
    if len(outline) < 3 or abs(_shoelace(outline)) < 1e-3:
        return None
    return outline


def _assign_outline_fillets(outline, fillets, tol):
    """One shared radius when every corner owns one of ``fillets``.

    Returns ``(radius, region ids)`` or ``(None, empty)``. Radii that do not
    agree stay unnamed: there is no per-corner parameter.
    """
    if not fillets or len(fillets) != len(outline):
        return None, set()
    radii = [float(item["radius"]) for item in fillets]
    med = float(np.median(radii))
    if med <= max(2.0 * tol, 0.05):
        return None, set()
    if max(abs(radius - med) for radius in radii) > max(4.0 * tol, 0.35):
        return None, set()
    pts = [(float(p[0]), float(p[1])) for p in outline]
    unused = list(fillets)
    limit = max(8.0 * tol, 0.6)
    # Noise pulls the tangent triangles into the wall, so the fitted arc is
    # shorter than the corner turn. The centre still has to sit on that corner.
    span_tol = np.radians(55.0)
    for i, corner in enumerate(pts):
        got = _fillet_center_at(pts[i - 1], corner, pts[(i + 1) % len(pts)], med)
        if got is None:
            return None, set()
        center, theta = got
        best_i, best_d = None, None
        for j, cand in enumerate(unused):
            dist = math.hypot(cand["x"] - center[0], cand["y"] - center[1])
            if best_d is None or dist < best_d:
                best_i, best_d = j, dist
        if best_i is None or best_d > limit or abs(unused[best_i]["span"] - abs(theta)) > span_tol:
            return None, set()
        unused.pop(best_i)
    ids = set()
    for item in fillets:
        ids.update(item["ids"])
    return med, ids


def _matched_outline(m, low, high, axis, wall_planes, U, V, tol, round_corners, boss_preview):
    """Sharp outline, cutouts, and curve loops shared by the two end caps."""
    low_poly = _polygons_from_cap(
        m, _cap_faces(m, low, axis), wall_planes, U, V, tol, round_corners)
    high_poly = _polygons_from_cap(
        m, _cap_faces(m, high, axis), wall_planes, U, V, tol, round_corners)
    curve_extra: list = []
    if boss_preview:
        # Boss footprints and a hole that continues through a boss are curves
        # on one end only. The polygonal outline still has to match.
        low_p, low_c = _loops_of(low_poly)
        high_p, high_c = _loops_of(high_poly)
        matched = _same_loops(low_p, high_p, tol) if low_p and high_p else None
        if matched is None and low_p:
            outline_try = _nice_loop(_drop_collinear(low_p[0][2], tol), tol)
            cut_try = []
            for _area, _kind, pts in low_p[1:]:
                loop = _nice_loop(_drop_collinear(pts, tol), tol)
                if len(loop) >= 3 and abs(_shoelace(loop)) > 1e-3:
                    cut_try.append(loop)
            if len(outline_try) >= 3 and _bosses_outside(boss_preview, outline_try, cut_try, tol):
                raise _Reject("a boss is partly outside the extruded outline", 6)
        curve_extra = list(low_c or []) + list(high_c or [])
    else:
        matched = _same_loops(low_poly, high_poly, tol)
    if matched is None:
        # one end can be a single triangle fan that still yields the outline;
        # accept it when the other end agrees, otherwise say so
        raise _Reject(
            "the two ends do not have the same outline, so this is not one straight extrusion",
            5,
        )
    if matched[0][0] != "polygon":
        raise _Reject("the end face is not bounded by straight walls", 3)
    outline = _nice_loop(_drop_collinear(matched[0][1], tol), tol)
    if len(outline) < 3 or abs(_shoelace(outline)) < 1e-3:
        raise _Reject("the outline area is zero", 3)
    cutouts = []
    curve_loops = []
    for kind, pts in matched[1:]:
        if kind == "polygon":
            loop = _nice_loop(_drop_collinear(pts, tol), tol)
            if len(loop) >= 3 and abs(_shoelace(loop)) > 1e-3:
                cutouts.append(loop)
        else:
            curve_loops.append(pts)
    curve_loops.extend(curve_extra)
    return outline, cutouts, curve_loops


def _along(m, regions, tol, axis: int) -> dict:
    ang = 4.0
    planes = []
    cylinders = []
    for r in regions:
        s = r.surface
        if isinstance(s, Plane):
            if r.area < _min_area(m):
                continue
            planes.append(r)
        elif isinstance(s, Cylinder):
            cylinders.append(r)
        elif isinstance(s, Cone):
            continue
        elif isinstance(s, Sphere):
            raise _Reject("sphere faces are not supported yet", 0)
        elif s is None:
            raise _Reject("part has unfitted (faceted) regions", 0)
        else:
            raise _Reject(f"{s.kind} faces are not supported yet", 0)
    tilted = []
    caps = []
    walls = []
    for r in planes:
        kind = _attitude(r.surface.normal, axis, ang)
        if kind == "tilted":
            tilted.append(r)
        elif kind == "cap":
            caps.append(r)
        else:
            walls.append(r)
    if tilted:
        tilted.sort(key=lambda r: -r.area)
        largest = tilted[0].area
        raise _Reject(
            f"{len(tilted)} face{'s' if len(tilted) != 1 else ''} "
            f"{'are' if len(tilted) != 1 else 'is'} tilted "
            f"(largest {largest:.0f} mm²)",
            1,
        )
    boss_preview = []
    if axis == 2:
        # A drafted boss is an outward cone that is not the 45 degree chamfer
        # a combined extrusion already writes. Name it before the extra end
        # cap is read as a second extrusion. A 45 degree cone stays on the
        # path that already handles hole chamfers.
        for region in regions:
            cone = region.surface
            if not isinstance(cone, Cone) or region.area < _min_area(m):
                continue
            if _principal_axis(cone.axis) != 2 or not _cone_outward(m, region, cone):
                continue
            if abs(float(cone.half_angle) - math.pi / 4.0) <= math.radians(3.0):
                continue
            raise _Reject("a tapered boss on an extruded outline is not supported yet", 6)
        boss_preview, drop_caps = _z_boss_plan(m, caps, cylinders, tol)
        if drop_caps:
            drop_ids = {id(cap) for cap in drop_caps}
            caps = [cap for cap in caps if id(cap) not in drop_ids]
    if len(caps) < 2:
        raise _Reject("fewer than two end faces", 2)
    coords = [(_plane_coord(r.surface, axis), r) for r in caps]
    coords.sort(key=lambda t: t[0])
    gap = _eps(tol, 0.2)
    clusters: list[list] = []
    for coord, region in coords:
        if not clusters or coord - clusters[-1][-1][1] > gap:
            clusters.append([(region.surface, coord, region)])
        else:
            clusters[-1].append((region.surface, coord, region))
    if len(clusters) < 2:
        raise _Reject("the end faces are not two parallel planes", 2)
    # intermediate caps must be ends of a round hole, not a step in the outline
    if len(clusters) > 2:
        for cluster in clusters[1:-1]:
            for _pl, _c, region in cluster:
                verts = m.vertices[np.unique(m.faces[region.faces])]
                explained = False
                for cr in cylinders:
                    cyl = cr.surface
                    if _principal_axis(cyl.axis) != axis:
                        continue
                    rel = verts - np.asarray(cyl.point, float)
                    rel[:, axis] = 0.0
                    if float(np.percentile(np.linalg.norm(rel, axis=1), 90)) <= float(cyl.radius) + gap:
                        explained = True
                        break
                if not explained:
                    raise _Reject(
                        "a face parallel to the ends sits in the middle, so this is more than one extrusion",
                        4,
                    )
    low, high = clusters[0], clusters[-1]
    origin = float(np.average([c for _p, c, _r in low], weights=[r.area for _p, _c, r in low]))
    end = float(np.average([c for _p, c, _r in high], weights=[r.area for _p, _c, r in high]))
    length = end - origin
    if length <= gap:
        raise _Reject("the extrusion length is zero", 2)
    extent = m.vertices[:, axis]
    if abs(float(extent.min()) - origin) > gap:
        raise _Reject("the solid sticks out past the two end faces", 3)
    if boss_preview:
        if abs(end - boss_preview[0]["z0"]) > _eps(tol, 0.35):
            raise _Reject("a boss does not start on the top face of the extruded outline", 6)
        above = m.vertices[extent > end + gap]
        if not _within_bosses(above, boss_preview, _eps(tol, 0.35)):
            raise _Reject("the solid sticks out past the two end faces", 3)
    elif abs(float(extent.max()) - end) > gap:
        raise _Reject("the solid sticks out past the two end faces", 3)
    wall_planes = _merge_walls([r.surface for r in walls], tol)
    if len(wall_planes) < 3:
        raise _Reject("not enough side walls to form an outline", 2)
    U, V = _PROFILE_AXES[axis]
    # Vertical rounds are cylinder patches, not arcs stored in the outline.
    # The walls still meet at the sharp corner. One shared radius on every
    # corner of a Z extrusion becomes outline_fillet_r. Anything else keeps
    # the old refusal.
    fillets = _vertical_fillets(m, cylinders, tol, origin, length) if axis == 2 else []
    modes = (True, False) if fillets else (False,)
    chosen = None
    fail = None
    for round_corners in modes:
        try:
            outline, cutouts, curve_loops = _matched_outline(
                m, low, high, axis, wall_planes, U, V, tol, round_corners, boss_preview)
        except _Reject as exc:
            fail = exc
            outline = cutouts = curve_loops = None
        radius, fillet_ids = (None, set())
        if round_corners:
            # Prefer the cap polygon. When its corners do not match the rounds
            # (noise pulls a wall intersection), rebuild each corner from the
            # fillet centre, which stays put.
            if outline is not None:
                radius, fillet_ids = _assign_outline_fillets(outline, fillets, tol)
            if radius is None:
                rebuilt = _rebuild_fillet_outline(
                    m, low, high, axis, wall_planes, fillets, U, V, tol)
                if rebuilt is not None:
                    radius, fillet_ids = _assign_outline_fillets(rebuilt, fillets, tol)
                    if radius is not None:
                        if outline is None:
                            outline, cutouts, curve_loops = rebuilt, [], []
                        else:
                            outline = rebuilt
            if radius is None or outline is None:
                continue
        elif outline is None:
            continue
        chosen = (outline, cutouts, curve_loops, radius, fillet_ids)
        break
    if chosen is None:
        if fillets and (fail is None or fail.rank < 6):
            raise _Reject("a rounded edge is not part of an editable outline", 6)
        if fail is not None:
            raise fail
        raise _Reject("a rounded edge is not part of an editable outline", 6)
    outline, cutouts, curve_loops, fillet_radius, fillet_ids = chosen
    if boss_preview and _bosses_outside(boss_preview, outline, cutouts, tol):
        raise _Reject("a boss is partly outside the extruded outline", 6)
    # a round hole is a non-polygonal loop on a cylinder along the extrusion
    axial = [r for r in cylinders if _principal_axis(r.surface.axis) == axis]
    axis_vec = np.zeros(3)
    axis_vec[axis] = 1.0
    cone_regions = [r for r in regions if isinstance(r.surface, Cone)]
    for pts in curve_loops:
        # A countersink mouth lies on an end face, not on the start plane.
        # A hole wall is a cylinder, so the axial position does not matter.
        on_wall = False
        on_mouth = False
        for z_face in (float(origin), float(origin + length)):
            cloud = np.array([p[0] * U + p[1] * V + z_face * axis_vec for p in pts])
            if any(_on_cylinder(cloud, r.surface, axis, tol) for r in axial):
                on_wall = True
            if any(_on_cone(cloud, r.surface, axis, tol) for r in cone_regions):
                on_mouth = True
        if not on_wall and not on_mouth:
            raise _Reject(
                "an opening in the outline is neither a straight cutout nor a round hole",
                5,
            )
    # every significant wall should be an edge of the outline or a cutout
    tol_edge = _eps(tol, 0.25)
    for w in wall_planes:
        ln = _wall_line(w, U, V)
        if ln is None:
            continue
        a, b, c = ln
        used = False
        for loop in [outline, *cutouts]:
            for i in range(len(loop)):
                p, q = loop[i], loop[(i + 1) % len(loop)]
                if abs(a * p[0] + b * p[1] - c) <= tol_edge and abs(a * q[0] + b * q[1] - c) <= tol_edge:
                    used = True
                    break
            if used:
                break
        if not used and fillet_radius is not None:
            # A rebuilt corner sits a few tenths off the noisy plane, and a
            # chord of a round can be fitted as its own vertical plane.
            loose = max(8.0 * tol, 0.6)
            for loop in [outline, *cutouts]:
                for i in range(len(loop)):
                    p, q = loop[i], loop[(i + 1) % len(loop)]
                    if (abs(a * p[0] + b * p[1] - c) <= loose
                            and abs(a * q[0] + b * q[1] - c) <= loose):
                        used = True
                        break
                if used:
                    break
            if not used:
                for fil in fillets:
                    if abs(a * fil["x"] + b * fil["y"] - c) <= float(fil["radius"]) + loose:
                        used = True
                        break
        if not used:
            raise _Reject("a side face is not part of the extruded outline", 4)
    # walls must run the full length (a step would end early). Coplanar regions share a plane.
    for w in wall_planes:
        n = _unit(w.normal)
        d = float(w.d)
        verts = []
        for r in walls:
            rn = _unit(r.surface.normal)
            rd = float(r.surface.d)
            dot = float(n @ rn)
            if dot > 0.99 and abs(d - rd) <= gap:
                verts.append(m.vertices[np.unique(m.faces[r.faces])])
            elif dot < -0.99 and abs(d + rd) <= gap:
                verts.append(m.vertices[np.unique(m.faces[r.faces])])
        if not verts:
            continue
        vv = np.vstack(verts)
        if float(vv[:, axis].max() - vv[:, axis].min()) < length - gap:
            raise _Reject(
                "a side face does not run the whole length, so the outline is not constant",
                4,
            )
    # A vertical countersink on a Z extrusion is a cone. A draft, a boss
    # chamfer, and a cone whose axis is not Z stay refused, with the same
    # message as before. An X or Y extrusion never takes the countersink path,
    # so a 45 degree chamfer on a side hole still falls through to the
    # combined-extrusion script.
    hole_cones = []
    if cone_regions:
        if axis != 2:
            raise _Reject(
                "a countersink, chamfer or draft on an extruded outline is not supported yet",
                6,
            )
        for region in cone_regions:
            cone = region.surface
            if (region.area < _min_area(m) or _principal_axis(cone.axis) != 2
                    or _cone_outward(m, region, cone)):
                raise _Reject(
                    "a countersink, chamfer or draft on an extruded outline is not supported yet",
                    6,
                )
            hole_cones.append(region)
    try:
        holes, bosses = _cylinder_features(
            m, regions, tol, axis, origin, length, outline, cutouts, hole_cones, fillet_ids)
    except Unsupported as e:
        raise _Reject(str(e), 6) from e
    profile = {
        "axis": _AXIS_NAME[axis],
        "origin": float(origin),
        "length": float(length),
        "outline": outline,
        "cutouts": cutouts,
    }
    if fillet_radius is not None:
        profile["fillet_r"] = float(fillet_radius)
    return {
        "profile": profile,
        "holes": holes,
        "bosses": bosses,
        "pads": [],
    }


def _summarize(failures: list[tuple[int, _Reject]]) -> str:
    if not failures:
        return "not a straight extrusion along X, Y, or Z"
    best = max(failures, key=lambda t: t[1].rank)
    if best[1].rank >= 4:
        axis = _AXIS_NAME[best[0]]
        return f"not a straight extrusion along {axis}: {best[1]}"
    # every axis failed while classifying faces. Say what each one ran into.
    bits = []
    for axis, err in failures:
        if err.rank <= 1:
            bits.append(f"along {_AXIS_NAME[axis]}, {err}")
        else:
            bits.append(f"along {_AXIS_NAME[axis]}, {err}")
    if best[1].rank <= 1 and all(err.rank <= 1 for _a, err in failures):
        return ("not a straight extrusion along X, Y, or Z. "
                + "; ".join(bits)
                + ". The cross-section changes along the part")
    return "not a straight extrusion. " + "; ".join(bits)


def extract_extrusion(m: trimesh.Trimesh, regions, tol: float) -> dict:
    """Recognise one straight extrusion. Raises :class:`Unsupported` otherwise."""
    if any(r.surface is None for r in regions):
        raise Unsupported("part has unfitted (faceted) regions")
    # A noisy fillet can come back as a tilted cone. If the vertices lie on a
    # principal cylinder, use that for this fit and put the original surface
    # back afterwards. A curve that is not that cylinder still refuses.
    replaced = []
    try:
        for r in regions:
            s = r.surface
            if isinstance(s, Sphere):
                raise Unsupported("sphere faces are not supported yet")
            if not isinstance(s, (Plane, Cylinder, Cone)):
                raise Unsupported(f"{s.kind} faces are not supported yet")
            if isinstance(s, (Cylinder, Cone)) and _principal_axis(s.axis) is None:
                pts = m.vertices[np.unique(m.faces[r.faces])] if len(r.faces) else np.zeros((0, 3))
                cyl = fit_principal_cylinder(pts, tol)
                if cyl is None:
                    raise Unsupported("curved surface whose axis is not aligned with X, Y or Z")
                replaced.append((r, s))
                r.surface = cyl
        failures = []
        for axis in (0, 1, 2):
            try:
                return _along(m, regions, tol, axis)
            except _Reject as e:
                failures.append((axis, e))
        raise Unsupported(_summarize(failures))
    finally:
        for region, old in replaced:
            region.surface = old
