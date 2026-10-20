"""Straight extrusions of a 2D outline, for parametric CadQuery and OpenSCAD.

A part qualifies when it is a prism along X, Y, or Z: two end caps, side
walls parallel to that axis, optional polygonal through-openings in the
outline, and optional round through-holes. The outline vertices and the
extrusion length become named parameters. The box/pad vocabulary stays in
charge for parts it already explains; this module is the fallback.
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
import trimesh

from .fit import Cone, Cylinder, Plane, Sphere
from .parametric import Unsupported, _nice, _principal_axis, _snap
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


def _loop_polygon(mesh, loop, walls, U, V, tol):
    """Corners of one boundary loop from the wall planes it runs along.

    Returns (points, "polygon") or (points2d, "curve") when the loop does not
    lie on walls (a round hole). None when the loop is too broken to use.
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


def _polygons_from_cap(mesh, face_ids, walls, U, V, tol):
    loops = _loops(mesh, face_ids)
    if not loops:
        return None
    built = []
    for loop in loops:
        got = _loop_polygon(mesh, loop, walls, U, V, tol)
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


def _same_loops(low, high, tol):
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


def _cylinder_features(m, regions, tol, axis, origin, length, outline, cutouts):
    """Through-holes whose axis is X, Y, or Z. Bosses and partial rounds are refused."""
    U, V = _PROFILE_AXES[axis]
    holes = []
    for r in regions:
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
            raise Unsupported("a boss on an extruded outline is not supported yet")
        verts = m.vertices[np.unique(m.faces[r.faces])]
        along = verts[:, ax]
        za, zb = float(along.min()), float(along.max())
        name = _AXIS_NAME[ax]
        pt = np.asarray(s.point, float)
        if ax == axis:
            uv = (float(pt @ U), float(pt @ V))
            if not _inside_profile(uv, outline, cutouts):
                raise Unsupported("hole is outside the extruded outline")
            covered = abs(za - origin) <= _eps(tol) and abs(zb - (origin + length)) <= _eps(tol)
            if not covered:
                raise Unsupported("hole along the extrusion does not pass through both ends")
            span_kind = "extrude"
            feat = {"d": 2.0 * float(s.radius), "axis": name, "span": span_kind}
            for i, coord in enumerate("xyz"):
                if i != ax:
                    feat[coord] = float(pt[i])
            holes.append(feat)
            continue
        # perpendicular to the extrusion: the profile coordinate that is not
        # the hole axis is fixed, and the hole has to cross every solid interval
        profile_axes = [i for i in range(3) if i != axis]
        # the hole varies along `ax`, so the fixed profile axis is the other one
        fixed_world = profile_axes[0] if profile_axes[1] == ax else profile_axes[1]
        fixed_index = profile_axes.index(fixed_world)  # 0 or 1 in (u, v)
        value = float(pt[fixed_world])
        intervals = _material_intervals(outline, cutouts, fixed_index, value)
        if not intervals:
            raise Unsupported("hole is outside the extruded outline")
        lim = _eps(tol, 0.35)
        if not all(za - lim <= a and b <= zb + lim for a, b in intervals):
            raise Unsupported("hole does not pass through the extruded profile")
        # also the cylinder itself must lie inside the extrusion's length
        # (its axis point's extrusion coordinate is between the caps)
        if not (origin - lim <= float(pt[axis]) <= origin + length + lim):
            raise Unsupported("hole is outside the extruded profile")
        vary_index = 1 - fixed_index
        span_kind = ("outline", vary_index)
        feat = {"d": 2.0 * float(s.radius), "axis": name, "span": span_kind}
        for i, coord in enumerate("xyz"):
            if i != ax:
                feat[coord] = float(pt[i])
        holes.append(feat)
    holes.sort(key=lambda h: (h["axis"], h.get("x", 0.0), h.get("y", 0.0), h.get("z", 0.0)))
    return holes


def profile_cutter(hole: dict, lang: str) -> tuple[str, str]:
    """``(start, length)`` of a through-hole cutter along its axis, past the ends by 1."""
    span = hole["span"]
    if span == "extrude":
        return "extrude_0 - 1", "extrude_len + 2"
    idx = span[1]
    if lang == "cadquery":
        lo = f"min(p[{idx}] for p in outline)"
        hi = f"max(p[{idx}] for p in outline)"
    else:
        lo = f"min([for (p = outline) p[{idx}]])"
        hi = f"max([for (p = outline) p[{idx}]])"
    return f"({lo} - 1)", f"({hi} - ({lo}) + 2)"


def profile_parameters(feats: dict, tol: float) -> list[tuple[str, float, str]]:
    """Float parameters of an extruded profile, in emitter order."""
    pr = feats["profile"]
    ax = pr["axis"]
    P = [
        ("extrude_len", _nice(pr["length"], tol), f"length along {ax}"),
        ("extrude_0", _nice(pr["origin"], tol), f"start along {ax}"),
    ]
    for i, d in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        P.append((f"{n}_d", _nice(d["d"], tol), f"diameter (axis {d['axis']})"))
        for coord in "xyz":
            if coord != d["axis"].lower() and coord in d:
                P.append((f"{n}_{coord}", _nice(d[coord], tol), ""))
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
    if abs(float(extent.min()) - origin) > gap or abs(float(extent.max()) - end) > gap:
        raise _Reject("the solid sticks out past the two end faces", 3)
    wall_planes = _merge_walls([r.surface for r in walls], tol)
    if len(wall_planes) < 3:
        raise _Reject("not enough side walls to form an outline", 2)
    U, V = _PROFILE_AXES[axis]
    low_poly = _polygons_from_cap(m, _cap_faces(m, low, axis), wall_planes, U, V, tol)
    high_poly = _polygons_from_cap(m, _cap_faces(m, high, axis), wall_planes, U, V, tol)
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
    # a round hole is a non-polygonal loop on a cylinder along the extrusion
    axial = [r for r in cylinders if _principal_axis(r.surface.axis) == axis]
    axis_vec = np.zeros(3)
    axis_vec[axis] = 1.0
    for pts in curve_loops:
        cloud = np.array([p[0] * U + p[1] * V + origin * axis_vec for p in pts])
        if not any(_on_cylinder(cloud, r.surface, axis, tol) for r in axial):
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
    if any(isinstance(r.surface, Cone) for r in regions):
        raise _Reject(
            "a countersink, chamfer or draft on an extruded outline is not supported yet",
            6,
        )
    try:
        holes = _cylinder_features(m, regions, tol, axis, origin, length, outline, cutouts)
    except Unsupported as e:
        raise _Reject(str(e), 6) from e
    return {
        "profile": {
            "axis": _AXIS_NAME[axis],
            "origin": float(origin),
            "length": float(length),
            "outline": outline,
            "cutouts": cutouts,
        },
        "holes": holes,
        "bosses": [],
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
    for r in regions:
        s = r.surface
        if isinstance(s, Sphere):
            raise Unsupported("sphere faces are not supported yet")
        if not isinstance(s, (Plane, Cylinder, Cone)):
            raise Unsupported(f"{s.kind} faces are not supported yet")
        if isinstance(s, (Cylinder, Cone)) and _principal_axis(s.axis) is None:
            raise Unsupported("curved surface whose axis is not aligned with X, Y or Z")
    failures = []
    for axis in (0, 1, 2):
        try:
            return _along(m, regions, tol, axis)
        except _Reject as e:
            failures.append((axis, e))
    raise Unsupported(_summarize(failures))
