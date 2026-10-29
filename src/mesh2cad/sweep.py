"""Segmented extrusions for a cross-section that changes along one axis.

Used only after the silhouette intersection, the partial-depth union, and the
wall-snapped rebuild have all missed the mesh. A full-length outline is the
union of every station, so a step, a moving tooth, or a pocket that occupies
only part of the length is filled in. This rebuilds the solid as one bounding
extrusion with a cut for each span where the true section stays still. A thin
gap between two parallel faces is its own polygonal cut, and so is an inner
loop that keeps its shape for a run of stations (a hex pocket, a hole).

The gate is unchanged. A part the earlier solids already accept never reaches
this module.
"""
from __future__ import annotations

import numpy as np

from .csg import (
    _ok, _point_in_ring, _score, _shoelace, _simplify_ring, _solid,
)
from .fit import Plane
from .parametric import _nice
from .stack import _polyline_distance

# How far a section may drift inside one segment, and how often stations are
# taken. Both sit under the export gate (max(4*tol, 0.1)) with room for the
# 0.02 mm overlap the cut solids use.
_DRIFT = 0.05
_PITCH = 0.12
_SIMPLIFY = 0.015


def _profile_axes(axis: int) -> tuple[int, int]:
    return {0: (1, 2), 1: (0, 2), 2: (0, 1)}[axis]


def _section_loops(mesh, axis: int, coord: float, lim: float):
    origin = [0.0, 0.0, 0.0]
    normal = [0.0, 0.0, 0.0]
    origin[axis] = float(coord)
    normal[axis] = 1.0
    sec = mesh.section(plane_origin=origin, plane_normal=normal)
    if sec is None:
        return []
    u, v = _profile_axes(axis)
    loops = []
    for path in sec.discrete:
        uv = [(float(q[u]), float(q[v])) for q in path]
        if len(uv) >= 2 and abs(uv[0][0] - uv[-1][0]) < 1e-3 and abs(uv[0][1] - uv[-1][1]) < 1e-3:
            uv = uv[:-1]
        simp = _simplify_ring(uv, lim)
        if len(simp) < 3 or abs(_shoelace(simp)) < 0.3:
            continue
        if _shoelace(simp) < 0:
            simp = list(reversed(simp))
        loops.append(simp)
    loops.sort(key=lambda s: -abs(_shoelace(s)))
    return loops


def _nest(loops):
    """Outers with the inner loops that sit inside them. Largest first."""
    groups = []
    for i, outer in enumerate(loops):
        cx = sum(p[0] for p in outer) / len(outer)
        cy = sum(p[1] for p in outer) / len(outer)
        if any(
            j != i and abs(_shoelace(other)) > abs(_shoelace(outer)) and _point_in_ring(cx, cy, other)
            for j, other in enumerate(loops)
        ):
            continue
        holes = []
        for j, other in enumerate(loops):
            if j == i:
                continue
            hx = sum(p[0] for p in other) / len(other)
            hy = sum(p[1] for p in other) / len(other)
            if _point_in_ring(hx, hy, outer) and abs(_shoelace(other)) < abs(_shoelace(outer)):
                holes.append(other)
        groups.append((outer, holes))
    return groups


def _center(loop):
    return (sum(p[0] for p in loop) / len(loop), sum(p[1] for p in loop) / len(loop))


def _plane_coord(region, axis: int) -> float | None:
    n = np.asarray(region.surface.normal, float)
    length = float(np.linalg.norm(n))
    if length < 1e-9 or abs(n[axis]) / length < 0.98:
        return None
    return float(region.surface.d) / length / (n[axis] / length)


def _slits(mesh, regions, axis: int):
    """Voids a fraction of a millimetre thick, outlined by the section between the walls.

    A 0.2 mm gap is narrower than the stations above, so a sample can land on
    the wall to either side and never see it. The two faces name the gap.
    """
    u, v = _profile_axes(axis)
    faces = []
    for region in regions:
        if not isinstance(region.surface, Plane) or region.area < 15:
            continue
        coord = _plane_coord(region, axis)
        if coord is None:
            continue
        verts = mesh.vertices[np.unique(mesh.faces[region.faces])]
        faces.append((
            float(coord), float(region.area),
            float(verts[:, u].min()), float(verts[:, u].max()),
            float(verts[:, v].min()), float(verts[:, v].max()),
        ))
    faces.sort()
    found = []
    for i, a in enumerate(faces):
        for b in faces[i + 1:]:
            gap = b[0] - a[0]
            if gap < 0.12:
                continue
            if gap > 0.55:
                break
            if min(a[3], b[3]) - max(a[2], b[2]) < 2 or min(a[5], b[5]) - max(a[4], b[4]) < 2:
                continue
            mid = 0.5 * (a[0] + b[0])
            loops = _section_loops(mesh, axis, mid, 0.03)
            target = min(a[1], b[1])
            best = None
            for loop in loops[1:]:
                area = abs(_shoelace(loop))
                if best is None or abs(area - target) < abs(abs(_shoelace(best)) - target):
                    best = loop
            if best is None or abs(abs(_shoelace(best)) - target) > 0.25 * target:
                continue
            cu, cv = _center(best)
            probe = [0.0, 0.0, 0.0]
            probe[axis] = mid
            probe[u], probe[v] = cu, cv
            if bool(mesh.contains([probe])[0]):
                continue
            found.append((a[0], b[0], best))
            break
    return found


def _in_slit(x, loop, slits, pitch):
    cu, cv = _center(loop)
    area = abs(_shoelace(loop))
    for a0, b0, outline in slits:
        if not (a0 - pitch <= x <= b0 + pitch):
            continue
        if abs(area - abs(_shoelace(outline))) > 8:
            continue
        ou, ov = _center(outline)
        if abs(cu - ou) < 2 and abs(cv - ov) < 2:
            return True
    return False


def _void_tracks(xs, samples, slits, pitch):
    """Inner loops that keep their place from one station to the next."""
    tracks = []
    for x, groups in zip(xs, samples):
        holes = []
        for _outer, inners in groups:
            for loop in inners:
                if not _in_slit(x, loop, slits, pitch):
                    holes.append(loop)
        for loop in holes:
            cu, cv = _center(loop)
            area = abs(_shoelace(loop))
            placed = False
            for track in tracks:
                if x - track["xs"][-1] > pitch * 1.6:
                    continue
                pu, pv = _center(track["loops"][-1])
                prev = abs(_shoelace(track["loops"][-1]))
                close = abs(cu - pu) < 1.2 and abs(cv - pv) < 1.2
                similar = abs(area - prev) < max(0.35 * prev, 2.0)
                if close and similar:
                    track["xs"].append(x)
                    track["loops"].append(loop)
                    placed = True
                    break
            if not placed:
                tracks.append({"xs": [x], "loops": [loop]})
    cuts = []
    for track in tracks:
        span = track["xs"][-1] - track["xs"][0]
        if len(track["xs"]) < 2 and span < 0.4 and abs(_shoelace(track["loops"][0])) < 3:
            continue
        loop = track["loops"][len(track["loops"]) // 2]
        cuts.append((track["xs"][0] - pitch / 2, track["xs"][-1] + pitch / 2, loop))
    return cuts


def _axis_order(mesh) -> list[int]:
    """The axis whose section area moves the most, first."""
    scores = []
    for axis in range(3):
        lo = float(mesh.bounds[0, axis])
        hi = float(mesh.bounds[1, axis])
        areas = []
        for t in (0.2, 0.4, 0.6, 0.8):
            loops = _section_loops(mesh, axis, lo + t * (hi - lo), 0.05)
            if loops:
                areas.append(abs(_shoelace(loops[0])))
        if len(areas) < 2:
            scores.append((0.0, axis))
            continue
        mean = sum(areas) / len(areas)
        var = sum((a - mean) ** 2 for a in areas) / len(areas)
        scores.append((var / max(mean, 1.0), axis))
    scores.sort(reverse=True)
    return [axis for _score, axis in scores]


def _round_loop(loop, step):
    return [(_nice(p[0], step), _nice(p[1], step)) for p in loop]


def _sweep_axis(mesh, regions, axis: int, tol: float):
    lo = float(mesh.bounds[0, axis])
    hi = float(mesh.bounds[1, axis])
    length = hi - lo
    if length < 0.5:
        return None
    pitch = _PITCH if length / _PITCH <= 200 else length / 200.0
    xs = list(np.arange(lo + pitch / 2.0, hi, pitch))
    if len(xs) < 2:
        return None
    samples = [_nest(_section_loops(mesh, axis, float(x), _SIMPLIFY)) for x in xs]
    if not any(samples):
        return None
    slits = _slits(mesh, regions, axis)
    voids = _void_tracks(xs, samples, slits, pitch)

    u, v = _profile_axes(axis)
    umin = float(mesh.bounds[0, u]) - 0.4
    umax = float(mesh.bounds[1, u]) + 0.4
    vmin = float(mesh.bounds[0, v]) - 0.4
    vmax = float(mesh.bounds[1, v]) + 0.4
    # a section vertex past the mesh bounds (a chord) still has to sit inside
    for groups in samples:
        for outer, _holes in groups:
            for p, q in outer:
                umin = min(umin, p - 0.4)
                umax = max(umax, p + 0.4)
                vmin = min(vmin, q - 0.4)
                vmax = max(vmax, q + 0.4)
    box = [(umin, vmin), (umax, vmin), (umax, vmax), (umin, vmax)]

    step = min(max(tol, 1e-4), 0.01)
    cuts = []
    for a0, b0, loop in slits:
        cuts.append({
            "axis": "XYZ"[axis],
            "axis_i": axis,
            "origin": _nice(a0, step),
            "length": _nice(b0 - a0, step),
            "outline": _round_loop(loop, step),
            "cutouts": [],
        })
    for a0, b0, loop in voids:
        area = abs(_shoelace(loop))
        # a pocket smaller than the slit faces stops just short of a tab that
        # sets the end of the span; the slice overlap would otherwise leave a skin
        if area < 40:
            a0 -= 0.04
            b0 += 0.04
        cuts.append({
            "axis": "XYZ"[axis],
            "axis_i": axis,
            "origin": _nice(a0, step),
            "length": _nice(max(b0 - a0, 0.05), step),
            "outline": _round_loop(loop, step),
            "cutouts": [],
        })

    outers = [[o for o, _h in groups] for groups in samples]
    i = 0
    segments = 0
    while i < len(xs):
        j = i + 1
        while j < len(xs):
            a, b = outers[i], outers[j]
            if len(a) != len(b) or not a:
                break
            if any(_polyline_distance(np.asarray(p, float), np.asarray(q, float)) > _DRIFT
                   for p, q in zip(a, b)):
                break
            j += 1
        if not outers[i]:
            i = j
            continue
        mid = outers[(i + j - 1) // 2]
        origin = lo if i == 0 else float(xs[i] - pitch / 2.0)
        end = hi if j == len(xs) else float(xs[j - 1] + pitch / 2.0)
        cuts.append({
            "axis": "XYZ"[axis],
            "axis_i": axis,
            "origin": _nice(origin, step),
            "length": _nice(max(end - origin, pitch), step),
            "outline": _round_loop(box, step),
            "cutouts": [_round_loop(loop, step) for loop in mid],
        })
        segments += 1
        i = j
    if segments < 2:
        return None
    positive = [{
        "axis": "XYZ"[axis],
        "axis_i": axis,
        "origin": _nice(lo, step),
        "length": _nice(hi - lo, step),
        "outline": _round_loop(box, step),
        "cutouts": [],
    }]
    try:
        shape = _solid(positive, cuts, [], [], mode="union")
    except Exception:
        return None
    score = _score(shape, mesh, tol, 4000)
    if score["max_deviation"] >= 1e8:
        return None
    result = {
        "csg": {
            "positive": [{k: p[k] for k in ("axis", "origin", "length", "outline", "cutouts")}
                         for p in positive],
            "cuts": [{k: c[k] for k in ("axis", "origin", "length", "outline", "cutouts")}
                     for c in cuts],
            "clips": [],
            "combine": "union",
        },
        "holes": [],
    }
    return score, result


def try_sweep(mesh, regions, tol: float):
    """The best segmented sweep, or None when no axis produces a solid.

    The caller keeps it only when it passes the same gate as every other
    script. A closer miss replaces the reported deviation.
    """
    best = None
    for axis in _axis_order(mesh):
        got = _sweep_axis(mesh, regions, axis, tol)
        if got is None:
            continue
        score, result = got
        if _ok(score, tol):
            return score, result
        if best is None or score["max_deviation"] < best[0]["max_deviation"]:
            best = (score, result)
        # one axis already explains the part to a millimetre; the others are
        # the long way around the same solid
        if score["max_deviation"] < 1.0:
            break
    return best
