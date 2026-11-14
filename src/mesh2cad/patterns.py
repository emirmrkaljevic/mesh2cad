"""Regular patterns of plain Z through-holes for the parametric exporters.

A row is three or more holes of one diameter whose centers lie on a straight
line and on one pitch. Two holes are left as ``hole_N_*`` (or ``base_hole_N_*``
/ ``leg_hole_N_*``). They have a single gap, so equal spacing is automatic,
and moving one of them should not drag the other.

A rectangular grid is an axis-aligned block, at least two by two, same
diameter, constant pitch along X and along Y. A bolt circle is three or more
holes of one diameter whose centers lie on one circle and are equally spaced
in angle. ``hole_circle_N_r`` is the pitch-circle radius, not the diameter.
Blind holes, countersinks, counterbores, and holes whose axis is X or Y stay
individual parameters.

:func:`apply_hole_patterns` rewrites a feature dict. The exporters try that
script first. If it misses the usual deviation and volume gate, they write
the individual holes instead.
"""
from __future__ import annotations

import math
from itertools import combinations

import numpy as np

_PLAIN = {"x", "y", "d", "axis"}


def _pos_eps(tol: float) -> float:
    """How far a center may sit from the ideal row, grid, or bolt circle."""
    return max(4.0 * float(tol), 0.25)


def _diam_eps(tol: float) -> float:
    return max(4.0 * float(tol), 0.25)


def _is_plain_z(hole: dict) -> bool:
    """A through-hole along Z with no countersink, counterbore, or fillet."""
    if hole.get("axis", "Z") != "Z":
        return False
    if "x" not in hole or "y" not in hole or "d" not in hole:
        return False
    return set(hole).issubset(_PLAIN)


def _by_diameter(holes: list[dict], tol: float) -> list[list[dict]]:
    eps = _diam_eps(tol)
    ordered = sorted(holes, key=lambda h: float(h["d"]))
    groups: list[list[dict]] = []
    for hole in ordered:
        if groups and abs(float(hole["d"]) - float(groups[-1][0]["d"])) <= eps:
            groups[-1].append(hole)
        else:
            groups.append([hole])
    return groups


def _lattice_fit(coords, eps: float):
    """``(start, pitch, residual)`` when ``coords[i]`` is ``start + i * pitch``.

    ``coords`` must already be sorted. ``None`` when the pitch is tiny or any
    point misses the line by more than ``eps``.
    """
    values = np.asarray(coords, float)
    n = int(values.size)
    if n < 2:
        return None
    idx = np.arange(n, dtype=float)
    i_mean = float(idx.mean())
    t_mean = float(values.mean())
    var = float(np.sum((idx - i_mean) ** 2))
    if var <= 0.0:
        return None
    pitch = float(np.sum((idx - i_mean) * (values - t_mean)) / var)
    if pitch <= eps:
        return None
    start = t_mean - pitch * i_mean
    resid = float(np.max(np.abs(values - (start + idx * pitch))))
    if resid > eps:
        return None
    return start, pitch, resid


def _clusters(values: dict[int, float], ids: list[int], eps: float) -> list[dict]:
    """Group ``ids`` whose coordinates fall within ``eps`` of the running end."""
    groups: list[dict] = []
    for i in sorted(ids, key=lambda k: values[k]):
        if groups and values[i] - groups[-1]["hi"] <= eps:
            groups[-1]["ids"].append(i)
            groups[-1]["hi"] = values[i]
        else:
            groups.append({"ids": [i], "hi": values[i]})
    for group in groups:
        group["mean"] = sum(values[i] for i in group["ids"]) / len(group["ids"])
    return groups


def _best_grid(holes: list[dict], eps: float):
    """Largest full rectangle in ``holes``, or None.

    Columns share an X coordinate and rows share a Y. Both pitches are
    constant. Every cell of the block has exactly one hole.
    """
    n = len(holes)
    if n < 4:
        return None
    ids = list(range(n))
    xs = {i: float(holes[i]["x"]) for i in ids}
    ys = {i: float(holes[i]["y"]) for i in ids}
    x_groups = _clusters(xs, ids, eps)
    y_groups = _clusters(ys, ids, eps)
    if len(x_groups) < 2 or len(y_groups) < 2:
        return None
    best = None
    for a in range(len(x_groups)):
        for b in range(a + 1, len(x_groups)):
            cols = x_groups[a:b + 1]
            fit_x = _lattice_fit([col["mean"] for col in cols], eps)
            if fit_x is None:
                continue
            for c in range(len(y_groups)):
                for d in range(c + 1, len(y_groups)):
                    rows = y_groups[c:d + 1]
                    fit_y = _lattice_fit([row["mean"] for row in rows], eps)
                    if fit_y is None:
                        continue
                    y_sets = [set(row["ids"]) for row in rows]
                    cells: list[int] = []
                    ok = True
                    for col in cols:
                        col_ids = set(col["ids"])
                        for y_ids in y_sets:
                            both = col_ids & y_ids
                            if len(both) != 1:
                                ok = False
                                break
                            cells.append(next(iter(both)))
                        if not ok:
                            break
                    if not ok:
                        continue
                    start_x, pitch_x, _rx = fit_x
                    start_y, pitch_y, _ry = fit_y
                    diam = float(np.median([holes[i]["d"] for i in cells]))
                    if pitch_x <= max(0.5 * diam, eps) or pitch_y <= max(0.5 * diam, eps):
                        continue
                    resid = 0.0
                    ix = 0
                    for col in cols:
                        iy = 0
                        col_ids = set(col["ids"])
                        for y_ids in y_sets:
                            i = next(iter(col_ids & y_ids))
                            err = math.hypot(
                                float(holes[i]["x"]) - (start_x + ix * pitch_x),
                                float(holes[i]["y"]) - (start_y + iy * pitch_y),
                            )
                            resid = max(resid, err)
                            iy += 1
                        ix += 1
                    if resid > eps:
                        continue
                    spec = {
                        "nx": len(cols),
                        "ny": len(rows),
                        "pitch_x": float(pitch_x),
                        "pitch_y": float(pitch_y),
                        "d": diam,
                        "x": float(start_x),
                        "y": float(start_y),
                        "angle": 0.0,
                    }
                    score = (len(cells), -resid)
                    if best is None or score > best[0]:
                        best = (score, spec, {id(holes[i]) for i in cells})
    if best is None:
        return None
    return best[1], best[2]


def _row_from_selection(pts: np.ndarray, holes: list[dict], eps: float):
    """One equal-pitch line through every point, or None."""
    center = pts.mean(axis=0)
    centered = pts - center
    if float(np.linalg.norm(centered)) <= eps:
        return None
    _u, _s, vt = np.linalg.svd(centered, full_matrices=False)
    direction = np.asarray(vt[0], float)
    rel = pts - center
    along = rel @ direction
    perp = np.linalg.norm(rel - along[:, None] * direction, axis=1)
    if float(perp.max()) > eps:
        return None
    order = np.argsort(along)
    fit = _lattice_fit(along[order], eps)
    if fit is None:
        return None
    start_t, pitch, resid_t = fit
    diam = float(np.median([h["d"] for h in holes]))
    if pitch <= max(0.5 * diam, eps):
        return None
    origin = center + direction * start_t
    far = center + direction * (start_t + pitch * (len(pts) - 1))
    ordered = [holes[int(i)] for i in order]
    # Travel toward +X, or toward +Y when the row is vertical, so a line
    # along Y is 90 degrees and not -90 from the other end.
    dx, dy = float(direction[0]), float(direction[1])
    if abs(dx) >= abs(dy):
        flip = dx < 0.0
    else:
        flip = dy < 0.0
    if flip:
        direction = -direction
        origin = far
        ordered = list(reversed(ordered))
    angle = math.degrees(math.atan2(float(direction[1]), float(direction[0])))
    ca = math.cos(math.radians(angle))
    sa = math.sin(math.radians(angle))
    resid = max(float(perp.max()), float(resid_t))
    ox, oy = float(origin[0]), float(origin[1])
    for i, hole in enumerate(ordered):
        x = ox + i * pitch * ca
        y = oy + i * pitch * sa
        resid = max(resid, math.hypot(x - float(hole["x"]), y - float(hole["y"])))
    if resid > eps:
        return None
    return {
        "n": len(holes),
        "pitch": float(pitch),
        "d": diam,
        "x": ox,
        "y": oy,
        "angle": float(angle),
        "_resid": resid,
    }


def _best_row(holes: list[dict], eps: float):
    """Longest equal-pitch line of at least three holes, or None."""
    n = len(holes)
    if n < 3:
        return None
    pts = np.array([[float(h["x"]), float(h["y"])] for h in holes], float)
    best = None
    for a, b in combinations(range(n), 2):
        delta = pts[b] - pts[a]
        span = float(np.hypot(delta[0], delta[1]))
        if span <= eps:
            continue
        direction = delta / span
        rel = pts - pts[a]
        along = rel @ direction
        perp = np.linalg.norm(rel - along[:, None] * direction, axis=1)
        on = [i for i in range(n) if float(perp[i]) <= eps]
        if len(on) < 3:
            continue
        on.sort(key=lambda i: float(along[i]))
        for i in range(len(on)):
            for j in range(i + 2, len(on)):
                sel = on[i:j + 1]
                spec = _row_from_selection(pts[sel], [holes[k] for k in sel], eps)
                if spec is None:
                    continue
                score = (spec["n"], -spec["_resid"])
                if best is None or score > best[0]:
                    chosen = dict(spec)
                    chosen.pop("_resid", None)
                    best = (score, chosen, {id(holes[k]) for k in sel})
    if best is None:
        return None
    return best[1], best[2]


def _take_grids(holes: list[dict], eps: float):
    remaining = list(holes)
    found = []
    used: set[int] = set()
    while len(remaining) >= 4:
        got = _best_grid(remaining, eps)
        if got is None:
            break
        spec, idset = got
        found.append(spec)
        used |= idset
        remaining = [h for h in remaining if id(h) not in idset]
    return found, used


def _take_rows(holes: list[dict], eps: float):
    remaining = list(holes)
    found = []
    used: set[int] = set()
    while len(remaining) >= 3:
        got = _best_row(remaining, eps)
        if got is None:
            break
        spec, idset = got
        found.append(spec)
        used |= idset
        remaining = [h for h in remaining if id(h) not in idset]
    return found, used


def _circumcircle(p, q, s):
    """``(cx, cy, r)`` through three points, or None when they are collinear."""
    d = 2.0 * (p[0] * (q[1] - s[1]) + q[0] * (s[1] - p[1]) + s[0] * (p[1] - q[1]))
    if abs(d) < 1e-9:
        return None
    p2 = p[0] * p[0] + p[1] * p[1]
    q2 = q[0] * q[0] + q[1] * q[1]
    s2 = s[0] * s[0] + s[1] * s[1]
    cx = (p2 * (q[1] - s[1]) + q2 * (s[1] - p[1]) + s2 * (p[1] - q[1])) / d
    cy = (p2 * (s[0] - q[0]) + q2 * (p[0] - s[0]) + s2 * (q[0] - p[0])) / d
    radius = math.hypot(cx - p[0], cy - p[1])
    if not math.isfinite(radius) or radius < 1e-6 or radius > 1e6:
        return None
    return float(cx), float(cy), float(radius)


def _fit_circle(pts: np.ndarray):
    """Least-squares circle ``(cx, cy, r)``, or None when the fit is degenerate."""
    x = np.asarray(pts[:, 0], float)
    y = np.asarray(pts[:, 1], float)
    design = np.column_stack((2.0 * x, 2.0 * y, np.ones(len(pts))))
    target = x * x + y * y
    sol, *_ = np.linalg.lstsq(design, target, rcond=None)
    cx, cy, c = float(sol[0]), float(sol[1]), float(sol[2])
    r2 = c + cx * cx + cy * cy
    if r2 <= 1e-12 or not math.isfinite(r2):
        return None
    return cx, cy, math.sqrt(r2)


def _fold_start(angle: float, step: float) -> float:
    """Fold ``angle`` into ``(-step/2, step/2]`` by whole steps."""
    half = step / 2.0
    start = float(angle)
    for _ in range(8):
        if start > half:
            start -= step
        elif start <= -half:
            start += step
        else:
            break
    if abs(start) < 1e-10:
        return 0.0
    return start


def _regular_circle(pts: np.ndarray, holes: list[dict], eps: float):
    """One equally spaced bolt circle, or None.

    The start angle is the polar angle of one hole, in degrees, 0 along +X,
    folded into ``(-180/n, 180/n]`` so the same holes have one value.
    """
    n = int(len(pts))
    if n < 3:
        return None
    fit = _fit_circle(pts)
    if fit is None:
        return None
    cx, cy, radius = fit
    if radius <= eps:
        return None
    rel = pts - np.array([cx, cy])
    dist = np.hypot(rel[:, 0], rel[:, 1])
    if float(np.max(np.abs(dist - radius))) > eps:
        return None
    ang = np.arctan2(rel[:, 1], rel[:, 0])
    order = np.argsort(ang)
    ang_s = ang[order]
    step = 2.0 * math.pi / n
    gaps = np.diff(ang_s)
    wrap = (ang_s[0] + 2.0 * math.pi) - ang_s[-1]
    gaps = np.append(gaps, wrap)
    # A missing hole leaves one double gap. The Euclidean check below is the
    # real gate; this only rejects a spacing that cannot land inside ``eps``.
    if float(np.max(np.abs(gaps - step))) > max(4.0 * eps / radius, 0.15):
        return None
    start = _fold_start(float(ang_s[0]), step)
    resid = 0.0
    used: set[int] = set()
    for point, theta in zip(pts, ang):
        slot = int(round((float(theta) - start) / step)) % n
        if slot in used:
            return None
        used.add(slot)
        ideal = start + slot * step
        err = math.hypot(
            float(point[0]) - (cx + radius * math.cos(ideal)),
            float(point[1]) - (cy + radius * math.sin(ideal)),
        )
        if err > eps:
            return None
        resid = max(resid, err)
    if len(used) != n:
        return None
    diam = float(np.median([h["d"] for h in holes]))
    chord = 2.0 * radius * math.sin(math.pi / n)
    if chord <= max(0.5 * diam, eps):
        return None
    angle = math.degrees(start)
    if abs(angle) < 1e-8:
        angle = 0.0
    return {
        "n": n,
        "d": diam,
        "r": float(radius),
        "cx": float(cx),
        "cy": float(cy),
        "angle": float(angle),
        "_resid": float(resid),
    }


def _circle_inliers(pts: np.ndarray, seed, eps: float):
    """Indices of every point within ``eps`` of the circle refit from ``seed``.

    The seed circle can be a noisy triple. A wider first pass collects the
    pattern, then the fit is tightened to ``eps``.
    """
    cx, cy, radius = seed
    wide = max(3.0 * eps, eps)
    radial = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - radius)
    idx = [i for i in range(len(pts)) if float(radial[i]) <= wide]
    if len(idx) < 3:
        return None
    fit = _fit_circle(pts[idx])
    if fit is None:
        return None
    cx, cy, radius = fit
    radial = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - radius)
    strict = tuple(i for i in range(len(pts)) if float(radial[i]) <= eps)
    if len(strict) < 3:
        return None
    fit2 = _fit_circle(pts[list(strict)])
    if fit2 is None:
        return None
    cx, cy, radius = fit2
    radial = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - radius)
    final = tuple(i for i in range(len(pts)) if float(radial[i]) <= eps)
    if final != strict:
        return None
    return final


def _best_circle(holes: list[dict], eps: float):
    """Largest equally spaced circle of at least three holes, or None.

    Every hole that lies on that circle has to be in the pattern. A subset is
    not peeled off when another hole of the same diameter sits on the circle,
    so five holes of a six-hole pattern stay individual.
    """
    n = len(holes)
    if n < 3:
        return None
    pts = np.array([[float(h["x"]), float(h["y"])] for h in holes], float)
    best = None
    seen: set[tuple[int, ...]] = set()
    for comb in combinations(range(n), 3):
        circ = _circumcircle(pts[comb[0]], pts[comb[1]], pts[comb[2]])
        if circ is None:
            continue
        idx = _circle_inliers(pts, circ, eps)
        if idx is None or idx in seen:
            continue
        seen.add(idx)
        spec = _regular_circle(pts[list(idx)], [holes[i] for i in idx], eps)
        if spec is None:
            continue
        cx, cy, radius = spec["cx"], spec["cy"], spec["r"]
        radial = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - radius)
        chosen = set(idx)
        if any(i not in chosen and float(radial[i]) <= eps for i in range(n)):
            continue
        score = (spec["n"], -spec["_resid"])
        if best is None or score > best[0]:
            chosen_spec = dict(spec)
            chosen_spec.pop("_resid", None)
            best = (score, chosen_spec, {id(holes[i]) for i in idx})
    if best is None:
        return None
    return best[1], best[2]


def _take_circles(holes: list[dict], eps: float):
    remaining = list(holes)
    found = []
    used: set[int] = set()
    while len(remaining) >= 3:
        got = _best_circle(remaining, eps)
        if got is None:
            break
        spec, idset = got
        found.append(spec)
        used |= idset
        remaining = [h for h in remaining if id(h) not in idset]
    return found, used


def partition_holes(holes: list[dict], tol: float):
    """``(kept, rows, grids, circles)``. ``kept`` preserves the original order.

    Grids are taken first, so a 2 by 3 block is one grid and not two rows.
    Rows are next, so a line is not also treated as a very large circle.
    A bolt circle is only what is left.
    """
    consumed: set[int] = set()
    rows: list[dict] = []
    grids: list[dict] = []
    circles: list[dict] = []
    plain = [h for h in holes if _is_plain_z(h)]
    eps = _pos_eps(tol)
    for group in _by_diameter(plain, tol):
        found_grids, grid_ids = _take_grids(group, eps)
        grids.extend(found_grids)
        consumed |= grid_ids
        left = [h for h in group if id(h) not in grid_ids]
        found_rows, row_ids = _take_rows(left, eps)
        rows.extend(found_rows)
        consumed |= row_ids
        left = [h for h in left if id(h) not in row_ids]
        found_circles, circle_ids = _take_circles(left, eps)
        circles.extend(found_circles)
        consumed |= circle_ids
    kept = [h for h in holes if id(h) not in consumed]
    rows.sort(key=lambda r: (r["y"], r["x"], r["angle"], r["d"]))
    grids.sort(key=lambda g: (g["y"], g["x"], g["d"]))
    circles.sort(key=lambda c: (c["cy"], c["cx"], c["r"], c["angle"], c["d"]))
    return kept, rows, grids, circles


def apply_hole_patterns(feats: dict, tol: float) -> dict:
    """A copy of ``feats`` with plain Z through-holes grouped into patterns.

    Profiles and combined extrusions are returned unchanged: their exporters
    do not write pattern loops. The input dict is not modified.
    """
    if "profile" in feats or "csg" in feats:
        return feats
    out = None

    def copy() -> dict:
        nonlocal out
        if out is None:
            out = dict(feats)
        return out

    changed = False
    if "holes" in feats:
        kept, rows, grids, circles = partition_holes(list(feats["holes"]), tol)
        if rows or grids or circles:
            changed = True
            copied = copy()
            copied["holes"] = kept
            if rows:
                copied["hole_rows"] = rows
            if grids:
                copied["hole_grids"] = grids
            if circles:
                copied["hole_circles"] = circles
    bracket = feats.get("bracket")
    if bracket:
        rewritten = None
        for member in ("base", "leg"):
            key = f"{member}_holes"
            if not bracket.get(key):
                continue
            kept, rows, grids, circles = partition_holes(list(bracket[key]), tol)
            if not rows and not grids and not circles:
                continue
            changed = True
            if rewritten is None:
                rewritten = dict(bracket)
            rewritten[key] = kept
            if rows:
                rewritten[f"{member}_hole_rows"] = rows
            if grids:
                rewritten[f"{member}_hole_grids"] = grids
            if circles:
                rewritten[f"{member}_hole_circles"] = circles
        if rewritten is not None:
            copy()["bracket"] = rewritten
    if not changed:
        return feats
    return out


def has_hole_patterns(feats: dict) -> bool:
    if feats.get("hole_rows") or feats.get("hole_grids") or feats.get("hole_circles"):
        return True
    bracket = feats.get("bracket") or {}
    return any(bracket.get(name) for name in (
        "base_hole_rows", "base_hole_grids", "base_hole_circles",
        "leg_hole_rows", "leg_hole_grids", "leg_hole_circles",
    ))


def pattern_attempts(feats: dict, tol: float):
    """Patterned feature dict first, then the original holes.

    One yield when nothing grouped. Callers that can measure the solid keep
    the first attempt that passes the deviation gate.
    """
    patterned = apply_hole_patterns(feats, tol)
    if has_hole_patterns(patterned):
        yield patterned
    yield feats


def pattern_parameters(rows, grids, stem: str, where: str = "", circles=()) -> list[tuple[str, float, str]]:
    """``(name, value, note)`` before rounding.

    ``stem`` is ``hole``, ``base_hole``, or ``leg_hole``. ``where`` is
    ``foot`` or ``upright`` on a bracket, and empty on a plate.
    ``circles`` use ``stem_circle_N_r`` for the pitch-circle radius.
    """
    diameter = "diameter" if not where else f"diameter, through the {where}"
    out: list[tuple[str, float, str]] = []
    for i, row in enumerate(rows, 1):
        name = f"{stem}_row_{i}"
        out.extend([
            (f"{name}_n", float(row["n"]), "count of holes"),
            (f"{name}_pitch", float(row["pitch"]), "center spacing"),
            (f"{name}_d", float(row["d"]), diameter),
            (f"{name}_x", float(row["x"]), "start center"),
            (f"{name}_y", float(row["y"]), "start center"),
            (f"{name}_angle", float(row["angle"]), "direction from the start, degrees"),
        ])
    for i, grid in enumerate(grids, 1):
        name = f"{stem}_grid_{i}"
        out.extend([
            (f"{name}_nx", float(grid["nx"]), "count along X"),
            (f"{name}_ny", float(grid["ny"]), "count along Y"),
            (f"{name}_pitch_x", float(grid["pitch_x"]), "center spacing along X"),
            (f"{name}_pitch_y", float(grid["pitch_y"]), "center spacing along Y"),
            (f"{name}_d", float(grid["d"]), diameter),
            (f"{name}_x", float(grid["x"]), "corner center"),
            (f"{name}_y", float(grid["y"]), "corner center"),
            (f"{name}_angle", float(grid["angle"]), "direction of the X axis, degrees"),
        ])
    for i, circle in enumerate(circles, 1):
        name = f"{stem}_circle_{i}"
        out.extend([
            (f"{name}_n", float(circle["n"]), "count of holes"),
            (f"{name}_d", float(circle["d"]), diameter),
            (f"{name}_r", float(circle["r"]), "pitch circle radius"),
            (f"{name}_cx", float(circle["cx"]), "pitch circle center"),
            (f"{name}_cy", float(circle["cy"]), "pitch circle center"),
            (f"{name}_angle", float(circle["angle"]), "start angle from +X, degrees"),
        ])
    return out


def box_pattern_parameters(feats: dict) -> list[tuple[str, float, str]]:
    return pattern_parameters(
        feats.get("hole_rows") or [],
        feats.get("hole_grids") or [],
        "hole",
        circles=feats.get("hole_circles") or [],
    )


def member_pattern_parameters(bracket: dict, member: str) -> list[tuple[str, float, str]]:
    where = "foot" if member == "base" else "upright"
    return pattern_parameters(
        bracket.get(f"{member}_hole_rows") or [],
        bracket.get(f"{member}_hole_grids") or [],
        f"{member}_hole",
        where,
        circles=bracket.get(f"{member}_hole_circles") or [],
    )


def iter_box_patterns(feats: dict):
    """``(prefix, kind)`` with ``kind`` ``row``, ``grid``, or ``circle``."""
    for i, _row in enumerate(feats.get("hole_rows") or (), 1):
        yield f"hole_row_{i}", "row"
    for i, _grid in enumerate(feats.get("hole_grids") or (), 1):
        yield f"hole_grid_{i}", "grid"
    for i, _circle in enumerate(feats.get("hole_circles") or (), 1):
        yield f"hole_circle_{i}", "circle"


def iter_bracket_patterns(feats: dict):
    """Foot patterns first, then the upright. ``kind`` is ``row``, ``grid``, or ``circle``."""
    bracket = feats.get("bracket") or {}
    for member in ("base", "leg"):
        for i, _row in enumerate(bracket.get(f"{member}_hole_rows") or (), 1):
            yield f"{member}_hole_row_{i}", "row"
        for i, _grid in enumerate(bracket.get(f"{member}_hole_grids") or (), 1):
            yield f"{member}_hole_grid_{i}", "grid"
        for i, _circle in enumerate(bracket.get(f"{member}_hole_circles") or (), 1):
            yield f"{member}_hole_circle_{i}", "circle"


def python_pattern_cuts(prefix: str, kind: str, style: str) -> list[str]:
    """A Python loop that cuts one cylinder per hole. ``style`` is the backend."""
    if style == "cadquery":
        cut = (f"result = result.cut(cq.Workplane(\"XY\", origin=(_x, _y, base_z0 - 1))"
               f".circle({prefix}_d / 2).extrude(base_t + 2))")
    elif style == "freecad":
        cut = (f"result = result.cut(Part.makeCylinder({prefix}_d / 2, base_t + 2, "
               f"App.Vector(_x, _y, base_z0 - 1), App.Vector(0, 0, 1)))")
    else:
        raise ValueError(style)
    if kind == "row":
        x = f"{prefix}_x + _i * {prefix}_pitch * math.cos(math.radians({prefix}_angle))"
        y = f"{prefix}_y + _i * {prefix}_pitch * math.sin(math.radians({prefix}_angle))"
        return [
            f"for _i in range(int(round({prefix}_n))):",
            f"    _x = {x}",
            f"    _y = {y}",
            f"    {cut}",
        ]
    if kind == "circle":
        x = (f"{prefix}_cx + {prefix}_r * math.cos(math.radians("
             f"{prefix}_angle + _i * 360.0 / {prefix}_n))")
        y = (f"{prefix}_cy + {prefix}_r * math.sin(math.radians("
             f"{prefix}_angle + _i * 360.0 / {prefix}_n))")
        return [
            f"for _i in range(int(round({prefix}_n))):",
            f"    _x = {x}",
            f"    _y = {y}",
            f"    {cut}",
        ]
    x = (f"{prefix}_x + _i * {prefix}_pitch_x * math.cos(math.radians({prefix}_angle))"
         f" - _j * {prefix}_pitch_y * math.sin(math.radians({prefix}_angle))")
    y = (f"{prefix}_y + _i * {prefix}_pitch_x * math.sin(math.radians({prefix}_angle))"
         f" + _j * {prefix}_pitch_y * math.cos(math.radians({prefix}_angle))")
    return [
        f"for _j in range(int(round({prefix}_ny))):",
        f"    for _i in range(int(round({prefix}_nx))):",
        f"        _x = {x}",
        f"        _y = {y}",
        f"        {cut}",
    ]


def openscad_pattern_cut(prefix: str, kind: str) -> str:
    """One ``for`` loop of cylinders, subtracted with the other holes."""
    if kind == "row":
        return (
            f"for (_i = [0 : {prefix}_n - 1])\n"
            f"  translate([{prefix}_x + _i * {prefix}_pitch * cos({prefix}_angle),\n"
            f"             {prefix}_y + _i * {prefix}_pitch * sin({prefix}_angle),\n"
            f"             base_z0 - 1])\n"
            f"    cylinder(h=base_t + 2, d={prefix}_d);"
        )
    if kind == "circle":
        return (
            f"for (_i = [0 : {prefix}_n - 1])\n"
            f"  translate([{prefix}_cx + {prefix}_r * cos({prefix}_angle + _i * 360 / {prefix}_n),\n"
            f"             {prefix}_cy + {prefix}_r * sin({prefix}_angle + _i * 360 / {prefix}_n),\n"
            f"             base_z0 - 1])\n"
            f"    cylinder(h=base_t + 2, d={prefix}_d);"
        )
    return (
        f"for (_j = [0 : {prefix}_ny - 1])\n"
        f"  for (_i = [0 : {prefix}_nx - 1])\n"
        f"    translate([{prefix}_x + _i * {prefix}_pitch_x * cos({prefix}_angle)"
        f" - _j * {prefix}_pitch_y * sin({prefix}_angle),\n"
        f"               {prefix}_y + _i * {prefix}_pitch_x * sin({prefix}_angle)"
        f" + _j * {prefix}_pitch_y * cos({prefix}_angle),\n"
        f"               base_z0 - 1])\n"
        f"      cylinder(h=base_t + 2, d={prefix}_d);"
    )
