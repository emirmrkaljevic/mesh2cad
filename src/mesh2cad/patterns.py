"""Regular patterns of plain Z through-holes for the parametric exporters.

A row is three or more holes of one diameter whose centers lie on a straight
line and on one pitch. Two holes are left as ``hole_N_*`` (or ``base_hole_N_*``
/ ``leg_hole_N_*``). They have a single gap, so equal spacing is automatic,
and moving one of them should not drag the other.

A rectangular grid is an axis-aligned block, at least two by two, same
diameter, constant pitch along X and along Y. A bolt circle is not detected.
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
    """How far a center may sit from the ideal row or grid, in millimetres."""
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


def partition_holes(holes: list[dict], tol: float):
    """``(kept, rows, grids)``. ``kept`` preserves the original order.

    Grids are taken first, so a 2 by 3 block is one grid and not two rows.
    """
    consumed: set[int] = set()
    rows: list[dict] = []
    grids: list[dict] = []
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
    kept = [h for h in holes if id(h) not in consumed]
    rows.sort(key=lambda r: (r["y"], r["x"], r["angle"], r["d"]))
    grids.sort(key=lambda g: (g["y"], g["x"], g["d"]))
    return kept, rows, grids


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
        kept, rows, grids = partition_holes(list(feats["holes"]), tol)
        if rows or grids:
            changed = True
            copied = copy()
            copied["holes"] = kept
            if rows:
                copied["hole_rows"] = rows
            if grids:
                copied["hole_grids"] = grids
    bracket = feats.get("bracket")
    if bracket:
        rewritten = None
        for member in ("base", "leg"):
            key = f"{member}_holes"
            if not bracket.get(key):
                continue
            kept, rows, grids = partition_holes(list(bracket[key]), tol)
            if not rows and not grids:
                continue
            changed = True
            if rewritten is None:
                rewritten = dict(bracket)
            rewritten[key] = kept
            if rows:
                rewritten[f"{member}_hole_rows"] = rows
            if grids:
                rewritten[f"{member}_hole_grids"] = grids
        if rewritten is not None:
            copy()["bracket"] = rewritten
    if not changed:
        return feats
    return out


def has_hole_patterns(feats: dict) -> bool:
    if feats.get("hole_rows") or feats.get("hole_grids"):
        return True
    bracket = feats.get("bracket") or {}
    return any(bracket.get(name) for name in (
        "base_hole_rows", "base_hole_grids", "leg_hole_rows", "leg_hole_grids",
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


def pattern_parameters(rows, grids, stem: str, where: str = "") -> list[tuple[str, float, str]]:
    """``(name, value, note)`` before rounding.

    ``stem`` is ``hole``, ``base_hole``, or ``leg_hole``. ``where`` is
    ``foot`` or ``upright`` on a bracket, and empty on a plate.
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
    return out


def box_pattern_parameters(feats: dict) -> list[tuple[str, float, str]]:
    return pattern_parameters(feats.get("hole_rows") or [], feats.get("hole_grids") or [], "hole")


def member_pattern_parameters(bracket: dict, member: str) -> list[tuple[str, float, str]]:
    where = "foot" if member == "base" else "upright"
    return pattern_parameters(
        bracket.get(f"{member}_hole_rows") or [],
        bracket.get(f"{member}_hole_grids") or [],
        f"{member}_hole",
        where,
    )


def iter_box_patterns(feats: dict):
    """``(prefix, kind)`` with ``kind`` ``row`` or ``grid``."""
    for i, _row in enumerate(feats.get("hole_rows") or (), 1):
        yield f"hole_row_{i}", "row"
    for i, _grid in enumerate(feats.get("hole_grids") or (), 1):
        yield f"hole_grid_{i}", "grid"


def iter_bracket_patterns(feats: dict):
    """Foot patterns first, then the upright. ``kind`` is ``row`` or ``grid``."""
    bracket = feats.get("bracket") or {}
    for member in ("base", "leg"):
        for i, _row in enumerate(bracket.get(f"{member}_hole_rows") or (), 1):
            yield f"{member}_hole_row_{i}", "row"
        for i, _grid in enumerate(bracket.get(f"{member}_hole_grids") or (), 1):
            yield f"{member}_hole_grid_{i}", "grid"


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
