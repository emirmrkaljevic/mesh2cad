"""Emit an editable, parametric OpenSCAD script from a fitted model.

Uses the same feature dict and parameter names as the CadQuery export
(:mod:`mesh2cad.parametric`). The script is plain CSG: ``cube`` / ``cylinder``
(with ``r1``/``r2`` for a draft or a chamfer; ``rotate`` when the axis is X or
Y), ``difference`` / ``union``, and ``offset`` for the one vertical corner
radius. A straight extrusion is ``linear_extrude`` of ``polygon(outline)``.
Several of those, cut, intersected, or stacked as partial-depth extrusions,
are the same polygons under ``intersection`` / ``union`` / ``difference``,
and a counterbored hole is a second cylinder. A face tilted to X, Y, and Z
is one planar clip.
Customizer annotations (``/* [Group] */``, ``// [min:max]``) are added
only where the part has a real bound.

When the ``openscad`` program is on ``PATH``, the script is rendered to STL and
compared with the input mesh. Without that program the script is still
returned, with no check. OpenSCAD is not a package dependency.
"""
from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import trimesh

from .extrude import (profile_cbore_cutter, profile_cutter, profile_parameters,
                      profile_point_lists)
from .parametric import (Unsupported, axis_point_expr, boss_span, cbore_cut_span,
                         chamfer_is_square, deviation_problem, hole_cut_span,
                         named_parameters, openscad_axis_rotate, recognize)
from .report import deviation

# Chord error of a radius-8 circle at this count is about 0.004 mm, under the
# fidelity limits used for these parts. Hidden from the Customizer.
CIRCLE_SEGMENTS = 96


def openscad_bin() -> str | None:
    return shutil.which("openscad")


def _fmt_bound(v: float) -> str:
    r = round(float(v), 2)
    if abs(r - round(r)) < 1e-9:
        return str(int(round(r)))
    if abs(r - round(r, 1)) < 1e-9:
        return f"{r:.1f}"
    return f"{r:.2f}".rstrip("0").rstrip(".")


def _slider(lo: float, hi: float, step: float, value: float, floor_hi: bool = False) -> str:
    """Customizer ``[min:step:max]``, widened so ``value`` sits inside it.

    ``floor_hi`` keeps the max from stepping past ``hi``. A draft slider that
    rounds up past the tip-closing angle, or a fillet slider that rounds up to
    half the short side, asks OpenSCAD for an empty solid.
    """
    lo = min(float(lo), float(value))
    hi = max(float(hi), float(value))
    lo_s = math.floor((lo / step) + 1e-9) * step
    if floor_hi:
        hi_s = math.floor((hi / step) + 1e-9) * step
        if hi_s + 1e-9 < value:
            hi_s = value
    else:
        hi_s = math.ceil((hi / step) - 1e-9) * step
    if hi_s - lo_s < step / 2:
        return ""
    if abs(step - 1.0) < 1e-9:
        return f"[{_fmt_bound(lo_s)}:{_fmt_bound(hi_s)}]"
    return f"[{_fmt_bound(lo_s)}:{_fmt_bound(step)}:{_fmt_bound(hi_s)}]"


def _feature_axis(by_name: dict[str, float], kind: str, idx: str) -> str:
    """X, Y, or Z from which of ``{kind}_{idx}_{x,y,z}`` the script actually has."""
    has = {c: f"{kind}_{idx}_{c}" in by_name for c in "xyz"}
    if has["y"] and has["z"] and not has["x"]:
        return "X"
    if has["x"] and has["z"] and not has["y"]:
        return "Y"
    return "Z"


def _face_limit(by_name: dict[str, float], axis: str) -> float:
    """Shorter side of the face perpendicular to ``axis``."""
    if axis == "X":
        return min(by_name["base_y"], by_name["base_t"])
    if axis == "Y":
        return min(by_name["base_x"], by_name["base_t"])
    return min(by_name["base_x"], by_name["base_y"])


def _customizer_range(name: str, value: float, by_name: dict[str, float]) -> str:
    """A slider bound only when the geometry actually limits the number."""
    plate = min(by_name["base_x"], by_name["base_y"])
    thick = by_name["base_t"]
    if name == "fillet_r":
        # offset() of a zero-size square is an empty solid, so stop short of half the short side
        return _slider(0.0, plate / 2.0 - 0.1, 0.1, value, floor_hi=True)
    hole_d = _indexed(name, "hole_", "_d")
    if hole_d:
        return _slider(0.2, _face_limit(by_name, _feature_axis(by_name, "hole", hole_d)), 0.1, value)
    hole_depth = _indexed(name, "hole_", "_depth")
    if hole_depth:
        span = {"X": by_name["base_x"], "Y": by_name["base_y"], "Z": thick}
        return _slider(0.2, span[_feature_axis(by_name, "hole", hole_depth)], 0.1, value)
    csk_d = _indexed(name, "hole_", "_csk_d")
    if csk_d:
        return _slider(by_name[f"hole_{csk_d}_d"], plate, 0.1, value)
    csk_a = _indexed(name, "hole_", "_csk_angle")
    if csk_a and 30.0 <= value <= 150.0:
        lo, hi = (60.0, 120.0) if 60.0 <= value <= 120.0 else (30.0, 150.0)
        return _slider(lo, hi, 1.0, value)
    cbore_d = _indexed(name, "hole_", "_cbore_d")
    if cbore_d:
        return _slider(by_name.get(f"hole_{cbore_d}_d", 0.2),
                       _face_limit(by_name, _feature_axis(by_name, "hole", cbore_d)), 0.1, value)
    cbore_depth = _indexed(name, "hole_", "_cbore_depth")
    if cbore_depth:
        span = {"X": by_name["base_x"], "Y": by_name["base_y"], "Z": thick}
        return _slider(0.2, span[_feature_axis(by_name, "hole", cbore_depth)], 0.1, value)
    boss = _indexed(name, "boss_", "_d")
    if boss:
        return _slider(0.2, _face_limit(by_name, _feature_axis(by_name, "boss", boss)), 0.1, value)
    draft = _indexed(name, "boss_", "_draft_deg")
    if draft:
        radius = by_name[f"boss_{draft}_d"] / 2.0
        height = by_name[f"boss_{draft}_h"]
        cap = math.degrees(math.atan2(radius, height)) if height > 1e-9 else value
        return _slider(0.0, cap, 0.5, value, floor_hi=True)
    chamfer = _indexed(name, "boss_", "_chamfer")
    if chamfer:
        cap = min(by_name[f"boss_{chamfer}_h"], by_name[f"boss_{chamfer}_d"] / 2.0)
        return _slider(0.0, cap, 0.1, value)
    ch_h = _indexed(name, "boss_", "_chamfer_h")
    if ch_h:
        return _slider(0.0, by_name[f"boss_{ch_h}_h"], 0.1, value)
    ch_w = _indexed(name, "boss_", "_chamfer_w")
    if ch_w:
        return _slider(0.0, by_name[f"boss_{ch_w}_d"] / 2.0, 0.1, value)
    pad_w = _indexed(name, "pad_", "_w")
    if pad_w:
        return _slider(0.2, by_name["base_x"], 0.1, value)
    pad_l = _indexed(name, "pad_", "_l")
    if pad_l:
        return _slider(0.2, by_name["base_y"], 0.1, value)
    return ""


def _indexed(name: str, prefix: str, suffix: str) -> str | None:
    if name.startswith(prefix) and name.endswith(suffix):
        mid = name[len(prefix):-len(suffix)]
        if mid.isdigit():
            return mid
    return None


def _param_lines(params: list[tuple[str, float, str]]) -> list[str]:
    by_name = {n: v for n, v, _ in params}
    width = max(len(n) for n, _, _ in params)

    def group(name: str) -> str:
        if name.startswith("pad_"):
            return "Pads"
        if name.startswith("boss_"):
            return "Bosses"
        if name.startswith("hole_"):
            return "Holes"
        return "Base"

    lines = ["// --- parameters ---"]
    current = None
    for name, value, note in params:
        g = group(name)
        if g != current:
            if current is not None:
                lines.append("")
            lines.append(f"/* [{g}] */")
            current = g
        if note:
            lines.append(f"// {note}")
        line = f"{name.ljust(width)} = {value!r};"
        slider = _customizer_range(name, value, by_name)
        if slider:
            line += f" // {slider}"
        lines.append(line)
    lines += ["", "/* [Hidden] */",
              f"$fn = {CIRCLE_SEGMENTS}; // circle segments"]
    return lines


def _needs_top(feats: dict) -> bool:
    if feats["pads"] or any(b.get("axis", "Z") == "Z" for b in feats["bosses"]):
        return True
    return any(h.get("axis", "Z") == "Z" and (
        "depth" in h or "csk_d" in h or ("cbore_d" in h and h.get("from", 1) > 0)
    ) for h in feats["holes"])


def _side_cylinder(prefix: str, feat: dict, diam: str, start: str, length: str) -> str:
    point = axis_point_expr(feat["axis"], prefix, start)
    return (f"translate([{point}])\n"
            f"  rotate({openscad_axis_rotate(feat['axis'])})\n"
            f"    cylinder(h={length}, d={diam});")


def _solids(feats: dict, tol: float) -> list[str]:
    if "fillet_r" in feats:
        base = (
            "translate([base_cx, base_cy, base_z0])\n"
            "  linear_extrude(base_t)\n"
            "    offset(r=fillet_r) // four vertical corners, one radius\n"
            "      square([base_x - 2 * fillet_r, base_y - 2 * fillet_r], center=true);"
        )
    else:
        base = (
            "translate([base_cx, base_cy, base_z0 + base_t / 2])\n"
            "  cube([base_x, base_y, base_t], center=true);"
        )
    solids = [base]
    for i, _d in enumerate(feats["pads"], 1):
        n = f"pad_{i}"
        solids.append(
            f"translate([{n}_x, {n}_y, top + {n}_h / 2])\n"
            f"  cube([{n}_w, {n}_l, {n}_h], center=true);"
        )
    for i, d in enumerate(feats["bosses"], 1):
        n = f"boss_{i}"
        if d.get("axis", "Z") != "Z":
            start, length = boss_span(n, d)
            solids.append(_side_cylinder(n, d, f"{n}_d", start, length))
            continue
        if "draft_deg" in d:
            solids.append(
                f"translate([{n}_x, {n}_y, top])\n"
                f"  cylinder(h={n}_h, r1={n}_d / 2, r2={n}_d / 2 - {n}_h * tan({n}_draft_deg));"
            )
        elif "chamfer" in d:
            if chamfer_is_square(d, tol):
                ch = cw = f"{n}_chamfer"
            else:
                ch, cw = f"{n}_chamfer_h", f"{n}_chamfer_w"
            solids.append(
                f"translate([{n}_x, {n}_y, top]) {{\n"
                f"  cylinder(h={n}_h - {ch}, d={n}_d);\n"
                f"  translate([0, 0, {n}_h - {ch}])\n"
                f"    cylinder(h={ch}, r1={n}_d / 2, r2={n}_d / 2 - {cw});\n"
                f"}}"
            )
        else:
            solids.append(
                f"translate([{n}_x, {n}_y, top])\n"
                f"  cylinder(h={n}_h, d={n}_d);"
            )
    return solids


def _cuts(feats: dict) -> list[str]:
    cuts = []
    for i, d in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        if d.get("axis", "Z") != "Z":
            start, length = hole_cut_span(n, d)
            cuts.append(_side_cylinder(n, d, f"{n}_d", start, length))
            if "cbore_d" in d:
                c0, clen = cbore_cut_span(n, d)
                cuts.append(_side_cylinder(n, d, f"{n}_cbore_d", c0, clen))
            continue
        if "depth" in d:
            cuts.append(
                f"translate([{n}_x, {n}_y, top - {n}_depth])\n"
                f"  cylinder(h={n}_depth + 1, d={n}_d);"
            )
        else:
            cuts.append(
                f"translate([{n}_x, {n}_y, base_z0 - 1])\n"
                f"  cylinder(h=base_t + 2, d={n}_d);"
            )
        if "cbore_d" in d:
            if d.get("from", 1) > 0:
                cuts.append(
                    f"translate([{n}_x, {n}_y, top - {n}_cbore_depth])\n"
                    f"  cylinder(h={n}_cbore_depth + 1, d={n}_cbore_d);"
                )
            else:
                cuts.append(
                    f"translate([{n}_x, {n}_y, base_z0 - 1])\n"
                    f"  cylinder(h={n}_cbore_depth + 1, d={n}_cbore_d);"
                )
        if "csk_d" in d:
            cuts.append(
                f"translate([{n}_x, {n}_y, top - {n}_csk_h])\n"
                f"  cylinder(h={n}_csk_h + 1, r1={n}_d / 2,\n"
                f"           r2={n}_csk_d / 2 + tan({n}_csk_angle / 2));"
            )
    return cuts


def _indent(text: str, levels: int) -> str:
    pad = "  " * levels
    return "\n".join(pad + ln if ln else ln for ln in text.splitlines())


def _model(feats: dict, tol: float) -> list[str]:
    solids = _solids(feats, tol)
    cuts = _cuts(feats)
    lines = ["// --- model ---", "module mesh2cad_part() {"]
    if _needs_top(feats):
        lines.append("  top = base_z0 + base_t;")
    for i, d in enumerate(feats["holes"], 1):
        if "csk_d" in d:
            n = f"hole_{i}"
            lines.append(
                f"  {n}_csk_h = ({n}_csk_d - {n}_d) / 2 / tan({n}_csk_angle / 2);"
            )
    body_at = 2 if cuts else 1
    if cuts:
        lines.append("  difference() {")
    if len(solids) == 1:
        lines.append(_indent(solids[0], body_at))
    else:
        lines.append(f"{'  ' * body_at}union() {{")
        for s in solids:
            lines.append(_indent(s, body_at + 1))
        lines.append(f"{'  ' * body_at}}}")
    if cuts:
        for c in cuts:
            lines.append(_indent(c, 2))
        lines.append("  }")
    lines += ["}", "", "mesh2cad_part();", ""]
    return lines


def _fmt_points_scad(pts) -> str:
    rows = ["["]
    for x, y in pts:
        rows.append(f"  [{x!r}, {y!r}],")
    rows.append("]")
    return "\n".join(rows)


def _profile_param_lines(feats: dict, tol: float) -> list[str]:
    numbers = profile_parameters(feats, tol)
    lists = profile_point_lists(feats)
    width = max([len(n) for n, _, _ in numbers] + [len(n) for n, _, _ in lists])
    # hole diameter slider: up to the shorter side of the outline's bounding box
    outline = feats["profile"]["outline"]
    span_u = max(p[0] for p in outline) - min(p[0] for p in outline)
    span_v = max(p[1] for p in outline) - min(p[1] for p in outline)
    face = min(span_u, span_v, feats["profile"]["length"])
    lines = ["// --- parameters ---", "/* [Outline] */"]
    for name, value, note in numbers:
        if name.startswith("hole_") or name.startswith("boss_"):
            continue
        if note:
            lines.append(f"// {note}")
        lines.append(f"{name.ljust(width)} = {value!r};")
    for name, pts, note in lists:
        if note:
            lines.append(f"// {note}")
        lines.append(f"{name} = {_fmt_points_scad(pts)};")
    hole_rows = [(n, v, note) for n, v, note in numbers if n.startswith("hole_")]
    if hole_rows:
        lines += ["", "/* [Holes] */"]
        for name, value, note in hole_rows:
            if note:
                lines.append(f"// {note}")
            line = f"{name.ljust(width)} = {value!r};"
            if name.endswith("_cbore_d"):
                stem = name[: -len("_cbore_d")]
                slider = _slider(by_name.get(f"{stem}_d", 0.2), max(face, value), 0.1, value)
            elif name.endswith("_cbore_depth"):
                slider = _slider(0.2, feats["profile"]["length"], 0.1, value)
            elif name.endswith("_d"):
                slider = _slider(0.2, face, 0.1, value)
            else:
                slider = ""
            if slider:
                line += f" // {slider}"
            lines.append(line)
    lines += ["", "/* [Hidden] */",
              f"$fn = {CIRCLE_SEGMENTS}; // circle segments"]
    return lines


def _profile_model(feats: dict) -> list[str]:
    pr = feats["profile"]
    axis = pr["axis"]
    cut_lines = []
    if pr["cutouts"]:
        cut_lines.append("difference() {")
        cut_lines.append("  polygon(outline);")
        for i in range(1, len(pr["cutouts"]) + 1):
            cut_lines.append(f"  polygon(cutout_{i});")
        cut_lines.append("}")
    else:
        cut_lines.append("polygon(outline);")
    profile_2d = "\n".join(cut_lines)
    if axis == "Z":
        placed = "translate([0, 0, extrude_0])\n  linear_extrude(extrude_len)\n" + _indent(profile_2d, 2)
    elif axis == "Y":
        # outline is (x, z). rotate 90° about X turns +Y of the polygon into +Z
        # and the extrusion (local +Z) into -Y, starting from the far end.
        placed = ("translate([0, extrude_0 + extrude_len, 0])\n"
                  "  rotate([90, 0, 0])\n"
                  "    linear_extrude(extrude_len)\n" + _indent(profile_2d, 3))
    else:
        # outline is (y, z), extruded along +X
        placed = ("translate([extrude_0, 0, 0])\n"
                  "  rotate([90, 0, 90])\n"
                  "    linear_extrude(extrude_len)\n" + _indent(profile_2d, 3))
    holes = []
    for i, hole in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        start, length = profile_cutter(hole, "openscad")
        holes.append(_side_cylinder(n, hole, f"{n}_d", start, length))
        if "cbore_d" in hole:
            c0, clen = profile_cbore_cutter(hole, n, "openscad")
            holes.append(_side_cylinder(n, hole, f"{n}_cbore_d", c0, clen))
    lines = ["// --- model ---", "module mesh2cad_part() {"]
    if holes:
        lines.append("  difference() {")
        lines.append(_indent(placed, 2))
        for h in holes:
            lines.append(_indent(h, 2))
        lines.append("  }")
    else:
        lines.append(_indent(placed, 1))
    lines += ["}", "", "mesh2cad_part();", ""]
    return lines


def _render_profile(feats: dict, tol: float, source: str | None) -> str:
    head = ["// Parametric OpenSCAD model generated by mesh2cad.",
            "//",
            f"// Source: {source}" if source else "// Edit the parameters below and re-render.",
            "// Units are the mesh units (usually mm). Render (F6), then export STL.",
            "//"]
    return "\n".join(head + _profile_param_lines(feats, tol) + [""] + _profile_model(feats))


def render(feats: dict, tol: float = 0.01, source: str | None = None) -> str:
    """Turn a feature dict into an OpenSCAD script with named parameters."""
    if "csg" in feats:
        from .csg import render_openscad
        return render_openscad(feats, tol, source)
    if "profile" in feats:
        return _render_profile(feats, tol, source)
    head = ["// Parametric OpenSCAD model generated by mesh2cad.",
            "//",
            f"// Source: {source}" if source else "// Edit the parameters below and re-render.",
            "// Units are the mesh units (usually mm). Render (F6), then export STL.",
            "//"]
    return "\n".join(head + _param_lines(named_parameters(feats, tol)) + [""] + _model(feats, tol))


def render_mesh(code: str, work: Path, timeout: float = 120.0) -> trimesh.Trimesh:
    """Render ``code`` with the OpenSCAD CLI and load the STL."""
    bin = openscad_bin()
    if not bin:
        raise FileNotFoundError("openscad is not on PATH")
    work.mkdir(parents=True, exist_ok=True)
    scad = work / "part.scad"
    stl = work / "part.stl"
    scad.write_text(code)
    proc = subprocess.run([bin, "-o", str(stl), str(scad)],
                          capture_output=True, text=True, timeout=timeout)
    err = (proc.stderr or proc.stdout or "").strip()
    if proc.returncode != 0 or not stl.is_file() or stl.stat().st_size == 0:
        tail = err.splitlines()[-1] if err else f"exit {proc.returncode}"
        raise RuntimeError(tail)
    out = trimesh.load(str(stl), force="mesh", process=True)
    if not isinstance(out, trimesh.Trimesh) or len(out.faces) == 0:
        raise RuntimeError("OpenSCAD wrote an empty mesh")
    if not out.is_watertight:
        raise RuntimeError("OpenSCAD mesh is not a closed solid")
    return out


def check(code: str, m: trimesh.Trimesh, work: Path, n: int = 4000) -> dict:
    """Render the script and measure how well it reproduces mesh ``m``."""
    out = render_mesh(code, work)
    dev = deviation(m, out, n=n)
    vm = abs(float(m.volume))
    vo = abs(float(out.volume))
    return {"max_deviation": float(dev["hausdorff"]),
            "hausdorff": float(dev["hausdorff"]),
            "rms": float(dev["rms"]),
            "volume_error_pct": float(100.0 * abs(vo - vm) / vm) if vm else float("inf")}


def to_openscad(m: trimesh.Trimesh, regions, tol: float, source: str | None = None,
                verify: bool = True) -> tuple[str, dict | None]:
    """Generate the script. When ``verify`` and ``openscad`` is installed, render
    it and raise :class:`Unsupported` if it does not reproduce the mesh.

    Returns ``(code, None)`` when OpenSCAD is not available: the script is
    still produced, and the caller reports it as unverified.
    """
    code = render(recognize(m, regions, tol), tol=tol, source=source)
    if not verify or openscad_bin() is None:
        return code, None
    try:
        with tempfile.TemporaryDirectory(prefix="mesh2cad-scad-") as d:
            q = check(code, m, Path(d))
    except Exception as e:  # noqa: BLE001
        raise Unsupported(f"generated script failed to render: {e}") from e
    problem = deviation_problem(q, tol)
    if problem:
        raise Unsupported(problem)
    return code, q
