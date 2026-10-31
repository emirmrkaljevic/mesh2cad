"""Emit an editable, parametric FreeCAD macro from a fitted model.

Uses the same feature dict and parameter names as the CadQuery and OpenSCAD
exports (:mod:`mesh2cad.parametric`). The macro is plain Python for FreeCAD's
Part API: ``makeBox`` / ``makeCylinder`` / ``makeCone``, a polygon extruded
along X, Y, or Z, ``fuse`` / ``cut`` / ``common``, and ``makeFillet`` for the
one vertical corner radius. A part that is a few extrusions (their
intersection, a union of partial-depth extrusions, or one minus cuts along
the other axes, plus a counterbored hole or a planar clip) is the same
boolean. Every dimension is a named parameter at the top of the file
(``hole_1_d = 6.0``).

When ``freecadcmd`` is on ``PATH``, the macro is executed and the solid is
compared with the input mesh. Without that program the macro is still
returned, with no check. FreeCAD is not a package dependency.
"""
from __future__ import annotations

import os
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path

import trimesh

from .csg import _OVERLAP, _fmt_py, _grown_expr, _hole_limits, _names, _param_rows
from .extrude import (profile_cbore_cutter, profile_cutter, profile_parameters,
                      profile_point_lists)
from .parametric import (Unsupported, axis_point_expr, boss_span, cbore_cut_span,
                         chamfer_is_square, deviation_problem, hole_cut_span,
                         named_parameters, recognize)
from .report import deviation

__all__ = ["Unsupported", "check", "freecad_bin", "render", "render_mesh", "to_freecad"]

_DIR = {
    "X": "App.Vector(1, 0, 0)",
    "Y": "App.Vector(0, 1, 0)",
    "Z": "App.Vector(0, 0, 1)",
}

# Chord error stays under the fidelity limits used for these parts.
_TESSELLATE = 0.01

_EXTRUDE_FN = '''
def _extrude(pts, axis, origin, length):
    """Extrude a closed outline along +X, +Y, or +Z.

    ``pts`` are (u, v) in the plane perpendicular to ``axis``:
    (y, z) for X, (x, z) for Y, (x, y) for Z.
    """
    pts = [(float(u), float(v)) for u, v in pts]
    area = 0.0
    for i, (u, v) in enumerate(pts):
        u2, v2 = pts[(i + 1) % len(pts)]
        area += u * v2 - u2 * v
    if area < 0.0:
        pts = list(reversed(pts))
    vecs = []
    for u, v in pts:
        if axis == "Z":
            vecs.append(App.Vector(u, v, origin))
        elif axis == "Y":
            vecs.append(App.Vector(u, origin, v))
        else:
            vecs.append(App.Vector(origin, u, v))
    vecs.append(vecs[0])
    if axis == "Z":
        direction = App.Vector(0, 0, length)
    elif axis == "Y":
        direction = App.Vector(0, length, 0)
    else:
        direction = App.Vector(length, 0, 0)
    return Part.Face(Part.makePolygon(vecs)).extrude(direction)
'''.strip("\n")

_CLIP_FN = '''
def _clip(nx, ny, nz, d):
    """Solid on the outward side of the plane nx*x + ny*y + nz*z = d."""
    nn = nx * nx + ny * ny + nz * nz
    origin = App.Vector(nx * d / nn, ny * d / nn, nz * d / nn)
    raw = App.Vector(0, -nz, ny) if abs(nx) < 0.9 else App.Vector(-nz, 0, nx)
    scale = raw.Length
    xdir = App.Vector(raw.x / scale, raw.y / scale, raw.z / scale)
    zscale = math.sqrt(nn)
    zdir = App.Vector(nx / zscale, ny / zscale, nz / zscale)
    yraw = zdir.cross(xdir)
    yscale = yraw.Length
    ydir = App.Vector(yraw.x / yscale, yraw.y / yscale, yraw.z / yscale)
    loop = []
    for su, sv in ((-200.0, -200.0), (200.0, -200.0), (200.0, 200.0), (-200.0, 200.0)):
        loop.append(App.Vector(
            origin.x + xdir.x * su + ydir.x * sv,
            origin.y + xdir.y * su + ydir.y * sv,
            origin.z + xdir.z * su + ydir.z * sv,
        ))
    loop.append(loop[0])
    return Part.Face(Part.makePolygon(loop)).extrude(
        App.Vector(zdir.x * 200.0, zdir.y * 200.0, zdir.z * 200.0))
'''.strip("\n")

_CHAMFER_FN = '''
def _edge_chamfer(shape, outline, z0, dist, edges):
    """Cut 45 degree wedges off the listed edges of a Z extrusion's start face."""
    count = len(outline)
    for index in edges:
        a = outline[index % count]
        b = outline[(index + 1) % count]
        dx = b[0] - a[0]
        dy = b[1] - a[1]
        length = math.hypot(dx, dy)
        if length < 1e-4:
            continue
        ox, oy = dy / length, -dx / length
        xdir = App.Vector(-ox, -oy, 0.0)
        zdir = App.Vector(dx / length, dy / length, 0.0)
        origin = App.Vector(a[0], a[1], z0)
        p0 = origin
        p1 = App.Vector(origin.x + xdir.x * dist, origin.y + xdir.y * dist, origin.z)
        p2 = App.Vector(origin.x, origin.y, origin.z + dist)
        wedge = Part.Face(Part.makePolygon([p0, p1, p2, p0])).extrude(
            App.Vector(zdir.x * length, zdir.y * length, zdir.z * length))
        shape = shape.cut(wedge)
    return shape
'''.strip("\n")

_DRIVER_BODY = r'''
import math
import struct
import sys
import traceback

def _xyz(p):
    if hasattr(p, "x"):
        return float(p.x), float(p.y), float(p.z)
    return float(p[0]), float(p[1]), float(p[2])

def _write_stl(path, pts, tris):
    xyz = [_xyz(p) for p in pts]
    with open(path, "wb") as handle:
        handle.write(b"mesh2cad".ljust(80, b"\0"))
        handle.write(struct.pack("<I", len(tris)))
        for tri in tris:
            a = xyz[tri[0]]
            b = xyz[tri[1]]
            c = xyz[tri[2]]
            ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
            vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
            nx = uy * vz - uz * vy
            ny = uz * vx - ux * vz
            nz = ux * vy - uy * vx
            length = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
            handle.write(struct.pack(
                "<12fH",
                nx / length, ny / length, nz / length,
                a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2],
                0,
            ))

def _main():
    ns = {"__name__": "mesh2cad_freecad", "__file__": MACRO}
    with open(MACRO, encoding="utf-8") as handle:
        code = handle.read()
    exec(compile(code, MACRO, "exec"), ns)
    shape = ns["result"]
    if shape is None or shape.isNull():
        raise RuntimeError("FreeCAD shape is null")
    solids = shape.Solids
    if len(solids) != 1:
        raise RuntimeError("expected one solid, got %d" % len(solids))
    pts, tris = solids[0].tessellate(TESSELLATE)
    if not tris:
        raise RuntimeError("tessellation is empty")
    _write_stl(STL, pts, tris)
    with open(MARKER, "w", encoding="utf-8") as handle:
        handle.write("ok\n")

try:
    _main()
except Exception:
    traceback.print_exc()
    sys.exit(1)
'''


def freecad_bin() -> str | None:
    for name in ("freecadcmd", "freecad-cmd", "FreeCADCmd"):
        found = shutil.which(name)
        if found:
            return found
    return None


def _head(source: str | None) -> list[str]:
    return [
        '"""Parametric FreeCAD model generated by mesh2cad.',
        "",
        f"Source: {source}" if source else "Edit the parameters below and run the macro again.",
        "Units are the mesh units (usually mm).",
        "In FreeCAD: Macro -> Macros... -> Execute.",
        "From a shell: freecadcmd this_file.FCMacro",
        '"""',
        "import math",
        "",
        "import FreeCAD as App",
        "import Part",
        "",
        "# --- parameters ---",
    ]


def _tail() -> list[str]:
    return [
        "",
        'if __name__ == "__main__":',
        '    _doc = App.newDocument("mesh2cad")',
        '    _obj = _doc.addObject("Part::Feature", "Part")',
        "    _obj.Shape = result",
        "    _doc.recompute()",
        "    try:",
        "        from pathlib import Path",
        '        _step = str(Path(__file__).with_suffix(".step"))',
        "    except NameError:",
        "        _step = None",
        "    if _step:",
        "        result.exportStep(_step)",
        '        print("wrote", _step)',
        "",
    ]


def _number_lines(numbers, lists) -> list[str]:
    names = [n for n, _, _ in numbers] + [n for n, _, _ in lists]
    width = max([len(n) for n in names] or [1])
    lines = []
    for name, value, note in numbers:
        line = f"{name.ljust(width)} = {value!r}"
        if note:
            line += f"  # {note}"
        lines.append(line)
    for name, pts, note in lists:
        if note:
            lines.append(f"# {note}")
        lines.append(f"{name} = {_fmt_py(pts)}")
    return lines


def _cylinder(prefix: str, feat: dict, radius: str, start: str, length: str) -> str:
    point = axis_point_expr(feat["axis"], prefix, start)
    return (f"Part.makeCylinder({radius}, {length}, "
            f"App.Vector({point}), {_DIR[feat['axis']]})")


def _cone(prefix: str, feat: dict, r1: str, r2: str, start: str, length: str) -> str:
    point = axis_point_expr(feat["axis"], prefix, start)
    return (f"Part.makeCone({r1}, {r2}, {length}, "
            f"App.Vector({point}), {_DIR[feat['axis']]})")


def _finish(head, params, helpers, body) -> str:
    chunks = head + params + ["", "# --- model ---"]
    if helpers:
        chunks.append("\n\n".join(helpers))
        chunks.append("")
    chunks.extend(body)
    chunks.extend(_tail())
    return "\n".join(chunks)


def _fillet_lines() -> list[str]:
    return [
        "_fillet_edges = []",
        "for _edge in result.Edges:",
        "    if len(_edge.Vertexes) != 2:",
        "        continue",
        "    _p0 = _edge.Vertexes[0].Point",
        "    _p1 = _edge.Vertexes[1].Point",
        "    if (abs(_p0.x - _p1.x) < 1e-7 and abs(_p0.y - _p1.y) < 1e-7 "
        "and abs(_p0.z - _p1.z) > 1e-7):",
        "        _fillet_edges.append(_edge)",
        "result = result.makeFillet(fillet_r, _fillet_edges)",
    ]


def _render_box(feats: dict, tol: float, source: str | None) -> str:
    body = [
        "result = Part.makeBox(base_x, base_y, base_t, "
        "App.Vector(base_cx - base_x / 2, base_cy - base_y / 2, base_z0))",
        "top = base_z0 + base_t",
    ]
    if "fillet_r" in feats:
        body.extend(_fillet_lines())
    for i, _pad in enumerate(feats["pads"], 1):
        n = f"pad_{i}"
        body.append(
            f"result = result.fuse(Part.makeBox({n}_w, {n}_l, {n}_h, "
            f"App.Vector({n}_x - {n}_w / 2, {n}_y - {n}_l / 2, top)))"
        )
    for i, d in enumerate(feats["bosses"], 1):
        n = f"boss_{i}"
        if d.get("axis", "Z") != "Z":
            start, length = boss_span(n, d)
            body.append(f"result = result.fuse({_cylinder(n, d, f'{n}_d / 2', start, length)})")
            continue
        if "draft_deg" in d:
            body.append(
                f"result = result.fuse(Part.makeCone({n}_d / 2, "
                f"{n}_d / 2 - {n}_h * math.tan(math.radians({n}_draft_deg)), "
                f"{n}_h, App.Vector({n}_x, {n}_y, top), App.Vector(0, 0, 1)))"
            )
        elif "chamfer" in d:
            if chamfer_is_square(d, tol):
                ch = cw = f"{n}_chamfer"
            else:
                ch, cw = f"{n}_chamfer_h", f"{n}_chamfer_w"
            body.append(
                f"result = result.fuse(Part.makeCylinder({n}_d / 2, {n}_h - {ch}, "
                f"App.Vector({n}_x, {n}_y, top), App.Vector(0, 0, 1)))"
            )
            body.append(
                f"result = result.fuse(Part.makeCone({n}_d / 2, {n}_d / 2 - {cw}, {ch}, "
                f"App.Vector({n}_x, {n}_y, top + {n}_h - {ch}), App.Vector(0, 0, 1)))"
            )
        else:
            body.append(
                f"result = result.fuse(Part.makeCylinder({n}_d / 2, {n}_h, "
                f"App.Vector({n}_x, {n}_y, top), App.Vector(0, 0, 1)))"
            )
    for i, d in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        if d.get("axis", "Z") != "Z":
            start, length = hole_cut_span(n, d)
            body.append(f"result = result.cut({_cylinder(n, d, f'{n}_d / 2', start, length)})")
            if "cbore_d" in d:
                c0, clen = cbore_cut_span(n, d)
                body.append(
                    f"result = result.cut({_cylinder(n, d, f'{n}_cbore_d / 2', c0, clen)})"
                )
            continue
        if "depth" in d:
            body.append(
                f"result = result.cut(Part.makeCylinder({n}_d / 2, {n}_depth + 1, "
                f"App.Vector({n}_x, {n}_y, top - {n}_depth), App.Vector(0, 0, 1)))"
            )
        else:
            body.append(
                f"result = result.cut(Part.makeCylinder({n}_d / 2, base_t + 2, "
                f"App.Vector({n}_x, {n}_y, base_z0 - 1), App.Vector(0, 0, 1)))"
            )
        if "csk_d" in d:
            body.append(
                f"_h = ({n}_csk_d - {n}_d) / 2 / math.tan(math.radians({n}_csk_angle / 2))"
            )
            body.append(
                f"result = result.cut(Part.makeCone({n}_d / 2, "
                f"{n}_csk_d / 2 + math.tan(math.radians({n}_csk_angle / 2)), "
                f"_h + 1, App.Vector({n}_x, {n}_y, top - _h), App.Vector(0, 0, 1)))"
            )
        if "cbore_d" in d:
            if d.get("from", 1) > 0:
                body.append(
                    f"result = result.cut(Part.makeCylinder({n}_cbore_d / 2, {n}_cbore_depth + 1, "
                    f"App.Vector({n}_x, {n}_y, top - {n}_cbore_depth), App.Vector(0, 0, 1)))"
                )
            else:
                body.append(
                    f"result = result.cut(Part.makeCylinder({n}_cbore_d / 2, {n}_cbore_depth + 1, "
                    f"App.Vector({n}_x, {n}_y, base_z0 - 1), App.Vector(0, 0, 1)))"
                )
    params = _number_lines(named_parameters(feats, tol), [])
    return _finish(_head(source), params, [], body)


def _extrude_lines(var, axis, origin, length, outline, cutouts, cutout_prefix) -> list[str]:
    lines = [f'{var} = _extrude({outline}, "{axis}", {origin}, {length})']
    for i in range(1, len(cutouts) + 1):
        lines.append(
            f'{var} = {var}.cut(_extrude({cutout_prefix}_{i}, "{axis}", {origin}, {length}))'
        )
    return lines


def _span(origin: str, length: str, overlap: bool) -> tuple[str, str]:
    if not overlap:
        return origin, length
    return f"({origin} - {_OVERLAP})", f"({length} + {2 * _OVERLAP})"


def _hole_lines(n: str, hole: dict, lo: str, hi: str) -> list[str]:
    """Cuts for one hole. Placement matches the CadQuery exporter."""
    lines = []
    if "depth" in hole:
        if hole.get("from", 1) < 0:
            start, length = f"({lo} - 1)", f"{n}_depth + 1"
        else:
            start, length = f"(({hi}) - {n}_depth)", f"{n}_depth + 1"
    else:
        start, length = f"({lo} - 1)", f"(({hi}) - ({lo}) + 2)"
    lines.append(f"result = result.cut({_cylinder(n, hole, f'{n}_d / 2', start, length)})")
    if "cbore_d" in hole:
        if hole.get("from", 1) < 0:
            c0, clen = f"({lo} - 1)", f"{n}_cbore_depth + 1"
        else:
            c0, clen = f"(({hi}) - {n}_cbore_depth)", f"{n}_cbore_depth + 1"
        lines.append(f"result = result.cut({_cylinder(n, hole, f'{n}_cbore_d / 2', c0, clen)})")
    if "chamfer" not in hole:
        return lines
    if hole.get("chamfer_in") and "depth" in hole:
        if hole.get("from", 1) < 0:
            c0 = f"(({lo}) + {n}_depth)"
            r1, r2 = f"{n}_d / 2", f"{n}_d / 2 - {n}_chamfer"
        else:
            c0 = f"(({hi}) - {n}_depth - {n}_chamfer)"
            r1, r2 = f"{n}_d / 2 - {n}_chamfer", f"{n}_d / 2"
    elif hole.get("chamfer_from", 1) < 0:
        c0 = f"({lo})"
        r1, r2 = f"{n}_d / 2 + {n}_chamfer", f"{n}_d / 2"
    else:
        c0 = f"(({hi}) - {n}_chamfer)"
        r1, r2 = f"{n}_d / 2", f"{n}_d / 2 + {n}_chamfer"
    lines.append(f"result = result.cut({_cone(n, hole, r1, r2, c0, f'{n}_chamfer')})")
    return lines


def _render_profile(feats: dict, tol: float, source: str | None) -> str:
    pr = feats["profile"]
    numbers = profile_parameters(feats, tol)
    lists = profile_point_lists(feats)
    origin, length = "extrude_0", "extrude_len"
    body = _extrude_lines("result", pr["axis"], origin, length, "outline",
                          pr["cutouts"], "cutout")
    for i, hole in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        start, hlen = profile_cutter(hole, "cadquery")
        body.append(f"result = result.cut({_cylinder(n, hole, f'{n}_d / 2', start, hlen)})")
        if "cbore_d" in hole:
            c0, clen = profile_cbore_cutter(hole, n, "cadquery")
            body.append(
                f"result = result.cut({_cylinder(n, hole, f'{n}_cbore_d / 2', c0, clen)})"
            )
    return _finish(_head(source), _number_lines(numbers, lists), [_EXTRUDE_FN], body)


def _render_csg(feats: dict, tol: float, source: str | None) -> str:
    del tol
    csg = feats["csg"]
    positives = csg["positive"]
    names = _names(positives)
    numbers, lists = _param_rows(feats, 0.01)
    union = csg.get("combine") == "union"
    overlap = len(positives) > 1 and not union
    body = []
    vars_ = []
    for spec, nm in zip(positives, names):
        var = nm["var"]
        vars_.append(var)
        origin, length = _grown_expr(nm["origin"], nm["len"], spec, union)
        origin, length = _span(origin, length, overlap)
        body.extend(_extrude_lines(var, spec["axis"], origin, length, nm["outline"],
                                   spec.get("cutouts") or [], nm["cutout"]))
        if union and spec.get("chamfer_at") and spec.get("axis") == "Z":
            body.append(
                f"{var} = _edge_chamfer({var}, {nm['outline']}, {nm['origin']}, "
                f"{nm['stem']}_chamfer, {nm['stem']}_chamfer_at)"
            )
        if union:
            for j, _clip in enumerate(spec.get("clips") or [], 1):
                pre = f"{nm['stem']}_clip{j}"
                body.append(
                    f"{var} = {var}.cut(_clip({pre}_x, {pre}_y, {pre}_z, {pre}_d))"
                )
    if len(vars_) > 1:
        body.append(f"result = {vars_[0]}")
        op = "fuse" if union else "common"
        for var in vars_[1:]:
            body.append(f"result = result.{op}({var})")
    for i, cut in enumerate(csg["cuts"], 1):
        origin, length = _span(f"cut_{i}_0", f"cut_{i}_len", True)
        body.extend(_extrude_lines(
            f"cut_{i}_solid", cut["axis"], origin, length, f"cut_{i}",
            cut.get("cutouts") or [], f"cut_{i}_cutout"))
        body.append(f"result = result.cut(cut_{i}_solid)")
    for i, hole in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        lo, hi = _hole_limits(hole, i, positives, names, "cadquery", union)
        body.extend(_hole_lines(n, hole, lo, hi))
    for i, cone in enumerate(csg.get("cones") or [], 1):
        n = f"cone_{i}"
        body.append(
            f"result = result.cut(Part.makeCone({n}_r0, {n}_r1, {n}_h, "
            f"App.Vector({n}_x, {n}_y, {n}_z), {_DIR[cone['axis']]}))"
        )
    for i, _clip in enumerate(csg["clips"], 1):
        body.append(f"result = result.cut(_clip(clip_{i}_x, clip_{i}_y, clip_{i}_z, clip_{i}_d))")
    helpers = [_EXTRUDE_FN]
    if csg["clips"] or any(p.get("clips") for p in positives):
        helpers.append(_CLIP_FN)
    if any(p.get("chamfer_at") and p.get("axis") == "Z" for p in positives):
        helpers.append(_CHAMFER_FN)
    return _finish(_head(source), _number_lines(numbers, lists), helpers, body)


def render(feats: dict, tol: float = 0.01, source: str | None = None) -> str:
    """Turn a feature dict into a FreeCAD macro with named parameters."""
    if "csg" in feats:
        return _render_csg(feats, tol, source)
    if "profile" in feats:
        return _render_profile(feats, tol, source)
    return _render_box(feats, tol, source)


def _driver_source(macro: Path, stl: Path, marker: Path) -> str:
    head = (
        f"MACRO = {str(macro)!r}\n"
        f"STL = {str(stl)!r}\n"
        f"MARKER = {str(marker)!r}\n"
        f"TESSELLATE = {_TESSELLATE!r}\n"
    )
    return head + _DRIVER_BODY


def render_mesh(code: str, work: Path, timeout: float = 300.0) -> trimesh.Trimesh:
    """Run ``code`` with ``freecadcmd`` and load the tessellated solid."""
    binary = freecad_bin()
    if not binary:
        raise FileNotFoundError("freecadcmd is not on PATH")
    work.mkdir(parents=True, exist_ok=True)
    macro = work / "part.FCMacro"
    driver = work / "drive.py"
    stl = work / "part.stl"
    marker = work / "part.ok"
    macro.write_text(code)
    driver.write_text(_driver_source(macro, stl, marker))
    env = os.environ.copy()
    env.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        proc = subprocess.run(
            [binary, str(driver)],
            capture_output=True, text=True, timeout=timeout, cwd=str(work), env=env,
        )
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(f"freecadcmd timed out after {timeout:.0f}s") from e
    err = "\n".join(part for part in ((proc.stderr or ""), (proc.stdout or "")) if part).strip()
    if proc.returncode != 0 or not marker.is_file() or not stl.is_file() or stl.stat().st_size == 0:
        tail = err.splitlines()[-1] if err else f"exit {proc.returncode}"
        raise RuntimeError(tail)
    out = trimesh.load(str(stl), force="mesh", process=True)
    if not isinstance(out, trimesh.Trimesh) or len(out.faces) == 0:
        raise RuntimeError("FreeCAD wrote an empty mesh")
    if not out.is_watertight:
        raise RuntimeError("FreeCAD mesh is not a closed solid")
    return out


def check(code: str, m: trimesh.Trimesh, work: Path, n: int = 4000) -> dict:
    """Run the macro and measure how well it reproduces mesh ``m``."""
    out = render_mesh(code, work)
    dev = deviation(m, out, n=n)
    vm = abs(float(m.volume))
    vo = abs(float(out.volume))
    return {"max_deviation": float(dev["hausdorff"]),
            "hausdorff": float(dev["hausdorff"]),
            "rms": float(dev["rms"]),
            "volume_error_pct": float(100.0 * abs(vo - vm) / vm) if vm else float("inf")}


def to_freecad(m: trimesh.Trimesh, regions, tol: float, source: str | None = None,
               verify: bool = True) -> tuple[str, dict | None]:
    """Generate the macro. When ``verify`` and ``freecadcmd`` is installed, run
    it and raise :class:`Unsupported` if it does not reproduce the mesh.

    Returns ``(code, None)`` when FreeCAD is not available: the macro is still
    produced, and the caller reports it as unverified.
    """
    code = render(recognize(m, regions, tol), tol=tol, source=source)
    if not verify or freecad_bin() is None:
        return code, None
    try:
        with tempfile.TemporaryDirectory(prefix="mesh2cad-fc-") as d:
            q = check(code, m, Path(d))
    except Exception as e:  # noqa: BLE001
        raise Unsupported(f"generated macro failed to run: {e}") from e
    problem = deviation_problem(q, tol)
    if problem:
        raise Unsupported(problem)
    return code, q
