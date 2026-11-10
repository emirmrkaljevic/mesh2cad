"""Emit an editable, parametric CadQuery script from a fitted model.

The vocabulary and the parameter list live in :mod:`mesh2cad.parametric`.
Every dimension is a named parameter at the top of the script
(``hole_1_d = 6.0``). Holes and bosses whose axis is X or Y are cylinders
along that axis; +Z features stay on the XY workplane. The four vertical
corners of a rectangular base, when they share one radius, become
``edges("|Z").fillet(fillet_r)``. A box with all twelve edges filleted at one
radius becomes ``edges().fillet(edge_fillet_r)``. A through-hole that clears
those rounds is cut after the fillet. A straight boss that clears them is
fused after the fillet. An L-bracket whose vertical edges, including the
inner corner, share one radius is a polyline extruded and then
``edges("|Z").fillet(fillet_r)``, with ``leg_x`` and ``leg_y`` for the two
thicknesses. A quarter-torus fillet on a vertical boss
or hole becomes ``boss_N_fillet``, ``hole_N_rim_r``, ``hole_N_bot_rim_r``, or
``hole_N_floor_r``, applied with ``Solid.fillet`` on that circular edge. A
straight extrusion of a 2D outline is a
``polyline`` extruded by ``extrude_len``. A part that is a few of those
extrusions (their intersection, a union of partial-depth extrusions, or one
minus cuts along the other axes, plus a counterbored hole or a planar clip)
is a short boolean of the same pieces. Anything outside that vocabulary
raises :class:`Unsupported`; the generated script is checked against the input
mesh before it is returned.
"""
from __future__ import annotations

import numpy as np
import trimesh

from .extrude import (profile_cbore_cutter, profile_cutter, profile_parameters,
                      profile_point_lists)
from .parametric import (Unsupported, axis_point_expr, boss_span, bracket_parameters,
                         bracket_profile_exprs, cadquery_axis_dir, cbore_cut_span,
                         chamfer_is_square, deviation_problem, extract, hole_cut_span,
                         named_parameters, recognize)

__all__ = ["Unsupported", "check", "extract", "render", "run_script", "to_cadquery"]

_CIRCLE_EDGES = '''
def _circle_edges(shape, x, y, z, radius):
    """Circular edges of ``shape`` whose centre and radius match."""
    found = []
    for edge in shape.Edges():
        if edge.geomType() != "CIRCLE":
            continue
        center = edge.arcCenter()
        if (abs(center.x - x) < 1e-4 and abs(center.y - y) < 1e-4
                and abs(center.z - z) < 1e-4 and abs(edge.radius() - radius) < 1e-3):
            found.append(edge)
    if not found:
        raise RuntimeError("fillet edge not found")
    return found
'''.strip("\n")


def _needs_circle_edges(feats: dict) -> bool:
    if any("fillet" in b for b in feats.get("bosses", ())):
        return True
    return any("rim_r" in h or "rim_lo" in h or "floor_r" in h for h in feats.get("holes", ()))


def _cq_fillet(point: str, radius: str, fillet: str) -> str:
    return (f"result = result.newObject([result.val().fillet({fillet}, "
            f"_circle_edges(result.val(), {point}, {radius}))])")


def _fmt_points_py(pts) -> str:
    rows = ["["]
    for x, y in pts:
        rows.append(f"    ({x!r}, {y!r}),")
    rows.append("]")
    return "\n".join(rows)


def _profile_workplane(axis: str) -> str:
    """Workplane whose +x/+y are the outline axes and whose normal is the extrusion."""
    if axis == "Z":
        return 'cq.Workplane("XY", origin=(0, 0, extrude_0))'
    if axis == "Y":
        # normal -Y, origin at the far end, so extrude(extrude_len) runs back to extrude_0.
        # local +x is world X and local +y is world Z, which is the (x, z) outline.
        return ("cq.Workplane(cq.Plane(origin=cq.Vector(0, extrude_0 + extrude_len, 0),"
                " xDir=cq.Vector(1, 0, 0), normal=cq.Vector(0, -1, 0)))")
    return ("cq.Workplane(cq.Plane(origin=cq.Vector(extrude_0, 0, 0),"
            " xDir=cq.Vector(0, 1, 0), normal=cq.Vector(1, 0, 0)))")


def _render_profile(feats: dict, tol: float, source: str | None) -> str:
    pr = feats["profile"]
    numbers = profile_parameters(feats, tol)
    lists = profile_point_lists(feats)
    names = [n for n, _, _ in numbers] + [n for n, _, _ in lists]
    width = max(len(n) for n in names)
    wp = _profile_workplane(pr["axis"])
    body = [f"result = ({wp}",
            "          .polyline(outline).close().extrude(extrude_len))"]
    for i, _loop in enumerate(pr["cutouts"], 1):
        body.append(f"result = result.cut({wp}")
        body.append(f"                     .polyline(cutout_{i}).close().extrude(extrude_len))")
    for i, hole in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        start, length = profile_cutter(hole, "cadquery")
        body.append(f"result = result.cut(cq.Workplane().add(cq.Solid.makeCylinder("
                    f"{n}_d / 2, {length}, cq.Vector({axis_point_expr(hole['axis'], n, start)}), "
                    f"cq.Vector({cadquery_axis_dir(hole['axis'])}))))")
        if "cbore_d" in hole:
            c0, clen = profile_cbore_cutter(hole, n, "cadquery")
            body.append(f"result = result.cut(cq.Workplane().add(cq.Solid.makeCylinder("
                        f"{n}_cbore_d / 2, {clen}, cq.Vector({axis_point_expr(hole['axis'], n, c0)}), "
                        f"cq.Vector({cadquery_axis_dir(hole['axis'])}))))")
    head = ['"""Parametric CadQuery model generated by mesh2cad.',
            "",
            f"Source: {source}" if source else "Edit the parameters below and re-run.",
            "Units are the mesh units (usually mm). Run: python this_file.py",
            '"""',
            "import math", "", "import cadquery as cq", "", "# --- parameters ---"]
    params = [f"{n.ljust(width)} = {v!r}" + (f"  # {note}" if note else "") for n, v, note in numbers]
    for name, pts, note in lists:
        if note:
            params.append(f"# {note}")
        params.append(f"{name} = {_fmt_points_py(pts)}")
    tail = ["", "if __name__ == \"__main__\":",
            "    from pathlib import Path",
            "    cq.exporters.export(result, str(Path(__file__).with_suffix(\".step\")))",
            "    print(\"wrote\", Path(__file__).with_suffix(\".step\"))", ""]
    return "\n".join(head + params + ["", "# --- model ---"] + body + tail)


def _render_bracket(feats: dict, tol: float, source: str | None) -> str:
    b = feats["bracket"]
    pts = bracket_profile_exprs(int(b["sx"]), int(b["sy"]))
    # two lines keeps the sharp outline readable next to the parameters
    row_a = ", ".join(f"({x}, {y})" for x, y in pts[:3])
    row_b = ", ".join(f"({x}, {y})" for x, y in pts[3:])
    P = bracket_parameters(feats, tol)
    w = max(len(n) for n, _, _ in P)
    head = ['"""Parametric CadQuery model generated by mesh2cad.',
            "",
            f"Source: {source}" if source else "Edit the parameters below and re-run.",
            "Units are the mesh units (usually mm). Run: python this_file.py",
            '"""',
            "import math", "", "import cadquery as cq", "", "# --- parameters ---"]
    params = [f"{n.ljust(w)} = {v!r}" + (f"  # {note}" if note else "") for n, v, note in P]
    body = [
        "result = (cq.Workplane(\"XY\", origin=(base_x0, base_y0, base_z0))",
        f"          .polyline([{row_a},",
        f"                     {row_b}])",
        "          .close().extrude(base_t)",
        "          .edges(\"|Z\").fillet(fillet_r))",
    ]
    tail = ["", "if __name__ == \"__main__\":",
            "    from pathlib import Path",
            "    cq.exporters.export(result, str(Path(__file__).with_suffix(\".step\")))",
            "    print(\"wrote\", Path(__file__).with_suffix(\".step\"))", ""]
    return "\n".join(head + params + ["", "# --- model ---"] + body + tail)


def render(feats: dict, tol: float = 0.01, source: str | None = None) -> str:
    """Turn a feature dict into a CadQuery script with named parameters."""
    if "bracket" in feats:
        return _render_bracket(feats, tol, source)
    if "csg" in feats:
        from .csg import render_cadquery
        return render_cadquery(feats, tol, source)
    if "profile" in feats:
        return _render_profile(feats, tol, source)
    P = named_parameters(feats, tol)
    if "edge_fillet_r" in feats:
        # every edge of the box, including the horizontal ones, at one radius
        body = ["result = (cq.Workplane(\"XY\", origin=(base_cx, base_cy, base_z0))",
                "          .rect(base_x, base_y).extrude(base_t)",
                "          .edges().fillet(edge_fillet_r))"]
    else:
        body = ["result = (cq.Workplane(\"XY\", origin=(base_cx, base_cy, base_z0))",
                "          .rect(base_x, base_y).extrude(base_t))"]
        if "fillet_r" in feats:
            body.append("result = result.edges(\"|Z\").fillet(fillet_r)")
    body.append("top = base_z0 + base_t")
    for i, d in enumerate(feats["pads"], 1):
        n = f"pad_{i}"
        body.append(f"result = result.union(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top))"
                    f".rect({n}_w, {n}_l).extrude({n}_h))")
    for i, d in enumerate(feats["bosses"], 1):
        n = f"boss_{i}"
        if d.get("axis", "Z") != "Z":
            start, length = boss_span(n, d)
            body.append(f"result = result.union(cq.Workplane().add(cq.Solid.makeCylinder("
                        f"{n}_d / 2, {length}, cq.Vector({axis_point_expr(d['axis'], n, start)}), "
                        f"cq.Vector({cadquery_axis_dir(d['axis'])}))))")
            continue
        if "draft_deg" in d:
            body.append(f"result = result.union(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top))"
                        f".circle({n}_d / 2).extrude({n}_h, taper={n}_draft_deg))")
        elif "chamfer" in d:
            if chamfer_is_square(d, tol):
                cw = ch = f"{n}_chamfer"
            else:
                ch, cw = f"{n}_chamfer_h", f"{n}_chamfer_w"
            body += [f"result = result.union(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top))"
                     f".circle({n}_d / 2).extrude({n}_h - {ch}))",
                     f"result = result.union(cq.Workplane().add(cq.Solid.makeCone({n}_d / 2, "
                     f"{n}_d / 2 - {cw}, {ch}, cq.Vector({n}_x, {n}_y, top + {n}_h - {ch}))))"]
        else:
            body.append(f"result = result.union(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top))"
                        f".circle({n}_d / 2).extrude({n}_h))")
            if "fillet" in d:
                body.append(_cq_fillet(f"{n}_x, {n}_y, top", f"{n}_d / 2", f"{n}_fillet"))
    for i, d in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        if d.get("axis", "Z") != "Z":
            start, length = hole_cut_span(n, d)
            body.append(f"result = result.cut(cq.Workplane().add(cq.Solid.makeCylinder("
                        f"{n}_d / 2, {length}, cq.Vector({axis_point_expr(d['axis'], n, start)}), "
                        f"cq.Vector({cadquery_axis_dir(d['axis'])}))))")
            if "cbore_d" in d:
                c0, clen = cbore_cut_span(n, d)
                body.append(f"result = result.cut(cq.Workplane().add(cq.Solid.makeCylinder("
                            f"{n}_cbore_d / 2, {clen}, cq.Vector({axis_point_expr(d['axis'], n, c0)}), "
                            f"cq.Vector({cadquery_axis_dir(d['axis'])}))))")
            continue
        if "depth" in d:
            body.append(f"result = result.cut(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top - {n}_depth))"
                        f".circle({n}_d / 2).extrude({n}_depth + 1))")
        else:
            body.append(f"result = result.cut(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, base_z0 - 1))"
                        f".circle({n}_d / 2).extrude(base_t + 2))")
        if "csk_d" in d:
            body += [f"_h = ({n}_csk_d - {n}_d) / 2 / math.tan(math.radians({n}_csk_angle / 2))",
                     f"result = result.cut(cq.Workplane().add(cq.Solid.makeCone({n}_d / 2, "
                     f"{n}_csk_d / 2 + 1 * math.tan(math.radians({n}_csk_angle / 2)), _h + 1, "
                     f"cq.Vector({n}_x, {n}_y, top - _h))))"]
        if "cbore_d" in d:
            if d.get("from", 1) > 0:
                body.append(f"result = result.cut(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, "
                            f"top - {n}_cbore_depth)).circle({n}_cbore_d / 2).extrude({n}_cbore_depth + 1))")
            else:
                body.append(f"result = result.cut(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, "
                            f"base_z0 - 1)).circle({n}_cbore_d / 2).extrude({n}_cbore_depth + 1))")
        if "rim_r" in d:
            body.append(_cq_fillet(f"{n}_x, {n}_y, top", f"{n}_d / 2", f"{n}_rim_r"))
        if "rim_lo" in d:
            body.append(_cq_fillet(f"{n}_x, {n}_y, base_z0", f"{n}_d / 2", f"{n}_bot_rim_r"))
        if "floor_r" in d:
            floor_z = f"top - {n}_depth" if d.get("floor_end", "low") == "low" else f"base_z0 + {n}_depth"
            body.append(_cq_fillet(f"{n}_x, {n}_y, {floor_z}", f"{n}_d / 2", f"{n}_floor_r"))
    w = max(len(n) for n, _, _ in P)
    head = ['"""Parametric CadQuery model generated by mesh2cad.',
            "",
            f"Source: {source}" if source else "Edit the parameters below and re-run.",
            "Units are the mesh units (usually mm). Run: python this_file.py",
            '"""',
            "import math", "", "import cadquery as cq", "", "# --- parameters ---"]
    params = [f"{n.ljust(w)} = {v!r}" + (f"  # {note}" if note else "") for n, v, note in P]
    model = ["", "# --- model ---"]
    if _needs_circle_edges(feats):
        model += [_CIRCLE_EDGES, ""]
    model += body
    tail = ["", "if __name__ == \"__main__\":",
            "    from pathlib import Path",
            "    cq.exporters.export(result, str(Path(__file__).with_suffix(\".step\")))",
            "    print(\"wrote\", Path(__file__).with_suffix(\".step\"))", ""]
    return "\n".join(head + params + model + tail)


def run_script(code: str):
    """Execute a generated script and return its ``result`` as a trimesh."""
    ns: dict = {"__name__": "mesh2cad_cq"}
    exec(compile(code, "<mesh2cad-cadquery>", "exec"), ns)
    v, f = ns["result"].val().tessellate(0.01, 0.1)
    return trimesh.Trimesh(np.array([q.toTuple() for q in v]), np.array(f), process=True)


def check(code: str, m: trimesh.Trimesh, n: int = 4000) -> dict:
    """Run the script and measure how well it reproduces mesh ``m``."""
    out = run_script(code)
    a, _ = trimesh.sample.sample_surface(out, n, seed=0)
    b, _ = trimesh.sample.sample_surface(m, n, seed=0)
    d1 = np.abs(trimesh.proximity.signed_distance(m, a))
    d2 = np.abs(trimesh.proximity.signed_distance(out, b))
    vm = abs(m.volume)
    return {"max_deviation": float(max(d1.max(), d2.max())),
            "rms": float(np.sqrt(np.mean(np.r_[d1, d2] ** 2))),
            "volume_error_pct": float(100 * abs(abs(out.volume) - vm) / vm)}


def to_cadquery(m: trimesh.Trimesh, regions, tol: float, source: str | None = None,
                verify: bool = True) -> tuple[str, dict | None]:
    """Generate the script; when ``verify`` and cadquery is installed, run it and
    raise :class:`Unsupported` if it does not reproduce the mesh."""
    code = render(recognize(m, regions, tol), tol=tol, source=source)
    if not verify:
        return code, None
    try:
        import cadquery  # noqa: F401
    except ImportError:
        return code, None
    try:
        q = check(code, m)
    except Exception as e:  # noqa: BLE001
        raise Unsupported(f"generated script failed to run: {e}") from e
    problem = deviation_problem(q, tol)
    if problem:
        raise Unsupported(problem)
    return code, q
