"""mesh2cad command line."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__


def _script(kind: str, path: str, a) -> int:
    from .parametric import Unsupported
    kw = dict(tol=a.tol, cylinders=not a.no_cylinders, spheres=not a.no_spheres, cones=not a.no_cones)
    if kind == "cadquery":
        from .convert import cadquery_file as write
        label, missing = "CadQuery", " (install cadquery to verify)"
    elif kind == "openscad":
        from .convert import openscad_file as write
        label, missing = "OpenSCAD", " (unverified: openscad not available)"
    else:
        from .convert import freecad_file as write
        label, missing = "FreeCAD", " (unverified: freecadcmd not available)"
    try:
        info = write(a.input, path, **kw)
    except Unsupported as e:
        print(f"cannot write {label} script: {e}", file=sys.stderr)
        return 2
    q = info["check"]
    chk = f" max_dev={q['max_deviation']:.4f}" if q else missing
    print(f"wrote {path}  params={len(info['parameters'])}{chk}")
    return 0


def _convert(a) -> int:
    from .convert import convert_file
    if a.cadquery:
        rc = _script("cadquery", a.cadquery, a)
        if rc:
            return rc
        if not a.output and not a.openscad and not a.freecad:
            return 0
    if a.openscad:
        rc = _script("openscad", a.openscad, a)
        if rc:
            return rc
        if not a.output and not a.freecad:
            return 0
    if a.freecad:
        rc = _script("freecad", a.freecad, a)
        if rc:
            return rc
        if not a.output:
            return 0
    from .report import dumps, to_markdown
    out = a.output or str(Path(a.input).with_suffix(".step"))
    res = convert_file(a.input, out, tol=a.tol, cylinders=not a.no_cylinders,
                       spheres=not a.no_spheres, cones=not a.no_cones)
    r = res.report
    if a.report:
        p = Path(a.report)
        p.write_text(to_markdown(r) if p.suffix.lower() == ".md" else dumps(r))
    if a.json:
        print(dumps(r))
    else:
        d = r["deviation"]
        print(f"wrote {out}  [{r['method']}] faces={r['output']['faces']} "
              f"valid={r['output']['valid_solid']} hausdorff={d['hausdorff']:.4f} rms={d['rms']:.4f}")
    return 0


def _report(a) -> int:
    from .convert import convert_mesh
    from .io import load_mesh
    from .report import dumps, to_markdown
    r = convert_mesh(load_mesh(a.input), tol=a.tol, cylinders=not a.no_cylinders,
                     spheres=not a.no_spheres, cones=not a.no_cones).report
    print(dumps(r) if a.json else to_markdown(r), end="")
    return 0


def _params(a) -> int:
    from .params import ParamError, format_parameters, list_parameters
    try:
        info = list_parameters(a.script)
    except ParamError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(format_parameters(info), end="")
    return 0


def _parse_set(item: str) -> tuple[str, str]:
    if "=" not in item:
        raise ValueError(f"expected name=value, got {item!r}")
    name, raw = item.split("=", 1)
    name, raw = name.strip(), raw.strip()
    if not name.isidentifier():
        raise ValueError(f"not a parameter name: {name!r}")
    if not raw:
        raise ValueError(f"{name} needs a value")
    return name, raw


def _edit(a) -> int:
    from .params import GeometryError, ParamError, edit_script, format_report
    if a.in_place and a.output_script:
        print("pass only one of --in-place and --output-script", file=sys.stderr)
        return 2
    changes: dict[str, str] = {}
    try:
        for item in a.set or []:
            name, raw = _parse_set(item)
            if name in changes:
                print(f"{name} was passed more than once", file=sys.stderr)
                return 2
            changes[name] = raw
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if not changes:
        print("pass at least one --set name=value", file=sys.stderr)
        return 2
    try:
        info = edit_script(
            a.script, changes, output_path=a.output,
            in_place=a.in_place, output_script=a.output_script,
        )
    except (ParamError, GeometryError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(format_report(info), end="")
    return 0


def _mcp(a) -> int:
    try:
        from .mcp_server import run
    except ImportError:
        print("MCP support needs: pip install 'mesh2cad[mcp]'", file=sys.stderr)
        return 1
    run()
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mesh2cad", description="Convert STL/OBJ/PLY meshes to STEP B-rep CAD")
    ap.add_argument("--version", action="version", version=f"mesh2cad {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("convert", help="convert a mesh to STEP")
    c.add_argument("input")
    c.add_argument("-o", "--output")
    c.add_argument("--tol", type=float, default=None,
                   help="fit tolerance in model units (default: 0.1%% of bbox diagonal); set to scan noise level")
    c.add_argument("--report", help="write fidelity report (.json or .md)")
    c.add_argument("--json", action="store_true", help="print full report as JSON")
    c.add_argument("--cadquery", metavar="OUT.py",
                   help="write an editable parametric CadQuery script instead of STEP "
                        "(add -o to also write STEP); boxes, straight extrusions, "
                        "and a small boolean of those")
    c.add_argument("--openscad", metavar="OUT.scad",
                   help="write an editable parametric OpenSCAD script instead of STEP "
                        "(add -o to also write STEP); same parts as --cadquery")
    c.add_argument("--freecad", metavar="OUT.FCMacro",
                   help="write an editable parametric FreeCAD macro instead of STEP "
                        "(add -o to also write STEP); same parts as --cadquery")
    c.add_argument("--no-cylinders", action="store_true")
    c.add_argument("--no-spheres", action="store_true")
    c.add_argument("--no-cones", action="store_true")
    c.set_defaults(fn=_convert)
    r = sub.add_parser("report", help="analyze a mesh and print a fidelity report without writing STEP")
    r.add_argument("input")
    r.add_argument("--tol", type=float, default=None)
    r.add_argument("--json", action="store_true")
    r.add_argument("--no-cylinders", action="store_true")
    r.add_argument("--no-spheres", action="store_true")
    r.add_argument("--no-cones", action="store_true")
    r.set_defaults(fn=_report)
    p = sub.add_parser("params", help="list editable parameters in a script mesh2cad generated")
    p.add_argument("script", help="CadQuery .py, OpenSCAD .scad, or FreeCAD .FCMacro")
    p.set_defaults(fn=_params)
    e = sub.add_parser("edit", help="change parameters in a generated script and re-render it")
    e.add_argument("script", help="CadQuery .py, OpenSCAD .scad, or FreeCAD .FCMacro")
    e.add_argument("--set", action="append", required=True, metavar="NAME=VALUE",
                   help="absolute value, or a delta when it starts with + or - (repeatable)")
    e.add_argument("-o", "--output", help="write the re-rendered solid (.step or .stl)")
    e.add_argument("--in-place", action="store_true", help="rewrite this script")
    e.add_argument("--output-script", metavar="PATH", help="write the edited script to PATH")
    e.set_defaults(fn=_edit)
    mm = sub.add_parser("mcp", help="run the MCP server (stdio) for AI agents")
    mm.set_defaults(fn=_mcp)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
