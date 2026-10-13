"""mesh2cad command line."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__


def _convert(a) -> int:
    from .convert import convert_file
    from .report import dumps, to_markdown
    out = a.output or str(Path(a.input).with_suffix(".step"))
    res = convert_file(a.input, out, tol=a.tol, cylinders=not a.no_cylinders)
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
    r = convert_mesh(load_mesh(a.input), tol=a.tol).report
    print(dumps(r) if a.json else to_markdown(r), end="")
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
    c.add_argument("--no-cylinders", action="store_true")
    c.set_defaults(fn=_convert)
    r = sub.add_parser("report", help="analyze a mesh and print a fidelity report without writing STEP")
    r.add_argument("input")
    r.add_argument("--tol", type=float, default=None)
    r.add_argument("--json", action="store_true")
    r.set_defaults(fn=_report)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
