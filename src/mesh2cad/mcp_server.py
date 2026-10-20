"""MCP server: lets AI agents convert meshes to STEP and inspect fidelity.

Run with `mesh2cad mcp` (stdio transport).
"""
from __future__ import annotations

from pathlib import Path

try:  # mcp >= 2
    from mcp.server.mcpserver import MCPServer as _Server
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP as _Server

mcp = _Server("mesh2cad")


@mcp.tool()
def convert(input_path: str, output_path: str | None = None, tol: float | None = None,
            cylinders: bool = True, spheres: bool = True, cones: bool = True) -> dict:
    """Convert an STL/OBJ/PLY mesh into a STEP B-rep with fitted planes, cylinders, spheres, cones, and constant-radius fillets.

    Args:
        input_path: path to the mesh file.
        output_path: STEP path to write (default: input with .step suffix).
        tol: fit tolerance in model units (~ mesh noise); default 0.1% of bbox diagonal.
        cylinders: also fit cylinders (holes, bosses, plane-plane fillets).
        spheres: also fit spheres (domes, spherical bosses).
        cones: also fit cones (chamfers, countersinks, tapered bosses).
    Returns the fidelity report (method, faces, deviation, volume error) and output path.
    """
    from .convert import convert_file
    out = output_path or str(Path(input_path).with_suffix(".step"))
    r = convert_file(input_path, out, tol=tol, cylinders=cylinders, spheres=spheres, cones=cones).report
    return {"output_path": out, **r}


@mcp.tool()
def analyze(input_path: str, tol: float | None = None, cylinders: bool = True,
            spheres: bool = True, cones: bool = True) -> dict:
    """Segment a mesh and report detected surfaces without writing STEP."""
    from .convert import convert_mesh
    from .io import load_mesh
    r = convert_mesh(load_mesh(input_path), tol=tol, cylinders=cylinders, spheres=spheres, cones=cones).report
    return {"regions": r["regions"], "area_fraction": r["area_fraction"],
            "method": r["method"], "input": r["input"]}


@mcp.tool()
def report(input_path: str, tol: float | None = None, cylinders: bool = True,
           spheres: bool = True, cones: bool = True) -> str:
    """Full Markdown fidelity report for converting this mesh."""
    from .convert import convert_mesh
    from .io import load_mesh
    from .report import to_markdown
    return to_markdown(convert_mesh(load_mesh(input_path), tol=tol, cylinders=cylinders,
                                    spheres=spheres, cones=cones).report)


@mcp.tool()
def to_cadquery(input_path: str, output_path: str | None = None, tol: float | None = None) -> dict:
    """Write an editable parametric CadQuery script (.py) for a simple part.

    Works for a box/plate base with rectangular pads, holes along X, Y, or Z
    (through, blind, and countersunk when the axis is vertical), bosses
    along X, Y, or Z (straight, and drafted or chamfered when the axis is
    vertical), and a straight extrusion of a 2D outline along X, Y, or Z.
    Each dimension is a named parameter (e.g. hole_1_d = 6.0, extrude_len): edit
    it and run the script with Python to regenerate the part. Returns the
    parameters and the measured deviation of the script's output from the mesh.
    Returns an "error" key when the part is outside that vocabulary.
    """
    from .parametric import Unsupported
    from .convert import cadquery_file
    out = output_path or str(Path(input_path).with_suffix(".py"))
    try:
        return cadquery_file(input_path, out, tol=tol)
    except Unsupported as e:
        return {"error": str(e)}


@mcp.tool()
def to_openscad(input_path: str, output_path: str | None = None, tol: float | None = None) -> dict:
    """Write an editable parametric OpenSCAD script (.scad) for a simple part.

    Same parts as to_cadquery: a box/plate base with rectangular pads, holes and
    bosses along X, Y, or Z (countersink, draft, and chamfer only when the axis
    is vertical), one radius on the four vertical corners, and a straight
    extrusion of a 2D outline. Parameter names match the CadQuery script
    (hole_1_d, boss_1_h, fillet_r, extrude_len, ...). Edit one and
    render in OpenSCAD. Returns the parameters, whether the openscad CLI
    verified the script against the mesh, and the measured deviation when it
    did. Without the openscad program the script is still written and verified
    is false. Returns an "error" key when the part is outside that vocabulary.
    """
    from .convert import openscad_file
    from .parametric import Unsupported
    out = output_path or str(Path(input_path).with_suffix(".scad"))
    try:
        return openscad_file(input_path, out, tol=tol)
    except Unsupported as e:
        return {"error": str(e)}


def run() -> None:
    mcp.run()
