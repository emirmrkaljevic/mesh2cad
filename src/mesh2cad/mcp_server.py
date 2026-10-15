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
    """Convert an STL/OBJ/PLY mesh into a STEP B-rep with fitted planes, cylinders, spheres, and cones.

    Args:
        input_path: path to the mesh file.
        output_path: STEP path to write (default: input with .step suffix).
        tol: fit tolerance in model units (~ mesh noise); default 0.1% of bbox diagonal.
        cylinders: also fit cylinders (holes, bosses).
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


def run() -> None:
    mcp.run()
