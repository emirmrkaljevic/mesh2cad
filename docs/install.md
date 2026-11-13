# Install

mesh2cad needs Python 3.10 or newer. The library is Python, on top of OpenCascade (`cadquery-ocp`), and pip installs that with the package. The command line tool and STEP export are in the base install. The MCP server and the CadQuery runner are extras. OpenSCAD and FreeCAD are separate programs, not Python packages.

```
pip install mesh2cad              # CLI + STEP
pip install "mesh2cad[mcp]"       # + MCP server
pip install "mesh2cad[cadquery]"  # + run/verify generated CadQuery scripts
```

OpenSCAD scripts are plain text. Verifying one needs the `openscad` program ([openscad.org](https://openscad.org)); it is not a Python dependency. Without it the `.scad` file is still written.

FreeCAD macros are plain text too. Measuring the macro itself needs `freecadcmd` (the FreeCAD console, from [freecad.org](https://www.freecad.org)). It is not a Python dependency. Without it, the same deviation and volume check is run with CadQuery, or with OpenSCAD when CadQuery is not installed. The `.FCMacro` file is written when that check passes, and it is not written when the check refuses. Those numbers are from the other backend, not from FreeCAD. Ubuntu 24.04 does not ship a `freecad` apt package, so install the console build and put `freecadcmd` on `PATH` when you want the macro measured directly.

Command flags are in the [CLI reference](cli.md). The server launch is in [MCP](mcp.md).

Back to the [README](../README.md).
