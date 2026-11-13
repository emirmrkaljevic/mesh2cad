# MCP server

mesh2cad speaks MCP over stdio. Install the extra, then launch it:

```
pip install "mesh2cad[mcp]"
mesh2cad mcp
```

Claude Desktop, Cursor, and any other client that starts a local command can use this config:

```json
{ "mcpServers": { "mesh2cad": { "command": "mesh2cad", "args": ["mcp"] } } }
```

`spheres` and `cones` default on for the tools that take them. The same fits are what the [command line](cli.md) runs.

## convert

`convert(input_path, output_path?, tol?, cylinders?, spheres?, cones?)`

Converts an STL, OBJ, or PLY mesh into a STEP file with fitted planes, cylinders, spheres, cones, and constant-radius fillets. `output_path` is the STEP path. The default is the input path with a `.step` suffix. `tol` is the fit tolerance in model units, about the mesh noise. The default is 0.1% of the bounding-box diagonal. `cylinders`, `spheres`, and `cones` turn those fits on or off. The result is the fidelity report (method, faces, deviation, volume error) and the output path.

## analyze

`analyze(input_path, tol?, cylinders?, spheres?, cones?)`

Segments the mesh and reports the detected surfaces. It does not write STEP. The result has the regions, the area fraction, the method, and the input summary.

## report

`report(input_path, tol?, cylinders?, spheres?, cones?)`

Returns the full Markdown fidelity report for this mesh. It does not write a file.

## to_cadquery

`to_cadquery(input_path, output_path?, tol?)`

Writes an editable parametric CadQuery script. The default path is the input with a `.py` suffix. The script tools return the parameters and, when the script was run or rendered, the measured deviation. A part outside that vocabulary comes back under `error`. The vocabulary is the one in [CadQuery](cadquery.md).

## to_openscad

`to_openscad(input_path, output_path?, tol?)`

Writes an OpenSCAD script for the same parts as `to_cadquery`. The default path ends in `.scad`. `to_openscad` returns `verified: false` when the `openscad` program is absent. The file is still written in that case. A part outside the vocabulary, or a render that does not match, comes back under `error` and nothing is written when the check itself refuses.

## to_freecad

`to_freecad(input_path, output_path?, tol?)`

Writes a FreeCAD macro for the same parts. The default path ends in `.FCMacro`. `to_freecad` returns `verified: false` when `freecadcmd` is absent. Without `freecadcmd`, the same deviation and volume check is run with CadQuery, or with OpenSCAD when CadQuery is not installed. The macro is written only when that check passes. A refusal comes back under `error`, and nothing is written.

## list_params

`list_params(script_path)`

Takes a script mesh2cad already wrote: a CadQuery `.py`, an OpenSCAD `.scad`, or a FreeCAD `.FCMacro`. Returns the name, value, kind, and unit for each parameter. Outline and cutout point lists are included, and they are not single numbers. A file mesh2cad did not generate comes back under `error`.

## edit_params

`edit_params(script_path, changes, output_path?)`

Takes a script mesh2cad already wrote and updates that script in place. A change is a number, or a string that starts with `+` or `-` for a delta. For example, `{"hole_1_d": 8, "base_x": "+2"}` sets the hole diameter to 8 and adds 2 mm to `base_x`.

`output_path` is the STEP or STL. When it is omitted, a `.step` (`.stl` for OpenSCAD) is written next to the script if that backend can run. The result includes the old and new values, the rendered path, validity, the volume and bounding-box delta, and warnings. An unknown name, or an edit that is not one valid solid, comes back under `error`, and nothing is written. The same wall checks as `mesh2cad edit` apply. They are listed in [Parameters](parameters.md).

Listed on [Glama](https://glama.ai/mcp/servers/emirmrkaljevic/mesh2cad) and the [TensorBlock MCP catalog](https://tensorblock.co/mcp/servers/github-emirmrkaljevic-mesh2cad-0f1c0e66).

Back to the [README](../README.md).
