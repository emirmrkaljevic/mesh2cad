# Command line

`mesh2cad` converts a mesh, prints a fidelity report, lists parameters in a script it wrote, or edits one of those parameters and rebuilds the solid. `mesh2cad --version` prints the version.

Input meshes are STL, OBJ, or PLY. Open the STEP file in FreeCAD, Fusion, Onshape, or any other CAD tool. When the fit is exact, the faces are planes, cylinders, and fillets you can select and dimension. `report.md` says how close the result is to the mesh. When a fit is bad, STEP output falls back to a faceted solid instead of guessing.

## convert

```
mesh2cad convert part.stl -o part.step --tol 0.1 --report report.md
mesh2cad convert bracket.stl --cadquery bracket.py
mesh2cad convert plate.stl --openscad plate.scad
mesh2cad convert plate.stl --freecad plate.FCMacro
```

`part.stl` is the mesh. `-o` / `--output` is the STEP path. If you leave it off, the STEP file is written next to the mesh with a `.step` suffix. If you ask only for a script and do not pass `-o`, no STEP file is written.

`--tol` is about the noise or accuracy of the mesh, in model units. The default is 0.1% of the bounding-box diagonal.

`--report` writes a fidelity report. A `.md` path is Markdown. Any other suffix is JSON. `--json` prints the full report as JSON. Without `--json`, the command prints one line: the output path, the method, the face count, whether the solid is valid, and the Hausdorff and RMS deviation.

The report records the method (`exact`, `faceted-hybrid`, or `faceted`), face types, solid validity, two-sided Hausdorff and RMS deviation, volume error, and the percent of area per surface type.

`--cadquery`, `--openscad`, and `--freecad` write an editable script. You can pass more than one, and `-o` still writes STEP. The script is checked against the mesh before it is written. If the part is outside the vocabulary, or the check fails, nothing is written for that script and the command exits with code 2 and the reason. What the scripts contain is in [CadQuery](cadquery.md), [OpenSCAD](openscad.md), and [FreeCAD](freecad.md). The refusal list is in [Limits](limits.md).

`--no-cylinders`, `--no-spheres`, and `--no-cones` turn those fits off. They are on by default.

## report

```
mesh2cad report part.stl --json
```

This runs the same fit and prints the fidelity report. It does not write STEP. `--tol`, `--json`, `--no-cylinders`, `--no-spheres`, and `--no-cones` mean the same thing as on `convert`. Without `--json` the report is Markdown on stdout.

## params

```
mesh2cad params bracket.py
```

Lists every named value in a CadQuery `.py`, an OpenSCAD `.scad`, or a FreeCAD `.FCMacro` that mesh2cad wrote: name, current value, kind, and unit when one is known. A hand-written script is refused, and the command exits with code 2. Kinds, units, and the names are in [Parameters](parameters.md).

## edit

```
mesh2cad edit bracket.py --set hole_1_d=8 --set base_x=+2 --in-place -o bracket.step
```

`--set` is required and can be repeated. `hole_1_d=8` sets that value. `base_x=+2` adds 2 mm, and `base_x=-1` subtracts 1 mm. A leading `+` or `-` is always a delta, so it is not how you type an absolute negative. Move a value that is already negative by a delta from its current value.

`--in-place` rewrites the script. `--output-script new.py` writes a copy. Pass only one of those two. With neither, the script on disk stays as it is. `-o` / `--output` writes the solid. CadQuery and FreeCAD write `.step` or `.stl`. OpenSCAD writes `.stl`.

The command prints volume and the bounding box before and after, and whether the result is one valid solid. An unknown name, a hole that does not fit the wall, or a render that is not one valid solid is refused. Nothing is written in those cases, and the command exits with code 2. The full list of those limits is in [Parameters](parameters.md).

CadQuery runs in-process (`pip install "mesh2cad[cadquery]"`). OpenSCAD and FreeCAD use the `openscad` and `freecadcmd` programs on `PATH`. If the program is missing, the script can still be rewritten and the command says the solid was not rendered. `-o` in that case writes nothing.

## mcp

```
mesh2cad mcp
```

Runs the MCP server on stdio. If the `mcp` package is missing, the command exits with code 1 and tells you to install `mesh2cad[mcp]`. The tools are in [MCP](mcp.md).

Back to the [README](../README.md).
