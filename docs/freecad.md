# FreeCAD macros

The same parts as the [CadQuery scripts](cadquery.md) can be written as a FreeCAD macro. Names and values match. The solid is analytic (cylinders and cones). The macro is a Part script, not a PartDesign body.

The same parts can be written as a FreeCAD macro. Names and values match the CadQuery script. `--freecad` can be combined with `--cadquery`, `--openscad`, and `-o`.
```
mesh2cad convert plate.stl --freecad plate.FCMacro
freecadcmd plate.FCMacro
```
A 40×30×6 plate with a Ø10 hole comes out as below (the file also has a short header and, when you run it, creates a `Part::Feature` and writes `plate.step` next to the macro):

```python
# --- parameters ---
base_x   = 40.0  # base length (X)
base_y   = 30.0  # base width (Y)
base_t   = 6.0  # base thickness (Z)
base_cx  = 0.0
base_cy  = 0.0
base_z0  = -3.0  # bottom face height
hole_1_d = 10.0  # diameter
hole_1_x = 0.0
hole_1_y = 0.0

# --- model ---
result = Part.makeBox(base_x, base_y, base_t, App.Vector(base_cx - base_x / 2, base_cy - base_y / 2, base_z0))
top = base_z0 + base_t
result = result.cut(Part.makeCylinder(hole_1_d / 2, base_t + 2, App.Vector(hole_1_x, hole_1_y, base_z0 - 1), App.Vector(0, 0, 1)))
```

Make the hole 2 mm bigger: set `hole_1_d = 12.0` and run the macro again. A drafted or chamfered boss is `Part.makeCone`. The four vertical corners of a rectangular base are `result.makeFillet(fillet_r, ...)` on the edges parallel to Z. All twelve edges of that box, at one radius, are `result.makeFillet(edge_fillet_r, ...)` on every edge. A through-hole that clears those rounds is cut after that fillet. A boss that clears those rounds is fused after that fillet. The same `makeFillet` rounds the circular edge of a vertical boss base, a hole rim, or a blind-hole floor, with `boss_1_fillet`, `hole_1_rim_r`, or `hole_1_floor_r`. An L-bracket is a polygon extruded along Z, then `makeFillet(fillet_r, ...)` on the edges parallel to Z, which includes the inner corner. `leg_x` and `leg_y` are the two thicknesses. A through-hole in the foot or the upright is `Part.makeCylinder` cut after that fillet (`base_hole_1_d` or `leg_hole_1_d`). A row is a Python loop of that cylinder, with the same `hole_row_1_*` numbers as the CadQuery script. A grid is two nested loops. An obround slot is two cylinders and the bar between them, with the same `slot_1_*` numbers. A hole or boss whose axis is X or Y is the same cylinder, with `App.Vector(1, 0, 0)` or `App.Vector(0, 1, 0)`. A straight extrusion is `_extrude(outline, "Z", extrude_0, extrude_len)` (or `"X"` / `"Y"`). A polygonal slot is a second extrusion subtracted from the first. Several extrusions are `.common()` (intersection) or `.fuse()` (a partial-depth union). A cut extrusion, a hole, a counterbore, a 45 degree chamfer, and a planar clip are `.cut()`. Names match the other scripts (`extrude_len`, `outline`, `hole_1_d`, `hole_1_cbore_d`, `cut_1`, `clip_1_d`).

The solid in FreeCAD is analytic: cylinders and cones, not a `$fn` polygon. The macro is a Part script with named numbers at the top. It does not build a PartDesign body or a constrained sketch, so the edit is the number, then run the macro again.

Same honesty as the other two. Before writing, mesh2cad runs the macro with `freecadcmd` and compares a fine tessellation of the solid to the mesh (sampled Hausdorff and volume, same limits). If `freecadcmd` is not installed, CadQuery or OpenSCAD applies that same check. The file is written only when the check passes, and the command says the FreeCAD solid itself was unverified. If the part is outside this vocabulary, or the check does not match, nothing is written and the command exits with code 2 and the reason. STEP remains the path for everything else.

See [CadQuery](cadquery.md), [OpenSCAD](openscad.md), [editing a parameter](parameters.md), and [limits](limits.md).

Back to the [README](../README.md).
