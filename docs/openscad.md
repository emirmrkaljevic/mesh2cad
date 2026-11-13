# OpenSCAD scripts

The same parts as the [CadQuery scripts](cadquery.md) can be written as OpenSCAD. Names and values match. This page is the generated script, what each feature becomes in OpenSCAD, and the check against the mesh.

The same parts can be written as OpenSCAD, which is what the 3D-printing Customizer on Thingiverse and Printables edits. Names and values match the CadQuery script. `--openscad` and `--cadquery` can be combined, and `-o` still writes STEP.
```
mesh2cad convert plate.stl --openscad plate.scad
openscad -o plate.stl plate.scad
```
A 40×30×6 plate with a Ø10 hole comes out as below, plus a short header. Before:

```openscad
// --- parameters ---
/* [Base] */
// base length (X)
base_x   = 40.0;
// base width (Y)
base_y   = 30.0;
// base thickness (Z)
base_t   = 6.0;
base_cx  = 0.0;
base_cy  = 0.0;
// bottom face height
base_z0  = -3.0;

/* [Holes] */
// diameter
hole_1_d = 10.0; // [0.2:0.1:30]
hole_1_x = 0.0;
hole_1_y = 0.0;

/* [Hidden] */
$fn = 96; // circle segments

// --- model ---
module mesh2cad_part() {
  difference() {
    translate([base_cx, base_cy, base_z0 + base_t / 2])
      cube([base_x, base_y, base_t], center=true);
    translate([hole_1_x, hole_1_y, base_z0 - 1])
      cylinder(h=base_t + 2, d=hole_1_d);
  }
}

mesh2cad_part();
```

After, the hole is 2 mm bigger: one number changes, then re-render:

```openscad
hole_1_d = 12.0; // [0.2:0.1:30]
```

Re-render (F6) and export STL. A drafted boss is a `cylinder` with `r1`/`r2`, a chamfer or countersink is the same, and the four vertical corners of a rectangular base are `offset(r=fillet_r)` on a smaller square, then extruded. A box with every edge rounded is `edge_fillet_r`: three cubes for the flat faces, twelve cylinders for the edges, and eight spheres for the corners. A through-hole that clears those rounds is a `cylinder` subtracted after that union. A boss that clears those rounds is a `cylinder` added after that union. An L-bracket is `offset` around the sharp L, so the outer corners and the inner corner all use `fillet_r`. The slider stops short of half the shortest edge, where those rounds would meet. A through-hole in the foot or the upright is a `cylinder` subtracted after that extrusion (`base_hole_1_d` or `leg_hole_1_d`). Its diameter slider stops at the shorter side of that leg. A row is a `for` loop of that cylinder, using `hole_row_1_n`, `hole_row_1_pitch`, and `hole_row_1_angle`. A grid is two nested loops (`hole_grid_1_nx`, `hole_grid_1_ny`). An obround slot is a `hull()` of two cylinders (`slot_1_len`, `slot_1_w`, `slot_1_angle`). A boss-base fillet, a hole-rim round, and a blind-floor fillet are a `rotate_extrude()` of a square with a quarter circle taken out of it, using `boss_1_fillet`, `hole_1_rim_r`, or `hole_1_floor_r`. The profile overlaps the sharp solid by 0.02 mm so the union stays closed. A hole or boss whose axis is X or Y is the same `cylinder`, turned with `rotate` so it points along that axis; its parameters are the diameter and the two coordinates of the axis (`hole_1_y` / `hole_1_z` for an X hole). A straight extrusion is `linear_extrude(extrude_len)` of `polygon(outline)`, with `polygon(cutout_1)` subtracted in 2D when the outline has a slot. An extrusion along Y or X is that solid turned so the outline sits in the XZ or YZ plane. A counterbored hole is a second `cylinder` of diameter `hole_1_cbore_d` and height `hole_1_cbore_depth`, from the face it opens on. Several extrusions are an `intersection()`, and a cut extrusion or a planar clip is subtracted. Names match the CadQuery script (`extrude_len`, `outline`, `hole_1_d`, `hole_1_cbore_d`, `cut_1`, `clip_1_d`). A slider (`// [min:max]`) is added only where the part gives a bound: hole diameter up to the shorter side of the face it is cut through (or of the outline, on an extrusion), counterbore diameter and depth on that same axis, fillet up to half the short side of the plate, or half the shortest edge of an L-bracket, draft until the tip would close, countersink angle in 60 to 120°, chamfer up to the boss, pad size up to the base, blind depth up to the length of the hole's axis. Positions (`hole_1_x`, `base_z0`, `extrude_0`, …) stay plain numbers. The outline is a point list, not a slider.

Same honesty as CadQuery. Before writing, mesh2cad renders the script with the `openscad` CLI and compares the STL to the mesh (sampled Hausdorff and volume, same limits as the CadQuery check). If `openscad` is not installed the file is still written and the command says it was unverified. If the part is outside this vocabulary, or the render does not match, nothing is written and the command exits with code 2 and the reason. Circles are `$fn = 96` polygons; on these parts that is a few thousandths of a millimetre, not an analytic solid. STEP remains the path for everything else.

See [CadQuery](cadquery.md), [FreeCAD](freecad.md), [editing a parameter](parameters.md), and [limits](limits.md).

Back to the [README](../README.md).
