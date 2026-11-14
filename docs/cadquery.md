# CadQuery scripts

For a part mesh2cad understands, `mesh2cad convert part.stl --cadquery part.py` writes a Python script whose dimensions are named parameters. This page is the parameter list, the measured examples, and the check that runs before the file is written. OpenSCAD and FreeCAD use the same names.

For simple prismatic parts mesh2cad can also write a **CadQuery Python script** where every dimension is a named parameter, so you change a number and re-run instead of editing B-rep faces:
```
mesh2cad convert bracket.stl --cadquery bracket.py   # add -o bracket.step to also write STEP
python bracket.py                                    # -> bracket.step
```
```python
# --- parameters ---
base_x           = 40.0  # base length (X)
base_y           = 30.0  # base width (Y)
base_t           = 6.0  # base thickness (Z)
...
hole_1_d         = 6.0  # diameter
hole_1_x         = 0.0
hole_1_y         = 0.0
hole_1_csk_d     = 12.0  # countersink diameter
hole_1_csk_angle = 90.0  # countersink angle
```
Make the hole 2 mm bigger: `mesh2cad edit bracket.py --set hole_1_d=8 --in-place -o bracket.step`, or set `hole_1_d = 8.0` and run the script yourself. A plate or block whose four vertical corners share one radius also gets `fillet_r` and `.edges("|Z").fillet(fillet_r)`, so the same edit grows the round. A quarter-torus fillet on a vertical boss or hole is a parameter of that feature: `boss_1_fillet` where a straight boss meets the top face, `hole_1_rim_r` on the top rim of a hole, `hole_1_bot_rim_r` on the bottom rim, and `hole_1_floor_r` where a blind hole's floor meets the wall. `mesh2cad edit part.py --set boss_1_fillet=3` rebuilds that round. A box with all twelve edges filleted at one radius gets `edge_fillet_r` and `.edges().fillet(edge_fillet_r)`. The box size stays `base_x`, `base_y`, and `base_t`. `mesh2cad edit part.py --set edge_fillet_r=4` rebuilds every round. An Ø8 through-hole whose wall clears those rounds is cut after the fillet (`hole_1_d` and its position). A straight boss on the top face whose wall clears those rounds is fused after the fillet (`boss_1_d`, `boss_1_h`, and its position). An L-bracket extruded along Z, with every vertical edge filleted including the inner corner, gets `fillet_r` together with the outer sizes `base_x`, `base_y`, `base_t` and the leg thicknesses `leg_x` and `leg_y`. The script draws that L and calls `.edges("|Z").fillet(fillet_r)`. `mesh2cad edit part.py --set fillet_r=3.5` grows every vertical round, and `--set base_x=44` lengthens the foot. A through-hole along Z in the foot is `base_hole_1_d`, `base_hole_1_x`, and `base_hole_1_y`. One in the upright is `leg_hole_1_d`, `leg_hole_1_x`, and `leg_hole_1_y`. A second hole in the same leg is `base_hole_2` or `leg_hole_2`. The script fillets the L and then cuts the cylinder. `mesh2cad edit part.py --set base_hole_1_d=5` rebuilds that hole.

A hole or boss on a side face uses the same kind of parameter. This is the block below (60×40×36, an Ø8 hole through X, an Ø6 hole 12 deep opening on +Y, an Ø10 boss 8 tall on +X):

```python
boss_1_d     = 10.0  # diameter
boss_1_h     = 8.0  # height from the +X face
boss_1_y     = 10.0
boss_1_z     = 8.0
hole_1_d     = 8.0  # diameter (axis X)
hole_1_y     = 0.0
hole_1_z     = -6.0
hole_2_d     = 6.0  # diameter (axis Y)
hole_2_x     = -12.0
hole_2_z     = 8.0
hole_2_depth = 12.0  # blind depth from the +Y face
```

Set `hole_1_d = 10.0` and run the script again. The volume drops by π(5²−4²)·60 mm³, the analytic change for that through hole.

Three or more plain through-holes along Z, one diameter, centers on a straight line at one pitch, are a single row. A 100×40×6 plate with four Ø6 holes at x = -18, -6, 6 and 18 comes out as:

```python
hole_row_1_n     = 4.0  # count of holes
hole_row_1_pitch = 12.0  # center spacing
hole_row_1_d     = 6.0  # diameter
hole_row_1_x     = -18.0  # start center
hole_row_1_y     = 0.0  # start center
hole_row_1_angle = 0.0  # direction from the start, degrees
```

`mesh2cad edit plate.py --set hole_row_1_n=5` cuts a fifth hole one pitch further along that line. `--set hole_row_1_pitch=+5` keeps four holes and spreads them to 17 mm. An axis-aligned block of at least two by two is `hole_grid_1_nx`, `hole_grid_1_ny`, `hole_grid_1_pitch_x`, `hole_grid_1_pitch_y`, the corner `hole_grid_1_x` / `hole_grid_1_y`, and `hole_grid_1_angle`. A hole that is not part of the pattern stays `hole_1_*`. The same grouping is used on a fully filleted box and on the foot and upright of an L-bracket (`base_hole_row_1_*`, `leg_hole_row_1_*`, and `base_hole_grid_1_*` / `leg_hole_grid_1_*`).

Two holes stay `hole_1` and `hole_2`. They have one gap, so equal spacing says nothing, and moving one should not drag the other. Three or more plain Z through-holes of one diameter, equally spaced on one circle, are `hole_circle_1_n`, `hole_circle_1_d`, `hole_circle_1_r` (pitch circle radius), `hole_circle_1_cx`, `hole_circle_1_cy`, and `hole_circle_1_angle` (degrees, 0 along +X). A center bore stays `hole_1_*`. A square that is also an axis-aligned grid stays `hole_grid_*`. A blind hole, a countersink, a counterbore, and a hole along X or Y stay `hole_N_*`. Spacing is equal when each center is within max(4·tol, 0.25) mm of the ideal position. If the patterned solid misses the usual gate, max deviation max(4·tol, 0.1) and volume error 1%, the script is written with one parameter set per hole.

A straight obround (stadium) through-slot along Z is the same kind of parameter. `slot_N_len` is the distance between the centers of the round ends, `slot_N_w` is the width (that diameter), `slot_N_x` and `slot_N_y` are the center, and `slot_N_angle` is the direction in degrees. Angle 0 runs along +X and angle 90 along +Y. Another angle is kept when the same two ends and two walls still fit. CadQuery cuts `slot2D(slot_N_len + slot_N_w, slot_N_w, slot_N_angle)`, because that call's length is end to end and the named length stays center to center. OpenSCAD hulls the two end cylinders. FreeCAD fuses those two cylinders with the bar between them. `mesh2cad edit plate.py --set slot_1_len=22` lengthens the straight section by 2 mm. The round ends stay put, and the volume drops by 2 times the width times the thickness.

On a 60×40×6 mm plate with one slot (center distance 20 mm, width 8 mm, center at the origin, angle 0) the script is written on 3 of 3 seeds at 0, 0.05, and 0.2 mm noise, for CadQuery and OpenSCAD. The parameters stay 20, 8, and 0. The gate is unchanged. With no noise the limit is 0.100 mm. CadQuery deviates by 0.0047 mm (volume error 0.0027%) and OpenSCAD by 0.0043 mm (0.0020%). At 0.05 mm the limit is 0.200 mm. The worst is seed 1: CadQuery 0.1609 mm (volume error 0.0356%), OpenSCAD 0.1609 mm (0.0349%). At 0.2 mm the limit is 0.800 mm. The worst is seed 0: CadQuery 0.6311 mm (volume error 0.0212%), OpenSCAD 0.6311 mm (0.0205%).

A second plate, same size, has two slots. One runs along X at (-10, 8), center distance 16 mm, width 6 mm. One runs along Y at (12, -6), center distance 14 mm, width 5 mm. With no noise CadQuery deviates by 0.0035 mm (volume error 0.0026%) and OpenSCAD by 0.0032 mm (0.0020%). At 0.05 mm all 3 seeds are written. Seed 1 fits the X slot length as 15.9. The worst is seed 2: CadQuery 0.1676 mm (volume error 0.0084%), OpenSCAD 0.1676 mm (0.0091%). At 0.2 mm none of the three seeds is written. Seeds 0 and 1 are refused with "fillet is not the four vertical corners of a rectangular base". Seed 2 is refused with "vertical fillet is not a base corner". The round ends were not paired, so the miss is not forced into a slot.

A slot at 30 degrees (center distance 18 mm, width 6 mm, centered) is written on the same gate. With no noise the angle is 30.0. CadQuery deviates by 0.0037 mm (volume error 0.0015%) and OpenSCAD by 0.0036 mm (0.0011%). At 0.05 mm all 3 seeds are written. Seed 0 fits the angle as 29.9. The worst is seed 1: CadQuery 0.1610 mm (volume error 0.0014%), OpenSCAD 0.1610 mm (0.0010%), angle 30.0. At 0.2 mm seed 0 is refused with "a rounded edge is not part of an editable outline". Seeds 1 and 2 are written. Seed 1 is length 18.2 and angle 30.3: CadQuery 0.5217 mm (volume error 0.0510%), OpenSCAD 0.5217 mm (0.0506%). Seed 2 is length 18.4 and angle 30.3: CadQuery 0.6870 mm (volume error 0.0640%), OpenSCAD 0.6870 mm (0.0636%).

FreeCAD writes the macro when the same CadQuery check passes. `freecadcmd` was not on PATH for this measurement, so the FreeCAD deviation is not reported. The macro text uses the same `slot_N_*` parameters.

A rectangular through-opening stays `cutout_1` on the extrusion path. It is not renamed `slot_1`. A round hole stays `hole_1_d`. A blind stadium, open on the top and stopped halfway through, is refused with "combined extrusion deviates by 0.162 mm (limit 0.100), volume error 0.00%". Two ends of different diameter (Ø6 and Ø10, 16 mm apart) are refused with "combined extrusion deviates by 0.160 mm (limit 0.100), volume error 0.01%". They are not averaged into one `slot_1_w`.

A part that is one straight extrusion of a 2D outline is the same kind of script. An L, a notched plate, or a block with a slanted face, extruded along X, Y, or Z, becomes `extrude_len` plus an `outline` list. plug-aligner (the printed part in the [real parts](real-parts.md) table, 7.3×17.7×5.5 mm) is a 5.5 mm extrusion along Z of an eight-point outline with a rectangular slot:

```python
extrude_len = 5.5  # length along Z
extrude_0   = 0.0  # start along Z
# outline in the XY plane, (x, y)
outline = [
    (0.0, 0.0),
    (6.1, 0.0),
    (6.1, 6.95),
    (7.3, 6.95),
    (7.3, 10.75),
    (6.1, 10.75),
    (6.1, 17.7),
    (0.0, 17.7),
]
# through opening in the outline
cutout_1 = [
    (0.45, 1.85),
    (5.65, 1.85),
    (5.65, 15.85),
    (0.45, 15.85),
]
```

Set `extrude_len = 7.5` and the part is 2 mm thicker. A corner moves when one pair in `outline` changes. A round through-hole on an extrusion is still `hole_1_d`, with `hole_1_y` and `hole_1_z` when the hole runs along X. `cutout_1` is a polygonal slot that goes all the way through. An extrusion along Y stores the outline as `(x, z)`; along X, as `(y, z)`.

Supported: a Z-up box/plate base, rectangular pads on the top face, holes along X, Y, or Z (through, blind, countersunk when the axis is vertical, and counterbored on any of those axes), bosses along X, Y, or Z (straight, and drafted or chamfered when the axis is vertical), that one vertical corner radius, one radius on all twelve edges of a rectangular box (`edge_fillet_r`), a through-hole cut after that fillet when its wall clears the rounds, a straight boss fused after that fillet when its wall clears the rounds, a quarter-torus fillet on a vertical boss or hole (`boss_N_fillet`, `hole_N_rim_r`, `hole_N_bot_rim_r`, `hole_N_floor_r`), an L-bracket extruded along Z whose vertical edges, including the inner corner, share `fillet_r` (`leg_x` and `leg_y` are the two thicknesses; a through-hole in the foot is `base_hole_N_d` and one in the upright is `leg_hole_N_d`), and one straight extrusion of a 2D outline along X, Y, or Z, including a polygonal through-slot and round through-holes. A straight obround through-slot along Z on a plate is `slot_N_len`, `slot_N_w`, `slot_N_x`, `slot_N_y`, and `slot_N_angle`. A counterbored hole is `hole_N_d` for the through diameter, `hole_N_cbore_d` for the large diameter, and `hole_N_cbore_depth` for how far that large diameter goes in from the face it opens on. A part that is a few of those extrusions is their intersection, or one of them minus straight cuts along the other axes. Two or three outlines use `extrude_z_len` / `outline_z` (and the same for X and Y); a single outline keeps `extrude_len` and `outline`. A cut is `cut_1_len`, `cut_1_0`, and `cut_1`. Up to four faces tilted to X, Y, and Z are planar clips (`clip_1_x`, `clip_1_y`, `clip_1_z`, `clip_1_d`). The box script is the one that is written when it already reproduces the mesh, so a plain plate stays `base_x` / `base_y` / `base_t` and does not become an outline list. A rectangular plate with a counterbored hole stays on those base parameters. Values are rounded to clean numbers (5.997 → 6.0) when that stays within `--tol`. Before writing, the script is executed and compared with the mesh (needs `cadquery` installed); if the part is outside this vocabulary or the script deviates, nothing is written and the command exits with code 2 and the reason. Use STEP output for everything else.

See [OpenSCAD](openscad.md), [FreeCAD](freecad.md), [editing a parameter](parameters.md), [limits](limits.md), and [real parts](real-parts.md).

Back to the [README](../README.md).
