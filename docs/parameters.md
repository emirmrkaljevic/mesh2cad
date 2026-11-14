# Editing a parameter

These are the `params` and `edit` commands, the names they accept, and the edits they refuse.

`params` and `edit` close the loop on a script mesh2cad already wrote. They do not fit a mesh.

```
mesh2cad params bracket.py
mesh2cad edit bracket.py --set hole_1_d=8 --set base_x=+2 --in-place -o bracket.step
```

`params` lists every named value in a CadQuery `.py`, an OpenSCAD `.scad`, or a FreeCAD `.FCMacro`: name, current value, kind, and unit when one is known. Kinds are the ones those scripts use (`diameter`, `length`, `depth`, `radius`, `angle`, `position`, `count` for `hole_row_1_n`, `hole_grid_1_nx`, `hole_grid_1_ny`, and `hole_circle_1_n`, and `points` for an outline or a cutout). Units are millimetres, or degrees for an angle (`hole_1_csk_angle`, `hole_row_1_angle`, `hole_circle_1_angle`, `boss_1_draft_deg`). A count has no unit. Names are the ones in the [CadQuery](cadquery.md), [OpenSCAD](openscad.md), and [FreeCAD](freecad.md) scripts: `base_x`, `base_y`, `base_t`, `leg_x`, `leg_y`, `hole_1_d`, `hole_row_1_n`, `hole_row_1_pitch`, `hole_grid_1_pitch_x`, `hole_circle_1_n`, `hole_circle_1_r`, `hole_circle_1_angle`, `slot_1_len`, `slot_1_w`, `slot_1_angle`, `base_hole_1_d`, `base_hole_row_1_n`, `leg_hole_1_d`, `hole_1_depth`, `hole_1_cbore_d`, `fillet_r`, `edge_fillet_r`, `extrude_len`, `outline`, and the same names with an axis or an index (`extrude_z_len`, `outline_z`, `cut_1`).

`edit` replaces the numeric token and leaves every other byte of the file as it was. `hole_1_d=8` sets that value. `base_x=+2` adds 2 mm, and `base_x=-1` subtracts 1 mm. A leading `+` or `-` is always a delta, so it is not how you type an absolute negative: move `base_z0` by a delta from its current value. `--in-place` rewrites the script. `--output-script new.py` writes a copy. With neither flag the script on disk stays as it is. `-o` writes the solid. CadQuery and FreeCAD write `.step` or `.stl`. OpenSCAD writes `.stl`.

The command prints volume and the bounding box before and after, and whether the result is one valid solid.

Limits, plainly:

- Only scripts mesh2cad generated, the ones with the parameter block it writes. A hand-written CadQuery file is refused.
- Outline and cutout point lists are listed and are not edited as a single number. Changing a corner of `outline` is still a text edit.
- An unknown name is refused, and the message lists the names that exist.
- A hole wider than the wall it is cut through, a mounting hole that leaves the foot or the upright of an L-bracket, a pattern count, pitch, or bolt-circle radius that puts a hole outside that wall, two holes that overlap, a slot whose round end leaves the plate, a count that is not a whole number from 1 to 200, a counterbore or countersink that no longer clears the through-hole, a fillet that eats the short side of the base, a fillet that meets on the shortest edge of an L-bracket, a leg thickness that is no longer smaller than the outer size, or a draft that closes a boss, is refused. So is a render that is empty or is not one valid solid. Nothing is written in those cases.
- The wall check uses the rectangular base, the L profile the bracket script extruded, or the outline of another extrusion. It is the obvious case (the hole does not fit in that face), not a full constraint solver.
- CadQuery runs in-process (`pip install 'mesh2cad[cadquery]'`). OpenSCAD and FreeCAD use the `openscad` and `freecadcmd` programs on `PATH`. If the program is missing, the script can still be rewritten and the command says the solid was not rendered. `-o` in that case writes nothing.
- Values from `--set` are numbers substituted into the existing literal. The script is executed only the same way a render already executes it.
- This does not build a PartDesign body or a constrained sketch. The number changes, then the same script runs again.

The same edits are available from the MCP server as `list_params` and `edit_params`. See [MCP](mcp.md).

Back to the [README](../README.md).
