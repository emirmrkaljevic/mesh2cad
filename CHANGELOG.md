# Changelog

## 0.10.0 - 2026-10-21
- Parametric CadQuery and OpenSCAD export now write a small boolean of straight extrusions along X, Y, and Z. Two or three outlines are intersected (`extrude_z_len`, `outline_z`, and the same for X and Y). One outline keeps `extrude_len` and `outline`. A straight cut along another axis is `cut_1_len`, `cut_1_0`, and `cut_1`. A round through-hole stays `hole_N_d`. A counterbored through-hole adds `hole_N_cbore_d` and `hole_N_cbore_depth`, measured from the face the large diameter opens on. Up to four faces tilted to all three axes are planar clips (`clip_1_x`, `clip_1_y`, `clip_1_z`, `clip_1_d`).
- The box script and the single-extrusion script are still written when they already match the mesh, so those parameter names are unchanged. A rectangular plate with a counterbored hole stays `base_x`, `hole_1_d`, `hole_1_cbore_d`, and `hole_1_cbore_depth`.
- Extrusions that are intersected overlap by 0.02 mm in the script so the shared end faces stay a manifold. The other outline stops at the true end, so the solid does not grow.
- Still refused, and the message says which: a boss on an extruded outline; a countersink, chamfer, or draft on an extrusion; a rounded edge that is not the four vertical corners of a rectangular base; a hole that does not pass through; a stepped hole that is not a counterbored through-hole; more than four faces tilted to X, Y, and Z; a recess that is not a straight cut along X, Y, or Z; a sphere; an unfitted region. A flat face that the fitter called a sphere, with most of its vertices on one axis-aligned plane, is read as that plane.
- endstop-block (308 triangles, 13.6×9×14 mm) exports both scripts. The solid is the intersection of the Z, X, and Y silhouettes, an Ø3.3 hole along X with an Ø6.2 counterbore 8 mm deep from the −X end, and one planar clip. The self-check deviation is 0.0421 mm, under the limit of 0.1 mm. plug-aligner is unchanged at 0.0000 mm. y-belt-idler, y-belt-holder, y-motor-holder, and y-rod-holder are still refused: "part has unfitted (faceted) regions".
- Tests: the intersection of two orthogonal extrusions, a block minus a cut extrusion along another axis, and a counterbored through-hole, clean and at 0.05 mm noise, for both exporters. Editing a hole diameter and an extrusion length changes the volume by the analytic amount.

## 0.9.0 - 2026-10-20
- Parametric CadQuery and OpenSCAD export now write a straight extrusion of a 2D outline along X, Y, or Z. `extrude_len` and `extrude_0` are the length and the start. `outline` is the point list in the plane perpendicular to that axis. A polygonal through-slot is `cutout_N`. A round through-hole stays `hole_N_d` with the same position parameters as on a box.
- The box, pad, fillet, countersink, and side-hole script is the one that is written when it already matches the mesh, so those parameters are unchanged.
- Still refused, and the message says which: a boss on an extruded outline; a countersink, chamfer, or draft on that outline; a rounded edge; two extrusions stacked on the same axis; a hole that does not pass through; a sphere; an unfitted region; a cross-section that changes along the part.
- plug-aligner (48 triangles, 7.3×17.7×5.5 mm) exports both scripts. The self-check deviation is 0.0000 mm. endstop-block does not: along X, 4 faces are tilted (largest 65 mm²); along Y, 4 faces are tilted (largest 73 mm²); along Z, 3 faces are tilted (largest 73 mm²). The other four Prusa parts still have unfitted regions. `benchmarks/real_parts.py` prints this next to the STEP table. A tighter tolerance that turns the endstop's short counterbore into a sphere is not the reason that gets reported.
- Tests: an L extruded along Y with a through hole, a chamfered outline extruded along Z with a hole, and a plate with a rectangular slot, clean and at 0.05 mm noise, for both exporters. Editing one hole diameter and `extrude_len` changes the volume by the analytic amount.

## 0.8.0 - 2026-10-19
- Clean CAD tessellations are segmented on the raw normals. A mesh whose dihedral angles sit on a few exact values (a CAD export with no noise, including faces that are a single quad) is no longer Taubin-smoothed, and a planar face of one or two triangles is kept. A boundary run of a single mesh edge is kept when the whole loop is coarse; a one-edge flicker between two long runs of the same region is still dropped.
- Two coplanar patches merge only when they share an edge, so disconnected faces on one plane stay separate and the exact builder does not see a non-manifold vertex.
- A sphere fit is rejected when the patch has fewer than eight normal directions. The side walls of a rectangular step lie on a circumsphere; they stay planes.
- On a clean mesh, cylinder strips are chained only across a shallow facet step, and a sphere has to beat that cylinder on the vertices by a real margin. Caps are not absorbed into a cylinder they only touch at the rim.
- `benchmarks/real_parts.py` downloads six Original Prusa i3 MK3S printed parts (Prusa Research, GPL-2.0, not vendored) and prints a fidelity table. CI skips it. `tests/test_coarse.py` converts cadquery solids tessellated once, with no subdivision.
- plug-aligner (48 triangles, 14 planes) is exact, Hausdorff 0, 0.6 s. endstop-block is exact: 15 planes and two cylinders (R 1.65 mm and 3.10 mm), Hausdorff 0.017 mm, volume error 0.20%, 2.0 s. On 0.7.0 both were faceted-hybrid (the aligner had no fitted region; the endstop was 75 planar faces).
- y-belt-idler, y-belt-holder, y-motor-holder and y-rod-holder stay faceted-hybrid or faceted. Each has unfitted patches (edge rounds and belt profiles; 0.5% to 6.6% of the area) and cylinders split by a groove, which is outside the vocabulary, so the exact solid is refused. They finish in 12–39 s. On 0.7.0 the holder, motor holder and rod holder did not finish within 300 s: segmentation was a few seconds, and the time was UnifySameDomain on the faceted fallback after the fit had given up.

## 0.7.0 - 2026-10-19
- Parametric CadQuery and OpenSCAD export now include through-holes, blind holes, and straight cylindrical bosses whose axis is the box X or Y axis, as well as +Z. A side hole is `hole_N_d` plus the two coordinates of its axis (`hole_N_y` and `hole_N_z` for an X hole; `hole_N_x` and `hole_N_z` for a Y hole). A blind hole keeps `hole_N_depth`, measured from the face it opens on. A side boss is `boss_N_d` and `boss_N_h` on that same pattern.
- Countersinks, drafts, and chamfers still have to be vertical. An axis that is not aligned with X, Y, or Z is refused. The script is still executed or rendered and compared with the mesh before it is written (same deviation limits); a mismatch is not written.
- Tests: a 60×40×36 block with an Ø8 through hole along X, an Ø6 blind hole along Y, and an Ø10 boss on the +X face, clean and at 0.05 mm noise, plus a side-hole diameter edit whose volume change matches the analytic cylinder, for both exporters. With no noise the scripts match the mesh to 0.006 mm (CadQuery volume error 0.004%, OpenSCAD 0.003%). At 0.05 mm noise the check against the noisy mesh is 0.139 mm / 0.015% (CadQuery) and 0.149 mm / 0.016% (OpenSCAD), under the existing limit max(4·tol, 0.1) = 0.20 mm. The exact STEP builder does not sew this blind hole (the end disk stays in the cylinder region); STEP falls back to a valid faceted solid, Hausdorff 0.004 mm.

## 0.6.0 - 2026-10-18
- Parametric OpenSCAD export: `mesh2cad convert in.stl --openscad out.scad` writes a Customizer-friendly `.scad` file for the same Z-up parts as the CadQuery script. Parameter names and values match (`base_t`, `hole_1_d`, `hole_1_csk_d`, `boss_1_draft_deg`, `boss_1_chamfer`, `pad_1_w`, `fillet_r`, …).
- Plain CSG: `cube` and `cylinder` (`r1`/`r2` for a draft, chamfer, or countersink), `union` / `difference`, and `offset(r=fillet_r)` for the four vertical corners of a rectangular base. A `// [min:max]` slider is written only where the part has a real bound (hole diameter, fillet, draft, countersink angle, chamfer, pad size). Positions stay plain numbers.
- The script is rendered with the `openscad` CLI and compared to the mesh (sampled Hausdorff and volume) before it is written. The same deviation limits as CadQuery apply. If `openscad` is not installed the file is still written and the command reports it unverified; OpenSCAD is not a Python dependency. Parts outside the vocabulary exit with code 2 and a reason. STEP output is unchanged.
- MCP tool `to_openscad(input_path, output_path?, tol?)`.
- Tests cover each vocabulary element on the synthetic parts (clean and 0.05 mm noise), a hole-diameter edit re-rendered with OpenSCAD, and the unsupported path.

## 0.5.0 - 2026-10-17
- Constant-radius fillets between two planes are rebuilt as one partial cylinder, tangent to both walls and trimmed by the two end caps. Covers rounded vertical corners of a plate or block, and the inner corner of an L-bracket as well as the outer ones.
- On the synthetic plate (r=5 mm), block (r=2 mm) and L-bracket (r=3 mm) this is exact at 0 and 0.05 mm noise (3/3 seeds). The radius is taken from the cylinder tangent to both walls that fits the arc, so wall triangles left in the patch do not move it. With no noise the radius error is 0.009 / 0.005 / 0.006 mm. At 0.05 mm the worst error is 0.041 mm on the r=2 corners, 0.100 mm on the plate and 0.131 mm on the bracket. At 0.2 mm the bracket falls back on every seed; the plate is exact on 2 of 3 seeds and one of those solids is missing a corner; the block stays valid but one seed keeps only three of the four rounds. The fallback stays a valid solid.
- CadQuery export writes `fillet_r` and `.edges("|Z").fillet(fillet_r)` when the four vertical corners of a rectangular base share one radius, so the radius is a named parameter. An L-bracket fillet stays on the STEP path.
- A fillet that runs into a spherical corner, and a fillet between a plane and a coaxial cylinder (a torus: boss base, rounded hole edge), are not rebuilt.

## 0.4.2 - 2026-10-16
- Narrow chamfers and countersinks that noise breaks into planar strips are gathered back onto the neighbouring cylinder. A chamfered boss stays exact through 0.2 mm noise; a countersink stays exact through 0.05 mm, and on one of three seeds at 0.2 mm.
- When exact sewing fails, the faceted and faceted-hybrid fallbacks are sewn at a tight tolerance and forced onto one parameter range, so the STEP is a valid solid. Coplanar triangles are merged only when that merge stays valid. UnifySameDomain still runs in an isolated process, so a crash keeps the un-merged solid.

## 0.4.1 - 2026-10-16
- Fix: `convert` segfaulted on some parts at the default `--tol` (e.g. a countersunk plate) inside OCCT's UnifySameDomain. Without an explicit `--tol`, tighter fits (0.4×, 0.2×) are now tried before falling back, and the fallback runs face merging in an isolated child process, keeping the un-merged valid solid if OCCT crashes.
- Fix: a fresh `pip install` pulled cadquery-ocp 8.x, whose API breaks the fallback path (`TopoDS.Shell_s`); pinned `cadquery-ocp>=7.8,<8`.
- Regression test for the countersunk plate at default and loose tolerance.

## 0.4.0 - 2026-10-16
- Parametric CadQuery export: `mesh2cad convert in.stl --cadquery out.py` writes an editable script with named parameters (`base_t`, `hole_1_d`, `hole_1_csk_d`, `boss_1_h`, `boss_1_draft_deg`, `boss_1_chamfer`, `pad_1_w`, …). Covers a box/plate base, rectangular pads, through/blind/countersunk holes and straight/drafted/chamfered bosses.
- The generated script is executed and checked against the mesh before it is written (max deviation ≤ max(4·tol, 0.1), volume error ≤ 1 %); unsupported parts exit with code 2 and a reason.
- MCP tool `to_cadquery(input_path, output_path?, tol?)`.
- New optional extra `mesh2cad[cadquery]`; README install now points at PyPI.
- Tests: script reproduces 8 reference parts (clean and 0.05 mm noise) against ground truth, parameter edits change the part as expected, unsupported parts are rejected.
- Benchmark (`benchmarks/run.py`) now includes the cone parts (countersunk plate, tapered boss, chamfered boss) and a cone face column; README table refreshed.
- Docs: Limits and `server.json` description now list cones; Glama score badge and MCP directory links in the README.

## 0.3.0 - 2026-10-15
- Cone primitive: `Cone`, `fit_cone` (constant normal-to-axis angle + tangent planes through the apex), `refine_cone`.
- Segmentation fits cones for chamfers, countersinks and tapered bosses; a coaxial cylinder + cone component is split by normal-to-axis angle before fitting. Cones win over sphere/cylinder fits only when clearly better.
- Exact B-rep conical frustum faces between two planes perpendicular to the axis, with plane ∩ cone circles on the planar faces. Hybrid fallback snaps vertices onto the cone.
- Fidelity report counts `cone` faces. CLI `--no-cones` and MCP `cones` flag (default on).
- Synthetic roundtrip: countersunk plate, 10° tapered boss, chamfered boss.

## 0.2.0 - 2026-10-14
- Sphere primitive: `Sphere`, `fit_sphere` (normal rays + algebraic fit), `refine_sphere`.
- Segmentation emits `Region(surface=Sphere)` for leftover caps and for clusters of small patches that tile a sphere (noisy or faceted ball ends).
- Exact B-rep spherical faces when the region is bounded by plane ∩ sphere circles (cap or zone). Hybrid fallback snaps vertices onto the sphere.
- Fidelity report counts `sphere` faces. CLI `--no-spheres` and MCP `spheres` flag (default on).
- Synthetic roundtrip: hemispherical dome, spherical pocket, and a sphere clipped to a zone between two planes.

## 0.1.0 - 2026-10-14
First release.
- Mesh loading/cleaning for STL, OBJ, PLY (trimesh).
- Segmentation into planes and cylinders (region growing, strip-chain cylinder detection, robust cylinder refinement, region merge/absorb).
- Exact B-rep STEP export via OpenCascade for planar parts with through-holes and bosses; hybrid/faceted fallbacks.
- Fidelity report (Hausdorff, RMS, volume error, face types, validity) as JSON/Markdown.
- CLI: `mesh2cad convert`, `mesh2cad report`, `mesh2cad mcp`.
- MCP server with `convert`, `analyze`, `report` tools.
- Synthetic self-test benchmark (6 parts × 3 noise levels × 3 seeds, all exact) in CI.
