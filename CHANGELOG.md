# Changelog

## 0.7.0 — 2026-10-19
- Parametric CadQuery and OpenSCAD export now include through-holes, blind holes, and straight cylindrical bosses whose axis is the box X or Y axis, as well as +Z. A side hole is `hole_N_d` plus the two coordinates of its axis (`hole_N_y` and `hole_N_z` for an X hole; `hole_N_x` and `hole_N_z` for a Y hole). A blind hole keeps `hole_N_depth`, measured from the face it opens on. A side boss is `boss_N_d` and `boss_N_h` on that same pattern.
- Countersinks, drafts, and chamfers still have to be vertical. An axis that is not aligned with X, Y, or Z is refused. The script is still executed or rendered and compared with the mesh before it is written (same deviation limits); a mismatch is not written.
- Tests: a 60×40×36 block with an Ø8 through hole along X, an Ø6 blind hole along Y, and an Ø10 boss on the +X face, clean and at 0.05 mm noise, plus a side-hole diameter edit whose volume change matches the analytic cylinder, for both exporters. With no noise the scripts match the mesh to 0.006 mm (CadQuery volume error 0.004%, OpenSCAD 0.003%). At 0.05 mm noise the check against the noisy mesh is 0.139 mm / 0.015% (CadQuery) and 0.149 mm / 0.016% (OpenSCAD), under the existing limit max(4·tol, 0.1) = 0.20 mm. The exact STEP builder does not sew this blind hole (the end disk stays in the cylinder region); STEP falls back to a valid faceted solid, Hausdorff 0.004 mm.

## 0.6.0 — 2026-10-18
- Parametric OpenSCAD export: `mesh2cad convert in.stl --openscad out.scad` writes a Customizer-friendly `.scad` file for the same Z-up parts as the CadQuery script. Parameter names and values match (`base_t`, `hole_1_d`, `hole_1_csk_d`, `boss_1_draft_deg`, `boss_1_chamfer`, `pad_1_w`, `fillet_r`, …).
- Plain CSG: `cube` and `cylinder` (`r1`/`r2` for a draft, chamfer, or countersink), `union` / `difference`, and `offset(r=fillet_r)` for the four vertical corners of a rectangular base. A `// [min:max]` slider is written only where the part has a real bound (hole diameter, fillet, draft, countersink angle, chamfer, pad size). Positions stay plain numbers.
- The script is rendered with the `openscad` CLI and compared to the mesh (sampled Hausdorff and volume) before it is written. The same deviation limits as CadQuery apply. If `openscad` is not installed the file is still written and the command reports it unverified; OpenSCAD is not a Python dependency. Parts outside the vocabulary exit with code 2 and a reason. STEP output is unchanged.
- MCP tool `to_openscad(input_path, output_path?, tol?)`.
- Tests cover each vocabulary element on the synthetic parts (clean and 0.05 mm noise), a hole-diameter edit re-rendered with OpenSCAD, and the unsupported path.

## 0.5.0 — 2026-10-17
- Constant-radius fillets between two planes are rebuilt as one partial cylinder, tangent to both walls and trimmed by the two end caps. Covers rounded vertical corners of a plate or block, and the inner corner of an L-bracket as well as the outer ones.
- On the synthetic plate (r=5 mm), block (r=2 mm) and L-bracket (r=3 mm) this is exact at 0 and 0.05 mm noise (3/3 seeds). The radius is taken from the cylinder tangent to both walls that fits the arc, so wall triangles left in the patch do not move it. With no noise the radius error is 0.009 / 0.005 / 0.006 mm. At 0.05 mm the worst error is 0.041 mm on the r=2 corners, 0.100 mm on the plate and 0.131 mm on the bracket. At 0.2 mm the bracket falls back on every seed; the plate is exact on 2 of 3 seeds and one of those solids is missing a corner; the block stays valid but one seed keeps only three of the four rounds. The fallback stays a valid solid.
- CadQuery export writes `fillet_r` and `.edges("|Z").fillet(fillet_r)` when the four vertical corners of a rectangular base share one radius, so the radius is a named parameter. An L-bracket fillet stays on the STEP path.
- A fillet that runs into a spherical corner, and a fillet between a plane and a coaxial cylinder (a torus: boss base, rounded hole edge), are not rebuilt.

## 0.4.2 — 2026-10-16
- Narrow chamfers and countersinks that noise breaks into planar strips are gathered back onto the neighbouring cylinder. A chamfered boss stays exact through 0.2 mm noise; a countersink stays exact through 0.05 mm, and on one of three seeds at 0.2 mm.
- When exact sewing fails, the faceted and faceted-hybrid fallbacks are sewn at a tight tolerance and forced onto one parameter range, so the STEP is a valid solid. Coplanar triangles are merged only when that merge stays valid. UnifySameDomain still runs in an isolated process, so a crash keeps the un-merged solid.

## 0.4.1 — 2026-10-16
- Fix: `convert` segfaulted on some parts at the default `--tol` (e.g. a countersunk plate) inside OCCT's UnifySameDomain. Without an explicit `--tol`, tighter fits (0.4×, 0.2×) are now tried before falling back, and the fallback runs face merging in an isolated child process, keeping the un-merged valid solid if OCCT crashes.
- Fix: a fresh `pip install` pulled cadquery-ocp 8.x, whose API breaks the fallback path (`TopoDS.Shell_s`); pinned `cadquery-ocp>=7.8,<8`.
- Regression test for the countersunk plate at default and loose tolerance.

## 0.4.0 — 2026-10-16
- Parametric CadQuery export: `mesh2cad convert in.stl --cadquery out.py` writes an editable script with named parameters (`base_t`, `hole_1_d`, `hole_1_csk_d`, `boss_1_h`, `boss_1_draft_deg`, `boss_1_chamfer`, `pad_1_w`, …). Covers a box/plate base, rectangular pads, through/blind/countersunk holes and straight/drafted/chamfered bosses.
- The generated script is executed and checked against the mesh before it is written (max deviation ≤ max(4·tol, 0.1), volume error ≤ 1 %); unsupported parts exit with code 2 and a reason.
- MCP tool `to_cadquery(input_path, output_path?, tol?)`.
- New optional extra `mesh2cad[cadquery]`; README install now points at PyPI.
- Tests: script reproduces 8 reference parts (clean and 0.05 mm noise) against ground truth, parameter edits change the part as expected, unsupported parts are rejected.
- Benchmark (`benchmarks/run.py`) now includes the cone parts (countersunk plate, tapered boss, chamfered boss) and a cone face column; README table refreshed.
- Docs: Limits and `server.json` description now list cones; Glama score badge and MCP directory links in the README.

## 0.3.0 — 2026-10-15
- Cone primitive: `Cone`, `fit_cone` (constant normal-to-axis angle + tangent planes through the apex), `refine_cone`.
- Segmentation fits cones for chamfers, countersinks and tapered bosses; a coaxial cylinder + cone component is split by normal-to-axis angle before fitting. Cones win over sphere/cylinder fits only when clearly better.
- Exact B-rep conical frustum faces between two planes perpendicular to the axis, with plane ∩ cone circles on the planar faces. Hybrid fallback snaps vertices onto the cone.
- Fidelity report counts `cone` faces. CLI `--no-cones` and MCP `cones` flag (default on).
- Synthetic roundtrip: countersunk plate, 10° tapered boss, chamfered boss.

## 0.2.0 — 2026-10-14
- Sphere primitive: `Sphere`, `fit_sphere` (normal rays + algebraic fit), `refine_sphere`.
- Segmentation emits `Region(surface=Sphere)` for leftover caps and for clusters of small patches that tile a sphere (noisy or faceted ball ends).
- Exact B-rep spherical faces when the region is bounded by plane ∩ sphere circles (cap or zone). Hybrid fallback snaps vertices onto the sphere.
- Fidelity report counts `sphere` faces. CLI `--no-spheres` and MCP `spheres` flag (default on).
- Synthetic roundtrip: hemispherical dome, spherical pocket, and a sphere clipped to a zone between two planes.

## 0.1.0 — 2026-10-14
First release.
- Mesh loading/cleaning for STL, OBJ, PLY (trimesh).
- Segmentation into planes and cylinders (region growing, strip-chain cylinder detection, robust cylinder refinement, region merge/absorb).
- Exact B-rep STEP export via OpenCascade for planar parts with through-holes and bosses; hybrid/faceted fallbacks.
- Fidelity report (Hausdorff, RMS, volume error, face types, validity) as JSON/Markdown.
- CLI: `mesh2cad convert`, `mesh2cad report`, `mesh2cad mcp`.
- MCP server with `convert`, `analyze`, `report` tools.
- Synthetic self-test benchmark (6 parts × 3 noise levels × 3 seeds, all exact) in CI.
