# Changelog

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
