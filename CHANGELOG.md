# Changelog

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
