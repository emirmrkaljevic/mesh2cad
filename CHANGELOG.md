# Changelog

## 0.1.0 — 2026-10-14
First release.
- Mesh loading/cleaning for STL, OBJ, PLY (trimesh).
- Segmentation into planes and cylinders (region growing, strip-chain cylinder detection, robust cylinder refinement, region merge/absorb).
- Exact B-rep STEP export via OpenCascade for planar parts with through-holes and bosses; hybrid/faceted fallbacks.
- Fidelity report (Hausdorff, RMS, volume error, face types, validity) as JSON/Markdown.
- CLI: `mesh2cad convert`, `mesh2cad report`, `mesh2cad mcp`.
- MCP server with `convert`, `analyze`, `report` tools.
- Synthetic self-test benchmark (6 parts × 3 noise levels × 3 seeds, all exact) in CI.
