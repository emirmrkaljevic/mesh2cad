# mesh2cad
[![ci](https://github.com/emirmrkaljevic/mesh2cad/actions/workflows/ci.yml/badge.svg)](https://github.com/emirmrkaljevic/mesh2cad/actions/workflows/ci.yml)

Turn 3D meshes (STL/OBJ/PLY) into **real, editable CAD (STEP)** — fitted planes, cylinders, spheres, and cones instead of triangle soup. CLI + MCP server for AI agents. MIT, pure Python on top of OpenCascade.

![demo](https://raw.githubusercontent.com/emirmrkaljevic/mesh2cad/main/docs/demo.gif)

![before/after](https://raw.githubusercontent.com/emirmrkaljevic/mesh2cad/main/docs/before_after.png)

*Left: 21k-triangle mesh with 0.1 mm noise. Right: the STEP mesh2cad writes — 10 planes + 2 cylinders, max deviation 0.05 mm from the true part.*

## Install
```
pip install "mesh2cad[mcp] @ git+https://github.com/emirmrkaljevic/mesh2cad"
```
(PyPI release coming.)

## Use
```
mesh2cad convert part.stl -o part.step --tol 0.1 --report report.md
mesh2cad report part.stl --json
```
`--tol` ≈ the noise / accuracy of your mesh in model units (default 0.1 % of the bounding-box diagonal).

Output: STEP + fidelity report — method (`exact` / `faceted-hybrid` / `faceted`), face types, solid validity, two-sided Hausdorff / RMS deviation, volume error, % area per surface type.

### MCP server (Claude Desktop, Cursor, …)
```json
{ "mcpServers": { "mesh2cad": { "command": "mesh2cad", "args": ["mcp"] } } }
```
Tools: `convert(input_path, output_path?, tol?, cylinders?, spheres?, cones?)`, `analyze(input_path, tol?, cylinders?, spheres?, cones?)`, `report(input_path, tol?, cylinders?, spheres?, cones?)`. `spheres` and `cones` default on.

## How it works
1. Load + clean mesh (trimesh); Taubin-smoothed copy for normals.
2. Region growing into planes; strip chains and leftovers are fitted as cylinders, spheres, or cones (normal-ray / algebraic / constant-normal-angle fits, robust nonlinear refinement). A coaxial cylinder + cone component (countersink, chamfered pin) is split by normal-to-axis angle first. Small noisy patches that tile a sphere are merged into one spherical region. Split/redundant regions are merged and absorbed.
3. Exact B-rep (OCP / OpenCascade): plane polygons from triple-plane intersections, circles from plane ∩ cylinder and plane ∩ sphere, plane ∩ cone, cylindrical faces, spherical caps and zones, conical frustums, sewn into a solid and validated.
4. Fallbacks: hybrid faceted solid (vertices snapped to fitted surfaces, coplanar triangles merged), then raw faceted solid.
5. Fidelity report against the input mesh.

## Self-test benchmark
Ground-truth cadquery parts → tessellated + subdivided (≤ 2 mm edges) → Gaussian vertex noise → mesh2cad → compared with the truth. 3 seeds, worst case shown (`python benchmarks/run.py`); the same cases run in CI.

| part | noise σ (mm) | exact | faces (plane/cyl[/sph[/cone]]) | Hausdorff (mm) | RMS (mm) | vol err % | time s |
|---|---|---|---|---|---|---|---|
| cube | 0.0 | 3/3 | 6/0 | 0.000 | 0.000 | 0.00 | 0.7 |
| cube | 0.05 | 3/3 | 6/0 | 0.019 | 0.006 | 0.07 | 0.6 |
| cube | 0.2 | 3/3 | 6/0 | 0.075 | 0.026 | 0.34 | 0.6 |
| l_bracket | 0.0 | 3/3 | 8/0 | 0.000 | 0.000 | 0.00 | 0.9 |
| l_bracket | 0.05 | 3/3 | 8/0 | 0.040 | 0.006 | 0.03 | 1.0 |
| l_bracket | 0.2 | 3/3 | 8/0 | 0.216 | 0.026 | 0.15 | 1.0 |
| stepped_block | 0.0 | 3/3 | 11/0 | 0.000 | 0.000 | 0.00 | 2.1 |
| stepped_block | 0.05 | 3/3 | 11/0 | 0.032 | 0.006 | 0.02 | 2.3 |
| stepped_block | 0.2 | 3/3 | 11/0 | 0.136 | 0.022 | 0.07 | 2.1 |
| plate_with_hole | 0.0 | 3/3 | 6/1 | 0.015 | 0.002 | 0.00 | 1.9 |
| plate_with_hole | 0.05 | 3/3 | 6/1 | 0.028 | 0.006 | 0.05 | 1.9 |
| plate_with_hole | 0.2 | 3/3 | 6/1 | 0.078 | 0.023 | 0.18 | 1.9 |
| boss | 0.0 | 3/3 | 7/1 | 0.027 | 0.006 | 0.01 | 1.9 |
| boss | 0.05 | 3/3 | 7/1 | 0.031 | 0.007 | 0.04 | 2.5 |
| boss | 0.2 | 3/3 | 7/1 | 0.090 | 0.016 | 0.07 | 3.5 |
| ribbed_plate | 0.0 | 3/3 | 10/2 | 0.013 | 0.002 | 0.00 | 2.3 |
| ribbed_plate | 0.05 | 3/3 | 10/2 | 0.028 | 0.005 | 0.03 | 2.4 |
| ribbed_plate | 0.2 | 3/3 | 10/2 | 0.108 | 0.020 | 0.12 | 2.3 |
| spherical_dome | 0.0 | 3/3 | 6/0/1 | 0.065 | 0.007 | 0.04 | 1.9 |
| spherical_dome | 0.05 | 3/3 | 6/0/1 | 0.067 | 0.008 | 0.06 | 2.0 |
| spherical_dome | 0.2 | 3/3 | 6/0/1 | 0.241 | 0.028 | 0.25 | 1.9 |
| spherical_pocket | 0.0 | 3/3 | 6/0/1 | 0.046 | 0.004 | 0.02 | 2.2 |
| spherical_pocket | 0.05 | 3/3 | 6/0/1 | 0.057 | 0.005 | 0.03 | 2.1 |
| spherical_pocket | 0.2 | 3/3 | 6/0/1 | 0.075 | 0.015 | 0.13 | 2.1 |
| countersunk_plate | 0.0 | 3/3 | 6/1/0/1 | 0.039 | 0.003 | 0.01 | 2.2 |
| countersunk_plate | 0.05 | 3/3 | 6/1/0/1 | 0.055 | 0.009 | 0.07 | 2.1 |
| countersunk_plate | 0.2 | 0/3 | 16322/0, 16324/0 | 0.745 | 0.130 | 0.17 | 41.6 |
| tapered_boss | 0.0 | 3/3 | 7/0/0/1 | 0.086 | 0.011 | 0.05 | 1.8 |
| tapered_boss | 0.05 | 3/3 | 7/0/0/1 | 0.095 | 0.012 | 0.05 | 1.7 |
| tapered_boss | 0.2 | 3/3 | 7/0/0/1 | 0.123 | 0.025 | 0.22 | 1.7 |
| chamfered_boss | 0.0 | 3/3 | 7/1/0/1 | 0.053 | 0.005 | 0.02 | 2.9 |
| chamfered_boss | 0.05 | 1/3 | 16461/0, 16473/0, 7/1/0/1 | 0.160 | 0.035 | 0.08 | 31.6 |
| chamfered_boss | 0.2 | 0/3 | 16565/0, 16569/0, 16570/0 | 0.703 | 0.133 | 0.26 | 34.8 |

"exact" = valid B-rep with the true face count. All rows are the worst of 3 seeds from `python benchmarks/run.py` (v0.3.0). Rows with fewer than 3/3 exact fell back to a faceted solid (face count = triangles): noisy narrow countersink/chamfer bands don't constrain a cone yet, so CI checks those parts only where they are exact.

## Limits (v0.3)
Planes, cylinders, and spheres. Cylinders must end on planes perpendicular to their axis (through-holes, bosses, pins). Spheres are exact when each boundary loop is a single plane ∩ sphere circle (domes, spherical bosses, spherical pockets, and zones between two parallel planes). A sphere that runs into an edge or a cylinder falls back to the hybrid solid — vertices are still snapped onto the fitted sphere, but the face stays faceted. Cones are exact when they run between two planes perpendicular to their axis: countersinks, chamfered holes and bosses, drafted/tapered bosses. Cones may also meet a coaxial cylinder (countersunk hole, chamfered pin). A cone that ends in its apex falls back to the hybrid solid, and countersinks are reliable up to ~0.05 mm noise, narrow chamfers on clean meshes only for now (wide drafted walls up to 0.2 mm). No fillets or free-form yet.

## Development
```
pip install -e .[dev]
pytest
python docs/make_demo.py   # regenerate images (needs matplotlib, pillow)
```

License: MIT

<!-- mcp-name: io.github.emirmrkaljevic/mesh2cad -->
