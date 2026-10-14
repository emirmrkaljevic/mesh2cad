# mesh2cad
[![ci](https://github.com/emirmrkaljevic/mesh2cad/actions/workflows/ci.yml/badge.svg)](https://github.com/emirmrkaljevic/mesh2cad/actions/workflows/ci.yml)

Turn 3D meshes (STL/OBJ/PLY) into **real, editable CAD (STEP)** — fitted planes, cylinders, and spheres instead of triangle soup. CLI + MCP server for AI agents. MIT, pure Python on top of OpenCascade.

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
Tools: `convert(input_path, output_path?, tol?, cylinders?, spheres?)`, `analyze(input_path, tol?, cylinders?, spheres?)`, `report(input_path, tol?, cylinders?, spheres?)`. `spheres` defaults on.

## How it works
1. Load + clean mesh (trimesh); Taubin-smoothed copy for normals.
2. Region growing into planes; strip chains and leftovers are fitted as cylinders or spheres (normal-ray / algebraic fit, robust nonlinear refinement). Small noisy patches that tile a sphere are merged into one spherical region. Split/redundant regions are merged and absorbed.
3. Exact B-rep (OCP / OpenCascade): plane polygons from triple-plane intersections, circles from plane ∩ cylinder and plane ∩ sphere, cylindrical faces, spherical caps and zones, sewn into a solid and validated.
4. Fallbacks: hybrid faceted solid (vertices snapped to fitted surfaces, coplanar triangles merged), then raw faceted solid.
5. Fidelity report against the input mesh.

## Self-test benchmark
Ground-truth cadquery parts → tessellated + subdivided (≤ 2 mm edges) → Gaussian vertex noise → mesh2cad → compared with the truth. 3 seeds, worst case shown (`python benchmarks/run.py`); the same cases run in CI.

| part | noise σ (mm) | exact | faces (plane/cyl[/sph]) | Hausdorff (mm) | RMS (mm) | vol err % | time s |
|---|---|---|---|---|---|---|---|
| cube | 0.0 | 3/3 | 6/0 | 0.000 | 0.000 | 0.00 | 0.9 |
| cube | 0.05 | 3/3 | 6/0 | 0.019 | 0.006 | 0.07 | 0.8 |
| cube | 0.2 | 3/3 | 6/0 | 0.075 | 0.026 | 0.34 | 0.9 |
| l_bracket | 0.0 | 3/3 | 8/0 | 0.000 | 0.000 | 0.00 | 1.3 |
| l_bracket | 0.05 | 3/3 | 8/0 | 0.040 | 0.006 | 0.03 | 1.2 |
| l_bracket | 0.2 | 3/3 | 8/0 | 0.216 | 0.026 | 0.15 | 1.2 |
| stepped_block | 0.0 | 3/3 | 11/0 | 0.000 | 0.000 | 0.00 | 1.2 |
| stepped_block | 0.05 | 3/3 | 11/0 | 0.032 | 0.006 | 0.02 | 1.6 |
| stepped_block | 0.2 | 3/3 | 11/0 | 0.136 | 0.022 | 0.07 | 1.6 |
| plate_with_hole | 0.0 | 3/3 | 6/1 | 0.015 | 0.002 | 0.00 | 2.4 |
| plate_with_hole | 0.05 | 3/3 | 6/1 | 0.028 | 0.006 | 0.05 | 2.3 |
| plate_with_hole | 0.2 | 3/3 | 6/1 | 0.078 | 0.023 | 0.18 | 2.7 |
| boss | 0.0 | 3/3 | 7/1 | 0.027 | 0.006 | 0.01 | 2.1 |
| boss | 0.05 | 3/3 | 7/1 | 0.031 | 0.007 | 0.04 | 2.7 |
| boss | 0.2 | 3/3 | 7/1 | 0.090 | 0.016 | 0.07 | 2.5 |
| ribbed_plate | 0.0 | 3/3 | 10/2 | 0.013 | 0.002 | 0.00 | 2.9 |
| ribbed_plate | 0.05 | 3/3 | 10/2 | 0.028 | 0.005 | 0.03 | 3.1 |
| ribbed_plate | 0.2 | 3/3 | 10/2 | 0.108 | 0.020 | 0.12 | 3.3 |
| spherical_dome | 0.0 | 3/3 | 6/0/1 | 0.065 | 0.007 | 0.04 | 1.7 |
| spherical_dome | 0.05 | 3/3 | 6/0/1 | 0.067 | 0.008 | 0.06 | 1.7 |
| spherical_dome | 0.2 | 3/3 | 6/0/1 | 0.241 | 0.028 | 0.25 | 1.6 |
| spherical_pocket | 0.0 | 3/3 | 6/0/1 | 0.046 | 0.004 | 0.02 | 2.0 |
| spherical_pocket | 0.05 | 3/3 | 6/0/1 | 0.057 | 0.005 | 0.03 | 1.9 |
| spherical_pocket | 0.2 | 3/3 | 6/0/1 | 0.075 | 0.015 | 0.13 | 1.9 |

"exact" = valid B-rep with the true face count. Plane and cylinder rows are the 0.1.0 measurements (accuracy unchanged). Sphere rows are the worst of 3 seeds from `python benchmarks/run.py`.

## Limits (v0.2)
Planes, cylinders, and spheres. Cylinders must end on planes perpendicular to their axis (through-holes, bosses, pins). Spheres are exact when each boundary loop is a single plane ∩ sphere circle (domes, spherical bosses, spherical pockets, and zones between two parallel planes). A sphere that runs into an edge or a cylinder falls back to the hybrid solid — vertices are still snapped onto the fitted sphere, but the face stays faceted. No cones, fillets, or free-form yet.

## Development
```
pip install -e .[dev]
pytest
python docs/make_demo.py   # regenerate images (needs matplotlib, pillow)
```

License: MIT

<!-- mcp-name: io.github.emirmrkaljevic/mesh2cad -->
