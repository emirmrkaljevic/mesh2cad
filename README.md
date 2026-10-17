# mesh2cad
[![ci](https://github.com/emirmrkaljevic/mesh2cad/actions/workflows/ci.yml/badge.svg)](https://github.com/emirmrkaljevic/mesh2cad/actions/workflows/ci.yml)
[![emirmrkaljevic/mesh2cad MCP server](https://glama.ai/mcp/servers/emirmrkaljevic/mesh2cad/badges/score.svg)](https://glama.ai/mcp/servers/emirmrkaljevic/mesh2cad)

<p align="center"><a href="https://skintoscan.com/#packages"><img src="https://img.shields.io/badge/%E2%98%95_Support_me-Get_a_SkinToScan_pack-ff5f5f?style=for-the-badge" alt="Support me: get a SkinToScan pack"></a><br><sub>If you find this useful, you can support me by getting a SkinToScan pack, my other project. It's a one-time purchase that directly funds my open-source work.</sub></p>

Turn 3D meshes (STL/OBJ/PLY) into **real, editable CAD (STEP)** — fitted planes, cylinders, spheres, and cones instead of triangle soup. CLI + MCP server for AI agents. MIT, pure Python on top of OpenCascade.

![demo](https://raw.githubusercontent.com/emirmrkaljevic/mesh2cad/main/docs/demo.gif)

![before/after](https://raw.githubusercontent.com/emirmrkaljevic/mesh2cad/main/docs/before_after.png)

*Left: 21k-triangle mesh with 0.1 mm noise. Right: the STEP mesh2cad writes — 10 planes + 2 cylinders, max deviation 0.05 mm from the true part.*

## Install
```
pip install mesh2cad              # CLI + STEP
pip install "mesh2cad[mcp]"       # + MCP server
pip install "mesh2cad[cadquery]"  # + run/verify generated CadQuery scripts
```

## Use
```
mesh2cad convert part.stl -o part.step --tol 0.1 --report report.md
mesh2cad report part.stl --json
```
`--tol` ≈ the noise / accuracy of your mesh in model units (default 0.1 % of the bounding-box diagonal).

Output: STEP + fidelity report — method (`exact` / `faceted-hybrid` / `faceted`), face types, solid validity, two-sided Hausdorff / RMS deviation, volume error, % area per surface type.

### Editable parametric CadQuery script
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
Make the hole 2 mm bigger: set `hole_1_d = 8.0` and run it again.

Supported: Z-up box/plate base, rectangular pads on the top face, vertical holes (through, blind, countersunk) and vertical bosses (straight, drafted, chamfered). Values are rounded to clean numbers (5.997 → 6.0) when that stays within `--tol`. Before writing, the script is executed and compared with the mesh (needs `cadquery` installed); if the part is outside this vocabulary or the script deviates, nothing is written and the command exits with code 2 and the reason. Use STEP output for everything else.

### MCP server (Claude Desktop, Cursor, …)
```json
{ "mcpServers": { "mesh2cad": { "command": "mesh2cad", "args": ["mcp"] } } }
```
Tools: `convert(input_path, output_path?, tol?, cylinders?, spheres?, cones?)`, `analyze(input_path, tol?, cylinders?, spheres?, cones?)`, `report(input_path, tol?, cylinders?, spheres?, cones?)`, `to_cadquery(input_path, output_path?, tol?)` (returns the parameters and measured deviation). `spheres` and `cones` default on.

Listed on [Glama](https://glama.ai/mcp/servers/emirmrkaljevic/mesh2cad) and the [TensorBlock MCP catalog](https://tensorblock.co/mcp/servers/github-emirmrkaljevic-mesh2cad-0f1c0e66).

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
| cube | 0.0 | 3/3 | 6/0 | 0.000 | 0.000 | 0.00 | 0.6 |
| cube | 0.05 | 3/3 | 6/0 | 0.019 | 0.006 | 0.07 | 0.6 |
| cube | 0.2 | 3/3 | 6/0 | 0.075 | 0.026 | 0.34 | 0.6 |
| l_bracket | 0.0 | 3/3 | 8/0 | 0.000 | 0.000 | 0.00 | 0.9 |
| l_bracket | 0.05 | 3/3 | 8/0 | 0.040 | 0.006 | 0.03 | 0.9 |
| l_bracket | 0.2 | 3/3 | 8/0 | 0.216 | 0.026 | 0.15 | 0.8 |
| stepped_block | 0.0 | 3/3 | 11/0 | 0.000 | 0.000 | 0.00 | 2.0 |
| stepped_block | 0.05 | 3/3 | 11/0 | 0.032 | 0.006 | 0.02 | 2.0 |
| stepped_block | 0.2 | 3/3 | 11/0 | 0.136 | 0.022 | 0.07 | 1.9 |
| plate_with_hole | 0.0 | 3/3 | 6/1 | 0.015 | 0.002 | 0.00 | 1.7 |
| plate_with_hole | 0.05 | 3/3 | 6/1 | 0.028 | 0.006 | 0.05 | 1.9 |
| plate_with_hole | 0.2 | 3/3 | 6/1 | 0.078 | 0.023 | 0.18 | 1.6 |
| boss | 0.0 | 3/3 | 7/1 | 0.027 | 0.006 | 0.01 | 1.8 |
| boss | 0.05 | 3/3 | 7/1 | 0.033 | 0.007 | 0.03 | 2.0 |
| boss | 0.2 | 3/3 | 7/1 | 0.090 | 0.017 | 0.07 | 3.3 |
| ribbed_plate | 0.0 | 3/3 | 10/2 | 0.013 | 0.002 | 0.00 | 2.2 |
| ribbed_plate | 0.05 | 3/3 | 10/2 | 0.028 | 0.005 | 0.03 | 2.1 |
| ribbed_plate | 0.2 | 3/3 | 10/2 | 0.108 | 0.020 | 0.12 | 2.1 |
| spherical_dome | 0.0 | 3/3 | 6/0/1 | 0.065 | 0.007 | 0.04 | 1.6 |
| spherical_dome | 0.05 | 3/3 | 6/0/1 | 0.067 | 0.008 | 0.06 | 1.6 |
| spherical_dome | 0.2 | 3/3 | 6/0/1 | 0.241 | 0.028 | 0.25 | 1.7 |
| spherical_pocket | 0.0 | 3/3 | 6/0/1 | 0.046 | 0.004 | 0.02 | 1.9 |
| spherical_pocket | 0.05 | 3/3 | 6/0/1 | 0.057 | 0.005 | 0.03 | 1.9 |
| spherical_pocket | 0.2 | 3/3 | 6/0/1 | 0.075 | 0.015 | 0.13 | 1.8 |
| countersunk_plate | 0.0 | 3/3 | 6/1/0/1 | 0.039 | 0.003 | 0.01 | 2.0 |
| countersunk_plate | 0.05 | 3/3 | 6/1/0/1 | 0.055 | 0.009 | 0.07 | 2.0 |
| countersunk_plate | 0.2 | 1/3 | 16342/0, 6/1/0/1 | 0.590 | 0.049 | 0.54 | 21.6 |
| tapered_boss | 0.0 | 3/3 | 7/0/0/1 | 0.086 | 0.011 | 0.05 | 1.6 |
| tapered_boss | 0.05 | 3/3 | 7/0/0/1 | 0.095 | 0.012 | 0.05 | 1.6 |
| tapered_boss | 0.2 | 3/3 | 7/0/0/1 | 0.123 | 0.025 | 0.22 | 1.5 |
| chamfered_boss | 0.0 | 3/3 | 7/1/0/1 | 0.053 | 0.005 | 0.02 | 2.4 |
| chamfered_boss | 0.05 | 3/3 | 7/1/0/1 | 0.096 | 0.012 | 0.08 | 2.3 |
| chamfered_boss | 0.2 | 3/3 | 7/1/0/1 | 0.511 | 0.059 | 0.49 | 2.1 |

"exact" = valid B-rep with the true face count. All rows are the worst of 3 seeds from `python benchmarks/run.py`. A row below 3/3 includes the faceted fallback (face count = triangles); that fallback is a valid solid. The 90° countersink (Ø6 hole, Ø12 head, 3 mm deep) at 0.2 mm noise is exact on 1 of 3 seeds. CI checks the cases that are exact, plus one fallback that must stay a valid solid.

## Limits (v0.4.2)
Planes, cylinders, spheres, and cones. Cylinders must end on planes perpendicular to their axis (through-holes, bosses, pins). Spheres are exact when each boundary loop is a single plane ∩ sphere circle (domes, spherical bosses, spherical pockets, and zones between two parallel planes). A sphere that runs into an edge or a cylinder falls back to the hybrid solid — vertices are still snapped onto the fitted sphere, but the face stays faceted. Cones are exact when they run between two planes perpendicular to their axis: countersinks, chamfered holes and bosses, drafted/tapered bosses. Cones may also meet a coaxial cylinder (countersunk hole, chamfered pin). A cone that ends in its apex falls back to the hybrid solid. Countersinks and narrow chamfers are exact through about 0.05 mm of noise; a 2 mm chamfer and a wide drafted wall stay exact at 0.2 mm. The countersunk plate's 90° countersink (Ø6 hole, Ø12 head, 3 mm deep) at 0.2 mm noise is exact on one of three seeds. On the other two the conical band and the hole wall meet at a single vertex, so those seeds fall back to a valid faceted solid. No fillets or free-form yet.

## Development
```
pip install -e .[dev]
pytest
python docs/make_demo.py   # regenerate images (needs matplotlib, pillow)
```

License: MIT

<!-- mcp-name: io.github.emirmrkaljevic/mesh2cad -->
