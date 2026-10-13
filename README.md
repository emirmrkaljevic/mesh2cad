# mesh2cad
Convert 3D meshes (STL/OBJ/PLY) into real CAD (STEP) with fitted planes and cylinders — CLI (MCP server for AI agents coming).

> Status: early development (v0.0.x). Not on PyPI yet — install from source.

```
pip install git+https://github.com/emirmrkaljevic/mesh2cad
mesh2cad convert part.stl -o part.step --tol 0.1 --report report.md
mesh2cad report part.stl --json
```

`--tol` should be about the noise level of your mesh (scan/print accuracy) in model units.

## How it works
1. Load + clean mesh (trimesh), Taubin-smoothed copy for normals.
2. Region growing into planes; strip chains / leftovers are fitted as cylinders; tiny or redundant regions are merged/absorbed.
3. Exact B-rep (OCP / OpenCascade): plane polygons from triple-plane intersections, circles from plane ∩ cylinder, cylindrical faces, sewn into a solid.
4. Fallback: hybrid faceted solid (vertices snapped to fitted surfaces, coplanar triangles merged), or raw faceted solid.
5. Fidelity report: two-sided sampled Hausdorff / RMS deviation, volume error, face types, solid validity.

## Self-test benchmark
Ground-truth cadquery parts → tessellated + subdivided (≤2 mm edges) → Gaussian vertex noise → mesh2cad → compared with the truth mesh. 3 seeds, worst case shown (`python benchmarks/run.py`).

| part | noise σ (mm) | exact | faces (plane/cyl) | Hausdorff (mm) | RMS (mm) | vol err % | time s |
|---|---|---|---|---|---|---|---|
| cube | 0.0 | 3/3 | 6/0 | 0.000 | 0.000 | 0.00 | 0.9 |
| cube | 0.05 | 3/3 | 6/0 | 0.019 | 0.006 | 0.07 | 0.8 |
| cube | 0.2 | 3/3 | 6/0 | 0.079 | 0.026 | 0.34 | 0.8 |
| l_bracket | 0.0 | 3/3 | 8/0 | 0.000 | 0.000 | 0.00 | 1.2 |
| l_bracket | 0.05 | 3/3 | 8/0 | 0.039 | 0.006 | 0.03 | 1.1 |
| l_bracket | 0.2 | 3/3 | 8/0 | 0.216 | 0.026 | 0.15 | 1.1 |
| stepped_block | 0.0 | 3/3 | 11/0 | 0.000 | 0.000 | 0.00 | 1.1 |
| stepped_block | 0.05 | 3/3 | 11/0 | 0.032 | 0.006 | 0.02 | 1.5 |
| stepped_block | 0.2 | 3/3 | 11/0 | 0.136 | 0.022 | 0.07 | 1.1 |
| plate_with_hole | 0.0 | 3/3 | 6/1 | 0.015 | 0.002 | 0.00 | 2.5 |
| plate_with_hole | 0.05 | 3/3 | 6/1 | 0.028 | 0.006 | 0.05 | 2.2 |
| plate_with_hole | 0.2 | 3/3 | 6/1 | 0.078 | 0.023 | 0.18 | 2.1 |
| boss | 0.0 | 3/3 | 7/1 | 0.027 | 0.006 | 0.01 | 2.1 |
| boss | 0.05 | 3/3 | 7/1 | 0.031 | 0.007 | 0.04 | 2.0 |
| boss | 0.2 | 3/3 | 7/1 | 0.090 | 0.016 | 0.07 | 2.3 |

"exact" = clean B-rep with the true face count. Known limits: no cones/spheres/fillets yet; cylinders must end on planes perpendicular to their axis.

## Development
```
pip install -e .[dev] cadquery
pytest
```

License: MIT
