# Benchmarks

These rows are STEP fidelity on synthetic parts. A cadquery solid is tessellated, noise is added, mesh2cad converts it, and the result is compared with the true solid. Parametric script checks are not in this table. They are in [CadQuery](cadquery.md), [OpenSCAD](openscad.md), [FreeCAD](freecad.md), and [Limits](limits.md). Real printed parts are in [Real parts](real-parts.md).

Ground-truth cadquery parts → tessellated + subdivided (≤ 2 mm edges) → Gaussian vertex noise → mesh2cad → compared with the truth. 3 seeds, worst case shown (`python benchmarks/run.py`); the same cases run in CI. Numbers are **STEP / B-rep fidelity** (Hausdorff, RMS, volume, face counts). Parametric OpenSCAD, CadQuery, and FreeCAD export are verified separately when you pass `--openscad` / `--cadquery` / `--freecad` (or the MCP `to_openscad` / `to_cadquery` / `to_freecad` tools); they are not rows in this table.

| part | noise σ (mm) | exact | faces (plane/cyl[/sph[/cone[/torus]]]) | radius err (mm) | Hausdorff (mm) | RMS (mm) | vol err % | time s |
|---|---|---|---|---|---|---|---|---|
| cube | 0.0 | 3/3 | 6/0 | - | 0.000 | 0.000 | 0.00 | 0.6 |
| cube | 0.05 | 3/3 | 6/0 | - | 0.019 | 0.006 | 0.07 | 0.6 |
| cube | 0.2 | 3/3 | 6/0 | - | 0.075 | 0.026 | 0.34 | 0.7 |
| l_bracket | 0.0 | 3/3 | 8/0 | - | 0.000 | 0.000 | 0.00 | 0.9 |
| l_bracket | 0.05 | 3/3 | 8/0 | - | 0.040 | 0.006 | 0.03 | 1.1 |
| l_bracket | 0.2 | 3/3 | 8/0 | - | 0.216 | 0.026 | 0.15 | 1.2 |
| stepped_block | 0.0 | 3/3 | 11/0 | - | 0.000 | 0.000 | 0.00 | 1.0 |
| stepped_block | 0.05 | 3/3 | 11/0 | - | 0.032 | 0.006 | 0.02 | 2.2 |
| stepped_block | 0.2 | 3/3 | 11/0 | - | 0.136 | 0.022 | 0.07 | 2.4 |
| plate_with_hole | 0.0 | 3/3 | 6/1 | - | 0.016 | 0.002 | 0.00 | 1.8 |
| plate_with_hole | 0.05 | 3/3 | 6/1 | - | 0.028 | 0.006 | 0.05 | 2.3 |
| plate_with_hole | 0.2 | 3/3 | 6/1 | - | 0.078 | 0.023 | 0.18 | 1.9 |
| boss | 0.0 | 3/3 | 7/1 | - | 0.027 | 0.006 | 0.01 | 1.8 |
| boss | 0.05 | 3/3 | 7/1 | - | 0.033 | 0.007 | 0.03 | 2.2 |
| boss | 0.2 | 3/3 | 7/1 | - | 0.091 | 0.017 | 0.09 | 4.6 |
| ribbed_plate | 0.0 | 3/3 | 10/2 | - | 0.012 | 0.001 | 0.00 | 2.5 |
| ribbed_plate | 0.05 | 3/3 | 10/2 | - | 0.028 | 0.005 | 0.03 | 2.5 |
| ribbed_plate | 0.2 | 3/3 | 10/2 | - | 0.108 | 0.020 | 0.12 | 2.6 |
| spherical_dome | 0.0 | 3/3 | 6/0/1 | - | 0.065 | 0.007 | 0.04 | 1.8 |
| spherical_dome | 0.05 | 3/3 | 6/0/1 | - | 0.067 | 0.008 | 0.06 | 1.8 |
| spherical_dome | 0.2 | 3/3 | 6/0/1 | - | 0.241 | 0.028 | 0.25 | 1.9 |
| spherical_pocket | 0.0 | 3/3 | 6/0/1 | - | 0.046 | 0.004 | 0.02 | 2.2 |
| spherical_pocket | 0.05 | 3/3 | 6/0/1 | - | 0.057 | 0.005 | 0.03 | 2.1 |
| spherical_pocket | 0.2 | 3/3 | 6/0/1 | - | 0.075 | 0.015 | 0.13 | 2.1 |
| countersunk_plate | 0.0 | 3/3 | 6/1/0/1 | - | 0.039 | 0.003 | 0.01 | 2.2 |
| countersunk_plate | 0.05 | 3/3 | 6/1/0/1 | - | 0.055 | 0.009 | 0.07 | 2.2 |
| countersunk_plate | 0.2 | 1/3 | 16342/0, 6/1/0/1 | - | 0.590 | 0.049 | 0.54 | 21.2 |
| tapered_boss | 0.0 | 3/3 | 7/0/0/1 | - | 0.086 | 0.011 | 0.05 | 1.9 |
| tapered_boss | 0.05 | 3/3 | 7/0/0/1 | - | 0.095 | 0.012 | 0.05 | 1.8 |
| tapered_boss | 0.2 | 3/3 | 7/0/0/1 | - | 0.123 | 0.025 | 0.22 | 2.1 |
| chamfered_boss | 0.0 | 3/3 | 7/1/0/1 | - | 0.053 | 0.005 | 0.02 | 2.7 |
| chamfered_boss | 0.05 | 3/3 | 7/1/0/1 | - | 0.096 | 0.012 | 0.08 | 2.5 |
| chamfered_boss | 0.2 | 3/3 | 7/1/0/1 | - | 0.511 | 0.059 | 0.49 | 2.3 |
| filleted_plate | 0.0 | 3/3 | 6/4 | 0.001 | 0.013 | 0.002 | 0.01 | 1.3 |
| filleted_plate | 0.05 | 3/3 | 6/4 | 0.100 | 0.038 | 0.008 | 0.07 | 1.3 |
| filleted_plate | 0.2 | 2/3 | 5796/0, 6/4, 7/3 | 0.447 | 0.857 | 0.086 | 0.33 | 6.8 |
| filleted_box | 0.0 | 3/3 | 6/4 | 0.003 | 0.006 | 0.001 | 0.00 | 1.8 |
| filleted_box | 0.05 | 3/3 | 6/4 | 0.041 | 0.028 | 0.006 | 0.03 | 1.9 |
| filleted_box | 0.2 | 3/3 | 6/3, 6/4 | 0.206 | 0.845 | 0.048 | 0.18 | 1.8 |
| filleted_bracket | 0.0 | 3/3 | 8/6 | 0.005 | 0.010 | 0.002 | 0.00 | 1.8 |
| filleted_bracket | 0.05 | 3/3 | 8/6 | 0.131 | 0.041 | 0.008 | 0.02 | 1.9 |
| filleted_bracket | 0.2 | 0/3 | 10452/0 | - | 0.958 | 0.130 | 0.50 | 13.7 |
| rounded_box | 0.0 | 3/3 | 6/12/8 | 0.006 | 0.023 | 0.003 | 0.01 | 12.6 |
| rounded_box | 0.05 | 3/3 | 6/12/8 | 0.007 | 0.035 | 0.007 | 0.07 | 8.4 |
| rounded_box | 0.2 | 3/3 | 6/12/8 | 0.112 | 0.116 | 0.026 | 0.18 | 7.4 |
| filleted_boss | 0.0 | 3/3 | 7/1/0/0/1 | 0.008 | 0.026 | 0.002 | 0.05 | 3.2 |
| filleted_boss | 0.05 | 3/3 | 7/1/0/0/1 | 0.279 | 0.150 | 0.021 | 0.10 | 3.2 |
| filleted_boss | 0.2 | 3/3 | 7/1 | - | 0.864 | 0.068 | 0.25 | 3.3 |
| rounded_hole | 0.0 | 3/3 | 6/1/0/0/1 | 0.068 | 0.026 | 0.004 | 0.05 | 3.9 |
| rounded_hole | 0.05 | 2/3 | 20908/0, 6/1/0/0/1 | 0.247 | 0.633 | 0.036 | 0.73 | 29.7 |
| rounded_hole | 0.2 | 3/3 | 6/1 | - | 0.681 | 0.041 | 0.09 | 2.8 |

"exact" = valid analytic B-rep. A face count equal to the triangle count is the faceted fallback, which is still a valid solid. When a sphere is present the face column is plane/cylinder/sphere, and when a torus is present it is plane/cylinder/sphere/cone/torus. radius err is the worst |fitted − true| cylinder radius among seeds that came back exact, or the torus minor radius for filleted_boss and rounded_hole, and - when the part has no fillet radius or no seed was exact. All rows are the worst of 3 seeds from `python benchmarks/run.py`. The 90° countersink (Ø6 hole, Ø12 head, 3 mm deep) at 0.2 mm noise is exact on 1 of 3 seeds. At 0.2 mm the filleted bracket falls back on every seed; the plate is exact on 2 of 3, and one of those solids is missing a corner (7 planes + 3 cylinders). The block stays a valid solid, but one seed keeps only three of the four rounds (Hausdorff 0.85 mm). At 0.2 mm the boss fillet and the hole-rim round are not recovered: the STEP is the sharp cylinder and planes, which is why the Hausdorff is the size of the missing round. The rounded box (40×30×20 mm, every edge at r=3 mm) is exact on 3 of 3 seeds at each noise level: 6 planes, 12 cylinders and 8 sphere octants. With no noise the radius error is 0.006 mm and the Hausdorff distance is 0.023 mm. At 0.05 mm the worst radius error is 0.007 mm (Hausdorff 0.035 mm). At 0.2 mm the worst radius error is 0.112 mm (Hausdorff 0.116 mm). The eight spheres share that cylinder radius. CI checks the cases that are exact, plus fallbacks that must stay valid solids.

Back to the [README](../README.md).
