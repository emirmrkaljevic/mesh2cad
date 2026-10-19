# mesh2cad
[![ci](https://github.com/emirmrkaljevic/mesh2cad/actions/workflows/ci.yml/badge.svg)](https://github.com/emirmrkaljevic/mesh2cad/actions/workflows/ci.yml)
[![emirmrkaljevic/mesh2cad MCP server](https://glama.ai/mcp/servers/emirmrkaljevic/mesh2cad/badges/score.svg)](https://glama.ai/mcp/servers/emirmrkaljevic/mesh2cad)

<p align="center"><a href="https://skintoscan.com/#packages"><img src=".github/assets/skintoscan-support.svg" width="684" alt="Support my open-source work: get a SkinToScan pack"></a><br><sub>If you find this useful, you can support me by getting a SkinToScan pack, my other project. It's a one-time purchase that directly funds my open-source work.</sub></p>
<p align="center"><a href="https://skintoscan.com/checkout/?add-to-cart=138"><img src=".github/assets/skintoscan-pack-1.svg" width="220" alt="Buy Single Picture Check, €4.49"></a> <a href="https://skintoscan.com/checkout/?add-to-cart=139"><img src=".github/assets/skintoscan-pack-5.svg" width="220" alt="Buy 5 Picture Pack, €19.99 (popular)"></a> <a href="https://skintoscan.com/checkout/?add-to-cart=140"><img src=".github/assets/skintoscan-pack-15.svg" width="220" alt="Buy 15 Picture Pack, €54.99"></a></p>

Turn 3D meshes (STL/OBJ/PLY) into **real, editable CAD (STEP)** — fitted planes, cylinders, spheres, cones, and constant-radius fillets instead of triangle soup. CLI + MCP server for AI agents. MIT, pure Python on top of OpenCascade.

![demo](https://raw.githubusercontent.com/emirmrkaljevic/mesh2cad/main/docs/demo.gif)

![before/after](https://raw.githubusercontent.com/emirmrkaljevic/mesh2cad/main/docs/before_after.png)

*Left: 21k-triangle mesh with 0.1 mm noise. Right: the STEP mesh2cad writes — 10 planes + 2 cylinders, max deviation 0.05 mm from the true part.*

## Install
```
pip install mesh2cad              # CLI + STEP
pip install "mesh2cad[mcp]"       # + MCP server
pip install "mesh2cad[cadquery]"  # + run/verify generated CadQuery scripts
```
OpenSCAD scripts are plain text. Verifying one needs the `openscad` program ([openscad.org](https://openscad.org)); it is not a Python dependency. Without it the `.scad` file is still written.

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
Make the hole 2 mm bigger: set `hole_1_d = 8.0` and run it again. A plate or block whose four vertical corners share one radius also gets `fillet_r` and `.edges("|Z").fillet(fillet_r)`, so the same edit grows the round.

A hole or boss on a side face uses the same kind of parameter. This is the block below (60×40×36, an Ø8 hole through X, an Ø6 hole 12 deep opening on +Y, an Ø10 boss 8 tall on +X):

```python
boss_1_d     = 10.0  # diameter
boss_1_h     = 8.0  # height from the +X face
boss_1_y     = 10.0
boss_1_z     = 8.0
hole_1_d     = 8.0  # diameter (axis X)
hole_1_y     = 0.0
hole_1_z     = -6.0
hole_2_d     = 6.0  # diameter (axis Y)
hole_2_x     = -12.0
hole_2_z     = 8.0
hole_2_depth = 12.0  # blind depth from the +Y face
```

Set `hole_1_d = 10.0` and run the script again. The volume drops by π(5²−4²)·60 mm³, the analytic change for that through hole.

Supported: a Z-up box/plate base, rectangular pads on the top face, holes along X, Y, or Z (through, blind, and countersunk when the axis is vertical), bosses along X, Y, or Z (straight, and drafted or chamfered when the axis is vertical), and that one vertical corner radius. Values are rounded to clean numbers (5.997 → 6.0) when that stays within `--tol`. Before writing, the script is executed and compared with the mesh (needs `cadquery` installed); if the part is outside this vocabulary or the script deviates, nothing is written and the command exits with code 2 and the reason. Use STEP output for everything else.

### Editable parametric OpenSCAD script
The same parts can be written as OpenSCAD, which is what the 3D-printing Customizer on Thingiverse and Printables edits. Names and values match the CadQuery script. `--openscad` and `--cadquery` can be combined, and `-o` still writes STEP.
```
mesh2cad convert plate.stl --openscad plate.scad
openscad -o plate.stl plate.scad
```
A 40×30×6 plate with a Ø10 hole comes out as below, plus a short header. Before:

```openscad
// --- parameters ---
/* [Base] */
// base length (X)
base_x   = 40.0;
// base width (Y)
base_y   = 30.0;
// base thickness (Z)
base_t   = 6.0;
base_cx  = 0.0;
base_cy  = 0.0;
// bottom face height
base_z0  = -3.0;

/* [Holes] */
// diameter
hole_1_d = 10.0; // [0.2:0.1:30]
hole_1_x = 0.0;
hole_1_y = 0.0;

/* [Hidden] */
$fn = 96; // circle segments

// --- model ---
module mesh2cad_part() {
  difference() {
    translate([base_cx, base_cy, base_z0 + base_t / 2])
      cube([base_x, base_y, base_t], center=true);
    translate([hole_1_x, hole_1_y, base_z0 - 1])
      cylinder(h=base_t + 2, d=hole_1_d);
  }
}

mesh2cad_part();
```

After, the hole is 2 mm bigger — one number changes, then re-render:

```openscad
hole_1_d = 12.0; // [0.2:0.1:30]
```

Re-render (F6) and export STL. A drafted boss is a `cylinder` with `r1`/`r2`, a chamfer or countersink is the same, and the four vertical corners of a rectangular base are `offset(r=fillet_r)` on a smaller square, then extruded. A hole or boss whose axis is X or Y is the same `cylinder`, turned with `rotate` so it points along that axis; its parameters are the diameter and the two coordinates of the axis (`hole_1_y` / `hole_1_z` for an X hole). A slider (`// [min:max]`) is added only where the part gives a bound: hole diameter up to the shorter side of the face it is cut through, fillet up to half the short side of the plate, draft until the tip would close, countersink angle in 60–120°, chamfer up to the boss, pad size up to the base, blind depth up to the length of the hole's axis. Positions (`hole_1_x`, `base_z0`, …) stay plain numbers.

Same honesty as CadQuery. Before writing, mesh2cad renders the script with the `openscad` CLI and compares the STL to the mesh (sampled Hausdorff and volume, same limits as the CadQuery check). If `openscad` is not installed the file is still written and the command says it was unverified. If the part is outside this vocabulary, or the render does not match, nothing is written and the command exits with code 2 and the reason. Circles are `$fn = 96` polygons; on these parts that is a few thousandths of a millimetre, not an analytic solid. STEP remains the path for everything else.

### MCP server (Claude Desktop, Cursor, …)
```json
{ "mcpServers": { "mesh2cad": { "command": "mesh2cad", "args": ["mcp"] } } }
```
Tools: `convert(input_path, output_path?, tol?, cylinders?, spheres?, cones?)`, `analyze(input_path, tol?, cylinders?, spheres?, cones?)`, `report(input_path, tol?, cylinders?, spheres?, cones?)`, `to_cadquery(input_path, output_path?, tol?)`, `to_openscad(input_path, output_path?, tol?)`. Both script tools return the parameters and, when the script was run or rendered, the measured deviation. `to_openscad` also returns `verified: false` when the `openscad` program is absent. `spheres` and `cones` default on.

Listed on [Glama](https://glama.ai/mcp/servers/emirmrkaljevic/mesh2cad) and the [TensorBlock MCP catalog](https://tensorblock.co/mcp/servers/github-emirmrkaljevic-mesh2cad-0f1c0e66).

## How it works
1. Load + clean mesh (trimesh); Taubin-smoothed copy for normals.
2. Region growing into planes; strip chains and leftovers are fitted as cylinders, spheres, or cones (normal-ray / algebraic / constant-normal-angle fits, robust nonlinear refinement). A narrow round between two planes is gathered back into one partial cylinder and snapped so its generators lie in those planes. A coaxial cylinder + cone component (countersink, chamfered pin) is split by normal-to-axis angle first. Small noisy patches that tile a sphere are merged into one spherical region. Split/redundant regions are merged and absorbed.
3. Exact B-rep (OCP / OpenCascade): plane polygons from triple-plane intersections, circles from plane ∩ cylinder and plane ∩ sphere, plane ∩ cone, full cylindrical walls, partial fillet cylinders, spherical caps and zones, conical frustums, sewn into a solid and validated.
4. Fallbacks: hybrid faceted solid (vertices snapped to fitted surfaces, coplanar triangles merged), then raw faceted solid.
5. Fidelity report against the input mesh.

## Self-test benchmark
Ground-truth cadquery parts → tessellated + subdivided (≤ 2 mm edges) → Gaussian vertex noise → mesh2cad → compared with the truth. 3 seeds, worst case shown (`python benchmarks/run.py`); the same cases run in CI. Numbers are **STEP / B-rep fidelity** (Hausdorff, RMS, volume, face counts). Parametric OpenSCAD and CadQuery export are verified separately when you pass `--openscad` / `--cadquery` (or the MCP `to_openscad` / `to_cadquery` tools); they are not rows in this table.

| part | noise σ (mm) | exact | faces (plane/cyl[/sph[/cone]]) | radius err (mm) | Hausdorff (mm) | RMS (mm) | vol err % | time s |
|---|---|---|---|---|---|---|---|---|
| cube | 0.0 | 3/3 | 6/0 | — | 0.000 | 0.000 | 0.00 | 0.6 |
| cube | 0.05 | 3/3 | 6/0 | — | 0.019 | 0.006 | 0.07 | 0.6 |
| cube | 0.2 | 3/3 | 6/0 | — | 0.075 | 0.026 | 0.34 | 0.7 |
| l_bracket | 0.0 | 3/3 | 8/0 | — | 0.000 | 0.000 | 0.00 | 0.9 |
| l_bracket | 0.05 | 3/3 | 8/0 | — | 0.040 | 0.006 | 0.03 | 1.1 |
| l_bracket | 0.2 | 3/3 | 8/0 | — | 0.216 | 0.026 | 0.15 | 1.2 |
| stepped_block | 0.0 | 3/3 | 11/0 | — | 0.000 | 0.000 | 0.00 | 2.2 |
| stepped_block | 0.05 | 3/3 | 11/0 | — | 0.032 | 0.006 | 0.02 | 2.2 |
| stepped_block | 0.2 | 3/3 | 11/0 | — | 0.136 | 0.022 | 0.07 | 2.4 |
| plate_with_hole | 0.0 | 3/3 | 6/1 | — | 0.015 | 0.002 | 0.00 | 1.8 |
| plate_with_hole | 0.05 | 3/3 | 6/1 | — | 0.028 | 0.006 | 0.05 | 2.3 |
| plate_with_hole | 0.2 | 3/3 | 6/1 | — | 0.078 | 0.023 | 0.18 | 1.9 |
| boss | 0.0 | 3/3 | 7/1 | — | 0.027 | 0.006 | 0.01 | 1.8 |
| boss | 0.05 | 3/3 | 7/1 | — | 0.033 | 0.007 | 0.03 | 2.2 |
| boss | 0.2 | 3/3 | 7/1 | — | 0.091 | 0.017 | 0.09 | 4.6 |
| ribbed_plate | 0.0 | 3/3 | 10/2 | — | 0.013 | 0.002 | 0.00 | 2.5 |
| ribbed_plate | 0.05 | 3/3 | 10/2 | — | 0.028 | 0.005 | 0.03 | 2.5 |
| ribbed_plate | 0.2 | 3/3 | 10/2 | — | 0.108 | 0.020 | 0.12 | 2.6 |
| spherical_dome | 0.0 | 3/3 | 6/0/1 | — | 0.065 | 0.007 | 0.04 | 1.8 |
| spherical_dome | 0.05 | 3/3 | 6/0/1 | — | 0.067 | 0.008 | 0.06 | 1.8 |
| spherical_dome | 0.2 | 3/3 | 6/0/1 | — | 0.241 | 0.028 | 0.25 | 1.9 |
| spherical_pocket | 0.0 | 3/3 | 6/0/1 | — | 0.046 | 0.004 | 0.02 | 2.2 |
| spherical_pocket | 0.05 | 3/3 | 6/0/1 | — | 0.057 | 0.005 | 0.03 | 2.1 |
| spherical_pocket | 0.2 | 3/3 | 6/0/1 | — | 0.075 | 0.015 | 0.13 | 2.1 |
| countersunk_plate | 0.0 | 3/3 | 6/1/0/1 | — | 0.039 | 0.003 | 0.01 | 2.2 |
| countersunk_plate | 0.05 | 3/3 | 6/1/0/1 | — | 0.055 | 0.009 | 0.07 | 2.2 |
| countersunk_plate | 0.2 | 1/3 | 16342/0, 6/1/0/1 | — | 0.590 | 0.049 | 0.54 | 21.2 |
| tapered_boss | 0.0 | 3/3 | 7/0/0/1 | — | 0.086 | 0.011 | 0.05 | 1.9 |
| tapered_boss | 0.05 | 3/3 | 7/0/0/1 | — | 0.095 | 0.012 | 0.05 | 1.8 |
| tapered_boss | 0.2 | 3/3 | 7/0/0/1 | — | 0.123 | 0.025 | 0.22 | 2.1 |
| chamfered_boss | 0.0 | 3/3 | 7/1/0/1 | — | 0.053 | 0.005 | 0.02 | 2.7 |
| chamfered_boss | 0.05 | 3/3 | 7/1/0/1 | — | 0.096 | 0.012 | 0.08 | 2.5 |
| chamfered_boss | 0.2 | 3/3 | 7/1/0/1 | — | 0.511 | 0.059 | 0.49 | 2.3 |
| filleted_plate | 0.0 | 3/3 | 6/4 | 0.009 | 0.017 | 0.002 | 0.00 | 1.3 |
| filleted_plate | 0.05 | 3/3 | 6/4 | 0.100 | 0.038 | 0.008 | 0.07 | 1.3 |
| filleted_plate | 0.2 | 2/3 | 5796/0, 6/4, 7/3 | 0.447 | 0.857 | 0.086 | 0.33 | 6.8 |
| filleted_box | 0.0 | 3/3 | 6/4 | 0.005 | 0.007 | 0.001 | 0.00 | 1.8 |
| filleted_box | 0.05 | 3/3 | 6/4 | 0.041 | 0.028 | 0.006 | 0.03 | 1.9 |
| filleted_box | 0.2 | 3/3 | 6/3, 6/4 | 0.206 | 0.845 | 0.048 | 0.18 | 1.8 |
| filleted_bracket | 0.0 | 3/3 | 8/6 | 0.006 | 0.011 | 0.002 | 0.00 | 1.8 |
| filleted_bracket | 0.05 | 3/3 | 8/6 | 0.131 | 0.041 | 0.008 | 0.02 | 1.9 |
| filleted_bracket | 0.2 | 0/3 | 10452/0 | — | 0.958 | 0.130 | 0.50 | 13.7 |

"exact" = valid analytic B-rep. A face count equal to the triangle count is the faceted fallback, which is still a valid solid. radius err is the worst |fitted − true| cylinder radius among seeds that came back exact, and — when the part has no fillet radius or no seed was exact. All rows are the worst of 3 seeds from `python benchmarks/run.py`. The 90° countersink (Ø6 hole, Ø12 head, 3 mm deep) at 0.2 mm noise is exact on 1 of 3 seeds. At 0.2 mm the filleted bracket falls back on every seed; the plate is exact on 2 of 3, and one of those solids is missing a corner (7 planes + 3 cylinders). The block stays a valid solid, but one seed keeps only three of the four rounds (Hausdorff 0.85 mm). CI checks the cases that are exact, plus fallbacks that must stay valid solids.

## Limits (v0.7.0)
Planes, cylinders, spheres, cones, and constant-radius fillets between two planes. Cylinders must end on planes perpendicular to their axis (through-holes, bosses, pins). Spheres are exact when each boundary loop is a single plane ∩ sphere circle (domes, spherical bosses, spherical pockets, and zones between two parallel planes). A sphere that runs into an edge or a cylinder falls back to the hybrid solid — vertices are still snapped onto the fitted sphere, but the face stays faceted. Cones are exact when they run between two planes perpendicular to their axis: countersinks, chamfered holes and bosses, drafted/tapered bosses. Cones may also meet a coaxial cylinder (countersunk hole, chamfered pin). A cone that ends in its apex falls back to the hybrid solid. Countersinks and narrow chamfers are exact through about 0.05 mm of noise; a 2 mm chamfer and a wide drafted wall stay exact at 0.2 mm. The countersunk plate's 90° countersink (Ø6 hole, Ø12 head, 3 mm deep) at 0.2 mm noise is exact on one of three seeds. On the other two the conical band and the hole wall meet at a single vertex, so those seeds fall back to a valid faceted solid.

A fillet between two planes is a partial cylinder. When the round runs from one cap to another — rounded vertical corners of a plate or block, a single edge, the inner corner of an L as well as the outer ones — it comes back as one cylinder, tangent to both walls, with the radius in the report's cylinder faces. The radius is the one whose axis is tangent to both walls and which fits the arc; wall triangles that remain in the patch do not pull it off. On the synthetic plate (r=5 mm), block (r=2 mm) and L-bracket (r=3 mm) that rebuild is exact at 0 and 0.05 mm of noise, all three seeds. With no noise the radius error is 0.009 mm on the plate, 0.005 mm on the block and 0.006 mm on the bracket. At 0.05 mm the worst radius error is 0.041 mm on the r=2 corners, 0.100 mm on the r=5 plate and 0.131 mm on the r=3 bracket. At 0.2 mm the bracket falls back on every seed. The plate is exact on 2 of 3 seeds, and one of those solids is missing a corner. The block is a valid solid on every seed, but one seed keeps only three of the four rounds and the Hausdorff rises to 0.85 mm. Whatever does not sew stays a valid faceted solid.

A fillet that runs into a spherical corner (every edge of a box filleted) is not rebuilt: the corner is a sphere octant bounded by three cylinders. A fillet between a plane and a coaxial cylinder — the base of a boss, the rim of a hole — is a torus, and is not rebuilt either. Free-form surfaces are still unsupported. CadQuery and OpenSCAD export write `fillet_r` only for the four vertical corners of a rectangular base (OpenSCAD as `offset(r=fillet_r)` on the base profile). An L-bracket fillet, a spherical corner, and a torus stay on the STEP path.

The parametric scripts are one Z-up box, pads on its top face, and holes and straight cylindrical bosses whose axis is X, Y, or Z. A through hole runs from one face of the base to the opposite face. A blind hole opens on one face; `hole_N_depth` is the distance from that face to the flat end. A side boss sits on one side face and grows straight out. Countersinks, drafts, and chamfers are still vertical only. An axis that is not aligned with X, Y, or Z is refused, and so is a side feature that does not open on a face of the base.

On the 60×40×36 block (Ø8 through hole along X at y=0, z=−6; Ø6 blind hole, 12 deep, opening on +Y at x=−12, z=8; Ø10 boss, 8 tall, on +X at y=10, z=8) both scripts recover those dimensions with no noise and at 0.05 mm noise. With no noise the CadQuery script matches the mesh to 0.006 mm (volume error 0.004%) and the OpenSCAD script to 0.006 mm (volume error 0.003%). At 0.05 mm noise the check against that noisy mesh is 0.139 mm / 0.015% for CadQuery and 0.149 mm / 0.016% for OpenSCAD, under the same limit as every other script, max(4·tol, 0.1) = 0.20 mm. The blind end is not its own plane on this part (the disk stays in the cylinder region), so the depth is the median of the vertices clustered at that end. The exact STEP builder does not sew that cylinder; STEP export of this block falls back to a valid faceted solid, Hausdorff 0.004 mm from the mesh. Anything else outside this vocabulary stays on the STEP path.

## Development
```
pip install -e .[dev]
pytest
python docs/make_demo.py   # regenerate images (needs matplotlib, pillow)
```

License: MIT

<!-- mcp-name: io.github.emirmrkaljevic/mesh2cad -->
