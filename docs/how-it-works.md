# How it works

The converter loads the mesh, groups triangles into primitive surfaces, builds an exact solid where those surfaces sew, and otherwise falls back to a faceted solid. A fidelity report is always computed against the input mesh.

1. Load + clean mesh (trimesh). Noisy meshes get a Taubin-smoothed copy for normals. A clean CAD tessellation (dihedral angles piled onto a few exact values, including an export whose faces are single quads) is segmented on the raw normals, and a face of one or two triangles is kept.
2. Region growing into planes; strip chains and leftovers are fitted as cylinders, spheres, or cones (normal-ray / algebraic / constant-normal-angle fits, robust nonlinear refinement). A narrow round between two planes is gathered back into one partial cylinder and snapped so its generators lie in those planes. A coaxial cylinder + cone component (countersink, chamfered pin) is split by normal-to-axis angle first. Small noisy patches that tile a sphere are merged into one spherical region. Split/redundant regions are merged and absorbed.
3. Exact B-rep (OCP / OpenCascade): plane polygons from triple-plane intersections, circles from plane ∩ cylinder and plane ∩ sphere, plane ∩ cone, full cylindrical walls, partial fillet cylinders, quarter tori tangent to a plane and a coaxial cylinder, spherical caps and zones, sphere octants where three edge fillets meet, conical frustums, sewn into a solid and validated.
4. Fallbacks: hybrid faceted solid (vertices snapped to fitted surfaces, coplanar triangles merged), then raw faceted solid.
5. Fidelity report against the input mesh.

What is refused, and the noise levels where a fit stops being exact, are in [Limits](limits.md). Measured numbers are in [Benchmarks](benchmarks.md) and [Real parts](real-parts.md).

Back to the [README](../README.md).
