# mesh2cad
Convert 3D meshes (STL/OBJ/PLY) into real CAD (STEP) with fitted planes and cylinders — CLI + MCP server for AI agents.

> Status: early development.

```
pip install mesh2cad
mesh2cad convert part.stl -o part.step --tol 0.1
```
Outputs a STEP file plus a fidelity report (Hausdorff distance, RMS deviation, % surface fitted).

License: MIT
