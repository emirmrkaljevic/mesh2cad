"""CAD-style tessellation: coarse triangles, no subdivision.

`parts.to_mesh` splits every edge to <= 2 mm, which hides a failure mode of
real STL export. A planar face there is often one quad. These meshes are the
cadquery solid tessellated once and left as-is. On the old segmenter the
stepped block came back faceted (its walls sit on a circumsphere, and a
two-triangle face was discarded).
"""
import numpy as np
import pytest
import trimesh

from mesh2cad.convert import convert_mesh
from parts import boss, cube, l_bracket, plate_with_hole, stepped_block, to_mesh

# (name, builder, output faces, hausdorff limit in mm)
CASES = [
    ("cube", cube, {"plane": 6}, 1e-6),
    ("l_bracket", l_bracket, {"plane": 8}, 1e-6),
    ("stepped_block", stepped_block, {"plane": 11}, 1e-6),
    ("plate_with_hole", plate_with_hole, {"plane": 6, "cylinder": 1}, 0.05),
    ("boss", boss, {"plane": 7, "cylinder": 1}, 0.1),
]


def coarse_mesh(wp, tol=0.1, ang=0.4):
    v, f = wp.val().tessellate(tol, ang)
    m = trimesh.Trimesh(np.array([p.toTuple() for p in v]), np.array(f), process=True)
    if len(m.faces) and (m.area_faces < 1e-12).any():
        m.update_faces(m.area_faces >= 1e-12)
        m.remove_unreferenced_vertices()
    return m


@pytest.mark.parametrize("name,builder,faces,haus", CASES)
def test_coarse_tessellation_is_exact(name, builder, faces, haus):
    wp = builder()
    m = coarse_mesh(wp)
    dense = to_mesh(wp)
    assert len(m.faces) < len(dense.faces) / 10
    r = convert_mesh(m, ref=m).report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == faces
    assert r["regions"]["by_type"].get("faceted", 0) == 0
    assert r["deviation"]["hausdorff"] < haus
    assert r["volume_error_pct"] < 0.5
