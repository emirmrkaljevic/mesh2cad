"""Sharp edges and short chamfers on a coarse CAD tessellation.

A ring of flat faces whose corners lie on a circle used to come back as one
cylinder or one circumscribed sphere. A 90 degree countersink only a few
tenths of a millimetre tall used to be dropped, so the hole was a cylinder
with no cone. Both are analytic surfaces on a clean mesh.
"""
import math

import cadquery as cq
import numpy as np
import trimesh

from mesh2cad.convert import convert_mesh


def _coarse(wp, tol, ang):
    v, f = wp.val().tessellate(tol, ang)
    m = trimesh.Trimesh(np.array([p.toTuple() for p in v]), np.array(f), process=True)
    if len(m.faces) and (m.area_faces < 1e-12).any():
        m.update_faces(m.area_faces >= 1e-12)
        m.remove_unreferenced_vertices()
    return m


def _octagon_boss():
    radius = 6.0
    pts = [(radius * math.cos(i * math.pi / 4), radius * math.sin(i * math.pi / 4))
           for i in range(8)]
    return (cq.Workplane().box(120, 120, 6).faces(">Z").workplane()
            .polyline(pts).close().extrude(8))


def _short_countersink():
    """80x80x10 plate, Ø3 hole, Ø4.2 head, 90 degree countersink.

    The conical band is about 0.6 mm tall. A coarse tessellation used to
    swallow it.
    """
    return cq.Workplane().box(80, 80, 10).faces(">Z").workplane().cskHole(3, 4.2, 90)


def test_octagon_boss_is_planes_not_a_round():
    mesh = _coarse(_octagon_boss(), tol=0.08, ang=0.35)
    res = convert_mesh(mesh, ref=mesh)
    report = res.report
    assert report["method"] == "exact"
    assert report["output"]["valid_solid"]
    assert report["output"]["faces"] == {"plane": 15}
    assert report["regions"]["by_type"] == {"plane": 15}
    assert report["deviation"]["hausdorff"] < 1e-6
    assert report["volume_error_pct"] < 0.01


def test_short_countersink_keeps_the_cone():
    mesh = _coarse(_short_countersink(), tol=0.05, ang=0.3)
    res = convert_mesh(mesh, ref=mesh)
    report = res.report
    assert report["method"] == "exact"
    assert report["output"]["valid_solid"]
    assert report["output"]["faces"] == {"plane": 6, "cylinder": 1, "cone": 1}
    cones = [r.surface for r in res.regions if r.surface is not None and r.surface.kind == "cone"]
    cylinders = [r.surface for r in res.regions if r.surface is not None and r.surface.kind == "cylinder"]
    assert len(cones) == 1 and len(cylinders) == 1
    assert abs(cones[0].half_angle - math.pi / 4) < 0.05
    assert abs(cylinders[0].radius - 1.5) < 0.05
    assert report["deviation"]["hausdorff"] < 0.05
    assert report["volume_error_pct"] < 0.05
