"""A cross-section that changes along one axis, cut back in segments.

The bar is full height for only part of its length. A hexagonal tunnel and a
pentagonal slit occupy part of that length; neither is a cylinder. The sweep
keeps the section on each span and cuts those openings as polygons. The gate
is the same one every other script uses.
"""
import shutil

import pytest
import trimesh

pytest.importorskip("cadquery")
import cadquery as cq  # noqa: E402

from mesh2cad.cadquery_code import check as cq_check  # noqa: E402
from mesh2cad.cadquery_code import to_cadquery  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.openscad_code import check as scad_check  # noqa: E402
from mesh2cad.openscad_code import to_openscad  # noqa: E402
from mesh2cad.sweep import try_sweep  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")


def _ramped_bar():
    """Full 16 x 12 section for x=0..16, then 16 x 5. A hex and a pentagon slit.

    The hex runs from x=4 to x=12. The slit is 0.3 mm along X at x=20, in the
    short section, so it is not a station of the outer profile.
    """
    tall = cq.Workplane("YZ").rect(16, 12).extrude(16)
    short = cq.Workplane("YZ").workplane(offset=16).rect(16, 5).extrude(12)
    bar = tall.union(short)
    pocket = cq.Workplane("YZ").workplane(offset=4).polygon(6, 5).extrude(8)
    slit = (
        cq.Workplane("YZ")
        .workplane(offset=20)
        .polyline([(2, -1.5), (4, -1.5), (5, -0.4), (3.2, 1.2), (2, 0.4)])
        .close()
        .extrude(0.3)
    )
    return bar.cut(pocket).cut(slit)


def _mesh(solid):
    import numpy as np

    verts, faces = solid.val().tessellate(0.05, 0.2)
    return trimesh.Trimesh(np.array([p.toTuple() for p in verts]), np.array(faces), process=True)


def test_sweep_follows_a_changing_section():
    mesh = _mesh(_ramped_bar())
    res = convert_mesh(mesh, tol=0.05, with_report=False)
    limit = max(4 * res.tol, 0.1)
    swept = try_sweep(res.mesh, res.regions, res.tol)
    assert swept is not None
    score, result = swept
    assert score["max_deviation"] <= limit
    assert score["volume_error_pct"] <= 1.0
    # the hex or the pentagon is a cut profile, not a circle and not a rectangle
    assert any(len(c["outline"]) >= 5 for c in result["csg"]["cuts"])
    assert any(c.get("cutouts") for c in result["csg"]["cuts"])
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= limit
    assert q["volume_error_pct"] <= 1.0
    against = cq_check(code, mesh)
    assert against["max_deviation"] <= limit
    assert against["volume_error_pct"] <= 1.0


@needs_openscad
def test_openscad_sweep_follows_a_changing_section(tmp_path):
    mesh = _mesh(_ramped_bar())
    res = convert_mesh(mesh, tol=0.05, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    limit = max(4 * res.tol, 0.1)
    assert q["max_deviation"] <= limit
    assert q["volume_error_pct"] <= 1.0
    against = scad_check(code, mesh, tmp_path)
    assert against["max_deviation"] <= limit
    assert against["volume_error_pct"] <= 1.0
