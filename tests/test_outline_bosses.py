"""Straight round bosses on the top face of one Z extrusion."""
import math
import re
import shutil

import pytest

pytest.importorskip("cadquery")
import cadquery as cq  # noqa: E402

from mesh2cad.cadquery_code import check as cq_check  # noqa: E402
from mesh2cad.cadquery_code import run_script, to_cadquery  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.freecad_code import to_freecad  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.openscad_code import check as scad_check  # noqa: E402
from mesh2cad.openscad_code import to_openscad  # noqa: E402
from mesh2cad.parametric import Unsupported  # noqa: E402
from mesh2cad.params import GeometryError, edit_script, list_parameters  # noqa: E402
from parts import l_plate_with_bosses, orthogonal_intersection, to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")

OUTLINE = [(0.0, 0.0), (48.0, 0.0), (48.0, 18.0), (18.0, 18.0), (18.0, 40.0), (0.0, 40.0)]


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _point_block(code, name):
    matched = re.search(rf"{name}\s*=\s*\[", code)
    assert matched, name
    depth = 1
    end = matched.end()
    while end < len(code) and depth:
        depth += (code[end] == "[") - (code[end] == "]")
        end += 1
    nums = [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?", code[matched.end():end - 1])]
    return list(zip(nums[0::2], nums[1::2]))


def _plate_mesh(sigma, seed):
    ref = to_mesh(l_plate_with_bosses())
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=seed) if sigma else ref
    return ref, mesh, tol


def _expect_bosses(code, tol):
    params = _params(code)
    assert params["extrude_len"] == pytest.approx(6.0, abs=tol)
    assert params["boss_1_d"] == pytest.approx(12.0, abs=tol)
    assert params["boss_1_h"] == pytest.approx(8.0, abs=tol)
    assert params["boss_1_x"] == pytest.approx(10.0, abs=tol)
    assert params["boss_1_y"] == pytest.approx(9.0, abs=tol)
    assert params["boss_2_d"] == pytest.approx(10.0, abs=tol)
    assert params["boss_2_h"] == pytest.approx(6.0, abs=tol)
    assert params["boss_2_x"] == pytest.approx(34.0, abs=tol)
    assert params["boss_2_y"] == pytest.approx(9.0, abs=tol)
    assert params["hole_1_d"] == pytest.approx(5.0, abs=tol)
    assert params["hole_1_x"] == pytest.approx(10.0, abs=tol)
    assert params["hole_1_y"] == pytest.approx(9.0, abs=tol)
    assert "boss_1_h" in code and "boss_2_h" in code
    # the through-hole continues out the top of the drilled boss
    assert "boss_1_h" in code.split("hole_1_d", 1)[-1] or "extrude_len + boss_1_h" in code
    return params


def _l():
    return cq.Workplane("XY").polyline(OUTLINE).close().extrude(6)


def _refuse(solid):
    mesh = to_mesh(solid)
    res = convert_mesh(mesh, tol=0.01, with_report=False)
    with pytest.raises(Unsupported) as caught:
        to_cadquery(res.mesh, res.regions, res.tol)
    return str(caught.value)


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_cadquery_bosses_on_an_extruded_outline(sigma, seed):
    ref, mesh, tol = _plate_mesh(sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    limit = max(4 * res.tol, 0.1)
    assert q["max_deviation"] <= limit
    assert q["volume_error_pct"] <= 1.0
    _expect_bosses(code, 0.05 if sigma == 0.0 else 0.15)
    got = _point_block(code, "outline")
    for point in OUTLINE:
        assert any(abs(point[0] - x) <= (0.05 if sigma == 0.0 else 0.25)
                   and abs(point[1] - y) <= (0.05 if sigma == 0.0 else 0.25) for x, y in got)
    assert ".circle(boss_1_d / 2).extrude(boss_1_h)" in code
    assert "extrude_len + boss_1_h" in code
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)
    assert against["volume_error_pct"] < 0.5


@needs_openscad
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_openscad_bosses_on_an_extruded_outline(sigma, seed, tmp_path):
    ref, mesh, tol = _plate_mesh(sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    _expect_bosses(code, 0.05 if sigma == 0.0 else 0.15)
    assert "/* [Bosses] */" in code
    assert "cylinder(h=boss_1_h, d=boss_1_d);" in code
    assert "extrude_len + boss_1_h" in code
    against = scad_check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)
    assert against["volume_error_pct"] < 0.5


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_freecad_macro_fuses_bosses_on_an_extruded_outline(sigma, seed):
    """The macro is written when the same gate passes.

    ``freecadcmd`` is not required: without it the check is CadQuery's, and
    the returned measurement is None.
    """
    _ref, mesh, tol = _plate_mesh(sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q is None or q["max_deviation"] <= max(4 * res.tol, 0.1)
    _expect_bosses(code, 0.05 if sigma == 0.0 else 0.15)
    assert "result.fuse(Part.makeCylinder(boss_1_d / 2, boss_1_h," in code
    assert "extrude_len + boss_1_h" in code
    assert code.index("result.fuse(Part.makeCylinder(boss_1_d / 2, boss_1_h,") < code.index("hole_1_d / 2")


def test_exporters_share_boss_parameters():
    _ref, mesh, tol = _plate_mesh(0.0, 0)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    scad, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    macro, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert _params(scad) == _params(cq) == _params(macro)
    slider = re.search(r"^boss_1_d\s*=\s*[\d.]+;\s*//\s*\[([^\]]+)\]", scad, re.M)
    assert slider
    assert float(slider.group(1).split(":")[-1]) >= _params(scad)["boss_1_d"]


def test_edit_lists_and_changes_an_outline_boss(tmp_path):
    _ref, mesh, tol = _plate_mesh(0.0, 0)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    script = tmp_path / "plate.py"
    script.write_text(code, encoding="utf-8")
    listed = {item["name"]: item for item in list_parameters(str(script))["parameters"]}
    assert listed["boss_2_d"]["kind"] == "diameter"
    assert listed["boss_2_h"]["kind"] == "length"
    assert listed["boss_2_x"]["kind"] == "position"
    assert listed["boss_1_h"]["value"] == pytest.approx(8.0)
    before = abs(run_script(code).volume)
    edit_script(str(script), {"boss_2_h": 9}, output_path=str(tmp_path / "out.step"), in_place=True)
    text = script.read_text(encoding="utf-8")
    assert re.search(r"^boss_2_h\s*=\s*9\.0", text, re.M)
    # the solid boss grows by a cylinder; the drilled boss is unchanged
    assert abs(run_script(text).volume) - before == pytest.approx(math.pi * 5.0 ** 2 * 3.0, rel=0.01)
    with pytest.raises(GeometryError, match="boss_2_d"):
        edit_script(str(script), {"boss_2_x": 47}, output_path=str(tmp_path / "no.step"))


def test_boss_partly_outside_the_outline_is_refused():
    # centered on the inner corner, so the circle leaves the L
    solid = _l().union(cq.Workplane().add(
        cq.Solid.makeCylinder(6.0, 8.0, cq.Vector(18.0, 18.0, 6.0), cq.Vector(0, 0, 1))))
    assert "a boss is partly outside the extruded outline" in _refuse(solid)


def test_boss_on_the_bottom_face_is_refused():
    solid = _l().union(cq.Workplane().add(
        cq.Solid.makeCylinder(5.0, 8.0, cq.Vector(10.0, 9.0, -8.0), cq.Vector(0, 0, 1))))
    assert "a boss on the bottom face of an extruded outline is not supported yet" in _refuse(solid)


def test_boss_on_a_side_face_is_refused():
    solid = _l().union(cq.Workplane().add(
        cq.Solid.makeCylinder(4.0, 8.0, cq.Vector(48.0, 9.0, 3.0), cq.Vector(1, 0, 0))))
    assert "a boss on a side face of an extruded outline is not supported yet" in _refuse(solid)


def test_tapered_boss_on_an_outline_is_refused():
    solid = _l().faces(">Z").workplane().center(10, 9).circle(6).extrude(8, taper=8)
    message = _refuse(solid)
    assert "a tapered boss on an extruded outline is not supported yet" in message


def test_hole_that_meets_a_boss_and_stops_is_refused():
    """A hole through the plate only, covered by a solid boss, is not rewritten."""
    plate = _l().cut(cq.Workplane().add(
        cq.Solid.makeCylinder(2.5, 20.0, cq.Vector(10.0, 9.0, -1.0), cq.Vector(0, 0, 1))))
    solid = plate.union(cq.Workplane().add(
        cq.Solid.makeCylinder(6.0, 8.0, cq.Vector(10.0, 9.0, 6.0), cq.Vector(0, 0, 1))))
    assert "a hole meets a boss but does not pass through it" in _refuse(solid)


def test_boss_on_a_combined_extrusion_is_refused():
    solid = orthogonal_intersection().union(cq.Workplane().add(
        cq.Solid.makeCylinder(4.0, 6.0, cq.Vector(10.0, 10.0, 12.0), cq.Vector(0, 0, 1))))
    assert "a boss on a combined extrusion is not supported yet" in _refuse(solid)


def test_boss_on_two_extrusions_is_refused():
    """A union of two Z outlines is not one extrusion, so the boss stays refused."""
    base = _l()
    step = (cq.Workplane("XY").polyline([(20, 22), (40, 22), (40, 36), (28, 36)])
            .close().extrude(4).translate((0, 0, 6)))
    boss = cq.Solid.makeCylinder(4.0, 5.0, cq.Vector(8.0, 8.0, 6.0), cq.Vector(0, 0, 1))
    message = _refuse(base.union(step).union(cq.Workplane().add(boss)))
    assert "a boss on a combined extrusion is not supported yet" in message
