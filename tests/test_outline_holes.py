"""Countersunk and counterbored through-holes on one extruded outline."""
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
from parts import countersunk_plate, hex_flange, to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")

# Flat-to-flat along X, same loop hex_flange extrudes.
_R = 28.0
HEX = [
    (_R * math.cos(math.radians(60 * i - 30)), _R * math.sin(math.radians(60 * i - 30)))
    for i in range(6)
]


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


def _mesh(face, sigma, seed):
    ref = to_mesh(hex_flange(face))
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=seed) if sigma else ref
    return ref, mesh, tol


def _expect(code, face, tol):
    params = _params(code)
    end = "+Z" if face == "top" else "-Z"
    assert "base_x" not in code
    assert "extrude_len" in code
    assert "hole_1_chamfer" not in code
    assert params["extrude_len"] == pytest.approx(6.0, abs=tol)
    assert params["hole_1_d"] == pytest.approx(6.0, abs=tol)
    assert params["hole_1_x"] == pytest.approx(-12.0, abs=tol)
    assert params["hole_1_csk_d"] == pytest.approx(12.0, abs=tol)
    assert params["hole_1_csk_angle"] == pytest.approx(90.0, abs=0.5)
    assert f"countersink diameter at the {end} end" in code
    assert params["hole_2_d"] == pytest.approx(5.5, abs=tol)
    assert params["hole_2_x"] == pytest.approx(12.0, abs=tol)
    assert params["hole_2_cbore_d"] == pytest.approx(10.0, abs=tol)
    assert params["hole_2_cbore_depth"] == pytest.approx(5.0, abs=tol)
    assert f"counterbore depth from the {end} end" in code
    return params


def _outline_ok(code, tol):
    got = _point_block(code, "outline")
    assert len(got) == len(HEX)
    for point in HEX:
        assert any(abs(point[0] - x) <= tol and abs(point[1] - y) <= tol for x, y in got)


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("face", ["top", "bottom"])
def test_cadquery_outline_countersink_and_counterbore(face, sigma, seed):
    ref, mesh, tol = _mesh(face, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    limit = max(4 * res.tol, 0.1)
    assert q["max_deviation"] <= limit
    assert q["volume_error_pct"] <= 1.0
    near = 0.05 if sigma == 0.0 else 0.2
    _expect(code, face, near)
    _outline_ok(code, 0.05 if sigma == 0.0 else 0.25)
    assert "Solid.makeCone(" in code
    assert "hole_1_csk_d / 2" in code
    assert "hole_2_cbore_d / 2" in code
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)
    assert against["volume_error_pct"] < 0.5


@needs_openscad
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("face", ["top", "bottom"])
def test_openscad_outline_countersink_and_counterbore(face, sigma, seed, tmp_path):
    ref, mesh, tol = _mesh(face, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    _expect(code, face, 0.05 if sigma == 0.0 else 0.2)
    assert "hole_1_csk_h + 1" in code
    assert "r1=" in code and "r2=" in code
    assert "d=hole_2_cbore_d" in code
    assert re.search(r"^hole_1_csk_angle\s*=\s*[\d.]+;\s*//\s*\[60:120\]", code, re.M)
    against = scad_check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)
    assert against["volume_error_pct"] < 0.5


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("face", ["top", "bottom"])
def test_freecad_macro_cuts_outline_countersink_and_counterbore(face, sigma, seed):
    """The macro is written when the same gate passes.

    ``freecadcmd`` is not required. Without it the check is CadQuery's.
    """
    _ref, mesh, tol = _mesh(face, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q is None or q["max_deviation"] <= max(4 * res.tol, 0.1)
    _expect(code, face, 0.05 if sigma == 0.0 else 0.2)
    assert "Part.makeCone(" in code
    assert "hole_1_csk_angle" in code
    assert "Part.makeCylinder(hole_2_cbore_d / 2" in code


def test_exporters_share_outline_hole_parameters():
    _ref, mesh, tol = _mesh("top", 0.0, 0)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    scad, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    macro, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert _params(scad) == _params(cq) == _params(macro)
    assert _point_block(scad, "outline") == _point_block(cq, "outline")


def test_rectangular_countersink_stays_a_box():
    """A rectangular plate keeps base_x. The outline vocabulary is not used."""
    ref = to_mesh(countersunk_plate())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    params = _params(code)
    assert params["base_x"] == pytest.approx(40.0, abs=0.05)
    assert params["hole_1_csk_d"] == pytest.approx(12.0, abs=0.05)
    assert params["hole_1_csk_angle"] == pytest.approx(90.0, abs=0.5)
    assert "extrude_len" not in code
    assert "outline" not in code


def _frustum(r1, r2, height):
    return math.pi * height * (r1 * r1 + r1 * r2 + r2 * r2) / 3.0


def test_edit_outline_countersink_and_counterbore(tmp_path):
    _ref, mesh, tol = _mesh("top", 0.0, 0)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    script = tmp_path / "flange.py"
    script.write_text(code, encoding="utf-8")
    listed = {item["name"]: item for item in list_parameters(str(script))["parameters"]}
    assert listed["hole_1_csk_d"]["kind"] == "diameter"
    assert listed["hole_1_csk_angle"]["kind"] == "angle"
    assert listed["hole_2_cbore_depth"]["kind"] == "depth"
    v0 = abs(run_script(code).volume)
    edit_script(str(script), {"hole_1_csk_d": 14.0}, in_place=True)
    # 90 degree countersink, Ø12 head (3 mm deep) to Ø14 (4 mm deep).
    old = _frustum(3.0, 6.0, 3.0) - math.pi * 9.0 * 3.0
    new = _frustum(3.0, 7.0, 4.0) - math.pi * 9.0 * 4.0
    assert v0 - abs(run_script(script.read_text(encoding="utf-8")).volume) == pytest.approx(new - old, rel=0.02)
    script.write_text(code, encoding="utf-8")
    edit_script(str(script), {"hole_2_cbore_d": 12.0}, in_place=True)
    # The counterbore is 5 mm deep. Only the annulus changes.
    expect = math.pi * (6.0 ** 2 - 5.0 ** 2) * 5.0
    assert v0 - abs(run_script(script.read_text(encoding="utf-8")).volume) == pytest.approx(expect, rel=0.02)
    script.write_text(code, encoding="utf-8")
    with pytest.raises(GeometryError, match="countersink"):
        edit_script(str(script), {"hole_1_d": 12.0})
    with pytest.raises(GeometryError, match="does not fit"):
        edit_script(str(script), {"hole_1_csk_d": 50.0})


def test_mouth_chamfer_is_a_90_degree_countersink():
    """A 45 degree mouth on a Z hole is the 90 degree countersink."""
    outline = HEX
    solid = cq.Workplane("XY").polyline(outline).close().extrude(6)
    hole = cq.Solid.makeCylinder(3.0, 10.0, cq.Vector(0.0, 0.0, -2.0), cq.Vector(0, 0, 1))
    # Axial length 1.5 mm, so the mouth diameter is 6 + 3 = 9.
    cone = cq.Solid.makeCone(3.0, 4.5, 1.5, cq.Vector(0.0, 0.0, 4.5), cq.Vector(0, 0, 1))
    part = solid.cut(cq.Workplane().add(hole)).cut(cq.Workplane().add(cone))
    ref = to_mesh(part)
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    params = _params(code)
    assert params["hole_1_d"] == pytest.approx(6.0, abs=0.05)
    assert params["hole_1_csk_d"] == pytest.approx(9.0, abs=0.05)
    assert params["hole_1_csk_angle"] == pytest.approx(90.0, abs=0.5)
    assert "hole_1_chamfer" not in code
    assert "extrude_len" in code


def test_non_45_countersink_on_an_outline():
    """82 degrees is not a chamfer. It stays hole_1_csk_angle."""
    part = cq.Workplane("XY").polyline(HEX).close().extrude(6)
    part = part.faces(">Z").workplane().pushPoints([(0.0, 0.0)]).cskHole(6, 12, 82)
    ref = to_mesh(part)
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    params = _params(code)
    assert params["hole_1_csk_d"] == pytest.approx(12.0, abs=0.05)
    assert params["hole_1_csk_angle"] == pytest.approx(82.0, abs=0.5)
    assert "hole_1_chamfer" not in code


def test_draft_on_an_outline_is_refused():
    """A cone that replaces the hole, rather than opening one mouth, is refused."""
    part = cq.Workplane("XY").polyline(HEX).close().extrude(6)
    cone = cq.Solid.makeCone(2.0, 5.0, 8.0, cq.Vector(0.0, 0.0, -1.0), cq.Vector(0, 0, 1))
    part = part.cut(cq.Workplane().add(cone))
    mesh = to_mesh(part)
    res = convert_mesh(mesh, tol=0.01, with_report=False)
    with pytest.raises(Unsupported, match="not supported yet"):
        to_cadquery(res.mesh, res.regions, res.tol)
