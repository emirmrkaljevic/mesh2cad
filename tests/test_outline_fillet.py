"""Vertical rounds on one Z extrusion of an outline, as outline_fillet_r."""
import math
import re
import shutil

import pytest

pytest.importorskip("cadquery")
import cadquery as cq  # noqa: E402

from mesh2cad.cadquery_code import check as cq_check  # noqa: E402
from mesh2cad.cadquery_code import run_script, to_cadquery  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.freecad_code import freecad_bin, to_freecad  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.openscad_code import check as scad_check  # noqa: E402
from mesh2cad.openscad_code import to_openscad  # noqa: E402
from mesh2cad.parametric import Unsupported  # noqa: E402
from mesh2cad.params import GeometryError, edit_script, list_parameters  # noqa: E402
from parts import filleted_t_plate, filleted_t_with_boss_and_countersink, to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")
needs_freecad = pytest.mark.skipif(freecad_bin() is None, reason="freecadcmd is not installed")

# Sharp corners of the T, before the r=3 rounds.
OUTLINE = [(16, 0), (34, 0), (34, 28), (50, 28), (50, 42), (0, 42), (0, 28), (16, 28)]
# (diameter, x, y). The bar hole is the wider one.
HOLES = [(4.0, 25.0, 14.0), (5.0, 25.0, 35.0)]



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


def _mesh(builder, sigma, seed):
    ref = to_mesh(builder())
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=seed) if sigma else ref
    return ref, mesh, tol


def _export(fn, res):
    """The script when the gate passes, or None when the fit is refused.

    At 0.05 mm one seed sits on the gate. Which tessellation crosses it is
    not stable, so a refusal is kept and a written script is checked.
    """
    try:
        return fn(res.mesh, res.regions, res.tol)
    except Unsupported:
        return None


def _expect_plate(code, near):
    params = _params(code)
    assert "base_x" not in params
    assert "fillet_r" not in params
    assert "edge_fillet_r" not in params
    assert "leg_x" not in params
    assert params["extrude_len"] == pytest.approx(6.0, abs=near)
    assert params["outline_fillet_r"] == pytest.approx(3.0, abs=near)
    got = _point_block(code, "outline")
    assert len(got) == len(OUTLINE)
    for point in OUTLINE:
        assert any(abs(point[0] - x) <= near and abs(point[1] - y) <= near for x, y in got)
    holes = [(params[f"hole_{i}_d"], params[f"hole_{i}_x"], params[f"hole_{i}_y"]) for i in (1, 2)]
    for expect in HOLES:
        assert any(abs(expect[0] - d) <= near and abs(expect[1] - x) <= near and abs(expect[2] - y) <= near
                   for d, x, y in holes)
    return params


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_cadquery_outline_fillet(sigma, seed):
    ref, mesh, tol = _mesh(filleted_t_plate, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    got = _export(to_cadquery, res)
    if sigma == 0.0:
        assert got is not None
    if got is None:
        return
    code, q = got
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    _expect_plate(code, 0.05 if sigma == 0.0 else 0.4)
    assert '.edges("|Z").fillet(outline_fillet_r)' in code
    assert code.index('.edges("|Z").fillet(outline_fillet_r)') < code.index("makeCylinder")
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.4)
    assert against["volume_error_pct"] < 0.5


@needs_openscad
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_openscad_outline_fillet(sigma, seed, tmp_path):
    ref, mesh, tol = _mesh(filleted_t_plate, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    got = _export(to_openscad, res)
    if sigma == 0.0:
        assert got is not None
    if got is None:
        return
    code, q = got
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    _expect_plate(code, 0.05 if sigma == 0.0 else 0.4)
    assert "offset(r=outline_fillet_r)" in code
    assert "offset(r=-outline_fillet_r)" in code
    assert code.index("offset(r=-outline_fillet_r)") < code.index("cylinder(")
    against = scad_check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.4)
    assert against["volume_error_pct"] < 0.5


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_freecad_macro_fillets_outline_corners(sigma, seed):
    """The macro is written when the same gate passes.

    ``freecadcmd`` is not required. Without it the check is CadQuery's.
    """
    _ref, mesh, tol = _mesh(filleted_t_plate, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    got = _export(to_freecad, res)
    if sigma == 0.0:
        assert got is not None
    if got is None:
        return
    code, q = got
    assert q is None or q["max_deviation"] <= max(4 * res.tol, 0.1)
    _expect_plate(code, 0.05 if sigma == 0.0 else 0.4)
    assert "result.makeFillet(outline_fillet_r, _fillet_edges)" in code
    assert code.index("result.makeFillet(outline_fillet_r, _fillet_edges)") < code.index("Part.makeCylinder")


@needs_freecad
@pytest.mark.parametrize("seed", [0])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_freecad_outline_fillet_deviation(sigma, seed, tmp_path):
    from mesh2cad.freecad_code import check as fc_check

    ref, mesh, tol = _mesh(filleted_t_plate, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    against = fc_check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)


def test_exporters_share_outline_fillet_parameters():
    _ref, mesh, tol = _mesh(filleted_t_plate, 0.0, 0)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    scad, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    macro, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert _params(scad) == _params(cq) == _params(macro)
    assert _point_block(scad, "outline") == _point_block(cq, "outline")


def test_edit_outline_fillet_radius(tmp_path):
    _ref, mesh, tol = _mesh(filleted_t_plate, 0.0, 0)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    script = tmp_path / "plate.py"
    script.write_text(code, encoding="utf-8")
    listed = {item["name"]: item for item in list_parameters(str(script))["parameters"]}
    assert listed["outline_fillet_r"]["kind"] == "radius"
    assert listed["outline_fillet_r"]["unit"] == "mm"
    assert listed["outline_fillet_r"]["value"] == pytest.approx(3.0)
    before = abs(run_script(code).volume)
    result = edit_script(str(script), {"outline_fillet_r": 4}, output_path=str(tmp_path / "out.step"),
                         in_place=True)
    text = script.read_text(encoding="utf-8")
    assert result["valid_solid"] is True
    assert re.search(r"^outline_fillet_r\s*=\s*4\.0", text, re.M)
    # Six outside corners and two inside corners. The net of a 90 degree round
    # is four pieces of r^2*(1-pi/4), times the 6 mm thickness.
    piece = (1.0 - math.pi / 4.0) * (16.0 - 9.0) * 6.0
    assert before - abs(run_script(text).volume) == pytest.approx(4.0 * piece, rel=0.02)
    with pytest.raises(GeometryError, match="outline_fillet_r"):
        edit_script(str(script), {"outline_fillet_r": 8}, output_path=str(tmp_path / "no.step"))
    assert re.search(r"^outline_fillet_r\s*=\s*4\.0", script.read_text(encoding="utf-8"), re.M)
    assert not (tmp_path / "no.step").exists()


def test_sharp_t_is_not_given_a_fillet():
    outline = list(OUTLINE)
    plate = cq.Workplane("XY").polyline(outline).close().extrude(6)
    bar = cq.Solid.makeCylinder(2.5, 20.0, cq.Vector(25.0, 35.0, -2.0), cq.Vector(0, 0, 1))
    stem = cq.Solid.makeCylinder(2.0, 20.0, cq.Vector(25.0, 14.0, -2.0), cq.Vector(0, 0, 1))
    solid = plate.cut(cq.Workplane().add(bar)).cut(cq.Workplane().add(stem))
    ref = to_mesh(solid)
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert "outline_fillet_r" not in code
    assert "fillet_r" not in code
    assert "edges(\"|Z\").fillet" not in code


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_cadquery_fillet_with_boss_and_countersink(sigma, seed):
    ref, mesh, tol = _mesh(filleted_t_with_boss_and_countersink, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    got = _export(to_cadquery, res)
    if sigma == 0.0:
        assert got is not None
    if got is None:
        return
    code, q = got
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    near = 0.05 if sigma == 0.0 else 0.4
    assert params["outline_fillet_r"] == pytest.approx(3.0, abs=near)
    assert params["hole_1_d"] == pytest.approx(5.0, abs=near)
    assert params["hole_1_csk_d"] == pytest.approx(9.0, abs=near)
    assert params["hole_1_csk_angle"] == pytest.approx(90.0, abs=0.5 if sigma == 0.0 else 2.0)
    assert params["boss_1_d"] == pytest.approx(8.0, abs=near)
    assert params["boss_1_h"] == pytest.approx(5.0, abs=near)
    assert params["boss_1_x"] == pytest.approx(25.0, abs=near)
    assert params["boss_1_y"] == pytest.approx(12.0, abs=near)
    fillet_at = code.index('.edges("|Z").fillet(outline_fillet_r)')
    assert fillet_at < code.index("boss_1_d / 2")
    assert fillet_at < code.index("makeCone")
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.4)


@needs_openscad
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_openscad_fillet_with_boss_and_countersink(sigma, seed, tmp_path):
    ref, mesh, tol = _mesh(filleted_t_with_boss_and_countersink, sigma, seed)
    res = convert_mesh(mesh, tol=tol, with_report=False)
    got = _export(to_openscad, res)
    if sigma == 0.0:
        assert got is not None
    if got is None:
        return
    code, q = got
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    near = 0.05 if sigma == 0.0 else 0.4
    assert params["outline_fillet_r"] == pytest.approx(3.0, abs=near)
    assert params["hole_1_csk_d"] == pytest.approx(9.0, abs=near)
    assert params["boss_1_d"] == pytest.approx(8.0, abs=near)
    assert "offset(r=outline_fillet_r)" in code
    against = scad_check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.4)
