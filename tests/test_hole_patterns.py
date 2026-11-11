"""Rows and grids of plain Z through-holes as pattern parameters.

Fixtures are built in this process. Nothing is loaded from disk. Two holes
stay hole_N_*: one gap is not a pattern. A line whose gaps are not equal
stays one parameter set per hole. A pattern the deviation gate rejects is
written the same way.
"""
import math
import re
import shutil

import pytest

pytest.importorskip("cadquery")
import cadquery as cq  # noqa: E402

from mesh2cad.cadquery_code import check, run_script, to_cadquery  # noqa: E402
from mesh2cad.cli import main  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.freecad_code import to_freecad  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.openscad_code import to_openscad  # noqa: E402
from mesh2cad.patterns import partition_holes  # noqa: E402
from mesh2cad.params import GeometryError, edit_script, list_parameters  # noqa: E402
from parts import _cut_z_hole, _filleted_bracket_solid, to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")


def _hole(x, y=0.0, d=6.0, **extra):
    hole = {"x": float(x), "y": float(y), "d": float(d), "axis": "Z"}
    hole.update(extra)
    return hole


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _plate(length, width, thick, points, diameter):
    part = cq.Workplane().box(length, width, thick)
    return part.faces(">Z").workplane().pushPoints(points).hole(diameter)


def plate_row():
    """100 x 40 x 6 plate, four Ø6 holes on y=0, pitch 12, from x=-18."""
    return _plate(100, 40, 6, [(-18, 0), (-6, 0), (6, 0), (18, 0)], 6)


def plate_grid():
    """100 x 70 x 6 plate, Ø5 holes on a 3 by 2 grid, pitch 20 and 24."""
    points = [(x, y) for x in (-20, 0, 20) for y in (-12, 12)]
    return _plate(100, 70, 6, points, 5)


def plate_near_miss():
    """Four collinear holes whose gaps are 12, 14 and 10 mm. Not one pitch."""
    return _plate(100, 40, 6, [(-18, 0), (-6, 0), (8, 0), (18, 0)], 6)


def plate_mild_miss():
    """Same row, last center 0.6 mm long. The detector still calls it a row.

    The snapped solid is farther from the mesh than the deviation gate, so
    the script that is written keeps one parameter set per hole.
    """
    return _plate(100, 40, 6, [(-18, 0), (-6, 0), (6, 0), (18.6, 0)], 6)


def plate_mixed():
    """The row of four, plus one Ø6 hole off the line at (0, 14)."""
    points = [(-18, 0), (-6, 0), (6, 0), (18, 0), (0, 14)]
    return _plate(100, 50, 6, points, 6)


def plate_vertical():
    """Three Ø6 holes along Y, pitch 16."""
    return _plate(40, 80, 6, [(0, -16), (0, 0), (0, 16)], 6)


def plate_filleted_row():
    """90 x 40 x 16 box, every edge r=2, then four Ø6 holes along X."""
    part = cq.Workplane().box(90, 40, 16).edges().fillet(2)
    return part.faces(">Z").workplane().pushPoints([(-24, 0), (-8, 0), (8, 0), (24, 0)]).hole(6)


def bracket_foot_row():
    """Filleted L with three Ø4 holes in the foot, pitch 8, at y=4."""
    part = _filleted_bracket_solid()
    for x in (16.0, 24.0, 32.0):
        part = _cut_z_hole(part, x, 4.0)
    return part


def _export(wp, sigma=0.0, seed=0):
    ref = to_mesh(wp)
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=seed) if sigma else ref
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    return ref, res, code, q


def test_detector_groups_rows_and_grids_and_leaves_the_rest():
    kept, rows, grids = partition_holes(
        [_hole(-18), _hole(-6), _hole(6), _hole(18)], 0.01)
    assert kept == []
    assert len(rows) == 1 and rows[0]["n"] == 4
    assert rows[0]["pitch"] == pytest.approx(12)
    assert rows[0]["angle"] == pytest.approx(0)
    assert rows[0]["x"] == pytest.approx(-18)

    kept, rows, grids = partition_holes([_hole(-18), _hole(18)], 0.01)
    assert rows == [] and grids == [] and len(kept) == 2

    kept, rows, grids = partition_holes(
        [_hole(-18), _hole(-6), _hole(8), _hole(18)], 0.01)
    assert rows == [] and grids == [] and len(kept) == 4

    # The detector accepts this row. The export gate is what refuses it.
    kept, rows, grids = partition_holes(
        [_hole(-18), _hole(-6), _hole(6), _hole(18.6)], 0.01)
    assert len(rows) == 1 and rows[0]["n"] == 4
    assert rows[0]["pitch"] == pytest.approx(12.18, abs=0.02)

    points = [(x, y) for x in (-20, 0, 20) for y in (-12, 12)]
    kept, rows, grids = partition_holes([_hole(x, y, d=5) for x, y in points], 0.01)
    assert rows == [] and kept == []
    assert grids[0]["nx"] == 3 and grids[0]["ny"] == 2
    assert grids[0]["pitch_x"] == pytest.approx(20)
    assert grids[0]["pitch_y"] == pytest.approx(24)
    assert grids[0]["x"] == pytest.approx(-20) and grids[0]["y"] == pytest.approx(-12)

    kept, rows, _grids = partition_holes(
        [_hole(0, 0), _hole(10, 10), _hole(20, 20), _hole(30, 30)], 0.01)
    assert kept == [] and rows[0]["n"] == 4
    assert rows[0]["angle"] == pytest.approx(45)
    assert rows[0]["pitch"] == pytest.approx(10 * math.sqrt(2))

    kept, rows, _grids = partition_holes(
        [_hole(0, -18), _hole(0, -6), _hole(0, 6), _hole(0, 18)], 0.01)
    assert rows[0]["angle"] == pytest.approx(90)
    assert rows[0]["y"] == pytest.approx(-18)

    kept, rows, _grids = partition_holes(
        [_hole(-18), _hole(-6), _hole(6), _hole(18), _hole(0, 14)], 0.01)
    assert len(rows) == 1 and len(kept) == 1 and kept[0]["y"] == pytest.approx(14)

    blind = [_hole(i, 0, depth=4) for i in (0, 12, 24, 36)]
    kept, rows, grids = partition_holes(blind, 0.01)
    assert rows == [] and grids == [] and len(kept) == 4


def test_row_of_four_is_one_pattern():
    _ref, _res, code, q = _export(plate_row())
    assert q["max_deviation"] <= 0.1
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["hole_row_1_n"] == pytest.approx(4)
    assert params["hole_row_1_pitch"] == pytest.approx(12)
    assert params["hole_row_1_d"] == pytest.approx(6, abs=0.05)
    assert params["hole_row_1_x"] == pytest.approx(-18, abs=0.05)
    assert params["hole_row_1_y"] == pytest.approx(0, abs=0.05)
    assert params["hole_row_1_angle"] == pytest.approx(0, abs=0.05)
    assert "hole_1_d" not in params
    assert "for _i in range(int(round(hole_row_1_n))):" in code
    solid = run_script(code)
    assert solid.is_watertight
    # four cylinders through a 6 mm plate
    block = 100 * 40 * 6
    assert abs(solid.volume) == pytest.approx(block - 4 * math.pi * 9 * 6, rel=0.01)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_row_survives_005_mm_noise(seed):
    ref, res, code, q = _export(plate_row(), sigma=0.05, seed=seed)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["hole_row_1_n"] == pytest.approx(4)
    assert params["hole_row_1_pitch"] == pytest.approx(12, abs=0.2)
    assert params["hole_row_1_d"] == pytest.approx(6, abs=0.1)
    assert "hole_1_d" not in params
    against = check(code, ref)
    assert against["max_deviation"] < 0.25
    assert against["volume_error_pct"] < 0.5


def test_grid_is_one_pattern():
    _ref, _res, code, q = _export(plate_grid())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert params["hole_grid_1_nx"] == pytest.approx(3)
    assert params["hole_grid_1_ny"] == pytest.approx(2)
    assert params["hole_grid_1_pitch_x"] == pytest.approx(20, abs=0.05)
    assert params["hole_grid_1_pitch_y"] == pytest.approx(24, abs=0.05)
    assert params["hole_grid_1_d"] == pytest.approx(5, abs=0.05)
    assert params["hole_grid_1_x"] == pytest.approx(-20, abs=0.05)
    assert params["hole_grid_1_y"] == pytest.approx(-12, abs=0.05)
    assert params["hole_grid_1_angle"] == pytest.approx(0, abs=0.05)
    assert "hole_row_1_n" not in params
    assert "hole_1_d" not in params
    assert "for _j in range(int(round(hole_grid_1_ny))):" in code
    solid = run_script(code)
    block = 100 * 70 * 6
    assert abs(solid.volume) == pytest.approx(block - 6 * math.pi * 2.5 ** 2 * 6, rel=0.01)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_grid_survives_005_mm_noise(seed):
    ref, res, code, q = _export(plate_grid(), sigma=0.05, seed=seed)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    params = _params(code)
    assert params["hole_grid_1_nx"] == pytest.approx(3)
    assert params["hole_grid_1_ny"] == pytest.approx(2)
    assert params["hole_grid_1_pitch_x"] == pytest.approx(20, abs=0.25)
    assert params["hole_grid_1_pitch_y"] == pytest.approx(24, abs=0.25)
    against = check(code, ref)
    assert against["max_deviation"] < 0.3


def test_near_miss_stays_individual_holes():
    _ref, _res, code, q = _export(plate_near_miss())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert "hole_row_1_n" not in params
    assert "hole_grid_1_nx" not in params
    xs = sorted(params[f"hole_{i}_x"] for i in range(1, 5))
    assert xs == pytest.approx([-18, -6, 8, 18], abs=0.05)
    assert all(params[f"hole_{i}_d"] == pytest.approx(6, abs=0.05) for i in range(1, 5))


def test_mild_spacing_error_falls_back_to_individual_holes():
    _ref, _res, code, q = _export(plate_mild_miss())
    assert q["max_deviation"] <= 0.1
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert "hole_row_1_n" not in params
    xs = sorted(params[f"hole_{i}_x"] for i in range(1, 5))
    assert xs == pytest.approx([-18, -6, 6, 18.6], abs=0.08)


def test_mixed_row_and_one_stray_hole():
    _ref, _res, code, q = _export(plate_mixed())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert params["hole_row_1_n"] == pytest.approx(4)
    assert params["hole_row_1_pitch"] == pytest.approx(12, abs=0.05)
    assert params["hole_1_d"] == pytest.approx(6, abs=0.05)
    assert params["hole_1_y"] == pytest.approx(14, abs=0.1)
    assert "hole_2_d" not in params


def test_vertical_row():
    _ref, _res, code, q = _export(plate_vertical())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert params["hole_row_1_n"] == pytest.approx(3)
    assert params["hole_row_1_pitch"] == pytest.approx(16, abs=0.05)
    assert params["hole_row_1_angle"] == pytest.approx(90, abs=0.5)
    assert params["hole_row_1_y"] == pytest.approx(-16, abs=0.1)


def test_edit_count_and_pitch(tmp_path):
    _ref, _res, code, _q = _export(plate_row())
    script = tmp_path / "row.py"
    script.write_text(code, encoding="utf-8")
    listed = list_parameters(str(script))
    by = {p["name"]: p for p in listed["parameters"]}
    assert by["hole_row_1_n"]["kind"] == "count" and by["hole_row_1_n"]["unit"] is None
    assert by["hole_row_1_pitch"]["kind"] == "length" and by["hole_row_1_pitch"]["unit"] == "mm"
    assert by["hole_row_1_d"]["kind"] == "diameter"
    assert by["hole_row_1_x"]["kind"] == "position"
    assert by["hole_row_1_angle"]["kind"] == "angle" and by["hole_row_1_angle"]["unit"] == "deg"
    assert by["hole_row_1_n"]["value"] == pytest.approx(4)

    five = edit_script(str(script), {"hole_row_1_n": 5}, output_path=str(tmp_path / "five.step"),
                       output_script=str(tmp_path / "five.py"))
    assert five["valid_solid"] is True
    assert five["rendered"] is True
    expect = math.pi * 3.0 ** 2 * 6.0
    assert five["before"]["volume"] - five["after"]["volume"] == pytest.approx(expect, rel=0.02)
    assert re.search(r"^hole_row_1_n\s*=\s*5\.0", (tmp_path / "five.py").read_text(encoding="utf-8"), re.M)
    # the script on disk is still the original
    assert re.search(r"^hole_row_1_n\s*=\s*4\.0", script.read_text(encoding="utf-8"), re.M)

    spread = edit_script(str(script), {"hole_row_1_pitch": "+5"}, output_path=str(tmp_path / "spread.step"),
                         output_script=str(tmp_path / "spread.py"))
    assert spread["valid_solid"] is True
    assert spread["delta"]["volume"] == pytest.approx(0, abs=2.0)
    assert re.search(r"^hole_row_1_pitch\s*=\s*17\.0", (tmp_path / "spread.py").read_text(encoding="utf-8"), re.M)
    spread_solid = run_script((tmp_path / "spread.py").read_text(encoding="utf-8"))
    assert spread_solid.is_watertight

    with pytest.raises(GeometryError, match="hole_row_1_d"):
        edit_script(str(script), {"hole_row_1_n": 8}, in_place=True)
    assert re.search(r"^hole_row_1_n\s*=\s*4\.0", script.read_text(encoding="utf-8"), re.M)

    assert main(["params", str(script)]) == 0
    assert main(["edit", str(script), "--set", "hole_row_1_n=5",
                 "--output-script", str(tmp_path / "cli.py"), "-o", str(tmp_path / "cli.step")]) == 0
    cli = (tmp_path / "cli.py").read_text(encoding="utf-8")
    assert re.search(r"^hole_row_1_n\s*=\s*5\.0", cli, re.M)
    assert run_script(cli).is_watertight

    pytest.importorskip("mcp")
    from mesh2cad.mcp_server import edit_params, list_params

    names = list_params(str(script))
    assert any(p["name"] == "hole_row_1_n" and p["kind"] == "count" for p in names["parameters"])
    copy = tmp_path / "mcp.py"
    copy.write_text(script.read_text(encoding="utf-8"), encoding="utf-8")
    result = edit_params(str(copy), {"hole_row_1_n": 5}, str(tmp_path / "mcp.step"))
    assert "error" not in result
    assert result["valid_solid"] is True
    text = copy.read_text(encoding="utf-8")
    assert re.search(r"^hole_row_1_n\s*=\s*5\.0", text, re.M)


def test_other_backends_use_the_same_row_and_grid_names():
    ref = to_mesh(plate_row())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    scad, _q = to_openscad(res.mesh, res.regions, res.tol)
    macro, _fq = to_freecad(res.mesh, res.regions, res.tol)
    for code in (scad, macro):
        params = _params(code)
        assert params["hole_row_1_n"] == pytest.approx(4)
        assert params["hole_row_1_pitch"] == pytest.approx(12, abs=0.05)
        assert "hole_1_d" not in params
    assert "for (_i = [0 : hole_row_1_n - 1])" in scad
    assert "cylinder(h=base_t + 2, d=hole_row_1_d);" in scad
    assert "for _i in range(int(round(hole_row_1_n))):" in macro
    assert "Part.makeCylinder(hole_row_1_d / 2, base_t + 2," in macro

    ref = to_mesh(plate_grid())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    scad, _q = to_openscad(res.mesh, res.regions, res.tol)
    macro, _fq = to_freecad(res.mesh, res.regions, res.tol)
    for code in (scad, macro):
        params = _params(code)
        assert params["hole_grid_1_nx"] == pytest.approx(3)
        assert params["hole_grid_1_ny"] == pytest.approx(2)
        assert params["hole_grid_1_pitch_x"] == pytest.approx(20, abs=0.05)
    assert "for (_j = [0 : hole_grid_1_ny - 1])" in scad
    assert "for _j in range(int(round(hole_grid_1_ny))):" in macro


@needs_openscad
def test_openscad_row_edit_renders():
    ref = to_mesh(plate_row())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q is not None
    assert q["max_deviation"] <= 0.1
    bigger = re.sub(r"^(hole_row_1_n\s*=\s*)[\d.]+", r"\g<1>5.0", code, count=1, flags=re.M)
    from mesh2cad.openscad_code import render_mesh
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        out = render_mesh(bigger, Path(d))
    assert out.is_watertight
    base = abs(run_script(to_cadquery(res.mesh, res.regions, res.tol)[0]).volume)
    # OpenSCAD volume is a $fn mesh. One extra hole is still about one cylinder.
    assert base - abs(out.volume) == pytest.approx(math.pi * 9 * 6, rel=0.05)


def test_filleted_box_row():
    # Default subdivision of every filleted edge is a couple of hundred thousand
    # triangles. The coarser mesh is the same kind the rounded-box tests use.
    ref = to_mesh(plate_filleted_row(), tol=0.1, ang=0.4, max_edge=4.0)
    res = convert_mesh(ref, tol=0.05, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    params = _params(code)
    assert "edge_fillet_r" in params
    assert params["hole_row_1_n"] == pytest.approx(4)
    assert params["hole_row_1_pitch"] == pytest.approx(16, abs=0.1)
    assert params["hole_row_1_d"] == pytest.approx(6, abs=0.1)
    assert ".edges().fillet(edge_fillet_r)" in code
    assert code.index(".edges().fillet(edge_fillet_r)") < code.index("for _i in range(int(round(hole_row_1_n))):")


def test_bracket_foot_row():
    _ref, _res, code, q = _export(bracket_foot_row())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert params["base_hole_row_1_n"] == pytest.approx(3)
    assert params["base_hole_row_1_pitch"] == pytest.approx(8, abs=0.1)
    assert params["base_hole_row_1_d"] == pytest.approx(4, abs=0.05)
    assert params["base_hole_row_1_x"] == pytest.approx(16, abs=0.15)
    assert params["base_hole_row_1_y"] == pytest.approx(4, abs=0.15)
    assert "base_hole_1_d" not in params
    assert "for _i in range(int(round(base_hole_row_1_n))):" in code
    fillet = code.index('.edges("|Z").fillet(fillet_r)')
    assert fillet < code.index("for _i in range(int(round(base_hole_row_1_n))):")
