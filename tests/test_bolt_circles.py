"""Bolt circles of plain Z through-holes as pattern parameters.

Fixtures are built in this process. Nothing is loaded from disk. A center
bore stays hole_N_*. A row or an axis-aligned grid keeps that pattern.
Five holes of a six-hole circle are not a circle. A pattern the deviation
gate rejects is written as one parameter set per hole.
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
from parts import to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")

PLATE = 80.0
THICK = 6.0
PCD_R = 28.0
HOLE_D = 6.0
BORE_D = 12.0


def _hole(x, y=0.0, d=HOLE_D, **extra):
    hole = {"x": float(x), "y": float(y), "d": float(d), "axis": "Z"}
    hole.update(extra)
    return hole


def _on_circle(count, radius=PCD_R, angle=0.0, diameter=HOLE_D, skip=()):
    holes = []
    for i in range(count):
        if i in skip:
            continue
        theta = math.radians(angle + i * 360.0 / count)
        holes.append(_hole(radius * math.cos(theta), radius * math.sin(theta), d=diameter))
    return holes


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _flange(points, diameter=HOLE_D, bore=BORE_D):
    part = cq.Workplane().box(PLATE, PLATE, THICK)
    part = part.faces(">Z").workplane().pushPoints([(0, 0)]).hole(bore)
    return part.faces(">Z").workplane().pushPoints(points).hole(diameter)


def plate_circle4():
    """80 x 80 x 6 plate, four Ø6 holes on a 28 mm radius, Ø12 center bore."""
    return _flange([(h["x"], h["y"]) for h in _on_circle(4)])


def plate_circle6():
    """Same plate, six Ø6 holes on the same radius, Ø12 center bore."""
    return _flange([(h["x"], h["y"]) for h in _on_circle(6)])


def plate_circle_angle():
    """Four Ø6 holes, start angle 20 degrees, Ø12 center bore."""
    return _flange([(h["x"], h["y"]) for h in _on_circle(4, angle=20.0)])


def plate_circle_missing():
    """Five of the six holes. One gap is double, so this is not a bolt circle."""
    return _flange([(h["x"], h["y"]) for h in _on_circle(6, skip=(0,))])


def plate_mild_miss():
    """Same four-hole circle, one center 0.2 mm long.

    The detector still calls it a circle. The snapped solid is farther from
    the mesh than the deviation gate, so the written script keeps one
    parameter set per hole.
    """
    points = [(h["x"], h["y"]) for h in _on_circle(4)]
    points[0] = (28.2, 0.0)
    return _flange(points)


def _export(wp, sigma=0.0, seed=0):
    ref = to_mesh(wp)
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=seed) if sigma else ref
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    return ref, res, code, q


def _circle_volume(count, hole_d=HOLE_D, bore_d=BORE_D):
    block = PLATE * PLATE * THICK
    hole = math.pi * (hole_d / 2.0) ** 2 * THICK
    bore = math.pi * (bore_d / 2.0) ** 2 * THICK
    return block - count * hole - bore


def test_detector_finds_a_circle_and_leaves_the_rest():
    # The center bore is a different diameter, so it is not part of the circle.
    holes = _on_circle(4) + [_hole(0, 0, d=BORE_D)]
    kept, rows, grids, circles = partition_holes(holes, 0.01)
    assert rows == [] and grids == []
    assert len(circles) == 1 and len(kept) == 1
    assert circles[0]["n"] == 4
    assert circles[0]["r"] == pytest.approx(PCD_R)
    assert circles[0]["cx"] == pytest.approx(0) and circles[0]["cy"] == pytest.approx(0)
    assert circles[0]["angle"] == pytest.approx(0)
    assert circles[0]["d"] == pytest.approx(HOLE_D)
    assert kept[0]["x"] == pytest.approx(0) and kept[0]["y"] == pytest.approx(0)
    assert kept[0]["d"] == pytest.approx(BORE_D)

    # A same-diameter center sits on the line through two opposite holes.
    # That line is a row, and rows are taken before circles.
    kept, rows, grids, circles = partition_holes(_on_circle(4) + [_hole(0, 0, d=HOLE_D)], 0.01)
    assert len(rows) == 1 and rows[0]["n"] == 3
    assert circles == [] and len(kept) == 2

    kept, rows, grids, circles = partition_holes(_on_circle(6), 0.01)
    assert rows == [] and grids == [] and kept == []
    assert circles[0]["n"] == 6 and circles[0]["angle"] == pytest.approx(0)
    assert circles[0]["r"] == pytest.approx(PCD_R)

    kept, rows, grids, circles = partition_holes(_on_circle(4, angle=20.0), 0.01)
    assert len(circles) == 1 and kept == []
    assert circles[0]["angle"] == pytest.approx(20)
    assert circles[0]["n"] == 4

    kept, rows, grids, circles = partition_holes(_on_circle(3, radius=20.0), 0.01)
    assert len(circles) == 1 and circles[0]["n"] == 3
    assert circles[0]["r"] == pytest.approx(20)
    assert circles[0]["angle"] == pytest.approx(0)

    # One of six is missing. The five that remain are not equally spaced,
    # and a regular subset is not split off while the others sit on the circle.
    kept, rows, grids, circles = partition_holes(_on_circle(6, skip=(0,)), 0.01)
    assert circles == [] and rows == [] and grids == [] and len(kept) == 5

    kept, rows, grids, circles = partition_holes(_on_circle(2), 0.01)
    assert circles == [] and len(kept) == 2

    # An axis-aligned square is a grid. It is not also written as a circle.
    square = [_hole(x, y) for x in (-10, 10) for y in (-10, 10)]
    kept, rows, grids, circles = partition_holes(square, 0.01)
    assert circles == [] and rows == [] and len(grids) == 1
    assert grids[0]["nx"] == 2 and grids[0]["ny"] == 2

    kept, rows, grids, circles = partition_holes(
        [_hole(-18), _hole(-6), _hole(6), _hole(18)], 0.01)
    assert circles == [] and len(rows) == 1

    # A different diameter, off the pitch circle, stays an individual hole.
    mixed = _on_circle(4) + [_hole(20, 0, d=10)]
    kept, rows, grids, circles = partition_holes(mixed, 0.01)
    assert len(circles) == 1 and circles[0]["d"] == pytest.approx(HOLE_D)
    assert len(kept) == 1 and kept[0]["d"] == pytest.approx(10)

    blind = [_hole(h["x"], h["y"], depth=4) for h in _on_circle(4)]
    kept, rows, grids, circles = partition_holes(blind, 0.01)
    assert circles == [] and len(kept) == 4


def test_four_holes_on_a_circle_are_one_pattern():
    _ref, _res, code, q = _export(plate_circle4())
    assert q["max_deviation"] <= 0.1
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["hole_circle_1_n"] == pytest.approx(4)
    assert params["hole_circle_1_d"] == pytest.approx(HOLE_D, abs=0.05)
    assert params["hole_circle_1_r"] == pytest.approx(PCD_R, abs=0.05)
    assert params["hole_circle_1_cx"] == pytest.approx(0, abs=0.05)
    assert params["hole_circle_1_cy"] == pytest.approx(0, abs=0.05)
    assert params["hole_circle_1_angle"] == pytest.approx(0, abs=0.05)
    assert params["hole_1_d"] == pytest.approx(BORE_D, abs=0.05)
    assert params["hole_1_x"] == pytest.approx(0, abs=0.05)
    assert params["hole_1_y"] == pytest.approx(0, abs=0.05)
    assert "hole_2_d" not in params
    assert "hole_row_1_n" not in params
    assert "hole_grid_1_nx" not in params
    assert "for _i in range(int(round(hole_circle_1_n))):" in code
    assert "hole_circle_1_r * math.cos" in code
    solid = run_script(code)
    assert solid.is_watertight
    assert abs(solid.volume) == pytest.approx(_circle_volume(4), rel=0.01)


def test_six_holes_and_a_turned_circle():
    _ref, _res, code, q = _export(plate_circle6())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert params["hole_circle_1_n"] == pytest.approx(6)
    assert params["hole_circle_1_r"] == pytest.approx(PCD_R, abs=0.05)
    assert params["hole_circle_1_angle"] == pytest.approx(0, abs=0.5)
    assert params["hole_1_d"] == pytest.approx(BORE_D, abs=0.05)
    assert "hole_2_d" not in params
    solid = run_script(code)
    assert abs(solid.volume) == pytest.approx(_circle_volume(6), rel=0.01)

    _ref, _res, code, q = _export(plate_circle_angle())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert params["hole_circle_1_n"] == pytest.approx(4)
    assert params["hole_circle_1_angle"] == pytest.approx(20, abs=0.5)
    assert params["hole_circle_1_r"] == pytest.approx(PCD_R, abs=0.05)
    assert params["hole_1_d"] == pytest.approx(BORE_D, abs=0.05)


def test_mild_radius_error_falls_back_to_individual_holes():
    _ref, _res, code, q = _export(plate_mild_miss())
    assert q["max_deviation"] <= 0.1
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert "hole_circle_1_n" not in params
    xs = sorted(params[f"hole_{i}_x"] for i in range(1, 6))
    assert xs == pytest.approx([-28, 0, 0, 0, 28.2], abs=0.08)


def test_missing_hole_stays_individual():
    _ref, _res, code, q = _export(plate_circle_missing())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert "hole_circle_1_n" not in params
    assert "hole_row_1_n" not in params
    assert "hole_grid_1_nx" not in params
    diameters = [params[f"hole_{i}_d"] for i in range(1, 7)]
    assert sorted(diameters) == pytest.approx([HOLE_D] * 5 + [BORE_D], abs=0.05)


@pytest.mark.parametrize("sigma", [0.05, 0.2])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_circle_survives_noise(sigma, seed):
    ref, res, code, q = _export(plate_circle4(), sigma=sigma, seed=seed)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["hole_circle_1_n"] == pytest.approx(4)
    assert params["hole_circle_1_r"] == pytest.approx(PCD_R, abs=0.5)
    assert params["hole_circle_1_d"] == pytest.approx(HOLE_D, abs=0.25)
    assert params["hole_circle_1_angle"] == pytest.approx(0, abs=2.0)
    assert "hole_row_1_n" not in params
    against = check(code, ref)
    # The noisy mesh is the one the gate used. Against the clean plate the
    # script can sit about as far away as the noise, plus a little tessellation.
    assert against["max_deviation"] < 0.35 + sigma
    assert against["volume_error_pct"] < 1.0


def test_edit_count_radius_and_refuses_a_bad_edit(tmp_path):
    _ref, _res, code, _q = _export(plate_circle4())
    script = tmp_path / "circle.py"
    script.write_text(code, encoding="utf-8")
    listed = list_parameters(str(script))
    by = {p["name"]: p for p in listed["parameters"]}
    assert by["hole_circle_1_n"]["kind"] == "count" and by["hole_circle_1_n"]["unit"] is None
    assert by["hole_circle_1_d"]["kind"] == "diameter" and by["hole_circle_1_d"]["unit"] == "mm"
    assert by["hole_circle_1_r"]["kind"] == "radius" and by["hole_circle_1_r"]["unit"] == "mm"
    assert by["hole_circle_1_cx"]["kind"] == "position"
    assert by["hole_circle_1_cy"]["kind"] == "position"
    assert by["hole_circle_1_angle"]["kind"] == "angle" and by["hole_circle_1_angle"]["unit"] == "deg"
    assert by["hole_circle_1_n"]["value"] == pytest.approx(4)

    eight = edit_script(str(script), {"hole_circle_1_n": 8}, output_path=str(tmp_path / "eight.step"),
                        output_script=str(tmp_path / "eight.py"))
    assert eight["valid_solid"] is True
    assert eight["rendered"] is True
    expect = 4.0 * math.pi * (HOLE_D / 2.0) ** 2 * THICK
    assert eight["before"]["volume"] - eight["after"]["volume"] == pytest.approx(expect, rel=0.02)
    assert re.search(r"^hole_circle_1_n\s*=\s*8\.0", (tmp_path / "eight.py").read_text(encoding="utf-8"), re.M)
    assert re.search(r"^hole_circle_1_n\s*=\s*4\.0", script.read_text(encoding="utf-8"), re.M)

    turned = edit_script(str(script), {"hole_circle_1_angle": 15}, output_path=str(tmp_path / "turn.step"),
                         output_script=str(tmp_path / "turn.py"))
    assert turned["valid_solid"] is True
    assert re.search(r"^hole_circle_1_angle\s*=\s*15\.0", (tmp_path / "turn.py").read_text(encoding="utf-8"), re.M)

    raw = script.read_bytes()
    with pytest.raises(GeometryError, match="whole number"):
        edit_script(str(script), {"hole_circle_1_n": 4.5}, in_place=True)
    assert script.read_bytes() == raw

    with pytest.raises(GeometryError, match="does not fit"):
        edit_script(str(script), {"hole_circle_1_r": 42}, in_place=True)
    assert script.read_bytes() == raw

    with pytest.raises(GeometryError, match="overlaps another hole of the same pattern"):
        edit_script(str(script), {"hole_circle_1_n": 36}, in_place=True)
    assert script.read_bytes() == raw

    with pytest.raises(GeometryError, match="overlaps"):
        edit_script(str(script), {"hole_circle_1_r": 8}, in_place=True)
    assert script.read_bytes() == raw

    assert main(["params", str(script)]) == 0
    assert main(["edit", str(script), "--set", "hole_circle_1_n=8",
                 "--output-script", str(tmp_path / "cli.py"), "-o", str(tmp_path / "cli.step")]) == 0
    cli = (tmp_path / "cli.py").read_text(encoding="utf-8")
    assert re.search(r"^hole_circle_1_n\s*=\s*8\.0", cli, re.M)
    assert run_script(cli).is_watertight

    pytest.importorskip("mcp")
    from mesh2cad.mcp_server import edit_params, list_params

    names = list_params(str(script))
    assert any(p["name"] == "hole_circle_1_n" and p["kind"] == "count" for p in names["parameters"])
    assert any(p["name"] == "hole_circle_1_r" and p["kind"] == "radius" for p in names["parameters"])
    copy = tmp_path / "mcp.py"
    copy.write_text(script.read_text(encoding="utf-8"), encoding="utf-8")
    result = edit_params(str(copy), {"hole_circle_1_n": 8}, str(tmp_path / "mcp.step"))
    assert "error" not in result
    assert result["valid_solid"] is True
    text = copy.read_text(encoding="utf-8")
    assert re.search(r"^hole_circle_1_n\s*=\s*8\.0", text, re.M)


def test_other_backends_use_the_same_circle_names():
    ref = to_mesh(plate_circle4())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    scad, _q = to_openscad(res.mesh, res.regions, res.tol)
    macro, _fq = to_freecad(res.mesh, res.regions, res.tol)
    for code in (scad, macro):
        params = _params(code)
        assert params["hole_circle_1_n"] == pytest.approx(4)
        assert params["hole_circle_1_r"] == pytest.approx(PCD_R, abs=0.05)
        assert params["hole_circle_1_angle"] == pytest.approx(0, abs=0.05)
        assert params["hole_1_d"] == pytest.approx(BORE_D, abs=0.05)
        assert "hole_2_d" not in params
    assert "for (_i = [0 : hole_circle_1_n - 1])" in scad
    assert "cylinder(h=base_t + 2, d=hole_circle_1_d);" in scad
    assert "cos(hole_circle_1_angle + _i * 360 / hole_circle_1_n)" in scad
    assert "for _i in range(int(round(hole_circle_1_n))):" in macro
    assert "Part.makeCylinder(hole_circle_1_d / 2, base_t + 2," in macro


@needs_openscad
def test_openscad_circle_edit_renders():
    ref = to_mesh(plate_circle4())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q is not None
    assert q["max_deviation"] <= 0.1
    bigger = re.sub(r"^(hole_circle_1_n\s*=\s*)[\d.]+", r"\g<1>8.0", code, count=1, flags=re.M)
    from mesh2cad.openscad_code import render_mesh
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        out = render_mesh(bigger, Path(d))
    assert out.is_watertight
    base = abs(run_script(to_cadquery(res.mesh, res.regions, res.tol)[0]).volume)
    extra = 4.0 * math.pi * (HOLE_D / 2.0) ** 2 * THICK
    assert base - abs(out.volume) == pytest.approx(extra, rel=0.08)
