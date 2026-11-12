"""Obround through-slots on a plate, as slot_N_* parameters.

Fixtures are built here. A rectangular opening stays a cutout. A round hole
stays hole_N_*. A blind stadium is not written as a through-slot.
"""
import math
import re

import pytest

pytest.importorskip("cadquery")
import cadquery as cq  # noqa: E402

from mesh2cad.cadquery_code import run_script, to_cadquery  # noqa: E402
from mesh2cad.cli import main  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.freecad_code import to_freecad  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.openscad_code import to_openscad  # noqa: E402
from mesh2cad.parametric import Unsupported  # noqa: E402
from mesh2cad.params import GeometryError, edit_script, list_parameters  # noqa: E402
from parts import plate_with_hole, slotted_plate, to_mesh  # noqa: E402


def _obround(length, width, angle, x, y, depth, z0):
    """Cutter for one stadium. ``length`` is the centre distance."""
    radius = width / 2.0
    ang = math.radians(angle)
    dx = (length / 2.0) * math.cos(ang)
    dy = (length / 2.0) * math.sin(ang)
    c1 = cq.Solid.makeCylinder(radius, depth, cq.Vector(x - dx, y - dy, z0), cq.Vector(0, 0, 1))
    c2 = cq.Solid.makeCylinder(radius, depth, cq.Vector(x + dx, y + dy, z0), cq.Vector(0, 0, 1))
    ca, sa = math.cos(ang), math.sin(ang)
    px, py = -sa * radius, ca * radius
    pts = [
        (x - dx + px, y - dy + py),
        (x + dx + px, y + dy + py),
        (x + dx - px, y + dy - py),
        (x - dx - px, y - dy - py),
    ]
    bar = cq.Workplane("XY", origin=(0, 0, z0)).polyline(pts).close().extrude(depth).val()
    return c1, c2, bar


def _plate(slots, thick=6.0, blind=False):
    """60 x 40 x thick plate. Each slot is (len, width, angle, x, y)."""
    part = cq.Workplane().box(60, 40, thick).val()
    if blind:
        z0, depth = 0.0, thick / 2.0 + 1.0
    else:
        z0, depth = -thick / 2.0 - 2.0, thick + 4.0
    for spec in slots:
        for solid in _obround(*spec, depth, z0):
            part = part.cut(solid)
    return cq.Workplane(obj=part)


def plate_one_slot():
    """60 x 40 x 6, one horizontal slot, centre distance 20, width 8, at the origin."""
    return _plate([(20.0, 8.0, 0.0, 0.0, 0.0)])


def plate_two_slots():
    """One slot along X and one along Y, different sizes, not overlapping."""
    return _plate([
        (16.0, 6.0, 0.0, -10.0, 8.0),
        (14.0, 5.0, 90.0, 12.0, -6.0),
    ])


def plate_angled_slot():
    """Same plate, one slot at 30 degrees."""
    return _plate([(18.0, 6.0, 30.0, 0.0, 0.0)])


def plate_blind_slot():
    """The stadium opens on the top face and stops halfway through."""
    return _plate([(16.0, 6.0, 0.0, 0.0, 0.0)], blind=True)


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _export(wp, sigma=0.0, seed=0):
    ref = to_mesh(wp)
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=seed) if sigma else ref
    res = convert_mesh(mesh, tol=tol, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    return ref, res, code, q


def _slot_area(length, width):
    radius = width / 2.0
    return length * width + math.pi * radius * radius


def test_one_axis_aligned_slot_is_named():
    _ref, _res, code, q = _export(plate_one_slot())
    assert q["max_deviation"] <= 0.1
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["slot_1_len"] == pytest.approx(20.0)
    assert params["slot_1_w"] == pytest.approx(8.0)
    assert params["slot_1_x"] == pytest.approx(0.0, abs=0.05)
    assert params["slot_1_y"] == pytest.approx(0.0, abs=0.05)
    assert params["slot_1_angle"] == pytest.approx(0.0)
    assert "hole_1_d" not in params
    assert "slot2D(slot_1_len + slot_1_w, slot_1_w, slot_1_angle)" in code
    solid = run_script(code)
    assert solid.is_watertight
    removed = _slot_area(20.0, 8.0) * 6.0
    assert abs(solid.volume) == pytest.approx(60 * 40 * 6 - removed, rel=0.002)


def test_vertical_slot_is_ninety_degrees():
    _ref, _res, code, q = _export(_plate([(16.0, 6.0, 90.0, -6.0, 2.0)]))
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert params["slot_1_len"] == pytest.approx(16.0)
    assert params["slot_1_w"] == pytest.approx(6.0)
    assert params["slot_1_x"] == pytest.approx(-6.0, abs=0.1)
    assert params["slot_1_y"] == pytest.approx(2.0, abs=0.1)
    assert params["slot_1_angle"] == pytest.approx(90.0)


def test_two_slots_of_different_orientation():
    _ref, _res, code, q = _export(plate_two_slots())
    assert q["max_deviation"] <= 0.1
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    found = []
    for i in (1, 2):
        found.append((
            params[f"slot_{i}_angle"],
            params[f"slot_{i}_len"],
            params[f"slot_{i}_w"],
            params[f"slot_{i}_x"],
            params[f"slot_{i}_y"],
        ))
    along_x = [row for row in found if row[0] == pytest.approx(0.0)]
    along_y = [row for row in found if row[0] == pytest.approx(90.0)]
    assert len(along_x) == 1 and len(along_y) == 1
    assert along_x[0][1] == pytest.approx(16.0)
    assert along_x[0][2] == pytest.approx(6.0)
    assert along_x[0][3] == pytest.approx(-10.0, abs=0.1)
    assert along_x[0][4] == pytest.approx(8.0, abs=0.1)
    assert along_y[0][1] == pytest.approx(14.0)
    assert along_y[0][2] == pytest.approx(5.0)
    assert along_y[0][3] == pytest.approx(12.0, abs=0.1)
    assert along_y[0][4] == pytest.approx(-6.0, abs=0.1)
    assert "slot_3_len" not in params
    solid = run_script(code)
    removed = (_slot_area(16.0, 6.0) + _slot_area(14.0, 5.0)) * 6.0
    assert abs(solid.volume) == pytest.approx(60 * 40 * 6 - removed, rel=0.002)


def test_angled_slot_keeps_its_angle():
    _ref, _res, code, q = _export(plate_angled_slot())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert params["slot_1_len"] == pytest.approx(18.0, abs=0.1)
    assert params["slot_1_w"] == pytest.approx(6.0, abs=0.1)
    assert params["slot_1_angle"] == pytest.approx(30.0, abs=0.5)
    assert params["slot_1_x"] == pytest.approx(0.0, abs=0.15)
    assert params["slot_1_y"] == pytest.approx(0.0, abs=0.15)


def test_other_backends_use_the_same_slot_names():
    ref = to_mesh(plate_one_slot())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    scad, sq = to_openscad(res.mesh, res.regions, res.tol)
    macro, _fq = to_freecad(res.mesh, res.regions, res.tol)
    for code in (scad, macro):
        params = _params(code)
        assert params["slot_1_len"] == pytest.approx(20.0)
        assert params["slot_1_w"] == pytest.approx(8.0)
        assert params["slot_1_angle"] == pytest.approx(0.0)
        assert "hole_1_d" not in params
    assert "/* [Slots] */" in scad
    assert "hull()" in scad
    assert "cylinder(h=base_t + 2, d=slot_1_w);" in scad
    assert sq is not None and sq["max_deviation"] <= 0.1
    assert "def _slot(length, width, angle, x, y, z0, height):" in macro
    assert "result = result.cut(_slot(slot_1_len, slot_1_w, slot_1_angle," in macro


def test_edit_lengthens_the_slot_and_refuses_a_breakout(tmp_path):
    _ref, _res, code, _q = _export(plate_one_slot())
    script = tmp_path / "slot.py"
    script.write_text(code, encoding="utf-8")
    listed = list_parameters(str(script))
    by = {p["name"]: p for p in listed["parameters"]}
    assert by["slot_1_len"]["kind"] == "length" and by["slot_1_len"]["unit"] == "mm"
    assert by["slot_1_w"]["kind"] == "diameter" and by["slot_1_w"]["unit"] == "mm"
    assert by["slot_1_x"]["kind"] == "position"
    assert by["slot_1_angle"]["kind"] == "angle" and by["slot_1_angle"]["unit"] == "deg"
    assert by["slot_1_len"]["value"] == pytest.approx(20.0)

    longer = edit_script(
        str(script), {"slot_1_len": 22},
        output_path=str(tmp_path / "longer.step"),
        output_script=str(tmp_path / "longer.py"),
    )
    assert longer["valid_solid"] is True
    assert longer["rendered"] is True
    # Two millimetres of straight section, the round ends stay put.
    assert longer["before"]["volume"] - longer["after"]["volume"] == pytest.approx(2 * 8 * 6, rel=0.02)
    text = (tmp_path / "longer.py").read_text(encoding="utf-8")
    assert re.search(r"^slot_1_len\s*=\s*22\.0", text, re.M)
    assert re.search(r"^slot_1_len\s*=\s*20\.0", script.read_text(encoding="utf-8"), re.M)
    assert run_script(text).is_watertight

    with pytest.raises(GeometryError, match="slot_1_w"):
        edit_script(str(script), {"slot_1_len": 70}, in_place=True)
    assert re.search(r"^slot_1_len\s*=\s*20\.0", script.read_text(encoding="utf-8"), re.M)

    with pytest.raises(GeometryError, match="slot_1_w"):
        edit_script(str(script), {"slot_1_w": 42})

    assert main(["params", str(script)]) == 0
    assert main(["edit", str(script), "--set", "slot_1_len=22",
                 "--output-script", str(tmp_path / "cli.py"),
                 "-o", str(tmp_path / "cli.step")]) == 0
    cli = (tmp_path / "cli.py").read_text(encoding="utf-8")
    assert re.search(r"^slot_1_len\s*=\s*22\.0", cli, re.M)

    pytest.importorskip("mcp")
    from mesh2cad.mcp_server import edit_params, list_params

    names = list_params(str(script))
    assert any(p["name"] == "slot_1_len" and p["kind"] == "length" for p in names["parameters"])
    copy = tmp_path / "mcp.py"
    copy.write_text(script.read_text(encoding="utf-8"), encoding="utf-8")
    result = edit_params(str(copy), {"slot_1_len": "+2"}, str(tmp_path / "mcp.step"))
    assert "error" not in result
    assert result["valid_solid"] is True
    assert result["before"]["volume"] - result["after"]["volume"] == pytest.approx(2 * 8 * 6, rel=0.02)
    assert re.search(r"^slot_1_len\s*=\s*22\.0", copy.read_text(encoding="utf-8"), re.M)


def test_round_hole_and_rectangular_slot_are_not_obrounds():
    _ref, _res, code, q = _export(plate_with_hole())
    assert q["max_deviation"] <= 0.1
    params = _params(code)
    assert "hole_1_d" in params
    assert not any(name.startswith("slot_") for name in params)

    ref = to_mesh(slotted_plate())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= 0.1
    assert "cutout_1" in code
    assert "slot_1_len" not in code
    assert "slot2D" not in code


def test_blind_slot_is_not_forced():
    ref = to_mesh(plate_blind_slot())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    with pytest.raises(Unsupported) as raised:
        to_cadquery(res.mesh, res.regions, res.tol)
    message = str(raised.value)
    assert "slot_" not in message
    assert message


def test_unequal_ends_are_not_one_slot():
    """A stadium whose ends have different diameters is not averaged into slot_N_w."""
    part = cq.Workplane().box(60, 40, 6).val()
    z0, depth = -5.0, 10.0
    # Centres 16 mm apart. Left end Ø6, right end Ø10. The walls are not parallel.
    left = cq.Solid.makeCylinder(3.0, depth, cq.Vector(-8.0, 0.0, z0), cq.Vector(0, 0, 1))
    right = cq.Solid.makeCylinder(5.0, depth, cq.Vector(8.0, 0.0, z0), cq.Vector(0, 0, 1))
    bar = (cq.Workplane("XY", origin=(0, 0, z0))
           .polyline([(-8, 3), (8, 5), (8, -5), (-8, -3)]).close().extrude(depth).val())
    part = part.cut(left).cut(right).cut(bar)
    ref = to_mesh(cq.Workplane(obj=part))
    res = convert_mesh(ref, tol=0.01, with_report=False)
    with pytest.raises(Unsupported) as raised:
        to_cadquery(res.mesh, res.regions, res.tol)
    assert "slot_" not in str(raised.value)


@pytest.mark.parametrize("sigma", [0.05, 0.2])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_one_slot_survives_noise(sigma, seed):
    ref, res, code, q = _export(plate_one_slot(), sigma=sigma, seed=seed)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["slot_1_len"] == pytest.approx(20.0, abs=0.25)
    assert params["slot_1_w"] == pytest.approx(8.0, abs=0.25)
    assert params["slot_1_angle"] == pytest.approx(0.0, abs=1.0)
    assert "hole_1_d" not in params
    solid = run_script(code)
    assert abs(solid.volume) == pytest.approx(abs(ref.volume), rel=0.01)
