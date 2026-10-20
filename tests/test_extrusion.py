"""Parametric export of a straight extrusion: outline, length, holes."""
import math
import re
import shutil

import pytest

pytest.importorskip("cadquery")
from mesh2cad.cadquery_code import check as cq_check  # noqa: E402
from mesh2cad.cadquery_code import run_script, to_cadquery  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.openscad_code import check as scad_check  # noqa: E402
from mesh2cad.openscad_code import render_mesh, to_openscad  # noqa: E402
from parts import l_profile_y, slanted_outline, slotted_plate, to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")

# L in XZ, extruded along Y. Area = 14*28 + 22*12.
L_OUTLINE = [(0.0, 0.0), (14.0, 0.0), (14.0, 16.0), (36.0, 16.0), (36.0, 28.0), (0.0, 28.0)]
L_AREA = 656.0
# 50×32 with the corner past (50, 20)–(34, 32) cut off. Area = 50*32 - 16*12/2.
SLANT_OUTLINE = [(0.0, 0.0), (50.0, 0.0), (50.0, 20.0), (34.0, 32.0), (0.0, 32.0)]
SLANT_AREA = 1504.0
SLOT_OUTLINE = [(0.0, 0.0), (40.0, 0.0), (40.0, 24.0), (0.0, 24.0)]
SLOT_CUTOUT = [(8.0, 8.0), (24.0, 8.0), (24.0, 16.0), (8.0, 16.0)]


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _point_block(code, name):
    m = re.search(rf"{name}\s*=\s*\[", code)
    assert m, name
    depth = 1
    j = m.end()
    while j < len(code) and depth:
        depth += (code[j] == "[") - (code[j] == "]")
        j += 1
    nums = [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?", code[m.end():j - 1])]
    return list(zip(nums[0::2], nums[1::2]))


def _same_loop(got, expect, tol):
    assert len(got) == len(expect)
    for px, py in expect:
        assert any(abs(px - gx) <= tol and abs(py - gy) <= tol for gx, gy in got)


def _fit(fn, sigma):
    ref = to_mesh(fn())
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=0) if sigma else ref
    res = convert_mesh(mesh, tol=tol, with_report=False)
    return ref, res


@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("fn,length,diameter,axis,outline", [
    (l_profile_y, 18.0, 6.0, "Y", L_OUTLINE),
    (slanted_outline, 12.0, 8.0, "Z", SLANT_OUTLINE),
])
def test_cadquery_extrusion(fn, length, diameter, axis, outline, sigma):
    ref, res = _fit(fn, sigma)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    limit = max(4 * res.tol, 0.1)
    assert q["max_deviation"] <= limit
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["extrude_len"] == pytest.approx(length, abs=0.05)
    assert params["hole_1_d"] == pytest.approx(diameter, abs=0.05)
    assert f"length along {axis}" in code
    assert "polyline(outline)" in code
    _same_loop(_point_block(code, "outline"), outline, 0.05 if sigma == 0.0 else 0.25)
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)
    assert against["volume_error_pct"] < 0.5


@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_rectangular_cutout(sigma):
    ref, res = _fit(slotted_plate, sigma)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    assert _params(code)["extrude_len"] == pytest.approx(6.0, abs=0.05)
    assert "hole_1_d" not in code
    tol = 0.05 if sigma == 0.0 else 0.25
    _same_loop(_point_block(code, "outline"), SLOT_OUTLINE, tol)
    _same_loop(_point_block(code, "cutout_1"), SLOT_CUTOUT, tol)
    assert "polyline(cutout_1)" in code
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)


def test_edit_hole_diameter_and_extrusion_length():
    """One hole diameter and one extrusion length, each checked against the analytic delta."""
    ref, res = _fit(l_profile_y, 0.0)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    v0 = abs(run_script(code).volume)
    wider = re.sub(r"^hole_1_d\s*=\s*[\d.]+", "hole_1_d = 8.0", code, flags=re.M)
    # Ø6 → Ø8 through the 36 mm foot, perpendicular to the extrusion
    expect_hole = math.pi * (4.0 ** 2 - 3.0 ** 2) * 36.0
    assert v0 - abs(run_script(wider).volume) == pytest.approx(expect_hole, rel=0.01)
    longer = re.sub(r"^extrude_len\s*=\s*[\d.]+", "extrude_len = 22.0", code, flags=re.M)
    assert abs(run_script(longer).volume) - v0 == pytest.approx(L_AREA * 4.0, rel=0.01)

    ref, res = _fit(slanted_outline, 0.0)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    v0 = abs(run_script(code).volume)
    wider = re.sub(r"^hole_1_d\s*=\s*[\d.]+", "hole_1_d = 10.0", code, flags=re.M)
    # Ø8 → Ø10 through the 12 mm extrusion
    expect_hole = math.pi * (5.0 ** 2 - 4.0 ** 2) * 12.0
    assert v0 - abs(run_script(wider).volume) == pytest.approx(expect_hole, rel=0.01)
    longer = re.sub(r"^extrude_len\s*=\s*[\d.]+", "extrude_len = 16.0", code, flags=re.M)
    # the hole runs along the extrusion, so it grows with the length
    expect_len = (SLANT_AREA - math.pi * 4.0 ** 2) * 4.0
    assert abs(run_script(longer).volume) - v0 == pytest.approx(expect_len, rel=0.01)


def test_exporters_share_extrusion_parameters():
    for fn in (l_profile_y, slanted_outline, slotted_plate):
        _ref, res = _fit(fn, 0.0)
        cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
        sc, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
        assert _params(sc) == _params(cq)
        assert _point_block(sc, "outline") == _point_block(cq, "outline")


@needs_openscad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("fn,length,diameter", [
    (l_profile_y, 18.0, 6.0),
    (slanted_outline, 12.0, 8.0),
])
def test_openscad_extrusion(fn, length, diameter, sigma, tmp_path):
    ref, res = _fit(fn, sigma)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["extrude_len"] == pytest.approx(length, abs=0.05)
    assert params["hole_1_d"] == pytest.approx(diameter, abs=0.05)
    assert "polygon(outline)" in code
    assert "linear_extrude(extrude_len)" in code
    against = scad_check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)
    assert against["volume_error_pct"] < 0.5


@needs_openscad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_openscad_rectangular_cutout(sigma, tmp_path):
    ref, res = _fit(slotted_plate, sigma)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert "polygon(cutout_1)" in code
    against = scad_check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.25)


@needs_openscad
def test_openscad_edit_hole_and_length(tmp_path):
    ref, res = _fit(l_profile_y, 0.0)
    code, _ = to_openscad(res.mesh, res.regions, res.tol)
    before = render_mesh(code, tmp_path / "l0")
    wider = re.sub(r"^(hole_1_d\s*=\s*)[\d.]+", r"\g<1>8.0", code, flags=re.M)
    after = render_mesh(wider, tmp_path / "l1")
    expect_hole = math.pi * (4.0 ** 2 - 3.0 ** 2) * 36.0
    assert before.volume - after.volume == pytest.approx(expect_hole, rel=0.02)
    longer = re.sub(r"^(extrude_len\s*=\s*)[\d.]+", r"\g<1>22.0", code, flags=re.M)
    grown = render_mesh(longer, tmp_path / "l2")
    assert grown.volume - before.volume == pytest.approx(L_AREA * 4.0, rel=0.02)

    ref, res = _fit(slanted_outline, 0.0)
    code, _ = to_openscad(res.mesh, res.regions, res.tol)
    before = render_mesh(code, tmp_path / "s0")
    wider = re.sub(r"^(hole_1_d\s*=\s*)[\d.]+", r"\g<1>10.0", code, flags=re.M)
    after = render_mesh(wider, tmp_path / "s1")
    expect_hole = math.pi * (5.0 ** 2 - 4.0 ** 2) * 12.0
    assert before.volume - after.volume == pytest.approx(expect_hole, rel=0.02)
    longer = re.sub(r"^(extrude_len\s*=\s*)[\d.]+", r"\g<1>16.0", code, flags=re.M)
    grown = render_mesh(longer, tmp_path / "s2")
    expect_len = (SLANT_AREA - math.pi * 4.0 ** 2) * 4.0
    assert grown.volume - before.volume == pytest.approx(expect_len, rel=0.02)
