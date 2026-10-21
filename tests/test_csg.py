"""Parametric export of a small boolean of straight extrusions."""
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
from mesh2cad.openscad_code import to_openscad  # noqa: E402
from parts import counterbored_plate, cut_block, orthogonal_intersection, to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")

# Z extrusion of the pentagon, clipped in Y by the X extrusion above z=6.
Z_OUTLINE = [(0.0, 0.0), (40.0, 0.0), (40.0, 16.0), (28.0, 28.0), (0.0, 28.0)]
# Projection of that intersection onto YZ. The upper step is y<=18, z<=12.
X_OUTLINE = [(0.0, 0.0), (28.0, 0.0), (28.0, 6.0), (18.0, 6.0), (18.0, 12.0), (0.0, 12.0)]
# Pentagon clipped to y<=18 (the X profile above the step). Area of z=8..12.
UPPER_SLICE = 718.0
BLOCK_OUTLINE = [(0.0, 0.0), (56.0, 0.0), (56.0, 36.0), (0.0, 36.0)]
SLOT_PROFILE = [(10.0, 8.0), (24.0, 8.0), (24.0, 18.0), (10.0, 18.0)]


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
    assert len(got) == len(expect), (got, expect)
    for px, py in expect:
        assert any(abs(px - gx) <= tol and abs(py - gy) <= tol for gx, gy in got), (px, py, got)


def _fit(fn, sigma):
    ref = to_mesh(fn())
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=0) if sigma else ref
    res = convert_mesh(mesh, tol=tol, with_report=False)
    return ref, res


@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_cadquery_intersection(sigma):
    ref, res = _fit(orthogonal_intersection, sigma)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    tol = 0.05 if sigma == 0.0 else 0.2
    assert params["extrude_z_len"] == pytest.approx(12.0, abs=tol)
    assert params["extrude_x_len"] == pytest.approx(40.0, abs=tol)
    assert params["hole_1_d"] == pytest.approx(6.0, abs=tol)
    assert "base_x" not in code
    assert "intersection" in code or ".intersect(" in code
    pt = 0.05 if sigma == 0.0 else 0.35
    _same_loop(_point_block(code, "outline_z"), Z_OUTLINE, pt)
    _same_loop(_point_block(code, "outline_x"), X_OUTLINE, pt)
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.3)
    assert against["volume_error_pct"] < 0.5


@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_cadquery_cut_extrusion(sigma):
    ref, res = _fit(cut_block, sigma)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    tol = 0.05 if sigma == 0.0 else 0.2
    assert params["extrude_len"] == pytest.approx(18.0, abs=tol)
    assert params["cut_1_len"] == pytest.approx(28.0, abs=tol)
    assert params["cut_1_0"] == pytest.approx(12.0, abs=tol)
    assert params["hole_1_d"] == pytest.approx(6.0, abs=tol)
    assert "length along Z" in code
    assert "length along X" in code
    pt = 0.05 if sigma == 0.0 else 0.35
    _same_loop(_point_block(code, "outline"), BLOCK_OUTLINE, pt)
    _same_loop(_point_block(code, "cut_1"), SLOT_PROFILE, pt)
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.3)


@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_cadquery_counterbore(sigma):
    ref, res = _fit(counterbored_plate, sigma)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    tol = 0.05 if sigma == 0.0 else 0.2
    assert params["base_x"] == pytest.approx(40.0, abs=tol)
    assert params["base_y"] == pytest.approx(28.0, abs=tol)
    assert params["base_t"] == pytest.approx(14.0, abs=tol)
    assert params["base_z0"] == pytest.approx(-7.0, abs=tol)
    assert params["hole_1_d"] == pytest.approx(6.0, abs=tol)
    assert params["hole_1_cbore_d"] == pytest.approx(12.0, abs=tol)
    assert params["hole_1_cbore_depth"] == pytest.approx(4.0, abs=tol)
    assert "extrude_len" not in code
    assert "counterbore depth from the top face" in code
    against = cq_check(code, ref)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.3)


def test_edit_hole_diameter_and_extrusion_length():
    """Hole diameter and extrusion length, each checked against the analytic delta."""
    _ref, res = _fit(orthogonal_intersection, 0.0)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    v0 = abs(run_script(code).volume)
    wider = re.sub(r"^hole_1_d\s*=\s*[\d.]+", "hole_1_d = 8.0", code, flags=re.M)
    # Ø6 → Ø8 along X, through the 40 mm intersection. The disk stays inside.
    assert v0 - abs(run_script(wider).volume) == pytest.approx(math.pi * 7.0 * 40.0, rel=0.01)
    shorter = re.sub(r"^extrude_z_len\s*=\s*[\d.]+", "extrude_z_len = 8.0", code, flags=re.M)
    # z=8..12 is the pentagon clipped to y<=18. The intersections overlap by
    # 0.02 mm so shared faces stay manifold; that trims this delta by that much.
    assert v0 - abs(run_script(shorter).volume) == pytest.approx(UPPER_SLICE * 4.0, rel=0.01)

    _ref, res = _fit(cut_block, 0.0)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    v0 = abs(run_script(code).volume)
    wider = re.sub(r"^hole_1_d\s*=\s*[\d.]+", "hole_1_d = 8.0", code, flags=re.M)
    assert v0 - abs(run_script(wider).volume) == pytest.approx(math.pi * 7.0 * 36.0, rel=0.01)
    taller = re.sub(r"^extrude_len\s*=\s*[\d.]+", "extrude_len = 22.0", code, flags=re.M)
    # the hole (z=4) and the slot (opens on z=18) do not meet the added slab
    assert abs(run_script(taller).volume) - v0 == pytest.approx(56.0 * 36.0 * 4.0, rel=0.01)

    _ref, res = _fit(counterbored_plate, 0.0)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    v0 = abs(run_script(code).volume)
    wider = re.sub(r"^hole_1_d\s*=\s*[\d.]+", "hole_1_d = 8.0", code, flags=re.M)
    # the counterbore already removed the annulus, so only the remaining depth changes
    assert v0 - abs(run_script(wider).volume) == pytest.approx(math.pi * 7.0 * 10.0, rel=0.01)


def test_exporters_share_csg_parameters():
    for fn in (orthogonal_intersection, cut_block, counterbored_plate):
        _ref, res = _fit(fn, 0.0)
        cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
        sc, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
        assert _params(sc) == _params(cq)


@needs_openscad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("fn", [orthogonal_intersection, cut_block, counterbored_plate])
def test_openscad_csg(fn, sigma, tmp_path):
    ref, res = _fit(fn, sigma)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    assert "linear_extrude" in code or "cube" in code
    against = scad_check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.05 if sigma == 0.0 else 0.3)
    assert against["volume_error_pct"] < 0.5
