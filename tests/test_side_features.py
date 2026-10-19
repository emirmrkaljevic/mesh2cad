"""Holes and bosses whose axis is X or Y, in both parametric exporters."""
import math
import re
import shutil

import pytest

pytest.importorskip("cadquery")
from mesh2cad.cadquery_code import check as cq_check  # noqa: E402
from mesh2cad.cadquery_code import run_script, to_cadquery  # noqa: E402
from mesh2cad.io import add_noise, clean  # noqa: E402
from mesh2cad.openscad_code import check as scad_check  # noqa: E402
from mesh2cad.openscad_code import render_mesh, to_openscad  # noqa: E402
from mesh2cad.parametric import Unsupported  # noqa: E402
from mesh2cad.segment import segment  # noqa: E402
from parts import angled_hole, side_hole_boss, to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _ids(params, kind):
    found = []
    for key in params:
        head, _, rest = key.partition("_")
        number, sep, _tail = rest.partition("_")
        if head == kind and sep and number.isdigit():
            i = int(number)
            if i not in found:
                found.append(i)
    return sorted(found)


class _Fit:
    def __init__(self, mesh, regions, tol):
        self.mesh, self.regions, self.tol = mesh, regions, tol


def _convert(sigma=0.0, seed=0, fn=side_hole_boss):
    """Segment the way ``--cadquery`` / ``--openscad`` do.

    Those commands fit the mesh and emit a script; they do not build the STEP
    solid first. The blind end of the side hole is a small disk that the
    segmenter keeps on the cylinder, so the exact B-rep builder declines the
    part and the STEP path falls back to a faceted solid.
    """
    ref = to_mesh(fn())
    tol = max(sigma, 0.01)
    mesh = clean(add_noise(ref, sigma, seed=seed))
    return ref, _Fit(mesh, segment(mesh, tol), tol)


def _expect_side_features(params, tol=0.05):
    """The Ø8 X through-hole, the Ø6 Y blind hole, and the Ø10 +X boss."""
    assert params["base_x"] == pytest.approx(60.0, abs=tol)
    assert params["base_y"] == pytest.approx(40.0, abs=tol)
    assert params["base_t"] == pytest.approx(36.0, abs=tol)
    holes = _ids(params, "hole")
    bosses = _ids(params, "boss")
    assert len(holes) == 2
    assert len(bosses) == 1
    through = [i for i in holes if f"hole_{i}_depth" not in params][0]
    blind = [i for i in holes if f"hole_{i}_depth" in params][0]
    assert f"hole_{through}_x" not in params
    assert params[f"hole_{through}_d"] == pytest.approx(8.0, abs=tol)
    assert params[f"hole_{through}_y"] == pytest.approx(0.0, abs=tol)
    assert params[f"hole_{through}_z"] == pytest.approx(-6.0, abs=tol)
    assert f"hole_{blind}_y" not in params
    assert params[f"hole_{blind}_d"] == pytest.approx(6.0, abs=tol)
    assert params[f"hole_{blind}_x"] == pytest.approx(-12.0, abs=tol)
    assert params[f"hole_{blind}_z"] == pytest.approx(8.0, abs=tol)
    assert params[f"hole_{blind}_depth"] == pytest.approx(12.0, abs=tol)
    boss = bosses[0]
    assert f"boss_{boss}_x" not in params
    assert params[f"boss_{boss}_d"] == pytest.approx(10.0, abs=tol)
    assert params[f"boss_{boss}_h"] == pytest.approx(8.0, abs=tol)
    assert params[f"boss_{boss}_y"] == pytest.approx(10.0, abs=tol)
    assert params[f"boss_{boss}_z"] == pytest.approx(8.0, abs=tol)
    return through, blind, boss


@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_cadquery_side_features(sigma):
    ref, res = _convert(sigma)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    through, blind, boss = _expect_side_features(_params(code))
    assert f"cq.Vector({cad_dir('X')})" in code or "cq.Vector(1, 0, 0)" in code
    assert "cq.Vector(0, 1, 0)" in code
    assert f"hole_{through}_d / 2" in code
    assert f"hole_{blind}_depth" in code
    assert f"boss_{boss}_h" in code
    against = cq_check(code, ref)
    assert against["max_deviation"] < 0.08
    assert against["volume_error_pct"] < 0.5


def cad_dir(axis):
    return {"X": "1, 0, 0", "Y": "0, 1, 0"}[axis]


def test_cadquery_edit_side_hole_diameter():
    ref, res = _convert(0.0)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    params = _params(code)
    through, _, _ = _expect_side_features(params)
    name = f"hole_{through}_d"
    d0 = params[name]
    v0 = abs(run_script(code).volume)
    bigger = re.sub(rf"^({name}\s*=\s*)[-\d.]+", rf"\g<1>{d0 + 2.0:.1f}", code, flags=re.M)
    v1 = abs(run_script(bigger).volume)
    expect = math.pi * (((d0 + 2.0) / 2.0) ** 2 - (d0 / 2.0) ** 2) * params["base_x"]
    assert v0 - v1 == pytest.approx(expect, rel=0.01)


def test_parameters_match_between_exporters():
    _ref, res = _convert(0.0)
    cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    sc, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    assert _params(sc) == _params(cq)


def test_angled_axis_is_refused():
    _ref, res = _convert(fn=angled_hole)
    with pytest.raises(Unsupported) as ei:
        to_cadquery(res.mesh, res.regions, res.tol)
    assert "aligned" in str(ei.value).lower() or "not vertical" in str(ei.value).lower()
    with pytest.raises(Unsupported):
        to_openscad(res.mesh, res.regions, res.tol, verify=False)


@needs_openscad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_openscad_side_features(sigma, tmp_path):
    ref, res = _convert(sigma)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    through, blind, boss = _expect_side_features(_params(code))
    assert "rotate([0, 90, 0])" in code
    assert "rotate([-90, 0, 0])" in code
    assert f"cylinder(h=base_x + 2, d=hole_{through}_d);" in code
    assert f"cylinder(h=hole_{blind}_depth + 1, d=hole_{blind}_d);" in code
    assert f"cylinder(h=boss_{boss}_h, d=boss_{boss}_d);" in code
    against = scad_check(code, ref, tmp_path)
    assert against["max_deviation"] < 0.08
    assert against["volume_error_pct"] < 0.5


def _slider_hi(code, name):
    m = re.search(rf"^{name}\s*=\s*[-\d.]+;\s*//\s*\[([^\]]+)\]", code, re.M)
    assert m, name
    return float(m.group(1).split(":")[-1])


def test_openscad_side_slider_uses_the_face():
    _ref, res = _convert(0.0)
    code, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    params = _params(code)
    through, blind, boss = _expect_side_features(params)
    # the X hole is limited by the YZ face, not by the length of the hole
    assert _slider_hi(code, f"hole_{through}_d") <= min(params["base_y"], params["base_t"]) + 1e-6
    assert _slider_hi(code, f"hole_{blind}_depth") == pytest.approx(params["base_y"], abs=0.2)
    assert _slider_hi(code, f"boss_{boss}_d") <= min(params["base_y"], params["base_t"]) + 1e-6
    assert re.search(rf"^hole_{through}_y\s*=\s*[-\d.]+;\s*$", code, re.M)


@needs_openscad
def test_openscad_edit_side_hole_diameter(tmp_path):
    ref, res = _convert(0.0)
    code, _ = to_openscad(res.mesh, res.regions, res.tol)
    params = _params(code)
    through, _, _ = _expect_side_features(params)
    name = f"hole_{through}_d"
    d0 = params[name]
    before = render_mesh(code, tmp_path / "before")
    bigger = re.sub(rf"^({name}\s*=\s*)[-\d.]+", rf"\g<1>{d0 + 2.0:.1f}", code, flags=re.M)
    after = render_mesh(bigger, tmp_path / "after")
    expect = math.pi * (((d0 + 2.0) / 2.0) ** 2 - (d0 / 2.0) ** 2) * params["base_x"]
    assert abs(before.volume) - abs(after.volume) == pytest.approx(expect, rel=0.02)
