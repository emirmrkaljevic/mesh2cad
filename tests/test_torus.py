"""Fillets between a plane and a coaxial cylinder come back as torus faces."""
import re

import pytest

from mesh2cad.cadquery_code import check as cq_check
from mesh2cad.cadquery_code import run_script, to_cadquery
from mesh2cad.convert import convert_mesh
from mesh2cad.freecad_code import freecad_bin, to_freecad
from mesh2cad.io import add_noise
from mesh2cad.openscad_code import openscad_bin, to_openscad
from mesh2cad.parametric import Unsupported
from mesh2cad.params import edit_script, list_parameters
from parts import TORI, TORUS_MAJOR, TORUS_MINOR, filleted_blind_hole, rounded_hole_both, to_mesh

needs_openscad = pytest.mark.skipif(openscad_bin() is None, reason="openscad is not installed")
needs_freecad = pytest.mark.skipif(freecad_bin() is None, reason="freecadcmd is not installed")

# (part, parameter, true radius). The fitted radius is the torus minor, which
# on a clean mesh sits a few hundredths under the nominal value.
TORUS_SCRIPTS = (
    ("filleted_boss", TORI["filleted_boss"], "boss_1_fillet", 2.0),
    ("rounded_hole", TORI["rounded_hole"], "hole_1_rim_r", 1.5),
    ("rounded_hole_both", rounded_hole_both, "hole_1_rim_r", 1.2),
    ("filleted_blind_hole", filleted_blind_hole, "hole_1_floor_r", 1.5),
)


def _minor(res):
    return [t["minor"] for t in res.report.get("tori", [])]


def _major(res):
    return [t["major"] for t in res.report.get("tori", [])]


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_boss_fillet_roundtrip(sigma, seed):
    ref = to_mesh(TORI["filleted_boss"]())
    tol = max(sigma, 0.01)
    res = convert_mesh(add_noise(ref, sigma, seed=seed), tol=tol, ref=ref)
    r = res.report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": 7, "cylinder": 1, "torus": 1}
    assert max(abs(v - TORUS_MINOR["filleted_boss"]) for v in _minor(res)) < (0.05 if sigma == 0 else 0.40)
    assert max(abs(v - TORUS_MAJOR["filleted_boss"]) for v in _major(res)) < (0.05 if sigma == 0 else 0.50)
    assert r["deviation"]["hausdorff"] < (0.08 if sigma == 0 else 0.25)
    assert r["volume_error_pct"] < 0.5


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_rounded_hole_clean(seed):
    ref = to_mesh(TORI["rounded_hole"]())
    res = convert_mesh(add_noise(ref, 0.0, seed=seed), tol=0.01, ref=ref)
    r = res.report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": 6, "cylinder": 1, "torus": 1}
    assert max(abs(v - 1.5) for v in _minor(res)) < 0.12
    assert max(abs(v - 7.5) for v in _major(res)) < 0.12
    assert r["deviation"]["hausdorff"] < 0.08
    assert r["volume_error_pct"] < 0.5


def test_rounded_hole_moderate_noise():
    """0.05 mm keeps a torus on at least two seeds. A miss stays a valid solid."""
    ref = to_mesh(TORI["rounded_hole"]())
    exact = 0
    for seed in (0, 1, 2):
        res = convert_mesh(add_noise(ref, 0.05, seed=seed), tol=0.05, ref=ref)
        r = res.report
        assert r["output"]["valid_solid"]
        assert r["volume_error_pct"] < 1.5
        if r["method"] == "exact" and r["output"]["faces"].get("torus") == 1:
            exact += 1
            assert max(abs(v - 1.5) for v in _minor(res)) < 0.40
            assert r["deviation"]["hausdorff"] < 0.20
        else:
            assert r["method"] in ("faceted-hybrid", "faceted")
            assert r["deviation"]["hausdorff"] < 1.0
    assert exact >= 2


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", ["filleted_boss", "rounded_hole"])
def test_torus_high_noise_stays_valid(name, seed):
    """0.2 mm loses the round. The solid is still valid, with or without a torus."""
    ref = to_mesh(TORI[name]())
    r = convert_mesh(add_noise(ref, 0.2, seed=seed), tol=0.2, ref=ref).report
    assert r["output"]["valid_solid"]
    assert r["method"] in ("exact", "faceted-hybrid", "faceted")
    assert r["deviation"]["hausdorff"] < 1.2
    assert r["volume_error_pct"] < 1.5


def test_blind_hole_floor_fillet():
    """Concave fillet where a blind-hole floor meets the wall (major < hole radius)."""
    ref = to_mesh(filleted_blind_hole())
    res = convert_mesh(ref, tol=0.01, ref=ref)
    r = res.report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": 7, "cylinder": 1, "torus": 1}
    assert _minor(res)[0] == pytest.approx(1.5, abs=0.1)
    assert _major(res)[0] == pytest.approx(3.5, abs=0.1)
    assert r["deviation"]["hausdorff"] < 0.08


def test_both_hole_rims_clean():
    ref = to_mesh(rounded_hole_both())
    res = convert_mesh(ref, tol=0.01, ref=ref)
    r = res.report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": 6, "cylinder": 1, "torus": 2}
    assert max(abs(v - 1.2) for v in _minor(res)) < 0.12
    assert r["deviation"]["hausdorff"] < 0.08


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", ["filleted_blind_hole", "rounded_hole_both"])
def test_other_tori_moderate_noise_stays_valid(name, seed):
    fn = filleted_blind_hole if name == "filleted_blind_hole" else rounded_hole_both
    ref = to_mesh(fn())
    r = convert_mesh(add_noise(ref, 0.05, seed=seed), tol=0.05, ref=ref).report
    assert r["output"]["valid_solid"]
    assert r["method"] in ("exact", "faceted-hybrid", "faceted")
    assert r["deviation"]["hausdorff"] < 2.0
    assert r["volume_error_pct"] < 1.5


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _convert(fn, sigma=0.0, seed=0):
    ref = to_mesh(fn())
    tol = max(sigma, 0.01)
    res = convert_mesh(add_noise(ref, sigma, seed=seed), tol=tol, with_report=False)
    return ref, res


@pytest.mark.parametrize("name,fn,key,radius", TORUS_SCRIPTS)
def test_cadquery_recovers_torus_fillet(name, fn, key, radius):
    """Clean mesh: the script names the fillet and stays inside the usual gate."""
    ref, res = _convert(fn)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params[key] == pytest.approx(radius, abs=0.12)
    if name == "rounded_hole_both":
        assert params["hole_1_bot_rim_r"] == pytest.approx(radius, abs=0.12)
    if name == "filleted_boss":
        assert "boss_1_fillet" in code and "_circle_edges" in code
    against = cq_check(code, ref)
    assert against["max_deviation"] < 0.08
    assert against["volume_error_pct"] < 0.5


@needs_openscad
@pytest.mark.parametrize("name,fn,key,radius", TORUS_SCRIPTS)
def test_openscad_recovers_torus_fillet(name, fn, key, radius, tmp_path):
    ref, res = _convert(fn)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params[key] == pytest.approx(radius, abs=0.12)
    assert "rotate_extrude()" in code
    if name == "rounded_hole_both":
        assert params["hole_1_bot_rim_r"] == pytest.approx(radius, abs=0.12)
    from mesh2cad.openscad_code import check
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < 0.08
    assert against["volume_error_pct"] < 0.5


@needs_freecad
@pytest.mark.parametrize("name,fn,key,radius", TORUS_SCRIPTS)
def test_freecad_recovers_torus_fillet(name, fn, key, radius, tmp_path):
    ref, res = _convert(fn)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params[key] == pytest.approx(radius, abs=0.12)
    assert "makeFillet" in code and "_circle_edges" in code
    from mesh2cad.freecad_code import check
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < 0.08
    assert against["volume_error_pct"] < 0.5


def test_freecad_script_names_the_fillet_without_freecad():
    """The macro text is checked even when freecadcmd is absent."""
    _ref, res = _convert(TORI["filleted_boss"])
    code, _q = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert _params(code)["boss_1_fillet"] == pytest.approx(2.0, abs=0.05)
    assert "result.makeFillet(boss_1_fillet," in code


def test_edit_boss_fillet_changes_the_solid(tmp_path):
    ref = to_mesh(TORI["filleted_boss"]())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    script = tmp_path / "boss.py"
    script.write_text(code, encoding="utf-8")
    listed = {p["name"]: p for p in list_parameters(str(script))["parameters"]}
    assert listed["boss_1_fillet"]["kind"] == "radius"
    assert listed["boss_1_fillet"]["unit"] == "mm"
    before = run_script(code)
    result = edit_script(str(script), {"boss_1_fillet": 3}, output_path=str(tmp_path / "boss.step"),
                         output_script=str(tmp_path / "edited.py"))
    after = run_script((tmp_path / "edited.py").read_text(encoding="utf-8"))
    # A point in the extra ring of a 3 mm fillet, outside the 2 mm fillet.
    assert not bool(before.contains([[9.2, 0.0, 3.5]])[0])
    assert bool(after.contains([[9.2, 0.0, 3.5]])[0])
    assert after.volume > before.volume
    assert result["valid_solid"] is True
    assert result["changes"][0]["old"] == pytest.approx(2.0, abs=0.05)
    assert result["changes"][0]["new"] == pytest.approx(3.0)


def test_edit_hole_rim_removes_material(tmp_path):
    ref = to_mesh(TORI["rounded_hole"]())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    script = tmp_path / "hole.py"
    script.write_text(code, encoding="utf-8")
    before = run_script(code)
    edit_script(str(script), {"hole_1_rim_r": 2.5}, output_script=str(tmp_path / "edited.py"))
    after = run_script((tmp_path / "edited.py").read_text(encoding="utf-8"))
    # Near the top face, outside the hole: material at the fitted rim, void at 2.5 mm.
    assert bool(before.contains([[7.0, 0.0, 3.7]])[0])
    assert not bool(after.contains([[7.0, 0.0, 3.7]])[0])
    assert after.volume < before.volume


def test_moderate_noise_can_keep_the_boss_fillet():
    """0.05 mm, seed 0: the fillet is a parameter and the usual gate still holds."""
    _ref, res = _convert(TORI["filleted_boss"], sigma=0.05, seed=0)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    assert _params(code)["boss_1_fillet"] == pytest.approx(1.7, abs=0.05)


def test_high_noise_boss_is_not_given_a_fillet():
    """0.2 mm, seed 0: the torus is gone and the sharp solid misses the gate."""
    _ref, res = _convert(TORI["filleted_boss"], sigma=0.2, seed=0)
    assert not any(getattr(r.surface, "kind", None) == "torus" for r in res.regions)
    with pytest.raises(Unsupported, match="deviates by"):
        to_cadquery(res.mesh, res.regions, res.tol)


def test_edit_floor_fillet_adds_material(tmp_path):
    ref = to_mesh(filleted_blind_hole())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    script = tmp_path / "blind.py"
    script.write_text(code, encoding="utf-8")
    before = run_script(code)
    edit_script(str(script), {"hole_1_floor_r": 2.4}, output_script=str(tmp_path / "edited.py"))
    after = run_script((tmp_path / "edited.py").read_text(encoding="utf-8"))
    # Inside the sharp corner of the floor, still void at 1.5 mm, filled at 2.4 mm.
    assert not bool(before.contains([[4.2, 0.0, -1.6]])[0])
    assert bool(after.contains([[4.2, 0.0, -1.6]])[0])
    assert after.volume > before.volume
