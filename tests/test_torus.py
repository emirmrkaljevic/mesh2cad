"""Fillets between a plane and a coaxial cylinder come back as torus faces."""
import pytest

from mesh2cad.cadquery_code import to_cadquery
from mesh2cad.convert import convert_mesh
from mesh2cad.io import add_noise
from mesh2cad.parametric import Unsupported
from parts import TORI, TORUS_MAJOR, TORUS_MINOR, filleted_blind_hole, rounded_hole_both, to_mesh


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


def test_torus_is_not_a_parametric_script():
    ref = to_mesh(TORI["filleted_boss"]())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    with pytest.raises(Unsupported) as ei:
        to_cadquery(res.mesh, res.regions, res.tol)
    assert "torus" in str(ei.value)
