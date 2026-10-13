"""Self-test: cadquery part -> STL -> noise -> mesh2cad -> compare to truth."""
import pytest

from mesh2cad.convert import convert_mesh
from mesh2cad.io import add_noise
from parts import PLANAR, to_mesh

EXPECTED_FACES = {"cube": 6, "l_bracket": 8, "stepped_block": 11}


@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.2])
@pytest.mark.parametrize("name", list(PLANAR))
def test_planar_roundtrip(name, sigma):
    ref = to_mesh(PLANAR[name]())
    tol = max(sigma, 0.01)
    r = convert_mesh(add_noise(ref, sigma, seed=1), tol=tol, ref=ref).report
    assert r["method"] == "planar-exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": EXPECTED_FACES[name]}
    assert r["deviation"]["hausdorff"] < max(3 * tol, 0.05)
    assert r["deviation"]["rms"] < max(tol, 0.01)
    assert r["volume_error_pct"] < 0.5
