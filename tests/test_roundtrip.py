"""Self-test: cadquery part -> STL -> noise -> mesh2cad -> compare to truth."""
import pytest

from mesh2cad.convert import convert_mesh
from mesh2cad.io import add_noise
from parts import PLANAR, to_mesh

EXPECTED_FACES = {"cube": 6, "l_bracket": 8, "stepped_block": 11}


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.2])
@pytest.mark.parametrize("name", list(PLANAR))
def test_planar_roundtrip(name, sigma, seed):
    ref = to_mesh(PLANAR[name]())
    tol = max(sigma, 0.01)
    r = convert_mesh(add_noise(ref, sigma, seed=seed), tol=tol, ref=ref).report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": EXPECTED_FACES[name]}
    assert r["deviation"]["hausdorff"] < max(3 * tol, 0.05)
    assert r["deviation"]["rms"] < max(tol, 0.01)
    assert r["volume_error_pct"] < 0.5


from parts import CURVED

CURVED_FACES = {"plate_with_hole": {"plane": 6, "cylinder": 1}, "boss": {"plane": 7, "cylinder": 1},
                "ribbed_plate": {"plane": 10, "cylinder": 2}}


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.2])
@pytest.mark.parametrize("name", list(CURVED))
def test_cylinder_roundtrip(name, sigma, seed):
    ref = to_mesh(CURVED[name]())
    tol = max(sigma, 0.01)
    r = convert_mesh(add_noise(ref, sigma, seed=seed), tol=tol, ref=ref).report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == CURVED_FACES[name]
    assert r["deviation"]["hausdorff"] < max(3 * tol, 0.1)
    assert r["volume_error_pct"] < 0.5


from parts import SPHERES, spherical_zone

SPHERE_FACES = {"spherical_dome": {"plane": 6, "sphere": 1},
                "spherical_pocket": {"plane": 6, "sphere": 1}}


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.2])
@pytest.mark.parametrize("name", list(SPHERES))
def test_sphere_roundtrip(name, sigma, seed):
    ref = to_mesh(SPHERES[name]())
    tol = max(sigma, 0.01)
    res = convert_mesh(add_noise(ref, sigma, seed=seed), tol=tol, ref=ref)
    r = res.report
    assert r["output"]["valid_solid"]
    assert r["regions"]["by_type"].get("sphere", 0) >= 1
    assert any(getattr(reg.surface, "kind", None) == "sphere" for reg in res.regions)
    assert r["method"] == "exact"
    assert r["output"]["faces"] == SPHERE_FACES[name]
    # zero-noise Hausdorff floor covers chordal error of the sphere tessellation
    assert r["deviation"]["hausdorff"] < max(3 * tol, 0.15)
    assert r["volume_error_pct"] < 0.5


@pytest.mark.parametrize("sigma", [0.0, 0.2])
def test_spherical_zone(sigma):
    """Sphere clipped by two parallel planes (two circular boundaries)."""
    ref = to_mesh(spherical_zone())
    tol = max(sigma, 0.01)
    r = convert_mesh(add_noise(ref, sigma, seed=0), tol=tol, ref=ref).report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": 2, "sphere": 1}
    assert r["deviation"]["hausdorff"] < max(3 * tol, 0.15)
    # ref mesh is chordal, so its volume sits ~0.6% under the exact zone
    assert r["volume_error_pct"] < 1.0
