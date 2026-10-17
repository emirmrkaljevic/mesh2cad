"""Constant-radius fillets: plane-plane rounds come back as partial cylinders."""
import cadquery as cq
import pytest

from mesh2cad.convert import convert_mesh
from mesh2cad.io import add_noise
from parts import FILLET_RADIUS, FILLETS, to_mesh

FACES = {"filleted_plate": {"plane": 6, "cylinder": 4},
         "filleted_box": {"plane": 6, "cylinder": 4},
         "filleted_bracket": {"plane": 8, "cylinder": 6}}


def _radii(res):
    return [r.surface.radius for r in res.regions
            if getattr(r.surface, "kind", None) == "cylinder"]


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("name", list(FILLETS))
def test_fillet_roundtrip(name, sigma, seed):
    ref = to_mesh(FILLETS[name]())
    tol = max(sigma, 0.01)
    res = convert_mesh(add_noise(ref, sigma, seed=seed), tol=tol, ref=ref)
    r = res.report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == FACES[name]
    truth = FILLET_RADIUS[name]
    radii = _radii(res)
    assert len(radii) == FACES[name]["cylinder"]
    assert max(abs(rad - truth) for rad in radii) < (0.08 if sigma == 0.0 else 0.35)
    assert r["deviation"]["hausdorff"] < (0.08 if sigma == 0.0 else 0.2)
    assert r["volume_error_pct"] < 0.5


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", list(FILLETS))
def test_fillet_high_noise_stays_valid(name, seed):
    """0.2 mm is past the noise these rounds survive reliably.

    Some seeds still sew an exact solid, and some of those miss a corner
    (Hausdorff well above the noise). The rest fall back. Either way the
    solid has to be valid and close to the mesh.
    """
    ref = to_mesh(FILLETS[name]())
    r = convert_mesh(add_noise(ref, 0.2, seed=seed), tol=0.2, ref=ref).report
    assert r["output"]["valid_solid"]
    assert r["method"] in ("exact", "faceted-hybrid", "faceted")
    assert r["deviation"]["hausdorff"] < 1.5
    assert r["volume_error_pct"] < 1.5


def test_spherical_corner_falls_back_valid():
    """Filleting every edge of a box adds sphere octants. Those stay faceted."""
    part = cq.Workplane().box(16, 14, 10).edges().fillet(1.2)
    m = to_mesh(part, tol=0.3, ang=0.6, max_edge=3.5)
    r = convert_mesh(m, tol=0.05).report
    assert r["output"]["valid_solid"]
    assert r["method"] in ("faceted-hybrid", "faceted")


def test_boss_base_fillet_falls_back_valid():
    """A plane-cylinder fillet is a torus and is not rebuilt."""
    part = (cq.Workplane().box(30, 30, 6).faces(">Z").workplane().circle(6).extrude(8)
            .edges(cq.NearestToPointSelector((0, 0, 3))).fillet(1.0))
    m = to_mesh(part, tol=0.2, ang=0.5, max_edge=3.0)
    r = convert_mesh(m, tol=0.05).report
    assert r["output"]["valid_solid"]
    assert r["method"] in ("faceted-hybrid", "faceted")
