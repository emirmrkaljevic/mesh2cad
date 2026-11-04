"""Constant-radius fillets: plane-plane rounds come back as partial cylinders."""
import numpy as np
import cadquery as cq
import pytest

from mesh2cad.convert import convert_mesh
from mesh2cad.io import add_noise
from parts import FILLET_RADIUS, FILLETS, ROUNDED, ROUNDED_RADIUS, to_mesh

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


def _unit_plane(plane):
    n = np.asarray(plane.normal, float)
    ln = float(np.linalg.norm(n))
    return n / ln, float(plane.d) / ln


def _tangent(regions, tol=1e-3):
    """Every cylinder meets two planes, every sphere three, and the fillet
    axis runs through the two corner centres."""
    planes = [r.surface for r in regions if r.surface.kind == "plane"]
    cyls = [r.surface for r in regions if r.surface.kind == "cylinder"]
    sphs = [r.surface for r in regions if r.surface.kind == "sphere"]
    for cyl in cyls:
        axis = np.asarray(cyl.axis, float)
        axis = axis / np.linalg.norm(axis)
        hits = 0
        for plane in planes:
            n, d = _unit_plane(plane)
            if abs(float(n @ axis)) > 0.05:
                continue
            if abs(abs(float(n @ cyl.point - d)) - cyl.radius) <= tol:
                hits += 1
        assert hits == 2
        near = 0
        for sph in sphs:
            rel = np.asarray(sph.center, float) - np.asarray(cyl.point, float)
            off = float(np.linalg.norm(rel - axis * float(rel @ axis)))
            if off <= tol and abs(sph.radius - cyl.radius) <= tol:
                near += 1
        assert near == 2
    for sph in sphs:
        hits = 0
        for plane in planes:
            n, d = _unit_plane(plane)
            if abs(abs(float(n @ sph.center - d)) - sph.radius) <= tol:
                hits += 1
        assert hits == 3


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.2])
def test_rounded_box_roundtrip(sigma, seed):
    """Every edge filleted: 6 planes, 12 cylinders, 8 sphere octants, all tangent."""
    ref = to_mesh(ROUNDED["rounded_box"]())
    tol = max(sigma, 0.01)
    res = convert_mesh(add_noise(ref, sigma, seed=seed), tol=tol, ref=ref)
    r = res.report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": 6, "cylinder": 12, "sphere": 8}
    truth = ROUNDED_RADIUS["rounded_box"]
    radii = [reg.surface.radius for reg in res.regions
             if reg.surface.kind in ("cylinder", "sphere")]
    assert len(radii) == 20
    assert max(abs(rad - truth) for rad in radii) < (0.05 if sigma < 0.2 else 0.2)
    assert r["deviation"]["hausdorff"] < (0.08 if sigma < 0.2 else 0.25)
    assert r["volume_error_pct"] < 0.5
    _tangent(res.regions)


def test_coarse_rounded_box_is_exact():
    """A coarser tessellation of the same corner is still one analytic solid."""
    part = cq.Workplane().box(16, 14, 10).edges().fillet(1.2)
    m = to_mesh(part, tol=0.3, ang=0.6, max_edge=3.5)
    res = convert_mesh(m, tol=0.05)
    r = res.report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"] == {"plane": 6, "cylinder": 12, "sphere": 8}
    _tangent(res.regions, tol=1e-3)


def test_coarse_boss_base_fillet_is_a_torus():
    """A coarser tessellation of a boss-base fillet is still one torus."""
    part = (cq.Workplane().box(30, 30, 6).faces(">Z").workplane().circle(6).extrude(8)
            .edges(cq.NearestToPointSelector((0, 0, 3))).fillet(1.0))
    m = to_mesh(part, tol=0.2, ang=0.5, max_edge=3.0)
    r = convert_mesh(m, tol=0.05).report
    assert r["method"] == "exact"
    assert r["output"]["valid_solid"]
    assert r["output"]["faces"].get("torus") == 1
    assert r["tori"][0]["minor"] == pytest.approx(1.0, abs=0.08)
    assert r["tori"][0]["major"] == pytest.approx(7.0, abs=0.08)
