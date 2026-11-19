import json
import urllib.request
from pathlib import Path

import pytest

from mesh2cad import brep
from mesh2cad.cli import main
from mesh2cad.convert import convert_mesh
from mesh2cad.fit import fit_cylinder, fit_plane, fit_sphere, refine_sphere
from mesh2cad.io import add_noise
from parts import cube, to_mesh
import numpy as np


def test_fit_plane():
    rng = np.random.default_rng(0)
    p = np.c_[rng.uniform(-5, 5, (200, 2)), np.full(200, 3.0)]
    pl = fit_plane(p)
    assert abs(abs(pl.normal[2]) - 1) < 1e-9 and abs(abs(pl.d) - 3) < 1e-9


def test_fit_cylinder():
    t = np.linspace(0, 2 * np.pi, 100, endpoint=False)
    z = np.tile([0.0, 5.0], 50)
    pts = np.c_[2 + 4 * np.cos(t), -1 + 4 * np.sin(t), z]
    nrm = np.c_[np.cos(t), np.sin(t), np.zeros_like(t)]
    c = fit_cylinder(pts, nrm)
    assert abs(c.radius - 4) < 1e-6 and abs(abs(c.axis[2]) - 1) < 1e-6


def test_fit_sphere():
    rng = np.random.default_rng(0)
    phi = rng.uniform(0.2, 1.2, 300)
    th = rng.uniform(0, 2 * np.pi, 300)
    nrm = np.c_[np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th), np.cos(phi)]
    center = np.array([1.0, -2.0, 0.5])
    pts = center + 5.0 * nrm
    sph = refine_sphere(fit_sphere(pts, nrm), pts)
    assert np.linalg.norm(sph.center - center) < 1e-6
    assert abs(sph.radius - 5.0) < 1e-6
    assert sph.distance(pts).max() < 1e-6
    assert np.linalg.norm(sph.project(pts) - pts, axis=1).max() < 1e-6


def test_cli_convert(tmp_path):
    stl = tmp_path / "cube.stl"
    to_mesh(cube()).export(stl)
    out = tmp_path / "cube.step"
    rep = tmp_path / "r.json"
    assert main(["convert", str(stl), "-o", str(out), "--report", str(rep)]) == 0
    assert out.read_text().startswith("ISO-10303-21")
    r = json.loads(rep.read_text())
    assert r["output"]["faces"] == {"plane": 6} and r["output"]["valid_solid"]


def test_snap_sphere():
    from mesh2cad import brep
    from mesh2cad.fit import Sphere
    from mesh2cad.io import clean
    from mesh2cad.segment import segment
    from parts import spherical_dome
    m = clean(to_mesh(spherical_dome()))
    regs = segment(m, 0.01)
    assert any(isinstance(r.surface, Sphere) for r in regs)
    x = brep.snap_vertices(m, regs, 0.01)
    for r in regs:
        if isinstance(r.surface, Sphere):
            vs = np.unique(m.faces[r.faces])
            assert r.surface.distance(x[vs]).max() < 1e-5


def test_snap_and_fallback():
    from mesh2cad import brep
    from mesh2cad.convert import convert_mesh
    from mesh2cad.io import add_noise, clean
    from mesh2cad.segment import segment
    from parts import plate_with_hole
    m = clean(add_noise(to_mesh(cube()), 0.05))
    regs = segment(m, 0.05)
    x = brep.snap_vertices(m, regs, 0.05)
    import numpy as np
    for r in regs:
        if r.surface is None:
            continue
        vs = np.unique(m.faces[r.faces])
        assert r.surface.distance(x[vs]).max() < 1e-6
    ref = to_mesh(plate_with_hole())
    r = convert_mesh(ref, tol=0.01, cylinders=False, ref=ref).report
    assert r["method"] in ("faceted-hybrid", "faceted")
    assert r["output"]["valid_solid"]
    assert r["deviation"]["hausdorff"] < 0.5


def test_countersunk_default_tol_no_crash(tmp_path):
    """Regression (0.4.1): default --tol lost the countersink cone and the hybrid
    fallback segfaulted in UnifySameDomain. Now tighter tols are retried and
    unify only runs on valid shells."""
    from parts import countersunk_plate
    stl = tmp_path / "csk.stl"
    to_mesh(countersunk_plate()).export(stl)
    out = tmp_path / "csk.json"
    assert main(["convert", str(stl), "--report", str(out)]) == 0
    r = json.loads(out.read_text())
    assert r["method"] == "exact" and r["output"]["valid_solid"]
    assert r["output"]["faces"].get("cone") == 1
    # the loose tol on its own must still produce a valid solid via the guarded fallback
    assert main(["convert", str(stl), "--tol", "0.0504", "-o", str(tmp_path / "loose.step")]) == 0


def _step_volume(path):
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Reader
    reader = STEPControl_Reader()
    assert reader.ReadFile(str(path)) == IFSelect_RetDone
    reader.TransferRoots()
    return brep.signed_volume(reader.OneShape())


def _coarse_plate(points):
    """CAD tessellation of an 80x50x6 plate, holes of radius 3.3 mm."""
    import cadquery as cq
    import trimesh
    part = (cq.Workplane().box(80, 50, 6).faces(">Z").workplane()
            .pushPoints(points).hole(6.6))
    verts, faces = part.val().tessellate(0.1, 0.4)
    mesh = trimesh.Trimesh(np.array([p.toTuple() for p in verts]), np.array(faces), process=True)
    if len(mesh.faces) and (mesh.area_faces < 1e-12).any():
        mesh.update_faces(mesh.area_faces >= 1e-12)
        mesh.remove_unreferenced_vertices()
    return mesh


_GRID = [(-20.0, -12.0), (20.0, -12.0), (-20.0, 12.0), (20.0, 12.0)]
_ROW = [(-30.0, 0.0), (-10.0, 0.0), (10.0, 0.0), (30.0, 0.0)]


def test_clean_four_hole_plate_is_exact():
    ref = _coarse_plate(_GRID)
    res = convert_mesh(ref, ref=ref)
    report = res.report
    assert report["method"] == "exact"
    assert report["output"]["valid_solid"]
    assert report["output"]["faces"] == {"plane": 6, "cylinder": 4}
    radii = sorted(float(region.surface.radius) for region in res.regions
                   if getattr(region.surface, "kind", None) == "cylinder")
    assert radii == pytest.approx([3.3, 3.3, 3.3, 3.3], abs=1e-3)
    assert brep.signed_volume(res.shape) > 0


@pytest.mark.parametrize("points,seed", [(_GRID, 0), (_ROW, 1)])
def test_noisy_plate_is_not_written_inside_out(points, seed, tmp_path):
    """sigma 0.02 mm at tol 0.02 used to sew a hybrid about 6 mm off.

    One grid seed measured Hausdorff 6.038. One row seed sewed a negative
    volume. Both have to come back as an outward solid close to the mesh,
    and the STEP file has to agree.
    """
    ref = _coarse_plate(points)
    noisy = add_noise(ref, 0.02, seed=seed)
    res = convert_mesh(noisy, tol=0.02, ref=ref)
    report = res.report
    assert report["output"]["valid_solid"]
    assert report["method"] in ("exact", "faceted-hybrid", "faceted")
    assert report["deviation"]["hausdorff"] < 1.0
    assert report["volume_error_pct"] < 2.0
    signed = brep.signed_volume(res.shape)
    assert signed > 0.95 * abs(float(ref.volume))
    out = tmp_path / "plate.step"
    brep.write_step(res.shape, str(out))
    assert _step_volume(out) > 0.95 * abs(float(ref.volume))


def test_write_step_turns_an_inside_out_solid_outward(tmp_path):
    solid = convert_mesh(to_mesh(cube()), with_report=False).shape
    out = tmp_path / "cube.step"
    brep.write_step(solid.Reversed(), str(out))
    assert _step_volume(out) == pytest.approx(8000.0, rel=1e-6)


_PRUSA_PIN = "529e5ad6892ce3b25015b7cca7f0481b8335a76c"
_PRUSA_FACES = {
    "plug-aligner": {"plane": 14},
    "endstop-block": {"plane": 15, "cylinder": 2},
    "y-belt-idler": {"plane": 197, "cylinder": 6, "cone": 3},
    "y-belt-holder": {"plane": 210, "cylinder": 4},
    "y-motor-holder": {"plane": 396, "cylinder": 6, "cone": 2},
    "y-rod-holder": {"plane": 158, "cylinder": 8, "cone": 11},
}


def _prusa_stl(name):
    cache = Path(__file__).resolve().parents[1] / "benchmarks" / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{name}.stl"
    if path.is_file() and path.stat().st_size > 0:
        return path
    url = (
        "https://raw.githubusercontent.com/prusa3d/Original-Prusa-i3/"
        f"{_PRUSA_PIN}/Printed-Parts/STL/{name}.stl"
    )
    urllib.request.urlretrieve(url, path)
    return path


@pytest.mark.slow
@pytest.mark.parametrize("name", list(_PRUSA_FACES))
def test_prusa_step_stays_exact(name):
    from mesh2cad.io import load_mesh
    mesh = load_mesh(str(_prusa_stl(name)))
    report = convert_mesh(mesh).report
    assert report["method"] == "exact"
    assert report["output"]["valid_solid"]
    assert report["output"]["faces"] == _PRUSA_FACES[name]
    assert report["volume_error_pct"] < 1.0
    assert report["deviation"]["hausdorff"] < 0.5
