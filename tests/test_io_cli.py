import json

from mesh2cad.cli import main
from mesh2cad.fit import fit_cylinder, fit_plane, fit_sphere, refine_sphere
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
