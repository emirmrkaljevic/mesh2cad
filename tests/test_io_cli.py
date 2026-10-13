import json

from mesh2cad.cli import main
from mesh2cad.fit import fit_cylinder, fit_plane
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


def test_cli_convert(tmp_path):
    stl = tmp_path / "cube.stl"
    to_mesh(cube()).export(stl)
    out = tmp_path / "cube.step"
    rep = tmp_path / "r.json"
    assert main(["convert", str(stl), "-o", str(out), "--report", str(rep)]) == 0
    assert out.read_text().startswith("ISO-10303-21")
    r = json.loads(rep.read_text())
    assert r["output"]["faces"] == {"plane": 6} and r["output"]["valid_solid"]
