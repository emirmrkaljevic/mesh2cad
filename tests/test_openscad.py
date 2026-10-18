"""Parametric OpenSCAD export: the script renders, matches the part, and edits apply."""
import re
import shutil

import pytest

pytest.importorskip("cadquery")
from mesh2cad.cadquery_code import to_cadquery  # noqa: E402
from mesh2cad.cli import main  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.openscad_code import check, render_mesh, to_openscad  # noqa: E402
from mesh2cad.parametric import Unsupported  # noqa: E402
from parts import (boss, chamfered_boss, countersunk_plate, cube, filleted_box,  # noqa: E402
                   filleted_bracket, filleted_plate, plate_with_hole, ribbed_plate,
                   spherical_dome, stepped_block, tapered_boss, to_mesh)

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")

PARTS = {"cube": cube, "stepped_block": stepped_block, "plate_with_hole": plate_with_hole,
         "boss": boss, "ribbed_plate": ribbed_plate, "countersunk_plate": countersunk_plate,
         "tapered_boss": tapered_boss, "chamfered_boss": chamfered_boss}

EXPECT = {"plate_with_hole": {"hole_1_d": 10.0, "base_t": 6.0},
          "countersunk_plate": {"hole_1_d": 6.0, "hole_1_csk_d": 12.0, "hole_1_csk_angle": 90.0},
          "boss": {"boss_1_d": 16.0, "boss_1_h": 15.0},
          "chamfered_boss": {"boss_1_d": 16.0, "boss_1_chamfer": 2.0},
          "tapered_boss": {"boss_1_d": 16.0, "boss_1_draft_deg": 10.0},
          "ribbed_plate": {"pad_1_w": 14.0, "hole_2_d": 8.0},
          "stepped_block": {"pad_1_w": 20.0, "pad_1_h": 8.0}}

SHAPE = {
    "cube": ("cube([base_x, base_y, base_t], center=true);",),
    "stepped_block": ("cube([pad_1_w, pad_1_l, pad_1_h], center=true);", "union()"),
    "plate_with_hole": ("difference()", "cylinder(h=base_t + 2, d=hole_1_d);"),
    "boss": ("cylinder(h=boss_1_h, d=boss_1_d);", "union()"),
    "ribbed_plate": ("difference()", "union()", "cylinder(h=base_t + 2, d=hole_1_d);"),
    "countersunk_plate": ("difference()", "r1=hole_1_d / 2", "tan(hole_1_csk_angle / 2)"),
    "tapered_boss": ("r1=boss_1_d / 2", "r2=boss_1_d / 2 - boss_1_h * tan(boss_1_draft_deg)"),
    "chamfered_boss": ("r1=boss_1_d / 2", "r2=boss_1_d / 2 - boss_1_chamfer"),
}


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _gen(name, sigma=0.0, fn=None):
    fn = fn or PARTS[name]
    ref = to_mesh(fn())
    tol = max(sigma, 0.01)
    res = convert_mesh(add_noise(ref, sigma, seed=0), tol=tol, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    return ref, code, q


@needs_openscad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("name", list(PARTS))
def test_script_reproduces_part(name, sigma, tmp_path):
    ref, code, q = _gen(name, sigma)
    assert q["max_deviation"] <= max(4 * max(sigma, 0.01), 0.1)
    assert q["volume_error_pct"] <= 1.0
    for bit in SHAPE[name]:
        assert bit in code
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < max(3 * sigma, 0.05)
    assert against["volume_error_pct"] < 0.5
    for k, v in EXPECT.get(name, {}).items():
        assert _params(code)[k] == pytest.approx(v, abs=0.5 if ("angle" in k or "deg" in k) else 0.05)


@pytest.mark.parametrize("name", list(PARTS) + ["filleted_plate", "filleted_box"])
def test_parameters_match_cadquery(name):
    fn = {**PARTS, "filleted_plate": filleted_plate, "filleted_box": filleted_box}[name]
    ref = to_mesh(fn())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    sc, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    assert _params(sc) == _params(cq)


@needs_openscad
@pytest.mark.parametrize("seed", [1, 2])
def test_recovered_chamfer_exports_openscad(seed, tmp_path):
    ref = to_mesh(chamfered_boss())
    res = convert_mesh(add_noise(ref, 0.05, seed=seed), tol=0.05, with_report=False)
    assert any(getattr(r.surface, "kind", None) == "cone" for r in res.regions)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < 0.15
    assert against["volume_error_pct"] < 0.5
    params = _params(code)
    assert params["boss_1_d"] == pytest.approx(16.0, abs=0.05)
    assert params["boss_1_chamfer"] == pytest.approx(2.0, abs=0.2)


@needs_openscad
def test_edit_hole_diameter_changes_the_hole(tmp_path):
    ref, code, _ = _gen("plate_with_hole")
    before = render_mesh(code, tmp_path / "before")
    bigger = re.sub(r"^(hole_1_d\s*=\s*)[\d.]+", r"\g<1>12.0", code, flags=re.M)
    after = render_mesh(bigger, tmp_path / "after")
    import math
    expect = math.pi * (6.0 ** 2 - 5.0 ** 2) * 6.0
    assert before.volume - after.volume == pytest.approx(expect, rel=0.02)
    params = _params(code)
    z = params["base_z0"] + params["base_t"] / 2.0
    ring = [[params["hole_1_x"] + 5.5, params["hole_1_y"], z]]
    assert bool(before.contains(ring)[0])
    assert not bool(after.contains(ring)[0])


@needs_openscad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("name,fn,radius", [
    ("filleted_plate", filleted_plate, 5.0),
    ("filleted_box", filleted_box, 2.0),
])
def test_vertical_fillet_is_a_parameter(name, fn, radius, sigma, tmp_path):
    ref = to_mesh(fn())
    tol = max(sigma, 0.01)
    res = convert_mesh(add_noise(ref, sigma, seed=0), tol=tol, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * tol, 0.1)
    params = _params(code)
    assert params["fillet_r"] == pytest.approx(radius, abs=0.15)
    assert "offset(r=fillet_r)" in code
    assert re.search(r"^fillet_r\s*=\s*[\d.]+;\s*//\s*\[", code, re.M)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < (0.08 if sigma == 0.0 else 0.25)
    assert against["volume_error_pct"] < 0.5


def _slider_hi(code, name):
    m = re.search(rf"^{name}\s*=\s*[\d.]+;\s*//\s*\[([^\]]+)\]", code, re.M)
    assert m, name
    return float(m.group(1).split(":")[-1])


def test_slider_only_where_a_bound_exists():
    ref = to_mesh(countersunk_plate())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    assert re.search(r"^hole_1_d\s*=\s*[\d.]+;\s*//\s*\[", code, re.M)
    assert re.search(r"^hole_1_csk_angle\s*=\s*[\d.]+;\s*//\s*\[60:120\]", code, re.M)
    # position and overall size stay plain numbers
    assert re.search(r"^base_cx\s*=\s*[-\d.]+;\s*$", code, re.M)
    assert re.search(r"^base_x\s*=\s*[\d.]+;\s*$", code, re.M)
    assert "/* [Holes] */" in code
    assert "/* [Hidden] */" in code

    ref = to_mesh(filleted_plate())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    params = _params(code)
    assert _slider_hi(code, "fillet_r") < min(params["base_x"], params["base_y"]) / 2

    ref = to_mesh(tapered_boss())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    import math
    params = _params(code)
    hi = _slider_hi(code, "boss_1_draft_deg")
    r2 = params["boss_1_d"] / 2 - params["boss_1_h"] * math.tan(math.radians(hi))
    assert r2 >= -1e-6


@needs_openscad
def test_blind_hole(tmp_path):
    import cadquery as cq
    wp = cq.Workplane().box(40, 30, 10).faces(">Z").workplane().hole(8, 4)
    ref = to_mesh(wp)
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["hole_1_d"] == pytest.approx(8.0, abs=0.05)
    assert params["hole_1_depth"] == pytest.approx(4.0, abs=0.05)
    assert "cylinder(h=hole_1_depth + 1, d=hole_1_d);" in code
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < 0.05
    assert against["volume_error_pct"] < 0.5


def test_bracket_fillet_is_not_a_box_script():
    ref = to_mesh(filleted_bracket())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    with pytest.raises(Unsupported):
        to_openscad(res.mesh, res.regions, res.tol)


def test_unsupported_part():
    res = convert_mesh(to_mesh(spherical_dome()), tol=0.01, with_report=False)
    with pytest.raises(Unsupported):
        to_openscad(res.mesh, res.regions, res.tol)


def test_unverified_when_openscad_missing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("mesh2cad.openscad_code.openscad_bin", lambda: None)
    stl = tmp_path / "p.stl"
    to_mesh(plate_with_hole()).export(stl)
    out = tmp_path / "p.scad"
    assert main(["convert", str(stl), "--openscad", str(out)]) == 0
    text = out.read_text()
    assert "hole_1_d" in text and "cube([base_x, base_y, base_t], center=true);" in text
    assert "unverified" in capsys.readouterr().out


@needs_openscad
def test_cli_openscad(tmp_path, capsys):
    stl = tmp_path / "p.stl"
    to_mesh(countersunk_plate()).export(stl)
    out = tmp_path / "p.scad"
    step = tmp_path / "p.step"
    assert main(["convert", str(stl), "--openscad", str(out), "-o", str(step)]) == 0
    assert "hole_1_csk_d" in out.read_text()
    assert step.is_file() and step.read_text().startswith("ISO-10303-21")
    assert "max_dev=" in capsys.readouterr().out
    dome = tmp_path / "dome.scad"
    to_mesh(spherical_dome()).export(stl)
    assert main(["convert", str(stl), "--openscad", str(dome), "-o", str(tmp_path / "nope.step")]) == 2
    assert not dome.exists()
    assert not (tmp_path / "nope.step").exists()
    # STEP on its own still accepts a part the script vocabulary refuses
    assert main(["convert", str(stl), "-o", str(tmp_path / "dome.step")]) == 0
    assert (tmp_path / "dome.step").is_file()


def test_mcp_to_openscad(tmp_path):
    pytest.importorskip("mcp")
    import asyncio

    from mesh2cad.mcp_server import mcp
    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert {"convert", "analyze", "report", "to_cadquery", "to_openscad"} <= names
    from mesh2cad.mcp_server import to_openscad as tool
    stl = tmp_path / "b.stl"
    to_mesh(boss()).export(stl)
    r = tool(str(stl))
    assert "error" not in r
    assert r["parameters"]["boss_1_d"] == pytest.approx(16.0, abs=0.05)
    assert (tmp_path / "b.scad").is_file()
    if shutil.which("openscad"):
        assert r["verified"] is True and r["check"]["max_deviation"] <= 0.1
    else:
        assert r["verified"] is False and r["check"] is None
    bad = tmp_path / "dome.stl"
    to_mesh(spherical_dome()).export(bad)
    err = tool(str(bad), str(tmp_path / "dome.scad"))
    assert "error" in err
    assert not (tmp_path / "dome.scad").exists()
