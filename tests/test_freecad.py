"""Parametric FreeCAD export: the macro runs, matches the part, and edits apply.

Text checks run without FreeCAD. Headless checks skip when ``freecadcmd`` is
not on PATH.
"""
import math
import re
import shutil

import pytest

pytest.importorskip("cadquery")
from mesh2cad.cadquery_code import to_cadquery  # noqa: E402
from mesh2cad.cli import main  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.freecad_code import check, freecad_bin, render_mesh, to_freecad  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.parametric import Unsupported  # noqa: E402
from parts import (  # noqa: E402
    angled_hole, boss, chamfered_boss, counterbored_plate, countersunk_plate,
    cube, cut_block, filleted_box, filleted_bracket, filleted_plate,
    l_profile_y, orthogonal_intersection, partial_edge_chamfer, plate_with_hole,
    ribbed_plate, side_hole_boss, slanted_outline, slotted_plate, spherical_dome,
    stepped_block, tapered_boss, to_mesh,
)
from test_csg import _chamfered_side_hole, _stepped_shelf  # noqa: E402
from test_partial_arc import _breakout_part, _mesh_of  # noqa: E402
from test_sweep import _mesh as _sweep_mesh  # noqa: E402
from test_sweep import _ramped_bar  # noqa: E402

needs_freecad = pytest.mark.skipif(freecad_bin() is None, reason="freecadcmd is not installed")

BOX = {
    "cube": cube,
    "stepped_block": stepped_block,
    "plate_with_hole": plate_with_hole,
    "boss": boss,
    "ribbed_plate": ribbed_plate,
    "countersunk_plate": countersunk_plate,
    "tapered_boss": tapered_boss,
    "chamfered_boss": chamfered_boss,
}

EXPECT = {
    "plate_with_hole": {"hole_1_d": 10.0, "base_t": 6.0},
    "countersunk_plate": {"hole_1_d": 6.0, "hole_1_csk_d": 12.0, "hole_1_csk_angle": 90.0},
    "boss": {"boss_1_d": 16.0, "boss_1_h": 15.0},
    "chamfered_boss": {"boss_1_d": 16.0, "boss_1_chamfer": 2.0},
    "tapered_boss": {"boss_1_d": 16.0, "boss_1_draft_deg": 10.0},
    "ribbed_plate": {"pad_1_w": 14.0, "hole_2_d": 8.0},
    "stepped_block": {"pad_1_w": 20.0, "pad_1_h": 8.0},
}

SHAPE = {
    "cube": ("Part.makeBox(base_x, base_y, base_t,",),
    "stepped_block": ("Part.makeBox(pad_1_w, pad_1_l, pad_1_h,", ".fuse("),
    "plate_with_hole": (".cut(", "Part.makeCylinder(hole_1_d / 2, base_t + 2,"),
    "boss": ("Part.makeCylinder(boss_1_d / 2, boss_1_h,", ".fuse("),
    "ribbed_plate": (".cut(", ".fuse(", "Part.makeCylinder(hole_1_d / 2, base_t + 2,"),
    "countersunk_plate": (".cut(", "Part.makeCone(hole_1_d / 2,", "hole_1_csk_angle"),
    "tapered_boss": ("Part.makeCone(boss_1_d / 2,", "boss_1_draft_deg"),
    "chamfered_boss": ("Part.makeCone(boss_1_d / 2, boss_1_d / 2 - boss_1_chamfer,",),
}


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _fit(fn, sigma=0.0, tol=None):
    ref = to_mesh(fn())
    tol = max(sigma, 0.01) if tol is None else tol
    mesh = add_noise(ref, sigma, seed=0) if sigma else ref
    res = convert_mesh(mesh, tol=tol, with_report=False)
    return ref, res


def _against(kind, sigma):
    """Hausdorff and volume limits the other exporters use against the clean mesh."""
    if kind == "fillet":
        return (0.08 if sigma == 0.0 else 0.25), 0.5
    if kind == "side":
        return 0.08, 0.5
    if kind == "extrude":
        return (0.05 if sigma == 0.0 else 0.25), 0.5
    if kind == "csg":
        return (0.05 if sigma == 0.0 else 0.3), 0.5
    if kind == "stack":
        return (0.08 if sigma == 0.0 else 0.35), 0.5
    if kind == "tight":
        return 0.08, 0.5
    return max(3 * sigma, 0.05), 0.5


@needs_freecad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("name", list(BOX))
def test_script_reproduces_part(name, sigma, tmp_path):
    ref, res = _fit(BOX[name], sigma)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * max(sigma, 0.01), 0.1)
    assert q["volume_error_pct"] <= 1.0
    for bit in SHAPE[name]:
        assert bit in code
    dev, vol = _against("box", sigma)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < dev
    assert against["volume_error_pct"] < vol
    for k, v in EXPECT.get(name, {}).items():
        assert _params(code)[k] == pytest.approx(v, abs=0.5 if ("angle" in k or "deg" in k) else 0.05)


@pytest.mark.parametrize("name", list(BOX) + ["filleted_plate", "filleted_box"])
def test_parameters_match_cadquery(name):
    fn = {**BOX, "filleted_plate": filleted_plate, "filleted_box": filleted_box}[name]
    ref = to_mesh(fn())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    fc, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert _params(fc) == _params(cq)
    compile(fc, "part.FCMacro", "exec")
    assert "# --- parameters ---" in fc and "# --- model ---" in fc
    assert "import FreeCAD as App" in fc and "import Part" in fc


def test_parameters_match_on_extrusions_and_booleans():
    for fn in (l_profile_y, slanted_outline, slotted_plate, orthogonal_intersection,
               cut_block, counterbored_plate, _stepped_shelf, side_hole_boss):
        ref = to_mesh(fn())
        res = convert_mesh(ref, tol=0.01, with_report=False)
        cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
        fc, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
        assert _params(fc) == _params(cq)
        compile(fc, "part.FCMacro", "exec")


def test_profile_and_boolean_text():
    ref = to_mesh(l_profile_y())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert '_extrude(outline, "Y", extrude_0, extrude_len)' in code
    # the through hole on this L runs along X, not along the extrusion
    assert "App.Vector(1, 0, 0)" in code

    ref = to_mesh(slanted_outline())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert '_extrude(outline, "Z", extrude_0, extrude_len)' in code

    ref = to_mesh(slotted_plate())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert "cutout_1" in code and ".cut(_extrude(cutout_1," in code

    ref = to_mesh(orthogonal_intersection())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert ".common(" in code
    assert "extrude_z_len" in code and "outline_z" in code

    ref = to_mesh(cut_block())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert "cut_1_len" in code and "result = result.cut(cut_1_solid)" in code

    ref = to_mesh(filleted_plate())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert "result.makeFillet(fillet_r, _fillet_edges)" in code

    ref = to_mesh(side_hole_boss())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, _ = to_freecad(res.mesh, res.regions, res.tol, verify=False)
    assert "App.Vector(1, 0, 0)" in code
    assert "App.Vector(0, 1, 0)" in code


@needs_freecad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("name,fn,radius", [
    ("filleted_plate", filleted_plate, 5.0),
    ("filleted_box", filleted_box, 2.0),
])
def test_vertical_fillet_is_a_parameter(name, fn, radius, sigma, tmp_path):
    ref, res = _fit(fn, sigma)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert _params(code)["fillet_r"] == pytest.approx(radius, abs=0.15)
    assert "makeFillet(fillet_r, _fillet_edges)" in code
    dev, vol = _against("fillet", sigma)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < dev
    assert against["volume_error_pct"] < vol


@needs_freecad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
def test_side_features(sigma, tmp_path):
    ref, res = _fit(side_hole_boss, sigma)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    dev, vol = _against("side", sigma)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < dev
    assert against["volume_error_pct"] < vol


@needs_freecad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("fn", [l_profile_y, slanted_outline, slotted_plate])
def test_extrusion(fn, sigma, tmp_path):
    ref, res = _fit(fn, sigma)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    dev, vol = _against("extrude", sigma)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < dev
    assert against["volume_error_pct"] < vol


@needs_freecad
@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("fn,kind", [
    (orthogonal_intersection, "csg"),
    (cut_block, "csg"),
    (counterbored_plate, "csg"),
    (_stepped_shelf, "stack"),
])
def test_csg(fn, kind, sigma, tmp_path):
    ref, res = _fit(fn, sigma)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    dev, vol = _against(kind, sigma)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < dev
    assert against["volume_error_pct"] < vol


@needs_freecad
def test_partial_edge_chamfer(tmp_path):
    ref = _mesh_of(partial_edge_chamfer())
    res = convert_mesh(ref, tol=0.05, with_report=False)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert "cut_1" in code
    dev, vol = _against("tight", 0.0)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < dev
    assert against["volume_error_pct"] < vol


@needs_freecad
def test_side_chamfer(tmp_path):
    ref = _mesh_of(_chamfered_side_hole())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert _params(code)["hole_1_chamfer"] == pytest.approx(0.8, abs=0.05)
    assert "Part.makeCone" in code
    dev, vol = _against("tight", 0.0)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < dev
    assert against["volume_error_pct"] < vol


@needs_freecad
def test_breakout_counterbore(tmp_path):
    ref = _mesh_of(_breakout_part())
    res = convert_mesh(ref, tol=0.05, with_report=False)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert "cbore_d" in "".join(params)
    dev, vol = _against("tight", 0.0)
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] < dev
    assert against["volume_error_pct"] < vol


@needs_freecad
def test_sweep_segmented_section(tmp_path):
    ref = _sweep_mesh(_ramped_bar())
    res = convert_mesh(ref, tol=0.05, with_report=False)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    limit = max(4 * res.tol, 0.1)
    assert q["max_deviation"] <= limit
    assert q["volume_error_pct"] <= 1.0
    assert "cut_" in code and "_extrude(" in code
    against = check(code, ref, tmp_path)
    assert against["max_deviation"] <= limit
    assert against["volume_error_pct"] <= 1.0


@needs_freecad
def test_edit_hole_diameter(tmp_path):
    ref, res = _fit(plate_with_hole)
    code, _ = to_freecad(res.mesh, res.regions, res.tol)
    before = render_mesh(code, tmp_path / "before")
    bigger = re.sub(r"^(hole_1_d\s*=\s*)[\d.]+", r"\g<1>12.0", code, flags=re.M)
    after = render_mesh(bigger, tmp_path / "after")
    expect = math.pi * (6.0 ** 2 - 5.0 ** 2) * 6.0
    assert abs(before.volume) - abs(after.volume) == pytest.approx(expect, rel=0.02)


def test_bracket_fillet_is_not_a_box_script():
    ref = to_mesh(filleted_bracket())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    with pytest.raises(Unsupported):
        to_freecad(res.mesh, res.regions, res.tol, verify=False)


def test_unsupported_part():
    res = convert_mesh(to_mesh(spherical_dome()), tol=0.01, with_report=False)
    with pytest.raises(Unsupported):
        to_freecad(res.mesh, res.regions, res.tol, verify=False)


def test_angled_axis_is_refused():
    ref = to_mesh(angled_hole())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    with pytest.raises(Unsupported):
        to_freecad(res.mesh, res.regions, res.tol, verify=False)


def test_unverified_when_freecad_missing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("mesh2cad.freecad_code.freecad_bin", lambda: None)
    stl = tmp_path / "p.stl"
    to_mesh(plate_with_hole()).export(stl)
    out = tmp_path / "p.FCMacro"
    assert main(["convert", str(stl), "--freecad", str(out)]) == 0
    text = out.read_text()
    assert "hole_1_d" in text and "Part.makeBox(base_x, base_y, base_t," in text
    assert "unverified" in capsys.readouterr().out


@needs_freecad
def test_cli_freecad(tmp_path, capsys):
    stl = tmp_path / "p.stl"
    to_mesh(countersunk_plate()).export(stl)
    out = tmp_path / "p.FCMacro"
    step = tmp_path / "p.step"
    assert main(["convert", str(stl), "--freecad", str(out), "-o", str(step)]) == 0
    assert "hole_1_csk_d" in out.read_text()
    assert step.is_file() and step.read_text().startswith("ISO-10303-21")
    assert "max_dev=" in capsys.readouterr().out
    dome = tmp_path / "dome.FCMacro"
    to_mesh(spherical_dome()).export(stl)
    assert main(["convert", str(stl), "--freecad", str(dome), "-o", str(tmp_path / "nope.step")]) == 2
    assert not dome.exists()
    assert not (tmp_path / "nope.step").exists()
    assert main(["convert", str(stl), "-o", str(tmp_path / "dome.step")]) == 0
    assert (tmp_path / "dome.step").is_file()


def test_mcp_to_freecad(tmp_path):
    pytest.importorskip("mcp")
    import asyncio

    from mesh2cad.mcp_server import mcp
    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert {"convert", "analyze", "report", "to_cadquery", "to_openscad", "to_freecad"} <= names
    from mesh2cad.mcp_server import to_freecad as tool
    stl = tmp_path / "b.stl"
    to_mesh(boss()).export(stl)
    r = tool(str(stl))
    assert "error" not in r
    assert r["parameters"]["boss_1_d"] == pytest.approx(16.0, abs=0.05)
    assert (tmp_path / "b.FCMacro").is_file()
    if shutil.which("freecadcmd"):
        assert r["verified"] is True and r["check"]["max_deviation"] <= 0.1
    else:
        assert r["verified"] is False and r["check"] is None
    bad = tmp_path / "dome.stl"
    to_mesh(spherical_dome()).export(bad)
    err = tool(str(bad), str(tmp_path / "dome.FCMacro"))
    assert "error" in err
    assert not (tmp_path / "dome.FCMacro").exists()
