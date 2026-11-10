"""Through-holes in the foot and upright of a filleted L-bracket.

Fixtures are built in this process. Nothing is loaded from disk. The pass
sets are the seeds the deviation gate accepts. A miss is a refusal, not a
script with the hole left out.
"""
import math
import re
import shutil

import pytest

pytest.importorskip("cadquery")
from mesh2cad.cadquery_code import to_cadquery  # noqa: E402
from mesh2cad.cli import main  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.freecad_code import freecad_bin, to_freecad  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.openscad_code import to_openscad  # noqa: E402
from mesh2cad.parametric import Unsupported  # noqa: E402
from mesh2cad.params import GeometryError, edit_script, list_parameters  # noqa: E402
from parts import filleted_bracket_holes, filleted_bracket_two_base_holes, to_mesh  # noqa: E402

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")
needs_freecad = pytest.mark.skipif(freecad_bin() is None, reason="freecadcmd is not installed")

# Foot hole at (24, 4), upright hole at (4, 19), both Ø4. Same bracket as the
# filleted L without holes: 40×30×12, legs 8, r=3.
_LAYOUT = {"both": (True, True), "base": (True, False), "leg": (False, True)}

# CadQuery and the FreeCAD proxy (CadQuery, when freecadcmd is absent).
_CQ_PASS = {
    "both": {(0.0, 0), (0.0, 1), (0.0, 2), (0.05, 1), (0.05, 2)},
    "base": {(0.0, 0), (0.0, 1), (0.0, 2), (0.05, 0)},
    "leg": {(0.0, 0), (0.0, 1), (0.0, 2), (0.05, 2)},
}
# OpenSCAD agrees, except the foot-only part at 0.05 mm seed 2: CadQuery is
# 0.211 mm (limit 0.200) and OpenSCAD is 0.188 mm.
_OS_PASS = {
    "both": _CQ_PASS["both"],
    "base": _CQ_PASS["base"] | {(0.05, 2)},
    "leg": _CQ_PASS["leg"],
}
# A script was produced and the deviation gate refused it.
_CQ_GATE = {"both": set(), "base": {(0.05, 1), (0.05, 2)}, "leg": {(0.05, 1)}}
_OS_GATE = {"both": set(), "base": {(0.05, 1)}, "leg": {(0.05, 1)}}

_MESH = {}


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _reference(layout):
    if layout not in _MESH:
        _MESH[layout] = to_mesh(filleted_bracket_holes(*_LAYOUT[layout]))
    return _MESH[layout]


def _case(layout, sigma, seed):
    ref = _reference(layout)
    tol = max(sigma, 0.01)
    mesh = add_noise(ref, sigma, seed=seed) if sigma else ref
    res = convert_mesh(mesh, tol=tol, with_report=False)
    return res, tol


def _expect_holes(params, layout):
    assert params["fillet_r"] == pytest.approx(3.0, abs=0.05)
    assert params["base_x"] == pytest.approx(40.0, abs=0.15)
    assert params["base_y"] == pytest.approx(30.0, abs=0.15)
    assert params["base_t"] == pytest.approx(12.0, abs=0.15)
    assert params["leg_x"] == pytest.approx(8.0, abs=0.15)
    assert params["leg_y"] == pytest.approx(8.0, abs=0.15)
    assert "edge_fillet_r" not in params
    if layout in ("both", "base"):
        assert params["base_hole_1_d"] == pytest.approx(4.0, abs=0.05)
        assert params["base_hole_1_x"] == pytest.approx(24.0, abs=0.15)
        assert params["base_hole_1_y"] == pytest.approx(4.0, abs=0.15)
    else:
        assert not any(name.startswith("base_hole_") for name in params)
    if layout in ("both", "leg"):
        assert params["leg_hole_1_d"] == pytest.approx(4.0, abs=0.05)
        assert params["leg_hole_1_x"] == pytest.approx(4.0, abs=0.15)
        assert params["leg_hole_1_y"] == pytest.approx(19.0, abs=0.15)
    else:
        assert not any(name.startswith("leg_hole_") for name in params)


def _run(export, res, passed, gate, layout, sigma, seed):
    key = (sigma, seed)
    if key not in passed:
        with pytest.raises(Unsupported) as exc:
            export(res.mesh, res.regions, res.tol)
        if key in gate:
            assert "deviates by" in str(exc.value)
        return None
    code, q = export(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    _expect_holes(_params(code), layout)
    return code


@pytest.mark.parametrize("layout", ["both", "base", "leg"])
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.2])
def test_cadquery_bracket_holes(layout, sigma, seed):
    res, _tol = _case(layout, sigma, seed)
    code = _run(to_cadquery, res, _CQ_PASS[layout], _CQ_GATE[layout], layout, sigma, seed)
    if code is None:
        return
    assert ".edges(\"|Z\").fillet(fillet_r)" in code
    hole = "base_hole_1_d" if layout != "leg" else "leg_hole_1_d"
    assert code.index(".edges(\"|Z\").fillet(fillet_r)") < code.index(f".circle({hole} / 2)")


@needs_openscad
@pytest.mark.parametrize("layout", ["both", "base", "leg"])
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.2])
def test_openscad_bracket_holes(layout, sigma, seed):
    res, _tol = _case(layout, sigma, seed)
    code = _run(to_openscad, res, _OS_PASS[layout], _OS_GATE[layout], layout, sigma, seed)
    if code is None:
        return
    assert "difference()" in code
    assert "offset(r=-fillet_r)" in code
    hole = "base_hole_1_d" if layout != "leg" else "leg_hole_1_d"
    assert code.index("offset(r=-fillet_r)") < code.index(f"cylinder(h=base_t + 2, d={hole});")
    slider = re.search(rf"^{hole}\s*=\s*[\d.]+;\s*//\s*\[([^\]]+)\]", code, re.M)
    assert slider
    assert float(slider.group(1).split(":")[-1]) <= min(8.0, 40.0)


@needs_freecad
@pytest.mark.parametrize("layout", ["both", "base", "leg"])
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("sigma", [0.0, 0.05, 0.2])
def test_freecad_bracket_holes(layout, sigma, seed):
    res, _tol = _case(layout, sigma, seed)
    code = _run(to_freecad, res, _CQ_PASS[layout], _CQ_GATE[layout], layout, sigma, seed)
    if code is None:
        return
    assert "result.makeFillet(fillet_r, _fillet_edges)" in code
    hole = "base_hole_1_d" if layout != "leg" else "leg_hole_1_d"
    assert code.index("result.makeFillet(fillet_r, _fillet_edges)") < code.index(
        f"Part.makeCylinder({hole} / 2, base_t + 2,")


def test_two_foot_holes_are_indexed():
    ref = to_mesh(filleted_bracket_two_base_holes())
    res = convert_mesh(ref, tol=0.01, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["base_hole_1_d"] == pytest.approx(4.0, abs=0.05)
    assert params["base_hole_1_x"] == pytest.approx(18.0, abs=0.15)
    assert params["base_hole_2_d"] == pytest.approx(4.0, abs=0.05)
    assert params["base_hole_2_x"] == pytest.approx(30.0, abs=0.15)
    assert params["base_hole_1_y"] == pytest.approx(4.0, abs=0.15)
    assert params["base_hole_2_y"] == pytest.approx(4.0, abs=0.15)
    assert not any(name.startswith("leg_hole_") for name in params)


def test_freecad_macro_text_and_proxy_gate(tmp_path, monkeypatch, capsys):
    """Without freecadcmd the macro is gated on the CadQuery check."""
    monkeypatch.setattr("mesh2cad.freecad_code.freecad_bin", lambda: None)
    res, tol = _case("both", 0.0, 0)
    code, q = to_freecad(res.mesh, res.regions, res.tol)
    assert q is None
    params = _params(code)
    _expect_holes(params, "both")
    assert code.index("result.makeFillet(fillet_r, _fillet_edges)") < code.index(
        "Part.makeCylinder(base_hole_1_d / 2, base_t + 2,")
    stl = tmp_path / "ok.stl"
    res.mesh.export(stl)
    out = tmp_path / "ok.FCMacro"
    assert main(["convert", str(stl), "--freecad", str(out), "--tol", str(tol)]) == 0
    assert "base_hole_1_d" in out.read_text(encoding="utf-8")
    assert "unverified" in capsys.readouterr().out

    # Foot-only, 0.05 mm, seed 2: OpenSCAD would pass. CadQuery is over the limit,
    # so the macro is not written.
    refused, rtol = _case("base", 0.05, 2)
    with pytest.raises(Unsupported, match="deviates by"):
        to_freecad(refused.mesh, refused.regions, refused.tol)
    bad = tmp_path / "bad.stl"
    refused.mesh.export(bad)
    macro = tmp_path / "bad.FCMacro"
    assert main(["convert", str(bad), "--freecad", str(macro), "--tol", str(rtol)]) == 2
    assert not macro.exists()


def test_edit_hole_diameter_and_position(tmp_path, capsys):
    ref = _reference("both")
    stl = tmp_path / "b.stl"
    ref.export(stl)
    script = tmp_path / "b.py"
    assert main(["convert", str(stl), "--cadquery", str(script)]) == 0
    assert main(["params", str(script)]) == 0
    listed = capsys.readouterr().out
    assert "base_hole_1_d" in listed and "diameter" in listed
    assert "leg_hole_1_d" in listed
    info = list_parameters(str(script))
    names = {p["name"]: p for p in info["parameters"]}
    assert names["base_hole_1_d"]["kind"] == "diameter" and names["base_hole_1_d"]["unit"] == "mm"
    assert names["base_hole_1_x"]["kind"] == "position"
    assert names["leg_hole_1_y"]["kind"] == "position"
    assert names["base_hole_1_d"]["value"] == pytest.approx(4.0)
    assert names["leg_hole_1_d"]["value"] == pytest.approx(4.0)
    raw = script.read_bytes()
    with pytest.raises(GeometryError, match="base_hole_1_d"):
        edit_script(str(script), {"base_hole_1_d": 12}, in_place=True)
    assert script.read_bytes() == raw
    with pytest.raises(GeometryError, match="does not fit"):
        edit_script(str(script), {"base_hole_1_y": 20}, in_place=True)
    assert script.read_bytes() == raw
    out = tmp_path / "out.step"
    assert main(["edit", str(script), "--set", "base_hole_1_d=5", "--set", "leg_hole_1_d=+1",
                 "--in-place", "-o", str(out)]) == 0
    report = capsys.readouterr().out
    assert "base_hole_1_d: 4.0 -> 5.0" in report
    assert "leg_hole_1_d: 4.0 -> 5.0" in report
    text = script.read_text(encoding="utf-8")
    assert re.search(r"^base_hole_1_d\s*=\s*5\.0", text, re.M)
    assert re.search(r"^leg_hole_1_d\s*=\s*5\.0", text, re.M)
    assert out.is_file()
    # Both holes grow from Ø4 to Ø5 through the 12 mm height.
    expect = 2.0 * math.pi * (2.5 ** 2 - 2.0 ** 2) * 12.0
    delta = float(re.search(r"delta ([+-][\d.]+)", report).group(1))
    assert delta == pytest.approx(-expect, rel=0.02)
    moved = edit_script(str(script), {"base_hole_1_x": "+2"}, output_path=str(tmp_path / "moved.step"),
                        in_place=True)
    assert moved["valid_solid"] is True
    assert re.search(r"^base_hole_1_x\s*=\s*26\.0", script.read_text(encoding="utf-8"), re.M)


def test_mcp_edits_bracket_hole(tmp_path):
    pytest.importorskip("mcp")
    from mesh2cad.mcp_server import edit_params, list_params
    ref = _reference("base")
    stl = tmp_path / "b.stl"
    ref.export(stl)
    script = tmp_path / "b.py"
    assert main(["convert", str(stl), "--cadquery", str(script)]) == 0
    listed = list_params(str(script))
    names = {p["name"]: p for p in listed["parameters"]}
    assert names["base_hole_1_d"]["kind"] == "diameter"
    assert names["base_hole_1_d"]["value"] == pytest.approx(4.0)
    assert "leg_hole_1_d" not in names
    out = tmp_path / "edited.step"
    result = edit_params(str(script), {"base_hole_1_d": 5}, str(out))
    assert "error" not in result
    assert result["valid_solid"] is True
    assert out.is_file()
    assert re.search(r"^base_hole_1_d\s*=\s*5\.0", script.read_text(encoding="utf-8"), re.M)
    refused = edit_params(str(script), {"base_hole_1_d": 12})
    assert "base_hole_1_d" in refused["error"]
    assert re.search(r"^base_hole_1_d\s*=\s*5\.0", script.read_text(encoding="utf-8"), re.M)
