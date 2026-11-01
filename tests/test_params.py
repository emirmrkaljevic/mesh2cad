"""Edit named parameters on scripts mesh2cad generated, then re-render."""
import json
import os
import re
import shutil
import urllib.request
from pathlib import Path

import pytest

pytest.importorskip("cadquery")
from mesh2cad.cadquery_code import run_script  # noqa: E402
from mesh2cad.cli import main  # noqa: E402
from mesh2cad.convert import cadquery_file, freecad_file, openscad_file  # noqa: E402
from mesh2cad.fit import Cylinder  # noqa: E402
from mesh2cad.freecad_code import freecad_bin  # noqa: E402
from mesh2cad.io import clean  # noqa: E402
from mesh2cad.params import (GeometryError, ParamError, check_edit,  # noqa: E402
                             edit_script, list_parameters)
from mesh2cad.segment import segment  # noqa: E402
from parts import (boss, chamfered_boss, counterbored_plate, countersunk_plate,  # noqa: E402
                   cube, cut_block, filleted_plate, l_profile_y, orthogonal_intersection,
                   plate_with_hole, ribbed_plate, side_hole_boss, slanted_outline,
                   slotted_plate, stepped_block, tapered_boss, to_mesh)

needs_openscad = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad is not installed")
needs_freecad = pytest.mark.skipif(freecad_bin() is None, reason="freecadcmd is not installed")

SYNTHETIC = {
    "cube": cube,
    "stepped_block": stepped_block,
    "plate_with_hole": plate_with_hole,
    "boss": boss,
    "ribbed_plate": ribbed_plate,
    "countersunk_plate": countersunk_plate,
    "tapered_boss": tapered_boss,
    "chamfered_boss": chamfered_boss,
    "filleted_plate": filleted_plate,
    "side_hole_boss": side_hole_boss,
    "l_profile_y": l_profile_y,
    "slanted_outline": slanted_outline,
    "slotted_plate": slotted_plate,
    "orthogonal_intersection": orthogonal_intersection,
    "cut_block": cut_block,
    "counterbored_plate": counterbored_plate,
}

PRUSA = (
    "plug-aligner",
    "endstop-block",
    "y-belt-idler",
    "y-belt-holder",
    "y-motor-holder",
    "y-rod-holder",
)
_PRUSA_PIN = "529e5ad6892ce3b25015b7cca7f0481b8335a76c"
_NUMBER = r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?"


def _mesh(tmp, fn, name):
    path = Path(tmp) / f"{name}.stl"
    to_mesh(fn()).export(path)
    return path


def _script(tmp, mesh, kind):
    suffix = {"cadquery": ".py", "openscad": ".scad", "freecad": ".FCMacro"}[kind]
    dest = Path(tmp) / f"part{suffix}"
    {"cadquery": cadquery_file, "openscad": openscad_file, "freecad": freecad_file}[kind](str(mesh), str(dest))
    return dest


def _prusa_stl(name):
    cache = Path(__file__).resolve().parents[1] / "benchmarks" / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{name}.stl"
    if not path.is_file() or path.stat().st_size == 0:
        url = (
            "https://raw.githubusercontent.com/prusa3d/Original-Prusa-i3/"
            f"{_PRUSA_PIN}/Printed-Parts/STL/{name}.stl"
        )
        urllib.request.urlretrieve(url, path)
    return path


def _blank_value(line):
    return re.sub(rf"(=\s*){_NUMBER}", r"\1#", line, count=1)


def _one_number_changed(before, after, name):
    b_lines, a_lines = before.splitlines(), after.splitlines()
    assert b_lines != a_lines
    assert len(b_lines) == len(a_lines)
    diffs = [(b, a) for b, a in zip(b_lines, a_lines) if b != a]
    assert len(diffs) == 1, diffs
    old, new = diffs[0]
    assert name in old and name in new
    assert _blank_value(old) == _blank_value(new)


def _candidates(info):
    params = {p["name"]: p for p in info["parameters"]}
    holes = [name for name in params if re.fullmatch(r"hole_\d+_d", name)]
    for name in holes:
        old = float(params[name]["value"])
        for new in (round(old + 1.0, 4), round(old + 0.5, 4), round(old - 0.5, 4), round(old - 1.0, 4)):
            if new > 0.4 and abs(new - old) >= 0.2:
                yield name, new
    for name in ("base_x", "base_y", "base_t", "extrude_len"):
        if name in params and isinstance(params[name]["value"], float):
            yield name, "+2"
            yield name, "+1"
    for name, param in params.items():
        if param["kind"] == "length" and name.endswith("_len") and name not in ("extrude_len",):
            yield name, "+1"


def _choose(path):
    info = list_parameters(str(path))
    for name, new in _candidates(info):
        try:
            rows = check_edit(str(path), {name: new})
        except (GeometryError, ParamError):
            continue
        return info, name, new, rows[0]["old"], rows[0]["new"]
    raise AssertionError(f"no safe edit for {path}")


def _length_axis(name, info):
    params = {p["name"]: p for p in info["parameters"]}
    if name == "base_x":
        return 0
    if name == "base_y":
        return 1
    if name == "base_t":
        return 2
    if name.startswith("extrude_x") and name.endswith("_len"):
        return 0
    if name.startswith("extrude_y") and name.endswith("_len"):
        return 1
    if name.startswith("extrude_z") and name.endswith("_len"):
        return 2
    if name == "extrude_len":
        note = (params.get("outline") or {}).get("note") or ""
        if "YZ" in note:
            return 0
        if "XZ" in note:
            return 1
        return 2
    return None


def _cylinder_radii(mesh):
    cleaned = clean(mesh)
    tol = max(0.05, 0.002 * float(cleaned.scale))
    regions = segment(cleaned, tol, cylinders=True, spheres=False, cones=False)
    return sorted(float(r.surface.radius) for r in regions if isinstance(r.surface, Cylinder))


def _assert_geometry(result, info, name, old, new, mesh):
    assert result["valid_solid"] is True
    if re.fullmatch(r"hole_\d+_d", name):
        assert abs(result["delta"]["volume"]) > 0.5
        if new > old:
            assert result["delta"]["volume"] < 0
        else:
            assert result["delta"]["volume"] > 0
        radii = _cylinder_radii(mesh)
        assert any(abs(radius - new / 2.0) <= 0.45 for radius in radii), (name, new, radii)
        return {"radii": radii}
    axis = _length_axis(name, info)
    assert axis is not None, name
    before = result["before"]["bbox_max"][axis] - result["before"]["bbox_min"][axis]
    after = result["after"]["bbox_max"][axis] - result["after"]["bbox_min"][axis]
    assert after - before == pytest.approx(new - old, abs=0.35)
    return {"bbox_span": [before, after], "axis": "XYZ"[axis]}


def _record(row):
    if not os.environ.get("MESH2CAD_EDIT_REPORT"):
        return
    path = Path(__file__).resolve().parents[1] / "benchmarks" / ".cache" / "param-edit-report.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row) + "\n")


def _run_edit(tmp, script, kind, part):
    info, name, raw, old, new = _choose(script)
    if any(p["name"] == "outline" for p in info["parameters"]):
        with pytest.raises(ParamError):
            check_edit(str(script), {"outline": 1})
    original = script.read_bytes()
    edited = Path(tmp) / f"edited{script.suffix}"
    # CadQuery writes an analytic STEP. OpenSCAD and FreeCAD are measured from STL.
    out = Path(tmp) / ("out.step" if kind == "cadquery" else "out.stl")
    result = edit_script(str(script), {name: raw}, output_path=str(out), output_script=str(edited))
    assert script.read_bytes() == original
    _one_number_changed(original.decode(), edited.read_text(encoding="utf-8"), name)
    assert out.is_file() and out.stat().st_size > 0
    if kind == "cadquery":
        mesh = run_script(edited.read_text(encoding="utf-8"))
        assert out.read_text(encoding="utf-8", errors="ignore").startswith("ISO-10303-21")
    else:
        import trimesh
        mesh = trimesh.load(str(out), force="mesh", process=True)
    measured = _assert_geometry(result, info, name, old, new, mesh)
    _record({
        "part": part,
        "backend": kind,
        "parameter": name,
        "old": old,
        "new": new,
        "volume_before": result["before"]["volume"],
        "volume_after": result["after"]["volume"],
        "bbox_before": [result["before"]["bbox_min"], result["before"]["bbox_max"]],
        "bbox_after": [result["after"]["bbox_min"], result["after"]["bbox_max"]],
        "measured": measured,
    })
    return result


@pytest.mark.parametrize("name", list(SYNTHETIC))
def test_cadquery_edit_changes_geometry(name, tmp_path):
    mesh = _mesh(tmp_path, SYNTHETIC[name], name)
    script = _script(tmp_path, mesh, "cadquery")
    _run_edit(tmp_path, script, "cadquery", name)


@pytest.mark.parametrize("name", list(SYNTHETIC))
@needs_openscad
def test_openscad_edit_changes_geometry(name, tmp_path):
    mesh = _mesh(tmp_path, SYNTHETIC[name], name)
    script = _script(tmp_path, mesh, "openscad")
    _run_edit(tmp_path, script, "openscad", name)


@pytest.mark.parametrize("name", list(SYNTHETIC))
@needs_freecad
def test_freecad_edit_changes_geometry(name, tmp_path):
    mesh = _mesh(tmp_path, SYNTHETIC[name], name)
    script = _script(tmp_path, mesh, "freecad")
    _run_edit(tmp_path, script, "freecad", name)


@pytest.mark.parametrize("name", PRUSA)
def test_prusa_cadquery_edit_changes_geometry(name, tmp_path):
    script = _script(tmp_path, _prusa_stl(name), "cadquery")
    _run_edit(tmp_path, script, "cadquery", name)


@pytest.mark.parametrize("name", PRUSA)
@needs_openscad
def test_prusa_openscad_edit_changes_geometry(name, tmp_path):
    script = _script(tmp_path, _prusa_stl(name), "openscad")
    _run_edit(tmp_path, script, "openscad", name)


@pytest.mark.parametrize("name", PRUSA)
@needs_freecad
def test_prusa_freecad_edit_changes_geometry(name, tmp_path):
    script = _script(tmp_path, _prusa_stl(name), "freecad")
    _run_edit(tmp_path, script, "freecad", name)


def test_params_lists_kind_and_edit_is_local(tmp_path, capsys):
    mesh = _mesh(tmp_path, plate_with_hole, "plate")
    script = _script(tmp_path, mesh, "cadquery")
    assert main(["params", str(script)]) == 0
    listed = capsys.readouterr().out
    assert "hole_1_d" in listed and "diameter" in listed and "mm" in listed
    assert "base_x" in listed and "length" in listed
    out = tmp_path / "out.step"
    edited = tmp_path / "edited.py"
    assert main(["edit", str(script), "--set", "hole_1_d=12", "--set", "base_x=+2",
                 "-o", str(out), "--output-script", str(edited)]) == 0
    report = capsys.readouterr().out
    assert "hole_1_d: 10.0 -> 12.0" in report
    assert "base_x: 40.0 -> 42.0" in report
    assert "valid solid: yes" in report
    assert "volume:" in report and "bbox min:" in report
    assert out.read_text(encoding="utf-8", errors="ignore").startswith("ISO-10303-21")
    text = edited.read_text(encoding="utf-8")
    assert "hole_1_d" in text and "12.0" in text
    original = script.read_text(encoding="utf-8")
    assert "hole_1_d" in original and "10.0" in original
    # two numbers change, everything else stays
    diffs = [(a, b) for a, b in zip(original.splitlines(), text.splitlines()) if a != b]
    assert len(diffs) == 2
    assert all(_blank_value(a) == _blank_value(b) for a, b in diffs)


def test_relative_base_changes_bbox(tmp_path):
    mesh = _mesh(tmp_path, cube, "cube")
    script = _script(tmp_path, mesh, "cadquery")
    result = edit_script(str(script), {"base_x": "+2"}, output_path=str(tmp_path / "out.step"),
                         output_script=str(tmp_path / "edited.py"))
    before = result["before"]["bbox_max"][0] - result["before"]["bbox_min"][0]
    after = result["after"]["bbox_max"][0] - result["after"]["bbox_min"][0]
    assert after - before == pytest.approx(2.0, abs=1e-6)
    assert result["valid_solid"] is True


def test_unknown_parameter_and_oversized_hole_write_nothing(tmp_path, capsys):
    mesh = _mesh(tmp_path, plate_with_hole, "plate")
    script = _script(tmp_path, mesh, "cadquery")
    raw = script.read_bytes()
    out = tmp_path / "no.step"
    assert main(["edit", str(script), "--set", "not_a_param=1", "--in-place", "-o", str(out)]) == 2
    err = capsys.readouterr().err
    assert "unknown parameter" in err and "hole_1_d" in err and "base_x" in err
    assert script.read_bytes() == raw and not out.exists()
    assert main(["edit", str(script), "--set", "hole_1_d=80", "--in-place", "-o", str(out)]) == 2
    err = capsys.readouterr().err
    assert "larger than the wall" in err
    assert script.read_bytes() == raw and not out.exists()


def test_fillet_that_collapses_the_base_is_refused(tmp_path):
    mesh = _mesh(tmp_path, filleted_plate, "fillet")
    script = _script(tmp_path, mesh, "cadquery")
    raw = script.read_bytes()
    with pytest.raises(GeometryError, match="fillet_r"):
        edit_script(str(script), {"fillet_r": 20}, in_place=True)
    assert script.read_bytes() == raw


def test_refuses_a_file_that_mesh2cad_did_not_write(tmp_path, capsys):
    script = tmp_path / "hand.py"
    script.write_text("hole_1_d = 4\n", encoding="utf-8")
    assert main(["params", str(script)]) == 2
    assert "not a script mesh2cad generated" in capsys.readouterr().err
    assert main(["edit", str(script), "--set", "hole_1_d=5", "--in-place"]) == 2
    assert script.read_text(encoding="utf-8") == "hole_1_d = 4\n"


def test_outline_is_listed_and_not_a_number(tmp_path):
    mesh = _mesh(tmp_path, slotted_plate, "slot")
    script = _script(tmp_path, mesh, "cadquery")
    info = list_parameters(str(script))
    outline = next(p for p in info["parameters"] if p["name"] == "outline")
    assert outline["kind"] == "points" and outline["count"] >= 4
    assert outline["unit"] is None
    with pytest.raises(ParamError, match="not a single number"):
        check_edit(str(script), {"outline": 1})


def test_missing_openscad_does_not_pretend_to_render(tmp_path, monkeypatch):
    mesh = _mesh(tmp_path, plate_with_hole, "plate")
    script = _script(tmp_path, mesh, "openscad")
    monkeypatch.setattr("mesh2cad.openscad_code.openscad_bin", lambda: None)
    edited = tmp_path / "edited.scad"
    result = edit_script(str(script), {"hole_1_d": 12}, output_script=str(edited))
    assert result["rendered"] is False
    assert result["valid_solid"] is None
    assert result["output_path"] is None
    assert any("not on PATH" in warning for warning in result["warnings"])
    assert "12.0" in edited.read_text(encoding="utf-8")
    out = tmp_path / "out.stl"
    with pytest.raises(ParamError, match="Nothing was written"):
        edit_script(str(script), {"hole_1_d": 12}, output_path=str(out), output_script=str(tmp_path / "no.scad"))
    assert not out.exists()
    assert not (tmp_path / "no.scad").exists()


def test_missing_freecad_is_reported(tmp_path, monkeypatch):
    mesh = _mesh(tmp_path, cube, "cube")
    script = _script(tmp_path, mesh, "freecad")
    monkeypatch.setattr("mesh2cad.freecad_code.freecad_bin", lambda: None)
    result = edit_script(str(script), {"base_x": "+2"}, output_script=str(tmp_path / "edited.FCMacro"))
    assert result["rendered"] is False
    assert any("freecadcmd" in warning for warning in result["warnings"])
    assert "base_x" in (tmp_path / "edited.FCMacro").read_text(encoding="utf-8")


def test_mcp_list_and_edit(tmp_path):
    pytest.importorskip("mcp")
    import asyncio

    from mesh2cad.mcp_server import edit_params, list_params, mcp
    names = {tool.name for tool in asyncio.run(mcp.list_tools())}
    assert {"list_params", "edit_params"} <= names
    mesh = _mesh(tmp_path, plate_with_hole, "plate")
    script = _script(tmp_path, mesh, "cadquery")
    listed = list_params(str(script))
    assert listed["backend"] == "cadquery"
    hole = next(p for p in listed["parameters"] if p["name"] == "hole_1_d")
    assert hole["value"] == pytest.approx(10.0)
    assert hole["kind"] == "diameter" and hole["unit"] == "mm"
    out = tmp_path / "edited.step"
    result = edit_params(str(script), {"hole_1_d": "+2"}, str(out))
    assert "error" not in result
    assert result["valid_solid"] is True
    assert result["changes"][0]["old"] == pytest.approx(10.0)
    assert result["changes"][0]["new"] == pytest.approx(12.0)
    assert result["delta"]["volume"] < 0
    assert out.is_file()
    assert "12.0" in script.read_text(encoding="utf-8")
    refused = edit_params(str(script), {"nope": 3})
    assert "unknown parameter" in refused["error"]
    assert "hole_1_d" in refused["valid_parameters"]
    untouched = script.read_text(encoding="utf-8")
    refused = edit_params(str(script), {"hole_1_d": 80})
    assert "larger than the wall" in refused["error"]
    assert script.read_text(encoding="utf-8") == untouched
