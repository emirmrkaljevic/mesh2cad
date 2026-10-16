"""Parametric CadQuery export: the script runs, matches the part, and edits apply."""
import re

import pytest

pytest.importorskip("cadquery")
from mesh2cad.cadquery_code import Unsupported, check, run_script, to_cadquery  # noqa: E402
from mesh2cad.cli import main  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from parts import (boss, chamfered_boss, countersunk_plate, cube, plate_with_hole,  # noqa: E402
                   ribbed_plate, spherical_dome, stepped_block, tapered_boss, to_mesh)

PARTS = {"cube": cube, "stepped_block": stepped_block, "plate_with_hole": plate_with_hole,
         "boss": boss, "ribbed_plate": ribbed_plate, "countersunk_plate": countersunk_plate,
         "tapered_boss": tapered_boss, "chamfered_boss": chamfered_boss}

EXPECT = {"plate_with_hole": {"hole_1_d": 10.0, "base_t": 6.0},
          "countersunk_plate": {"hole_1_d": 6.0, "hole_1_csk_d": 12.0, "hole_1_csk_angle": 90.0},
          "boss": {"boss_1_d": 16.0, "boss_1_h": 15.0},
          "chamfered_boss": {"boss_1_d": 16.0, "boss_1_chamfer": 2.0},
          "tapered_boss": {"boss_1_d": 16.0, "boss_1_draft_deg": 10.0},
          "ribbed_plate": {"pad_1_w": 14.0, "hole_2_d": 8.0}}


def _params(code):
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _gen(name, sigma=0.0):
    ref = to_mesh(PARTS[name]())
    tol = max(sigma, 0.01)
    res = convert_mesh(add_noise(ref, sigma, seed=0), tol=tol, with_report=False)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    return ref, code


@pytest.mark.parametrize("sigma", [0.0, 0.05])
@pytest.mark.parametrize("name", list(PARTS))
def test_script_reproduces_part(name, sigma):
    ref, code = _gen(name, sigma)
    q = check(code, ref)  # compare against the clean ground truth
    assert q["max_deviation"] < max(3 * sigma, 0.05)
    assert q["volume_error_pct"] < 0.5
    for k, v in EXPECT.get(name, {}).items():
        assert _params(code)[k] == pytest.approx(v, abs=0.5 if ("angle" in k or "deg" in k) else 0.05)


def test_edit_parameter_changes_part():
    ref, code = _gen("plate_with_hole")
    v0 = run_script(code).volume
    bigger = re.sub(r"^hole_1_d\s*=\s*[\d.]+", "hole_1_d = 12.0", code, flags=re.M)
    v1 = run_script(bigger).volume
    expect = 3.141592653589793 * (6.0 ** 2 - 5.0 ** 2) * 6.0
    assert v0 - v1 == pytest.approx(expect, rel=0.01)


def test_unsupported_part():
    res = convert_mesh(to_mesh(spherical_dome()), tol=0.01, with_report=False)
    with pytest.raises(Unsupported):
        to_cadquery(res.mesh, res.regions, res.tol)


def test_cli_cadquery(tmp_path, capsys):
    stl = tmp_path / "p.stl"
    to_mesh(countersunk_plate()).export(stl)
    out = tmp_path / "p.py"
    assert main(["convert", str(stl), "--cadquery", str(out)]) == 0
    assert "hole_1_csk_d" in out.read_text()
    assert not (tmp_path / "p.step").exists()
    to_mesh(spherical_dome()).export(stl)
    assert main(["convert", str(stl), "--cadquery", str(out)]) == 2


def test_mcp_to_cadquery(tmp_path):
    pytest.importorskip("mcp")
    import asyncio

    from mesh2cad.mcp_server import mcp
    assert "to_cadquery" in {t.name for t in asyncio.run(mcp.list_tools())}
    from mesh2cad.mcp_server import to_cadquery as tool
    stl = tmp_path / "b.stl"
    to_mesh(boss()).export(stl)
    r = tool(str(stl))
    assert r["parameters"]["boss_1_d"] == pytest.approx(16.0, abs=0.05)
    assert (tmp_path / "b.py").exists()
