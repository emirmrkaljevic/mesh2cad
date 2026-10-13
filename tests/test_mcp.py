import asyncio
import json

import pytest

pytest.importorskip("mcp")
from mesh2cad.mcp_server import mcp  # noqa: E402
from parts import cube, to_mesh  # noqa: E402


def test_mcp_tools(tmp_path):
    stl = tmp_path / "c.stl"
    to_mesh(cube()).export(stl)
    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert {"convert", "analyze", "report"} <= names
    res = asyncio.run(mcp.call_tool("convert", {"input_path": str(stl)}))
    if isinstance(res, tuple):
        res = res[0]
    content = getattr(res, "content", res)
    data = getattr(res, "structured_content", None) or getattr(res, "structuredContent", None)
    if not data or "output" not in data:
        data = json.loads(content[0].text)
    assert data["output"]["faces"] == {"plane": 6}
    assert (tmp_path / "c.step").exists()
