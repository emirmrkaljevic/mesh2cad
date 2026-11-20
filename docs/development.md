# Development

Part of [mesh2cad](../README.md). How to run the tests and rebuild the README images from a checkout of the repository.

```
pip install -e .[dev]
pytest
python benchmarks/real_parts.py   # Prusa STLs, needs the network; not in CI
python scripts/make_demo.py   # regenerate images (needs matplotlib, pillow)
```

`scripts/make_demo.py` writes [before_after.png](before_after.png) and [demo.gif](demo.gif) into this folder.

The project site is built from `site/` plus these pages by the `pages` workflow.

The license is MIT. See [LICENSE](../LICENSE).
