# Development

From a checkout of the repository:

```
pip install -e .[dev]
pytest
python benchmarks/real_parts.py   # Prusa STLs, needs the network; not in CI
python docs/make_demo.py   # regenerate images (needs matplotlib, pillow)
```

The license is MIT. See [LICENSE](../LICENSE).

Image files written by `docs/make_demo.py` are [before_after.png](before_after.png) and [demo.gif](demo.gif).

Back to the [README](../README.md).
