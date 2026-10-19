"""CAD-exported STLs that were never subdivided.

The synthetic benchmark tessellates a cadquery solid and then splits every
edge to <= 2 mm. A mesh written by a CAD tool is the opposite: a few long
triangles per face, no noise. These six files are printed parts from the
Original Prusa i3 MK3S (Prusa Research, GPL-2.0). They are downloaded from a
pinned commit and cached under benchmarks/.cache, which is gitignored — the
STLs are not part of this MIT tree.

    python benchmarks/real_parts.py

CI does not run this (it needs the network).
"""
import pathlib
import sys
import time
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.io import load_mesh  # noqa: E402

PIN = "529e5ad6892ce3b25015b7cca7f0481b8335a76c"
BASE = (
    "https://raw.githubusercontent.com/prusa3d/Original-Prusa-i3/"
    f"{PIN}/Printed-Parts/STL/"
)
PARTS = (
    "plug-aligner",
    "endstop-block",
    "y-belt-idler",
    "y-belt-holder",
    "y-motor-holder",
    "y-rod-holder",
)
CACHE = pathlib.Path(__file__).resolve().parent / ".cache"


def _counts(d, order):
    keys = [k for k in order if d.get(k)]
    extra = [k for k in d if k not in order and d.get(k)]
    if not keys and not extra:
        return "—"
    return ", ".join(f"{d[k]} {k}" for k in keys + extra)


def _faces(f):
    return _counts(f, ("plane", "cylinder", "sphere", "cone", "other"))


def _regions(by_type):
    return _counts(by_type, ("plane", "cylinder", "sphere", "cone", "faceted"))


def _stl(name):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{name}.stl"
    if not path.is_file() or path.stat().st_size == 0:
        urllib.request.urlretrieve(BASE + f"{name}.stl", path)
    return path


def main():
    print(
        "| part | triangles | method | output faces | regions | valid | "
        "Hausdorff | RMS | vol err % | time s |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|")
    for name in PARTS:
        mesh = load_mesh(str(_stl(name)))
        t = time.perf_counter()
        res = convert_mesh(mesh)
        dt = time.perf_counter() - t
        r = res.report
        dev = r["deviation"]
        vol = r["volume_error_pct"]
        vol_s = f"{vol:.3f}" if vol is not None else "—"
        print(
            f"| {name} | {r['input']['faces']} | {r['method']} | "
            f"{_faces(r['output']['faces'])} | {_regions(r['regions']['by_type'])} | "
            f"{'yes' if r['output']['valid_solid'] else 'no'} | "
            f"{dev['hausdorff']:.4f} | {dev['rms']:.4f} | {vol_s} | {dt:.1f} |"
        )


if __name__ == "__main__":
    main()
