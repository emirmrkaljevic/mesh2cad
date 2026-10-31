"""CAD-exported STLs that were never subdivided.

The synthetic benchmark tessellates a cadquery solid and then splits every
edge to <= 2 mm. A mesh written by a CAD tool is the opposite: a few long
triangles per face, no noise. These six files are printed parts from the
Original Prusa i3 MK3S (Prusa Research, GPL-2.0). They are downloaded from a
pinned commit and cached under benchmarks/.cache, which is gitignored. The
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
from mesh2cad.io import clean, load_mesh  # noqa: E402
from mesh2cad.parametric import Unsupported  # noqa: E402
from mesh2cad.segment import segment  # noqa: E402

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
        return "-"
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


def _script_attempt(emit, mesh, regions, tol):
    try:
        _code, q = emit(mesh, regions, tol)
    except Unsupported as e:
        return {"ok": False, "error": str(e), "dev": None}
    dev = None if q is None else q["max_deviation"]
    return {"ok": True, "error": None, "dev": dev}


def _fit_rank(regions):
    """Prefer the segmentation the STEP path keeps.

    A tighter tolerance can turn a short cylinder into a circumscribed sphere.
    That fit is a worse explanation of the part than the looser one, so its
    refusal (``sphere faces...``) must not hide the real reason.
    """
    penalty = 0
    for r in regions:
        kind = r.surface.kind if r.surface is not None else "faceted"
        if kind == "sphere":
            penalty += 10
        elif kind == "faceted":
            penalty += 5
    return penalty


def _parametric(mesh):
    """Same tolerance retries as ``mesh2cad convert --cadquery/--openscad/--freecad``.

    One segmentation is shared by the exporters. A success is kept. When
    every tolerance fails, the reported reason is the one from the
    segmentation with the fewest spheres, then the fewest unfitted regions,
    which is the same choice ``convert`` makes.
    """
    from mesh2cad.cadquery_code import to_cadquery
    from mesh2cad.freecad_code import to_freecad
    from mesh2cad.openscad_code import to_openscad

    tols = [f * 1e-3 * mesh.scale for f in (1.0, 0.4, 0.2)]
    out = {}
    best_fail = {}
    kinds = (("cadquery", to_cadquery), ("openscad", to_openscad), ("freecad", to_freecad))
    for t in tols:
        regions = segment(mesh, t)
        rank = (_fit_rank(regions), -t)
        for kind, emit in kinds:
            if out.get(kind, {}).get("ok"):
                continue
            attempt = _script_attempt(emit, mesh, regions, t)
            if attempt["ok"]:
                out[kind] = attempt
                continue
            prev = best_fail.get(kind)
            if prev is None or rank < prev[0]:
                best_fail[kind] = (rank, attempt)
        if all(v.get("ok") for v in out.values()) and len(out) == len(kinds):
            break
    for kind, (_rank, attempt) in best_fail.items():
        out.setdefault(kind, attempt)
    return out


def _cell(info):
    if info.get("ok"):
        dev = info["dev"]
        return "yes", (f"{dev:.4f}" if dev is not None else "unverified")
    return "no", "-"


def main():
    print(
        "| part | triangles | method | output faces | regions | valid | "
        "Hausdorff | RMS | vol err % | time s |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|")
    parametric = []
    for name in PARTS:
        mesh = load_mesh(str(_stl(name)))
        t = time.perf_counter()
        res = convert_mesh(mesh)
        dt = time.perf_counter() - t
        r = res.report
        dev = r["deviation"]
        vol = r["volume_error_pct"]
        vol_s = f"{vol:.3f}" if vol is not None else "-"
        print(
            f"| {name} | {r['input']['faces']} | {r['method']} | "
            f"{_faces(r['output']['faces'])} | {_regions(r['regions']['by_type'])} | "
            f"{'yes' if r['output']['valid_solid'] else 'no'} | "
            f"{dev['hausdorff']:.4f} | {dev['rms']:.4f} | {vol_s} | {dt:.1f} |"
        )
        parametric.append((name, _parametric(clean(mesh))))
    print()
    print("| part | cadquery | dev mm | openscad | dev mm | freecad | dev mm |")
    print("|---|---|---|---|---|---|---|")
    reasons = []
    for name, info in parametric:
        cq, cq_dev = _cell(info["cadquery"])
        sc, sc_dev = _cell(info["openscad"])
        fc, fc_dev = _cell(info["freecad"])
        print(f"| {name} | {cq} | {cq_dev} | {sc} | {sc_dev} | {fc} | {fc_dev} |")
        for kind in ("cadquery", "openscad", "freecad"):
            if not info[kind]["ok"]:
                reasons.append(f"{name} {kind}: {info[kind]['error']}")
    for line in reasons:
        print(line)


if __name__ == "__main__":
    main()
