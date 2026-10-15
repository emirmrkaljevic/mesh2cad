"""Synthetic benchmark: cadquery parts -> mesh -> noise -> mesh2cad -> compare to truth.

Usage: python benchmarks/run.py [--seeds 3]
"""
import argparse
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tests"))
from parts import CONES, CURVED, PLANAR, SPHERES, to_mesh  # noqa: E402

from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, default=3)
a = ap.parse_args()
print("| part | noise σ (mm) | exact | faces (plane/cyl[/sph[/cone]]) | Hausdorff (mm) | RMS (mm) | vol err % | time s |")
print("|---|---|---|---|---|---|---|---|")


def _faces(f):
    parts = [str(f.get("plane", 0)), str(f.get("cylinder", 0))]
    if f.get("sphere") or f.get("cone") or f.get("other"):
        parts.append(str(f.get("sphere", 0)))
    if f.get("cone"):
        parts.append(str(f.get("cone", 0)))
    s = "/".join(parts)
    return s + ("+" if f.get("other") else "")


for name, fn in {**PLANAR, **CURVED, **SPHERES, **CONES}.items():
    ref = to_mesh(fn())
    for sigma in (0.0, 0.05, 0.2):
        ok, H, R, V, T, F = 0, [], [], [], [], set()
        for seed in range(a.seeds):
            t = time.time()
            r = convert_mesh(add_noise(ref, sigma, seed=seed), tol=max(sigma, 0.01), ref=ref).report
            T.append(time.time() - t)
            ok += r["method"] == "exact" and r["output"]["valid_solid"]
            F.add(_faces(r["output"]["faces"]))
            H.append(r["deviation"]["hausdorff"]); R.append(r["deviation"]["rms"]); V.append(r["volume_error_pct"])
        print(f"| {name} | {sigma} | {ok}/{a.seeds} | {', '.join(sorted(F))} | {max(H):.3f} | {max(R):.3f} | {max(V):.2f} | {max(T):.1f} |")
