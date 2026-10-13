"""Synthetic benchmark: cadquery parts -> mesh -> noise -> mesh2cad -> compare to truth.

Usage: python benchmarks/run.py [--seeds 3]
"""
import argparse
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tests"))
from parts import CURVED, PLANAR, to_mesh  # noqa: E402

from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, default=3)
a = ap.parse_args()
print("| part | noise σ (mm) | exact | faces (plane/cyl) | Hausdorff (mm) | RMS (mm) | vol err % | time s |")
print("|---|---|---|---|---|---|---|---|")
for name, fn in {**PLANAR, **CURVED}.items():
    ref = to_mesh(fn())
    for sigma in (0.0, 0.05, 0.2):
        ok, H, R, V, T, F = 0, [], [], [], [], set()
        for seed in range(a.seeds):
            t = time.time()
            r = convert_mesh(add_noise(ref, sigma, seed=seed), tol=max(sigma, 0.01), ref=ref).report
            T.append(time.time() - t)
            ok += r["method"] == "exact" and r["output"]["valid_solid"]
            f = r["output"]["faces"]
            F.add(f"{f.get('plane', 0)}/{f.get('cylinder', 0)}" + ("+" if f.get("other") else ""))
            H.append(r["deviation"]["hausdorff"]); R.append(r["deviation"]["rms"]); V.append(r["volume_error_pct"])
        print(f"| {name} | {sigma} | {ok}/{a.seeds} | {', '.join(sorted(F))} | {max(H):.3f} | {max(R):.3f} | {max(V):.2f} | {max(T):.1f} |")
