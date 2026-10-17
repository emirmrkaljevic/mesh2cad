"""Synthetic benchmark: cadquery parts -> mesh -> noise -> mesh2cad -> compare to truth.

Usage: python benchmarks/run.py [--seeds 3]
"""
import argparse
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tests"))
from parts import CONES, CURVED, FILLET_RADIUS, FILLETS, PLANAR, SPHERES, to_mesh  # noqa: E402

from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, default=3)
a = ap.parse_args()
print("| part | noise σ (mm) | exact | faces (plane/cyl[/sph[/cone]]) | radius err (mm) | Hausdorff (mm) | RMS (mm) | vol err % | time s |")
print("|---|---|---|---|---|---|---|---|---|")


def _faces(f):
    parts = [str(f.get("plane", 0)), str(f.get("cylinder", 0))]
    if f.get("sphere") or f.get("cone") or f.get("other"):
        parts.append(str(f.get("sphere", 0)))
    if f.get("cone"):
        parts.append(str(f.get("cone", 0)))
    s = "/".join(parts)
    return s + ("+" if f.get("other") else "")


def _radius_err(res, truth):
    if truth is None or res.report["method"] != "exact":
        return None
    radii = [r.surface.radius for r in res.regions if getattr(r.surface, "kind", None) == "cylinder"]
    if not radii:
        return None
    return max(abs(rad - truth) for rad in radii)


parts = [(name, fn, None) for name, fn in {**PLANAR, **CURVED, **SPHERES, **CONES}.items()]
parts += [(name, fn, FILLET_RADIUS[name]) for name, fn in FILLETS.items()]
for name, fn, truth in parts:
    ref = to_mesh(fn())
    for sigma in (0.0, 0.05, 0.2):
        ok, H, R, V, T, F, RE = 0, [], [], [], [], set(), []
        for seed in range(a.seeds):
            t = time.time()
            res = convert_mesh(add_noise(ref, sigma, seed=seed), tol=max(sigma, 0.01), ref=ref)
            r = res.report
            T.append(time.time() - t)
            ok += r["method"] == "exact" and r["output"]["valid_solid"]
            F.add(_faces(r["output"]["faces"]))
            H.append(r["deviation"]["hausdorff"]); R.append(r["deviation"]["rms"]); V.append(r["volume_error_pct"])
            err = _radius_err(res, truth)
            if err is not None:
                RE.append(err)
        re_s = f"{max(RE):.3f}" if RE else "—"
        print(f"| {name} | {sigma} | {ok}/{a.seeds} | {', '.join(sorted(F))} | {re_s} | {max(H):.3f} | {max(R):.3f} | {max(V):.2f} | {max(T):.1f} |")
