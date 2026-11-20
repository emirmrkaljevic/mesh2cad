"""Render README before/after images and demo GIF (needs matplotlib, pillow)."""
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from PIL import Image  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
import cadquery as cq  # noqa: E402
from parts import ribbed_plate as demo_part, to_mesh  # noqa: E402

from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.io import add_noise  # noqa: E402
from mesh2cad.report import shape_to_mesh  # noqa: E402
from OCP.TopAbs import TopAbs_FACE  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402

OUT = ROOT / "docs"


def faces_of(shape):
    exp, out = TopExp_Explorer(shape, TopAbs_FACE), []
    while exp.More():
        out.append(shape_to_mesh(exp.Current(), 0.05))
        exp.Next()
    return out


def draw(ax, polys, colors, edges, title, azim):
    ax.clear()
    pc = Poly3DCollection(polys, facecolors=colors, edgecolors=colors if edges == "same" else edges,
                          linewidths=0.15 if edges != "same" else 0.4)
    ax.add_collection3d(pc)
    ax.set_xlim(-20, 20); ax.set_ylim(-20, 20); ax.set_zlim(-10, 14)
    ax.set_box_aspect((1, 1, 0.65))
    ax.view_init(28, azim)
    ax.set_axis_off()
    ax.set_title(title, fontsize=11)


def shade(normals, base):
    light = np.array([0.4, -0.5, 0.75]); light /= np.linalg.norm(light)
    k = 0.45 + 0.55 * np.clip(normals @ light, 0, 1)
    return np.clip(np.asarray(base)[None, :3] * k[:, None], 0, 1)


def main():
    ref = to_mesh(demo_part())
    noisy = add_noise(ref, 0.1, seed=0)
    res = convert_mesh(noisy, tol=0.1, ref=ref)
    r = res.report
    print(r["method"], r["output"]["faces"], r["deviation"])
    before = (noisy.triangles, shade(noisy.face_normals, (0.75, 0.75, 0.8)), (0.2, 0.2, 0.25, 0.6))
    cmap = plt.get_cmap("tab20")
    P, C = [], []
    for i, fm in enumerate(faces_of(res.shape)):
        P += list(fm.triangles)
        C += list(shade(fm.face_normals, cmap(i % 20)))
    after = (np.array(P), np.array(C), "same")
    f = r["output"]["faces"]
    t_before = f"input: {len(noisy.faces)} noisy triangles"
    t_after = f"mesh2cad STEP: {f.get('plane', 0)} planes + {f.get('cylinder', 0)} cylinders"
    fig = plt.figure(figsize=(10, 4.4), dpi=90)
    a1, a2 = fig.add_subplot(121, projection="3d"), fig.add_subplot(122, projection="3d")
    draw(a1, *before, t_before, -60); draw(a2, *after, t_after, -60)
    fig.tight_layout(); fig.savefig(OUT / "before_after.png"); plt.close(fig)
    frames = []
    for az in range(-60, 300, 15):
        fig = plt.figure(figsize=(8, 3.6), dpi=70)
        a1, a2 = fig.add_subplot(121, projection="3d"), fig.add_subplot(122, projection="3d")
        draw(a1, *before, t_before, az); draw(a2, *after, t_after, az)
        fig.tight_layout(); fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3]).convert("P", palette=Image.ADAPTIVE))
        plt.close(fig)
    frames[0].save(OUT / "demo.gif", save_all=True, append_images=frames[1:], duration=120, loop=0, optimize=True)


if __name__ == "__main__":
    main()
