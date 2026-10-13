"""Mesh loading and preprocessing."""
from __future__ import annotations

import numpy as np
import trimesh


def load_mesh(path: str) -> trimesh.Trimesh:
    """Load STL/OBJ/PLY into a single cleaned Trimesh."""
    m = trimesh.load(path, force="mesh", process=True)
    if isinstance(m, trimesh.Scene):
        m = trimesh.util.concatenate(list(m.geometry.values()))
    return clean(m)


def clean(m: trimesh.Trimesh) -> trimesh.Trimesh:
    m = m.copy()
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces())
    m.update_faces(m.unique_faces())
    m.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(m)
    return m


def add_noise(m: trimesh.Trimesh, sigma: float, seed: int = 0) -> trimesh.Trimesh:
    """Gaussian vertex noise (for self-tests / benchmarks)."""
    rng = np.random.default_rng(seed)
    out = m.copy()
    out.vertices = out.vertices + rng.normal(0.0, sigma, out.vertices.shape)
    return out
