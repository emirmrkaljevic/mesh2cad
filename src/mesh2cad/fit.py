"""Primitive fitting: planes and cylinders."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Plane:
    normal: np.ndarray
    d: float  # normal . x = d
    kind: str = field(default="plane", init=False)

    def distance(self, pts):
        return np.abs(pts @ self.normal - self.d)

    def project(self, pts):
        return pts - np.outer(pts @ self.normal - self.d, self.normal)


@dataclass
class Cylinder:
    axis: np.ndarray
    point: np.ndarray
    radius: float
    kind: str = field(default="cylinder", init=False)

    def _radial(self, pts):
        v = pts - self.point
        return v - np.outer(v @ self.axis, self.axis)

    def distance(self, pts):
        return np.abs(np.linalg.norm(self._radial(pts), axis=1) - self.radius)

    def project(self, pts):
        r = self._radial(pts)
        n = np.linalg.norm(r, axis=1, keepdims=True)
        n[n < 1e-12] = 1e-12
        return pts - r + r / n * self.radius


def fit_plane(pts, weights=None) -> Plane:
    w = np.ones(len(pts)) if weights is None else np.asarray(weights, float)
    c = (pts * w[:, None]).sum(0) / w.sum()
    q = (pts - c) * np.sqrt(w)[:, None]
    _, _, vt = np.linalg.svd(q, full_matrices=False)
    n = vt[-1]
    return Plane(n, float(n @ c))


def fit_cylinder(pts, normals) -> Cylinder | None:
    """Axis = direction most orthogonal to all normals; then circle fit."""
    if len(pts) < 6:
        return None
    _, s, vt = np.linalg.svd(normals, full_matrices=False)
    axis = vt[-1] / np.linalg.norm(vt[-1])
    # normals must span a 2D fan around the axis
    if s[1] < 0.15 * s[0]:
        return None
    u = np.cross(axis, [1.0, 0, 0])
    if np.linalg.norm(u) < 0.1:
        u = np.cross(axis, [0, 1.0, 0])
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    x, y = pts @ u, pts @ v
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]
    sol, *_ = np.linalg.lstsq(A, x * x + y * y, rcond=None)
    cx, cy, k = sol
    r2 = k + cx * cx + cy * cy
    if r2 <= 0:
        return None
    return Cylinder(axis, cx * u + cy * v, float(np.sqrt(r2)))


def refine_cylinder(cyl: Cylinder, pts) -> Cylinder:
    """Nonlinear least-squares refinement of a cylinder on points."""
    from scipy.optimize import least_squares
    a = cyl.axis / np.linalg.norm(cyl.axis)
    p0 = cyl.point - a * (cyl.point @ a)  # closest axis point to origin-ish
    if len(pts) > 4000:
        pts = pts[np.random.default_rng(0).choice(len(pts), 4000, replace=False)]

    def unpack(x):
        ax = x[:3] / np.linalg.norm(x[:3])
        return ax, x[3:6], x[6]

    def res(x):
        ax, p, r = unpack(x)
        v = pts - p
        rad = v - np.outer(v @ ax, ax)
        return np.linalg.norm(rad, axis=1) - r

    try:
        sol = least_squares(res, np.r_[a, p0, cyl.radius], loss="soft_l1",
                            f_scale=max(1e-3, 0.02 * cyl.radius), max_nfev=200)
    except Exception:
        return cyl
    ax, p, r = unpack(sol.x)
    if not np.isfinite(r) or r <= 0:
        return cyl
    return Cylinder(ax, p, float(abs(r)))
