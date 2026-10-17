"""Primitive fitting: planes, cylinders, spheres, and cones."""
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


@dataclass
class Sphere:
    center: np.ndarray
    radius: float
    kind: str = field(default="sphere", init=False)

    def distance(self, pts):
        return np.abs(np.linalg.norm(pts - self.center, axis=1) - self.radius)

    def project(self, pts):
        v = pts - self.center
        n = np.linalg.norm(v, axis=1, keepdims=True)
        n = np.maximum(n, 1e-12)
        return self.center + v / n * self.radius


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


def _unit_normal(plane: Plane):
    n = np.asarray(plane.normal, float)
    ln = float(np.linalg.norm(n))
    if ln < 1e-12:
        return None, None
    return n / ln, float(plane.d) / ln


def _tangent_radius(c0, c1, axis, pts, radius0, lo, hi):
    """Radius of the cylinder ``c0 + R c1`` that best fits ``pts``.

    A free circle fit is pulled off when the patch also contains wall
    triangles. Those points lie on the tangent planes, so a median residual
    against the tangent-constrained axis stays on the arc.
    """
    if pts is None or len(pts) < 8 or not np.isfinite(radius0) or hi <= lo:
        return radius0
    sample = np.asarray(pts, float)
    if len(sample) > 2000:
        sample = sample[::int(np.ceil(len(sample) / 2000.0))]
    axis = np.asarray(axis, float)
    axis = axis / max(float(np.linalg.norm(axis)), 1e-12)

    def mad(radius):
        centre = c0 + radius * c1
        rel = sample - centre
        rel = rel - np.outer(rel @ axis, axis)
        return float(np.median(np.abs(np.linalg.norm(rel, axis=1) - radius)))

    grid = np.linspace(lo, hi, 33)
    errs = np.empty(len(grid))
    for i, radius in enumerate(grid):
        errs[i] = mad(float(radius))
    best = int(np.argmin(errs))
    step = float(grid[1] - grid[0])
    grid = np.linspace(max(lo, float(grid[best]) - step), min(hi, float(grid[best]) + step), 17)
    errs = np.empty(len(grid))
    for i, radius in enumerate(grid):
        errs[i] = mad(float(radius))
    best = int(np.argmin(errs))
    # a flat residual means the points do not pick out a radius; keep the hint
    if errs[best] > mad(radius0) * 0.99 + 1e-6:
        return radius0
    return float(grid[best])


def snap_tangent_cylinder(cyl: Cylinder, planes, tol: float, pts=None) -> Cylinder | None:
    """Cylinder of ``cyl`` made exactly tangent to two neighbouring planes.

    A constant-radius fillet between two planes is a partial cylinder whose
    generators lie in those planes. The fitted axis is only a hint: the
    returned axis is the cross product of the two tangent normals, and the
    radius is the distance from that axis to each plane. When ``pts`` is
    given, that radius is the one whose tangent cylinder actually passes
    through the arc, not the free fit. Returns None when no pair of
    ``planes`` is tangent within ``tol``.
    """
    if cyl.radius <= 0 or not planes:
        return None
    axis0 = np.asarray(cyl.axis, float)
    axis0 = axis0 / max(float(np.linalg.norm(axis0)), 1e-12)
    point = np.asarray(cyl.point, float)
    limit = max(8.0 * tol, 0.15 * cyl.radius)
    cands = []
    for plane in planes:
        if not isinstance(plane, Plane):
            continue
        n, d = _unit_normal(plane)
        if n is None or abs(float(n @ axis0)) > 0.30:
            continue
        signed = float(n @ point - d)
        err = abs(abs(signed) - cyl.radius)
        if err <= limit:
            cands.append((err, n, d, signed))
    if len(cands) < 2:
        return None
    cands.sort(key=lambda t: t[0])
    pair = None
    for i in range(len(cands)):
        for j in range(i + 1, len(cands)):
            if np.linalg.norm(np.cross(cands[i][1], cands[j][1])) >= 0.30:
                pair = (cands[i], cands[j])
                break
        if pair is not None:
            break
    if pair is None:
        return None
    # a third plane equally tangent means the strip is not one edge fillet
    third = max(4.0 * tol, 0.08 * cyl.radius)
    for cand in cands:
        if cand is pair[0] or cand is pair[1] or cand[0] > third:
            continue
        extra = cand[1]
        if min(np.linalg.norm(np.cross(extra, pair[0][1])),
               np.linalg.norm(np.cross(extra, pair[1][1]))) >= 0.30:
            return None
    (_, n1, d1, s1), (_, n2, d2, s2) = pair
    axis = np.cross(n1, n2)
    axis = axis / np.linalg.norm(axis)
    if float(axis @ axis0) < 0.0:
        axis = -axis
    radius0 = float(np.median([cyl.radius, abs(s1), abs(s2)]))
    if not np.isfinite(radius0) or radius0 <= 0:
        return None
    # exactly r on the same side of each plane as the fitted centre.
    # c(R) = c0 + R c1 is the axis at distance R from both planes.
    sigma1 = np.sign(s1) if s1 != 0 else 1.0
    sigma2 = np.sign(s2) if s2 != 0 else 1.0
    mat = np.stack([n1, n2, axis])
    try:
        c0 = np.linalg.solve(mat, np.array([d1, d2, float(axis @ point)]))
        c1 = np.linalg.solve(mat, np.array([sigma1, sigma2, 0.0]))
    except np.linalg.LinAlgError:
        return None
    if not np.isfinite(c0).all() or not np.isfinite(c1).all():
        return None
    half = max(6.0 * tol, 0.2 * cyl.radius)
    radius = _tangent_radius(c0, c1, axis, pts, radius0, max(1e-6, radius0 - half), radius0 + half)
    if not np.isfinite(radius) or radius <= 0:
        return None
    origin = c0 + radius * c1
    if not np.isfinite(origin).all():
        return None
    shift = float(np.linalg.norm(np.cross(origin - point, axis)))
    if shift > max(6.0 * tol, 0.25 * cyl.radius) or abs(radius - cyl.radius) > half:
        return None
    return Cylinder(axis, origin, radius)


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


def _sphere_from_normals(pts, normals):
    """Center = least-squares intersection of normal rays; radius = median distance."""
    n = np.asarray(normals, float)
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.maximum(ln, 1e-12)
    P = np.eye(3) - n[:, :, None] * n[:, None, :]
    A = P.sum(axis=0)
    ev = np.linalg.eigvalsh(A)
    if ev[0] <= 1e-8 * max(ev[-1], 1e-12):
        return None
    b = np.einsum("mij,mj->i", P, pts)
    try:
        center = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return None
    if not np.isfinite(center).all():
        return None
    radius = float(np.median(np.linalg.norm(pts - center, axis=1)))
    if not np.isfinite(radius) or radius <= 0:
        return None
    return center, radius


def _sphere_algebraic(pts):
    """Linearized (x-c)^2 = r^2 fit. Unstable on tiny caps; caller compares candidates."""
    A = np.c_[pts, np.ones(len(pts))]
    rhs = -np.sum(pts * pts, axis=1)
    sol, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    center = -0.5 * sol[:3]
    r2 = float(center @ center - sol[3])
    if not np.isfinite(r2) or r2 <= 0:
        return None
    return center, float(np.sqrt(r2))


def fit_sphere(pts, normals=None) -> Sphere | None:
    """Fit a sphere from points, using normal rays when normals are available.

    Both the normal-ray center and the algebraic fit are tried; the candidate
    with the smaller median residual is kept. Call refine_sphere afterwards.
    """
    pts = np.asarray(pts, float)
    if len(pts) < 6:
        return None
    cands = []
    if normals is not None:
        hit = _sphere_from_normals(pts, normals)
        if hit is not None:
            cands.append(hit)
    hit = _sphere_algebraic(pts)
    if hit is not None:
        cands.append(hit)
    ext = float(np.linalg.norm(pts - pts.mean(0), axis=1).max())
    best = None
    for center, radius in cands:
        if ext > 0 and radius > max(1e6, 1e4 * ext):
            continue
        err = float(np.median(np.abs(np.linalg.norm(pts - center, axis=1) - radius)))
        if best is None or err < best[0]:
            best = (err, center, radius)
    if best is None:
        return None
    return Sphere(np.asarray(best[1], float), float(best[2]))


def refine_sphere(sph: Sphere, pts) -> Sphere:
    """Nonlinear least-squares refinement of a sphere on points."""
    from scipy.optimize import least_squares
    pts = np.asarray(pts, float)
    if len(pts) < 4:
        return sph
    if len(pts) > 4000:
        pts = pts[np.random.default_rng(0).choice(len(pts), 4000, replace=False)]

    def res(x):
        return np.linalg.norm(pts - x[:3], axis=1) - x[3]

    try:
        sol = least_squares(res, np.r_[sph.center, sph.radius], loss="soft_l1",
                            f_scale=max(1e-3, 0.02 * sph.radius), max_nfev=200)
    except Exception:
        return sph
    center, radius = sol.x[:3], float(abs(sol.x[3]))
    if not np.isfinite(center).all() or not np.isfinite(radius) or radius <= 0:
        return sph
    return Sphere(center, radius)


@dataclass
class Cone:
    """Right circular cone. ``axis`` points from the apex into the cone;
    ``half_angle`` is in radians, 0 < half_angle < pi/2."""
    apex: np.ndarray
    axis: np.ndarray
    half_angle: float
    kind: str = field(default="cone", init=False)

    def _hr(self, pts):
        v = pts - self.apex
        h = v @ self.axis
        rad = v - np.outer(h, self.axis)
        return h, rad, np.linalg.norm(rad, axis=1)

    def radius_at(self, h):
        return h * np.tan(self.half_angle)

    def distance(self, pts):
        h, _, rho = self._hr(pts)
        c, s = np.cos(self.half_angle), np.sin(self.half_angle)
        return np.abs(rho * c - h * s)

    def normal(self, pts):
        """Unit outward normal (away from the axis)."""
        _, rad, rho = self._hr(pts)
        u = rad / np.maximum(rho, 1e-12)[:, None]
        return u * np.cos(self.half_angle) - self.axis * np.sin(self.half_angle)

    def project(self, pts):
        n = self.normal(pts)
        h, _, rho = self._hr(pts)
        sd = rho * np.cos(self.half_angle) - h * np.sin(self.half_angle)
        return pts - n * sd[:, None]


def fit_cone(pts, normals) -> Cone | None:
    """Linear cone fit from oriented normals.

    Every normal of a cone makes the same angle with the axis, so the axis is
    the direction minimising the variance of n.a; every tangent plane passes
    through the apex, so the apex solves n_i . x = n_i . p_i in least squares.
    """
    pts = np.asarray(pts, float)
    n = np.asarray(normals, float)
    if len(pts) < 8:
        return None
    n = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    w, V = np.linalg.eigh(np.cov(n.T))
    axis = V[:, 0]
    if w[0] > 0.35 * max(w[1], 1e-12):  # noisy normals: loose, residual checks decide
        return None
    c = float(np.mean(n @ axis))
    if abs(c) < 0.03 or abs(c) > 0.97:  # cylinder-like or plane-like
        return None
    A = n.T @ n
    if np.linalg.eigvalsh(A)[0] < 1e-6 * len(n):
        return None
    apex = np.linalg.solve(A, n.T @ np.einsum("ij,ij->i", n, pts))
    # outward normal has n.a = -sin(alpha) when a points into the cone
    if np.median((pts - apex) @ axis) < 0:
        axis = -axis
    c = float(np.mean(n @ axis))
    alpha = float(np.arcsin(np.clip(abs(c), 0.0, 1.0)))
    if not np.isfinite(apex).all() or not 0.02 < alpha < np.pi / 2 - 0.02:
        return None
    return Cone(apex, axis, alpha)


def refine_cone(cone: Cone, pts) -> Cone:
    from scipy.optimize import least_squares
    pts = np.asarray(pts, float)
    if len(pts) < 6:
        return cone
    if len(pts) > 4000:
        pts = pts[np.random.default_rng(0).choice(len(pts), 4000, replace=False)]
    scale = float(np.linalg.norm(pts - pts.mean(0), axis=1).max()) or 1.0

    def unpack(x):
        return x[:3], x[3:6] / np.linalg.norm(x[3:6]), x[6]

    def res(x):
        apex, ax, al = unpack(x)
        v = pts - apex
        h = v @ ax
        rho = np.linalg.norm(v - np.outer(h, ax), axis=1)
        return rho * np.cos(al) - h * np.sin(al)

    try:
        sol = least_squares(res, np.r_[cone.apex, cone.axis, cone.half_angle], loss="soft_l1",
                            f_scale=max(1e-3, 0.01 * scale), max_nfev=300)
    except Exception:
        return cone
    apex, ax, al = unpack(sol.x)
    if not np.isfinite(sol.x).all() or not 0.01 < al < np.pi / 2 - 0.01:
        return cone
    return Cone(apex, ax, float(al))


def fit_cone_about_axis(origin, axis, pts, normals=None) -> Cone | None:
    """Right circular cone locked to ``axis`` (a line through ``origin``).

    Narrow chamfers and countersinks do not span enough of the normal circle
    for a free axis estimate. The neighbouring cylinder or cap plane supplies
    the axis; half-angle comes from ``|n·axis| = sin(alpha)`` when normals
    agree, otherwise from ``rho = tan(alpha) * h``.
    """
    pts = np.asarray(pts, float)
    origin = np.asarray(origin, float)
    axis = np.asarray(axis, float)
    nrm = float(np.linalg.norm(axis))
    if nrm < 1e-12 or len(pts) < 6:
        return None
    axis = axis / nrm

    def frame(ax):
        v = pts - origin
        h = v @ ax
        rho = np.linalg.norm(v - np.outer(h, ax), axis=1)
        return h, rho

    alpha = None
    if normals is not None:
        n = np.asarray(normals, float)
        ln = np.linalg.norm(n, axis=1)
        ok = ln > 1e-12
        if int(ok.sum()) >= 6:
            n = n[ok] / ln[ok, None]
            s = np.abs(n @ axis)
            band = (s > 0.12) & (s < 0.96)
            if int(band.sum()) >= 6:
                alpha = float(np.arcsin(np.clip(np.median(s[band]), 0.0, 1.0)))
    h, rho = frame(axis)
    # normals give |sin(alpha)| but not which way the cone opens. Radius has
    # to increase along the axis; a negative slope means the axis is backwards.
    if len(h) >= 6 and float(np.ptp(h)) > 1e-4:
        sol, *_ = np.linalg.lstsq(np.c_[h, np.ones_like(h)], rho, rcond=None)
        if float(sol[0]) < -0.05:
            axis = -axis
            h, rho = frame(axis)
    if alpha is None:
        if len(h) < 6:
            return None
        sol, *_ = np.linalg.lstsq(np.c_[h, np.ones_like(h)], rho, rcond=None)
        slope = float(sol[0])
        if slope <= 0.05:
            return None
        alpha = float(np.arctan(slope))
    if not np.isfinite(alpha) or not 0.08 < alpha < np.pi / 2 - 0.08:
        return None
    tana = float(np.tan(alpha))
    h_i = h - rho / tana
    h_apex = float(np.median(h_i))
    if not np.isfinite(h_apex):
        return None
    resid = np.abs(rho - (h - h_apex) * tana)
    scale = float(np.median(resid))
    inl = resid <= max(3.0 * scale, 1e-9)
    if int(inl.sum()) >= 6:
        h_apex = float(np.median(h_i[inl]))
    apex = origin + axis * h_apex
    along = (pts - apex) @ axis
    if float(np.median(along)) < 0:
        axis = -axis
        along = -along
    # a handful of boundary vertices may cross the apex under noise
    if float(np.percentile(along, 8)) <= 0:
        return None
    return Cone(np.asarray(apex, float), np.asarray(axis, float), float(alpha))
