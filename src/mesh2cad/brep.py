"""Build a B-rep solid from segmented mesh and export STEP.

Strategy: snap every mesh vertex onto the fitted surfaces of all regions it
touches (corner -> plane intersection point, edge -> intersection curve),
rebuild triangles as OCC faces, sew, then merge co-surface faces with
ShapeUpgrade_UnifySameDomain so each fitted plane becomes one clean face.
Unfitted regions stay faceted (hybrid fallback).

Exact path (build_planar): plane polygons from triple-plane intersections,
circles from plane ∩ cylinder and plane ∩ sphere, full cylindrical faces,
partial fillet cylinders tangent to two planes, and spherical caps/zones.
"""
from __future__ import annotations

import numpy as np
import trimesh
from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon,
                                BRepBuilderAPI_MakeSolid, BRepBuilderAPI_Sewing)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepLib import BRepLib
from OCP.gp import gp_Pnt
from OCP.IFSelect import IFSelect_RetDone
from OCP.ShapeFix import ShapeFix_Shape, ShapeFix_Solid
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer
from OCP.TopAbs import TopAbs_SHELL
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Shape

from .fit import Cone, Cylinder, Plane, Sphere
from .segment import Region, _analytic_mesh


def _curved_normal(surface, pts):
    """Unit surface normal used to linearise a cylinder or sphere constraint."""
    if isinstance(surface, Cone):
        return surface.normal(pts)
    if isinstance(surface, Cylinder):
        rad = surface._radial(pts)
    elif isinstance(surface, Sphere):
        rad = pts - surface.center
    else:
        raise TypeError(type(surface))
    n = np.linalg.norm(rad, axis=1, keepdims=True)
    n = np.maximum(n, 1e-12)
    return rad / n


def snap_vertices(m: trimesh.Trimesh, regions: list[Region], tol: float, iters: int = 10):
    """Move vertices onto the fitted surfaces of every region they touch.
    Plane-only vertices are solved in one batched least-squares step;
    vertices touching a cylinder or sphere use a few linearised iterations."""
    nv = len(m.vertices)
    x0 = m.vertices.copy()
    lam = 1e-6
    rows_v, rows_r = [], []
    for i, r in enumerate(regions):
        if r.surface is None:
            continue
        vs = np.unique(m.faces[r.faces])
        rows_v.append(vs)
        rows_r.append(np.full(len(vs), i))
    if not rows_v:
        return x0
    rv, rr = np.concatenate(rows_v), np.concatenate(rows_r)
    is_plane = np.array([isinstance(r.surface, Plane) for r in regions])
    curved = np.zeros(nv, bool)
    curved[rv[~is_plane[rr]]] = True
    touched = np.zeros(nv, bool)
    touched[rv] = True

    def solve(p, N, D, idx):
        # accumulate N^T N and N^T D per vertex
        A = np.zeros((nv, 3, 3))
        b = np.zeros((nv, 3))
        np.add.at(A, idx, N[:, :, None] * N[:, None, :])
        np.add.at(b, idx, N * D[:, None])
        A += lam * np.eye(3)
        b += lam * p
        return np.linalg.solve(A, b[..., None])[..., 0]

    pl = is_plane[rr]
    Np = np.array([regions[i].surface.normal if is_plane[i] else np.zeros(3) for i in range(len(regions))])
    Dp = np.array([regions[i].surface.d if is_plane[i] else 0.0 for i in range(len(regions))])
    x = x0.copy()
    for _ in range(iters if curved.any() else 1):
        N = np.zeros((len(rv), 3))
        D = np.zeros(len(rv))
        N[pl], D[pl] = Np[rr[pl]], Dp[rr[pl]]
        for i in np.unique(rr[~pl]):
            sel = (rr == i)
            c = regions[i].surface
            p = x[rv[sel]]
            q = c.project(p)
            n = _curved_normal(c, p)
            N[sel], D[sel] = n, np.einsum("ij,ij->i", n, q)
        xn = solve(x0, N, D, rv)
        xn[~touched] = x0[~touched]
        # guard against ill-conditioned (near-parallel) solves
        bad = np.linalg.norm(xn - x0, axis=1) > 4 * tol + 1e-9
        if bad.any():
            first = {}
            for v, i in zip(rv, rr):
                first.setdefault(v, i)
            for v in np.where(bad)[0]:
                xn[v] = regions[first[v]].surface.project(x0[v][None])[0]
        done = np.abs(xn - x).max() < 1e-10
        x = xn
        if done:
            break
    return x


def build_solid(m: trimesh.Trimesh, verts: np.ndarray, sew_tol: float) -> TopoDS_Shape:
    sew = BRepBuilderAPI_Sewing(sew_tol)
    pts = [gp_Pnt(*map(float, v)) for v in verts]
    for f in m.faces:
        a, b, c = f
        if len({a, b, c}) < 3:
            continue
        poly = BRepBuilderAPI_MakePolygon(pts[a], pts[b], pts[c], True)
        if not poly.IsDone():
            continue
        mf = BRepBuilderAPI_MakeFace(poly.Wire(), True)
        if mf.IsDone():
            sew.Add(mf.Face())
    sew.Perform()
    shape = sew.SewedShape()
    exp = TopExp_Explorer(shape, TopAbs_SHELL)
    if exp.More():
        ms = BRepBuilderAPI_MakeSolid()
        while exp.More():
            ms.Add(TopoDS.Shell_s(exp.Current()))
            exp.Next()
        fix = ShapeFix_Solid(ms.Solid())
        fix.Perform()
        shape = fix.Solid()
    # sewing leaves a few edges without SameParameter; BRepCheck then
    # rejects an otherwise closed solid
    same_parameter(shape, tol=max(min(sew_tol, 1e-4), 1e-7))
    return shape


def same_parameter(shape: TopoDS_Shape, tol: float = 1e-5) -> TopoDS_Shape:
    """Force edge 3D curves and pcurves onto one parameter range, in place."""
    BRepLib.SameParameter_s(shape, float(tol), True)
    return shape


def fix_shape(shape: TopoDS_Shape, prec: float = 1e-6, max_tol: float = 1e-3) -> TopoDS_Shape:
    """Run ShapeFix. Used when unify flips a handful of faces inside-out."""
    fx = ShapeFix_Shape(shape)
    fx.SetPrecision(float(prec))
    fx.SetMaxTolerance(float(max_tol))
    fx.Perform()
    return fx.Shape()


def unify(shape: TopoDS_Shape, lin_tol: float, ang_tol: float = 1e-3) -> TopoDS_Shape:
    u = ShapeUpgrade_UnifySameDomain(shape, True, True, True)
    u.SetLinearTolerance(lin_tol)
    u.SetAngularTolerance(ang_tol)
    u.Build()
    return u.Shape()


_UNIFY_CHILD = """
import sys
from OCP.BRep import BRep_Builder
from OCP.BRepTools import BRepTools
from OCP.TopoDS import TopoDS_Shape
from mesh2cad.brep import unify
s = TopoDS_Shape()
BRepTools.Read_s(s, sys.argv[1], BRep_Builder())
BRepTools.Write_s(unify(s, float(sys.argv[3]), float(sys.argv[4])), sys.argv[2])
"""


def safe_unify(shape: TopoDS_Shape, lin_tol: float, ang_tol: float = 1e-3,
               timeout: float = 120.0, keep_invalid: bool = False) -> TopoDS_Shape:
    """UnifySameDomain in a child process. OCCT can segfault inside it on some
    snapped shells (valid per BRepCheck); then the shape is returned un-merged
    instead of taking the whole process down.

    With ``keep_invalid``, a unify that finishes but fails BRepCheck is returned
    anyway so the caller can run SameParameter or ShapeFix. A crash, timeout, or
    null result still returns ``shape`` itself.
    """
    import os
    import subprocess
    import sys
    import tempfile

    from OCP.BRepTools import BRepTools
    with tempfile.TemporaryDirectory() as d:
        src, dst = os.path.join(d, "in.brep"), os.path.join(d, "out.brep")
        BRepTools.Write_s(shape, src)
        try:
            p = subprocess.run([sys.executable, "-c", _UNIFY_CHILD, src, dst, repr(lin_tol), repr(ang_tol)],
                               capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return shape
        if p.returncode != 0 or not os.path.exists(dst):
            return shape
        out = TopoDS_Shape()
        BRepTools.Read_s(out, dst, BRep_Builder())
        if out.IsNull():
            return shape
        if keep_invalid or is_valid(out):
            return out
        return shape


def is_valid(shape: TopoDS_Shape) -> bool:
    return bool(BRepCheck_Analyzer(shape).IsValid())


def write_step(shape: TopoDS_Shape, path: str) -> None:
    w = STEPControl_Writer()
    w.Transfer(shape, STEPControl_AsIs)
    if w.Write(str(path)) != IFSelect_RetDone:
        raise RuntimeError(f"STEP write failed: {path}")


# ---------------------------------------------------------------- planar B-rep
def _plane_basis(n):
    u = np.cross(n, [1.0, 0, 0])
    if np.linalg.norm(u) < 0.1:
        u = np.cross(n, [0, 1.0, 0])
    u /= np.linalg.norm(u)
    return u, np.cross(n, u)


def _intersect3(planes):
    N = np.array([p.normal for p in planes])
    D = np.array([p.d for p in planes])
    if abs(np.linalg.det(N)) < 1e-6:
        return None
    return np.linalg.solve(N, D)


def _loops(m, label, rid):
    """Directed boundary edges of region rid chained into loops of (a, b, nbr)."""
    faces = np.where(label == rid)[0]
    adj = {}
    for (fa, fb), (ea, eb) in zip(m.face_adjacency, m.face_adjacency_edges):
        adj[(min(ea, eb), max(ea, eb))] = (fa, fb)
    nxt = {}
    for f in faces:
        tri = m.faces[f]
        for k in range(3):
            a, b = int(tri[k]), int(tri[(k + 1) % 3])
            fa, fb = adj.get((min(a, b), max(a, b)), (f, f))
            other = fb if fa == f else fa
            if label[other] != rid:
                if a in nxt:
                    return None  # non-manifold boundary vertex
                nxt[a] = (b, int(label[other]))
    loops = []
    while nxt:
        start = next(iter(nxt))
        loop, a = [], start
        while True:
            if a not in nxt:
                return None
            b, s = nxt.pop(a)
            loop.append((a, b, s))
            a = b
            if a == start:
                break
        loops.append(loop)
    return loops


def _collapse_runs(runs):
    """Merge adjacent runs of the same region, including across the loop seam."""
    if runs and runs[0][0] == runs[-1][0] and len(runs) > 1:
        runs[0][1] += runs.pop()[1]
    j = 0
    while j < len(runs) and len(runs) > 1:
        k = (j + 1) % len(runs)
        if j != k and runs[j][0] == runs[k][0]:
            runs[j][1] += runs[k][1]
            runs.pop(k)
            if runs and runs[0][0] == runs[-1][0] and len(runs) > 1:
                runs[0][1] += runs.pop()[1]
            j = 0
        else:
            j += 1
    return runs


def _transition_vertex(loop, prev_sid, this_sid):
    """Mesh vertex where neighbour ``prev_sid`` hands off to ``this_sid``."""
    for i in range(len(loop)):
        if loop[i][2] == this_sid and loop[i - 1][2] == prev_sid:
            return int(loop[i][0])
    return None


def _corner_point(m, face, prev, this, loop, prev_sid, this_sid, tol, analytic):
    """Corner of three planes.

    A CAD vertex can sit on three planes whose normals are linearly dependent
    (two boundary edges of the face are parallel, so there is no unique
    triple intersection). Nearly parallel planes also have an algebraic
    intersection that lands far from the mesh. On a clean tessellation the
    mesh vertex is that corner when it lies on all three planes. On a noisy
    mesh the fitted intersection is the corner; swapping in the noisy vertex
    pulls the solid off the true planes.
    """
    corner = _intersect3([face, prev, this])
    if not analytic:
        return corner
    vert = _transition_vertex(loop, prev_sid, this_sid)
    point = None if vert is None else np.asarray(m.vertices[vert], float)
    near = max(4.0 * tol, 1e-3)
    if corner is not None and (point is None or float(np.linalg.norm(np.asarray(corner, float) - point)) <= near):
        return corner
    if point is None:
        return None
    # the vertex is a tessellation point of a fitted plane, so it can sit a
    # little further off a neighbour than the intersection test allows
    limit = max(8.0 * tol, 1e-3)
    for plane in (face, prev, this):
        if float(plane.distance(point[None])[0]) > limit:
            return None
    return point


def _runs(loop, min_run):
    labs = [s for _, _, s in loop]
    if not labs:
        return []
    # rotate so we start at a label change
    for i in range(len(labs)):
        if labs[i] != labs[i - 1]:
            labs = labs[i:] + labs[:i]
            break
    runs = []
    for s in labs:
        if runs and runs[-1][0] == s:
            runs[-1][1] += 1
        else:
            runs.append([s, 1])
    runs = _collapse_runs(runs)
    # A dense noisy boundary flickers: one or two edges of a foreign region
    # between two long runs of the same neighbour. A CAD edge is often a
    # single mesh edge between two different neighbours, and every run on
    # that loop is short — dropping the shortest one deletes the corner.
    if len(runs) > 1 and max(r[1] for r in runs) >= min_run:
        changed = True
        while changed and len(runs) > 1:
            changed = False
            victim = None
            for i in range(len(runs)):
                if runs[i][1] >= min_run:
                    continue
                prev = runs[(i - 1) % len(runs)]
                nxt = runs[(i + 1) % len(runs)]
                if prev[0] == nxt[0] and prev[1] >= runs[i][1] and nxt[1] >= runs[i][1]:
                    victim = i
                    break
            if victim is None:
                break
            prev_i = (victim - 1) % len(runs)
            runs[prev_i][1] += runs[victim][1]
            del runs[victim]
            runs = _collapse_runs(runs)
            changed = True
    return [r[0] for r in runs]


def _inside(poly, p):
    """Even-odd point-in-polygon test (2D)."""
    x, y = poly[:, 0], poly[:, 1]
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    c = ((y > p[1]) != (y2 > p[1])) & (p[0] < (x2 - x) * (p[1] - y) / (y2 - y + 1e-300) + x)
    return bool(c.sum() % 2)


def _unit_plane(plane):
    n = np.asarray(plane.normal, float)
    ln = float(np.linalg.norm(n))
    if ln < 1e-12:
        return None, None
    return n / ln, float(plane.d) / ln


def _same_plane(a, b, tol):
    """True when two fitted planes are the same geometric plane."""
    na, da = _unit_plane(a)
    nb, db = _unit_plane(b)
    if na is None or nb is None:
        return False
    if float(na @ nb) < 0.0:
        nb = -nb
        db = -db
    if float(na @ nb) < 0.9999:
        return False
    return abs(float(da - db)) <= max(4.0 * tol, 1e-4)


def _plane_sphere_circle(plane, sphere):
    """plane ∩ sphere → (center, unit normal, radius), or None if they miss."""
    n, d = _unit_plane(plane)
    if n is None or sphere.radius <= 0:
        return None
    h = float(n @ sphere.center - d)
    if not np.isfinite(h) or abs(h) >= sphere.radius * (1.0 - 1e-9):
        return None
    rho2 = sphere.radius ** 2 - h * h
    if rho2 <= 0:
        return None
    rho = float(np.sqrt(rho2))
    if rho <= 1e-8 * sphere.radius:
        return None
    return sphere.center - h * n, n, rho


def _make_sphere_face(sphere, planes, pts, normals):
    """Spherical cap (1 plane) or zone (2 parallel planes) bounded by circles."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.Geom import Geom_SphericalSurface
    from OCP.gp import gp_Ax3, gp_Dir
    R = float(sphere.radius)
    C = np.asarray(sphere.center, float)
    try:
        if len(planes) == 1:
            n, d = _unit_plane(planes[0])
            if n is None:
                return None
            h = float(n @ C - d)
            if abs(h) >= R * (1.0 - 1e-9):
                return None
            side = 1.0 if float(np.median(pts @ n - d)) >= 0.0 else -1.0
            axis = side * n
            v0 = float(np.arcsin(np.clip(-h * side / R, -1.0, 1.0)))
            if np.pi / 2 - v0 < 1e-5:
                return None
            ax = gp_Ax3(gp_Pnt(*map(float, C)), gp_Dir(*map(float, axis)))
            surf = Geom_SphericalSurface(ax, R)
            mf = BRepBuilderAPI_MakeFace(surf, 0.0, 2 * np.pi, v0, np.pi / 2, 1e-7)
        elif len(planes) == 2:
            n1, _d1 = _unit_plane(planes[0])
            n2, _d2 = _unit_plane(planes[1])
            if n1 is None or n2 is None or abs(float(n1 @ n2)) < 0.995:
                return None

            def latitude(plane):
                n, d = _unit_plane(plane)
                h = float(n @ C - d)
                if abs(h) >= R * (1.0 - 1e-9):
                    return None
                return float(np.arcsin(np.clip(-h * float(n @ n1) / R, -1.0, 1.0)))

            v1, v2 = latitude(planes[0]), latitude(planes[1])
            if v1 is None or v2 is None or abs(v1 - v2) < 1e-5:
                return None
            ax = gp_Ax3(gp_Pnt(*map(float, C)), gp_Dir(*map(float, n1)))
            surf = Geom_SphericalSurface(ax, R)
            mf = BRepBuilderAPI_MakeFace(surf, 0.0, 2 * np.pi, min(v1, v2), max(v1, v2), 1e-7)
        else:
            return None
        if not mf.IsDone():
            return None
        face = mf.Face()
    except Exception:
        return None
    radial = pts - C
    radial = radial / np.maximum(np.linalg.norm(radial, axis=1, keepdims=True), 1e-12)
    if float(np.mean(np.sum(normals * radial, axis=1))) < 0.0:
        face = face.Reversed()
    return face


def _plane_cone_circle(plane, cone):
    """plane ∩ cone for a plane perpendicular to the axis → (center, normal, radius)."""
    n, d = _unit_plane(plane)
    if n is None or abs(float(n @ cone.axis)) < 0.995:
        return None
    h = (d - float(n @ cone.apex)) / float(n @ cone.axis)
    if h <= 0:
        return None
    rho = float(cone.radius_at(h))
    if rho <= 1e-9:
        return None
    return cone.apex + cone.axis * h, n, rho


def _cyl_cone_circle(cyl, cone):
    """Coaxial cylinder ∩ cone → (center, axis, radius)."""
    if abs(float(cyl.axis @ cone.axis)) < 0.995:
        return None
    off = cone.apex - cyl.point
    if np.linalg.norm(off - cyl.axis * (off @ cyl.axis)) > 0.05 * cyl.radius:
        return None
    h = float(cyl.radius / np.tan(cone.half_angle))
    return cone.apex + cone.axis * h, cone.axis, float(cyl.radius)


def _ring(surface, other):
    """Circle where a cone/cylinder region meets a neighbouring region."""
    if isinstance(surface, Cone):
        if isinstance(other, Plane):
            return _plane_cone_circle(other, surface)
        if isinstance(other, Cylinder):
            return _cyl_cone_circle(other, surface)
    if isinstance(surface, Cylinder) and isinstance(other, Cone):
        return _cyl_cone_circle(surface, other)
    return None


def _make_cone_face(cone, hits):
    """Conical frustum face between two coaxial circles."""
    from OCP.Geom import Geom_ConicalSurface
    from OCP.gp import gp_Ax3, gp_Dir
    if len(hits) != 2 or any(h is None for h in hits):
        return None
    (c1, _, r1), (c2, _, r2) = hits
    z = c2 - c1
    H = float(np.linalg.norm(z))
    if H < 1e-9:
        return None
    z = z / H
    semi = float(np.arctan2(r2 - r1, H))
    try:
        surf = Geom_ConicalSurface(gp_Ax3(gp_Pnt(*map(float, c1)), gp_Dir(*map(float, z))), semi, float(r1))
        mf = BRepBuilderAPI_MakeFace(surf, 0.0, 2 * np.pi, 0.0, H / np.cos(semi), 1e-7)
        return mf.Face() if mf.IsDone() else None
    except Exception:
        return None


def _circle_plane_hits(center, normal, radius, plane_n, plane_d):
    """Points where a circle meets a plane. None when the plane misses it."""
    normal = np.asarray(normal, float)
    plane_n = np.asarray(plane_n, float)
    center = np.asarray(center, float)
    ln = float(np.linalg.norm(normal))
    pn = float(np.linalg.norm(plane_n))
    if ln < 1e-12 or pn < 1e-12 or radius <= 0:
        return None
    normal, plane_n = normal / ln, plane_n / pn
    in_plane = plane_n - normal * float(plane_n @ normal)
    ln = float(np.linalg.norm(in_plane))
    if ln < 1e-8:
        return None
    in_plane /= ln
    denom = float(plane_n @ in_plane)
    if abs(denom) < 1e-8:
        return None
    closest = center + in_plane * ((float(plane_d) - float(plane_n @ center)) / denom)
    closest = closest - normal * float(normal @ (closest - center))
    delta = float(np.linalg.norm(closest - center))
    if delta > radius + max(1e-7, 1e-6 * radius):
        return None
    half = float(np.sqrt(max(radius * radius - min(delta, radius) ** 2, 0.0)))
    if half <= max(1e-8, 1e-6 * radius):
        return [closest], 0.0
    line_dir = np.cross(normal, plane_n)
    line_dir = line_dir / np.linalg.norm(line_dir)
    return [closest + line_dir * half, closest - line_dir * half], half


def _plane_cyl_circle(plane, cyl):
    """plane ∩ cylinder when the plane is nearly perpendicular to the axis."""
    n, d = _unit_plane(plane)
    if n is None or cyl.radius <= 0 or abs(float(n @ cyl.axis)) < 0.85:
        return None
    center = cyl.point + cyl.axis * ((d - float(n @ cyl.point)) / float(n @ cyl.axis))
    return center, n, float(cyl.radius)


def _fillet_vertex(plane_a, plane_b, cyl, tol):
    """Corner where a fillet cylinder meets two planes.

    The corner is the cap plane meeting the generator that lies in the tangent
    wall. A circular section of radius ``cyl.radius`` misses that wall by
    about ``r α²/2`` once noise tilts the axis off the cap normal, which is
    enough for the circle test to report a miss.
    """
    na, da = _unit_plane(plane_a)
    nb, db = _unit_plane(plane_b)
    if na is None or nb is None or cyl.radius <= 0:
        return None
    axis = np.asarray(cyl.axis, float)
    axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
    aa, ab = abs(float(na @ axis)), abs(float(nb @ axis))
    if max(aa, ab) < 0.85:
        return None
    if aa >= ab:
        cap_n, cap_d, side_n, side_d = na, da, nb, db
    else:
        cap_n, cap_d, side_n, side_d = nb, db, na, da
    signed = float(side_n @ np.asarray(cyl.point, float) - side_d)
    gap = abs(abs(signed) - float(cyl.radius))
    # sewing uses 2*tol; a wall this close is the tangent plane, not a cut
    if (abs(signed) > 1e-12 and gap <= max(2.0 * tol, 1e-4)
            and abs(float(side_n @ axis)) <= 0.20 and abs(float(cap_n @ axis)) >= 0.85):
        radial = -side_n * np.sign(signed)
        radial = radial - axis * float(radial @ axis)
        ln = float(np.linalg.norm(radial))
        if ln >= 1e-8:
            radial /= ln
            origin = np.asarray(cyl.point, float) + radial * float(cyl.radius)
            origin = origin + axis * ((cap_d - float(cap_n @ origin)) / float(cap_n @ axis))
            if np.isfinite(origin).all():
                return origin
    ca = _plane_cyl_circle(plane_a, cyl)
    cb = _plane_cyl_circle(plane_b, cyl)
    if ca is None and cb is None:
        return None
    if ca is not None and (cb is None or abs(float(ca[1] @ axis)) >= abs(float(cb[1] @ axis))):
        center, normal, radius = ca
        sn, sd = nb, db
    else:
        center, normal, radius = cb
        sn, sd = na, da
    hit = _circle_plane_hits(center, normal, radius, sn, sd)
    if hit is None:
        return None
    pts, half = hit
    if half > max(8.0 * tol, 0.35 * radius):
        return None
    if abs(signed) < 1e-12:
        return None
    radial = -sn * np.sign(float(sn @ np.asarray(cyl.point, float) - sd))
    radial = radial - normal * float(radial @ normal)
    ln = float(np.linalg.norm(radial))
    if ln < 1e-8:
        return None
    expect = center + radial / ln * radius
    return min(pts, key=lambda p: float(np.linalg.norm(p - expect)))


def _line_edge(p0, p1):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
    from OCP.gp import gp_Pnt
    if float(np.linalg.norm(p0 - p1)) < 1e-8:
        return None
    edge = BRepBuilderAPI_MakeEdge(gp_Pnt(*map(float, p0)), gp_Pnt(*map(float, p1)))
    return edge.Edge() if edge.IsDone() else None


def _arc_edge(plane, cyl, p0, p1, hint):
    """Circular arc of plane ∩ cylinder from p0 to p1 through ``hint``.

    ``GC_MakeArcOfCircle(circle, p0, p1, sense)`` ignores the sense flag in
    this OCCT build and always emits the same sweep, which is the long arc
    when the fillet is concave. A point on the desired arc fixes the sweep.
    """
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
    from OCP.GC import GC_MakeArcOfCircle
    from OCP.gp import gp_Pnt
    circ = _plane_cyl_circle(plane, cyl)
    if circ is None:
        return None
    center, normal, radius = circ
    if float(np.linalg.norm(p0 - p1)) < 1e-8:
        return None
    hint = np.asarray(hint, float)
    hint = hint - normal * (float(normal @ hint) - float(normal @ center))
    radial = hint - center
    ln = float(np.linalg.norm(radial))
    if ln < 1e-9:
        return None
    mid = center + radial / ln * radius
    try:
        arc = GC_MakeArcOfCircle(
            gp_Pnt(*map(float, p0)), gp_Pnt(*map(float, mid)), gp_Pnt(*map(float, p1)))
    except Exception:
        return None
    if not arc.IsDone():
        return None
    edge = BRepBuilderAPI_MakeEdge(arc.Value())
    return edge.Edge() if edge.IsDone() else None


def _mixed_fillet_wire(m, face, seq, regions, tol, loop=None):
    """Planar wire whose runs alternate straight edges and fillet arcs."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeWire
    if len(seq) < 3:
        return None
    verts = []
    for i, sid in enumerate(seq):
        prev, this = regions[seq[i - 1]].surface, regions[sid].surface
        if isinstance(prev, Plane) and isinstance(this, Plane):
            vert = _intersect3([face, prev, this])
        elif isinstance(prev, Cylinder) and isinstance(this, Plane):
            vert = _fillet_vertex(face, this, prev, tol)
        elif isinstance(this, Cylinder) and isinstance(prev, Plane):
            vert = _fillet_vertex(face, prev, this, tol)
        else:
            return None
        if vert is None:
            return None
        verts.append(np.asarray(vert, float))
    edges = []
    for i, sid in enumerate(seq):
        surf = regions[sid].surface
        p0, p1 = verts[i], verts[(i + 1) % len(seq)]
        if isinstance(surf, Cylinder) and _plane_cyl_circle(face, surf) is not None:
            # a fillet's triangles sit out in the arc, so their centroid picks
            # the sweep. A hole wall's centroid sits on the axis and does not,
            # and the shared mesh edges do.
            cent = m.triangles_center[regions[sid].faces].mean(axis=0)
            axis = np.asarray(surf.axis, float)
            axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
            rel = cent - np.asarray(surf.point, float)
            radial = rel - axis * float(rel @ axis)
            sample = [m.vertices[a] for a, _, s in loop if s == sid] if loop is not None else []
            if float(np.linalg.norm(radial)) > 0.3 * float(surf.radius) or not sample:
                hint = cent
            else:
                hint = np.mean(sample, axis=0)
            edge = _arc_edge(face, surf, p0, p1, hint)
        else:
            edge = _line_edge(p0, p1)
        if edge is None:
            return None
        edges.append(edge)
    wire = BRepBuilderAPI_MakeWire()
    for edge in edges:
        wire.Add(edge)
    if not wire.IsDone():
        return None
    return wire.Wire()


def _rotate_about(axis, radial, angle):
    y = np.cross(axis, radial)
    y = y / np.linalg.norm(y)
    return np.cos(angle) * radial + np.sin(angle) * y


def _radial_from(cyl, point):
    d = np.asarray(point, float) - cyl.point
    d = d - cyl.axis * float(d @ cyl.axis)
    ln = float(np.linalg.norm(d))
    if ln < 1e-8:
        return None
    return d / ln


def _make_fillet_cylinder(m, region, loops, regions, tol):
    """Partial cylinder bounded by two tangent walls and two end planes."""
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.Geom import Geom_CylindricalSurface
    from OCP.gp import gp_Ax3, gp_Dir, gp_Pnt
    if len(loops) != 1:
        return None
    seq = _runs(loops[0], min_run=3)
    if len(seq) != 4:
        seq = _runs(loops[0], min_run=2)
    if len(seq) != 4 or any(not isinstance(regions[s].surface, Plane) for s in seq):
        return None
    cyl = region.surface
    planes = [regions[s].surface for s in seq]
    caps = []
    for plane in planes:
        n, _d = _unit_plane(plane)
        caps.append(n is not None and abs(float(n @ cyl.axis)) > 0.70)
    if caps.count(True) != 2 or any(caps[i] == caps[i - 1] for i in range(4)):
        return None
    verts = []
    for i in range(4):
        vert = _fillet_vertex(planes[i - 1], planes[i], cyl, tol)
        if vert is None:
            return None
        verts.append(np.asarray(vert, float))
    if float(cyl.distance(np.vstack(verts)).max()) > max(8.0 * tol, 0.2 * cyl.radius):
        return None
    rads = []
    for i, is_cap in enumerate(caps):
        if is_cap:
            continue
        rad = _radial_from(cyl, verts[i])
        if rad is None:
            return None
        rads.append(rad)
    if len(rads) != 2 or float(np.linalg.norm(np.cross(rads[0], rads[1]))) < 0.15:
        return None
    axis = cyl.axis / max(float(np.linalg.norm(cyl.axis)), 1e-12)
    y = np.cross(axis, rads[0])
    y = y / np.linalg.norm(y)

    def _ang(rad):
        a = float(np.arctan2(float(rad @ y), float(rad @ rads[0])))
        return a if a >= 0.0 else a + 2.0 * np.pi

    hint = _radial_from(cyl, m.triangles_center[region.faces].mean(axis=0))
    if hint is None:
        return None
    u_b, u_h = _ang(rads[1]), _ang(hint)
    sweep = u_b if u_h <= u_b + 1e-4 else u_b - 2.0 * np.pi
    if not 0.20 <= abs(sweep) <= np.pi * 1.15:
        return None
    hs = [float((v - cyl.point) @ axis) for v in verts]
    h0, h1 = min(hs), max(hs)
    if h1 - h0 < 1e-6:
        return None
    origin = cyl.point + axis * h0
    try:
        ax = gp_Ax3(gp_Pnt(*map(float, origin)), gp_Dir(*map(float, axis)), gp_Dir(*map(float, rads[0])))
        surf = Geom_CylindricalSurface(ax, float(cyl.radius))
        made = BRepBuilderAPI_MakeFace(surf, 0.0 if sweep > 0.0 else float(sweep),
                                       float(sweep) if sweep > 0.0 else 0.0, 0.0, float(h1 - h0), 1e-7)
    except Exception:
        return None
    if not made.IsDone():
        return None
    face = made.Face()
    mid = _rotate_about(axis, rads[0], sweep / 2.0)
    if float(np.mean(m.face_normals[region.faces], axis=0) @ mid) < 0.0:
        face = face.Reversed()
    return face


def _full_circle_wire(center, normal, radius):
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
    from OCP.gp import gp_Ax2, gp_Circ, gp_Dir, gp_Pnt
    n = np.asarray(normal, float)
    ln = float(np.linalg.norm(n))
    if ln < 1e-12 or radius <= 0:
        return None
    circ = gp_Circ(gp_Ax2(gp_Pnt(*map(float, center)), gp_Dir(*map(float, n / ln))), float(radius))
    edge = BRepBuilderAPI_MakeEdge(circ)
    if not edge.IsDone():
        return None
    wire = BRepBuilderAPI_MakeWire(edge.Edge())
    return wire.Wire() if wire.IsDone() else None


def _split_cap_wire(m, label, regions, loop, seq, cap, cyl, tol):
    """One circular wire split into an arc per coplanar cap fragment.

    The arc ends are the same fillet vertices the cap faces use, so sewing
    can fuse the cylinder to pads that do not share an edge with each other.
    """
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeWire
    points, hints = [], []
    for i, sid in enumerate(seq):
        prev = seq[i - 1]
        vert = _transition_vertex(loop, prev, sid)
        side = None
        if vert is not None:
            for f in np.asarray(m.vertex_faces[vert]).ravel():
                f = int(f)
                if f < 0:
                    continue
                rid = int(label[f])
                if rid < 0:
                    continue
                surf = regions[rid].surface
                if isinstance(surf, Plane) and not _same_plane(surf, cap, tol):
                    side = surf
                    break
        pt = _fillet_vertex(cap, side, cyl, tol) if side is not None else None
        if pt is None:
            return None
        points.append(np.asarray(pt, float))
        sample = [m.vertices[a] for a, _, s in loop if s == sid]
        hints.append(np.mean(sample, axis=0) if sample else points[-1])
    edges = []
    n = len(points)
    for i in range(n):
        edge = _arc_edge(cap, cyl, points[i], points[(i + 1) % n], hints[i])
        if edge is None:
            return None
        edges.append(edge)
    wire = BRepBuilderAPI_MakeWire()
    for edge in edges:
        wire.Add(edge)
    return wire.Wire() if wire.IsDone() else None


def build_planar(m: trimesh.Trimesh, regions: list[Region], tol: float):
    """Exact B-rep: plane polygons, cylindrical walls, spherical caps/zones.
    Returns None when the topology can't be resolved (caller falls back)."""
    from OCP.Geom import Geom_Plane
    from OCP.gp import gp_Dir, gp_Pln
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
    from OCP.Geom import Geom_CylindricalSurface
    from OCP.gp import gp_Ax2, gp_Ax3, gp_Circ
    if any(r.surface is None for r in regions):
        return None
    analytic = _analytic_mesh(m)
    label = -np.ones(len(m.faces), int)
    for i, r in enumerate(regions):
        label[r.faces] = i
    if (label < 0).any():
        return None
    sew = BRepBuilderAPI_Sewing(max(tol, 1e-6) * 2)
    def circle_on(plane, cyl):
        """plane ∩ cylinder when the axis is parallel to the plane normal."""
        if abs(plane.normal @ cyl.axis) < 0.995:
            return None
        n = plane.normal
        c = cyl.point + cyl.axis * ((plane.d - n @ cyl.point) / (n @ cyl.axis))
        return c, n

    for rid, r in enumerate(regions):
        loops = _loops(m, label, rid)
        if not loops:
            return None
        if isinstance(r.surface, Cylinder):
            # a fillet is one loop: two tangent walls and two end caps
            if len(loops) == 1:
                face = _make_fillet_cylinder(m, r, loops, regions, tol)
                if face is None:
                    return None
                sew.Add(face)
                continue
            # cylinder bounded by exactly two planar circular loops.
            # A circle may be split across several coplanar cap fragments
            # that do not touch; those ends are sewn as matching arcs.
            if len(loops) != 2:
                return None
            ends = []
            for loop in loops:
                seq = _runs(loop, min_run=3)
                if not seq:
                    return None
                surfs = [regions[s].surface for s in seq]
                if all(isinstance(s, Plane) for s in surfs) and all(
                        _same_plane(surfs[0], s, tol) for s in surfs):
                    o = surfs[0]
                elif len(seq) == 1:
                    o = surfs[0]
                else:
                    return None
                cc = (circle_on(o, r.surface) if isinstance(o, Plane) else _ring(r.surface, o))
                if cc is None:
                    return None
                ends.append((loop, seq, o, cc))
            split = any(len(seq) > 1 for _, seq, _, _ in ends)
            hs = [cc[0] for _, _, _, cc in ends]
            a = hs[1] - hs[0]
            h = np.linalg.norm(a)
            if h < 1e-9:
                return None
            a = a / h
            ax = gp_Ax3(gp_Pnt(*map(float, hs[0])), gp_Dir(*map(float, a)))
            surf = Geom_CylindricalSurface(ax, float(r.surface.radius))
            if not split:
                mf = BRepBuilderAPI_MakeFace(surf, 0.0, 2 * np.pi, 0.0, float(h), 1e-7)
                if not mf.IsDone():
                    return None
                sew.Add(mf.Face())
                continue
            if any(not isinstance(o, Plane) for _, _, o, _ in ends):
                return None
            wires = []
            for loop, seq, cap, cc in ends:
                if len(seq) == 1:
                    wire = _full_circle_wire(cc[0], cc[1], r.surface.radius)
                else:
                    wire = _split_cap_wire(m, label, regions, loop, seq, cap, r.surface, tol)
                if wire is None:
                    return None
                wires.append(wire)
            mf = BRepBuilderAPI_MakeFace(surf, wires[0], False)
            mf.Add(wires[1])
            if not mf.IsDone():
                mf = BRepBuilderAPI_MakeFace(surf, wires[0].Reversed(), False)
                mf.Add(wires[1])
            if not mf.IsDone():
                return None
            face = mf.Face()
            # the two wires can come back unorientable until the face is
            # rebuilt on the cylinder; keep the fix only when the area matches
            if not BRepCheck_Analyzer(face).IsValid():
                from OCP.BRepGProp import BRepGProp
                from OCP.GProp import GProp_GProps
                from OCP.ShapeFix import ShapeFix_Face
                fixer = ShapeFix_Face(face)
                fixer.SetPrecision(max(tol, 1e-6))
                fixer.SetMaxTolerance(max(5.0 * tol, 1e-4))
                fixer.Perform()
                fixed = fixer.Face()
                props = GProp_GProps()
                BRepGProp.SurfaceProperties_s(fixed, props)
                if (not BRepCheck_Analyzer(fixed).IsValid()
                        or abs(props.Mass() - r.area) > 0.25 * max(r.area, 1e-9)):
                    return None
                face = fixed
            sew.Add(face)
            continue
        if isinstance(r.surface, Sphere):
            # sphere bounded by planar circles (cap, or a zone between two planes)
            planes = []
            for loop in loops:
                seq = _runs(loop, min_run=3)
                if len(seq) != 1 or not isinstance(regions[seq[0]].surface, Plane):
                    return None
                if _plane_sphere_circle(regions[seq[0]].surface, r.surface) is None:
                    return None
                planes.append(regions[seq[0]].surface)
            pts = m.triangles_center[r.faces]
            face = _make_sphere_face(r.surface, planes, pts, m.face_normals[r.faces])
            if face is None:
                return None
            sew.Add(face)
            continue
        if isinstance(r.surface, Cone):
            hits = []
            for loop in loops:
                seq = _runs(loop, min_run=3)
                if len(seq) != 1:
                    return None
                hits.append(_ring(r.surface, regions[seq[0]].surface))
            face = _make_cone_face(r.surface, hits)
            if face is None:
                return None
            sew.Add(face)
            continue
        wires = []
        u, v = _plane_basis(r.surface.normal)
        for loop in loops:
            lp = m.vertices[[a for a, _, _ in loop]]
            q2 = np.c_[lp @ u, lp @ v]
            # signed area of the mesh loop: >0 outer boundary, <0 hole (u,v,n right-handed)
            sa = 0.5 * np.sum(q2[:, 0] * np.roll(q2[:, 1], -1) - np.roll(q2[:, 0], -1) * q2[:, 1])
            seq = _runs(loop, min_run=3)
            # a fillet run can be only a couple of edges once the wings have been
            # absorbed into the wall; min_run=3 would delete the arc
            seq_fine = _runs(loop, min_run=1)
            if any(isinstance(regions[k].surface, Cylinder) for k in seq_fine):
                seq = _runs(loop, min_run=2)
                fine_cyl = [k for k in seq_fine if isinstance(regions[k].surface, Cylinder)]
                if any(k not in seq for k in fine_cyl):
                    seq = seq_fine
            if len(seq) == 1 and isinstance(regions[seq[0]].surface, Cylinder):
                cyl = regions[seq[0]].surface
                cc = circle_on(r.surface, cyl)
                if cc is None:
                    return None
                circ = gp_Circ(gp_Ax2(gp_Pnt(*map(float, cc[0])), gp_Dir(*map(float, cc[1]))), float(cyl.radius))
                w = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire()
                wires.append((sa, q2, w))
                continue
            if len(seq) == 1 and isinstance(regions[seq[0]].surface, (Sphere, Cone)):
                o = regions[seq[0]].surface
                hit = (_plane_sphere_circle(r.surface, o) if isinstance(o, Sphere)
                       else _plane_cone_circle(r.surface, o))
                if hit is None:
                    return None
                cc, nrm, rho = hit
                circ = gp_Circ(gp_Ax2(gp_Pnt(*map(float, cc)), gp_Dir(*map(float, nrm))), float(rho))
                w = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire()
                wires.append((sa, q2, w))
                continue
            if any(isinstance(regions[k].surface, Cylinder) for k in seq):
                mixed = _mixed_fillet_wire(m, r.surface, seq, regions, tol, loop)
                if mixed is None:
                    return None
                wires.append((sa, q2, mixed))
                continue
            if len(seq) < 3 or any(not isinstance(regions[k].surface, Plane) for k in seq):
                return None
            pts = []
            for i in range(len(seq)):
                c = _corner_point(m, r.surface, regions[seq[i - 1]].surface,
                                  regions[seq[i]].surface, loop, seq[i - 1], seq[i], tol,
                                  analytic)
                if c is None:
                    return None
                c = np.asarray(c, float)
                if pts and float(np.linalg.norm(c - pts[-1])) <= 1e-7:
                    continue
                pts.append(c)
            if len(pts) >= 2 and float(np.linalg.norm(pts[0] - pts[-1])) <= 1e-7:
                pts.pop()
            if len(pts) < 3:
                return None
            poly = BRepBuilderAPI_MakePolygon()
            for p in pts:
                poly.Add(gp_Pnt(*map(float, p)))
            poly.Close()
            if not poly.IsDone():
                return None
            wires.append((sa, q2, poly.Wire()))
        from OCP.ShapeFix import ShapeFix_Face
        outers = [w for w in wires if w[0] > 0]
        holes = [w for w in wires if w[0] <= 0]
        if not outers:
            return None
        n = r.surface.normal
        o = n * r.surface.d
        pln = gp_Pln(gp_Pnt(*map(float, o)), gp_Dir(*map(float, n)))
        groups = {id(w): [] for w in outers}
        for h in holes:
            host = [w for w in outers if _inside(w[1], h[1][0])]
            if not host:
                return None
            groups[id(min(host, key=lambda w: w[0]))].append(h)
        for w in outers:
            mf = BRepBuilderAPI_MakeFace(pln, w[2], True)
            for h in groups[id(w)]:
                mf.Add(h[2])
            if not mf.IsDone():
                return None
            fx = ShapeFix_Face(mf.Face())
            fx.FixOrientation()
            fx.Perform()
            sew.Add(fx.Face())
    sew.Perform()
    shape = sew.SewedShape()
    exp = TopExp_Explorer(shape, TopAbs_SHELL)
    if not exp.More():
        return None
    ms = BRepBuilderAPI_MakeSolid()
    while exp.More():
        ms.Add(TopoDS.Shell_s(exp.Current()))
        exp.Next()
    fix = ShapeFix_Solid(ms.Solid())
    fix.Perform()
    return fix.Solid()
