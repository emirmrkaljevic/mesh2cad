"""Build a B-rep solid from segmented mesh and export STEP.

Strategy (v0.1): snap every mesh vertex onto the fitted surfaces of all
regions it touches (corner -> plane intersection point, edge -> intersection
line), rebuild triangles as OCC faces, sew, then merge co-surface faces with
ShapeUpgrade_UnifySameDomain so each fitted plane becomes one clean face.
Unfitted regions stay faceted (hybrid fallback).
"""
from __future__ import annotations

import numpy as np
import trimesh
from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon,
                                BRepBuilderAPI_MakeSolid, BRepBuilderAPI_Sewing)
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.gp import gp_Pnt
from OCP.IFSelect import IFSelect_RetDone
from OCP.ShapeFix import ShapeFix_Solid
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer
from OCP.TopAbs import TopAbs_SHELL
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Shape

from .fit import Cylinder, Plane
from .segment import Region


def snap_vertices(m: trimesh.Trimesh, regions: list[Region], tol: float, iters: int = 10):
    """Move vertices onto the fitted surfaces of every region they touch.
    Plane-only vertices are solved in one batched least-squares step;
    vertices touching a cylinder use a few linearised iterations."""
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
    has_cyl = np.zeros(nv, bool)
    has_cyl[rv[~is_plane[rr]]] = True
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
    for _ in range(iters if has_cyl.any() else 1):
        N = np.zeros((len(rv), 3))
        D = np.zeros(len(rv))
        N[pl], D[pl] = Np[rr[pl]], Dp[rr[pl]]
        for i in np.unique(rr[~pl]):
            sel = (rr == i)
            c = regions[i].surface
            p = x[rv[sel]]
            q = c.project(p)
            rad = c._radial(p)
            n = rad / np.maximum(np.linalg.norm(rad, axis=1, keepdims=True), 1e-12)
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
    return shape


def unify(shape: TopoDS_Shape, lin_tol: float, ang_tol: float = 1e-3) -> TopoDS_Shape:
    u = ShapeUpgrade_UnifySameDomain(shape, True, True, True)
    u.SetLinearTolerance(lin_tol)
    u.SetAngularTolerance(ang_tol)
    u.Build()
    return u.Shape()


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


def _runs(loop, min_run):
    labs = [s for _, _, s in loop]
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
    # drop flicker runs
    changed = True
    while changed and len(runs) > 1:
        changed = False
        i = min(range(len(runs)), key=lambda k: runs[k][1])
        if runs[i][1] < min_run:
            runs[i - 1][1] += runs[i][1]
            runs.pop(i)
            # merge equal neighbours
            j = 0
            while j < len(runs) and len(runs) > 1:
                k = (j + 1) % len(runs)
                if runs[j][0] == runs[k][0] and j != k:
                    runs[j][1] += runs[k][1]
                    runs.pop(k)
                    j = 0
                else:
                    j += 1
            changed = True
    if runs and runs[0][0] == runs[-1][0] and len(runs) > 1:
        runs[0][1] += runs.pop()[1]
    return [r[0] for r in runs]


def _inside(poly, p):
    """Even-odd point-in-polygon test (2D)."""
    x, y = poly[:, 0], poly[:, 1]
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    c = ((y > p[1]) != (y2 > p[1])) & (p[0] < (x2 - x) * (p[1] - y) / (y2 - y + 1e-300) + x)
    return bool(c.sum() % 2)


def build_planar(m: trimesh.Trimesh, regions: list[Region], tol: float):
    """Exact planar B-rep: polygon corners = triple-plane intersections.
    Returns None when the topology can't be resolved (caller falls back)."""
    from OCP.Geom import Geom_Plane
    from OCP.gp import gp_Dir, gp_Pln
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
    from OCP.Geom import Geom_CylindricalSurface
    from OCP.gp import gp_Ax2, gp_Ax3, gp_Circ
    if any(r.surface is None for r in regions):
        return None
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
            # cylinder bounded by exactly two planar circular loops
            if len(loops) != 2:
                return None
            hs = []
            for loop in loops:
                seq = _runs(loop, min_run=3)
                if len(seq) != 1 or not isinstance(regions[seq[0]].surface, Plane):
                    return None
                cc = circle_on(regions[seq[0]].surface, r.surface)
                if cc is None:
                    return None
                hs.append(cc[0])
            a = hs[1] - hs[0]
            h = np.linalg.norm(a)
            a = a / h
            ax = gp_Ax3(gp_Pnt(*map(float, hs[0])), gp_Dir(*map(float, a)))
            surf = Geom_CylindricalSurface(ax, float(r.surface.radius))
            mf = BRepBuilderAPI_MakeFace(surf, 0.0, 2 * np.pi, 0.0, float(h), 1e-7)
            if not mf.IsDone():
                return None
            sew.Add(mf.Face())
            continue
        wires = []
        u, v = _plane_basis(r.surface.normal)
        for loop in loops:
            lp = m.vertices[[a for a, _, _ in loop]]
            q2 = np.c_[lp @ u, lp @ v]
            # signed area of the mesh loop: >0 outer boundary, <0 hole (u,v,n right-handed)
            sa = 0.5 * np.sum(q2[:, 0] * np.roll(q2[:, 1], -1) - np.roll(q2[:, 0], -1) * q2[:, 1])
            seq = _runs(loop, min_run=3)
            if len(seq) == 1 and isinstance(regions[seq[0]].surface, Cylinder):
                cyl = regions[seq[0]].surface
                cc = circle_on(r.surface, cyl)
                if cc is None:
                    return None
                circ = gp_Circ(gp_Ax2(gp_Pnt(*map(float, cc[0])), gp_Dir(*map(float, cc[1]))), float(cyl.radius))
                w = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire()
                wires.append((sa, q2, w))
                continue
            if len(seq) < 3 or any(not isinstance(regions[k].surface, Plane) for k in seq):
                return None
            pts = []
            for i in range(len(seq)):
                c = _intersect3([r.surface, regions[seq[i - 1]].surface, regions[seq[i]].surface])
                if c is None:
                    return None
                pts.append(c)
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
