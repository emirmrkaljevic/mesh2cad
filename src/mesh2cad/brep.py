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

from .fit import Plane
from .segment import Region


def snap_vertices(m: trimesh.Trimesh, regions: list[Region], tol: float, iters: int = 10):
    nv = len(m.vertices)
    vregs = [[] for _ in range(nv)]
    for i, r in enumerate(regions):
        if r.surface is None:
            continue
        for v in np.unique(m.faces[r.faces]):
            vregs[v].append(i)
    x0 = m.vertices.copy()
    x = x0.copy()
    lam = 1e-6
    for v in range(nv):
        rs = vregs[v]
        if not rs:
            continue
        surfs = [regions[i].surface for i in rs]
        planes = [s for s in surfs if isinstance(s, Plane)]
        others = [s for s in surfs if not isinstance(s, Plane)]
        p = x0[v]
        for _ in range(iters if others else 1):
            # linearize every surface at the current point as a plane
            N, D = [], []
            for s in planes:
                N.append(s.normal); D.append(s.d)
            for s in others:
                q = s.project(p[None])[0]
                n = p - q
                ln = np.linalg.norm(n)
                if ln < 1e-12:
                    r = s._radial(p[None])[0]
                    n = r / max(np.linalg.norm(r), 1e-12)
                else:
                    n = n / ln
                N.append(n); D.append(n @ q)
            N, D = np.array(N), np.array(D)
            A = N.T @ N + lam * np.eye(3)
            b = N.T @ D + lam * p
            pn = np.linalg.solve(A, b)
            # guard against wildly ill-conditioned (near-parallel) solves
            if np.linalg.norm(pn - x0[v]) > 4 * tol + 1e-9 and len(N) > 1:
                pn = surfs[0].project(x0[v][None])[0]
            if np.linalg.norm(pn - p) < 1e-10:
                p = pn
                break
            p = pn
        x[v] = p
    return x


def relax(m, x, vregs, regions, tol, iters=60):
    """Untangle triangles: move vertices toward neighbour mean while keeping
    them on their constraint surfaces (linearised as planes at the snapped
    point). Corners (3+ independent planes) stay put; interior vertices are
    irrelevant to the final faces, so smoothing them is free."""
    nv = len(x)
    lam = 1.0
    P = np.zeros((nv, 3, 3))
    C = np.zeros((nv, 3))
    active = np.zeros(nv, bool)
    for v in range(nv):
        if not vregs[v]:
            continue
        N, D = [], []
        for i in vregs[v]:
            sfc = regions[i].surface
            if isinstance(sfc, Plane):
                n, d = sfc.normal, sfc.d
            else:
                r = sfc._radial(x[v][None])[0]
                n = r / max(np.linalg.norm(r), 1e-12)
                d = n @ x[v]
            N.append(n); D.append(d)
        N, D = np.array(N), np.array(D)
        # cap rank: drop near-duplicate constraints
        A = N.T @ N * 1e6 + lam * np.eye(3)
        Ai = np.linalg.inv(A)
        P[v] = Ai * lam
        C[v] = Ai @ (N.T @ D * 1e6)
        active[v] = True
    edges = m.edges_unique
    deg = np.bincount(edges.ravel(), minlength=nv).astype(float)
    deg[deg == 0] = 1
    idx = np.where(active)[0]
    for _ in range(iters):
        s = np.zeros_like(x)
        np.add.at(s, edges[:, 0], x[edges[:, 1]])
        np.add.at(s, edges[:, 1], x[edges[:, 0]])
        mean = s / deg[:, None]
        x[idx] = np.einsum("nij,nj->ni", P[idx], mean[idx]) + C[idx]
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


def build_planar(m: trimesh.Trimesh, regions: list[Region], tol: float):
    """Exact planar B-rep: polygon corners = triple-plane intersections.
    Returns None when the topology can't be resolved (caller falls back)."""
    from OCP.Geom import Geom_Plane
    from OCP.gp import gp_Dir, gp_Pln
    if any(not isinstance(r.surface, Plane) for r in regions):
        return None
    label = -np.ones(len(m.faces), int)
    for i, r in enumerate(regions):
        label[r.faces] = i
    if (label < 0).any():
        return None
    sew = BRepBuilderAPI_Sewing(max(tol, 1e-6) * 2)
    for rid, r in enumerate(regions):
        loops = _loops(m, label, rid)
        if not loops:
            return None
        wires = []
        for loop in loops:
            seq = _runs(loop, min_run=3)
            if len(seq) < 3:
                return None
            pts = []
            for i in range(len(seq)):
                c = _intersect3([r.surface, regions[seq[i - 1]].surface, regions[seq[i]].surface])
                if c is None:
                    return None
                pts.append(c)
            u, v = _plane_basis(r.surface.normal)
            q = np.array([[p @ u, p @ v] for p in pts])
            area = 0.5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - np.roll(q[:, 0], -1) * q[:, 1])
            poly = BRepBuilderAPI_MakePolygon()
            for p in pts:
                poly.Add(gp_Pnt(*map(float, p)))
            poly.Close()
            if not poly.IsDone():
                return None
            wires.append((abs(area), poly.Wire()))
        wires.sort(key=lambda w: -w[0])
        n = r.surface.normal
        o = n * r.surface.d
        pln = gp_Pln(gp_Pnt(*map(float, o)), gp_Dir(*map(float, n)))
        mf = BRepBuilderAPI_MakeFace(pln, wires[0][1], True)
        for _, w in wires[1:]:
            mf.Add(w)
        if not mf.IsDone():
            return None
        from OCP.ShapeFix import ShapeFix_Face
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
