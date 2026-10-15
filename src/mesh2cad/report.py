"""Fidelity report: compare output B-rep to input mesh."""
from __future__ import annotations

import json
from collections import Counter

import numpy as np
import trimesh
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.GeomAbs import GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Plane, GeomAbs_Sphere
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_FACE, TopAbs_REVERSED
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS


def face_types(shape) -> dict:
    c = Counter()
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        t = BRepAdaptor_Surface(TopoDS.Face_s(exp.Current())).GetType()
        kind = {GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder", GeomAbs_Sphere: "sphere",
                GeomAbs_Cone: "cone"}.get(t, "other")
        c[kind] += 1
        exp.Next()
    return dict(c)


def shape_to_mesh(shape, deflection: float) -> trimesh.Trimesh:
    BRepMesh_IncrementalMesh(shape, deflection, False, 0.3, True)
    V, F, off = [], [], 0
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        face = TopoDS.Face_s(exp.Current())
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(face, loc)
        if tri is not None:
            tr = loc.Transformation()
            for i in range(1, tri.NbNodes() + 1):
                p = tri.Node(i).Transformed(tr)
                V.append((p.X(), p.Y(), p.Z()))
            rev = face.Orientation() == TopAbs_REVERSED
            for i in range(1, tri.NbTriangles() + 1):
                a, b, c = tri.Triangle(i).Get()
                F.append((a - 1 + off, c - 1 + off, b - 1 + off) if rev else (a - 1 + off, b - 1 + off, c - 1 + off))
            off += tri.NbNodes()
        exp.Next()
    return trimesh.Trimesh(np.array(V), np.array(F), process=True)


def volume(shape) -> float:
    p = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, p)
    return abs(p.Mass())


def deviation(a: trimesh.Trimesh, b: trimesh.Trimesh, n: int = 20000, seed: int = 0) -> dict:
    """Two-sided sampled distances between surfaces."""
    rng = np.random.default_rng(seed)
    pa = trimesh.sample.sample_surface(a, n, seed=rng)[0]
    pb = trimesh.sample.sample_surface(b, n, seed=rng)[0]
    dab = np.abs(trimesh.proximity.closest_point(b, pa)[1])
    dba = np.abs(trimesh.proximity.closest_point(a, pb)[1])
    d = np.r_[dab, dba]
    return {"hausdorff": float(d.max()), "mean": float(d.mean()),
            "rms": float(np.sqrt((d ** 2).mean())), "p95": float(np.percentile(d, 95))}


def build_report(src: trimesh.Trimesh, shape, regions, valid: bool, ref=None) -> dict:
    """ref: optional ground-truth mesh to compare against instead of src."""
    out_mesh = shape_to_mesh(shape, 0.01 * src.scale)
    target = ref if ref is not None else src
    total = src.area
    fitted = Counter()
    for r in regions:
        fitted[r.surface.kind if r.surface is not None else "faceted"] += r.area
    vin = abs(target.volume) if target.is_volume else float("nan")
    vout = volume(shape)
    return {
        "input": {"vertices": len(src.vertices), "faces": len(src.faces),
                  "watertight": bool(src.is_watertight)},
        "regions": {"count": len(regions),
                    "by_type": dict(Counter(r.surface.kind if r.surface else "faceted" for r in regions))},
        "area_fraction": {k: v / total for k, v in fitted.items()},
        "output": {"faces": face_types(shape), "valid_solid": valid, "volume": vout},
        "volume_error_pct": float(abs(vout - vin) / vin * 100) if vin == vin and vin > 0 else None,
        "deviation": deviation(target, out_mesh),
    }


def to_markdown(r: dict) -> str:
    d = r["deviation"]
    lines = ["# mesh2cad fidelity report", "",
             f"- Input: {r['input']['faces']} triangles, watertight={r['input']['watertight']}",
             f"- Regions: {r['regions']['count']} {r['regions']['by_type']}",
             f"- Output faces: {r['output']['faces']}, valid solid: {r['output']['valid_solid']}",
             "- Area fitted: " + ", ".join(f"{k} {v:.1%}" for k, v in r["area_fraction"].items()),
             f"- Volume error: {r['volume_error_pct']:.3f}%" if r["volume_error_pct"] is not None else "- Volume error: n/a",
             f"- Deviation: Hausdorff {d['hausdorff']:.4f}, mean {d['mean']:.4f}, RMS {d['rms']:.4f}, p95 {d['p95']:.4f}"]
    return "\n".join(lines) + "\n"


def dumps(r: dict) -> str:
    return json.dumps(r, indent=2)
