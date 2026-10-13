"""High-level pipeline."""
from __future__ import annotations

from dataclasses import dataclass

import trimesh

from . import brep, report as rep
from .io import clean, load_mesh
from .segment import segment


@dataclass
class Result:
    shape: object
    regions: list
    report: dict


def convert_mesh(m: trimesh.Trimesh, tol: float | None = None, cylinders: bool = True,
                 angle_deg: float = 20.0, ref: trimesh.Trimesh | None = None,
                 with_report: bool = True) -> Result:
    m = clean(m)
    if tol is None:
        tol = 1e-3 * m.scale
    regions = segment(m, tol, angle_deg=angle_deg, cylinders=cylinders)
    shape, method = None, "faceted-hybrid"
    if all(r.surface is not None for r in regions):
        shape = brep.build_planar(m, regions, tol)
        if shape is not None:
            from .report import volume
            vm = abs(m.volume) if m.is_volume else None
            if not brep.is_valid(shape) or (vm and abs(volume(shape) - vm) > 0.05 * vm):
                shape = None
            else:
                method = "exact"
    if shape is None:
        verts = brep.snap_vertices(m, regions, tol)
        shape = brep.build_solid(m, verts, sew_tol=max(tol, 1e-6))
        shape = brep.unify(shape, lin_tol=max(tol * 0.05, 1e-6))
        if not brep.is_valid(shape):
            # snapping tangled triangles: export the raw mesh as a faithful faceted solid
            shape = brep.unify(brep.build_solid(m, m.vertices, sew_tol=1e-6), lin_tol=1e-6)
            method = "faceted"
    valid = brep.is_valid(shape)
    r = rep.build_report(m, shape, regions, valid, ref=ref) if with_report else {}
    if r:
        r["method"] = method
    return Result(shape, regions, r)


def convert_file(src: str, dst: str, tol: float | None = None, **kw) -> Result:
    res = convert_mesh(load_mesh(src), tol=tol, **kw)
    brep.write_step(res.shape, dst)
    return res
