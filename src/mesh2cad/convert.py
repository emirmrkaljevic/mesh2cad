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
                 spheres: bool = True, cones: bool = True, angle_deg: float = 20.0, ref: trimesh.Trimesh | None = None,
                 with_report: bool = True) -> Result:
    m = clean(m)
    if tol is None:
        tol = 1e-3 * m.scale
    regions = segment(m, tol, angle_deg=angle_deg, cylinders=cylinders, spheres=spheres, cones=cones)
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
        # Sew with a tight tolerance: the mesh is already closed, and a sew
        # tolerance on the order of the noise welds distinct vertices together.
        shape = _faceted_solid(m, brep.snap_vertices(m, regions, tol), lin_tol=max(tol * 0.05, 1e-6))
        method = "faceted-hybrid"
        if shape is None:
            # raw triangles are already coincident, so a broken unify is usually
            # just inside-out faces and ShapeFix can put them back
            shape = _faceted_solid(m, m.vertices, lin_tol=1e-6, repair_unify=True)
            method = "faceted"
        if shape is None:
            shape = brep.build_solid(m, m.vertices, sew_tol=1e-6)
            method = "faceted"
    valid = brep.is_valid(shape)
    r = rep.build_report(m, shape, regions, valid, ref=ref) if with_report else {}
    if r:
        r["method"] = method
    return Result(shape, regions, r)


def _faceted_solid(m: trimesh.Trimesh, verts, lin_tol: float, repair_unify: bool = False):
    """Sew ``verts`` into a solid. Unify coplanar triangles when that stays valid.

    ShapeUpgrade_UnifySameDomain can flip a few faces (BRepCheck
    UnorientableShape) or, after vertex snapping, merge triangles that are not
    quite coplanar. A broken unify is dropped. ``repair_unify`` retries ShapeFix,
    which recovers the raw mesh (inside-out faces) and not a snapped one.
    """
    sewn = brep.build_solid(m, verts, sew_tol=1e-6)
    if not brep.is_valid(sewn):
        sewn = brep.fix_shape(sewn)
    if not brep.is_valid(sewn):
        return None
    unified = brep.same_parameter(brep.unify(sewn, lin_tol=lin_tol), tol=1e-5)
    if brep.is_valid(unified):
        return unified
    if repair_unify:
        fixed = brep.same_parameter(brep.fix_shape(unified), tol=1e-5)
        if brep.is_valid(fixed):
            return fixed
    return sewn


def convert_file(src: str, dst: str, tol: float | None = None, **kw) -> Result:
    res = convert_mesh(load_mesh(src), tol=tol, **kw)
    brep.write_step(res.shape, dst)
    return res
