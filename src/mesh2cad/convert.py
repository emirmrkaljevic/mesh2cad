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
    mesh: object = None
    tol: float = 0.0


def _exact(m, regions, tol):
    if not all(r.surface is not None for r in regions):
        return None
    shape = brep.build_planar(m, regions, tol)
    if shape is None:
        return None
    from .report import volume
    vm = abs(m.volume) if m.is_volume else None
    if not brep.is_valid(shape) or (vm and abs(volume(shape) - vm) > 0.05 * vm):
        return None
    return shape


def convert_mesh(m: trimesh.Trimesh, tol: float | None = None, cylinders: bool = True,
                 spheres: bool = True, cones: bool = True, angle_deg: float = 20.0, ref: trimesh.Trimesh | None = None,
                 with_report: bool = True) -> Result:
    m = clean(m)
    seg = dict(angle_deg=angle_deg, cylinders=cylinders, spheres=spheres, cones=cones)
    # without an explicit tol, also try tighter fits before falling back: small
    # cones/holes on clean meshes are often lost at 0.1% of the bbox diagonal
    tols = [tol] if tol is not None else [f * 1e-3 * m.scale for f in (1.0, 0.4, 0.2)]
    shape, method, first = None, "faceted-hybrid", None
    for t in tols:
        regions = segment(m, t, **seg)
        first = first or (t, regions)
        shape = _exact(m, regions, t)
        if shape is not None:
            tol, method = t, "exact"
            break
    if shape is None:
        tol, regions = first
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
    return Result(shape, regions, r, m, tol)


def _faceted_solid(m: trimesh.Trimesh, verts, lin_tol: float, repair_unify: bool = False):
    """Sew ``verts`` into a solid. Unify coplanar triangles when that stays valid.

    ShapeUpgrade_UnifySameDomain can segfault, flip a few faces (BRepCheck
    UnorientableShape), or, after vertex snapping, merge triangles that are not
    quite coplanar. Unify runs in a child process so a crash keeps the sewn
    solid. A broken unify is dropped. ``repair_unify`` retries ShapeFix, which
    recovers the raw mesh (inside-out faces) and not a snapped one.
    """
    sewn = brep.build_solid(m, verts, sew_tol=1e-6)
    if not brep.is_valid(sewn):
        sewn = brep.fix_shape(sewn)
    if not brep.is_valid(sewn):
        return None
    # keep an invalid unify so SameParameter / ShapeFix can still recover it;
    # a crash or a null result comes back as ``sewn`` itself
    unified = brep.safe_unify(sewn, lin_tol=lin_tol, keep_invalid=True)
    if unified is not sewn:
        unified = brep.same_parameter(unified, tol=1e-5)
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


def cadquery_file(src: str, dst: str, tol: float | None = None, **kw) -> dict:
    """Write an editable parametric CadQuery script for a simple part.

    Raises ``mesh2cad.cadquery_code.Unsupported`` when the part is outside the
    box/plate + pads + holes + bosses vocabulary or the script does not match."""
    from pathlib import Path

    from .cadquery_code import to_cadquery
    from .cadquery_code import Unsupported
    m = clean(load_mesh(src))
    # without an explicit tol, also try tighter fits: small cones and holes on
    # clean meshes are often lost at the default 0.1% of the bbox diagonal
    tols = [tol] if tol is not None else [f * 1e-3 * m.scale for f in (1.0, 0.4, 0.2)]
    for i, t in enumerate(tols):
        regions = segment(m, t, cylinders=kw.get("cylinders", True), spheres=kw.get("spheres", True),
                          cones=kw.get("cones", True))
        try:
            code, q = to_cadquery(m, regions, t, source=Path(src).name)
            break
        except Unsupported:
            if i == len(tols) - 1:
                raise
    Path(dst).write_text(code)
    params = {}
    for line in code.split("# --- parameters ---")[1].split("# --- model ---")[0].splitlines():
        if "=" in line:
            k, v = line.split("#")[0].split("=")
            params[k.strip()] = float(v)
    return {"output_path": dst, "parameters": params, "check": q}
