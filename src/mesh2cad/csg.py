"""Small booleans of straight extrusions, for parametric CadQuery and OpenSCAD.

A part qualifies when every region is a plane, an axis-aligned cylinder, or a
45 degree axis-aligned cone, and the solid is one of:

* the intersection of outlines extruded along X, Y, or Z;
* one of those extrusions minus straight cut extrusions along the other axes,
  including a blind pocket that opens on one face;
* either of those, plus round through-holes, counterbored through-holes,
  blind holes, and a 45 degree chamfer on a hole mouth or a blind hole bottom;
* a planar clip where one face is tilted to all three axes (at most four).

The box vocabulary and a single straight extrusion are tried first. This
module is the fallback, and it keeps their parameter names when the part is
still one extrusion (``extrude_len``, ``outline``, ``hole_N_d``).
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
import trimesh
from scipy import ndimage
from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon,
                                BRepBuilderAPI_Sewing)
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import (BRepPrimAPI_MakeCone, BRepPrimAPI_MakeCylinder,
                             BRepPrimAPI_MakePrism)
from OCP.BRepTools import BRepTools, BRepTools_WireExplorer
from OCP.BOPAlgo import BOPAlgo_Builder
from OCP.GeomAbs import GeomAbs_Line
from OCP.GProp import GProp_GProps
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_SOLID, TopAbs_WIRE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt, gp_Vec

from .extrude import (
    _AXIS_NAME, _PROFILE_AXES, _attitude, _drop_collinear, _drop_near, _eps,
    _loop_polygon, _loops, _merge_walls, _min_area, _nice_loop, _outward, _partial_span,
    _perp_dist, _plane_coord, _polygons_from_cap, _project, _shoelace, _unit,
)
from .fit import Cone, Cylinder, Plane, Sphere
from .parametric import (
    Segment, Unsupported, _cluster_end, _nice, _plain_or_counterbore, _principal_axis,
    _snap, axis_point_expr, cadquery_axis_dir, openscad_axis_rotate,
)
from .report import deviation, shape_to_mesh, volume
from .segment import _analytic_mesh

# Outline axis of an extrusion: X stores (y, z), Y stores (x, z), Z stores (x, y).
_OUTLINE_OF = {0: {1: 0, 2: 1}, 1: {0: 0, 2: 1}, 2: {0: 0, 1: 1}}
_PLANE = {0: ("YZ", "y, z"), 1: ("XZ", "x, z"), 2: ("XY", "x, y")}
# Prefer a Z extrusion as the one that survives when the others are redundant.
_DROP_ORDER = (1, 0, 2)
# Extrusions that meet on a face are overlapped by this much in the script so
# the intersection stays a manifold. The other outline stops at the true end,
# so the solid does not grow. Small on purpose: a length edit should still
# change the volume by very nearly the outline area times the edit.
_OVERLAP = 0.02


class _Cloud:
    """Vertex array ``_loop_polygon`` can index, without building a mesh."""

    def __init__(self, vertices):
        self.vertices = np.asarray(vertices, float)


def _refuse_surfaces(m, regions):
    """Reject geometry this boolean cannot name. The message says which."""
    if any(r.surface is None for r in regions):
        raise Unsupported("part has unfitted (faceted) regions")
    for r in regions:
        s = r.surface
        if isinstance(s, Sphere):
            raise Unsupported("sphere faces are not supported yet")
        if isinstance(s, Cone):
            # A 45 degree cone on a principal axis is a hole chamfer. Anything
            # else (a draft, a countersink that is not 45 degrees) stays out.
            ang = abs(float(s.half_angle) - math.pi / 4.0)
            if _principal_axis(s.axis) is None or ang > math.radians(3.0):
                raise Unsupported(
                    "a countersink, chamfer or draft on a combined extrusion is not supported yet")
            continue
        if not isinstance(s, (Plane, Cylinder)):
            raise Unsupported(f"{s.kind} faces are not supported yet")
        if isinstance(s, Cylinder):
            if _principal_axis(s.axis) is None:
                raise Unsupported("curved surface whose axis is not aligned with X, Y or Z")
            if r.area < _min_area(m):
                continue
            span = _partial_span(m, s, r.faces)
            if span is not None and span < np.radians(250.0):
                raise Unsupported("a rounded edge is not part of an editable outline")
            if _outward(m, r, s):
                raise Unsupported("a boss on an extruded outline is not supported yet")


def _annulus(m, region, axis, point, radius, tol) -> bool:
    """A planar ring or disk around ``axis``: every vertex stays inside the bore."""
    if _attitude(region.surface.normal, axis, 4.0) != "cap":
        return False
    verts = m.vertices[np.unique(m.faces[region.faces])]
    rel = verts - np.asarray(point, float)
    rel[:, axis] = 0.0
    return float(np.max(np.linalg.norm(rel, axis=1))) <= float(radius) + _eps(tol, 0.35)


def _cap_coords(regions, axis, skip):
    coords = []
    for r in regions:
        if id(r) in skip or not isinstance(r.surface, Plane):
            continue
        if _attitude(r.surface.normal, axis, 4.0) != "cap":
            continue
        coords.append(_plane_coord(r.surface, axis))
    return coords


def _points_inside(mesh, pts):
    pts = np.asarray(pts, float)
    if len(pts) == 0:
        return np.zeros(0, dtype=bool)
    try:
        return np.asarray(mesh.contains(pts), dtype=bool)
    except Exception:
        return np.ones(len(pts), dtype=bool)


def _part_span(m, axis, exterior, tol):
    """Axial extent of the part, snapped to a cap that actually sits on the end."""
    coords = m.vertices[:, axis]
    lo, hi = float(coords.min()), float(coords.max())
    band = max(8.0 * tol, 0.4)
    if exterior:
        c = min(exterior, key=lambda q: abs(q - lo))
        if abs(c - lo) <= band:
            lo = float(c)
        c = min(exterior, key=lambda q: abs(q - hi))
        if abs(c - hi) <= band:
            hi = float(c)
    return lo, hi


def _split_chains(segs, split_gap):
    """Coaxial pieces with a wide empty gap are different holes, not one step."""
    if not segs:
        return []
    chains = [[segs[0]]]
    for a, b in zip(segs, segs[1:]):
        if b.za - a.zb > split_gap:
            chains.append([b])
        else:
            chains[-1].append(b)
    return chains


def _span_clear(mesh, point, axis, lo, hi, occupied):
    """True when the centerline between ``lo`` and ``hi`` misses the solid.

    Occupied intervals are the cylinders already accounted for. A hole that
    only shows up on either side of a channel still goes through when the
    channel itself is empty.
    """
    samples = []
    t = float(lo) + 0.25
    while t < float(hi) - 0.2:
        if not any(a - 0.2 <= t <= b + 0.2 for a, b in occupied):
            samples.append(t)
        t += 0.45
    if not samples:
        return True
    pts = np.repeat(np.asarray(point, float)[None, :], len(samples), axis=0)
    pts[:, axis] = samples
    inside = _points_inside(mesh, pts)
    return not bool(np.any(inside))


def _extend_small(chain, lo, hi, r_small, gap):
    """Run the smallest radius across the whole span, keeping larger steps."""
    big = []
    for s in chain:
        if float(s.ra) > r_small + gap:
            za, zb = float(s.za), float(s.zb)
            if abs(za - lo) <= gap:
                za = float(lo)
            if abs(zb - hi) <= gap:
                zb = float(hi)
            big.append(Segment("cylinder", za, zb, float(s.ra), float(s.ra)))
    big.sort(key=lambda s: s.za)
    spans = []
    cursor = float(lo)
    for s in big:
        if s.za > cursor + 1e-6:
            spans.append(Segment("cylinder", cursor, s.za, r_small, r_small))
        spans.append(s)
        cursor = float(s.zb)
    if float(hi) > cursor + 1e-6:
        spans.append(Segment("cylinder", cursor, float(hi), r_small, r_small))
    return spans


def _chain_opened(chain, exterior, gap):
    """The existing through-hole test, on one chain."""
    za = min(s.za for s in chain)
    zb = max(s.zb for s in chain)
    lo = _snap(za, exterior, gap)
    hi = _snap(zb, exterior, gap)
    opened = abs(lo - za) <= 1e-9 and abs(hi - zb) <= 1e-9 and hi > lo + gap
    if exterior:
        lo_hit = min(exterior, key=lambda c: abs(c - za))
        hi_hit = min(exterior, key=lambda c: abs(c - zb))
        opened = abs(lo_hit - za) <= gap and abs(hi_hit - zb) <= gap and hi_hit > lo_hit + gap
        lo, hi = (lo_hit, hi_hit) if opened else (za, zb)
    return opened, lo, hi


def _blind_hole(chain, exterior, gap):
    """A single radius that opens on exactly one exterior cap."""
    if any(abs(float(s.ra) - float(chain[0].ra)) > gap for s in chain):
        return None
    za = min(float(s.za) for s in chain)
    zb = max(float(s.zb) for s in chain)
    if not exterior or zb <= za + gap:
        return None
    lo_hit = min(exterior, key=lambda c: abs(c - za))
    hi_hit = min(exterior, key=lambda c: abs(c - zb))
    on_lo = abs(lo_hit - za) <= gap
    on_hi = abs(hi_hit - zb) <= gap
    radius = float(chain[0].ra)
    if on_hi and not on_lo:
        return {"d": 2.0 * radius, "depth": float(hi_hit - za), "from": 1}
    if on_lo and not on_hi:
        return {"d": 2.0 * radius, "depth": float(zb - lo_hit), "from": -1}
    return None


def _covered_by_through(chain, emitted, axis, point, gap):
    """A later piece of a hole that was already extended through the part."""
    rmax = max(float(s.ra) for s in chain)
    za = min(float(s.za) for s in chain)
    zb = max(float(s.zb) for s in chain)
    for prev in emitted:
        if prev["axis"] != axis or _perp_dist(axis, prev["point"], point) > gap:
            continue
        if prev["radius"] + gap < rmax:
            continue
        if za >= prev["lo"] - gap and zb <= prev["hi"] + gap:
            return True
    return False


def _cone_ends(m, region, axis):
    """Axial ends of a cone band and the radius at each end."""
    s = region.surface
    hat = _unit(s.axis)
    verts = m.vertices[np.unique(m.faces[region.faces])]
    along = verts[:, axis]
    height = (verts - np.asarray(s.apex, float)) @ hat
    rho = height * math.tan(float(s.half_angle))
    lo, hi = float(along.min()), float(along.max())
    band = max(0.05, 0.3 * (hi - lo))
    low = rho[along <= lo + band]
    high = rho[along >= hi - band]
    r_lo = float(np.median(low)) if len(low) else float(np.median(rho))
    r_hi = float(np.median(high)) if len(high) else float(np.median(rho))
    return lo, hi, r_lo, r_hi


def _attach_chamfer(m, feat, cones, used, axis, point, gap):
    """A 45 degree cone on this hole: a mouth flare, or a taper under a blind hole."""
    radius = float(feat["d"]) / 2.0
    best = None
    for region in cones:
        if id(region) in used:
            continue
        if _principal_axis(region.surface.axis) != axis:
            continue
        if _perp_dist(axis, point, region.surface.apex) > gap:
            continue
        lo, hi, r_lo, r_hi = _cone_ends(m, region, axis)
        axial = abs(hi - lo)
        if axial < max(0.05, 0.25 * gap):
            continue
        r_big, r_small = max(r_lo, r_hi), min(r_lo, r_hi)
        side = -1 if r_lo > r_hi else 1
        mouth = abs(r_small - radius) <= abs(r_big - radius) and r_big > radius + max(0.05, 0.2 * axial)
        taper = abs(r_big - radius) < abs(r_small - radius) and radius - r_small > max(0.05, 0.2 * axial)
        if mouth and abs(r_small - radius) <= gap:
            kind = "mouth"
        elif taper and abs(r_big - radius) <= gap and "depth" in feat:
            kind = "in"
        else:
            continue
        if best is None or axial > best[0]:
            best = (axial, kind, side, region)
    if best is None:
        return
    axial, kind, side, region = best
    used.add(id(region))
    feat["chamfer"] = float(axial)
    feat["chamfer_from"] = int(side)
    feat["chamfer_in"] = kind == "in"


def _feat_from_step(stepped, axis, point):
    feat = {"d": stepped["d"], "axis": _AXIS_NAME[axis]}
    for i, coord in enumerate("xyz"):
        if i != axis:
            feat[coord] = float(point[i])
    if "cbore_d" in stepped:
        feat["cbore_d"] = stepped["cbore_d"]
        feat["cbore_depth"] = stepped["cbore_depth"]
        feat["from"] = stepped["from"]
    if "depth" in stepped:
        feat["depth"] = stepped["depth"]
        feat["from"] = stepped["from"]
    return feat


def _holes(m, regions, tol):
    """Through-holes, counterbores, blind holes, and a 45 degree chamfer on them."""
    cyls = []
    cones = []
    for r in regions:
        s = r.surface
        if r.area < _min_area(m):
            continue
        if isinstance(s, Cone):
            if _principal_axis(s.axis) is not None and abs(float(s.half_angle) - math.pi / 4.0) <= math.radians(3.0):
                cones.append(r)
            continue
        if not isinstance(s, Cylinder):
            continue
        ax = _principal_axis(s.axis)
        if ax is None:
            continue
        span = _partial_span(m, s, r.faces)
        if span is not None and span < np.radians(250.0):
            continue
        if _outward(m, r, s):
            continue
        cyls.append(r)
    groups = []
    gtol = _eps(tol, 0.35)
    for r in cyls:
        s = r.surface
        ax = _principal_axis(s.axis)
        pt = np.asarray(s.point, float)
        for g in groups:
            if g["axis"] == ax and _perp_dist(ax, g["point"], pt) <= gtol:
                g["regs"].append(r)
                g["point"] = g["point"] + pt
                g["n"] += 1
                break
        else:
            groups.append({"axis": ax, "point": pt.copy(), "regs": [r], "n": 1})
    holes = []
    annuli = set()
    gap = _eps(tol)
    # A real break between two holes, not the zero gap where a counterbore steps.
    split_gap = max(1.0, 8.0 * tol)
    used_cones = set()
    for g in groups:
        ax = g["axis"]
        pt = g["point"] / g["n"]
        radii = [float(r.surface.radius) for r in g["regs"]]
        big = max(radii)
        for r in regions:
            if isinstance(r.surface, Plane) and r.area >= _min_area(m) and _annulus(m, r, ax, pt, big, tol):
                annuli.add(id(r))
        segs = []
        for r in g["regs"]:
            verts = m.vertices[np.unique(m.faces[r.faces])]
            along = verts[:, ax]
            rad = float(r.surface.radius)
            segs.append(Segment("cylinder", float(along.min()), float(along.max()), rad, rad))
        segs.sort(key=lambda s: s.za)
        exterior = []
        for r in regions:
            if not isinstance(r.surface, Plane) or r.area < _min_area(m):
                continue
            if _attitude(r.surface.normal, ax, 4.0) != "cap":
                continue
            if id(r) in annuli:
                continue
            exterior.append(_plane_coord(r.surface, ax))
        targets = exterior + _cap_coords([r for r in regions if id(r) in annuli], ax, set())
        for s in segs:
            s.za = _snap(s.za, targets, gap)
            s.zb = _snap(s.zb, targets, gap)
        for a, b in zip(segs, segs[1:]):
            if abs(a.zb - b.za) <= gap:
                mid = 0.5 * (a.zb + b.za)
                a.zb = b.za = mid
        span_lo, span_hi = _part_span(m, ax, exterior, tol)
        chains = _split_chains(segs, split_gap)
        chains.sort(key=lambda ch: -max(float(s.ra) for s in ch))
        emitted = []
        for chain in chains:
            if _covered_by_through(chain, emitted, ax, pt, gap):
                continue
            opened, lo, hi = _chain_opened(chain, exterior, gap)
            stepped = _plain_or_counterbore(chain, lo, hi, gap) if opened else None
            through = stepped is not None and "depth" not in (stepped or {})
            if stepped is None:
                r_small = min(float(s.ra) for s in chain)
                occupied = [(float(s.za), float(s.zb)) for s in chain]
                if span_hi > span_lo + gap and _span_clear(m, pt, ax, span_lo, span_hi, occupied):
                    rebuilt = _extend_small(chain, span_lo, span_hi, r_small, gap)
                    stepped = _plain_or_counterbore(rebuilt, span_lo, span_hi, gap)
                    through = stepped is not None
                    if through:
                        lo, hi = span_lo, span_hi
            if stepped is None:
                stepped = _blind_hole(chain, exterior, gap)
                through = False
            if stepped is None:
                if len(chain) == 1:
                    raise Unsupported("hole does not pass through")
                raise Unsupported("stepped hole is not a counterbored through-hole")
            feat = _feat_from_step(stepped, ax, pt)
            _attach_chamfer(m, feat, cones, used_cones, ax, pt, gap)
            holes.append(feat)
            if through:
                # the cutter runs the whole part, so a fragment on the far
                # side of a channel is the same hole and must not be repeated
                emitted.append({
                    "axis": ax,
                    "point": pt,
                    "radius": float(stepped["d"]) / 2.0,
                    "lo": float(span_lo),
                    "hi": float(span_hi),
                })
    holes.sort(key=lambda h: (h["axis"], h.get("x", 0.0), h.get("y", 0.0), h.get("z", 0.0)))
    return holes, annuli


def _clips(m, regions, tol):
    """Faces tilted to X, Y and Z. A handful become planar clips; more is a refusal."""
    found = []
    for r in regions:
        if not isinstance(r.surface, Plane) or r.area < _min_area(m):
            continue
        n = r.surface.normal
        if all(_attitude(n, ax, 4.0) == "tilted" for ax in range(3)):
            found.append(r)
    if len(found) > 4:
        found.sort(key=lambda r: -r.area)
        raise Unsupported(
            f"{len(found)} faces are tilted to X, Y and Z "
            f"(largest {found[0].area:.0f} mm²), more than a few planar clips")
    clips = []
    for r in sorted(found, key=lambda r: -r.area):
        n = _unit(r.surface.normal)
        d = float(r.surface.d) / max(float(np.linalg.norm(r.surface.normal)), 1e-12)
        clips.append({
            "x": float(f"{n[0]:.4f}"),
            "y": float(f"{n[1]:.4f}"),
            "z": float(f"{n[2]:.4f}"),
            "d": float(f"{d:.4f}"),
        })
    return clips, {id(r) for r in found}


def _make_face_2d(rings):
    """Planar face in z=0. ``rings[0]`` is the outer loop; the rest are holes."""
    def wire(pts, ccw):
        seq = list(pts)
        if (_shoelace(seq) > 0) != ccw:
            seq = list(reversed(seq))
        poly = BRepBuilderAPI_MakePolygon()
        for x, y in seq:
            poly.Add(gp_Pnt(float(x), float(y), 0.0))
        poly.Close()
        if not poly.IsDone():
            return None
        return poly.Wire()

    outer = wire(rings[0], True)
    if outer is None:
        return None
    mk = BRepBuilderAPI_MakeFace(outer)
    if not mk.IsDone():
        return None
    for inner in rings[1:]:
        w = wire(inner, False)
        if w is not None:
            mk.Add(w)
    if not mk.IsDone():
        mk = BRepBuilderAPI_MakeFace(outer)
    return mk.Face() if mk.IsDone() else None


def _region_rings(m, region, U, V, tol):
    loops = _loops(m, region.faces)
    rings = []
    for loop in loops:
        xy = _project(m.vertices[np.asarray(loop, int)], U, V)
        ring = _drop_near([(float(p[0]), float(p[1])) for p in xy], tol)
        if len(ring) >= 3 and abs(_shoelace(ring)) > 1e-3:
            rings.append(ring)
    rings.sort(key=lambda r: -abs(_shoelace(r)))
    return rings


def _fuse_2d(faces, fuzzy):
    if not faces:
        return None
    if len(faces) == 1:
        return faces[0]
    op = BOPAlgo_Builder()
    for face in faces:
        op.AddArgument(face)
    op.SetFuzzyValue(float(fuzzy))
    op.Perform()
    if op.HasErrors():
        return None
    sew = BRepBuilderAPI_Sewing(float(fuzzy))
    sew.Add(op.Shape())
    sew.Perform()
    sewn = sew.SewedShape()
    # UnifySameDomain segfaults on a large messy union (a fitted printed part
    # whose outline is hundreds of faces). A child process keeps the refusal
    # path alive; a small outline stays in-process, as the tests do.
    if len(faces) > 24:
        from .brep import safe_unify
        return safe_unify(sewn, float(fuzzy))
    uni = ShapeUpgrade_UnifySameDomain(sewn, True, True, False)
    uni.Build()
    return uni.Shape()


def _faces_of(shape):
    found = []
    if shape is None:
        return found
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        found.append(TopoDS.Face_s(exp.Current()))
        exp.Next()
    return found


def _largest_face(shape):
    best, best_a = None, -1.0
    for face in _faces_of(shape):
        props = GProp_GProps()
        BRepGProp.SurfaceProperties_s(face, props)
        area = float(props.Mass())
        if area > best_a:
            best, best_a = face, area
    return best


def _wire_points(wire, face):
    pts = []
    exp = BRepTools_WireExplorer(wire, face)
    while exp.More():
        edge = exp.Current()
        curve = BRepAdaptor_Curve(edge)
        if curve.GetType() == GeomAbs_Line:
            p = BRep_Tool.Pnt_s(exp.CurrentVertex())
            pts.append((p.X(), p.Y()))
        else:
            u0, u1 = curve.FirstParameter(), curve.LastParameter()
            for i in range(24):
                p = curve.Value(u0 + (u1 - u0) * i / 24.0)
                pts.append((p.X(), p.Y()))
        exp.Next()
    return _drop_near([(float(p[0]), float(p[1])) for p in pts], 1e-6)


def _face_loops(face):
    outer = BRepTools.OuterWire_s(face)
    out_pts = _wire_points(outer, face)
    inners = []
    exp = TopExp_Explorer(face, TopAbs_WIRE)
    while exp.More():
        wire = TopoDS.Wire_s(exp.Current())
        if not wire.IsSame(outer):
            pts = _wire_points(wire, face)
            if len(pts) >= 3:
                inners.append(pts)
        exp.Next()
    return out_pts, inners


def _as_circle(pts, tol):
    p = np.asarray(pts, float)
    if len(p) < 8:
        return None
    c = p.mean(axis=0)
    rho = np.linalg.norm(p - c, axis=1)
    radius = float(np.median(rho))
    if radius < 0.2:
        return None
    if float(np.max(np.abs(rho - radius))) > _eps(tol, 0.25):
        return None
    return (float(c[0]), float(c[1])), radius


def _hole_uv(hole, axis):
    U, V = _PROFILE_AXES[axis]
    pt = np.zeros(3)
    for coord, value in hole.items():
        if coord in "xyz":
            pt["xyz".index(coord)] = float(value)
    return (float(pt @ U), float(pt @ V))


def _is_hole_loop(pts, holes, axis, tol):
    got = _as_circle(pts, tol)
    if got is None:
        return False
    (u, v), radius = got
    for hole in holes:
        hu, hv = _hole_uv(hole, axis)
        if math.hypot(u - hu, v - hv) > _eps(tol, 0.5):
            continue
        for key in ("d", "cbore_d"):
            if key in hole and abs(radius - float(hole[key]) / 2.0) <= _eps(tol, 0.35):
                return True
    return False


def _snap_ring(pts, walls, U, V, tol):
    if len(pts) < 3:
        return None
    cloud = _Cloud([np.asarray(U) * p[0] + np.asarray(V) * p[1] for p in pts])
    got = _loop_polygon(cloud, list(range(len(pts))), walls, U, V, tol)
    if got is None or got[1] != "polygon":
        return None
    return got[0]


def _accept_ring(raw, snapped, tol):
    """Keep the wall-snapped ring when it is the same polygon, not a smaller one."""
    raw_c = _drop_collinear(list(raw), tol)
    if snapped is None:
        return raw_c
    raw_a = abs(_shoelace(raw_c))
    snap_a = abs(_shoelace(snapped))
    if raw_a < 1e-6:
        return snapped
    if snap_a < 0.97 * raw_a or snap_a > 1.08 * raw_a:
        return raw_c
    return snapped


def _walls_of(regions, axis, skip):
    walls = []
    for r in regions:
        if id(r) in skip or not isinstance(r.surface, Plane):
            continue
        if r.area <= 0:
            continue
        if _attitude(r.surface.normal, axis, 4.0) == "wall":
            walls.append(r.surface)
    return _merge_walls(walls, 0.05)


def _extent(m, axis, regions, skip, tol):
    coords = m.vertices[:, axis]
    lo = float(coords.min())
    hi = float(coords.max())
    caps = _cap_coords(regions, axis, skip)
    band = max(8.0 * tol, 0.4)

    def one(extreme, sign):
        if caps:
            c = min(caps, key=lambda q: abs(q - extreme))
            if abs(c - extreme) <= band:
                return float(c)
        return float(_cluster_end(coords, extreme, sign, tol))

    return one(lo, -1), one(hi, 1)


def _usable_face(face) -> bool:
    if face is None:
        return False
    try:
        return bool(BRepCheck_Analyzer(face).IsValid())
    except Exception:
        return False


def _marching_loops(mask, lo, pitch, eps):
    """Closed contours of a binary mask. Largest area first.

    The ring is split at its farthest point before simplification, because the
    closing chord of a loop has zero length and would otherwise collapse it.
    """
    grid = np.pad(mask.astype(np.uint8), 1)
    height, width = grid.shape
    edges = defaultdict(list)

    def add(p, q):
        edges[p].append(q)
        edges[q].append(p)

    def mid(x, y, edge):
        return ((x + 0.5, y + 1.0), (x + 1.0, y + 0.5), (x + 0.5, y), (x, y + 0.5))[edge]

    for y in range(height - 1):
        row = grid[y:y + 2]
        for x in range(width - 1):
            nw, ne = int(row[0, x]), int(row[0, x + 1])
            sw, se = int(row[1, x]), int(row[1, x + 1])
            pts = []
            if sw != se:
                pts.append(mid(x, y, 0))
            if se != ne:
                pts.append(mid(x, y, 1))
            if ne != nw:
                pts.append(mid(x, y, 2))
            if nw != sw:
                pts.append(mid(x, y, 3))
            if len(pts) == 2:
                add(*pts)
            elif len(pts) == 4:
                # a saddle: connect along the filled diagonal
                if nw == se:
                    add(mid(x, y, 2), mid(x, y, 3))
                    add(mid(x, y, 0), mid(x, y, 1))
                else:
                    add(mid(x, y, 3), mid(x, y, 0))
                    add(mid(x, y, 1), mid(x, y, 2))
    used = set()

    def key(a, b):
        return (a, b) if a <= b else (b, a)

    loops = []
    for start, nbs in edges.items():
        for nb in nbs:
            if key(start, nb) in used:
                continue
            loop = [start]
            cur = nb
            used.add(key(start, cur))
            while cur != start:
                loop.append(cur)
                nxts = [p for p in edges[cur] if key(cur, p) not in used]
                if not nxts:
                    break
                nxt = nxts[0]
                used.add(key(cur, nxt))
                cur = nxt
            if cur == start and len(loop) >= 4:
                loops.append(np.asarray(loop, float))
    out = []
    for raw in loops:
        world = np.stack([
            lo[0] + (raw[:, 0] - 1.0) * pitch,
            lo[1] + (raw[:, 1] - 1.0) * pitch,
        ], axis=1)
        ring = np.vstack([world, world[0]])

        def rec(a, b, keep):
            if b <= a + 1:
                return
            seg = ring[b] - ring[a]
            length = float(np.linalg.norm(seg))
            if length < 1e-9:
                return
            rel = ring[a + 1:b] - ring[a]
            dist = np.abs(rel[:, 0] * seg[1] - rel[:, 1] * seg[0]) / length
            i = int(np.argmax(dist))
            if dist[i] > eps:
                rec(a, a + 1 + i, keep)
                keep.append(a + 1 + i)
                rec(a + 1 + i, b, keep)

        rel = ring[1:-1] - ring[0]
        if len(rel) == 0:
            continue
        far = int(np.argmax(np.linalg.norm(rel, axis=1))) + 1
        keep = [0, far]
        rec(0, far, keep)
        rec(far, len(ring) - 1, keep)
        keep = sorted({k for k in keep if k < len(world)})
        if len(keep) >= 3:
            out.append([(float(p[0]), float(p[1])) for p in ring[keep]])
    out.sort(key=lambda loop: -abs(_shoelace(loop)))
    return out


def _raster_mask(mesh, face_ids, axis, pitch):
    """Pixels covered by ``face_ids`` projected along ``axis``. ``(u, v)`` matches the outline."""
    axes = [i for i in range(3) if i != axis]
    u, v = axes
    if face_ids is None:
        tris = mesh.triangles
        mins = mesh.vertices[:, [u, v]].min(0)
        maxs = mesh.vertices[:, [u, v]].max(0)
    else:
        tris = mesh.triangles[np.asarray(face_ids, int)]
        if len(tris) == 0:
            return None, None, None
        flat = tris.reshape(-1, 3)
        mins = flat[:, [u, v]].min(0)
        maxs = flat[:, [u, v]].max(0)
    nu = int(np.ceil((maxs[0] - mins[0]) / pitch)) + 3
    nv = int(np.ceil((maxs[1] - mins[1]) / pitch)) + 3
    if nu > 2500 or nv > 2500 or nu < 2 or nv < 2:
        return None, None, None
    mask = np.zeros((nv, nu), dtype=bool)
    lo = np.asarray(mins, float)
    for tri in tris:
        p = tri[:, [u, v]]
        b0 = np.clip(np.floor((p.min(0) - lo) / pitch).astype(int), 0, [nu - 1, nv - 1])
        b1 = np.clip(np.ceil((p.max(0) - lo) / pitch).astype(int), 0, [nu - 1, nv - 1])
        xs = np.arange(b0[0], b1[0] + 1)
        ys = np.arange(b0[1], b1[1] + 1)
        if xs.size == 0 or ys.size == 0:
            continue
        xx, yy = np.meshgrid(xs, ys)
        pts = np.stack([lo[0] + xx * pitch, lo[1] + yy * pitch], axis=-1)
        a, b, c = p
        v0, v1, v2 = c - a, b - a, pts - a
        den = v0[0] * v1[1] - v1[0] * v0[1]
        if abs(den) < 1e-14:
            continue
        uu = (v2[..., 0] * v1[1] - v1[0] * v2[..., 1]) / den
        vv = (v0[0] * v2[..., 1] - v2[..., 0] * v0[1]) / den
        hit = (uu >= -1e-6) & (vv >= -1e-6) & (uu + vv <= 1.0 + 1e-6)
        mask[yy[hit], xx[hit]] = True
    return mask, lo, (u, v)


def _raster_outline(mesh, axis, pitch, eps, fill=True, face_ids=None):
    got = _raster_mask(mesh, face_ids, axis, pitch)
    mask, lo, _axes = got
    if mask is None or not np.any(mask):
        return []
    if fill:
        mask = ndimage.binary_fill_holes(mask)
    return _marching_loops(mask, lo, pitch, eps)


def _simplify_ring(pts, lim):
    """Drop points that sit on the chord of their neighbours, down to ``lim``."""
    pts = [(float(p[0]), float(p[1])) for p in pts]
    if len(pts) <= 3:
        return pts
    changed = True
    while changed and len(pts) > 3:
        changed = False
        n = len(pts)
        for i in range(n):
            a = np.asarray(pts[(i - 1) % n], float)
            b = np.asarray(pts[i], float)
            c = np.asarray(pts[(i + 1) % n], float)
            ac = c - a
            norm = float(np.linalg.norm(ac))
            if norm < 1e-9:
                dist = float(np.linalg.norm(b - a))
            else:
                dist = abs(float(ac[0] * (a[1] - b[1]) - ac[1] * (a[0] - b[0]))) / norm
            if dist <= lim:
                del pts[i]
                changed = True
                break
    return pts


def _point_in_ring(x, y, loop):
    inside = False
    j = len(loop) - 1
    for i in range(len(loop)):
        xi, yi = loop[i]
        xj, yj = loop[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-30) + xi):
            inside = not inside
        j = i
    return inside


def _recess_is_void(mesh, region, axis, coord, exterior):
    """The slab between a floor and the face it opens on is empty in the mesh."""
    mid = 0.5 * (float(coord) + float(exterior))
    centers = mesh.triangles_center[region.faces]
    if len(centers) == 0:
        return False
    step = max(1, len(centers) // 40)
    pts = np.array(centers[::step], float, copy=True)
    pts[:, axis] = mid
    inside = _points_inside(mesh, pts)
    if len(inside) == 0:
        return False
    return float(np.mean(inside)) < 0.25


def _floor_ring(m, region, axis, tol):
    """Outline of one planar floor, outer loop first. Raster when the boundary is a scribble."""
    U, V = _PROFILE_AXES[axis]
    rings = _region_rings(m, region, U, V, tol)
    lim = max(tol, 0.04)

    def usable(loops):
        if not loops:
            return None
        outer = _simplify_ring(loops[0], lim)
        if len(outer) < 3 or abs(_shoelace(outer)) < 0.45 * float(region.area):
            return None
        return [outer, *[_simplify_ring(inner, lim) for inner in loops[1:]]]

    got = usable(rings)
    if got is not None:
        return got
    pitch = max(0.08, 1.2 * tol)
    eps = max(0.05, tol)
    return usable(_raster_outline(m, axis, pitch, eps, fill=False, face_ids=region.faces))


def _floor_cuts(m, regions, skip, tol):
    """Blind pockets whose floor is parallel to a larger face and opens onto it.

    A lettered recess is not a straight prism: the floor is flat, the walls
    are not. The outline is the floor itself, cut at constant depth.
    """
    caps = []
    for region in regions:
        if id(region) in skip or not isinstance(region.surface, Plane):
            continue
        if region.area < _min_area(m):
            continue
        normal = _unit(region.surface.normal)
        axis = int(np.argmax(np.abs(normal)))
        if _attitude(region.surface.normal, axis, 6.0) != "cap":
            continue
        caps.append((region, axis, normal, _plane_coord(region.surface, axis)))
    found = []
    for region, axis, normal, coord in caps:
        best = None
        for other, oaxis, onormal, ocoord in caps:
            if other is region or oaxis != axis or other.area <= region.area * 1.2:
                continue
            if float(normal @ onormal) < 0.98:
                continue
            if normal[axis] > 0.0 and ocoord <= coord:
                continue
            if normal[axis] < 0.0 and ocoord >= coord:
                continue
            depth = abs(float(ocoord) - float(coord))
            if depth < max(3.0 * tol, 0.08):
                continue
            if best is None or depth < best[0]:
                best = (depth, float(ocoord))
        if best is None:
            continue
        depth, exterior = best
        if not _recess_is_void(m, region, axis, coord, exterior):
            continue
        rings = _floor_ring(m, region, axis, tol)
        if not rings:
            continue
        outer = rings[0]
        # a floor that spans most of the part is a wall seen edge-on, not a pocket
        extent = float(m.vertices[:, axis].max() - m.vertices[:, axis].min())
        if extent <= 0 or depth > 0.45 * extent:
            continue
        others = [i for i in range(3) if i != axis]
        cross = 1.0
        for i in others:
            cross *= float(m.vertices[:, i].max() - m.vertices[:, i].min())
        if cross > 0 and abs(_shoelace(outer)) > 0.35 * cross:
            continue
        cutouts = []
        for inner in rings[1:]:
            if len(inner) < 3 or abs(_shoelace(inner)) >= abs(_shoelace(outer)):
                continue
            # the bore through a pocket floor is drilled as a hole, not kept as an island
            if _as_circle(inner, tol) is not None:
                continue
            cx = sum(p[0] for p in inner) / len(inner)
            cy = sum(p[1] for p in inner) / len(inner)
            if _point_in_ring(cx, cy, outer):
                cutouts.append(inner)
        found.append({
            "area": float(region.area),
            "spec": {
                "axis": _AXIS_NAME[axis],
                "axis_i": axis,
                "origin": float(min(coord, exterior)),
                "length": float(depth),
                "outline": outer,
                "cutouts": cutouts,
            },
        })
    found.sort(key=lambda item: -item["area"])
    return [item["spec"] for item in found]


def _silhouette(m, regions, axis, holes, skip, tol):
    """Outer profile of the part looking along ``axis``, plus polygonal openings."""
    U, V = _PROFILE_AXES[axis]
    # Wall planes snap a noisy boundary back to the straight outline. The faces
    # themselves are edge-on along this axis: a noisy sliver of one cuts the
    # union into pieces, and the largest piece is then only part of the profile.
    walls = _walls_of(regions, axis, skip)
    faces = []
    for r in regions:
        if id(r) in skip or not isinstance(r.surface, Plane):
            continue
        if r.area < _min_area(m):
            continue
        if _attitude(r.surface.normal, axis, 4.0) == "wall":
            continue
        # A corner clip stays in the union. Leaving it out bites a notch into
        # the outline, and the extrusion would carry that notch along the axis.
        rings = _region_rings(m, r, U, V, tol)
        if not rings:
            continue
        outer = _accept_ring(rings[0], _snap_ring(rings[0], walls, U, V, tol), tol)
        inners = [
            _accept_ring(inner, _snap_ring(inner, walls, U, V, tol), tol)
            for inner in rings[1:]
        ]
        face = _make_face_2d([outer, *inners])
        if not _usable_face(face):
            face = _make_face_2d([outer])
        if _usable_face(face):
            faces.append(face)
    if not faces:
        return None
    fuzzy = max(4.0 * tol, 0.15)
    fused_shape = _fuse_2d(faces, fuzzy)
    # A wider gap only when the faces did not sew into one. On a clean part the
    # wider gap bridges a real notch (a counterbore, a shallow step).
    if len(_faces_of(fused_shape)) != 1:
        wider = _fuse_2d(faces, max(fuzzy * 2.0, 0.3))
        if len(_faces_of(wider)) == 1 or fused_shape is None:
            fused_shape = wider
    fused = _largest_face(fused_shape)
    if fused is None:
        return None
    outer, inners = _face_loops(fused)
    if len(outer) < 3:
        return None
    walls = _walls_of(regions, axis, skip)
    outer = _accept_ring(outer, _snap_ring(outer, walls, U, V, tol), tol)
    outer = _nice_loop(_drop_collinear(outer, tol), tol)
    if len(outer) < 3 or abs(_shoelace(outer)) < 1e-3:
        return None
    cutouts = []
    for inner in inners:
        if _is_hole_loop(inner, holes, axis, tol):
            continue
        if _as_circle(inner, tol) is not None:
            continue
        ring = _accept_ring(inner, _snap_ring(inner, walls, U, V, tol), tol)
        ring = _nice_loop(_drop_collinear(ring, tol), tol)
        if len(ring) >= 3 and abs(_shoelace(ring)) > 1e-3:
            cutouts.append(ring)
    origin, end = _extent(m, axis, regions, skip, tol)
    length = end - origin
    if length <= _eps(tol, 0.2):
        return None
    return {
        "axis": _AXIS_NAME[axis],
        "axis_i": axis,
        "origin": float(origin),
        "length": float(length),
        "outline": outer,
        "cutouts": cutouts,
    }


def _profile_face(outline, cutouts, axis, origin):
    U, V = _PROFILE_AXES[axis]
    base = np.zeros(3)
    base[axis] = float(origin)
    # (U, V) is right-handed about +axis except for a Y extrusion, where U×V = -Y.
    flip = float(np.cross(U, V)[axis]) < 0

    def wire(pts, ccw_wanted):
        seq = [(float(p[0]), float(p[1])) for p in pts]
        if (_shoelace(seq) > 0) != ccw_wanted:
            seq = list(reversed(seq))
        poly = BRepBuilderAPI_MakePolygon()
        for u, v in seq:
            q = base + u * U + v * V
            poly.Add(gp_Pnt(float(q[0]), float(q[1]), float(q[2])))
        poly.Close()
        return poly.Wire()

    mk = BRepBuilderAPI_MakeFace(wire(outline, not flip))
    for cut in cutouts:
        mk.Add(wire(cut, flip))
    if not mk.IsDone():
        raise Unsupported("an extruded outline is not a closed polygon")
    return mk.Face()


def _prism(spec):
    face = _profile_face(spec["outline"], spec.get("cutouts") or [], spec["axis_i"], spec["origin"])
    vec = [0.0, 0.0, 0.0]
    vec[spec["axis_i"]] = float(spec["length"])
    return BRepPrimAPI_MakePrism(face, gp_Vec(*vec), True).Shape()


def _numeric_span(positives, axis):
    lows, highs = [], []
    for p in positives:
        if p["axis_i"] == axis:
            lows.append(float(p["origin"]))
            highs.append(float(p["origin"]) + float(p["length"]))
        else:
            idx = _OUTLINE_OF[p["axis_i"]][axis]
            vals = [float(pt[idx]) for pt in p["outline"]]
            lows.append(min(vals))
            highs.append(max(vals))
    return max(lows), min(highs)


def _hole_point(hole):
    pt = [0.0, 0.0, 0.0]
    for i, coord in enumerate("xyz"):
        if coord in hole:
            pt[i] = float(hole[coord])
    return pt


def _add_cut(shape, tool):
    op = BRepAlgoAPI_Cut(shape, tool)
    op.Build()
    return op.Shape()


def _cylinder(point, axis, radius, length):
    direction = [0.0, 0.0, 0.0]
    direction[axis] = 1.0
    ax = gp_Ax2(gp_Pnt(*[float(v) for v in point]), gp_Dir(*direction))
    return BRepPrimAPI_MakeCylinder(ax, float(radius), float(length)).Shape()


def _cone(point, axis, r1, r2, length):
    direction = [0.0, 0.0, 0.0]
    direction[axis] = 1.0
    ax = gp_Ax2(gp_Pnt(*[float(v) for v in point]), gp_Dir(*direction))
    return BRepPrimAPI_MakeCone(ax, max(float(r1), 1e-4), max(float(r2), 1e-4), float(length)).Shape()


def _hole_cutters(hole, lo, hi):
    """Cylinders and 45 degree cones that cut one hole out of the positive solid.

    ``lo`` and ``hi`` are the part extent along the hole axis. A through hole
    runs past both ends. A blind hole stops at ``depth`` from the mouth. A
    mouth chamfer widens the opening; a bottom chamfer tapers inward from the
    end of a blind cylinder. The script emitters use the same placement.
    """
    ax = "XYZ".index(hole["axis"])
    radius = float(hole["d"]) / 2.0
    tools = []
    if "depth" in hole:
        depth = float(hole["depth"])
        if hole.get("from", 1) < 0:
            start_c, length_c = lo - 1.0, depth + 1.0
            inner = lo + depth
        else:
            start_c, length_c = hi - depth, depth + 1.0
            inner = hi - depth
    else:
        start_c, length_c = lo - 1.0, (hi - lo) + 2.0
        inner = None
    origin = _hole_point(hole)
    origin[ax] = start_c
    tools.append((origin, ax, radius, radius, length_c))
    if "cbore_d" in hole:
        depth = float(hole["cbore_depth"])
        c0 = _hole_point(hole)
        if hole.get("from", 1) < 0:
            c0[ax] = lo - 1.0
            clen = depth + 1.0
        else:
            c0[ax] = hi - depth
            clen = depth + 1.0
        tools.append((c0, ax, float(hole["cbore_d"]) / 2.0, float(hole["cbore_d"]) / 2.0, clen))
    if "chamfer" in hole:
        h = float(hole["chamfer"])
        if hole.get("chamfer_in") and inner is not None:
            c0 = _hole_point(hole)
            if hole.get("from", 1) < 0:
                c0[ax] = inner
                tools.append((c0, ax, radius, max(radius - h, 1e-3), h))
            else:
                c0[ax] = inner - h
                tools.append((c0, ax, max(radius - h, 1e-3), radius, h))
        elif not hole.get("chamfer_in"):
            c0 = _hole_point(hole)
            if hole.get("chamfer_from", 1) < 0:
                c0[ax] = lo
                tools.append((c0, ax, radius + h, radius, h))
            else:
                c0[ax] = hi - h
                tools.append((c0, ax, radius, radius + h, h))
    return tools


def _clip_tool(clip):
    n = np.array([clip["x"], clip["y"], clip["z"]], float)
    nn = float(n @ n)
    if nn < 1e-12:
        raise Unsupported("a planar clip has a zero normal")
    origin = n * (float(clip["d"]) / nn)
    direction = n / math.sqrt(nn)
    helper = np.array([1.0, 0.0, 0.0]) if abs(direction[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    xdir = np.cross(direction, helper)
    xdir /= np.linalg.norm(xdir)
    ydir = np.cross(direction, xdir)
    size, depth = 400.0, 200.0
    poly = BRepBuilderAPI_MakePolygon()
    for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        q = origin + xdir * (sx * size / 2.0) + ydir * (sy * size / 2.0)
        poly.Add(gp_Pnt(float(q[0]), float(q[1]), float(q[2])))
    poly.Close()
    face = BRepBuilderAPI_MakeFace(poly.Wire()).Face()
    return BRepPrimAPI_MakePrism(face, gp_Vec(*[float(v) for v in direction * depth]), True).Shape()


def _solid(positives, cuts, holes, clips):
    if not positives:
        return None
    shape = _prism(positives[0])
    for spec in positives[1:]:
        op = BRepAlgoAPI_Common(shape, _prism(spec))
        op.Build()
        shape = op.Shape()
    for spec in cuts:
        # the same overlap the positive solids use, so a pocket that opens on
        # an outer face is not left with a skin where those solids were grown
        grown = dict(spec)
        grown["origin"] = float(spec["origin"]) - _OVERLAP
        grown["length"] = float(spec["length"]) + 2.0 * _OVERLAP
        shape = _add_cut(shape, _prism(grown))
    for hole in holes:
        ax = "XYZ".index(hole["axis"])
        lo, hi = _numeric_span(positives, ax)
        if hi <= lo:
            continue
        for start, _ax, r1, r2, length in _hole_cutters(hole, lo, hi):
            if abs(r1 - r2) <= 1e-9:
                tool = _cylinder(start, ax, r1, length)
            else:
                tool = _cone(start, ax, r1, r2, length)
            shape = _add_cut(shape, tool)
    for clip in clips:
        shape = _add_cut(shape, _clip_tool(clip))
    return shape


def _one_solid(shape):
    if shape is None:
        return None
    exp = TopExp_Explorer(shape, TopAbs_SOLID)
    found = []
    while exp.More():
        found.append(TopoDS.Solid_s(exp.Current()))
        exp.Next()
    if len(found) == 1:
        return found[0]
    return shape


def _score(shape, mesh, tol, samples):
    empty = {"max_deviation": 1e9, "volume_error_pct": 1e9}
    if shape is None:
        return empty
    try:
        out = shape_to_mesh(shape, 0.02)
    except Exception:
        return empty
    if out is None or len(getattr(out, "faces", [])) == 0:
        return empty
    dev = deviation(mesh, out, n=samples)
    vm = abs(float(mesh.volume))
    vo = volume(shape)
    vol_pct = 100.0 * abs(vo - vm) / vm if vm else 1e9
    return {"max_deviation": float(dev["hausdorff"]), "volume_error_pct": float(vol_pct)}


def _ok(score, tol) -> bool:
    if score is None:
        return False
    lim = max(4.0 * tol, 0.1)
    return score["max_deviation"] <= lim and score["volume_error_pct"] <= 1.0


def _fail(score, tol):
    lim = max(4.0 * tol, 0.1)
    raise Unsupported(
        f"combined extrusion deviates by {score['max_deviation']:.3f} mm "
        f"(limit {lim:.3f}), volume error {score['volume_error_pct']:.2f}%")


def _components(m, regs):
    parent = list(range(len(regs)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    edge_owner = {}
    for i, r in enumerate(regs):
        for fi in r.faces:
            tri = m.faces[int(fi)]
            for a, b in ((int(tri[0]), int(tri[1])), (int(tri[1]), int(tri[2])),
                         (int(tri[2]), int(tri[0]))):
                e = (a, b) if a < b else (b, a)
                j = edge_owner.get(e)
                if j is None:
                    edge_owner[e] = i
                else:
                    parent[find(i)] = find(j)
    groups = {}
    for i, r in enumerate(regs):
        groups.setdefault(find(i), []).append(r)
    return list(groups.values())


def _sample_out(m, region, delta):
    n = _unit(region.surface.normal)
    d = float(np.asarray(region.surface.normal, float) @ n)
    # Plane.d is n·x with the stored normal; rescale onto the unit normal.
    plane_d = float(region.surface.d) / max(d, 1e-12)
    faces = region.faces
    w = m.area_faces[faces]
    center = (m.triangles_center[faces] * w[:, None]).sum(0) / max(float(w.sum()), 1e-12)
    center = center - n * (float(center @ n) - plane_d)
    return center + n * delta


def _pocket_faces(m, regions, solid, skip, tol):
    if solid is None:
        return []
    delta = max(2.0 * tol, 0.15)
    clf = BRepClass3d_SolidClassifier(_one_solid(solid))
    found = []
    for r in regions:
        if id(r) in skip or not isinstance(r.surface, Plane) or r.area < _min_area(m):
            continue
        pt = _sample_out(m, r, delta)
        clf.Perform(gp_Pnt(float(pt[0]), float(pt[1]), float(pt[2])), 1e-6)
        if clf.State() == TopAbs_IN:
            found.append(r)
    return found


def _cut_from_component(m, comp, regions, skip, tol):
    """A pocket that is one straight extrusion. Longest axis wins when several fit.

    Two caps are a through cut. One cap is a blind pocket: the open end is the
    far side of the component, snapped to an exterior face when the pocket
    actually breaks that face. A lone plane is not a pocket. An axis whose
    caps do not make a polygon is skipped so a longer neighbour can still fit.
    """
    best = None
    gap = _eps(tol, 0.25)
    for axis in range(3):
        if any(_attitude(r.surface.normal, axis, 4.0) == "tilted" for r in comp):
            continue
        if not any(_attitude(r.surface.normal, axis, 4.0) == "wall" for r in comp):
            continue
        caps = []
        for r in comp:
            if _attitude(r.surface.normal, axis, 4.0) == "cap":
                caps.append((_plane_coord(r.surface, axis), r))
        caps.sort(key=lambda t: t[0])
        clusters = []
        for coord, region in caps:
            if not clusters or coord - clusters[-1][-1][0] > gap:
                clusters.append([(coord, region)])
            else:
                clusters[-1].append((coord, region))
        if len(clusters) not in (1, 2):
            continue
        verts = np.vstack([m.vertices[np.unique(m.faces[r.faces])] for r in comp])
        along = verts[:, axis]
        extent = float(along.max() - along.min())
        if len(clusters) == 2:
            origin = float(np.average([c for c, _r in clusters[0]]))
            end = float(np.average([c for c, _r in clusters[-1]]))
        else:
            cap_c = float(np.average([c for c, _r in clusters[0]]))
            exteriors = _cap_coords(regions, axis, skip)
            band = max(8.0 * tol, 0.4)
            if abs(cap_c - float(along.min())) <= abs(cap_c - float(along.max())):
                origin, end = cap_c, float(along.max())
                open_is_end = True
            else:
                origin, end = float(along.min()), cap_c
                open_is_end = False
            if exteriors:
                if open_is_end:
                    snapped = min(exteriors, key=lambda q: abs(q - end))
                    if abs(snapped - end) <= band:
                        end = float(snapped)
                else:
                    snapped = min(exteriors, key=lambda q: abs(q - origin))
                    if abs(snapped - origin) <= band:
                        origin = float(snapped)
        length = end - origin
        if length <= gap:
            continue
        U, V = _PROFILE_AXES[axis]
        walls = _walls_of(regions, axis, skip)
        loops = []
        for cluster in clusters:
            ids = []
            for _c, region in cluster:
                ids.extend(int(i) for i in region.faces)
            built = _polygons_from_cap(m, ids, walls, U, V, tol)
            if built and built[0][1] == "polygon":
                loops.append(built[0][2])
        if not loops:
            continue
        outline = _nice_loop(_drop_collinear(loops[0], tol), tol)
        if len(outline) < 3 or abs(_shoelace(outline)) < 1e-3:
            continue
        # the first axis wins a near tie, matching the previous longest-axis rule
        if best is None or extent > best[0] + 0.05:
            best = (extent, axis, origin, length, outline)
    if best is None:
        area = max(r.area for r in comp)
        raise Unsupported(
            f"a recess is not a straight cut along X, Y, or Z "
            f"({len(comp)} faces, largest {area:.0f} mm²)")
    _extent_len, axis, origin, length, outline = best
    return {
        "axis": _AXIS_NAME[axis],
        "axis_i": axis,
        "origin": float(origin),
        "length": float(length),
        "outline": outline,
        "cutouts": [],
    }


def _round_spec(spec, tol):
    spec = dict(spec)
    spec["origin"] = _nice(spec["origin"], tol)
    spec["length"] = _nice(spec["length"], tol)
    spec["outline"] = _nice_loop(spec["outline"], tol)
    spec["cutouts"] = [_nice_loop(c, tol) for c in spec.get("cutouts") or []]
    return spec


def _round_hole(hole, tol):
    out = dict(hole)
    out["d"] = _nice(out["d"], tol)
    for coord in "xyz":
        if coord in out:
            out[coord] = _nice(out[coord], tol)
    if "cbore_d" in out:
        out["cbore_d"] = _nice(out["cbore_d"], tol)
        out["cbore_depth"] = _nice(out["cbore_depth"], tol)
    if "depth" in out:
        out["depth"] = _nice(out["depth"], tol)
        out["from"] = 1 if out.get("from", 1) > 0 else -1
    if "chamfer" in out:
        out["chamfer"] = _nice(out["chamfer"], tol)
        out["chamfer_from"] = 1 if out.get("chamfer_from", 1) > 0 else -1
        out["chamfer_in"] = bool(out.get("chamfer_in"))
    return out


def _public(spec):
    return {k: spec[k] for k in ("axis", "origin", "length", "outline", "cutouts")}


def _with_index(spec):
    out = dict(spec)
    out["axis_i"] = "XYZ".index(spec["axis"])
    return out


def extract_csg(m: trimesh.Trimesh, regions, tol: float) -> dict:
    """Recognise a small boolean of extrusions. Raises :class:`Unsupported`."""
    _refuse_surfaces(m, regions)
    holes, annuli = _holes(m, regions, tol)
    clips, clip_ids = _clips(m, regions, tol)
    # Annulus faces would fill a counterbore. Clip faces stay: they are part of
    # the silhouette, and the half-space cut removes only the corner.
    skip = set(annuli)
    specs = {}
    for axis in range(3):
        sil = _silhouette(m, regions, axis, holes, skip, tol)
        if sil is not None:
            specs[axis] = _round_spec(sil, tol)
    if not specs:
        raise Unsupported("not a combined extrusion: no straight outline along X, Y, or Z")
    holes = [_round_hole(h, tol) for h in holes]
    samples = 4000

    def evaluate(axes, cuts, use_clips):
        positives = [specs[a] for a in axes if a in specs]
        try:
            shape = _solid(positives, cuts, holes, clips if use_clips else [])
        except Exception:
            return {"max_deviation": 1e9, "volume_error_pct": 1e9}
        return _score(shape, m, tol, samples)

    axes = [a for a in (2, 0, 1) if a in specs]
    cache = {}

    def score_of(key_axes, cuts, use_clips):
        key = (tuple(sorted(key_axes)), len(cuts), bool(use_clips))
        if key not in cache:
            cache[key] = evaluate(list(key_axes), cuts, use_clips)
        return cache[key]

    def shrink(chosen, cuts, use_clips):
        for drop in _DROP_ORDER:
            if drop not in chosen or len(chosen) == 1:
                continue
            trial = [a for a in chosen if a != drop]
            if _ok(score_of(trial, cuts, use_clips), tol):
                chosen = trial
        return chosen

    chosen = list(axes)
    use_clips = bool(clips)
    cuts = []
    best = score_of(chosen, cuts, use_clips)
    if not _ok(best, tol) and clips:
        bare = score_of(chosen, [], False)
        if _ok(bare, tol) or bare["max_deviation"] < best["max_deviation"]:
            use_clips = False
            best = bare
    if not _ok(best, tol) and _analytic_mesh(m):
        # Fusing the projected faces can notch a silhouette (a screw hole whose
        # wall sticks out of the union). A raster of the mesh is the true
        # outline. Only replace an outline the raster actually contains, and
        # only after the analytic solid has already missed the gate, so a part
        # whose fused outline is already right never changes.
        changed = False
        pitch = max(0.12, 2.0 * tol)
        eps = max(0.08, tol)
        for axis, spec in specs.items():
            loops = _raster_outline(m, axis, pitch, eps, fill=True)
            if not loops:
                continue
            if abs(_shoelace(loops[0])) <= abs(_shoelace(spec["outline"])) * 1.03:
                continue
            spec["outline"] = _nice_loop(_simplify_ring(loops[0], max(tol, 0.02)), tol)
            spec["cutouts"] = []
            changed = True
        if changed:
            cache.clear()
            best = score_of(chosen, cuts, use_clips)
    if not _ok(best, tol):
        # A silhouette is a superset. Extra material left inside it is a cut
        # when that material's faces form a straight extrusion. One pocket
        # that is not a straight cut does not throw the pockets that are.
        host = _solid([specs[a] for a in chosen], [], holes, clips if use_clips else [])
        built = []
        try:
            pockets = _pocket_faces(m, regions, host, skip, tol)
        except Exception:
            pockets = []
        for comp in _components(m, pockets):
            if sum(r.area for r in comp) < _min_area(m):
                continue
            try:
                built.append(_round_spec(_cut_from_component(m, comp, regions, skip, tol), tol))
            except Unsupported:
                continue
        if built:
            cache.clear()
            got = evaluate(chosen, built, use_clips)
            if _ok(got, tol) or got["max_deviation"] < best["max_deviation"]:
                cuts = built
                best = got
    if not _ok(best, tol) and _analytic_mesh(m):
        # A shallow floor with curved walls (a logo) is still a constant-depth
        # cut. Keep a floor only when it moves the solid closer to the mesh.
        for spec in _floor_cuts(m, regions, skip, tol):
            trial = cuts + [_round_spec(spec, tol)]
            cache.clear()
            got = evaluate(chosen, trial, use_clips)
            closer = got["max_deviation"] < best["max_deviation"] - 0.01
            same_dev = got["max_deviation"] <= best["max_deviation"] + 0.02
            if closer or (same_dev and got["volume_error_pct"] < best["volume_error_pct"] - 0.05):
                cuts = trial
                best = got
                if _ok(best, tol):
                    break
    if not _ok(best, tol):
        # one outline can be slightly tight and cut off real material; drop it
        for drop in _DROP_ORDER:
            if drop not in chosen or len(chosen) == 1:
                continue
            trial = [a for a in chosen if a != drop]
            cache.clear()
            got = evaluate(trial, cuts, use_clips)
            if not _ok(got, tol):
                got = evaluate(trial, [], use_clips)
                if _ok(got, tol):
                    chosen, cuts, best = trial, [], got
                    break
                continue
            chosen, best = trial, got
            break
    if not _ok(best, tol):
        _fail(best, tol)
    chosen = shrink(chosen, cuts, use_clips)
    # confirm the reduced solid, not a lucky sample of a larger one
    final = score_of(chosen, cuts, use_clips)
    if not _ok(final, tol):
        _fail(final, tol)
    positives = [_public(specs[a]) for a in (2, 0, 1) if a in chosen]
    return {
        "csg": {
            "positive": positives,
            "cuts": [_public(c) for c in cuts],
            "clips": clips if use_clips else [],
        },
        "holes": holes,
    }


def _names(positives):
    solo = len(positives) == 1
    named = []
    for p in positives:
        if solo:
            named.append({"len": "extrude_len", "origin": "extrude_0", "outline": "outline",
                          "cutout": "cutout"})
        else:
            a = p["axis"].lower()
            named.append({"len": f"extrude_{a}_len", "origin": f"extrude_{a}_0",
                          "outline": f"outline_{a}", "cutout": f"outline_{a}_cutout"})
    return named


def _span_expr(positives, names, hole_axis, lang):
    lows, highs = [], []
    for p, nm in zip(positives, names):
        ax = "XYZ".index(p["axis"])
        if ax == hole_axis:
            lows.append(nm["origin"])
            highs.append(f"({nm['origin']} + {nm['len']})")
        else:
            idx = _OUTLINE_OF[ax][hole_axis]
            ol = nm["outline"]
            if lang == "cadquery":
                lows.append(f"(min(p[{idx}] for p in {ol}))")
                highs.append(f"(max(p[{idx}] for p in {ol}))")
            else:
                lows.append(f"min([for (p = {ol}) p[{idx}]])")
                highs.append(f"max([for (p = {ol}) p[{idx}]])")
    if len(lows) == 1:
        return lows[0], highs[0]
    return f"max({', '.join(lows)})", f"min({', '.join(highs)})"


def _fmt_py(pts) -> str:
    rows = ["["]
    for x, y in pts:
        rows.append(f"    ({x!r}, {y!r}),")
    rows.append("]")
    return "\n".join(rows)


def _fmt_scad(pts) -> str:
    rows = ["["]
    for x, y in pts:
        rows.append(f"  [{x!r}, {y!r}],")
    rows.append("]")
    return "\n".join(rows)


def _indent(text: str, levels: int) -> str:
    pad = "  " * levels
    return "\n".join(pad + ln if ln else ln for ln in text.splitlines())


def _param_rows(feats, tol):
    """Float parameters, then point-list parameters, in script order."""
    csg = feats["csg"]
    positives = csg["positive"]
    names = _names(positives)
    numbers = []
    lists = []
    for p, nm in zip(positives, names):
        plane, coords = _PLANE["XYZ".index(p["axis"])]
        numbers.append((nm["len"], p["length"], f"length along {p['axis']}"))
        numbers.append((nm["origin"], p["origin"], f"start along {p['axis']}"))
        lists.append((nm["outline"], p["outline"],
                      f"outline in the {plane} plane, ({coords})"))
        for i, loop in enumerate(p.get("cutouts") or [], 1):
            lists.append((f"{nm['cutout']}_{i}", loop, "through opening in the outline"))
    for i, cut in enumerate(csg["cuts"], 1):
        plane, coords = _PLANE["XYZ".index(cut["axis"])]
        numbers.append((f"cut_{i}_len", cut["length"], f"length along {cut['axis']}"))
        numbers.append((f"cut_{i}_0", cut["origin"], f"start along {cut['axis']}"))
        lists.append((f"cut_{i}", cut["outline"],
                      f"cut profile in the {plane} plane, ({coords})"))
        for j, loop in enumerate(cut.get("cutouts") or [], 1):
            lists.append((f"cut_{i}_cutout_{j}", loop, "island left inside the cut"))
    for i, hole in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        numbers.append((f"{n}_d", hole["d"], f"diameter (axis {hole['axis']})"))
        for coord in "xyz":
            if coord in hole:
                numbers.append((f"{n}_{coord}", hole[coord], ""))
        if "depth" in hole:
            end = f"{'+' if hole.get('from', 1) > 0 else '-'}{hole['axis']}"
            numbers.append((f"{n}_depth", hole["depth"], f"blind depth from the {end} end"))
        if "cbore_d" in hole:
            end = f"{'+' if hole.get('from', 1) > 0 else '-'}{hole['axis']}"
            numbers.append((f"{n}_cbore_d", hole["cbore_d"], "counterbore diameter"))
            numbers.append((f"{n}_cbore_depth", hole["cbore_depth"],
                            f"counterbore depth from the {end} end"))
        if "chamfer" in hole:
            if hole.get("chamfer_in"):
                note = "45 degree taper at the inner end of the blind hole"
            else:
                end = f"{'+' if hole.get('chamfer_from', 1) > 0 else '-'}{hole['axis']}"
                note = f"45 degree chamfer from the {end} end"
            numbers.append((f"{n}_chamfer", hole["chamfer"], note))
    for i, clip in enumerate(csg["clips"], 1):
        n = f"clip_{i}"
        numbers.append((f"{n}_x", clip["x"], "outward normal"))
        numbers.append((f"{n}_y", clip["y"], ""))
        numbers.append((f"{n}_z", clip["z"], ""))
        numbers.append((f"{n}_d", clip["d"], "plane: normal · point = d"))
    return numbers, lists


def _workplane(axis: str, origin: str, length: str) -> str:
    if axis == "Z":
        return f'cq.Workplane("XY", origin=(0, 0, {origin}))'
    if axis == "Y":
        return (f"cq.Workplane(cq.Plane(origin=cq.Vector(0, {origin} + {length}, 0),"
                f" xDir=cq.Vector(1, 0, 0), normal=cq.Vector(0, -1, 0)))")
    return (f"cq.Workplane(cq.Plane(origin=cq.Vector({origin}, 0, 0),"
            f" xDir=cq.Vector(0, 1, 0), normal=cq.Vector(1, 0, 0)))")


def _cq_extrude(var, axis, origin, length, outline, cutouts, cutout_prefix, overlap=False):
    if overlap:
        origin = f"({origin} - {_OVERLAP})"
        length = f"({length} + {2 * _OVERLAP})"
    wp = _workplane(axis, origin, length)
    lines = [f"{var} = ({wp}",
             f"         .polyline({outline}).close().extrude({length}))"]
    for i in range(1, len(cutouts) + 1):
        lines.append(f"{var} = {var}.cut({wp}")
        lines.append(f"        .polyline({cutout_prefix}_{i}).close().extrude({length}))")
    return lines


def _cq_hole(n, hole, lo, hi):
    """CadQuery cuts for one hole. Placement matches :func:`_hole_cutters`."""
    lines = []
    direction = cadquery_axis_dir(hole["axis"])
    if "depth" in hole:
        if hole.get("from", 1) < 0:
            start, length = f"({lo} - 1)", f"{n}_depth + 1"
        else:
            start, length = f"(({hi}) - {n}_depth)", f"{n}_depth + 1"
    else:
        start, length = f"({lo} - 1)", f"(({hi}) - ({lo}) + 2)"
    lines.append(
        f"result = result.cut(cq.Workplane().add(cq.Solid.makeCylinder("
        f"{n}_d / 2, {length}, cq.Vector({axis_point_expr(hole['axis'], n, start)}), "
        f"cq.Vector({direction}))))")
    if "cbore_d" in hole:
        if hole.get("from", 1) < 0:
            c0, clen = f"({lo} - 1)", f"{n}_cbore_depth + 1"
        else:
            c0, clen = f"(({hi}) - {n}_cbore_depth)", f"{n}_cbore_depth + 1"
        lines.append(
            f"result = result.cut(cq.Workplane().add(cq.Solid.makeCylinder("
            f"{n}_cbore_d / 2, {clen}, cq.Vector({axis_point_expr(hole['axis'], n, c0)}), "
            f"cq.Vector({direction}))))")
    if "chamfer" not in hole:
        return lines
    if hole.get("chamfer_in") and "depth" in hole:
        if hole.get("from", 1) < 0:
            c0 = f"(({lo}) + {n}_depth)"
            r1, r2 = f"{n}_d / 2", f"{n}_d / 2 - {n}_chamfer"
        else:
            c0 = f"(({hi}) - {n}_depth - {n}_chamfer)"
            r1, r2 = f"{n}_d / 2 - {n}_chamfer", f"{n}_d / 2"
        clen = f"{n}_chamfer"
    elif hole.get("chamfer_from", 1) < 0:
        c0 = f"({lo})"
        r1, r2 = f"{n}_d / 2 + {n}_chamfer", f"{n}_d / 2"
        clen = f"{n}_chamfer"
    else:
        c0 = f"(({hi}) - {n}_chamfer)"
        r1, r2 = f"{n}_d / 2", f"{n}_d / 2 + {n}_chamfer"
        clen = f"{n}_chamfer"
    lines.append(
        f"result = result.cut(cq.Workplane().add(cq.Solid.makeCone("
        f"{r1}, {r2}, {clen}, cq.Vector({axis_point_expr(hole['axis'], n, c0)}), "
        f"cq.Vector({direction}))))")
    return lines


def render_cadquery(feats: dict, tol: float = 0.01, source: str | None = None) -> str:
    """CadQuery script for a boolean of extrusions."""
    del tol  # values are already rounded to the fit tolerance
    csg = feats["csg"]
    positives = csg["positive"]
    names = _names(positives)
    numbers, lists = _param_rows(feats, 0.01)
    body = []
    overlap = len(positives) > 1
    vars_ = []
    for p, nm in zip(positives, names):
        var = "result" if len(positives) == 1 else f"solid_{p['axis'].lower()}"
        vars_.append(var)
        body += _cq_extrude(var, p["axis"], nm["origin"], nm["len"], nm["outline"],
                            p.get("cutouts") or [], nm["cutout"], overlap=overlap)
    if len(vars_) > 1:
        body.append(f"result = {vars_[0]}")
        for var in vars_[1:]:
            body.append(f"result = result.intersect({var})")
    for i, cut in enumerate(csg["cuts"], 1):
        body += _cq_extrude(f"cut_{i}_solid", cut["axis"], f"cut_{i}_0", f"cut_{i}_len",
                            f"cut_{i}", cut.get("cutouts") or [], f"cut_{i}_cutout", overlap=True)
        body.append(f"result = result.cut(cut_{i}_solid)")
    for i, hole in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        ax = "XYZ".index(hole["axis"])
        lo, hi = _span_expr(positives, names, ax, "cadquery")
        body.extend(_cq_hole(n, hole, lo, hi))
    for i, _clip in enumerate(csg["clips"], 1):
        n = f"clip_{i}"
        body.append(
            f"_n{i} = cq.Vector({n}_x, {n}_y, {n}_z)\n"
            f"_nn{i} = {n}_x*{n}_x + {n}_y*{n}_y + {n}_z*{n}_z\n"
            f"_o{i} = _n{i} * ({n}_d / _nn{i})\n"
            f"_t{i} = cq.Vector(0, -{n}_z, {n}_y) if abs({n}_x) < 0.9 else cq.Vector(-{n}_z, 0, {n}_x)\n"
            f"result = result.cut(cq.Workplane(cq.Plane(origin=_o{i}, xDir=_t{i}, normal=_n{i}))"
            f".rect(400, 400).extrude(200))")
    width = max([len(n) for n, _, _ in numbers] + [len(n) for n, _, _ in lists] + [1])
    head = ['"""Parametric CadQuery model generated by mesh2cad.',
            "",
            f"Source: {source}" if source else "Edit the parameters below and re-run.",
            "Units are the mesh units (usually mm). Run: python this_file.py",
            '"""',
            "import math", "", "import cadquery as cq", "", "# --- parameters ---"]
    params = [f"{n.ljust(width)} = {v!r}" + (f"  # {note}" if note else "") for n, v, note in numbers]
    for name, pts, note in lists:
        if note:
            params.append(f"# {note}")
        params.append(f"{name} = {_fmt_py(pts)}")
    tail = ["", "if __name__ == \"__main__\":",
            "    from pathlib import Path",
            "    cq.exporters.export(result, str(Path(__file__).with_suffix(\".step\")))",
            "    print(\"wrote\", Path(__file__).with_suffix(\".step\"))", ""]
    return "\n".join(head + params + ["", "# --- model ---"] + body + tail)


def _scad_placed(axis, origin, length, polygon_lines, overlap=False) -> str:
    if overlap:
        origin = f"({origin} - {_OVERLAP})"
        length = f"({length} + {2 * _OVERLAP})"
    body = "\n".join(polygon_lines)
    if axis == "Z":
        return f"translate([0, 0, {origin}])\n  linear_extrude({length})\n" + _indent(body, 2)
    if axis == "Y":
        return (f"translate([0, {origin} + {length}, 0])\n"
                f"  rotate([90, 0, 0])\n"
                f"    linear_extrude({length})\n" + _indent(body, 3))
    return (f"translate([{origin}, 0, 0])\n"
            f"  rotate([90, 0, 90])\n"
            f"    linear_extrude({length})\n" + _indent(body, 3))


def _scad_profile(outline, cutouts, cutout_prefix) -> list[str]:
    if not cutouts:
        return [f"polygon({outline});"]
    lines = ["difference() {", f"  polygon({outline});"]
    for i in range(1, len(cutouts) + 1):
        lines.append(f"  polygon({cutout_prefix}_{i});")
    lines.append("}")
    return lines


def _scad_hole(prefix, hole, diam, start, length) -> str:
    point = axis_point_expr(hole["axis"], prefix, start)
    return (f"translate([{point}])\n"
            f"  rotate({openscad_axis_rotate(hole['axis'])})\n"
            f"    cylinder(h={length}, d={diam});")


def _scad_cone(prefix, hole, r1, r2, start, length) -> str:
    point = axis_point_expr(hole["axis"], prefix, start)
    return (f"translate([{point}])\n"
            f"  rotate({openscad_axis_rotate(hole['axis'])})\n"
            f"    cylinder(h={length}, r1={r1}, r2={r2});")


def _scad_hole_pieces(n, hole, lo, hi):
    """OpenSCAD cuts for one hole. Placement matches :func:`_hole_cutters`."""
    pieces = []
    if "depth" in hole:
        if hole.get("from", 1) < 0:
            start, length = f"({lo} - 1)", f"{n}_depth + 1"
        else:
            start, length = f"(({hi}) - {n}_depth)", f"{n}_depth + 1"
    else:
        start, length = f"({lo} - 1)", f"(({hi}) - ({lo}) + 2)"
    pieces.append(_scad_hole(n, hole, f"{n}_d", start, length))
    if "cbore_d" in hole:
        if hole.get("from", 1) < 0:
            c0, clen = f"({lo} - 1)", f"{n}_cbore_depth + 1"
        else:
            c0, clen = f"(({hi}) - {n}_cbore_depth)", f"{n}_cbore_depth + 1"
        pieces.append(_scad_hole(n, hole, f"{n}_cbore_d", c0, clen))
    if "chamfer" not in hole:
        return pieces
    if hole.get("chamfer_in") and "depth" in hole:
        if hole.get("from", 1) < 0:
            c0 = f"(({lo}) + {n}_depth)"
            r1, r2 = f"{n}_d / 2", f"{n}_d / 2 - {n}_chamfer"
        else:
            c0 = f"(({hi}) - {n}_depth - {n}_chamfer)"
            r1, r2 = f"{n}_d / 2 - {n}_chamfer", f"{n}_d / 2"
    elif hole.get("chamfer_from", 1) < 0:
        c0 = f"({lo})"
        r1, r2 = f"{n}_d / 2 + {n}_chamfer", f"{n}_d / 2"
    else:
        c0 = f"(({hi}) - {n}_chamfer)"
        r1, r2 = f"{n}_d / 2", f"{n}_d / 2 + {n}_chamfer"
    pieces.append(_scad_cone(n, hole, r1, r2, c0, f"{n}_chamfer"))
    return pieces


def _scad_clip(i) -> str:
    n = f"clip_{i}"
    return (
        f"translate([{n}_x, {n}_y, {n}_z] * {n}_d / ({n}_x*{n}_x + {n}_y*{n}_y + {n}_z*{n}_z))\n"
        f"  rotate(acos(max(-1, min(1, {n}_z / sqrt({n}_x*{n}_x + {n}_y*{n}_y + {n}_z*{n}_z)))),"
        f" [-{n}_y, {n}_x, 0])\n"
        f"    translate([-200, -200, 0])\n"
        f"      cube([400, 400, 200]);"
    )


def _hole_slider_limit(positives, hole) -> float:
    ax = "XYZ".index(hole["axis"])
    spans = []
    for i in range(3):
        if i == ax:
            continue
        lo, hi = _numeric_span([_with_index(p) for p in positives], i)
        spans.append(hi - lo)
    return max(min(spans), float(hole["d"]))


def render_openscad(feats: dict, tol: float = 0.01, source: str | None = None) -> str:
    """OpenSCAD script for a boolean of extrusions."""
    del tol
    from .openscad_code import CIRCLE_SEGMENTS, _slider

    csg = feats["csg"]
    positives = csg["positive"]
    names = _names(positives)
    numbers, lists = _param_rows(feats, 0.01)
    by_name = {n: v for n, v, _ in numbers}
    width = max([len(n) for n, _, _ in numbers] + [len(n) for n, _, _ in lists] + [1])
    head = ["// Parametric OpenSCAD model generated by mesh2cad.",
            "//",
            f"// Source: {source}" if source else "// Edit the parameters below and re-render.",
            "// Units are the mesh units (usually mm). Render (F6), then export STL.",
            "//",
            "// --- parameters ---"]
    lines = list(head)
    pending = list(lists)
    current = None

    def section(name: str) -> str:
        if name.startswith("cut_"):
            return "Cuts"
        if name.startswith("hole_"):
            return "Holes"
        if name.startswith("clip_"):
            return "Clips"
        return "Extrusions"

    def flush(sec: str):
        nonlocal pending
        keep = []
        for lname, pts, lnote in pending:
            if section(lname) != sec:
                keep.append((lname, pts, lnote))
                continue
            if lnote:
                lines.append(f"// {lnote}")
            lines.append(f"{lname} = {_fmt_scad(pts)};")
        pending = keep

    for name, value, note in numbers:
        sec = section(name)
        if sec != current:
            if current is not None:
                flush(current)
                lines.append("")
            lines.append(f"/* [{sec}] */")
            current = sec
        if note:
            lines.append(f"// {note}")
        line = f"{name.ljust(width)} = {value!r};"
        slider = ""
        if name.endswith("_cbore_d"):
            stem = name[: -len("_cbore_d")]
            slider = _slider(by_name.get(f"{stem}_d", 0.2),
                             max(_hole_slider_limit(positives, feats["holes"][int(stem.split("_")[1]) - 1]), value),
                             0.1, value)
        elif name.endswith("_cbore_depth"):
            stem = name[: -len("_cbore_depth")]
            hole = feats["holes"][int(stem.split("_")[1]) - 1]
            lo, hi = _numeric_span([_with_index(p) for p in positives], "XYZ".index(hole["axis"]))
            slider = _slider(0.2, max(hi - lo, value), 0.1, value)
        elif name.startswith("hole_") and name.endswith("_d"):
            hole = feats["holes"][int(name.split("_")[1]) - 1]
            slider = _slider(0.2, _hole_slider_limit(positives, hole), 0.1, value)
        if slider:
            line += f" // {slider}"
        lines.append(line)
    if current is not None:
        flush(current)
    lines += ["", "/* [Hidden] */", f"$fn = {CIRCLE_SEGMENTS}; // circle segments", ""]

    overlap = len(positives) > 1
    pieces = []
    for p, nm in zip(positives, names):
        pieces.append(_scad_placed(
            p["axis"], nm["origin"], nm["len"],
            _scad_profile(nm["outline"], p.get("cutouts") or [], nm["cutout"]),
            overlap=overlap))
    if len(pieces) == 1:
        positive = pieces[0]
    else:
        positive = "intersection() {\n" + "\n".join(_indent(p, 1) for p in pieces) + "\n}"
    sub = []
    for i, cut in enumerate(csg["cuts"], 1):
        sub.append(_scad_placed(cut["axis"], f"cut_{i}_0", f"cut_{i}_len",
                                _scad_profile(f"cut_{i}", cut.get("cutouts") or [], f"cut_{i}_cutout"),
                                overlap=True))
    for i, hole in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        ax = "XYZ".index(hole["axis"])
        lo, hi = _span_expr(positives, names, ax, "openscad")
        sub.extend(_scad_hole_pieces(n, hole, lo, hi))
    for i in range(1, len(csg["clips"]) + 1):
        sub.append(_scad_clip(i))
    model = ["// --- model ---", "module mesh2cad_part() {"]
    if sub:
        model.append("  difference() {")
        model.append(_indent(positive, 2))
        for piece in sub:
            model.append(_indent(piece, 2))
        model.append("  }")
    else:
        model.append(_indent(positive, 1))
    model += ["}", "", "mesh2cad_part();", ""]
    return "\n".join(lines + model)
