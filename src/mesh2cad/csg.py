"""Small booleans of straight extrusions, for parametric CadQuery and OpenSCAD.

A part qualifies when every region is a plane or an axis-aligned cylinder and
the solid is one of:

* the intersection of outlines extruded along X, Y, or Z;
* one of those extrusions minus straight cut extrusions along the other axes;
* either of those, plus round through-holes and counterbored through-holes;
* a planar clip where one face is tilted to all three axes (at most four).

The box vocabulary and a single straight extrusion are tried first. This
module is the fallback, and it keeps their parameter names when the part is
still one extrusion (``extrude_len``, ``outline``, ``hole_N_d``).
"""
from __future__ import annotations

import math

import numpy as np
import trimesh
from OCP.BRep import BRep_Tool
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon,
                                BRepBuilderAPI_Sewing)
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism
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
            raise Unsupported(
                "a countersink, chamfer or draft on a combined extrusion is not supported yet")
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


def _holes(m, regions, tol):
    """Through-holes and counterbored through-holes. Blind holes are refused."""
    cyls = []
    for r in regions:
        s = r.surface
        if not isinstance(s, Cylinder) or r.area < _min_area(m):
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
        za = min(s.za for s in segs)
        zb = max(s.zb for s in segs)
        lo = _snap(za, exterior, gap)
        hi = _snap(zb, exterior, gap)
        opened = abs(lo - za) <= 1e-9 and abs(hi - zb) <= 1e-9 and hi > lo + gap
        # _snap returns the original value when nothing is close, so require a real plane.
        if exterior:
            lo_hit = min(exterior, key=lambda c: abs(c - za))
            hi_hit = min(exterior, key=lambda c: abs(c - zb))
            opened = abs(lo_hit - za) <= gap and abs(hi_hit - zb) <= gap and hi_hit > lo_hit + gap
            lo, hi = (lo_hit, hi_hit) if opened else (za, zb)
        if not opened:
            raise Unsupported("hole does not pass through")
        stepped = _plain_or_counterbore(segs, lo, hi, gap)
        if stepped is None:
            if len(segs) == 1:
                raise Unsupported("hole does not pass through")
            raise Unsupported("stepped hole is not a counterbored through-hole")
        feat = {"d": stepped["d"], "axis": _AXIS_NAME[ax]}
        for i, coord in enumerate("xyz"):
            if i != ax:
                feat[coord] = float(pt[i])
        if "cbore_d" in stepped:
            feat["cbore_d"] = stepped["cbore_d"]
            feat["cbore_depth"] = stepped["cbore_depth"]
            feat["from"] = stepped["from"]
        holes.append(feat)
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
    uni = ShapeUpgrade_UnifySameDomain(sew.SewedShape(), True, True, False)
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
        shape = _add_cut(shape, _prism(spec))
    for hole in holes:
        ax = "XYZ".index(hole["axis"])
        lo, hi = _numeric_span(positives, ax)
        if hi <= lo:
            continue
        start = _hole_point(hole)
        start[ax] = lo - 1.0
        shape = _add_cut(shape, _cylinder(start, ax, float(hole["d"]) / 2.0, (hi - lo) + 2.0))
        if "cbore_d" in hole:
            depth = float(hole["cbore_depth"])
            c0 = _hole_point(hole)
            if hole.get("from", 1) < 0:
                c0[ax] = lo - 1.0
            else:
                c0[ax] = hi - depth
            shape = _add_cut(shape, _cylinder(c0, ax, float(hole["cbore_d"]) / 2.0, depth + 1.0))
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
    """A pocket that is one straight extrusion. Longest axis wins when several fit."""
    best = None
    for axis in range(3):
        if any(_attitude(r.surface.normal, axis, 4.0) == "tilted" for r in comp):
            continue
        caps = []
        for r in comp:
            if _attitude(r.surface.normal, axis, 4.0) == "cap":
                caps.append(( _plane_coord(r.surface, axis), r))
        caps.sort(key=lambda t: t[0])
        gap = _eps(tol, 0.25)
        clusters = []
        for coord, region in caps:
            if not clusters or coord - clusters[-1][-1][0] > gap:
                clusters.append([(coord, region)])
            else:
                clusters[-1].append((coord, region))
        if len(clusters) != 2:
            continue
        origin = float(np.average([c for c, _r in clusters[0]]))
        end = float(np.average([c for c, _r in clusters[-1]]))
        length = end - origin
        if length <= gap:
            continue
        verts = np.vstack([m.vertices[np.unique(m.faces[r.faces])] for r in comp])
        extent = float(verts[:, axis].max() - verts[:, axis].min())
        if best is None or extent > best[0] + 0.05:
            best = (extent, axis, origin, length, clusters)
    if best is None:
        area = max(r.area for r in comp)
        raise Unsupported(
            f"a recess is not a straight cut along X, Y, or Z "
            f"({len(comp)} faces, largest {area:.0f} mm²)")
    _extent_len, axis, origin, length, clusters = best
    U, V = _PROFILE_AXES[axis]
    walls = _walls_of(regions, axis, skip)
    # both end caps should agree; use the larger one
    loops = []
    for cluster in clusters:
        ids = []
        for _c, region in cluster:
            ids.extend(int(i) for i in region.faces)
        built = _polygons_from_cap(m, ids, walls, U, V, tol)
        if built and built[0][1] == "polygon":
            loops.append(built[0][2])
    if not loops:
        raise Unsupported(
            f"a cut along {_AXIS_NAME[axis]} has no straight end face")
    outline = _nice_loop(_drop_collinear(loops[0], tol), tol)
    if len(outline) < 3:
        raise Unsupported(f"a cut along {_AXIS_NAME[axis]} has an empty profile")
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
        shape = _solid(positives, cuts, holes, clips if use_clips else [])
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
    if not _ok(best, tol):
        # A silhouette is a superset. Extra material left inside it is a cut
        # when that material's faces form a straight extrusion.
        host = _solid([specs[a] for a in chosen], [], holes, clips if use_clips else [])
        try:
            pockets = _pocket_faces(m, regions, host, skip, tol)
            built = []
            for comp in _components(m, pockets):
                if sum(r.area for r in comp) < _min_area(m):
                    continue
                built.append(_round_spec(_cut_from_component(m, comp, regions, skip, tol), tol))
            if built:
                cuts = built
                best = score_of(chosen, cuts, use_clips)
        except Unsupported:
            best = best
    if not _ok(best, tol):
        # one outline can be slightly tight and cut off real material; drop it
        for drop in _DROP_ORDER:
            if drop not in chosen or len(chosen) == 1:
                continue
            trial = [a for a in chosen if a != drop]
            got = score_of(trial, [], use_clips)
            if _ok(got, tol):
                chosen, cuts, best = trial, [], got
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
    for i, hole in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        numbers.append((f"{n}_d", hole["d"], f"diameter (axis {hole['axis']})"))
        for coord in "xyz":
            if coord in hole:
                numbers.append((f"{n}_{coord}", hole[coord], ""))
        if "cbore_d" in hole:
            end = f"{'+' if hole.get('from', 1) > 0 else '-'}{hole['axis']}"
            numbers.append((f"{n}_cbore_d", hole["cbore_d"], "counterbore diameter"))
            numbers.append((f"{n}_cbore_depth", hole["cbore_depth"],
                            f"counterbore depth from the {end} end"))
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
                            f"cut_{i}", [], "cut")
        body.append(f"result = result.cut(cut_{i}_solid)")
    for i, hole in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        ax = "XYZ".index(hole["axis"])
        lo, hi = _span_expr(positives, names, ax, "cadquery")
        start, length = f"({lo} - 1)", f"(({hi}) - ({lo}) + 2)"
        body.append(
            f"result = result.cut(cq.Workplane().add(cq.Solid.makeCylinder("
            f"{n}_d / 2, {length}, cq.Vector({axis_point_expr(hole['axis'], n, start)}), "
            f"cq.Vector({cadquery_axis_dir(hole['axis'])}))))")
        if "cbore_d" in hole:
            if hole.get("from", 1) < 0:
                c0, clen = f"({lo} - 1)", f"{n}_cbore_depth + 1"
            else:
                c0, clen = f"(({hi}) - {n}_cbore_depth)", f"{n}_cbore_depth + 1"
            body.append(
                f"result = result.cut(cq.Workplane().add(cq.Solid.makeCylinder("
                f"{n}_cbore_d / 2, {clen}, cq.Vector({axis_point_expr(hole['axis'], n, c0)}), "
                f"cq.Vector({cadquery_axis_dir(hole['axis'])}))))")
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
                                _scad_profile(f"cut_{i}", [], "cut")))
    for i, hole in enumerate(feats.get("holes") or [], 1):
        n = f"hole_{i}"
        ax = "XYZ".index(hole["axis"])
        lo, hi = _span_expr(positives, names, ax, "openscad")
        sub.append(_scad_hole(n, hole, f"{n}_d", f"({lo} - 1)", f"(({hi}) - ({lo}) + 2)"))
        if "cbore_d" in hole:
            if hole.get("from", 1) < 0:
                c0, clen = f"({lo} - 1)", f"{n}_cbore_depth + 1"
            else:
                c0, clen = f"(({hi}) - {n}_cbore_depth)", f"{n}_cbore_depth + 1"
            sub.append(_scad_hole(n, hole, f"{n}_cbore_d", c0, clen))
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
