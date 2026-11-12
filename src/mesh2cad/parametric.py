"""Shared vocabulary for parametric CadQuery, OpenSCAD, and FreeCAD export.

Supported parts:

* a Z-up box or plate, rectangular pads on its top face, holes along X, Y, or
  Z (through, blind; countersunk when the axis is +Z), bosses along X, Y, or Z
  (straight; drafted or chamfered when the axis is +Z), one constant radius on
  all four vertical corners of the base (``fillet_r``), one constant radius on
  all twelve edges of a rectangular box (``edge_fillet_r``), a through-hole cut
  after that fillet when its wall clears the rounds, a straight boss fused
  after that fillet when its wall clears the rounds (the eight corner spheres
  stay part of ``edge_fillet_r``), a quarter-torus
  fillet on a vertical boss or hole: ``boss_N_fillet`` at the base of a straight
  boss, ``hole_N_rim_r`` / ``hole_N_bot_rim_r`` on a through-hole rim, and
  ``hole_N_floor_r`` where a blind-hole floor meets the wall, and an L-bracket
  extruded along Z whose vertical edges, including the inner corner, share one
  radius (``fillet_r``, with outer sizes ``base_x`` / ``base_y`` / ``base_t``
  and leg thicknesses ``leg_x`` / ``leg_y``), plus round through-holes in the
  foot (``base_hole_N_d``, ``base_hole_N_x``, ``base_hole_N_y``) and the
  upright (``leg_hole_N_d``, ``leg_hole_N_x``, ``leg_hole_N_y``) when the
  circle stays in that member and clears the vertical rounds. Three or more
  plain Z through-holes of one diameter, in a line at one pitch, are one row
  (``hole_row_1_n``, ``hole_row_1_pitch``, and the start and direction). An
  axis-aligned rectangle of them is one grid (``hole_grid_1_nx``,
  ``hole_grid_1_ny``, ``hole_grid_1_pitch_x``, ``hole_grid_1_pitch_y``). The
  same grouping is applied to ``base_hole_*`` and ``leg_hole_*``. Two holes,
  a bolt circle, and any hole that is blind, countersunk, counterbored, or
  not along Z stay individual parameters. A straight obround through-slot
  along Z is ``slot_N_len`` (centre to centre), ``slot_N_w`` (width, the end
  diameter), ``slot_N_x``, ``slot_N_y``, and ``slot_N_angle`` (0 along +X,
  90 along +Y);
* a straight extrusion of a 2D outline along X, Y, or Z, including a polygonal
  through-opening in that outline and round through-holes;
* a counterbored through-hole (a larger cylinder opening on one face, the
  smaller one continuing through);
* a small boolean of those extrusions: the intersection of outlines extruded
  along X, Y, or Z, or one extrusion minus cut extrusions along the other
  axes, plus a planar clip where a single face is tilted to all three axes.

:func:`extract` is the box vocabulary. :func:`recognize` uses it when the box
explains the solid, then one extrusion, then the combined extrusion.
:func:`named_parameters` is the box parameter list; extrusion parameters are
produced next to the outline. Anything outside these vocabularies raises
:class:`Unsupported`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import trimesh

from .fit import Cone, Cylinder, Plane, Sphere, Torus, fit_cylinder
from .segment import (
    _angular_gap, _outward_plane, _partial_span, _rounded_edge_radius, _rounded_primitives,
)

Z = np.array([0.0, 0.0, 1.0])


class Unsupported(ValueError):
    """The part is not in the simple parametric vocabulary."""


@dataclass
class Segment:
    kind: str  # "cylinder" | "cone"
    za: float
    zb: float
    ra: float  # radius at za
    rb: float  # radius at zb
    # end ("za" or "zb") -> (radius, kind) for a quarter-torus fillet.
    # kind is "boss" (base of an outward cylinder), "rim" (hole mouth), or
    # "floor" (blind-hole floor).
    fillets: dict = field(default_factory=dict)


@dataclass
class Stack:
    axis: int  # 0 = X, 1 = Y, 2 = Z
    point: np.ndarray  # a point on the axis
    outward: bool
    segs: list = field(default_factory=list)

    @property
    def x(self) -> float:
        return float(self.point[0])

    @property
    def y(self) -> float:
        return float(self.point[1])


def _nice(v: float, tol: float) -> float:
    for nd in (0, 1, 2):
        r = round(v, nd)
        if abs(r - v) <= max(tol, 0.005):
            return float(r) + 0.0
    return round(float(v), 3) + 0.0


def _verts(m, r):
    return m.vertices[np.unique(m.faces[r.faces])]


def _plane_pos(pl: Plane):
    """(axis index, signed outward sign, coordinate) of an axis-aligned plane."""
    i = int(np.argmax(np.abs(pl.normal)))
    return i, float(np.sign(pl.normal[i])), float(pl.d / pl.normal[i])


def _snap(v, cands, eps):
    if len(cands):
        c = min(cands, key=lambda c: abs(c - v))
        if abs(c - v) <= eps:
            return c
    return v


def _vertical_fillets(m, regions):
    """Partial vertical cylinders. Full holes and bosses cover ~360° and are not here."""
    found = []
    for r in regions:
        s = r.surface
        if not isinstance(s, Cylinder) or not len(r.faces):
            continue
        axis = np.asarray(s.axis, float)
        axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
        if abs(abs(float(axis @ Z)) - 1.0) > 0.02:
            continue
        span = _partial_span(m, s, r.faces)
        if span is None or not np.radians(25.0) <= span <= np.radians(190.0):
            continue
        found.append(r)
    return found


def _base_fillet_radius(m, fillet_regs, lo, hi, zmin, ztop, tol):
    """Radius when ``fillet_regs`` are the four vertical corners of the base box.

    The script emitters can only round that one case (CadQuery ``edges("|Z").fillet``,
    OpenSCAD ``offset``, FreeCAD ``makeFillet``). Any other fillet is left to STEP.
    """
    if len(fillet_regs) != 4:
        raise Unsupported("fillet is not the four vertical corners of a rectangular base")
    radii = [float(r.surface.radius) for r in fillet_regs]
    med = float(np.median(radii))
    if max(abs(r - med) for r in radii) > max(4.0 * tol, 0.25):
        raise Unsupported("vertical fillets do not share one radius")
    used = set()
    eps = max(4.0 * tol, 0.25)
    zeps = max(6.0 * tol, 0.5)
    for r in fillet_regs:
        rad = float(r.surface.radius)
        pt = np.asarray(r.surface.point, float)
        cands = [(lo[0] + rad, lo[1] + rad), (lo[0] + rad, hi[1] - rad),
                 (hi[0] - rad, lo[1] + rad), (hi[0] - rad, hi[1] - rad)]
        best_i, best_d = None, 1e9
        for i, (x, y) in enumerate(cands):
            d = float(np.hypot(pt[0] - x, pt[1] - y))
            if d < best_d:
                best_i, best_d = i, d
        if best_d > eps or best_i in used:
            raise Unsupported("vertical fillet is not a base corner")
        used.add(best_i)
        vz = _verts(m, r)[:, 2]
        if abs(float(vz.min()) - zmin) > zeps or abs(float(vz.max()) - ztop) > zeps:
            raise Unsupported("vertical fillet does not span the base")
    return med


def _edge_fillet_removed(x: float, y: float, z: float, radius: float) -> float:
    """Volume a sharp box loses when every edge is filleted at ``radius``.

    Each edge, away from the corners, loses a square minus a quarter circle.
    Each corner loses a cube minus a sphere octant.
    """
    r = float(radius)
    edges = 4.0 * r * r * (1.0 - math.pi / 4.0) * (x + y + z - 6.0 * r)
    corners = 8.0 * r * r * r * (1.0 - math.pi / 6.0)
    return edges + corners


def _axis_aligned_box(planes, tol: float):
    """``(lo, hi)`` of six axis-aligned planes, two opposite faces on each axis."""
    groups: dict[int, list[tuple[float, float]]] = {0: [], 1: [], 2: []}
    for plane in planes:
        n = np.abs(np.asarray(plane.normal, float))
        ln = float(np.linalg.norm(n))
        if ln < 1e-12 or float(n.max()) / ln < math.cos(math.radians(8.0)):
            return None
        i, sg, coord = _plane_pos(plane)
        groups[i].append((sg, coord))
    lo = np.zeros(3)
    hi = np.zeros(3)
    for i in range(3):
        items = groups[i]
        if len(items) != 2:
            return None
        signs = {sg for sg, _c in items}
        if signs != {-1.0, 1.0}:
            return None
        coords = [c for _sg, c in items]
        lo[i], hi[i] = min(coords), max(coords)
        if hi[i] - lo[i] <= max(8.0 * tol, 1.0):
            return None
    return lo, hi


def _inset_side(coord: float, lo: float, hi: float, radius: float, eps: float) -> int | None:
    """-1 on the low face, +1 on the high face, when ``coord`` is inset by ``radius``."""
    if abs(coord - (lo + radius)) <= eps:
        return -1
    if abs(coord - (hi - radius)) <= eps:
        return 1
    return None


def _edge_fillet_box(m, regions, tol: float) -> dict | None:
    """A rectangular box with all twelve edges filleted at one radius.

    The fitted solid is 6 planes, 12 quarter cylinders and 8 sphere octants.
    A through-hole, or a boss on the top face, is :func:`_edge_fillet_with_holes`
    when this box alone does not account for the mesh. The claim that built
    these regions can fold a small boss into the fillet, so the boss is read
    from the faces that solid does not cover.
    ``m`` is unused: the faces already carry the fitted surfaces.
    """
    del m
    planes, cyls, sphs = [], [], []
    for r in regions:
        if r.surface is None or len(r.faces) == 0:
            return None
        s = r.surface
        if isinstance(s, Plane):
            planes.append(s)
        elif isinstance(s, Cylinder):
            cyls.append(s)
        elif isinstance(s, Sphere):
            sphs.append(s)
        else:
            return None
    if len(planes) != 6 or len(cyls) != 12 or len(sphs) != 8:
        return None
    radii = [float(c.radius) for c in cyls] + [float(s.radius) for s in sphs]
    radius = float(np.median(radii))
    if radius <= 0 or max(abs(v - radius) for v in radii) > max(4.0 * tol, 0.25):
        return None
    framed = _axis_aligned_box(planes, tol)
    if framed is None:
        return None
    lo, hi = framed
    short = float(min(hi - lo))
    if radius >= 0.5 * short - 1e-6:
        return None
    eps = max(8.0 * tol, 0.5)
    corners = set()
    for sph in sphs:
        center = np.asarray(sph.center, float)
        signs = []
        for i in range(3):
            side = _inset_side(float(center[i]), float(lo[i]), float(hi[i]), radius, eps)
            if side is None:
                return None
            signs.append(side)
        key = tuple(signs)
        if key in corners:
            return None
        corners.add(key)
    if len(corners) != 8:
        return None
    edges = set()
    for cyl in cyls:
        axis = _principal_axis(cyl.axis)
        if axis is None:
            return None
        # a few degrees off the box axis is still that edge; a diagonal is not
        direction = np.asarray(cyl.axis, float)
        direction = direction / max(float(np.linalg.norm(direction)), 1e-12)
        if abs(float(direction[axis])) < math.cos(math.radians(8.0)):
            return None
        point = np.asarray(cyl.point, float)
        signs = []
        for i in range(3):
            if i == axis:
                continue
            side = _inset_side(float(point[i]), float(lo[i]), float(hi[i]), radius, eps)
            if side is None:
                return None
            signs.append((i, side))
        key = (axis, tuple(signs))
        if key in edges:
            return None
        edges.add(key)
    if len(edges) != 12:
        return None
    return {"base": dict(cx=(lo[0] + hi[0]) / 2, cy=(lo[1] + hi[1]) / 2, z0=float(lo[2]),
                         x=float(hi[0] - lo[0]), y=float(hi[1] - lo[1]), t=float(hi[2] - lo[2])),
            "pads": [], "holes": [], "bosses": [],
            "edge_fillet_r": radius}


def _end_cluster(coords, sign: int, tol: float) -> float | None:
    """Median coordinate of the vertices piled on one flat face.

    The band stays narrower than a fillet, so the round beside the face does
    not pull the measurement inward.
    """
    coords = np.asarray(coords, float)
    reach = max(0.45, 3.0 * tol)
    if sign > 0:
        tip = float(np.percentile(coords, 99.5))
        band = coords[coords >= tip - reach]
    else:
        tip = float(np.percentile(coords, 0.5))
        band = coords[coords <= tip + reach]
    if band.size < 12:
        return None
    return float(np.median(band))


def _filleted_box_frame(m, regions, tol: float):
    """``(lo, hi, flat)`` of a rectangular box, or None.

    Each side is the largest plane within 4° of that axis. A side the fitter
    folded into a fillet is the median of the vertices piled at that end.
    ``flat`` is six outward axis-aligned planes, in the pairs the fillet-radius
    fit expects.
    """
    origin = np.asarray(m.centroid, float)
    best: dict[tuple[int, int], tuple[float, float]] = {}
    min_area = max(0.02 * float(m.area), 40.0)
    align = math.cos(math.radians(4.0))
    for r in regions:
        if not isinstance(r.surface, Plane) or r.area < min_area or len(r.faces) == 0:
            continue
        got = _outward_plane(r.surface, origin)
        if got is None:
            continue
        n, d = got
        i = int(np.argmax(np.abs(n)))
        if abs(float(n[i])) < align:
            continue
        sign = 1 if float(n[i]) > 0.0 else -1
        coord = float(d / n[i])
        key = (i, sign)
        if key not in best or r.area > best[key][0]:
            best[key] = (float(r.area), coord)
    lo = np.zeros(3)
    hi = np.zeros(3)
    for i in range(3):
        pos = best.get((i, 1))
        neg = best.get((i, -1))
        high = pos[1] if pos is not None else _end_cluster(m.vertices[:, i], 1, tol)
        low = neg[1] if neg is not None else _end_cluster(m.vertices[:, i], -1, tol)
        if high is None or low is None:
            return None
        if high < low:
            low, high = high, low
        if float(high - low) <= max(8.0 * tol, 1.0):
            return None
        lo[i], hi[i] = low, high
    flat = []
    for i in range(3):
        positive = np.zeros(3)
        positive[i] = 1.0
        negative = np.zeros(3)
        negative[i] = -1.0
        flat.append((positive, float(hi[i])))
        flat.append((negative, float(-lo[i])))
    return lo, hi, flat


def _loose_axis(axis, degrees: float = 8.0) -> int | None:
    """0, 1, or 2 when ``axis`` is within ``degrees`` of X, Y, or Z."""
    a = np.asarray(axis, float)
    nrm = float(np.linalg.norm(a))
    if nrm < 1e-12:
        return None
    a = a / nrm
    i = int(np.argmax(np.abs(a)))
    if abs(float(a[i])) < math.cos(math.radians(degrees)):
        return None
    return i


def _fillet_distance(cents, flat, cyls, spheres) -> np.ndarray:
    """Distance from each face centre to the nearest filleted-box surface."""
    best = np.full(len(cents), np.inf)
    for n, d in flat:
        best = np.minimum(best, np.abs(cents @ n - d))
    for cyl in cyls:
        best = np.minimum(best, cyl.distance(cents))
    for sph in spheres:
        best = np.minimum(best, sph.distance(cents))
    return best


def _on_cylinder(pts, point, radius: float, axis: int, tol: float) -> np.ndarray:
    """Mask of points that sit on the cylinder wall, not on a neighbouring fillet."""
    rel = np.asarray(pts, float) - np.asarray(point, float)
    rel[:, axis] = 0.0
    rho = np.linalg.norm(rel, axis=1)
    return np.abs(rho - float(radius)) <= max(6.0 * tol, 0.45)


def _as_through_hole(m, group, lo, hi, fillet_r: float, tol: float) -> dict | None:
    """One coaxial group as a through-hole that clears the fillets, or None."""
    ax = int(group["axis"])
    point = np.asarray(group["point"], float)
    radius = float(group["radius"])
    if radius <= max(2.0 * tol, 0.2):
        return None
    faces = np.concatenate([r.faces for r in group["regions"]])
    pts = m.vertices[np.unique(m.faces[faces])]
    on = _on_cylinder(pts, point, radius, ax, tol)
    if int(on.sum()) < 12:
        return None
    wall = pts[on]
    along = wall[:, ax]
    eps = max(6.0 * tol, 0.75)
    if float(along.min()) > float(lo[ax]) + eps or float(along.max()) < float(hi[ax]) - eps:
        return None
    axis_vec = np.zeros(3)
    axis_vec[ax] = 1.0
    gap = _angular_gap(axis_vec, point, wall)
    if 2.0 * math.pi - gap < math.radians(270.0):
        return None
    centres = m.triangles_center[faces]
    near = _on_cylinder(centres, point, radius, ax, tol)
    if int(near.sum()) < 8:
        return None
    radial = centres[near] - point
    radial[:, ax] = 0.0
    ln = np.linalg.norm(radial, axis=1, keepdims=True)
    radial = radial / np.maximum(ln, 1e-12)
    # A hole wall points at its axis. A boss points away.
    if float(np.mean(np.sum(m.face_normals[faces][near] * radial, axis=1))) > -0.25:
        return None
    slack = max(4.0 * tol, 0.35)
    for i in range(3):
        if i == ax:
            continue
        inset = min(float(point[i] - lo[i]), float(hi[i] - point[i]))
        if inset < radius + float(fillet_r) - slack:
            return None
    feat = {"d": 2.0 * radius, "axis": "XYZ"[ax]}
    for i, name in enumerate("xyz"):
        if i != ax:
            feat[name] = float(point[i])
    return feat


def _through_holes_clear_of_fillet(m, regions, lo, hi, fillet_r: float, dist, tol: float):
    """Through-holes left over once the filleted box accounts for its faces.

    None when a leftover surface is not such a hole: a boss, a blind hole, a
    free sphere, or a round the fillet model does not cover. Regions smaller
    than half a percent of the mesh are segmentation crumbs and are ignored.
    """
    area = m.area_faces
    total = float(area.sum()) or 1.0
    crumb = 0.005 * total
    thresh = max(8.0 * tol, 0.15 * float(fillet_r))
    pieces = []
    for r in regions:
        if r.area <= 0 or len(r.faces) == 0:
            continue
        outside = dist[r.faces] > thresh
        frac = float(area[r.faces][outside].sum()) / float(r.area)
        if frac <= 0.5 or r.area < crumb:
            continue
        if not isinstance(r.surface, Cylinder):
            return None
        ax = _loose_axis(r.surface.axis)
        if ax is None:
            return None
        pieces.append((r, ax, np.asarray(r.surface.point, float), float(r.surface.radius)))
    if not pieces:
        return None
    groups = []
    gtol = max(5.0 * tol, 0.2)
    rtol = max(4.0 * tol, 0.25)
    for region, ax, point, radius in pieces:
        for g in groups:
            if g["axis"] != ax or abs(radius - g["radius"]) > rtol:
                continue
            if _perp_dist(ax, g["point"], point) > gtol:
                continue
            w = g["area"] + float(region.area)
            g["point"] = (g["point"] * g["area"] + point * float(region.area)) / w
            g["radius"] = (g["radius"] * g["area"] + radius * float(region.area)) / w
            g["area"] = w
            g["regions"].append(region)
            break
        else:
            groups.append({"axis": ax, "point": point.copy(), "radius": radius,
                           "area": float(region.area), "regions": [region]})
    holes = []
    for g in groups:
        hole = _as_through_hole(m, g, lo, hi, fillet_r, tol)
        if hole is None:
            return None
        holes.append(hole)
    holes.sort(key=lambda h: ("XYZ".index(h["axis"]), tuple(h.get(c, 0.0) for c in "xyz")))
    return holes


def _refit_boss_circle(m, point, radius, z0, z1, tol):
    """Center and radius from vertices on the wall, not from inset face centres.

    The leftover faces start above the top plane, so a fit on those centres
    alone drifts. Vertices on the whole wall, including the part the fillet
    model still claims, pull the axis back.
    """
    band = max(6.0 * tol, 0.45)
    verts = m.vertices
    rel = verts - np.asarray(point, float)
    rel[:, 2] = 0.0
    rho = np.linalg.norm(rel, axis=1)
    z = verts[:, 2]
    sel = (np.abs(rho - float(radius)) <= band) & (z >= float(z0) - band) & (z <= float(z1) + band)
    pts = verts[sel]
    if len(pts) < 40:
        return np.asarray(point, float), float(radius)
    x, y = pts[:, 0], pts[:, 1]
    fit = np.c_[2.0 * x, 2.0 * y, np.ones(len(x))]
    sol, *_ = np.linalg.lstsq(fit, x * x + y * y, rcond=None)
    cx, cy, k = float(sol[0]), float(sol[1]), float(sol[2])
    r2 = k + cx * cx + cy * cy
    if r2 <= 0:
        return np.asarray(point, float), float(radius)
    rad = math.sqrt(r2)
    if abs(rad - float(radius)) > max(4.0 * tol, 0.25):
        return np.asarray(point, float), float(radius)
    if math.hypot(cx - float(point[0]), cy - float(point[1])) > max(4.0 * tol, 0.4):
        return np.asarray(point, float), float(radius)
    return np.array([cx, cy, 0.0]), rad


def _top_bosses_clear_of_fillet(m, lo, hi, fillet_r: float, dist, tol: float):
    """Straight bosses on the top face left over once the filleted box is accounted for.

    None when nothing is left over, or what is left over is not such a boss.
    The rounded-box claim folds a small boss into the fillet regions, so the
    wall is the faces the fillet model does not cover, not a cylinder region.
    The script fillets the box first and fuses the boss after.
    """
    area = m.area_faces
    total = float(area.sum()) or 1.0
    crumb = 0.005 * total
    thresh = max(8.0 * tol, 0.15 * float(fillet_r))
    outside = np.flatnonzero(dist > thresh)
    if outside.size == 0 or float(area[outside].sum()) < crumb:
        return None
    norms = m.face_normals
    cents = m.triangles_center
    above = float(hi[2]) + max(4.0 * tol, 0.4)
    cap_sel = (norms[outside, 2] > 0.70) & (cents[outside, 2] > above)
    wall_sel = np.abs(norms[outside, 2]) < 0.50
    cap_ids = outside[cap_sel]
    wall_ids = outside[wall_sel]
    other = outside[~(cap_sel | wall_sel)]
    outside_area = float(area[outside].sum())
    if float(area[other].sum()) > max(2.0 * crumb, 0.20 * outside_area):
        return None
    if cap_ids.size < 8 or wall_ids.size < 24:
        return None
    cyl = fit_cylinder(cents[wall_ids], norms[wall_ids])
    if cyl is None or _loose_axis(cyl.axis) != 2:
        return None
    point = np.asarray(cyl.point, float)
    radius = float(cyl.radius)
    if radius <= max(2.0 * tol, 0.2):
        return None
    on_wall = _on_cylinder(cents[wall_ids], point, radius, 2, tol)
    if float(area[wall_ids][on_wall].sum()) < 0.70 * float(area[wall_ids].sum()):
        return None
    wall_ids = wall_ids[on_wall]
    if wall_ids.size < 24:
        return None
    z_top = float(np.percentile(cents[wall_ids, 2], 98))
    point, radius = _refit_boss_circle(m, point, radius, float(hi[2]), z_top, tol)
    radial = cents[wall_ids] - point
    radial[:, 2] = 0.0
    ln = np.linalg.norm(radial, axis=1, keepdims=True)
    radial = radial / np.maximum(ln, 1e-12)
    if float(np.mean(np.sum(norms[wall_ids] * radial, axis=1))) < 0.25:
        return None
    gap = _angular_gap(np.array([0.0, 0.0, 1.0]), point, cents[wall_ids])
    if 2.0 * math.pi - gap < math.radians(270.0):
        return None
    cap_rho = np.linalg.norm(cents[cap_ids, :2] - point[:2], axis=1)
    cap_near = cap_rho <= radius + max(4.0 * tol, 0.5)
    cap_ids = cap_ids[cap_near]
    if cap_ids.size < 8:
        return None
    cap_z = float(np.median(cents[cap_ids, 2]))
    face_z = float(hi[2])
    height = cap_z - face_z
    if height <= max(3.0 * tol, 0.5):
        return None
    if float(np.percentile(np.abs(cents[cap_ids, 2] - cap_z), 90)) > max(6.0 * tol, 0.45):
        return None
    disk = math.pi * radius * radius
    cap_area = float(area[cap_ids].sum())
    if cap_area < 0.45 * disk or cap_area > 1.45 * disk:
        return None
    wall_z = cents[wall_ids, 2]
    # Faces within ``thresh`` of the top plane stay with the fillet model, so
    # the leftover wall begins a little above the face it actually meets.
    if float(wall_z.min()) > face_z + thresh + max(6.0 * tol, 0.8):
        return None
    if float(wall_z.max()) < cap_z - max(6.0 * tol, 1.0):
        return None
    expect = 2.0 * math.pi * radius * height
    if float(area[wall_ids].sum()) < 0.35 * expect:
        return None
    slack = max(4.0 * tol, 0.35)
    for i in (0, 1):
        inset = min(float(point[i] - lo[i]), float(hi[i] - point[i]))
        if inset < radius + float(fillet_r) - slack:
            return None
    boss = {"d": 2.0 * radius, "h": float(height), "axis": "Z", "from": 1,
            "x": float(point[0]), "y": float(point[1])}
    return [boss]


def _edge_fillet_with_holes(m, regions, tol: float) -> dict | None:
    """A fully filleted box with a through-hole or a top boss that clears the rounds.

    The twelve edges and eight corner spheres are ``edge_fillet_r``. A fitted
    sphere or cone that lies on those surfaces is that fillet, not a free-form
    sphere. A leftover cylinder that runs from one face to the opposite face
    and points inward is a through-hole, cut after the fillet. Faces the fillet
    model does not cover, forming an outward cylinder on the top face, are a
    boss fused after the fillet. Either way the corner spheres stay part of
    the same radius.
    """
    if not any(isinstance(r.surface, Cylinder) and r.area > 0 for r in regions):
        return None
    framed = _filleted_box_frame(m, regions, tol)
    if framed is None:
        return None
    lo, hi, flat = framed
    separations = [float(hi[i] - lo[i]) for i in range(3)]
    radius = _rounded_edge_radius(m, flat, separations, tol)
    if radius is None or radius >= 0.5 * float(min(separations)) - 1e-6:
        return None
    built = _rounded_primitives(flat, radius)
    if built is None:
        return None
    spheres, cyls = built
    dist = _fillet_distance(m.triangles_center, flat, cyls, spheres)
    holes = _through_holes_clear_of_fillet(m, regions, lo, hi, radius, dist, tol)
    bosses = []
    if not holes:
        bosses = _top_bosses_clear_of_fillet(m, lo, hi, radius, dist, tol) or []
        if not bosses:
            return None
        holes = []
    feats = {
        "base": dict(cx=float(lo[0] + hi[0]) / 2, cy=float(lo[1] + hi[1]) / 2, z0=float(lo[2]),
                     x=float(hi[0] - lo[0]), y=float(hi[1] - lo[1]), t=float(hi[2] - lo[2])),
        "pads": [], "holes": holes, "bosses": bosses,
        "edge_fillet_r": float(radius),
    }
    if not _box_explains_mesh(m, feats):
        return None
    return feats


def _principal_axis(axis) -> int | None:
    """0, 1, or 2 when ``axis`` lies along X, Y, or Z. None otherwise."""
    a = np.asarray(axis, float)
    nrm = float(np.linalg.norm(a))
    if nrm < 1e-12:
        return None
    a = a / nrm
    i = int(np.argmax(np.abs(a)))
    if abs(abs(float(a[i])) - 1.0) > 1e-3:
        return None
    return i


def _cluster_end(coords, extreme: float, sign: int, tol: float) -> float:
    """Median axial coordinate of the vertices piled up at a free end.

    Used when that end did not snap to a plane. A blind cap absorbed into the
    cylinder region has no plane of its own, and under noise the single
    furthest vertex sits past the cap.
    """
    band = max(6.0 * tol, 0.3)
    c = np.asarray(coords, float)
    near = c[c <= extreme + band] if sign < 0 else c[c >= extreme - band]
    if near.size < 8:
        return float(extreme)
    return float(np.median(near))


def _perp_dist(axis: int, p, q) -> float:
    d = np.asarray(p, float) - np.asarray(q, float)
    d[axis] = 0.0
    return float(np.linalg.norm(d))


def _stacks(m, regions, tol, planes):
    """Group coaxial cylinder/cone pieces.

    ``planes`` maps an axis index to the plane coordinates along that axis, used
    to snap each piece's ends onto the faces it actually meets. A cone is kept
    only on a vertical axis (countersink, draft, chamfer); a cone on X or Y is
    outside the vocabulary.
    """
    stacks: list[Stack] = []
    gtol = max(5 * tol, 0.05)
    for r in regions:
        s = r.surface
        if not isinstance(s, (Cylinder, Cone)):
            continue
        ax = _principal_axis(s.axis)
        if ax is None:
            raise Unsupported("curved surface whose axis is not aligned with X, Y or Z")
        if isinstance(s, Cone) and ax != 2:
            raise Unsupported("countersink, chamfer or draft whose axis is not vertical")
        c = np.asarray(s.point if isinstance(s, Cylinder) else s.apex, float)
        v = _verts(m, r)
        rel = v - c
        rel[:, ax] = 0.0
        rho = np.linalg.norm(rel, axis=1)
        eps = max(4 * tol, 0.05)
        coords = planes.get(ax, ())
        along = v[:, ax]
        za = _snap(float(along.min()), coords, eps)
        zb = _snap(float(along.max()), coords, eps)
        # A blind end whose disk was absorbed into the cylinder has no plane to
        # snap to, so the furthest noisy vertex would set the depth. Use the
        # median of the vertices clustered at that free end instead. An end
        # that already snapped to a face is left alone.
        if abs(za - float(along.min())) <= 1e-9:
            za = _snap(_cluster_end(along, za, -1, tol), coords, eps)
        if abs(zb - float(along.max())) <= 1e-9:
            zb = _snap(_cluster_end(along, zb, 1, tol), coords, eps)
        if isinstance(s, Cylinder):
            ra = rb = float(s.radius)
            kind = "cylinder"
        else:
            b, a = np.polyfit(v[:, ax], rho, 1)
            ra, rb, kind = float(a + b * za), float(a + b * zb), "cone"
        fc = m.triangles_center[r.faces]
        radial = fc - c
        radial[:, ax] = 0.0
        radial /= np.maximum(np.linalg.norm(radial, axis=1, keepdims=True), 1e-12)
        out = float(np.mean(np.sum(m.face_normals[r.faces] * radial, axis=1))) > 0
        for st in stacks:
            if st.axis == ax and st.outward == out and _perp_dist(ax, st.point, c) < gtol:
                st.segs.append(Segment(kind, za, zb, ra, rb))
                break
        else:
            stacks.append(Stack(ax, c.copy(), out, [Segment(kind, za, zb, ra, rb)]))
    for st in stacks:
        st.segs.sort(key=lambda g: g.za)
        for g, h in zip(st.segs, st.segs[1:]):  # make junctions meet
            if abs(g.zb - h.za) <= max(4 * tol, 0.05):
                zj = (g.zb + h.za) / 2
                g.zb = h.za = zj
    for st in stacks:
        if st.axis != 2:
            continue
        for g in st.segs:
            if g.kind == "cone":
                k = next(r.surface for r in regions if isinstance(r.surface, Cone)
                         and np.hypot(*(r.surface.apex[:2] - st.point[:2])) < gtol)
                g.ra = abs(k.apex[2] - g.za) * math.tan(k.half_angle)
                g.rb = abs(k.apex[2] - g.zb) * math.tan(k.half_angle)
    return stacks


def _axis_name(axis: int) -> str:
    return "XYZ"[axis]


def _plain_or_counterbore(segs, lo: float, hi: float, eps: float) -> dict | None:
    """A through cylinder, or a counterbored through-hole, spanning ``lo``..``hi``.

    ``segs`` are sorted along the axis. A counterbore is the larger cylinder
    opening on one end and the smaller one continuing to the other end.
    Returns None when the stack is a blind hole or not this pattern; the
    caller still handles blind holes and countersinks.
    """
    if not segs or any(getattr(s, "kind", "cylinder") != "cylinder" for s in segs):
        return None
    if len(segs) == 1:
        c = segs[0]
        if abs(c.za - lo) <= eps and abs(c.zb - hi) <= eps:
            return {"d": 2.0 * float(c.ra)}
        return None
    if len(segs) != 2:
        return None
    a, b = segs
    if abs(a.zb - b.za) > eps or abs(a.za - lo) > eps or abs(b.zb - hi) > eps:
        return None
    if float(b.ra) > float(a.ra) + eps:
        return {"d": 2.0 * float(a.ra), "cbore_d": 2.0 * float(b.ra),
                "cbore_depth": float(b.zb - b.za), "from": 1}
    if float(a.ra) > float(b.ra) + eps:
        return {"d": 2.0 * float(b.ra), "cbore_d": 2.0 * float(a.ra),
                "cbore_depth": float(a.zb - a.za), "from": -1}
    return None


def _face_name(axis: str, sign: int) -> str:
    return f"{'+' if sign > 0 else '-'}{axis}"


def base_low_expr(axis: str) -> str:
    return {"X": "(base_cx - base_x / 2)", "Y": "(base_cy - base_y / 2)", "Z": "base_z0"}[axis]


def base_high_expr(axis: str) -> str:
    return {"X": "(base_cx + base_x / 2)", "Y": "(base_cy + base_y / 2)",
            "Z": "(base_z0 + base_t)"}[axis]


def base_size_expr(axis: str) -> str:
    return {"X": "base_x", "Y": "base_y", "Z": "base_t"}[axis]


def axis_point_expr(axis: str, prefix: str, along: str) -> str:
    """``x, y, z`` with the coordinate along ``axis`` replaced by ``along``."""
    coords = []
    for name in "xyz":
        coords.append(along if name == axis.lower() else f"{prefix}_{name}")
    return ", ".join(coords)


def hole_cut_span(prefix: str, hole: dict) -> tuple[str, str]:
    """Start and length of a side-hole cutter, along +axis, overlapping the mouth by 1."""
    ax = hole["axis"]
    lo, hi = base_low_expr(ax), base_high_expr(ax)
    if "depth" not in hole:
        return f"{lo} - 1", f"{base_size_expr(ax)} + 2"
    if hole["from"] < 0:
        return f"{lo} - 1", f"{prefix}_depth + 1"
    return f"{hi} - {prefix}_depth", f"{prefix}_depth + 1"


def cbore_cut_span(prefix: str, hole: dict) -> tuple[str, str]:
    """Start and length of a counterbore cutter, along +axis, overlapping the mouth by 1."""
    ax = hole.get("axis", "Z")
    lo, hi = base_low_expr(ax), base_high_expr(ax)
    if hole.get("from", 1) < 0:
        return f"{lo} - 1", f"{prefix}_cbore_depth + 1"
    return f"{hi} - {prefix}_cbore_depth", f"{prefix}_cbore_depth + 1"


def boss_span(prefix: str, boss: dict) -> tuple[str, str]:
    """Start and length of a side boss, growing outward from the face it sits on."""
    ax = boss["axis"]
    if boss["from"] < 0:
        return f"{base_low_expr(ax)} - {prefix}_h", f"{prefix}_h"
    return base_high_expr(ax), f"{prefix}_h"


def quarter_ring_volume(wall_r: float, minor: float, outside: bool) -> float:
    """Volume of one revolved square-minus-quarter-disk.

    ``outside`` is true when the ring sits outside the cylinder wall (a boss
    base fillet, or the round removed from a hole rim) and false when it sits
    inside the wall (the material a blind-hole floor fillet puts back).
    """
    if minor <= 0.0 or wall_r <= 0.0:
        return 0.0
    area = minor * minor * (1.0 - math.pi / 4.0)
    # distance from the sharp corner to the centroid, toward the square interior
    frac = (0.5 - (math.pi / 4.0) * (1.0 - 4.0 / (3.0 * math.pi))) / (1.0 - math.pi / 4.0)
    centre = wall_r + (minor * frac if outside else -minor * frac)
    if centre <= 0.0:
        return 0.0
    return 2.0 * math.pi * centre * area


def _apply_torus_fillets(stacks, tori, plane_coords, tol: float) -> None:
    """Extend each straight cylinder through a tangent quarter-torus fillet.

    The cylinder region stops where the fillet starts, so the nominal boss or
    hole would come out short. The end is moved out to the tangent plane and
    the minor radius is stored on that end. A fillet whose axis is not Z is
    refused: the emitters only write the vertical case.
    """
    if not tori:
        return
    slack = max(8.0 * tol, 0.4)
    snap = max(10.0 * tol, 0.55)
    for reg in tori:
        torus = reg.surface
        if not isinstance(torus, Torus) or torus.minor <= 0.0 or torus.major <= 0.0:
            raise Unsupported("torus fillet has no usable radius")
        ax = _principal_axis(torus.axis)
        if ax is None:
            raise Unsupported("torus fillet whose axis is not aligned with X, Y or Z")
        if ax != 2:
            raise Unsupported(
                "torus fillet on a hole or boss whose axis is not Z is not supported yet"
            )
        direction = np.asarray(torus.axis, float)
        direction = direction / max(float(np.linalg.norm(direction)), 1e-12)
        centre = np.asarray(torus.center, float)
        plane_c = float(centre[ax] + direction[ax] * float(torus.minor))
        plane_c = _snap(plane_c, plane_coords.get(ax, ()), snap)
        centre_c = float(centre[ax])
        minor = float(torus.minor)
        wall = (float(torus.major) - minor) if torus.inner else (float(torus.major) + minor)
        coaxial = [st for st in stacks
                   if st.axis == ax and _perp_dist(ax, st.point, centre) <= max(5.0 * tol, 0.25)]
        if not coaxial:
            raise Unsupported(
                "torus fillet is not the base of a straight boss, the rim of a hole, "
                "or the floor of a blind hole"
            )
        match = None
        radius_hit = False
        for st in coaxial:
            close = [s for s in st.segs if abs(float(s.ra) - wall) <= max(slack, 0.5)
                     or abs(float(s.rb) - wall) <= max(slack, 0.5)]
            if not close:
                continue
            radius_hit = True
            if len(st.segs) != 1 or st.segs[0].kind != "cylinder":
                raise Unsupported(
                    "torus fillet on a drafted, chamfered or stepped feature is not supported yet"
                )
            seg = st.segs[0]
            for end, end_c in (("za", float(seg.za)), ("zb", float(seg.zb))):
                if end in seg.fillets:
                    continue
                if abs(end_c - centre_c) > max(slack, 0.55 * minor):
                    continue
                other = float(seg.zb if end == "za" else seg.za)
                # the plane is beyond this end, away from the rest of the cylinder
                if (plane_c - end_c) * (other - end_c) > 0.0:
                    continue
                if abs(abs(plane_c - end_c) - minor) > max(slack, 0.45 * minor):
                    continue
                match = (seg, end, plane_c, st)
                break
            if match:
                break
        if match is None:
            if not radius_hit:
                raise Unsupported(
                    "torus fillet is not the base of a straight boss, the rim of a hole, "
                    "or the floor of a blind hole"
                )
            raise Unsupported(
                "torus fillet does not meet the end of its cylinder"
            )
        seg, end, plane_c, st = match
        if st.outward:
            if not torus.inner:
                raise Unsupported("torus on a boss is not a base fillet")
            kind = "boss"
        else:
            kind = "floor" if not torus.inner else "rim"
        if end == "za":
            seg.za = plane_c
        else:
            seg.zb = plane_c
        if seg.za > seg.zb:
            raise Unsupported("torus fillet reverses its cylinder")
        seg.fillets[end] = (minor, kind)


def _hole_fillets(seg, hole: dict) -> None:
    """Copy rim and floor fillets from ``seg`` onto a hole feature."""
    for end, rim_key in (("zb", "rim_r"), ("za", "rim_lo")):
        info = seg.fillets.get(end)
        if not info:
            continue
        radius, kind = info
        if kind == "rim":
            hole[rim_key] = radius
        elif kind == "floor":
            if "floor_r" in hole:
                raise Unsupported("a hole has more than one floor fillet")
            hole["floor_r"] = radius
            hole["floor_end"] = "high" if end == "zb" else "low"
        else:
            raise Unsupported("torus on a hole is not a rim or a floor fillet")


def _boss_fillet(seg, from_sign: int, boss: dict) -> None:
    """Copy a base fillet onto a boss. A round on the free end is refused."""
    base = "za" if from_sign > 0 else "zb"
    free = "zb" if base == "za" else "za"
    if seg.fillets.get(free):
        raise Unsupported("round on the free end of a boss is not supported yet")
    info = seg.fillets.get(base)
    if not info:
        return
    radius, kind = info
    if kind != "boss":
        raise Unsupported("torus on a boss is not a base fillet")
    boss["fillet"] = radius


def cadquery_axis_dir(axis: str) -> str:
    return {"X": "1, 0, 0", "Y": "0, 1, 0", "Z": "0, 0, 1"}[axis]


def openscad_axis_rotate(axis: str) -> str:
    """Rotation taking OpenSCAD's +Z cylinder onto +X or +Y."""
    return {"X": "[0, 90, 0]", "Y": "[-90, 0, 0]", "Z": "[0, 0, 0]"}[axis]


def extract(m: trimesh.Trimesh, regions, tol: float) -> dict:
    """Recognise base, pads, holes and bosses. Returns a feature dict."""
    if any(r.surface is None for r in regions):
        raise Unsupported("part has unfitted (faceted) regions")
    planes = [r for r in regions if isinstance(r.surface, Plane)]
    tori = [r for r in regions if isinstance(r.surface, Torus)]
    for r in regions:
        if not isinstance(r.surface, (Plane, Cylinder, Cone, Torus)):
            raise Unsupported(f"{r.surface.kind} faces are not supported yet")
    # Round ends of a slot are partial cylinders, or a ring of flat facets on
    # a clean tessellation. Claiming them first keeps those facets from being
    # read as a tilted plane, and keeps the ends from being read as fillets.
    from .slots import find_obround_slots
    slots, slot_planes, slot_cyls = find_obround_slots(m, regions, tol)
    pos = []
    for r in planes:
        if id(r) in slot_planes:
            continue
        n = np.abs(r.surface.normal)
        if n.max() < math.cos(math.radians(2)):
            raise Unsupported("plane that is not axis-aligned")
        pos.append((r, *_plane_pos(r.surface)))
    down = [(r, c) for r, i, sg, c in pos if i == 2 and sg < 0]
    if not down:
        raise Unsupported("no bottom face")
    zmin = min(down, key=lambda t: -t[0].area)[1]
    plane_coords = {0: [], 1: [], 2: []}
    for _, i, _, c in pos:
        plane_coords[i].append(c)
    plane_coords[2].append(zmin)
    # corner fillets are partial cylinders; leaving them in the stacks makes
    # them look like bosses that do not start on the top face
    fillet_regs = [r for r in _vertical_fillets(m, regions) if id(r) not in slot_cyls]
    fillet_ids = {id(r) for r in fillet_regs}
    stacks = _stacks(m, [r for r in regions
                         if id(r) not in fillet_ids and id(r) not in slot_cyls
                         and not isinstance(r.surface, Torus)],
                     tol, plane_coords)
    _apply_torus_fillets(stacks, tori, plane_coords, tol)

    def on_stack(r):
        v = _verts(m, r)
        for st in stacks:
            rmax = max(max(g.ra, g.rb) for g in st.segs)
            d = v - st.point
            d[:, st.axis] = 0.0
            if np.percentile(np.linalg.norm(d, axis=1), 95) <= rmax + max(3 * tol, 0.05):
                return True
        return False

    up = []
    for r in planes:
        nz = r.surface.normal @ Z
        z = _plane_pos(r.surface)[2] if abs(nz) > 0.99 else 0.0
        # a counterbore floor faces into the hole; it is on a stack, not an overhang
        if nz < -0.999 and z > zmin + max(3 * tol, 0.05) and not on_stack(r):
            raise Unsupported("downward-facing face above the bottom (overhang)")
        if nz > 0.999:
            up.append((r, z))
    if not up:
        raise Unsupported("no top face")
    base_r, ztop = max(up, key=lambda t: t[0].area)
    # Pad edges snap to every side plane. The base box does not: a boss end
    # cap sits outside the box and is on a stack, so it must not move the wall.
    xs = [c for _, i, _, c in pos if i == 0]
    ys = [c for _, i, _, c in pos if i == 1]
    xs_base = [c for r, i, _, c in pos if i == 0 and not on_stack(r)]
    ys_base = [c for r, i, _, c in pos if i == 1 and not on_stack(r)]
    if not xs_base or not ys_base:
        raise Unsupported("no side faces")
    lo, hi = np.array([min(xs_base), min(ys_base)]), np.array([max(xs_base), max(ys_base)])
    box_lo = (lo[0], lo[1], zmin)
    box_hi = (hi[0], hi[1], ztop)
    feats = {"base": dict(cx=(lo[0] + hi[0]) / 2, cy=(lo[1] + hi[1]) / 2, z0=zmin,
                          x=hi[0] - lo[0], y=hi[1] - lo[1], t=ztop - zmin),
             "pads": [], "holes": [], "bosses": []}
    eps = max(3 * tol, 0.05)
    for r, z in up:
        if r is base_r or on_stack(r):
            continue
        # a pad that crosses the top splits that face into coplanar pieces.
        # They share z with the base top; a pocket sits below it.
        if abs(z - ztop) <= eps:
            continue
        if z < ztop:
            raise Unsupported("pocket or step below the top face")
        pv = _verts(m, r)
        e2 = max(4 * tol, 0.1)
        plo = np.array([_snap(pv[:, 0].min(), xs, e2), _snap(pv[:, 1].min(), ys, e2)])
        phi = np.array([_snap(pv[:, 0].max(), xs, e2), _snap(pv[:, 1].max(), ys, e2)])
        feats["pads"].append(dict(cx=(plo[0] + phi[0]) / 2, cy=(plo[1] + phi[1]) / 2,
                                  w=phi[0] - plo[0], l=phi[1] - plo[1], h=z - ztop))
    def _cross_section_ok(st):
        for i in range(3):
            if i == st.axis:
                continue
            if not (box_lo[i] - eps <= float(st.point[i]) <= box_hi[i] + eps):
                return False
        return True

    for st in stacks:
        g = st.segs
        if st.axis != 2:
            if not _cross_section_ok(st):
                raise Unsupported("side feature axis is outside the base")
            if st.outward and (len(g) != 1 or g[0].kind != "cylinder"):
                raise Unsupported("side boss is not a plain cylinder")
            if not st.outward and any(s.kind != "cylinder" for s in g):
                raise Unsupported("side hole is not a plain cylinder")
            lo_a, hi_a = box_lo[st.axis], box_hi[st.axis]
            stepped = None if st.outward else _plain_or_counterbore(g, lo_a, hi_a, eps)
            if stepped and "cbore_d" in stepped:
                feat = {name: float(st.point[i]) for i, name in enumerate("xyz") if i != st.axis}
                feat["d"] = stepped["d"]
                feat["axis"] = _axis_name(st.axis)
                feat["cbore_d"] = stepped["cbore_d"]
                feat["cbore_depth"] = stepped["cbore_depth"]
                feat["from"] = stepped["from"]
                feats["holes"].append(feat)
                continue
            if len(g) != 1 or g[0].kind != "cylinder":
                raise Unsupported("side hole is not a plain or counterbored cylinder")
            c = g[0]
            feat = {name: float(st.point[i]) for i, name in enumerate("xyz") if i != st.axis}
            feat["d"] = 2 * c.ra
            feat["axis"] = _axis_name(st.axis)
            if not st.outward:
                on_lo = abs(c.za - lo_a) <= eps
                on_hi = abs(c.zb - hi_a) <= eps
                if on_lo and on_hi:
                    pass
                elif on_lo and c.zb < hi_a - eps:
                    feat["depth"] = c.zb - lo_a
                    feat["from"] = -1
                elif on_hi and c.za > lo_a + eps:
                    feat["depth"] = hi_a - c.za
                    feat["from"] = 1
                else:
                    raise Unsupported("side hole does not open on a face of the base")
                _hole_fillets(c, feat)
                feats["holes"].append(feat)
            else:
                # the wall starts on a side face and runs outward, past that face
                if abs(c.za - hi_a) <= eps and c.zb > hi_a + eps:
                    feat["h"] = c.zb - hi_a
                    feat["from"] = 1
                elif abs(c.zb - lo_a) <= eps and c.za < lo_a - eps:
                    feat["h"] = lo_a - c.za
                    feat["from"] = -1
                else:
                    raise Unsupported("side boss does not start on a side face of the base")
                _boss_fillet(c, feat["from"], feat)
                feats["bosses"].append(feat)
            continue
        if not st.outward:
            h = dict(x=st.x, y=st.y, axis="Z")
            stepped = _plain_or_counterbore(g, zmin, ztop, eps)
            if stepped and "cbore_d" in stepped:
                h["d"] = stepped["d"]
                h["cbore_d"] = stepped["cbore_d"]
                h["cbore_depth"] = stepped["cbore_depth"]
                h["from"] = stepped["from"]
                feats["holes"].append(h)
                continue
            cyl = [s for s in g if s.kind == "cylinder"]
            if len(cyl) != 1 or cyl[0] is not g[0]:
                raise Unsupported("hole profile is not a plain, counterbored or countersunk hole")
            c = cyl[0]
            h["d"] = 2 * c.ra
            top_of_hole = c.zb
            if len(g) == 2 and g[1].kind == "cone" and g[1].rb > g[1].ra:
                k = g[1]
                h["csk_d"] = 2 * k.rb
                h["csk_angle"] = 2 * math.degrees(math.atan2(k.rb - k.ra, k.zb - k.za))
                top_of_hole = k.zb
            elif len(g) != 1:
                raise Unsupported("hole profile is not a plain or countersunk hole")
            if abs(top_of_hole - ztop) > eps:
                raise Unsupported("hole does not start at the top face")
            if c.za > zmin + eps:
                h["depth"] = ztop - c.za
                h["from"] = 1
            _hole_fillets(c, h)
            feats["holes"].append(h)
        else:
            b = dict(x=st.x, y=st.y, axis="Z")
            if abs(g[0].za - ztop) > eps:
                raise Unsupported("boss does not start on the top face")
            b["d"] = 2 * g[0].ra
            b["h"] = g[-1].zb - ztop
            b["from"] = 1
            kinds = [s.kind for s in g]
            if kinds == ["cylinder"]:
                pass
            elif kinds == ["cone"]:
                b["draft_deg"] = math.degrees(math.atan2(g[0].ra - g[0].rb, g[0].zb - g[0].za))
            elif kinds == ["cylinder", "cone"] and g[1].rb < g[1].ra:
                b["chamfer"] = g[1].zb - g[1].za
                b["chamfer_w"] = g[1].ra - g[1].rb
            else:
                raise Unsupported("boss profile is not straight, drafted or chamfered")
            _boss_fillet(g[0], b["from"], b)
            feats["bosses"].append(b)
    if fillet_regs:
        feats["fillet_r"] = _base_fillet_radius(m, fillet_regs, lo, hi, zmin, ztop, tol)
    if slots:
        feats["slots"] = slots
    return feats


def _feature_volume(feats: dict) -> float:
    """Analytic volume of a box-vocabulary feature dict."""
    b = feats["base"]
    v = float(b["x"]) * float(b["y"]) * float(b["t"])
    if "edge_fillet_r" in feats:
        v -= _edge_fillet_removed(float(b["x"]), float(b["y"]), float(b["t"]), float(feats["edge_fillet_r"]))
    elif "fillet_r" in feats:
        r = float(feats["fillet_r"])
        v -= 4.0 * float(b["t"]) * r * r * (1.0 - math.pi / 4.0)

    def cyl(radius, height):
        return math.pi * radius * radius * height

    def frust(r1, r2, height):
        return math.pi * height * (r1 * r1 + r1 * r2 + r2 * r2) / 3.0

    for pad in feats["pads"]:
        v += float(pad["w"]) * float(pad["l"]) * float(pad["h"])
    for boss in feats["bosses"]:
        radius, height = float(boss["d"]) / 2.0, float(boss["h"])
        if "draft_deg" in boss:
            r2 = radius - height * math.tan(math.radians(float(boss["draft_deg"])))
            v += frust(radius, r2, height)
        elif "chamfer" in boss:
            ch = float(boss["chamfer"])
            cw = float(boss.get("chamfer_w", ch))
            v += cyl(radius, height - ch) + frust(radius, radius - cw, ch)
        else:
            v += cyl(radius, height)
        if "fillet" in boss:
            v += quarter_ring_volume(radius, float(boss["fillet"]), True)
    for hole in feats["holes"]:
        radius = float(hole["d"]) / 2.0
        ax = hole.get("axis", "Z")
        depth = float(hole.get("depth", {"X": b["x"], "Y": b["y"], "Z": b["t"]}[ax]))
        v -= cyl(radius, depth)
        if "csk_d" in hole:
            big = float(hole["csk_d"]) / 2.0
            csk_h = (big - radius) / math.tan(math.radians(float(hole["csk_angle"]) / 2.0))
            v -= frust(radius, big, csk_h) - cyl(radius, csk_h)
        if "cbore_d" in hole:
            big = float(hole["cbore_d"]) / 2.0
            v -= math.pi * (big * big - radius * radius) * float(hole["cbore_depth"])
        if "rim_r" in hole:
            v -= quarter_ring_volume(radius, float(hole["rim_r"]), True)
        if "rim_lo" in hole:
            v -= quarter_ring_volume(radius, float(hole["rim_lo"]), True)
        if "floor_r" in hole:
            v += quarter_ring_volume(radius, float(hole["floor_r"]), False)
    for slot in feats.get("slots") or ():
        radius = float(slot["w"]) / 2.0
        v -= (float(slot["len"]) * float(slot["w"]) + math.pi * radius * radius) * float(b["t"])
    return v


def _box_explains_mesh(m: trimesh.Trimesh, feats: dict) -> bool:
    """False when the box solid is a different shape from the mesh.

    An L or a notched plate still fits inside a bounding box, and
    :func:`extract` will return that box. The volume gap is how that case is
    told apart from a real box, pad, hole, or boss.
    """
    mv = abs(float(m.volume))
    if mv < 1e-6:
        return True
    rel = abs(_feature_volume(feats) - mv) / mv
    if m.is_watertight:
        return rel <= 0.01
    return rel <= 0.25


# Box failures that only mean "this is not a rectangular box". An extruded
# outline is allowed to have those, so the extrusion error replaces them.
_BOX_SUPERSEDED = (
    "plane that is not axis-aligned",
    "pocket or step below the top face",
    "downward-facing face above the bottom",
    "no side faces",
    "no bottom face",
    "no top face",
    "overhang",
)


def _flatten_planar_spheres(m, regions, tol):
    """Read a sphere as a plane when the patch really is an axis-aligned face.

    Noise on a flat face (a counterbore floor, a shallow step) can beat the
    plane fit and come back as a large sphere. A real spherical cap does not
    put most of its vertices on one axis-aligned plane, so it is left as a
    sphere and still refused by the vocabulary.
    """
    band = max(6.0 * tol, 0.25)
    res_lim = max(4.0 * tol, 0.15)
    for r in regions:
        if not isinstance(r.surface, Sphere) or r.area <= 0 or len(r.faces) == 0:
            continue
        verts = m.vertices[np.unique(m.faces[r.faces])]
        if len(verts) < 12:
            continue
        best = None
        for axis in range(3):
            coord = verts[:, axis]
            med = float(np.median(coord))
            mask = np.abs(coord - med) <= band
            if float(mask.mean()) < 0.7:
                continue
            inl = verts[mask]
            level = float(np.median(inl[:, axis]))
            if float(np.percentile(np.abs(inl[:, axis] - level), 95)) > res_lim:
                continue
            span = max(float(inl[:, i].max() - inl[:, i].min())
                       for i in range(3) if i != axis)
            if span < max(8.0 * tol, 1.0):
                continue
            frac = float(mask.mean())
            if best is None or frac > best[0]:
                best = (frac, axis, level)
        if best is None:
            continue
        _frac, axis, level = best
        sign = 1.0 if float(np.median(m.face_normals[r.faces][:, axis])) >= 0.0 else -1.0
        normal = np.zeros(3)
        normal[axis] = sign
        r.surface = Plane(normal, sign * level)


def bracket_volume(base_x: float, base_y: float, leg_x: float, leg_y: float,
                   height: float, radius: float) -> float:
    """Volume of an L extruded along Z with every vertical edge filleted.

    Five convex corners each lose a square minus a quarter circle. The inner
    corner is concave: its centre sits in the notch and that quarter circle is
    added, so the net change is four corners removed.
    """
    sharp = float(base_x) * float(leg_y) + float(leg_x) * (float(base_y) - float(leg_y))
    corner = float(radius) ** 2 * (1.0 - math.pi / 4.0)
    return (sharp - 4.0 * corner) * float(height)


def bracket_profile_exprs(sx: int, sy: int) -> list[tuple[str, str]]:
    """Sharp L profile as ``(x, y)`` expressions, heel at the origin, CCW.

    ``sx`` and ``sy`` are +1 or -1, the direction each leg runs from the heel.
    The parameters are ``base_x``, ``base_y``, ``leg_x``, and ``leg_y``.
    """
    def term(sign: int, name: str) -> str:
        return name if sign > 0 else f"-{name}"

    pts = [
        ("0", "0"),
        (term(sx, "base_x"), "0"),
        (term(sx, "base_x"), term(sy, "leg_y")),
        (term(sx, "leg_x"), term(sy, "leg_y")),
        (term(sx, "leg_x"), term(sy, "base_y")),
        ("0", term(sy, "base_y")),
    ]
    if sx * sy < 0:
        pts = list(reversed(pts))
    return pts


def iter_bracket_holes(feats: dict):
    """``(name, hole)`` for mounting holes, foot first, then the upright."""
    b = feats["bracket"]
    for i, hole in enumerate(b.get("base_holes") or (), 1):
        yield f"base_hole_{i}", hole
    for i, hole in enumerate(b.get("leg_holes") or (), 1):
        yield f"leg_hole_{i}", hole


def bracket_parameters(feats: dict, tol: float = 0.01) -> list[tuple[str, float, str]]:
    """``(name, value, note)`` for a filleted L-bracket."""
    from .patterns import member_pattern_parameters

    b = feats["bracket"]
    rows = [
        ("base_x", b["x"], "outer length along X"),
        ("base_y", b["y"], "outer length along Y"),
        ("base_t", b["t"], "height along Z"),
        ("leg_x", b["leg_x"], "upright thickness along X"),
        ("leg_y", b["leg_y"], "foot thickness along Y"),
        ("fillet_r", b["fillet_r"], "vertical edge radius, including the inner corner"),
        ("base_x0", b["x0"], "heel position along X"),
        ("base_y0", b["y0"], "heel position along Y"),
        ("base_z0", b["z0"], "bottom face height"),
    ]
    for member in ("base", "leg"):
        rows.extend(member_pattern_parameters(b, member))
        where = "foot" if member == "base" else "upright"
        stem = f"{member}_hole"
        for i, hole in enumerate(b.get(f"{member}_holes") or (), 1):
            rows.append((f"{stem}_{i}_d", hole["d"], f"diameter, through the {where}"))
            rows.append((f"{stem}_{i}_x", hole["x"], ""))
            rows.append((f"{stem}_{i}_y", hole["y"], ""))
    # A pattern direction is rounded like a length. Snapping it to the nearest
    # half degree would slide the far holes off the mesh.
    def _tol_for(name: str) -> float:
        if ("angle" in name or "deg" in name) and "_row_" not in name and "_grid_" not in name:
            return 0.5
        return tol

    return [(name, _nice(val, _tol_for(name)), note) for name, val, note in rows]


def _cluster_walls(items, eps: float) -> list[dict]:
    """Group wall fragments that share a coordinate. ``items`` are ``(coord, lo, hi)``."""
    groups: list[dict] = []
    for coord, lo, hi in sorted(items, key=lambda t: t[0]):
        if groups and abs(coord - groups[-1]["coord"]) <= eps:
            g = groups[-1]
            n = g["n"] + 1
            g["coord"] = (g["coord"] * g["n"] + coord) / n
            g["n"] = n
            g["lo"] = min(g["lo"], lo)
            g["hi"] = max(g["hi"], hi)
        else:
            groups.append({"coord": float(coord), "lo": float(lo), "hi": float(hi), "n": 1})
    return groups


def _axis_ok(normal) -> int | None:
    n = np.abs(np.asarray(normal, float))
    ln = float(np.linalg.norm(n))
    if ln < 1e-12 or float(n.max()) / ln < math.cos(math.radians(8.0)):
        return None
    return int(np.argmax(n))


def _split_bracket_cylinders(m, cyls):
    """Six vertical fillets, plus Z cylinders that could be through-holes.

    Returns ``(None, None)`` when a cylinder is neither a fillet arc nor a
    hole wall. The caller still requires exactly six fillets.
    """
    fillets = []
    holes = []
    for r in cyls:
        s = r.surface
        axis = np.asarray(s.axis, float)
        axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
        span = _partial_span(m, s, r.faces)
        if span is None:
            return None, None
        vertical = abs(abs(float(axis @ Z)) - 1.0) <= 0.02
        if vertical and np.radians(50.0) <= span <= np.radians(140.0):
            fillets.append(r)
            continue
        if _loose_axis(s.axis) == 2 and span > np.radians(140.0):
            holes.append(r)
            continue
        return None, None
    if len(fillets) != 6:
        return None, None
    return fillets, holes


def _ray_inside_poly(poly, x: float, y: float) -> bool:
    inside = False
    j = len(poly) - 1
    for i, (xi, yi) in enumerate(poly):
        xj, yj = poly[j]
        if (yi > y) != (yj > y):
            cross = (xj - xi) * (y - yi) / ((yj - yi) or 1e-30) + xi
            if x < cross:
                inside = not inside
        j = i
    return inside


def _bracket_member(lx, ly, radius, base_x, base_y, leg_x, leg_y, med, tol) -> str | None:
    """``base`` or ``leg`` when the circle sits in that member and clears the rounds.

    A circle in the corner square is inside both members. That hole is named
    with the foot. A circle that crosses from the exclusive foot into the
    exclusive upright is neither.
    """
    poly = (
        (0.0, 0.0),
        (float(base_x), 0.0),
        (float(base_x), float(leg_y)),
        (float(leg_x), float(leg_y)),
        (float(leg_x), float(base_y)),
        (0.0, float(base_y)),
    )
    pad = max(2.0 * tol, 0.15)
    test_r = max(0.0, float(radius) - pad)
    if not _ray_inside_poly(poly, lx, ly):
        return None
    for k in range(24):
        ang = 2.0 * math.pi * k / 24.0
        if not _ray_inside_poly(poly, lx + test_r * math.cos(ang), ly + test_r * math.sin(ang)):
            return None
    # Outer fillets sit inside the sharp L. The inner fillet's centre is in
    # the notch, so a circle inside the sharp outline does not meet it.
    need = float(med) + float(radius) - max(4.0 * tol, 0.35)
    outers = (
        (med, med),
        (base_x - med, med),
        (base_x - med, leg_y - med),
        (med, base_y - med),
        (leg_x - med, base_y - med),
    )
    for fx, fy in outers:
        if math.hypot(lx - fx, ly - fy) < need:
            return None

    def in_rect(x0, x1, y0, y1) -> bool:
        return (x0 - pad + radius <= lx <= x1 + pad - radius
                and y0 - pad + radius <= ly <= y1 + pad - radius)

    in_foot = in_rect(0.0, base_x, 0.0, leg_y)
    in_upright = in_rect(0.0, leg_x, 0.0, base_y)
    if in_foot and not in_upright:
        return "base"
    if in_upright and not in_foot:
        return "leg"
    if in_foot and in_upright:
        return "base"
    return None


def _refit_hole_circle(wall, point, radius, tol):
    """Centre and radius from wall vertices, kept only when close to the fit."""
    if len(wall) < 40:
        return np.asarray(point, float), float(radius)
    x, y = wall[:, 0], wall[:, 1]
    fit = np.c_[2.0 * x, 2.0 * y, np.ones(len(x))]
    sol, *_ = np.linalg.lstsq(fit, x * x + y * y, rcond=None)
    cx, cy, k = float(sol[0]), float(sol[1]), float(sol[2])
    r2 = k + cx * cx + cy * cy
    if r2 <= 0.0:
        return np.asarray(point, float), float(radius)
    rad = math.sqrt(r2)
    if abs(rad - float(radius)) > max(4.0 * tol, 0.25):
        return np.asarray(point, float), float(radius)
    if math.hypot(cx - float(point[0]), cy - float(point[1])) > max(4.0 * tol, 0.4):
        return np.asarray(point, float), float(radius)
    return np.array([cx, cy, 0.0]), rad


def _group_hole_regs(regs, tol):
    """Coaxial hole fragments. A noisy wall often arrives as more than one region."""
    groups = []
    gtol = max(5.0 * tol, 0.2)
    rtol = max(4.0 * tol, 0.25)
    for region in regs:
        s = region.surface
        point = np.asarray(s.point, float)
        radius = float(s.radius)
        for g in groups:
            if abs(radius - g["radius"]) > rtol:
                continue
            if _perp_dist(2, g["point"], point) > gtol:
                continue
            w = g["area"] + float(region.area)
            g["point"] = (g["point"] * g["area"] + point * float(region.area)) / w
            g["radius"] = (g["radius"] * g["area"] + radius * float(region.area)) / w
            g["area"] = w
            g["regions"].append(region)
            break
        else:
            groups.append({
                "point": point.copy(), "radius": radius,
                "area": float(region.area), "regions": [region],
            })
    return groups


def _as_bracket_hole(m, group, z0, z1, tol):
    """One coaxial group as a Z through-hole, or None."""
    point = np.asarray(group["point"], float)
    radius = float(group["radius"])
    if radius <= max(2.0 * tol, 0.2):
        return None
    faces = np.concatenate([r.faces for r in group["regions"]])
    pts = m.vertices[np.unique(m.faces[faces])]
    on = _on_cylinder(pts, point, radius, 2, tol)
    if int(on.sum()) < 12:
        return None
    wall = pts[on]
    along = wall[:, 2]
    eps = max(6.0 * tol, 0.75)
    if float(along.min()) > float(z0) + eps or float(along.max()) < float(z1) - eps:
        return None
    gap = _angular_gap(np.array([0.0, 0.0, 1.0]), point, wall)
    if 2.0 * math.pi - gap < math.radians(270.0):
        return None
    centres = m.triangles_center[faces]
    near = _on_cylinder(centres, point, radius, 2, tol)
    if int(near.sum()) < 8:
        return None
    radial = centres[near] - point
    radial[:, 2] = 0.0
    ln = np.linalg.norm(radial, axis=1, keepdims=True)
    radial = radial / np.maximum(ln, 1e-12)
    if float(np.mean(np.sum(m.face_normals[faces][near] * radial, axis=1))) > -0.25:
        return None
    point, radius = _refit_hole_circle(wall, point, radius, tol)
    return {"d": 2.0 * float(radius), "x": float(point[0]), "y": float(point[1])}


def _bracket_mounting_holes(m, regs, x0, y0, z0, z1, base_x, base_y, leg_x, leg_y,
                            med, sx, sy, tol):
    """Through-holes in the foot and the upright, or None if one does not fit.

    Overlapping holes are refused: the volume check would subtract the shared
    cylinder twice and the script would not match the mesh.
    """
    found = []
    for g in _group_hole_regs(regs, tol):
        hole = _as_bracket_hole(m, g, z0, z1, tol)
        if hole is None:
            return None
        lx = sx * (hole["x"] - x0)
        ly = sy * (hole["y"] - y0)
        member = _bracket_member(
            lx, ly, hole["d"] / 2.0, base_x, base_y, leg_x, leg_y, med, tol)
        if member is None:
            return None
        found.append((member, hole))
    slack = max(4.0 * tol, 0.35)
    for i, (_ma, a) in enumerate(found):
        for _mb, b in found[i + 1:]:
            if math.hypot(a["x"] - b["x"], a["y"] - b["y"]) < (a["d"] + b["d"]) / 2.0 - slack:
                return None
    base, leg = [], []
    for member, hole in found:
        row = {"d": hole["d"], "x": hole["x"], "y": hole["y"]}
        (base if member == "base" else leg).append(row)
    base.sort(key=lambda h: (h["x"], h["y"]))
    leg.sort(key=lambda h: (h["x"], h["y"]))
    if not base and not leg:
        return None
    return base, leg


def _filleted_bracket(m: trimesh.Trimesh, regions, tol: float) -> dict | None:
    """An L extruded along Z, every vertical edge filleted at one radius.

    The inner corner is included. A round through-hole along Z is kept when
    its circle lies in the foot or the upright and clears those rounds.
    Anything else returns None so the other vocabularies keep the part they
    already explain.
    """
    planes = []
    cyls = []
    for r in regions:
        s = r.surface
        if s is None or not len(r.faces):
            return None
        if isinstance(s, Plane):
            planes.append(r)
        elif isinstance(s, Cylinder):
            cyls.append(r)
        else:
            return None
    fillets, hole_regs = _split_bracket_cylinders(m, cyls)
    if fillets is None:
        return None
    radii = []
    centers = []
    z_ends = []
    for r in fillets:
        s = r.surface
        axis = np.asarray(s.axis, float)
        axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
        if abs(abs(float(axis @ Z)) - 1.0) > 0.02:
            return None
        span = _partial_span(m, s, r.faces)
        if span is None or not np.radians(50.0) <= span <= np.radians(140.0):
            return None
        radii.append(float(s.radius))
        pt = np.asarray(s.point, float)
        centers.append((float(pt[0]), float(pt[1])))
        vz = _verts(m, r)[:, 2]
        z_ends.append((float(vz.min()), float(vz.max())))
    med = float(np.median(radii))
    if med <= max(2.0 * tol, 0.05):
        return None
    if max(abs(r - med) for r in radii) > max(4.0 * tol, 0.25):
        return None

    walls = {0: [], 1: [], 2: []}
    for r in planes:
        axis = _axis_ok(r.surface.normal)
        if axis is None:
            return None
        _i, _sg, coord = _plane_pos(r.surface)
        if axis != _i:
            return None
        free = 1 if axis == 0 else 0 if axis == 1 else None
        if free is None:
            walls[2].append((coord, coord, coord))
            continue
        verts = _verts(m, r)
        walls[axis].append((coord, float(verts[:, free].min()), float(verts[:, free].max())))
    eps = max(4.0 * tol, 0.3)
    xs = _cluster_walls(walls[0], eps)
    ys = _cluster_walls(walls[1], eps)
    zs = _cluster_walls(walls[2], eps)
    if len(xs) != 3 or len(ys) != 3 or len(zs) != 2:
        return None

    def _orient(groups):
        a, mid, b = groups
        la, lb = a["hi"] - a["lo"], b["hi"] - b["lo"]
        if la >= lb:
            heel, toe, sign = a, b, 1
            long, short = la, lb
        else:
            heel, toe, sign = b, a, -1
            long, short = lb, la
        # The upright's outer face runs the long way. The foot's end is only
        # the leg thickness. A near-rectangle does not separate that cleanly.
        if long < 1.8 * short or long - short < max(8.0 * tol, 2.0):
            return None
        return heel, mid, toe, sign

    x_orient = _orient(xs)
    y_orient = _orient(ys)
    if x_orient is None or y_orient is None:
        return None
    x_heel, x_mid, x_toe, sx = x_orient
    y_heel, y_mid, y_toe, sy = y_orient
    base_x = abs(x_toe["coord"] - x_heel["coord"])
    base_y = abs(y_toe["coord"] - y_heel["coord"])
    leg_x = abs(x_mid["coord"] - x_heel["coord"])
    leg_y = abs(y_mid["coord"] - y_heel["coord"])
    if not (med < leg_x < base_x - med and med < leg_y < base_y - med):
        return None
    z0 = min(zs[0]["coord"], zs[1]["coord"])
    z1 = max(zs[0]["coord"], zs[1]["coord"])
    height = z1 - z0
    if height <= max(8.0 * tol, 1.0):
        return None
    zeps = max(6.0 * tol, 0.5)
    for lo, hi in z_ends:
        if abs(lo - z0) > zeps or abs(hi - z1) > zeps:
            return None

    # Centres: convex corners inset into the material, inner corner in the notch.
    expected = [
        (x_heel["coord"] + sx * med, y_heel["coord"] + sy * med),
        (x_heel["coord"] + sx * (base_x - med), y_heel["coord"] + sy * med),
        (x_heel["coord"] + sx * (base_x - med), y_heel["coord"] + sy * (leg_y - med)),
        (x_heel["coord"] + sx * med, y_heel["coord"] + sy * (base_y - med)),
        (x_heel["coord"] + sx * (leg_x - med), y_heel["coord"] + sy * (base_y - med)),
        (x_heel["coord"] + sx * (leg_x + med), y_heel["coord"] + sy * (leg_y + med)),
    ]
    ceps = max(8.0 * tol, 0.45)
    unused = list(centers)
    for ex, ey in expected:
        best = min(range(len(unused)), key=lambda i: math.hypot(unused[i][0] - ex, unused[i][1] - ey))
        if math.hypot(unused[best][0] - ex, unused[best][1] - ey) > ceps:
            return None
        unused.pop(best)
    base_holes: list[dict] = []
    leg_holes: list[dict] = []
    if hole_regs:
        mounted = _bracket_mounting_holes(
            m, hole_regs, x_heel["coord"], y_heel["coord"], z0, z1,
            base_x, base_y, leg_x, leg_y, med, sx, sy, tol)
        if mounted is None:
            return None
        base_holes, leg_holes = mounted
    vol = bracket_volume(base_x, base_y, leg_x, leg_y, height, med)
    for hole in base_holes + leg_holes:
        vol -= math.pi * (float(hole["d"]) / 2.0) ** 2 * height
    mv = abs(float(m.volume))
    if mv < 1e-6 or abs(vol - mv) / mv > 0.02:
        return None
    bracket = {
        "x0": float(x_heel["coord"]), "y0": float(y_heel["coord"]), "z0": float(z0),
        "x": float(base_x), "y": float(base_y), "t": float(height),
        "leg_x": float(leg_x), "leg_y": float(leg_y),
        "fillet_r": float(med), "sx": int(sx), "sy": int(sy),
    }
    if base_holes:
        bracket["base_holes"] = base_holes
    if leg_holes:
        bracket["leg_holes"] = leg_holes
    return {"bracket": bracket}


def recognize(m: trimesh.Trimesh, regions, tol: float) -> dict:
    """Box vocabulary when it matches, else one extrusion, else a small CSG."""
    from .csg import extract_csg
    from .extrude import extract_extrusion

    # A fully filleted box is six planes, twelve cylinders and eight spheres.
    # Claim it before the box extractor, which refuses every sphere. A
    # through-hole or a top boss on that box is the same solid plus that
    # feature; the corner spheres stay part of the fillet. The plain box is
    # kept when it already accounts for the mesh, and when the extra feature
    # does not, so the script gate still names the deviation.
    rounded = _edge_fillet_box(m, regions, tol)
    if rounded is not None and _box_explains_mesh(m, rounded):
        return rounded
    holed = _edge_fillet_with_holes(m, regions, tol)
    if holed is not None:
        return holed
    if rounded is not None:
        return rounded
    # An L with every vertical edge filleted, including the inner corner, is
    # not a rectangular box. Claim it before the box extractor refuses the
    # six fillets. Through-holes in the foot or the upright stay on this
    # claim. A plate, a plain extrusion, and a fully filleted box do not have
    # this wall pattern, so they stay on the path they already use.
    bracket = _filleted_bracket(m, regions, tol)
    if bracket is not None:
        return bracket
    _flatten_planar_spheres(m, regions, tol)

    box_err = None
    feats = None
    try:
        feats = extract(m, regions, tol)
    except Unsupported as e:
        box_err = e
    if feats is not None and _box_explains_mesh(m, feats):
        return feats
    ext = None
    try:
        return extract_extrusion(m, regions, tol)
    except Unsupported as e:
        ext = e
    # A hard box failure (a sphere, an angled axis, an unfitted patch) is the
    # reason to report. "Not a rectangular box" is not: the extrusion or the
    # CSG may still explain the part, and its message names the real obstacle.
    hard_box = box_err is not None and not any(s in str(box_err) for s in _BOX_SUPERSEDED)
    if hard_box:
        raise box_err
    try:
        return extract_csg(m, regions, tol)
    except Unsupported:
        # extract_csg names the obstacle (an unfitted patch, a rounded edge,
        # a recess that is not a straight cut, a measured deviation).
        raise


def chamfer_is_square(d: dict, tol: float) -> bool:
    """True when the top chamfer is 45°, so one length names both width and height."""
    return abs(d["chamfer"] - d["chamfer_w"]) <= max(2 * tol, 0.02)


def named_parameters(feats: dict, tol: float = 0.01) -> list[tuple[str, float, str]]:
    """``(name, value, note)`` in the order every emitter prints them."""
    from .patterns import box_pattern_parameters

    P: list[tuple[str, float, str]] = []

    def p(name, val, note=""):
        t = 0.5 if ("angle" in name or "deg" in name) else tol
        P.append((name, _nice(val, t), note))
        return name

    b = feats["base"]
    p("base_x", b["x"], "base length (X)")
    p("base_y", b["y"], "base width (Y)")
    p("base_t", b["t"], "base thickness (Z)")
    p("base_cx", b["cx"])
    p("base_cy", b["cy"])
    p("base_z0", b["z0"], "bottom face height")
    if "edge_fillet_r" in feats:
        p("edge_fillet_r", feats["edge_fillet_r"], "radius of every edge")
    if "fillet_r" in feats:
        p("fillet_r", feats["fillet_r"], "vertical corner radius")
    for i, d in enumerate(feats["pads"], 1):
        n = f"pad_{i}"
        for k, note in (("w", "X size"), ("l", "Y size"), ("h", "height above top face")):
            p(f"{n}_{k}", d[k], note)
        p(f"{n}_x", d["cx"])
        p(f"{n}_y", d["cy"])
    for i, d in enumerate(feats["bosses"], 1):
        n = f"boss_{i}"
        ax = d.get("axis", "Z")
        if ax == "Z":
            p(f"{n}_d", d["d"], "diameter at the top face")
            p(f"{n}_h", d["h"], "height")
            if "fillet" in d:
                p(f"{n}_fillet", d["fillet"], "fillet at the base")
            p(f"{n}_x", d["x"])
            p(f"{n}_y", d["y"])
            if "draft_deg" in d:
                p(f"{n}_draft_deg", d["draft_deg"], "taper angle")
            elif "chamfer" in d:
                if chamfer_is_square(d, tol):
                    p(f"{n}_chamfer", d["chamfer"], "45° top-edge chamfer")
                else:
                    p(f"{n}_chamfer_h", d["chamfer"])
                    p(f"{n}_chamfer_w", d["chamfer_w"])
        else:
            face = _face_name(ax, d["from"])
            p(f"{n}_d", d["d"], "diameter")
            p(f"{n}_h", d["h"], f"height from the {face} face")
            if "fillet" in d:
                p(f"{n}_fillet", d["fillet"], f"fillet at the {face} face")
            for coord in "xyz":
                if coord != ax.lower() and coord in d:
                    p(f"{n}_{coord}", d[coord])
    for name, val, note in box_pattern_parameters(feats):
        # Direction uses the length tolerance. Half a degree would move the
        # last hole of a long row by more than the deviation gate.
        if name.endswith("_angle"):
            P.append((name, _nice(val, tol), note))
        else:
            p(name, val, note)
    for i, d in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        ax = d.get("axis", "Z")
        if ax == "Z":
            p(f"{n}_d", d["d"], "diameter")
            p(f"{n}_x", d["x"])
            p(f"{n}_y", d["y"])
            if "depth" in d:
                p(f"{n}_depth", d["depth"], "blind depth from top face")
            if "floor_r" in d:
                p(f"{n}_floor_r", d["floor_r"], "fillet where the floor meets the wall")
            if "rim_r" in d:
                p(f"{n}_rim_r", d["rim_r"], "round on the top rim")
            if "rim_lo" in d:
                p(f"{n}_bot_rim_r", d["rim_lo"], "round on the bottom rim")
            if "csk_d" in d:
                p(f"{n}_csk_d", d["csk_d"], "countersink diameter")
                p(f"{n}_csk_angle", d["csk_angle"], "countersink angle")
            if "cbore_d" in d:
                face = "top face" if d.get("from", 1) > 0 else "bottom face"
                p(f"{n}_cbore_d", d["cbore_d"], "counterbore diameter")
                p(f"{n}_cbore_depth", d["cbore_depth"], f"counterbore depth from the {face}")
        else:
            p(f"{n}_d", d["d"], f"diameter (axis {ax})")
            for coord in "xyz":
                if coord != ax.lower() and coord in d:
                    p(f"{n}_{coord}", d[coord])
            if "depth" in d:
                p(f"{n}_depth", d["depth"], f"blind depth from the {_face_name(ax, d['from'])} face")
            if "floor_r" in d:
                p(f"{n}_floor_r", d["floor_r"], "fillet where the floor meets the wall")
            if "rim_r" in d:
                face = _face_name(ax, 1)
                p(f"{n}_rim_r", d["rim_r"], f"round on the {face} rim")
            if "rim_lo" in d:
                face = _face_name(ax, -1)
                p(f"{n}_bot_rim_r", d["rim_lo"], f"round on the {face} rim")
            if "cbore_d" in d:
                p(f"{n}_cbore_d", d["cbore_d"], "counterbore diameter")
                p(f"{n}_cbore_depth", d["cbore_depth"],
                  f"counterbore depth from the {_face_name(ax, d['from'])} face")
    for i, slot in enumerate(feats.get("slots") or (), 1):
        n = f"slot_{i}"
        p(f"{n}_len", slot["len"], "center-to-center length")
        p(f"{n}_w", slot["w"], "width, the end diameter")
        p(f"{n}_x", slot["x"], "center")
        p(f"{n}_y", slot["y"], "center")
        # Same tolerance as a length. Half a degree would slide the round ends.
        P.append((f"{n}_angle", _nice(float(slot["angle"]), tol),
                  "direction of the slot axis, degrees"))
    return P


def deviation_problem(q: dict, tol: float) -> str | None:
    """None when a generated solid is close enough to the mesh it came from."""
    lim = max(4 * tol, 0.1)
    if q["max_deviation"] > lim or q["volume_error_pct"] > 1.0:
        return (f"generated script deviates by {q['max_deviation']:.3f} "
                f"(limit {lim:.3f}), volume error {q['volume_error_pct']:.2f}%")
    return None
