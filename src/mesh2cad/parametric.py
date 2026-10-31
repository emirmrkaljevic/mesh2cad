"""Shared vocabulary for parametric CadQuery, OpenSCAD, and FreeCAD export.

Supported parts:

* a Z-up box or plate, rectangular pads on its top face, holes along X, Y, or
  Z (through, blind; countersunk when the axis is +Z), bosses along X, Y, or Z
  (straight; drafted or chamfered when the axis is +Z), and one constant
  radius on all four vertical corners of the base;
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

from .fit import Cone, Cylinder, Plane, Sphere
from .segment import _partial_span

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
    for r in regions:
        if not isinstance(r.surface, (Plane, Cylinder, Cone)):
            raise Unsupported(f"{r.surface.kind} faces are not supported yet")
    for r in planes:
        n = np.abs(r.surface.normal)
        if n.max() < math.cos(math.radians(2)):
            raise Unsupported("plane that is not axis-aligned")
    pos = [(r, *_plane_pos(r.surface)) for r in planes]
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
    fillet_regs = _vertical_fillets(m, regions)
    fillet_ids = {id(r) for r in fillet_regs}
    stacks = _stacks(m, [r for r in regions if id(r) not in fillet_ids], tol, plane_coords)

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
            feats["bosses"].append(b)
    if fillet_regs:
        feats["fillet_r"] = _base_fillet_radius(m, fillet_regs, lo, hi, zmin, ztop, tol)
    return feats


def _feature_volume(feats: dict) -> float:
    """Analytic volume of a box-vocabulary feature dict."""
    b = feats["base"]
    v = float(b["x"]) * float(b["y"]) * float(b["t"])
    if "fillet_r" in feats:
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


def recognize(m: trimesh.Trimesh, regions, tol: float) -> dict:
    """Box vocabulary when it matches, else one extrusion, else a small CSG."""
    from .csg import extract_csg
    from .extrude import extract_extrusion

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
            for coord in "xyz":
                if coord != ax.lower() and coord in d:
                    p(f"{n}_{coord}", d[coord])
    for i, d in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        ax = d.get("axis", "Z")
        if ax == "Z":
            p(f"{n}_d", d["d"], "diameter")
            p(f"{n}_x", d["x"])
            p(f"{n}_y", d["y"])
            if "depth" in d:
                p(f"{n}_depth", d["depth"], "blind depth from top face")
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
            if "cbore_d" in d:
                p(f"{n}_cbore_d", d["cbore_d"], "counterbore diameter")
                p(f"{n}_cbore_depth", d["cbore_depth"],
                  f"counterbore depth from the {_face_name(ax, d['from'])} face")
    return P


def deviation_problem(q: dict, tol: float) -> str | None:
    """None when a generated solid is close enough to the mesh it came from."""
    lim = max(4 * tol, 0.1)
    if q["max_deviation"] > lim or q["volume_error_pct"] > 1.0:
        return (f"generated script deviates by {q['max_deviation']:.3f} "
                f"(limit {lim:.3f}), volume error {q['volume_error_pct']:.2f}%")
    return None
