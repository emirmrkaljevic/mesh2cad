"""Straight obround (stadium) through-slots along Z, for the box vocabulary.

A slot is two round ends of one diameter and two straight walls between them.
``slot_N_len`` is the distance between the end centres, ``slot_N_w`` is the
width (that diameter), ``slot_N_x`` and ``slot_N_y`` are the centre, and
``slot_N_angle`` is the direction of the long axis in degrees. Axis-aligned
slots are 0 (along +X) and 90 (along +Y). Another angle is kept when the same
two ends and two walls still fit; the usual deviation gate decides whether
the script is written.

The round ends are often a partial cylinder. On a clean tessellation they are
a ring of flat facets instead, the same way a broken-out counterbore is. A
pair is claimed only when both ends, both walls, and a through span along Z
agree. Anything else is left alone, so a rectangular opening stays a cutout
and a round hole stays ``hole_N_*``.
"""
from __future__ import annotations

import math

import numpy as np

from .fit import Cylinder, Plane
from .segment import _partial_span, partial_arc_cylinders

_Z = np.array([0.0, 0.0, 1.0])


def _unit_plane(plane: Plane):
    raw = np.asarray(plane.normal, float)
    ln = float(np.linalg.norm(raw))
    if ln < 1e-12:
        return None
    return raw / ln, float(plane.d) / ln


def _plate_z(m, regions, tol: float):
    """``(zmin, ztop)`` from the horizontal faces, or None."""
    levels = []
    floor = max(0.012 * float(m.area), 2.0)
    for region in regions:
        surf = region.surface
        if not isinstance(surf, Plane) or region.area < floor:
            continue
        got = _unit_plane(surf)
        if got is None:
            continue
        normal, level = got
        if abs(float(normal[2])) < math.cos(math.radians(12.0)):
            continue
        levels.append(level / float(normal[2]))
    if len(levels) < 2:
        return None
    z0, z1 = float(min(levels)), float(max(levels))
    if z1 - z0 <= max(8.0 * tol, 1.0):
        return None
    return z0, z1


def _inward_bulge(m, faces, point, radius: float):
    """Unit direction from the axis toward the arc, when the wall faces the axis.

    None when the patch faces outward (a boss or a corner fillet) or the arc
    has no clear side.
    """
    centres = m.triangles_center[faces]
    radial = centres - np.asarray(point, float)
    radial[:, 2] = 0.0
    norm = np.linalg.norm(radial, axis=1)
    ok = norm > max(0.15 * radius, 1e-4)
    if int(ok.sum()) < 6:
        return None
    dirs = radial[ok] / norm[ok, None]
    facing = float(np.median(np.sum(m.face_normals[faces][ok] * dirs, axis=1)))
    if facing > -0.2:
        return None
    bulge = dirs.mean(axis=0)[:2]
    scale = float(np.linalg.norm(bulge))
    if scale < 0.12:
        return None
    return bulge / scale


def _as_cap(m, region, tol: float, z0: float, z1: float, owners: dict):
    """One round end, or None."""
    surf = region.surface
    if not isinstance(surf, Cylinder) or surf.radius <= max(2.0 * tol, 0.3):
        return None
    axis = np.asarray(surf.axis, float)
    axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
    if abs(abs(float(axis @ _Z)) - 1.0) > math.sin(math.radians(12.0)):
        return None
    faces = np.asarray(region.faces, int)
    span = _partial_span(m, surf, faces)
    if span is None or not math.radians(70.0) <= float(span) <= math.radians(250.0):
        return None
    point = np.asarray(surf.point, float)
    bulge = _inward_bulge(m, faces, point, float(surf.radius))
    if bulge is None:
        return None
    verts = m.vertices[np.unique(m.faces[faces])]
    zeps = max(6.0 * tol, 0.55)
    if float(verts[:, 2].min()) > z0 + zeps or float(verts[:, 2].max()) < z1 - zeps:
        return None
    src, planes = set(), set()
    for face in faces:
        owner = owners.get(int(face))
        if owner is None:
            continue
        if owner is region and isinstance(surf, Cylinder):
            src.add(id(owner))
        elif isinstance(owner.surface, Plane):
            planes.add(id(owner))
        elif isinstance(owner.surface, Cylinder):
            src.add(id(owner))
    return {
        "x": float(point[0]),
        "y": float(point[1]),
        "r": float(surf.radius),
        "bulge": bulge,
        "faces": faces,
        "area": float(region.area) if region.area else float(m.area_faces[faces].sum()),
        "src": src,
        "planes": planes,
    }


def _owners(regions) -> dict:
    out = {}
    for region in regions:
        for face in np.asarray(region.faces, int).ravel():
            out[int(face)] = region
    return out


def _merge_caps(caps: list[dict], tol: float) -> list[dict]:
    """Coaxial fragments of one end, area-weighted."""
    gtol = max(5.0 * tol, 0.45)
    rtol = max(4.0 * tol, 0.35)
    groups: list[dict] = []
    for cap in caps:
        for group in groups:
            if abs(cap["r"] - group["r"]) > rtol:
                continue
            if math.hypot(cap["x"] - group["x"], cap["y"] - group["y"]) > gtol:
                continue
            weight = group["area"] + cap["area"]
            group["x"] = (group["x"] * group["area"] + cap["x"] * cap["area"]) / weight
            group["y"] = (group["y"] * group["area"] + cap["y"] * cap["area"]) / weight
            group["r"] = (group["r"] * group["area"] + cap["r"] * cap["area"]) / weight
            group["bulge"] = group["bulge"] * group["area"] + cap["bulge"] * cap["area"]
            scale = float(np.linalg.norm(group["bulge"]))
            if scale > 1e-9:
                group["bulge"] = group["bulge"] / scale
            group["area"] = weight
            group["faces"] = np.concatenate([group["faces"], cap["faces"]])
            group["src"] |= cap["src"]
            group["planes"] |= cap["planes"]
            break
        else:
            groups.append({
                "x": cap["x"], "y": cap["y"], "r": cap["r"],
                "bulge": np.asarray(cap["bulge"], float),
                "faces": cap["faces"], "area": cap["area"],
                "src": set(cap["src"]), "planes": set(cap["planes"]),
            })
    return groups


def _caps(m, regions, tol: float, z0: float, z1: float) -> list[dict]:
    owners = _owners(regions)
    found = []
    for region in regions:
        cap = _as_cap(m, region, tol, z0, z1, owners)
        if cap is not None and cap["src"]:
            found.append(cap)
    if len(found) < 2:
        for arc in partial_arc_cylinders(m, regions, tol):
            cap = _as_cap(m, arc, tol, z0, z1, owners)
            if cap is not None:
                found.append(cap)
    return _merge_caps(found, tol)


def _wall_candidates(m, regions, direction, perp, mid, radius: float, tol: float):
    eps = max(4.0 * tol, 0.4)
    out = []
    for region in regions:
        surf = region.surface
        if not isinstance(surf, Plane) or region.area < max(6.0 * tol, 1.5):
            continue
        got = _unit_plane(surf)
        if got is None:
            continue
        normal, level = got
        if abs(float(normal[2])) > math.sin(math.radians(18.0)):
            continue
        horiz = np.array([normal[0], normal[1]], float)
        scale = float(np.linalg.norm(horiz))
        if scale < 0.8:
            continue
        horiz /= scale
        if abs(float(horiz @ direction)) > math.sin(math.radians(20.0)):
            continue
        if float(horiz @ perp) < 0.0:
            horiz = -horiz
            level = -level
            normal = -normal
        sep = float(horiz @ mid - level)
        if abs(abs(sep) - radius) > eps:
            continue
        verts = m.vertices[np.unique(m.faces[region.faces])]
        out.append({
            "sep": sep,
            "area": float(region.area),
            "region": region,
            "verts": verts,
        })
    return out


def _cluster_walls(cands, tol: float):
    if not cands:
        return []
    eps = max(4.0 * tol, 0.4)
    ordered = sorted(cands, key=lambda c: c["sep"])
    groups = []
    for cand in ordered:
        if groups and abs(cand["sep"] - groups[-1]["sep"]) <= eps:
            group = groups[-1]
            n = group["n"] + 1
            group["sep"] = (group["sep"] * group["n"] + cand["sep"]) / n
            group["n"] = n
            group["area"] += cand["area"]
            group["regions"].append(cand["region"])
            group["verts"].append(cand["verts"])
        else:
            groups.append({
                "sep": float(cand["sep"]),
                "n": 1,
                "area": float(cand["area"]),
                "regions": [cand["region"]],
                "verts": [cand["verts"]],
            })
    return groups


def _wall_ok(group, direction, mid, length: float, radius: float, z0: float, z1: float, tol: float) -> bool:
    verts = np.vstack(group["verts"])
    along = (verts[:, 0] - mid[0]) * direction[0] + (verts[:, 1] - mid[1]) * direction[1]
    lo, hi = np.percentile(along, [2.0, 98.0])
    margin = max(0.9 * radius, 1.0)
    half = 0.5 * length
    if float(lo) < -half - margin or float(hi) > half + margin:
        return False
    if float(hi - lo) < 0.45 * length:
        return False
    thick = z1 - z0
    if group["area"] < 0.22 * length * thick:
        return False
    zeps = max(6.0 * tol, 0.6)
    if float(verts[:, 2].min()) > z0 + zeps or float(verts[:, 2].max()) < z1 - zeps:
        return False
    return True


def _pair(m, a: dict, b: dict, regions, z0: float, z1: float, tol: float):
    dx, dy = b["x"] - a["x"], b["y"] - a["y"]
    length = float(math.hypot(dx, dy))
    radius = 0.5 * (a["r"] + b["r"])
    if length < max(1.0, 3.0 * tol):
        return None
    if abs(a["r"] - b["r"]) > max(4.0 * tol, 0.4):
        return None
    direction = np.array([dx / length, dy / length], float)
    # The arc sits on the outer side, so each bulge points away from the other end.
    if float(a["bulge"] @ direction) > -0.15:
        return None
    if float(b["bulge"] @ -direction) > -0.15:
        return None
    # Undirected: keep the axis in [0, 180) so a slot along -X is angle 0.
    if direction[0] < 0.0 or (abs(direction[0]) < 1e-9 and direction[1] < 0.0):
        direction = -direction
        a, b = b, a
        dx, dy = -dx, -dy
    perp = np.array([-direction[1], direction[0]], float)
    mid = np.array([(a["x"] + b["x"]) / 2.0, (a["y"] + b["y"]) / 2.0], float)
    groups = _cluster_walls(_wall_candidates(m, regions, direction, perp, mid, radius, tol), tol)
    pos = [g for g in groups if g["sep"] >= 0.0 and _wall_ok(g, direction, mid, length, radius, z0, z1, tol)]
    neg = [g for g in groups if g["sep"] < 0.0 and _wall_ok(g, direction, mid, length, radius, z0, z1, tol)]
    if not pos or not neg:
        return None
    right = min(pos, key=lambda g: abs(g["sep"] - radius))
    left = min(neg, key=lambda g: abs(-g["sep"] - radius))
    eps = max(4.0 * tol, 0.4)
    if abs(right["sep"] - radius) > eps or abs(-left["sep"] - radius) > eps:
        return None
    width = float(right["sep"] - left["sep"])
    if abs(width - 2.0 * radius) > max(6.0 * tol, 0.5):
        return None
    angle = math.degrees(math.atan2(float(direction[1]), float(direction[0])))
    if angle < 0.0:
        angle += 180.0
    if angle >= 180.0 - 1e-6:
        angle = 0.0
    # Snap to the axis when the ends already share a coordinate. A few tenths
    # of a millimetre of noise must not turn a horizontal slot into 1 degree,
    # and a real diagonal is left at the angle the ends measured.
    lateral = abs(float(a["y"] - b["y"])) if abs(dx) >= abs(dy) else abs(float(a["x"] - b["x"]))
    move = 0.5 * lateral
    if move <= max(2.0 * tol, 0.05) and lateral <= max(4.0 * tol, 0.25):
        if abs(dx) >= abs(dy):
            angle = 0.0
            length = abs(float(a["x"] - b["x"]))
        else:
            angle = 90.0
            length = abs(float(a["y"] - b["y"]))
    if length < max(1.0, 3.0 * tol):
        return None
    walls = right["regions"] + left["regions"]
    residual = abs(right["sep"] - radius) + abs(-left["sep"] - radius) + abs(a["r"] - b["r"])
    return {
        "len": float(length),
        "w": float(width),
        "x": float(mid[0]),
        "y": float(mid[1]),
        "angle": float(angle),
        "_residual": float(residual),
        "_caps": (a, b),
        "_walls": walls,
    }


def _nearby_planes(m, regions, cap: dict, tol: float) -> set[int]:
    """Facet planes of this end that the arc fit did not already consume.

    The tangent slivers where the semicircle meets the straight wall are still
    planes, and they are not axis-aligned. They belong to the end.
    """
    cyl = Cylinder(
        np.array([0.0, 0.0, 1.0]),
        np.array([cap["x"], cap["y"], 0.0]),
        float(cap["r"]),
    )
    thr = max(4.0 * tol, 0.15)
    cap_area = math.pi * float(cap["r"]) * max(float(cap.get("span_h", 1.0)), 1.0)
    found = set()
    for region in regions:
        if not isinstance(region.surface, Plane):
            continue
        if id(region) in cap["planes"]:
            found.add(id(region))
            continue
        if region.area > max(4.0 * cap_area, 80.0):
            continue
        verts = m.vertices[np.unique(m.faces[region.faces])]
        if len(verts) < 3:
            continue
        centre = verts.mean(axis=0)
        if math.hypot(float(centre[0] - cap["x"]), float(centre[1] - cap["y"])) > 1.35 * float(cap["r"]):
            continue
        if float(np.percentile(cyl.distance(verts), 90)) <= thr:
            found.add(id(region))
    return found


def find_obround_slots(m, regions, tol: float):
    """``(slots, exempt_plane_ids, cap_cylinder_ids)``.

    ``slots`` is a list of ``len, w, x, y, angle``. The id sets are empty
    when nothing was claimed.
    """
    span = _plate_z(m, regions, tol)
    if span is None:
        return [], set(), set()
    z0, z1 = span
    caps = _caps(m, regions, tol, z0, z1)
    if len(caps) < 2:
        return [], set(), set()
    pairs = []
    for i in range(len(caps)):
        for j in range(i + 1, len(caps)):
            got = _pair(m, caps[i], caps[j], regions, z0, z1, tol)
            if got is not None:
                pairs.append(got)
    pairs.sort(key=lambda p: (p["_residual"], p["y"], p["x"]))
    used_caps = set()
    used_walls = set()
    chosen = []
    for pair in pairs:
        cap_ids = {id(pair["_caps"][0]), id(pair["_caps"][1])}
        wall_ids = {id(r) for r in pair["_walls"]}
        if cap_ids & used_caps or wall_ids & used_walls:
            continue
        used_caps |= cap_ids
        used_walls |= wall_ids
        chosen.append(pair)
    if not chosen:
        return [], set(), set()
    exempt: set[int] = set()
    cyl_ids: set[int] = set()
    slots = []
    for pair in chosen:
        for cap in pair["_caps"]:
            cyl_ids |= cap["src"]
            exempt |= cap["planes"]
            exempt |= _nearby_planes(m, regions, cap, tol)
        exempt |= {id(r) for r in pair["_walls"]}
        slots.append({
            "len": pair["len"],
            "w": pair["w"],
            "x": pair["x"],
            "y": pair["y"],
            "angle": pair["angle"],
        })
    slots.sort(key=lambda s: (s["y"], s["x"], s["angle"], s["w"]))
    return slots, exempt, cyl_ids


def cadquery_slot_cut(name: str) -> str:
    """One obround cut. CadQuery's ``slot2D`` length is end to end.

    End to end is the centre distance plus the width, so ``name_len`` stays
    the distance between the end centres.
    """
    return (
        f"result = result.cut(cq.Workplane(\"XY\", origin=({name}_x, {name}_y, base_z0 - 1))"
        f".slot2D({name}_len + {name}_w, {name}_w, {name}_angle).extrude(base_t + 2))"
    )


def openscad_slot_cut(name: str) -> str:
    """Hull of the two end cylinders. ``name_len`` is the centre distance."""
    return (
        f"hull() {{\n"
        f"  translate([{name}_x - {name}_len / 2 * cos({name}_angle),\n"
        f"             {name}_y - {name}_len / 2 * sin({name}_angle),\n"
        f"             base_z0 - 1])\n"
        f"    cylinder(h=base_t + 2, d={name}_w);\n"
        f"  translate([{name}_x + {name}_len / 2 * cos({name}_angle),\n"
        f"             {name}_y + {name}_len / 2 * sin({name}_angle),\n"
        f"             base_z0 - 1])\n"
        f"    cylinder(h=base_t + 2, d={name}_w);\n"
        f"}}"
    )


FREECAD_SLOT_FN = '''
def _slot(length, width, angle, x, y, z0, height):
    """Obround cutter. ``length`` is the distance between the end centers."""
    ang = math.radians(angle)
    dx = length / 2.0 * math.cos(ang)
    dy = length / 2.0 * math.sin(ang)
    px = -width / 2.0 * math.sin(ang)
    py = width / 2.0 * math.cos(ang)
    pts = [
        App.Vector(x - dx + px, y - dy + py, z0),
        App.Vector(x + dx + px, y + dy + py, z0),
        App.Vector(x + dx - px, y + dy - py, z0),
        App.Vector(x - dx - px, y - dy - py, z0),
    ]
    area = 0.0
    for i, p in enumerate(pts):
        q = pts[(i + 1) % len(pts)]
        area += p.x * q.y - q.x * p.y
    if area < 0.0:
        pts.reverse()
    pts.append(pts[0])
    bar = Part.Face(Part.makePolygon(pts)).extrude(App.Vector(0, 0, height))
    c1 = Part.makeCylinder(width / 2.0, height, App.Vector(x - dx, y - dy, z0), App.Vector(0, 0, 1))
    c2 = Part.makeCylinder(width / 2.0, height, App.Vector(x + dx, y + dy, z0), App.Vector(0, 0, 1))
    return c1.fuse(c2).fuse(bar)
'''.strip("\n")


def freecad_slot_cut(name: str) -> str:
    return (
        f"result = result.cut(_slot({name}_len, {name}_w, {name}_angle, "
        f"{name}_x, {name}_y, base_z0 - 1, base_t + 2))"
    )
