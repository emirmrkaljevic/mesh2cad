"""Emit an editable, parametric CadQuery script from a fitted model.

Supported (Z-up, axis-aligned) parts: a box/plate base, rectangular pads on its
top face, vertical holes (through, blind, countersunk), vertical bosses
(straight, drafted, chamfered), and one constant radius on all four vertical
corners of the base (``.edges("|Z").fillet(fillet_r)``). Every dimension
becomes a named parameter at the top of the script (``hole_1_d = 6.0``) so
the part can be edited and re-run. Anything outside that vocabulary raises
:class:`Unsupported`; the generated script is checked against the input mesh
before it is returned.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import trimesh

from .fit import Cone, Cylinder, Plane
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
    x: float
    y: float
    outward: bool
    segs: list = field(default_factory=list)


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

    The script can only say ``edges("|Z").fillet(r)`` for that one case. Any
    other fillet (one edge, an L, a horizontal round) is left to STEP.
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


def _stacks(m, regions, tol, zs=()):
    stacks: list[Stack] = []
    gtol = max(5 * tol, 0.05)
    for r in regions:
        s = r.surface
        if not isinstance(s, (Cylinder, Cone)):
            continue
        if abs(abs(s.axis @ Z) - 1) > 1e-3:
            raise Unsupported("curved surface with a non-vertical axis")
        c = s.point if isinstance(s, Cylinder) else s.apex
        v = _verts(m, r)
        rho = np.linalg.norm(v[:, :2] - c[:2], axis=1)
        eps = max(4 * tol, 0.05)
        za, zb = _snap(v[:, 2].min(), zs, eps), _snap(v[:, 2].max(), zs, eps)
        if isinstance(s, Cylinder):
            ra = rb = s.radius
            kind = "cylinder"
        else:
            b, a = np.polyfit(v[:, 2], rho, 1)
            ra, rb, kind = a + b * za, a + b * zb, "cone"
        fc = m.triangles_center[r.faces]
        radial = fc[:, :2] - c[:2]
        radial /= np.maximum(np.linalg.norm(radial, axis=1, keepdims=True), 1e-12)
        out = float(np.mean(np.sum(m.face_normals[r.faces][:, :2] * radial, axis=1))) > 0
        for st in stacks:
            if np.hypot(st.x - c[0], st.y - c[1]) < gtol and st.outward == out:
                st.segs.append(Segment(kind, za, zb, ra, rb))
                break
        else:
            stacks.append(Stack(c[0], c[1], out, [Segment(kind, za, zb, ra, rb)]))
    for st in stacks:
        st.segs.sort(key=lambda g: g.za)
        for g, h in zip(st.segs, st.segs[1:]):  # make junctions meet
            if abs(g.zb - h.za) <= max(4 * tol, 0.05):
                zj = (g.zb + h.za) / 2
                g.zb = h.za = zj
    for st in stacks:
        for g in st.segs:
            if g.kind == "cone":
                k = next(r.surface for r in regions if isinstance(r.surface, Cone)
                         and np.hypot(*(r.surface.apex[:2] - [st.x, st.y])) < gtol)
                g.ra = abs(k.apex[2] - g.za) * math.tan(k.half_angle)
                g.rb = abs(k.apex[2] - g.zb) * math.tan(k.half_angle)
    return stacks


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
    zs = [c for _, i, _, c in pos if i == 2]
    # corner fillets are partial cylinders; leaving them in the stacks makes
    # them look like bosses that do not start on the top face
    fillet_regs = _vertical_fillets(m, regions)
    fillet_ids = {id(r) for r in fillet_regs}
    stacks = _stacks(m, [r for r in regions if id(r) not in fillet_ids], tol, zs + [zmin])

    def on_stack(r):
        v = _verts(m, r)
        for st in stacks:
            rmax = max(max(g.ra, g.rb) for g in st.segs)
            if np.percentile(np.hypot(v[:, 0] - st.x, v[:, 1] - st.y), 95) <= rmax + max(3 * tol, 0.05):
                return True
        return False

    up = []
    for r in planes:
        nz = r.surface.normal @ Z
        z = _plane_pos(r.surface)[2] if abs(nz) > 0.99 else 0.0
        if nz < -0.999 and z > zmin + max(3 * tol, 0.05):
            raise Unsupported("downward-facing face above the bottom (overhang)")
        if nz > 0.999:
            up.append((r, z))
    if not up:
        raise Unsupported("no top face")
    base_r, ztop = max(up, key=lambda t: t[0].area)
    xs = [c for _, i, _, c in pos if i == 0]
    ys = [c for _, i, _, c in pos if i == 1]
    if not xs or not ys:
        raise Unsupported("no side faces")
    lo, hi = np.array([min(xs), min(ys)]), np.array([max(xs), max(ys)])
    feats = {"base": dict(cx=(lo[0] + hi[0]) / 2, cy=(lo[1] + hi[1]) / 2, z0=zmin,
                          x=hi[0] - lo[0], y=hi[1] - lo[1], t=ztop - zmin),
             "pads": [], "holes": [], "bosses": []}
    eps = max(3 * tol, 0.05)
    for r, z in up:
        if r is base_r or on_stack(r) or z < ztop + eps:
            if r is not base_r and not on_stack(r):
                raise Unsupported("pocket or step below the top face")
            continue
        pv = _verts(m, r)
        e2 = max(4 * tol, 0.1)
        plo = np.array([_snap(pv[:, 0].min(), xs, e2), _snap(pv[:, 1].min(), ys, e2)])
        phi = np.array([_snap(pv[:, 0].max(), xs, e2), _snap(pv[:, 1].max(), ys, e2)])
        feats["pads"].append(dict(cx=(plo[0] + phi[0]) / 2, cy=(plo[1] + phi[1]) / 2,
                                  w=phi[0] - plo[0], l=phi[1] - plo[1], h=z - ztop))
    for st in stacks:
        g = st.segs
        if not st.outward:
            h = dict(x=st.x, y=st.y)
            cyl = [s for s in g if s.kind == "cylinder"]
            if len(cyl) != 1 or cyl[0] is not g[0]:
                raise Unsupported("hole profile is not a plain or countersunk hole")
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
            feats["holes"].append(h)
        else:
            b = dict(x=st.x, y=st.y)
            if abs(g[0].za - ztop) > eps:
                raise Unsupported("boss does not start on the top face")
            b["d"] = 2 * g[0].ra
            b["h"] = g[-1].zb - ztop
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


def render(feats: dict, tol: float = 0.01, source: str | None = None) -> str:
    """Turn a feature dict into a CadQuery script with named parameters."""
    P: list[tuple[str, float, str]] = []

    def p(name, val, note=""):
        t = 0.5 if ("angle" in name or "deg" in name) else tol  # angles in degrees
        P.append((name, _nice(val, t), note))
        return name

    b = feats["base"]
    p("base_x", b["x"], "base length (X)")
    p("base_y", b["y"], "base width (Y)")
    p("base_t", b["t"], "base thickness (Z)")
    p("base_cx", b["cx"]); p("base_cy", b["cy"]); p("base_z0", b["z0"], "bottom face height")
    body = ["result = (cq.Workplane(\"XY\", origin=(base_cx, base_cy, base_z0))",
            "          .rect(base_x, base_y).extrude(base_t))",
            "top = base_z0 + base_t"]
    if "fillet_r" in feats:
        p("fillet_r", feats["fillet_r"], "vertical corner radius")
        body.append("result = result.edges(\"|Z\").fillet(fillet_r)")
    for i, d in enumerate(feats["pads"], 1):
        n = f"pad_{i}"
        for k, note in (("w", "X size"), ("l", "Y size"), ("h", "height above top face")):
            p(f"{n}_{k}", d[k], note)
        p(f"{n}_x", d["cx"]); p(f"{n}_y", d["cy"])
        body.append(f"result = result.union(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top))"
                    f".rect({n}_w, {n}_l).extrude({n}_h))")
    for i, d in enumerate(feats["bosses"], 1):
        n = f"boss_{i}"
        p(f"{n}_d", d["d"], "diameter at the top face"); p(f"{n}_h", d["h"], "height")
        p(f"{n}_x", d["x"]); p(f"{n}_y", d["y"])
        if "draft_deg" in d:
            p(f"{n}_draft_deg", d["draft_deg"], "taper angle")
            body.append(f"result = result.union(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top))"
                        f".circle({n}_d / 2).extrude({n}_h, taper={n}_draft_deg))")
        elif "chamfer" in d:
            if abs(d["chamfer"] - d["chamfer_w"]) <= max(2 * tol, 0.02):
                p(f"{n}_chamfer", d["chamfer"], "45° top-edge chamfer")
                cw = ch = f"{n}_chamfer"
            else:
                p(f"{n}_chamfer_h", d["chamfer"]); p(f"{n}_chamfer_w", d["chamfer_w"])
                ch, cw = f"{n}_chamfer_h", f"{n}_chamfer_w"
            body += [f"result = result.union(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top))"
                     f".circle({n}_d / 2).extrude({n}_h - {ch}))",
                     f"result = result.union(cq.Workplane().add(cq.Solid.makeCone({n}_d / 2, "
                     f"{n}_d / 2 - {cw}, {ch}, cq.Vector({n}_x, {n}_y, top + {n}_h - {ch}))))"]
        else:
            body.append(f"result = result.union(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top))"
                        f".circle({n}_d / 2).extrude({n}_h))")
    for i, d in enumerate(feats["holes"], 1):
        n = f"hole_{i}"
        p(f"{n}_d", d["d"], "diameter"); p(f"{n}_x", d["x"]); p(f"{n}_y", d["y"])
        if "depth" in d:
            p(f"{n}_depth", d["depth"], "blind depth from top face")
            body.append(f"result = result.cut(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, top - {n}_depth))"
                        f".circle({n}_d / 2).extrude({n}_depth + 1))")
        else:
            body.append(f"result = result.cut(cq.Workplane(\"XY\", origin=({n}_x, {n}_y, base_z0 - 1))"
                        f".circle({n}_d / 2).extrude(base_t + 2))")
        if "csk_d" in d:
            p(f"{n}_csk_d", d["csk_d"], "countersink diameter")
            p(f"{n}_csk_angle", d["csk_angle"], "countersink angle")
            body += [f"_h = ({n}_csk_d - {n}_d) / 2 / math.tan(math.radians({n}_csk_angle / 2))",
                     f"result = result.cut(cq.Workplane().add(cq.Solid.makeCone({n}_d / 2, "
                     f"{n}_csk_d / 2 + 1 * math.tan(math.radians({n}_csk_angle / 2)), _h + 1, "
                     f"cq.Vector({n}_x, {n}_y, top - _h))))"]
    w = max(len(n) for n, _, _ in P)
    head = ['"""Parametric CadQuery model generated by mesh2cad.',
            "",
            f"Source: {source}" if source else "Edit the parameters below and re-run.",
            "Units are the mesh units (usually mm). Run: python this_file.py",
            '"""',
            "import math", "", "import cadquery as cq", "", "# --- parameters ---"]
    params = [f"{n.ljust(w)} = {v!r}" + (f"  # {note}" if note else "") for n, v, note in P]
    tail = ["", "if __name__ == \"__main__\":",
            "    from pathlib import Path",
            "    cq.exporters.export(result, str(Path(__file__).with_suffix(\".step\")))",
            "    print(\"wrote\", Path(__file__).with_suffix(\".step\"))", ""]
    return "\n".join(head + params + ["", "# --- model ---"] + body + tail)


def run_script(code: str):
    """Execute a generated script and return its ``result`` as a trimesh."""
    ns: dict = {"__name__": "mesh2cad_cq"}
    exec(compile(code, "<mesh2cad-cadquery>", "exec"), ns)
    v, f = ns["result"].val().tessellate(0.01, 0.1)
    return trimesh.Trimesh(np.array([q.toTuple() for q in v]), np.array(f), process=True)


def check(code: str, m: trimesh.Trimesh, n: int = 4000) -> dict:
    """Run the script and measure how well it reproduces mesh ``m``."""
    out = run_script(code)
    a, _ = trimesh.sample.sample_surface(out, n, seed=0)
    b, _ = trimesh.sample.sample_surface(m, n, seed=0)
    d1 = np.abs(trimesh.proximity.signed_distance(m, a))
    d2 = np.abs(trimesh.proximity.signed_distance(out, b))
    vm = abs(m.volume)
    return {"max_deviation": float(max(d1.max(), d2.max())),
            "rms": float(np.sqrt(np.mean(np.r_[d1, d2] ** 2))),
            "volume_error_pct": float(100 * abs(abs(out.volume) - vm) / vm)}


def to_cadquery(m: trimesh.Trimesh, regions, tol: float, source: str | None = None,
                verify: bool = True) -> tuple[str, dict | None]:
    """Generate the script; when ``verify`` and cadquery is installed, run it and
    raise :class:`Unsupported` if it does not reproduce the mesh."""
    code = render(extract(m, regions, tol), tol=tol, source=source)
    if not verify:
        return code, None
    try:
        import cadquery  # noqa: F401
    except ImportError:
        return code, None
    try:
        q = check(code, m)
    except Exception as e:  # noqa: BLE001
        raise Unsupported(f"generated script failed to run: {e}") from e
    lim = max(4 * tol, 0.1)
    if q["max_deviation"] > lim or q["volume_error_pct"] > 1.0:
        raise Unsupported(f"generated script deviates by {q['max_deviation']:.3f} "
                          f"(limit {lim:.3f}), volume error {q['volume_error_pct']:.2f}%")
    return code, q
