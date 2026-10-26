"""A counterbore that breaks out of a side face is only a partial ring of planes.

The STEP segmentation leaves that ring as planes. Parametric export fits the
arc and, when it is coaxial with a through hole and the gap out to the part
end is empty, writes it as a counterbored through-hole. A tangent fillet and a
flat polygon whose corners lie on a circle stay planes. The self-check gate
is the same one every other script uses.
"""
import math

import numpy as np
import pytest
import trimesh

pytest.importorskip("cadquery")
import cadquery as cq  # noqa: E402

from mesh2cad.cadquery_code import check as cq_check  # noqa: E402
from mesh2cad.cadquery_code import run_script, to_cadquery  # noqa: E402
from mesh2cad.convert import convert_mesh  # noqa: E402
from mesh2cad.fit import Cylinder, Plane  # noqa: E402
from mesh2cad.openscad_code import to_openscad  # noqa: E402
from mesh2cad.segment import Region, partial_arc_cylinders, segment  # noqa: E402

# Breakout counterbore: centre x=2.2, R=5, block starts at x=0.
# The wall that remains inside the block is about 232 degrees of arc.
_CBORE_R = 5.0
_HOLE_R = 1.5
_CBORE_DEPTH = 6.0


def _breakout_part(tab=False):
    """Block, a through hole, and a counterbore that opens through the x=0 face.

    With ``tab``, a pad on the far end sets the -Y bound past the counterbore
    mouth. The mouth itself is not on that pad, so the gap is empty.
    """
    block = cq.Workplane("XY").box(40, 24, 14, centered=(False, False, False))
    solid = block.val()
    if tab:
        pad = cq.Workplane("XY").box(10, 8, 14, centered=(False, False, False)).translate((30, -8, 0))
        solid = solid.fuse(pad.val())
    hole = cq.Solid.makeCylinder(_HOLE_R, 40.0, cq.Vector(2.2, -12.0, 7.0), cq.Vector(0, 1, 0))
    cbore = cq.Solid.makeCylinder(_CBORE_R, _CBORE_DEPTH, cq.Vector(2.2, 0.0, 7.0), cq.Vector(0, 1, 0))
    return cq.Workplane(obj=solid).cut(cq.Workplane().add(hole)).cut(cq.Workplane().add(cbore))


def _mesh_of(wp, tol=0.08, ang=0.25):
    v, f = wp.val().tessellate(tol, ang)
    return trimesh.Trimesh(np.array([p.toTuple() for p in v]), np.array(f), process=True)


def _params(code):
    import re
    return {k: float(v) for k, v in re.findall(r"^(\w+)\s*=\s*([-\d.]+)", code, re.M)}


def _cap_area(radius, distance):
    """Area of the disk of ``radius`` on the far side of a plane ``distance`` from the centre."""
    if distance >= radius:
        return 0.0
    if distance <= -radius:
        return math.pi * radius * radius
    height = radius - distance
    return (radius * radius * math.acos((radius - height) / radius)
            - (radius - height) * math.sqrt(max(0.0, 2.0 * radius * height - height * height)))


def _inside_area(radius, distance):
    """Area of the disk that lies on the material side of the breakout plane."""
    return math.pi * radius * radius - _cap_area(radius, distance)


def test_partial_arc_recovers_breakout_counterbore():
    mesh = _mesh_of(_breakout_part())
    tol = 0.05
    regions = segment(mesh, tol)
    fitted = [r.surface.radius for r in regions
              if isinstance(r.surface, Cylinder) and abs(r.surface.radius - _CBORE_R) < 0.4]
    assert fitted == []
    found = partial_arc_cylinders(mesh, regions, tol)
    arcs = [r for r in found if abs(r.surface.radius - _CBORE_R) < 0.3]
    assert len(arcs) == 1
    arc = arcs[0].surface
    assert abs(abs(arc.axis[1]) - 1.0) < 1e-6
    assert arc.point[0] == pytest.approx(2.2, abs=0.15)
    assert arc.point[2] == pytest.approx(7.0, abs=0.15)
    pts = mesh.vertices[np.unique(mesh.faces[arcs[0].faces])]
    from mesh2cad.segment import _arc_span, _cylinder_angles
    span, second, _mid = _arc_span(_cylinder_angles(arc.axis, arc.point, pts)[0])
    assert math.degrees(span) == pytest.approx(232.0, abs=20.0)
    assert math.degrees(second) < 40.0


def test_partial_arc_rejects_fillet_and_flat_polygon():
    """A 90 degree fillet and a coarse flat ring are not cylinders.

    The fillet is exploded back into one plane per triangle, which is the
    state a partial round has before this fitter runs. The flat ring is eight
    chords of a circle: the corners lie on the circle and the face centres do not.
    """
    bracket = (cq.Workplane("XY").polyline([(0, 0), (40, 0), (40, 8), (8, 8), (8, 30), (0, 30)])
               .close().extrude(12).edges("|Z").fillet(3))
    mesh = _mesh_of(bracket, tol=0.1, ang=0.3)
    tol = 0.05
    regions = segment(mesh, tol)
    exploded = []
    for r in regions:
        if isinstance(r.surface, Cylinder):
            for f in r.faces:
                f = int(f)
                pl = Plane(mesh.face_normals[f], float(mesh.face_normals[f] @ mesh.triangles_center[f]))
                exploded.append(Region(np.array([f]), pl, float(mesh.area_faces[f])))
        else:
            exploded.append(r)
    found = partial_arc_cylinders(mesh, exploded, tol)
    assert not any(isinstance(r.surface, Cylinder) and abs(r.surface.radius - 3.0) < 0.5 for r in found)

    flat = _flat_ring(radius=15.0, n=8, span=math.radians(200.0), height=4.0)
    found = partial_arc_cylinders(flat, _one_plane_per_face(flat), 0.05)
    assert found == []

    fillet = _tangent_arc(radius=5.0, a0=math.radians(180.0), a1=math.radians(300.0), n=10, height=8.0)
    found = partial_arc_cylinders(fillet, _one_plane_per_face(fillet), 0.05)
    assert found == []


def _one_plane_per_face(mesh):
    regs = []
    for f in range(len(mesh.faces)):
        n = mesh.face_normals[f]
        regs.append(Region(np.array([f]), Plane(n, float(n @ mesh.triangles_center[f])),
                           float(mesh.area_faces[f])))
    return regs


def _flat_ring(radius, n, span, height):
    """``n`` flat quads whose corners lie on a circle. Face centres sit inside it."""
    angles = np.linspace(-0.5 * span, 0.5 * span, n + 1)
    verts = []
    faces = []
    for i, a in enumerate(angles):
        verts.append((radius * math.cos(a), radius * math.sin(a), 0.0))
        verts.append((radius * math.cos(a), radius * math.sin(a), height))
    for i in range(n):
        a, b = 2 * i, 2 * (i + 1)
        c, d = a + 1, b + 1
        faces.append((a, b, d))
        faces.append((a, d, c))
    wall = trimesh.Trimesh(vertices=np.array(verts), faces=np.array(faces), process=False)
    # A large closed box keeps each chord under the small-face cut.
    box = trimesh.creation.box(extents=(200, 200, 200))
    box.apply_translation((400, 0, 0))
    return trimesh.util.concatenate([wall, box])


def _tangent_arc(radius, a0, a1, n, height):
    """A partial cylinder with a tangent wall at each end. The arc is a fillet.

    The walls share the arc's end vertices, so the fitter can see that both
    planes are tangent and leave the arc as planes.
    """
    angles = np.linspace(a0, a1, n + 1)
    verts = []
    for a in angles:
        verts.append([radius * math.cos(a), radius * math.sin(a), 0.0])
        verts.append([radius * math.cos(a), radius * math.sin(a), height])
    faces = []
    for i in range(n):
        a, b = 2 * i, 2 * (i + 1)
        c, d = a + 1, b + 1
        faces.append((a, b, d))
        faces.append((a, d, c))

    def add_wall(end, tangent):
        far = np.array(verts[end]) + np.array(tangent) * (4.0 * radius)
        i0 = len(verts)
        verts.append(far.tolist())
        verts.append((far + (0.0, 0.0, height)).tolist())
        # end is the z=0 vertex; end+1 is z=height. Both already exist.
        faces.append((end, i0, i0 + 1))
        faces.append((end, i0 + 1, end + 1))

    t0 = np.array((-math.sin(a0), math.cos(a0), 0.0))
    t1 = np.array((-math.sin(a1), math.cos(a1), 0.0))
    add_wall(0, -t0)
    add_wall(2 * n, t1)
    mesh = trimesh.Trimesh(vertices=np.array(verts, float), faces=np.array(faces), process=False)
    # Arc normals point at the axis, the way a concave fillet does.
    centres = mesh.triangles_center
    radial = centres.copy()
    radial[:, 2] = 0.0
    radial /= np.maximum(np.linalg.norm(radial, axis=1, keepdims=True), 1e-9)
    outward = np.sum(mesh.face_normals * radial, axis=1) > 0
    faces_a = np.array(mesh.faces)
    faces_a[outward] = faces_a[outward][:, ::-1]
    mesh = trimesh.Trimesh(vertices=mesh.vertices, faces=faces_a, process=False)
    box = trimesh.creation.box(extents=(200, 200, 200))
    box.apply_translation((400, 0, 0))
    return trimesh.util.concatenate([mesh, box])


def test_cadquery_breakout_counterbore():
    ref = _mesh_of(_breakout_part())
    res = convert_mesh(ref, tol=0.05, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["hole_1_d"] == pytest.approx(2.0 * _HOLE_R, abs=0.15)
    assert params["hole_1_cbore_d"] == pytest.approx(2.0 * _CBORE_R, abs=0.25)
    assert params["hole_1_cbore_depth"] == pytest.approx(_CBORE_DEPTH, abs=0.3)
    assert "counterbore depth" in code
    against = cq_check(code, ref)
    assert against["max_deviation"] < 0.15
    assert against["volume_error_pct"] < 0.5


def test_edit_breakout_counterbore_diameter():
    """Widening the counterbore removes the part of the annulus that lies in the block."""
    ref = _mesh_of(_breakout_part())
    res = convert_mesh(ref, tol=0.05, with_report=False)
    code, _ = to_cadquery(res.mesh, res.regions, res.tol)
    params = _params(code)
    v0 = abs(run_script(code).volume)
    import re
    wider = re.sub(r"^hole_1_cbore_d\s*=\s*[-\d.]+", "hole_1_cbore_d = 12.0", code, count=1, flags=re.M)
    # centre is 2.2 mm from the x=0 face. The cutter is a full cylinder; only
    # the portion with x>=0 meets the block. The through hole is already gone.
    depth = params["hole_1_cbore_depth"]
    old_r = params["hole_1_cbore_d"] / 2.0
    new_r = 6.0
    delta = (_inside_area(new_r, 2.2) - _inside_area(old_r, 2.2)) * depth
    assert v0 - abs(run_script(wider).volume) == pytest.approx(delta, rel=0.03)


def test_breakout_counterbore_reaches_across_an_empty_end():
    """A tab past the mouth does not shorten the counterbore or move its shoulder."""
    ref = _mesh_of(_breakout_part(tab=True))
    res = convert_mesh(ref, tol=0.05, with_report=False)
    code, q = to_cadquery(res.mesh, res.regions, res.tol)
    assert q["max_deviation"] <= max(4 * res.tol, 0.1)
    assert q["volume_error_pct"] <= 1.0
    params = _params(code)
    assert params["hole_1_cbore_d"] == pytest.approx(2.0 * _CBORE_R, abs=0.25)
    # The -Y bound is the tab (y=-8), the shoulder is at y=6, and the gap is empty.
    assert params["hole_1_cbore_depth"] == pytest.approx(14.0, abs=0.4)


def test_openscad_breakout_counterbore():
    ref = _mesh_of(_breakout_part())
    res = convert_mesh(ref, tol=0.05, with_report=False)
    code, q = to_openscad(res.mesh, res.regions, res.tol, verify=False)
    cq, _ = to_cadquery(res.mesh, res.regions, res.tol, verify=False)
    assert _params(code) == _params(cq)
    assert "hole_1_cbore_d" in code
    assert "hole_1_cbore_depth" in code
    if q is not None:
        assert q["max_deviation"] <= max(4 * res.tol, 0.1)
        assert q["volume_error_pct"] <= 1.0
