"""Synthetic ground-truth parts (cadquery) -> trimesh."""
import cadquery as cq
import numpy as np
import trimesh


def cube():
    return cq.Workplane().box(20, 20, 20)


def l_bracket():
    return (cq.Workplane("XY").polyline([(0, 0), (40, 0), (40, 5), (5, 5), (5, 30), (0, 30)])
            .close().extrude(20))


def stepped_block():
    return cq.Workplane().box(40, 30, 10).faces(">Z").workplane().rect(20, 15).extrude(8)


def plate_with_hole():
    return cq.Workplane().box(40, 30, 6).faces(">Z").workplane().hole(10)


def boss():
    return cq.Workplane().box(40, 40, 5).faces(">Z").workplane().circle(8).extrude(15)


def ribbed_plate():
    return (cq.Workplane().box(50, 30, 6).faces(">Z").workplane()
            .pushPoints([(-15, 0), (15, 0)]).hole(8)
            .faces(">Z").workplane().center(0, 0).rect(14, 30).extrude(10))


def spherical_dome():
    """40×40×10 block with a hemisphere (r=8) centered on the top face."""
    block = cq.Workplane().box(40, 40, 10)
    dome = cq.Workplane("XY").workplane(offset=5).sphere(8)
    return block.union(dome)


def spherical_pocket():
    """40×40×10 block with a spherical pocket (r=6) opening through the top face."""
    block = cq.Workplane().box(40, 40, 10)
    cutter = cq.Workplane("XY").workplane(offset=5).sphere(6)
    return block.cut(cutter)


def spherical_zone():
    """Sphere (r=10) clipped to a slab: one spherical zone between two planes."""
    return cq.Workplane().sphere(10).intersect(cq.Workplane().box(30, 30, 12))


def countersunk_plate():
    """40×30×6 plate with a Ø6 through hole and a 90° Ø12 countersink."""
    return cq.Workplane().box(40, 30, 6).faces(">Z").workplane().cskHole(6, 12, 90)


def tapered_boss():
    """40×40×5 plate with a Ø16 boss drafted 10° over 10 mm."""
    return cq.Workplane().box(40, 40, 5).faces(">Z").workplane().circle(8).extrude(10, taper=10)


def chamfered_boss():
    """Ø16 cylinder boss with a 2 mm chamfer on its top edge."""
    return (cq.Workplane().box(40, 40, 5).faces(">Z").workplane().circle(8).extrude(12)
            .faces(">Z").edges().chamfer(2))


def filleted_plate():
    """40×30×6 plate with the four vertical corners rounded at r=5."""
    return cq.Workplane().box(40, 30, 6).edges("|Z").fillet(5)


def filleted_box():
    """30×24×16 block with the four vertical edges rounded at r=2."""
    return cq.Workplane().box(30, 24, 16).edges("|Z").fillet(2)


def side_hole_boss():
    """60×40×36 block with three side features.

    Ø8 through hole along X, axis at y=0, z=-6.
    Ø6 blind hole, 12 deep, opening on the +Y face, axis at x=-12, z=8.
    Ø10 boss, 8 tall, on the +X face, axis at y=10, z=8.
    """
    block = cq.Workplane().box(60, 40, 36).val()
    # cutters run 1 mm past the opening so the boolean does not leave a coincident face
    through = cq.Solid.makeCylinder(4.0, 90.0, cq.Vector(-45.0, 0.0, -6.0), cq.Vector(1, 0, 0))
    blind = cq.Solid.makeCylinder(3.0, 13.0, cq.Vector(-12.0, 8.0, 8.0), cq.Vector(0, 1, 0))
    boss = cq.Solid.makeCylinder(5.0, 9.0, cq.Vector(29.0, 10.0, 8.0), cq.Vector(1, 0, 0))
    return cq.Workplane(obj=block.cut(through).cut(blind).fuse(boss))


def angled_hole():
    """40×40×20 block with a hole whose axis is in the XY plane but not on X or Y."""
    block = cq.Workplane().box(40, 40, 20).val()
    cut = cq.Solid.makeCylinder(4.0, 80.0, cq.Vector(-30.0, -14.0, 0.0), cq.Vector(1.0, 0.4, 0.0))
    return cq.Workplane(obj=block.cut(cut))


def filleted_bracket():
    """L-bracket extruded 12 mm; every vertical edge, including the inner corner, at r=3."""
    return (cq.Workplane("XY").polyline([(0, 0), (40, 0), (40, 8), (8, 8), (8, 30), (0, 30)])
            .close().extrude(12).edges("|Z").fillet(3))


PLANAR = {"cube": cube, "l_bracket": l_bracket, "stepped_block": stepped_block}
FILLETS = {"filleted_plate": filleted_plate, "filleted_box": filleted_box,
           "filleted_bracket": filleted_bracket}
FILLET_RADIUS = {"filleted_plate": 5.0, "filleted_box": 2.0, "filleted_bracket": 3.0}
CURVED = {"plate_with_hole": plate_with_hole, "boss": boss, "ribbed_plate": ribbed_plate}
CONES = {"countersunk_plate": countersunk_plate, "tapered_boss": tapered_boss,
         "chamfered_boss": chamfered_boss}
SPHERES = {"spherical_dome": spherical_dome, "spherical_pocket": spherical_pocket}


def to_mesh(wp, tol=0.05, ang=0.2, max_edge=2.0):
    v, f = wp.val().tessellate(tol, ang)
    m = trimesh.Trimesh(np.array([p.toTuple() for p in v]), np.array(f), process=True)
    # make it scan-like: dense, roughly uniform triangles
    v2, f2 = trimesh.remesh.subdivide_to_size(m.vertices, m.faces, max_edge)
    out = trimesh.Trimesh(v2, f2, process=True)
    # sphere tessellation can leave a zero-area triangle at the pole
    if len(out.faces) and (out.area_faces < 1e-12).any():
        out.update_faces(out.nondegenerate_faces())
        out.remove_unreferenced_vertices()
    return out
