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


def l_profile_y():
    """L in the XZ plane, extruded 18 mm along Y, with an Ø6 hole along X.

    The upright is the full height and the foot sticks out at the top, so the
    shape is not a Z-up box with a pad. Outline (x, z):
    (0, 0), (14, 0), (14, 16), (36, 16), (36, 28), (0, 28).
    The hole axis is y=9, z=22, through the foot.
    """
    outline = [(0, 0), (14, 0), (14, 16), (36, 16), (36, 28), (0, 28)]
    solid = (cq.Workplane(cq.Plane(origin=cq.Vector(0, 18, 0),
                                    xDir=cq.Vector(1, 0, 0),
                                    normal=cq.Vector(0, -1, 0)))
             .polyline(outline).close().extrude(18))
    hole = cq.Solid.makeCylinder(3.0, 50.0, cq.Vector(-5.0, 9.0, 22.0), cq.Vector(1, 0, 0))
    return solid.cut(cq.Workplane().add(hole))


def slanted_outline():
    """50×32 outline with one slanted corner, extruded 12 mm along Z from z=-6.

    Points (x, y): (0, 0), (50, 0), (50, 20), (34, 32), (0, 32).
    An Ø8 hole through Z at (14, 12).
    """
    outline = [(0, 0), (50, 0), (50, 20), (34, 32), (0, 32)]
    solid = cq.Workplane("XY").polyline(outline).close().extrude(12).translate((0, 0, -6))
    hole = cq.Solid.makeCylinder(4.0, 30.0, cq.Vector(14.0, 12.0, -20.0), cq.Vector(0, 0, 1))
    return solid.cut(cq.Workplane().add(hole))


def slotted_plate():
    """40×24×6 plate extruded along Z with a 16×8 rectangular through-slot."""
    outer = [(0, 0), (40, 0), (40, 24), (0, 24)]
    inner = [(8, 8), (24, 8), (24, 16), (8, 16)]
    return (cq.Workplane("XY").polyline(outer).close().extrude(6)
            .cut(cq.Workplane("XY").polyline(inner).close().extrude(6)))


def orthogonal_intersection():
    """Intersection of a Z extrusion and an X extrusion, with an Ø6 hole along Y.

    Z: z=0..12, outline (x, y) = (0, 0), (40, 0), (40, 16), (28, 28), (0, 28).
    X: x=-4, length 48, outline (y, z) = (0, -2), (28, -2), (28, 6), (18, 6),
    (18, 20), (0, 20).     The hole axis is y=10, z=4 and it passes through x. The Ø6 disk sits inside
    the profile, so Ø8 still does.
    """
    z_outline = [(0, 0), (40, 0), (40, 16), (28, 28), (0, 28)]
    z_solid = cq.Workplane("XY").polyline(z_outline).close().extrude(12)
    yz = [(0, -2), (28, -2), (28, 6), (18, 6), (18, 20), (0, 20)]
    x_solid = (cq.Workplane(cq.Plane(origin=cq.Vector(-4, 0, 0),
                                      xDir=cq.Vector(0, 1, 0),
                                      normal=cq.Vector(1, 0, 0)))
               .polyline(yz).close().extrude(48))
    hole = cq.Solid.makeCylinder(3.0, 60.0, cq.Vector(-10.0, 10.0, 4.0), cq.Vector(1, 0, 0))
    return z_solid.intersect(x_solid).cut(cq.Workplane().add(hole))


def cut_block():
    """56×36×18 block with a slot cut along X and an Ø6 through-hole along Y.

    The block is extruded along Z from z=0. The slot starts at x=12, runs 28 mm,
    and its profile in YZ is (10, 8), (24, 8), (24, 22), (10, 22), so it opens
    on the top face. The hole axis is x=6, z=4, outside the slot.
    """
    base = (cq.Workplane("XY").polyline([(0, 0), (56, 0), (56, 36), (0, 36)])
            .close().extrude(18))
    slot = (cq.Workplane(cq.Plane(origin=cq.Vector(12, 0, 0),
                                   xDir=cq.Vector(0, 1, 0),
                                   normal=cq.Vector(1, 0, 0)))
            .polyline([(10, 8), (24, 8), (24, 22), (10, 22)]).close().extrude(28))
    hole = cq.Solid.makeCylinder(3.0, 50.0, cq.Vector(6.0, -7.0, 4.0), cq.Vector(0, 1, 0))
    return base.cut(slot).cut(cq.Workplane().add(hole))


def counterbored_plate():
    """40×28×14 plate with an Ø6 through-hole and an Ø12 counterbore, 4 mm deep."""
    return cq.Workplane("XY").box(40, 28, 14).faces(">Z").workplane().cboreHole(6, 12, 4)


def partial_edge_chamfer():
    """40×30×16 block with a 2 mm 45° chamfer on only the middle of one top edge.

    The chamfer runs along X from x=-10 to x=10 on the +Y/+Z corner. The two
    ends of that edge stay square, so every full-length silhouette is still
    the box and the chamfer is a separate cut.
    """
    block = cq.Workplane("XY").box(40, 30, 16)
    cutter = (cq.Workplane(cq.Plane(origin=cq.Vector(-10, 0, 0),
                                     xDir=cq.Vector(0, 1, 0),
                                     normal=cq.Vector(1, 0, 0)))
              .polyline([(15, 8), (13, 8), (15, 6)]).close().extrude(20))
    return block.cut(cutter)


def breakout_slant():
    """Plate with one slanted wall and a slot that opens into a through hole.

    Outline (x, y): (-25, -18), (25, -18), (25, 10), (10, 16), (-25, 16),
    extruded 10 mm. Ø10 hole at the origin. The slot is x=-4..4 from the -Y
    edge into that hole, so the hole is part of the outer loop.
    """
    outline = [(-25, -18), (25, -18), (25, 10), (10, 16), (-25, 16)]
    plate = cq.Workplane("XY").polyline(outline).close().extrude(10)
    hole = cq.Solid.makeCylinder(5.0, 20.0, cq.Vector(0, 0, -5), cq.Vector(0, 0, 1))
    slot = (cq.Workplane("XY").polyline([(-4, -18), (4, -18), (4, 0), (-4, 0)])
            .close().extrude(10))
    return plate.cut(cq.Workplane().add(hole)).cut(slot)


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
