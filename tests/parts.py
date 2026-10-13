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


PLANAR = {"cube": cube, "l_bracket": l_bracket, "stepped_block": stepped_block}
CURVED = {"plate_with_hole": plate_with_hole, "boss": boss}


def to_mesh(wp, tol=0.05, ang=0.2, max_edge=2.0):
    v, f = wp.val().tessellate(tol, ang)
    m = trimesh.Trimesh(np.array([p.toTuple() for p in v]), np.array(f), process=True)
    # make it scan-like: dense, roughly uniform triangles
    v2, f2 = trimesh.remesh.subdivide_to_size(m.vertices, m.faces, max_edge)
    return trimesh.Trimesh(v2, f2, process=True)
