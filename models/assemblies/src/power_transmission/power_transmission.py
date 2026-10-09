"""Power-transmission assembly from the textbook section drawing.

Datum: shaft axis along +X, up +Z, base sitting on Z=0.
The drawing states the 70 mm center height. Bearing 6003 and the M2 Z34 spur
gear fix the mating sizes; every other length is an assumption.
"""

from __future__ import annotations

import math
from pathlib import Path

from cadgen import build123d as bd
from cadgen import srgb, step, stl

CENTER_Z = 70.0
HOUSING_LEN = 52.0
BARREL_R = 30.0
BASE_H = 12.0
BCD = 48.0
BOLT_ANGLES = (45.0, 135.0, 225.0, 315.0)

GEAR_MODULE = 2.0
GEAR_TEETH = 34
GEAR_PA = math.radians(20.0)
GEAR_FACE_X = 68.0
GEAR_FACE_W = 18.0
GEAR_NOSE = 5.0

CAST = srgb("#6E7378")
COVER_COLOR = srgb("#8B919A")
SHAFT_COLOR = srgb("#C5CAD1")
GEAR_COLOR = srgb("#3C4A57")
COUPLING_COLOR = srgb("#4E5964")
BEARING_COLOR = srgb("#E2E6EB")
BOLT_COLOR = srgb("#2A3036")
SEAL_COLOR = srgb("#1C1C1C")
KEY_COLOR = srgb("#8A7349")


def cyl_x(radius: float, x0: float, x1: float, y: float = 0.0, z: float = CENTER_Z) -> bd.Solid:
    """Closed cylinder along X, from x0 to x1, axis at (y, z)."""
    solid = bd.Rot(0, 90, 0) * bd.Cylinder(radius, x1 - x0)
    return bd.Pos((x0 + x1) / 2.0, y, z) * solid


def ring_x(outer_r: float, inner_r: float, x0: float, x1: float, y: float = 0.0, z: float = CENTER_Z) -> bd.Solid:
    return cyl_x(outer_r, x0, x1, y, z) - cyl_x(inner_r, x0 - 0.2, x1 + 0.2, y, z)


def bolt_center(angle_deg: float) -> tuple[float, float]:
    angle = math.radians(angle_deg)
    return (BCD / 2.0) * math.cos(angle), CENTER_Z + (BCD / 2.0) * math.sin(angle)


def _rot2(x: float, y: float, angle: float) -> tuple[float, float]:
    c, s = math.cos(angle), math.sin(angle)
    return (c * x - s * y, s * x + c * y)


def _involute(base_r: float, t: float) -> tuple[float, float]:
    return (
        base_r * (math.cos(t) + t * math.sin(t)),
        base_r * (math.sin(t) - t * math.cos(t)),
    )


def spur_gear_points(module: float, teeth: int, pressure: float, samples: int = 10) -> list[tuple[float, float]]:
    """Full involute outline, module in mm, centered on the origin. CCW."""
    pitch_r = module * teeth / 2.0
    base_r = pitch_r * math.cos(pressure)
    tip_r = pitch_r + module
    root_r = pitch_r - 1.25 * module
    t_tip = math.sqrt((tip_r / base_r) ** 2 - 1.0)
    t_pitch = math.tan(pressure)
    px, py = _involute(base_r, t_pitch)
    shift = (math.pi / (2.0 * teeth)) - math.atan2(py, px)

    flank: list[tuple[float, float]] = []
    for i in range(samples + 1):
        x, y = _involute(base_r, t_tip * i / samples)
        flank.append(_rot2(x, y, shift))
    start_r = math.hypot(*flank[0])
    if start_r > root_r:
        ang = math.atan2(flank[0][1], flank[0][0])
        flank.insert(0, (root_r * math.cos(ang), root_r * math.sin(ang)))

    pos = flank
    neg = [(x, -y) for x, y in flank]
    tooth = neg + list(reversed(pos))
    pitch = math.tau / teeth
    points: list[tuple[float, float]] = []
    for i in range(teeth):
        placed = [_rot2(x, y, i * pitch) for x, y in tooth]
        points.extend(placed)
        end = placed[-1]
        nxt = _rot2(tooth[0][0], tooth[0][1], (i + 1) * pitch)
        a0 = math.atan2(end[1], end[0])
        a1 = math.atan2(nxt[1], nxt[0])
        if a1 <= a0:
            a1 += math.tau
        for s in range(1, 4):
            a = a0 + (a1 - a0) * s / 4.0
            points.append((root_r * math.cos(a), root_r * math.sin(a)))
    return points


def face_from_points(points: list[tuple[float, float]]) -> bd.Face:
    area = 0.0
    for i, (x1, y1) in enumerate(points):
        x2, y2 = points[(i + 1) % len(points)]
        area += x1 * y2 - x2 * y1
    if area < 0.0:
        points = list(reversed(points))
    wire = bd.Polyline(*points, close=True)
    return bd.Face(wire.wire())


def extrude_yz(points_yz: list[tuple[float, float]], x0: float, x1: float) -> bd.Solid:
    """Extrude a YZ polygon from x0 to x1. points are (y, z)."""
    sketch = [(-z, y) for y, z in points_yz]
    solid = bd.Rot(0, 90, 0) * bd.extrude(face_from_points(sketch), amount=x1 - x0)
    return bd.Pos(x0, 0, 0) * solid


def extrude_xz(points_xz: list[tuple[float, float]], y0: float, y1: float) -> bd.Solid:
    """Extrude an XZ polygon from y0 to y1. points are (x, z)."""
    sketch = [(x, -z) for x, z in points_xz]
    solid = bd.Rot(-90, 0, 0) * bd.extrude(face_from_points(sketch), amount=y1 - y0)
    return bd.Pos(0, y0, 0) * solid


def cap_screw(x_head_face: float, y: float, z: float, *, head_d: float, head_h: float, shank_d: float, shank_len: float, hex_af: float, direction: float) -> bd.Solid:
    """Socket head along X. direction +1 puts the head toward -X of the face."""
    sign = -1.0 if direction > 0 else 1.0
    hx0 = x_head_face
    hx1 = x_head_face + sign * head_h
    head = cyl_x(head_d / 2.0, min(hx0, hx1), max(hx0, hx1), y, z)
    r = hex_af / math.sqrt(3.0)
    hex_pts = [(r * math.cos(math.radians(30 + i * 60)), r * math.sin(math.radians(30 + i * 60))) for i in range(6)]
    socket = bd.extrude(face_from_points(hex_pts), amount=head_h * 0.55)
    socket = bd.Rot(0, 90, 0) * socket
    socket = bd.Pos(x_head_face + sign * head_h * 0.15, y, z) * socket
    sx0 = x_head_face
    sx1 = x_head_face - sign * shank_len
    shank = cyl_x(shank_d / 2.0, min(sx0, sx1) - 0.3, max(sx0, sx1), y, z)
    return head.cut(socket).fuse(shank)


def set_screw(x: float, y_outer: float, z: float, *, radius: float, length: float) -> bd.Solid:
    """Headless screw along +Y, tip toward the shaft, outer end at y_outer."""
    screw = bd.Rot(90, 0, 0) * bd.Cylinder(radius, length)
    return bd.Pos(x, y_outer - length / 2.0, z) * screw


def key_block(x0: float, x1: float, shaft_r: float, *, width: float = 5.0, t1: float = 3.0, height: float = 5.0) -> bd.Solid:
    z0 = CENTER_Z + shaft_r - t1 + 0.08
    return bd.Pos((x0 + x1) / 2.0, 0, z0 + (height - 0.25) / 2.0) * bd.Box(x1 - x0, width - 0.1, height - 0.25)


def keyway_cut(x0: float, x1: float, shaft_r: float, *, width: float = 5.0, depth: float = 3.0) -> bd.Solid:
    z0 = CENTER_Z + shaft_r - depth
    return bd.Pos((x0 + x1) / 2.0, 0, z0 + (depth + 1.5) / 2.0) * bd.Box(x1 - x0, width, depth + 1.5)


def make_bearing() -> bd.Solid:
    """Simplified open 6003: 17 x 35 x 10, eight balls. Axis along Z, centered."""
    outer = bd.Cylinder(17.5, 10) - bd.Cylinder(14.15, 10.4)
    inner = bd.Cylinder(11.85, 10) - bd.Cylinder(8.5, 10.4)
    balls = []
    for i in range(8):
        angle = i * math.tau / 8.0
        balls.append(bd.Pos(13.0 * math.cos(angle), 13.0 * math.sin(angle), 0) * bd.Sphere(2.9))
    solid = outer.fuse(inner)
    for ball in balls:
        solid = solid.fuse(ball)
    return solid


def place_bearing(x0: float) -> bd.Solid:
    solid = bd.Rot(0, 90, 0) * make_bearing()
    return bd.Pos(x0 + 5.0, 0, CENTER_Z) * solid


def make_shaft() -> bd.Solid:
    parts = [
        cyl_x(8.0, -56.0, 0.0),
        cyl_x(8.5, 0.0, 10.0),
        cyl_x(10.0, 10.0, 32.0),
        cyl_x(8.5, 32.0, 63.0),
        cyl_x(8.0, 63.0, 86.0),
    ]
    shaft = parts[0]
    for part in parts[1:]:
        shaft = shaft.fuse(part)
    shaft = shaft.cut(keyway_cut(-50.0, -30.0, 8.0))
    shaft = shaft.cut(keyway_cut(70.0, 84.0, 8.0))
    return shaft


def make_housing() -> bd.Solid:
    barrel = cyl_x(BARREL_R, 0.0, HOUSING_LEN)
    base = bd.Pos(26.0, 0, BASE_H / 2.0) * bd.Box(70.0, 96.0, BASE_H)
    web = bd.Pos(26.0, 0, 28.0) * bd.Box(46.0, 12.0, 32.0)
    rib_a = extrude_yz([(-42, 12), (42, 12), (14, 46), (-14, 46)], 12.0, 20.0)
    rib_b = extrude_yz([(-42, 12), (42, 12), (14, 46), (-14, 46)], 34.0, 42.0)
    gusset_l = extrude_xz([(4, 14), (24, 14), (4, 42)], -4.0, 4.0)
    gusset_r = extrude_xz([(48, 14), (28, 14), (48, 42)], -4.0, 4.0)
    body = barrel
    for part in (base, web, rib_a, rib_b, gusset_l, gusset_r):
        body = body.fuse(part)

    body = body.cut(cyl_x(17.5, -1.0, 42.0))
    body = body.cut(cyl_x(11.0, 42.0, 45.0))
    body = body.cut(cyl_x(15.0, 45.0, HOUSING_LEN + 0.2))
    for angle in BOLT_ANGLES:
        y, z = bolt_center(angle)
        body = body.cut(cyl_x(2.6, -0.2, 14.0, y, z))
        body = body.cut(cyl_x(2.6, 40.0, HOUSING_LEN + 0.2, y, z))
    for x in (10.0, 42.0):
        for y in (-32.0, 32.0):
            body = body.cut(bd.Pos(x, y, BASE_H / 2.0) * bd.Cylinder(4.5, BASE_H + 2.0))
    return body


def make_cover() -> bd.Solid:
    flange = ring_x(BARREL_R, 9.25, -8.0, 0.0)
    boss = ring_x(17.0, 9.25, -18.0, -8.0)
    cover = flange.fuse(boss)
    for angle in BOLT_ANGLES:
        y, z = bolt_center(angle)
        cover = cover.cut(cyl_x(2.8, -8.2, 0.2, y, z))
    return cover


def make_seal_cover() -> bd.Solid:
    cover = ring_x(BARREL_R, 11.5, 52.0, 57.0)
    for angle in BOLT_ANGLES:
        y, z = bolt_center(angle)
        cover = cover.cut(cyl_x(2.8, 51.8, 57.2, y, z))
    return cover


def make_coupling() -> bd.Solid:
    flange = ring_x(28.0, 8.0, -56.0, -46.0)
    hub = ring_x(14.0, 8.0, -46.0, -22.0)
    coupling = flange.fuse(hub)
    for angle in BOLT_ANGLES:
        y = 21.0 * math.cos(math.radians(angle))
        z = CENTER_Z + 21.0 * math.sin(math.radians(angle))
        coupling = coupling.cut(cyl_x(3.3, -56.2, -45.8, y, z))
    # Radial M4 hole in the hub, along +Y, clear of the top keyway.
    coupling = coupling.cut(bd.Pos(-36.0, 11.0, CENTER_Z) * (bd.Rot(90, 0, 0) * bd.Cylinder(2.2, 8.0)))
    return coupling


def make_gear() -> bd.Solid:
    tip_r = GEAR_MODULE * (GEAR_TEETH + 2) / 2.0
    disk = bd.extrude(face_from_points(spur_gear_points(GEAR_MODULE, GEAR_TEETH, GEAR_PA)), amount=GEAR_FACE_W)
    hub = bd.Pos(0, 0, (GEAR_FACE_W - GEAR_NOSE) / 2.0) * bd.Cylinder(18.0, GEAR_FACE_W + GEAR_NOSE)
    gear = disk.fuse(hub)
    gear = gear.cut(bd.Cylinder(8.0, GEAR_FACE_W + GEAR_NOSE + 4.0))
    # Keyway on local -X so that Rot(0, 90, 0) carries it to world +Z.
    gear = gear.cut(bd.Pos(-9.2, 0, 9.0) * bd.Box(4.6, 5.0, 14.0))
    gear = gear.cut(bd.Pos(0, 12.0, -2.5) * (bd.Rot(90, 0, 0) * bd.Cylinder(2.6, 12.0)))
    gear = bd.Rot(0, 90, 0) * gear
    gear = bd.Pos(GEAR_FACE_X, 0, CENTER_Z) * gear
    # Tip diameter is the modeled requirement. The chord between flank ends is inside it.
    peak = 0.0
    for vertex in gear.vertices():
        peak = max(peak, math.hypot(vertex.Y, vertex.Z - CENTER_Z))
    if abs(peak - tip_r) > 0.05:
        raise RuntimeError(f"gear tip radius {peak:.3f} mm, expected {tip_r:.3f} mm")
    return gear


def paint(solid: bd.Solid, name: str, color) -> bd.Solid:
    solid.label = name
    solid.color = color
    return solid


@step(out="../../STEP/power_transmission/power_transmission.step")
@stl(
    out="../../STL/power_transmission/power_transmission.stl",
    mesh_tolerance=1.0e-3,
    mesh_angular_tolerance=0.20,
)
def power_transmission():
    housing = paint(make_housing(), "housing", CAST)
    cover = paint(make_cover(), "cover", COVER_COLOR)
    seal_cover = paint(make_seal_cover(), "seal_cover", COVER_COLOR)
    shaft = paint(make_shaft(), "shaft", SHAFT_COLOR)
    gear = paint(make_gear(), "spur_gear", GEAR_COLOR)
    coupling = paint(make_coupling(), "coupling", COUPLING_COLOR)
    bearing_l = paint(place_bearing(0.0), "bearing_6003_left", BEARING_COLOR)
    bearing_r = paint(place_bearing(32.0), "bearing_6003_right", BEARING_COLOR)
    seal = paint(ring_x(15.0, 8.6, 45.0, 51.8), "oil_seal", SEAL_COLOR)
    key_gear = paint(key_block(70.0, 84.0, 8.0), "key_gear", KEY_COLOR)
    key_coupling = paint(key_block(-50.0, -30.0, 8.0), "key_coupling", KEY_COLOR)

    bolts = []
    for index, angle in enumerate(BOLT_ANGLES):
        y, z = bolt_center(angle)
        bolts.append(paint(
            cap_screw(-8.0, y, z, head_d=8.5, head_h=5.0, shank_d=4.6, shank_len=12.0, hex_af=4.0, direction=1.0),
            f"cover_bolt_{index + 1}",
            BOLT_COLOR,
        ))
        bolts.append(paint(
            cap_screw(57.0, y, z, head_d=8.5, head_h=5.0, shank_d=4.6, shank_len=12.0, hex_af=4.0, direction=-1.0),
            f"seal_bolt_{index + 1}",
            BOLT_COLOR,
        ))
        fy = 21.0 * math.cos(math.radians(angle))
        fz = CENTER_Z + 21.0 * math.sin(math.radians(angle))
        bolts.append(paint(
            cap_screw(-56.0, fy, fz, head_d=10.0, head_h=6.0, shank_d=5.5, shank_len=10.0, hex_af=5.0, direction=1.0),
            f"coupling_bolt_{index + 1}",
            BOLT_COLOR,
        ))

    set_gear = paint(set_screw(65.5, 18.2, CENTER_Z, radius=2.4, length=9.5), "set_screw_gear", BOLT_COLOR)
    set_coupling = paint(set_screw(-36.0, 14.2, CENTER_Z, radius=1.9, length=5.6), "set_screw_coupling", BOLT_COLOR)

    parts = [
        housing, cover, seal_cover, shaft, gear, coupling,
        bearing_l, bearing_r, seal, key_gear, key_coupling,
        set_gear, set_coupling, *bolts,
    ]
    for part in parts:
        if part.volume <= 0:
            raise RuntimeError(f"{part.label} has non-positive volume")

    base = housing.bounding_box()
    if abs(base.min.Z) > 1e-3:
        raise RuntimeError(f"base bottom at Z={base.min.Z:.4f}, expected 0")
    shaft_bb = shaft.bounding_box()
    if abs(((shaft_bb.min.Z + shaft_bb.max.Z) / 2.0) - CENTER_Z) > 1e-3:
        raise RuntimeError("shaft axis is not at the 70 mm center height")
    for bearing, name in ((bearing_l, "left"), (bearing_r, "right")):
        box = bearing.bounding_box()
        od = box.max.Y - box.min.Y
        width = box.max.X - box.min.X
        if abs(od - 35.0) > 0.05 or abs(width - 10.0) > 0.05:
            raise RuntimeError(f"{name} bearing envelope {od:.2f} x {width:.2f}, expected 35 x 10")

    assembly = bd.Compound(children=parts, label="power_transmission")
    return assembly


if __name__ == "__main__":
    power_transmission()
    print(Path(__file__).resolve().parent.parent.parent / "STL" / "power_transmission" / "power_transmission.stl")
