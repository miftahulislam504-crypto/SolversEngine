"""
Mesh Generation — General Polygon → Quad Mesh
==================================================

Slab/Wall element গুলো frontend থেকে একটা general (possibly-concave)
polygon হিসেবে আসে (vertices: Point3D[], ন্যূনতম ৩টা, counter-clockwise)।
আমাদের shell element (cpp/src/shell.cpp) শুধু 4-node quadrilateral
সমর্থন করে, তাই এই polygon-কে quad-এ mesh করতে হয়।

পদ্ধতি (দুই ধাপ):
  ১. Ear-Clipping Triangulation — general simple polygon (concave
     সহ, কিন্তু self-intersecting না) কে triangle এ ভাগ করা। এটা একটা
     প্রতিষ্ঠিত, সরল O(n²) algorithm (Meister, 1975 থেকে উদ্ভূত ধারণা,
     কম্পিউটেশনাল জ্যামিতিতে ব্যাপকভাবে পড়ানো ও ব্যবহৃত) — n vertex
     পর্যন্ত ছোট/মাঝারি polygon (slab/wall boundary, সাধারণত <20
     vertex) এর জন্য যথেষ্ট দ্রুত, আরও efficient algorithm
     (Delaunay-based ইত্যাদি) এই স্কেলে প্রয়োজন নেই।
  ২. Triangle → Quad রূপান্তর — প্রতিটা triangle কে তার centroid ও
     তিনটা edge-midpoint ব্যবহার করে ৩টা quad এ ভাগ করা হয় (একটা
     সুপরিচিত, সহজ triangle-to-quad conversion কৌশল — J-এর প্রতিটা
     triangle থেকে "fan of quads around the centroid" তৈরি করে)। এটা
     Delaunay-based advancing-front quad meshing এর চেয়ে সরল কিন্তু
     robust (কখনো invalid/self-intersecting quad তৈরি করে না, যেকোনো
     non-degenerate triangle এর জন্য)।

সীমাবদ্ধতা (এই সংস্করণ):
  - Mesh density uniform/adaptive না — শুধু polygon boundary অনুযায়ী
    triangulate করে প্রতিটা triangle কে ৩-quad এ ভাগ করে, কোনো
    additional mesh refinement (target element size অনুযায়ী subdivide
    করা) নেই। ফলে বড় slab এ কম, uneven-size element হবে — Phase 4
    এর এই ধাপের জন্য এটা গ্রহণযোগ্য (displacement সমাধান পাওয়া মূল
    লক্ষ্য), কিন্তু আরও নিখুঁত stress recovery ভবিষ্যতে uniform mesh
    density দাবি করবে।
  - Self-intersecting (non-simple) polygon সমর্থিত না — ear-clipping
    ধরে নেয় polygon simple (edges একে অপরকে ছেদ করে না)।
  - Hole-সহ polygon (একটা slab-এর মধ্যে opening/void) সমর্থিত না এই
    সংস্করণে।
"""

import math
from typing import Any


class MeshGenerationError(Exception):
    """Polygon triangulation বা mesh generation ব্যর্থ হলে (invalid geometry)।"""


Point = tuple[float, float, float]


def _cross_2d(o: Point, a: Point, b: Point) -> float:
    """
    ২D cross product (XZ প্লেনে প্রজেক্ট করে — frontend এর
    computePolygonPlanArea() এর মতো একই convention, element.ts এ
    দেখুন) — polygon winding direction ও ear-validity যাচাই করতে
    ব্যবহৃত। XY না, XZ ব্যবহারের কারণ: frontend এ "plan" (horizontal)
    এরিয়া XZ প্লেনে প্রজেক্ট করে গণনা করা হয় (Y = উচ্চতা/vertical axis
    কনভেনশন অনুযায়ী, CivilOS এর পুরো geometry system এ)। কিন্তু একটা
    Wall vertical plane এ থাকতে পারে (XY বা YZ প্লেনে) — তাই এই ফাংশন
    আসলে polygon যে প্লেনেই থাকুক না কেন কাজ করার জন্য একটা local ২D
    projection দরকার, যা _project_to_2d() এ করা হয় triangulate() এর
    ভেতরে, এখানে সরাসরি না।
    """
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _project_polygon_to_2d(vertices: list[Point]) -> tuple[list[tuple[float, float]], Point, Point, Point]:
    """
    ৩D polygon vertex তালিকা (যা একটা সমতলে অবস্থিত বলে ধরা হয়, flat
    slab/wall — non-planar polygon সমর্থিত না) কে একটা local ২D
    coordinate system এ project করে, যাতে standard ২D ear-clipping
    algorithm প্রয়োগ করা যায়।

    রিটার্ন করে: (2d_points, origin, local_u_axis, local_v_axis) —
    পরবর্তীতে ২D triangulation ফলাফলকে আবার ৩D তে ফিরিয়ে আনতে এই
    axis গুলো লাগবে (unproject করতে)।
    """
    if len(vertices) < 3:
        raise MeshGenerationError(f"Polygon এ ন্যূনতম ৩টা vertex দরকার, পাওয়া গেছে {len(vertices)}টা")

    p0 = vertices[0]

    # local U axis = vertex0 → vertex1 দিকের unit vector
    p1 = vertices[1]
    u_raw = (p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2])
    u_len = math.sqrt(u_raw[0] ** 2 + u_raw[1] ** 2 + u_raw[2] ** 2)
    if u_len < 1e-9:
        raise MeshGenerationError("Polygon এর প্রথম দুইটা vertex প্রায় coincident — degenerate geometry")
    u_axis = (u_raw[0] / u_len, u_raw[1] / u_len, u_raw[2] / u_len)

    # normal (local W axis) = প্রথম non-degenerate cross product খুঁজে বের করা
    # (শুধু vertex0,1,2 ব্যবহার না করে — যদি সেগুলো collinear হয়, পরের
    # vertex চেষ্টা করা)
    normal = None
    for i in range(2, len(vertices)):
        pi = vertices[i]
        v_raw = (pi[0] - p0[0], pi[1] - p0[1], pi[2] - p0[2])
        cross = (
            u_axis[1] * v_raw[2] - u_axis[2] * v_raw[1],
            u_axis[2] * v_raw[0] - u_axis[0] * v_raw[2],
            u_axis[0] * v_raw[1] - u_axis[1] * v_raw[0],
        )
        cross_len = math.sqrt(cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2)
        if cross_len > 1e-9:
            normal = (cross[0] / cross_len, cross[1] / cross_len, cross[2] / cross_len)
            break

    if normal is None:
        raise MeshGenerationError("Polygon এর সব vertex collinear মনে হচ্ছে — degenerate (zero-area) geometry")

    # local V axis = normal × u_axis (orthogonal সম্পূর্ণ করতে)
    v_axis = (
        normal[1] * u_axis[2] - normal[2] * u_axis[1],
        normal[2] * u_axis[0] - normal[0] * u_axis[2],
        normal[0] * u_axis[1] - normal[1] * u_axis[0],
    )

    points_2d = []
    for p in vertices:
        rel = (p[0] - p0[0], p[1] - p0[1], p[2] - p0[2])
        u_coord = rel[0] * u_axis[0] + rel[1] * u_axis[1] + rel[2] * u_axis[2]
        v_coord = rel[0] * v_axis[0] + rel[1] * v_axis[1] + rel[2] * v_axis[2]
        points_2d.append((u_coord, v_coord))

    return points_2d, p0, u_axis, v_axis


def _unproject_2d_to_3d(
    point_2d: tuple[float, float], origin: Point, u_axis: Point, v_axis: Point
) -> Point:
    """local ২D (u,v) কোঅর্ডিনেটকে আবার গ্লোবাল ৩D তে ফিরিয়ে আনে।"""
    u, v = point_2d
    return (
        origin[0] + u * u_axis[0] + v * v_axis[0],
        origin[1] + u * u_axis[1] + v * v_axis[1],
        origin[2] + u * u_axis[2] + v * v_axis[2],
    )


def _polygon_signed_area_2d(points: list[tuple[float, float]]) -> float:
    """Shoelace formula — ধনাত্মক মানে counter-clockwise winding।"""
    area = 0.0
    n = len(points)
    for i in range(n):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    return area / 2.0


def _point_in_triangle_2d(
    p: tuple[float, float], a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]
) -> bool:
    """একটা point triangle (a,b,c) এর ভিতরে (বা edge-এ) আছে কিনা — barycentric sign test।"""
    def sign(p1, p2, p3):
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

    d1 = sign(p, a, b)
    d2 = sign(p, b, c)
    d3 = sign(p, c, a)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not (has_neg and has_pos)


def _ear_clip_triangulate_2d(points: list[tuple[float, float]]) -> list[tuple[int, int, int]]:
    """
    Ear-Clipping Triangulation — একটা simple (self-intersecting না)
    ২D polygon কে triangle এ ভাগ করে, vertex index এর tuple হিসেবে
    ফেরত দেয় (মূল points তালিকার index অনুযায়ী)।

    Algorithm: প্রতি iteration এ একটা "ear" (একটা vertex যেখানে তার
    দুই প্রতিবেশী দিয়ে গঠিত triangle সম্পূর্ণ polygon এর ভিতরে থাকে,
    এবং সেই triangle এ polygon এর অন্য কোনো vertex নেই) খুঁজে সেটা
    কেটে ফেলা (triangle হিসেবে রেকর্ড করে সেই vertex বাদ দেওয়া), যতক্ষণ
    না মাত্র ৩টা vertex অবশিষ্ট থাকে।
    """
    n = len(points)
    if n < 3:
        raise MeshGenerationError(f"Triangulation এর জন্য ন্যূনতম ৩টা vertex দরকার, পাওয়া গেছে {n}টা")
    if n == 3:
        return [(0, 1, 2)]

    # Winding direction ঠিক করা — ear-clipping algorithm ধরে নেয়
    # counter-clockwise polygon; clockwise হলে vertex index reverse
    # করে নেওয়া হচ্ছে (triangulation শেষে original index এ ফিরিয়ে
    # আনা হবে)।
    signed_area = _polygon_signed_area_2d(points)
    indices = list(range(n))
    if signed_area < 0:
        indices = list(reversed(indices))

    remaining = list(indices)
    triangles: list[tuple[int, int, int]] = []

    # Safety limit — infinite loop প্রতিরোধ করতে (কোনো valid ear খুঁজে
    # না পেলে, যেমন self-intersecting polygon এ ঘটতে পারে)
    max_iterations = n * n + 10
    iteration = 0

    while len(remaining) > 3 and iteration < max_iterations:
        iteration += 1
        ear_found = False

        for i in range(len(remaining)):
            prev_idx = remaining[(i - 1) % len(remaining)]
            curr_idx = remaining[i]
            next_idx = remaining[(i + 1) % len(remaining)]

            a, b, c = points[prev_idx], points[curr_idx], points[next_idx]

            # Convexity check: এই vertex এ polygon convex (ear হওয়ার
            # প্রথম শর্ত) — counter-clockwise polygon এ cross product
            # ধনাত্মক হওয়া উচিত
            cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            if cross <= 1e-12:
                continue  # reflex vertex, ear না

            # কোনো অন্য (non-adjacent) vertex এই triangle এর ভিতরে আছে কিনা
            contains_other_vertex = False
            for j in remaining:
                if j in (prev_idx, curr_idx, next_idx):
                    continue
                if _point_in_triangle_2d(points[j], a, b, c):
                    contains_other_vertex = True
                    break

            if contains_other_vertex:
                continue

            # এটা একটা বৈধ ear — কেটে ফেলা
            triangles.append((prev_idx, curr_idx, next_idx))
            remaining.pop(i)
            ear_found = True
            break

        if not ear_found:
            raise MeshGenerationError(
                "Polygon triangulation ব্যর্থ হয়েছে — কোনো বৈধ 'ear' খুঁজে পাওয়া যায়নি। এটা সাধারণত "
                "মানে polygon self-intersecting (edge গুলো একে অপরকে ছেদ করছে) অথবা degenerate "
                "(duplicate/collinear vertex সহ) — polygon geometry পুনর্বিবেচনা করুন।"
            )

    if len(remaining) == 3:
        triangles.append((remaining[0], remaining[1], remaining[2]))

    return triangles


def _triangle_to_quads_2d(
    a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]
) -> list[tuple[tuple[float, float], tuple[float, float], tuple[float, float], tuple[float, float]]]:
    """
    একটা triangle (a,b,c) কে centroid ও তিনটা edge-midpoint ব্যবহার
    করে ৩টা quad এ ভাগ করে। প্রতিটা quad = [vertex, edge-midpoint,
    centroid, edge-midpoint] (counter-clockwise)।

    এই কৌশল সবসময় valid, non-self-intersecting quad দেয় যেকোনো
    non-degenerate (নন-জিরো-এরিয়া) triangle এর জন্য — কারণ centroid
    সবসময় triangle এর ভিতরে থাকে (convex shape), ও edge-midpoint
    সংজ্ঞা অনুযায়ীই সংশ্লিষ্ট edge এর উপর থাকে।
    """
    centroid = ((a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0)
    mid_ab = ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)
    mid_bc = ((b[0] + c[0]) / 2.0, (b[1] + c[1]) / 2.0)
    mid_ca = ((c[0] + a[0]) / 2.0, (c[1] + a[1]) / 2.0)

    return [
        (a, mid_ab, centroid, mid_ca),
        (b, mid_bc, centroid, mid_ab),
        (c, mid_ca, centroid, mid_bc),
    ]


def generate_quad_mesh(vertices_3d: list[dict[str, float]]) -> list[list[Point]]:
    """
    একটা polygon (frontend Point3D dict এর তালিকা হিসেবে, {"x":.., "y":..,
    "z":..}) কে quad panel এর তালিকায় mesh করে — প্রতিটা quad ৪টা ৩D
    corner point (dict না, plain tuple) এর তালিকা, counter-clockwise।

    এই ফাংশন module-level docstring এ বর্ণিত দুই ধাপ (triangulate,
    triangle→quad) একত্রিত করে একটা পূর্ণ pipeline হিসেবে চালায়।
    """
    vertices: list[Point] = [(v["x"], v["y"], v["z"]) for v in vertices_3d]

    points_2d, origin, u_axis, v_axis = _project_polygon_to_2d(vertices)
    triangle_indices = _ear_clip_triangulate_2d(points_2d)

    quads_3d: list[list[Point]] = []
    for (i0, i1, i2) in triangle_indices:
        a2d, b2d, c2d = points_2d[i0], points_2d[i1], points_2d[i2]
        quads_2d = _triangle_to_quads_2d(a2d, b2d, c2d)
        for quad in quads_2d:
            quad_3d = [_unproject_2d_to_3d(pt, origin, u_axis, v_axis) for pt in quad]
            quads_3d.append(quad_3d)

    return quads_3d
