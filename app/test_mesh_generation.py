"""
test_mesh_generation.py — mesh_generation.py এর জন্য unit test।
চালানো: python3 test_mesh_generation.py
"""

import math
import sys

from mesh_generation import generate_quad_mesh, MeshGenerationError

testsPassed = 0
testsFailed = 0


def check(name: str, condition: bool) -> None:
    global testsPassed, testsFailed
    if condition:
        print(f"  [PASS] {name}")
        testsPassed += 1
    else:
        print(f"  [FAIL] {name}")
        testsFailed += 1


def quad_area_3d(q: list[tuple[float, float, float]]) -> float:
    """একটা planar quad এর area — cross-product-sum পদ্ধতি (Newell's method এর সরলীকৃত রূপ)।"""
    total = [0.0, 0.0, 0.0]
    for i in range(4):
        a, b = q[i], q[(i + 1) % 4]
        cross = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
        total[0] += cross[0]
        total[1] += cross[1]
        total[2] += cross[2]
    return math.sqrt(total[0] ** 2 + total[1] ** 2 + total[2] ** 2) / 2.0


print("=== Test 1: Rectangle (4 vertices, horizontal XZ plane) ===")
rect = [
    {"x": 0.0, "y": 3.0, "z": 0.0},
    {"x": 4.0, "y": 3.0, "z": 0.0},
    {"x": 4.0, "y": 3.0, "z": 3.0},
    {"x": 0.0, "y": 3.0, "z": 3.0},
]
quads = generate_quad_mesh(rect)
check("Rectangle produces 6 quads (2 triangles x 3)", len(quads) == 6)
total_area = sum(quad_area_3d(q) for q in quads)
check("Total mesh area matches rectangle area (4x3=12)", abs(total_area - 12.0) < 1e-6)

print("\n=== Test 2: Concave L-shape (6 vertices) ===")
l_shape = [
    {"x": 0.0, "y": 0.0, "z": 0.0},
    {"x": 4.0, "y": 0.0, "z": 0.0},
    {"x": 4.0, "y": 0.0, "z": 2.0},
    {"x": 2.0, "y": 0.0, "z": 2.0},
    {"x": 2.0, "y": 0.0, "z": 4.0},
    {"x": 0.0, "y": 0.0, "z": 4.0},
]
quads_l = generate_quad_mesh(l_shape)
check("L-shape produces 12 quads (4 triangles x 3)", len(quads_l) == 12)
total_area_l = sum(quad_area_3d(q) for q in quads_l)
check("Total mesh area matches L-shape area (16-4=12)", abs(total_area_l - 12.0) < 1e-6)

print("\n=== Test 3: Vertical Wall (XY plane) ===")
wall = [
    {"x": 0.0, "y": 0.0, "z": 0.0},
    {"x": 5.0, "y": 0.0, "z": 0.0},
    {"x": 5.0, "y": 3.0, "z": 0.0},
    {"x": 0.0, "y": 3.0, "z": 0.0},
]
quads_wall = generate_quad_mesh(wall)
total_area_wall = sum(quad_area_3d(q) for q in quads_wall)
check("Vertical wall (non-horizontal plane) mesh area matches (5x3=15)",
      abs(total_area_wall - 15.0) < 1e-6)

print("\n=== Test 4: Minimal Triangle (3 vertices) ===")
tri = [
    {"x": 0.0, "y": 0.0, "z": 0.0},
    {"x": 4.0, "y": 0.0, "z": 0.0},
    {"x": 2.0, "y": 0.0, "z": 3.0},
]
quads_tri = generate_quad_mesh(tri)
check("Triangle produces exactly 3 quads", len(quads_tri) == 3)
total_area_tri = sum(quad_area_3d(q) for q in quads_tri)
check("Triangle mesh area matches (0.5*4*3=6)", abs(total_area_tri - 6.0) < 1e-6)

print("\n=== Test 5: Degenerate (Collinear) Polygon Raises Error ===")
raised = False
try:
    generate_quad_mesh([
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 0.0, "z": 0.0},
        {"x": 2.0, "y": 0.0, "z": 0.0},
    ])
except MeshGenerationError:
    raised = True
check("Collinear (zero-area) polygon raises MeshGenerationError", raised)

print("\n=== Test 6: Too-Few-Vertices Polygon Raises Error ===")
raised2 = False
try:
    generate_quad_mesh([
        {"x": 0.0, "y": 0.0, "z": 0.0},
        {"x": 1.0, "y": 0.0, "z": 0.0},
    ])
except MeshGenerationError:
    raised2 = True
check("Two-vertex 'polygon' raises MeshGenerationError", raised2)

print("\n=== Test 7: All Quads Have Exactly 4 Distinct Corners ===")
all_valid = True
for q in quads_l:
    if len(q) != 4:
        all_valid = False
        break
    # check no two corners are coincident (degenerate quad)
    for i in range(4):
        for j in range(i + 1, 4):
            dist = math.sqrt(sum((q[i][k] - q[j][k]) ** 2 for k in range(3)))
            if dist < 1e-9:
                all_valid = False
check("All generated quads have 4 distinct, non-coincident corners", all_valid)

print("\n=== Test 8: Pentagon (5 vertices, convex) ===")
# Regular-ish pentagon in the XZ plane
pentagon = [
    {"x": 2.0, "y": 0.0, "z": 0.0},
    {"x": 4.0, "y": 0.0, "z": 1.5},
    {"x": 3.2, "y": 0.0, "z": 4.0},
    {"x": 0.8, "y": 0.0, "z": 4.0},
    {"x": 0.0, "y": 0.0, "z": 1.5},
]
quads_pentagon = generate_quad_mesh(pentagon)
check("Pentagon (5 vertices) produces 9 quads (3 triangles x 3)", len(quads_pentagon) == 9)

print(f"\n{'=' * 40}")
print(f"Results: {testsPassed} passed, {testsFailed} failed")
print("=" * 40)

sys.exit(1 if testsFailed > 0 else 0)
