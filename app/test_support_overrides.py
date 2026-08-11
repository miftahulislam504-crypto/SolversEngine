"""
test_support_overrides.py — analysis_orchestration.py এর Phase 5
support_overrides logic এর জন্য unit test।

চালানো: python3 test_support_overrides.py  (pytest ছাড়াই চলে, plain
assert — test_midspan_split.py এর একই কনভেনশন)।

এই টেস্ট C++ solver ছাড়াই চলে (pure Python, build_solver_model এর
boundary_conditions আউটপুট list শেপ ও মান যাচাই করে)।
"""

import sys

from analysis_orchestration import build_solver_model

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


def make_two_story_column() -> tuple[list, list, list]:
    # Y=0 (base, heuristic ধরবে) → Y=3 (mid, ধরবে না) → Y=6 (roof, ধরবে না)
    elements = [
        {
            "elementId": "col-1",
            "category": "column",
            "startPoint": {"x": 0.0, "y": 0.0, "z": 0.0},
            "endPoint": {"x": 0.0, "y": 3.0, "z": 0.0},
            "materialId": "mat-1",
            "sectionId": "sec-1",
            "connectionType": "moment",
        },
        {
            "elementId": "col-2",
            "category": "column",
            "startPoint": {"x": 0.0, "y": 3.0, "z": 0.0},
            "endPoint": {"x": 0.0, "y": 6.0, "z": 0.0},
            "materialId": "mat-1",
            "sectionId": "sec-1",
            "connectionType": "moment",
        },
    ]
    materials = [{"materialId": "mat-1", "type": "steel", "es": 200000, "unitWeight": 78.5}]
    sections = [{"sectionId": "sec-1", "properties": {"area": 10000, "ixx": 1e8, "iyy": 5e7, "j": 2e7}}]
    return elements, materials, sections


def find_bc_for_node(boundary_conditions: list, node_index: int) -> dict | None:
    for bc in boundary_conditions:
        if bc["nodeIndex"] == node_index:
            return bc
    return None


print("=== Test 1: No override — Y≈0 heuristic behaves as before (backward compatible) ===")
elements, materials, sections = make_two_story_column()
load_cases = [{"elementId": "col-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": 0, "forceZ": -10}]
solver_input, warnings, _ = build_solver_model(elements, materials, sections, load_cases, set())
base_node = next(i for i, n in enumerate(solver_input["nodes"]) if n["y"] == 0.0)
bc = find_bc_for_node(solver_input["boundaryConditions"], base_node)
check("Base node (Y=0) has a boundary condition", bc is not None)
check("Base node is fully fixed (all 6 DOF)", bc is not None and all(bc[k] for k in
      ["restrainX", "restrainY", "restrainZ", "restrainRx", "restrainRy", "restrainRz"]))
check("Exactly 1 boundary condition (no override, only base node)", len(solver_input["boundaryConditions"]) == 1)
check("No support_override warning present when support_overrides=None",
      not any("support_override" in w for w in warnings))

print("\n=== Test 2: supportType='pinned' override at base — translations fixed, rotations free ===")
elements, materials, sections = make_two_story_column()
load_cases = [{"elementId": "col-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": 0, "forceZ": -10}]
overrides = [{"x": 0.0, "y": 0.0, "z": 0.0, "supportType": "pinned"}]
solver_input, warnings, _ = build_solver_model(elements, materials, sections, load_cases, set(), overrides)
base_node = next(i for i, n in enumerate(solver_input["nodes"]) if n["y"] == 0.0)
bc = find_bc_for_node(solver_input["boundaryConditions"], base_node)
check("Base node still has a boundary condition after override", bc is not None)
check("Pinned: translations restrained", bc is not None and bc["restrainX"] and bc["restrainY"] and bc["restrainZ"])
check("Pinned: rotations free", bc is not None and not bc["restrainRx"] and not bc["restrainRy"] and not bc["restrainRz"])
check("Applied-override warning present", any("প্রয়োগ করা হয়েছে" in w for w in warnings))

print("\n=== Test 3: supportType='free' override removes an existing heuristic support ===")
elements, materials, sections = make_two_story_column()
load_cases = [{"elementId": "col-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": 0, "forceZ": -10}]
overrides = [{"x": 0.0, "y": 0.0, "z": 0.0, "supportType": "free"}]
solver_input, warnings, _ = build_solver_model(elements, materials, sections, load_cases, set(), overrides)
base_node = next(i for i, n in enumerate(solver_input["nodes"]) if n["y"] == 0.0)
bc = find_bc_for_node(solver_input["boundaryConditions"], base_node)
check("Base node boundary condition removed (supportType='free')", bc is None)
check("All-removed warning present (structure now unstable)",
      any("সব boundary condition override দিয়ে সরিয়ে" in w for w in warnings))

print("\n=== Test 4: supportType='custom' with explicit per-DOF flags ===")
elements, materials, sections = make_two_story_column()
load_cases = [{"elementId": "col-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": 0, "forceZ": -10}]
overrides = [{"x": 0.0, "y": 0.0, "z": 0.0, "supportType": "custom",
              "restrainX": True, "restrainY": True, "restrainZ": False,
              "restrainRx": False, "restrainRy": True, "restrainRz": False}]
solver_input, warnings, _ = build_solver_model(elements, materials, sections, load_cases, set(), overrides)
base_node = next(i for i, n in enumerate(solver_input["nodes"]) if n["y"] == 0.0)
bc = find_bc_for_node(solver_input["boundaryConditions"], base_node)
check("Custom DOF flags applied exactly as specified",
      bc is not None and bc["restrainX"] and bc["restrainY"] and not bc["restrainZ"]
      and not bc["restrainRx"] and bc["restrainRy"] and not bc["restrainRz"])

print("\n=== Test 5: Override at a non-base coordinate (mid-height, Y=3) — supplements heuristic, doesn't replace it ===")
elements, materials, sections = make_two_story_column()
load_cases = [{"elementId": "col-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": 0, "forceZ": -10}]
overrides = [{"x": 0.0, "y": 3.0, "z": 0.0, "supportType": "pinned"}]
solver_input, warnings, _ = build_solver_model(elements, materials, sections, load_cases, set(), overrides)
base_node = next(i for i, n in enumerate(solver_input["nodes"]) if n["y"] == 0.0)
mid_node = next(i for i, n in enumerate(solver_input["nodes"]) if n["y"] == 3.0)
base_bc = find_bc_for_node(solver_input["boundaryConditions"], base_node)
mid_bc = find_bc_for_node(solver_input["boundaryConditions"], mid_node)
check("Base node (Y=0) still fully fixed via heuristic (untouched)",
      base_bc is not None and all(base_bc[k] for k in ["restrainX", "restrainY", "restrainZ", "restrainRx", "restrainRy", "restrainRz"]))
check("Mid-height node (Y=3) now also has a boundary condition (pinned, via override)",
      mid_bc is not None and mid_bc["restrainX"] and not mid_bc["restrainRx"])
check("Total boundary conditions = 2 (base heuristic + mid override)", len(solver_input["boundaryConditions"]) == 2)

print("\n=== Test 6: Override at a coordinate with no matching node — unmatched, warned, not silently ignored ===")
elements, materials, sections = make_two_story_column()
load_cases = [{"elementId": "col-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": 0, "forceZ": -10}]
overrides = [{"x": 99.0, "y": 99.0, "z": 99.0, "supportType": "fixed"}]
solver_input, warnings, _ = build_solver_model(elements, materials, sections, load_cases, set(), overrides)
check("Unmatched-override warning present", any("মেলেনি" in w for w in warnings))
check("Base node heuristic still applies (unaffected by the unmatched override)",
      len(solver_input["boundaryConditions"]) == 1)

print("\n=== Test 7: Unknown supportType string — override skipped with warning, doesn't crash ===")
elements, materials, sections = make_two_story_column()
load_cases = [{"elementId": "col-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": 0, "forceZ": -10}]
overrides = [{"x": 0.0, "y": 0.0, "z": 0.0, "supportType": "spring"}]
solver_input, warnings, _ = build_solver_model(elements, materials, sections, load_cases, set(), overrides)
check("Unknown supportType warning present", any("অজানা supportType" in w for w in warnings))
base_node = next(i for i, n in enumerate(solver_input["nodes"]) if n["y"] == 0.0)
bc = find_bc_for_node(solver_input["boundaryConditions"], base_node)
check("Base node heuristic untouched by the rejected override",
      bc is not None and all(bc[k] for k in ["restrainX", "restrainY", "restrainZ", "restrainRx", "restrainRy", "restrainRz"]))

print("\n=== Test 8: Empty support_overrides list behaves identically to None ===")
elements, materials, sections = make_two_story_column()
load_cases = [{"elementId": "col-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": 0, "forceZ": -10}]
solver_input_none, warnings_none, _ = build_solver_model(elements, materials, sections, load_cases, set(), None)
solver_input_empty, warnings_empty, _ = build_solver_model(elements, materials, sections, load_cases, set(), [])
check("Same boundary conditions with None vs [] support_overrides",
      solver_input_none["boundaryConditions"] == solver_input_empty["boundaryConditions"])

print("\n" + "=" * 40)
print(f"Results: {testsPassed} passed, {testsFailed} failed")
print("=" * 40)

if testsFailed > 0:
    sys.exit(1)
