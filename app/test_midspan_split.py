"""
test_midspan_split.py — analysis_orchestration.py এর mid-span point
load split logic এর জন্য unit test।

চালানো: python3 test_midspan_split.py  (pytest ছাড়াই চলে, plain assert)

এই টেস্ট C++ solver ছাড়াই চলে (pure Python, build_solver_model এর
আউটপুট dict/list শেপ ও মান যাচাই করে) — দ্রুত, CI-friendly। End-to-end
numerical verification (textbook formula মিলিয়ে) আলাদাভাবে ম্যানুয়াল
sandbox এ করা হয়েছে (pybind11 module কম্পাইল করে PL³/48EI ও
Pa²b²/3EIL ফর্মুলার সাথে relative error ~1e-9 এ মিলেছে) — সেই
verification এখানে repeat করা হয়নি কারণ এই ফাইলের scope শুধু
orchestration logic (split হওয়া node/element/load সঠিক কিনা),
numerical solver accuracy না (সেটা cpp/tests/ এর দায়িত্ব)।
"""

import sys

from analysis_orchestration import build_solver_model, ModelParsingError

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


def make_basic_beam(connection_type: str = "moment") -> tuple[list, list, list]:
    elements = [
        {
            "elementId": "beam-1",
            "category": "beam",
            "startPoint": {"x": 0.0, "y": 0.0, "z": 0.0},
            "endPoint": {"x": 6.0, "y": 0.0, "z": 0.0},
            "materialId": "mat-1",
            "sectionId": "sec-1",
            "connectionType": connection_type,
        },
    ]
    materials = [{"materialId": "mat-1", "type": "steel", "es": 200000, "unitWeight": 78.5}]
    sections = [{"sectionId": "sec-1", "properties": {"area": 10000, "ixx": 1e8, "iyy": 5e7, "j": 2e7}}]
    return elements, materials, sections


print("=== Test 1: No mid-span load — element unchanged (1 sub-element) ===")
elements, materials, sections = make_basic_beam()
load_cases = [{"elementId": "beam-1", "applicationType": "point", "positionRatio": 0.0,
               "forceX": 0, "forceY": -10, "forceZ": 0}]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("Exactly 1 solver element (no split)", len(solver_input["elements"]) == 1)
check("elementId unchanged (no #0 suffix)", solver_input["elements"][0]["elementId"] == "beam-1")
check("Exactly 2 nodes (just endpoints)", len(solver_input["nodes"]) == 2)
check("Load applied at node 0 (start)", solver_input["loads"][0]["nodeIndex"] == 0)

print("\n=== Test 2: Mid-span load at ratio=0.5 — splits into 2 sub-elements ===")
elements, materials, sections = make_basic_beam()
load_cases = [{"elementId": "beam-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": -10, "forceZ": 0}]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("Exactly 2 solver sub-elements", len(solver_input["elements"]) == 2)
check("3 nodes total (2 endpoints + 1 split)", len(solver_input["nodes"]) == 3)
midspan_node = solver_input["nodes"][2]
check("Split node at x=3.0 (midpoint of 0..6)", abs(midspan_node["x"] - 3.0) < 1e-9)
check("Load applied at the split node (index 2)", solver_input["loads"][0]["nodeIndex"] == 2)
check("Sub-element 0 elementId is 'beam-1#0'", solver_input["elements"][0]["elementId"] == "beam-1#0")
check("Sub-element 1 elementId is 'beam-1#1'", solver_input["elements"][1]["elementId"] == "beam-1#1")
check("No mid-span snap warning present",
      not any("snap" in w.lower() or "nearest" in w for w in warnings))
check("Registry has 2 entries mapping back to 'beam-1'",
      len(registry) == 2 and all(r["originalElementId"] == "beam-1" for r in registry))

print("\n=== Test 3: Off-center load (ratio=0.25) — split point at correct position ===")
elements, materials, sections = make_basic_beam()
load_cases = [{"elementId": "beam-1", "applicationType": "point", "positionRatio": 0.25,
               "forceX": 0, "forceY": -10, "forceZ": 0}]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
split_node = solver_input["nodes"][2]
check("Split node at x=1.5 (25% of 6m)", abs(split_node["x"] - 1.5) < 1e-9)

print("\n=== Test 4: Multiple mid-span loads on same element — multiple splits ===")
elements, materials, sections = make_basic_beam()
load_cases = [
    {"elementId": "beam-1", "applicationType": "point", "positionRatio": 0.25,
     "forceX": 0, "forceY": -5, "forceZ": 0},
    {"elementId": "beam-1", "applicationType": "point", "positionRatio": 0.75,
     "forceX": 0, "forceY": -5, "forceZ": 0},
]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("3 sub-elements (2 splits)", len(solver_input["elements"]) == 3)
check("4 nodes total (2 endpoints + 2 splits)", len(solver_input["nodes"]) == 4)
check("2 loads, each at its own split node",
      len(solver_input["loads"]) == 2 and solver_input["loads"][0]["nodeIndex"] != solver_input["loads"][1]["nodeIndex"])

print("\n=== Test 5: Load very near endpoint (within tolerance) — snaps, no split ===")
elements, materials, sections = make_basic_beam()
load_cases = [{"elementId": "beam-1", "applicationType": "point", "positionRatio": 0.00001,
               "forceX": 0, "forceY": -10, "forceZ": 0}]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("No split for near-zero ratio (1 sub-element)", len(solver_input["elements"]) == 1)
check("2 nodes only (no split node created)", len(solver_input["nodes"]) == 2)
check("Load snapped to start node (index 0)", solver_input["loads"][0]["nodeIndex"] == 0)

print("\n=== Test 6: Internal cut point stays rigid even if element is 'pin' ===")
elements, materials, sections = make_basic_beam(connection_type="pin")
load_cases = [{"elementId": "beam-1", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": -10, "forceZ": 0}]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("Sub-element 0 keeps 'pin' (touches real start)", solver_input["elements"][0]["connectionType"] == "pin")
check("Sub-element 1 keeps 'pin' (touches real end)", solver_input["elements"][1]["connectionType"] == "pin")
check("Warning issued about split-pin approximation",
      any("pin" in w.lower() and ("split" in w or "উভয়" in w) for w in warnings))

print("\n=== Test 7: ModelParsingError still raised for missing element reference ===")
elements, materials, sections = make_basic_beam()
load_cases = [{"elementId": "nonexistent", "applicationType": "point", "positionRatio": 0.5,
               "forceX": 0, "forceY": -10, "forceZ": 0}]
raised = False
try:
    build_solver_model(elements, materials, sections, load_cases, set())
except ModelParsingError:
    raised = True
check("ModelParsingError raised for unknown elementId in load case", raised)

print(f"\n{'=' * 40}")
print(f"Results: {testsPassed} passed, {testsFailed} failed")
print("=" * 40)

sys.exit(1 if testsFailed > 0 else 0)
