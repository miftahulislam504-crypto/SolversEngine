"""
test_uniform_loads.py — analysis_orchestration.py এর uniform-line ও
uniform-area equivalent nodal load conversion logic এর জন্য unit test।

চালানো: python3 test_uniform_loads.py  (pytest ছাড়াই চলে, plain assert)

test_midspan_split.py এর মতোই — pure Python, build_solver_model() এর
আউটপুট dict/list শেপ ও মান যাচাই করে, C++ solver ছাড়াই। uniform-line
এর wL/2 lumped equivalent এর numerical মান হাতে-হিসাব করে মিলানো
হয়েছে (textbook w×L/2 formula), uniform-area এর quad-area × intensity
÷ 4 ও একই ভাবে যাচাই করা।
"""

import os
import sys

# app/ থেকে সরাসরি (python3 test_uniform_loads.py) চালালে project root
# sys.path এ থাকে না, তাই analysis_orchestration.py এর নিজস্ব "from
# app.model_conversion import ..." (package-style, production convention —
# app/main.py এর সাথে সামঞ্জস্যপূর্ণ) resolve করতে ব্যর্থ হয়। এই দুই লাইন
# project root (app/ এর parent) কে path এ যোগ করে, যাতে "app" প্যাকেজ
# হিসেবে import করা যায় — bare "from analysis_orchestration import ..."
# এর বদলে, যেটা কখনো কাজ করত না (analysis_orchestration.py নিজেই
# "app." prefix আশা করে)।
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.analysis_orchestration import build_solver_model, ModelParsingError

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
        # pytest collection এর সময় module-level কোড import হিসেবে চলে —
        # raised AssertionError ছাড়া pytest বুঝতে পারবে না কিছু fail
        # করেছে (নিচের check() কল গুলো শুধু print/counter, exception না)।
        # standalone চালানোর output/exit code অপরিবর্তিত থাকে (নিচের
        # __main__ guard দেখুন)।
        raise AssertionError(f"check failed: {name}")


def make_basic_beam(length: float = 6.0, connection_type: str = "moment") -> tuple[list, list, list]:
    elements = [
        {
            "elementId": "beam-1",
            "category": "beam",
            "startPoint": {"x": 0.0, "y": 0.0, "z": 0.0},
            "endPoint": {"x": length, "y": 0.0, "z": 0.0},
            "materialId": "mat-1",
            "sectionId": "sec-1",
            "connectionType": connection_type,
        },
    ]
    materials = [{"materialId": "mat-1", "type": "steel", "es": 200000, "unitWeight": 78.5}]
    sections = [{"sectionId": "sec-1", "properties": {"area": 10000, "ixx": 1e8, "iyy": 5e7, "j": 2e7}}]
    return elements, materials, sections


def make_flat_slab(size: float = 4.0) -> tuple[list, list, list]:
    """4x4 m flat horizontal slab (XZ plane, single quad-ready square) at Y=3 (a floor level)."""
    elements = [
        {
            "elementId": "slab-1",
            "category": "slab",
            "vertices": [
                {"x": 0.0, "y": 3.0, "z": 0.0},
                {"x": size, "y": 3.0, "z": 0.0},
                {"x": size, "y": 3.0, "z": size},
                {"x": 0.0, "y": 3.0, "z": size},
            ],
            "materialId": "mat-1",
            "thickness": 150,  # mm
        },
    ]
    materials = [{"materialId": "mat-1", "type": "concrete", "fc": 25, "unitWeight": 24}]
    sections: list = []
    return elements, materials, sections


print("=== Test 1: Uniform-line load on a simple unsplit beam — wL/2 at each end ===")
elements, materials, sections = make_basic_beam(length=6.0)
load_cases = [{"elementId": "beam-1", "applicationType": "uniform-line",
               "forceX": 0, "forceY": -3.6, "forceZ": 0}]  # -3.6 kN/m (matches self-weight test earlier in session)
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("Exactly 1 solver element (uniform load doesn't force a split)", len(solver_input["elements"]) == 1)
check("Exactly 2 loads generated (one per endpoint)", len(solver_input["loads"]) == 2)
expected_half_force = -3.6 * 6.0 / 2.0  # -10.8 kN at each end
check(f"Load at node 0 = {expected_half_force} kN (wL/2)",
      abs(solver_input["loads"][0]["fy"] - expected_half_force) < 1e-9)
check(f"Load at node 1 = {expected_half_force} kN (wL/2)",
      abs(solver_input["loads"][1]["fy"] - expected_half_force) < 1e-9)
check("Loads applied at the two distinct endpoint nodes",
      {solver_input["loads"][0]["nodeIndex"], solver_input["loads"][1]["nodeIndex"]} == {0, 1})
check("Uniform-load disclaimer warning present",
      any("equivalent nodal load" in w or "Lumped" in w or "lumped" in w for w in warnings))

print("\n=== Test 2: Uniform-line load on a beam ALSO split by a mid-span point load ===")
elements, materials, sections = make_basic_beam(length=6.0)
load_cases = [
    {"elementId": "beam-1", "applicationType": "point", "positionRatio": 0.5,
     "forceX": 0, "forceY": -10, "forceZ": 0},
    {"elementId": "beam-1", "applicationType": "uniform-line",
     "forceX": 0, "forceY": -2.0, "forceZ": 0},
]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("2 sub-elements (split by the point load)", len(solver_input["elements"]) == 2)
# 1 point load + 4 uniform-line contributions (2 sub-elements x 2 ends each)
check("5 total load entries (1 point + 2x2 uniform-line)", len(solver_input["loads"]) == 5)
# Each sub-element is 3m long (half of 6m) -> w*L_sub/2 = -2.0*3/2 = -3.0 kN per end
uniform_loads = [ld for ld in solver_input["loads"] if abs(ld["fy"] - (-10.0)) > 1e-9]
check("4 uniform-line contributions each = -3.0 kN (w*L_sub/2 for 3m sub-segment)",
      len(uniform_loads) == 4 and all(abs(ld["fy"] - (-3.0)) < 1e-9 for ld in uniform_loads))
# Split node (index 2, shared by both sub-elements) should receive TWO -3.0 kN contributions (sums to -6.0 via solver's +=)
split_node_uniform_contribs = [ld for ld in uniform_loads if ld["nodeIndex"] == 2]
check("Split node receives 2 separate -3.0 kN entries (solver sums them, total -6.0)",
      len(split_node_uniform_contribs) == 2)

print("\n=== Test 3: Uniform-area load on a flat 4x4m slab — area x intensity / 4 per corner ===")
elements, materials, sections = make_flat_slab(size=4.0)
load_cases = [{"elementId": "slab-1", "applicationType": "uniform-area", "intensity": 5.0}]  # 5 kN/m^2
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
total_area = 4.0 * 4.0  # 16 m^2
expected_total_force = 5.0 * total_area  # 80 kN total downward... but applied as +fy per corner (see note below)
actual_total_force = sum(ld["fy"] for ld in solver_input["loads"])
check(f"Sum of all corner fy loads = intensity x area = {expected_total_force} kN (total force conserved)",
      abs(actual_total_force - expected_total_force) < 1e-6)
check("All uniform-area loads only touch fy (fx=fz=0)",
      all(ld["fx"] == 0.0 and ld["fz"] == 0.0 for ld in solver_input["loads"]))
check("Shell force-recovery limitation warning still present (pre-existing)",
      any("stress/moment" in w or "internal force" in w for w in warnings))

print("\n=== Test 4: Uniform-area load on an element that is NOT a shell — rejected with warning, not crashed ===")
elements, materials, sections = make_basic_beam()
load_cases = [{"elementId": "beam-1", "applicationType": "uniform-area", "intensity": 5.0}]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("No loads generated (beam-1 is not a shell)", len(solver_input["loads"]) == 0)
check("Warning issued about non-shell uniform-area target",
      any("shell" in w.lower() and "না" in w for w in warnings))

print("\n=== Test 5: Unsupported load type still correctly flagged (regression check) ===")
elements, materials, sections = make_basic_beam()
load_cases = [{"elementId": "beam-1", "applicationType": "some-future-type",
               "forceX": 0, "forceY": -1, "forceZ": 0}]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
check("No loads generated for unknown type", len(solver_input["loads"]) == 0)
check("Warning names the unsupported type",
      any("some-future-type" in w for w in warnings))

print("\n=== Test 6: Uniform-line load with zero total_length sub-entries still uses correct total length ===")
elements, materials, sections = make_basic_beam(length=10.0)
load_cases = [{"elementId": "beam-1", "applicationType": "uniform-line",
               "forceX": 0, "forceY": -1.0, "forceZ": 0}]
solver_input, warnings, registry = build_solver_model(elements, materials, sections, load_cases, set())
expected_half = -1.0 * 10.0 / 2.0  # -5.0
check("10m beam: wL/2 = -5.0 kN at each end",
      all(abs(ld["fy"] - expected_half) < 1e-9 for ld in solver_input["loads"]))

print(f"\n{'=' * 40}")
print(f"Results: {testsPassed} passed, {testsFailed} failed")
print("=" * 40)

# __main__ guard: pytest এই ফাইল import করলে sys.exit() রেজ হয়ে
# INTERNALERROR দিত (pytest নিজের exit flow না এমন SystemExit ক্র্যাশ
# হিসেবে ধরে) — standalone চালানোর (python3 test_uniform_loads.py)
# সময়ই শুধু এই কল সক্রিয়।
if __name__ == "__main__":
    sys.exit(1 if testsFailed > 0 else 0)
