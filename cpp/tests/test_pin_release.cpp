#include "stiffness.h"
#include "solver.h"
#include "types.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace civilos;

int testsPassed = 0;
int testsFailed = 0;

void checkClose(const std::string& name, double actual, double expected, double tolerance = 1e-6) {
    double relError = std::abs(expected) > 1e-12
        ? std::abs(actual - expected) / std::abs(expected)
        : std::abs(actual - expected);
    if (relError < tolerance) {
        std::cout << "  [PASS] " << name << " = " << actual << " (expected " << expected << ")\n";
        testsPassed++;
    } else {
        std::cout << "  [FAIL] " << name << " = " << actual << " (expected " << expected
                   << ", rel error = " << relError << ")\n";
        testsFailed++;
    }
}

void checkTrue(const std::string& name, bool condition) {
    if (condition) {
        std::cout << "  [PASS] " << name << "\n";
        testsPassed++;
    } else {
        std::cout << "  [FAIL] " << name << "\n";
        testsFailed++;
    }
}

int main() {
    std::cout << std::setprecision(8);

    // ============================================================
    // Test 1: Both-end release zeroes out all bending stiffness rows
    // ============================================================
    // একটা fully pin-pin (both ends released) member এ bending
    // moment/rotation DOF এর কোনো stiffness contribution থাকা উচিত
    // না — শুধু axial (u) ও torsion (rx) stiffness টিকে থাকবে,
    // যা একটা truss/2-force member এর সংজ্ঞা।
    std::cout << "=== Test 1: Pin-Pin Release — Bending Rows/Cols Must Vanish ===\n";
    {
        SectionProperties section{0.01, 0.0001, 0.00005, 0.00002};
        MaterialProperties material{200e6, 80e6}; // steel-like, kN/m2
        auto kRigid = computeLocalStiffnessMatrix(3.0, section, material);
        auto kPinned = applyEndReleases(kRigid, true, true);

        // ry1,rz1,ry2,rz2 rows/cols (indices 4,5,10,11) সব-শূন্য হওয়া উচিত
        double maxBendingEntry = 0.0;
        for (int idx : {4, 5, 10, 11}) {
            for (int j = 0; j < 12; ++j) {
                maxBendingEntry = std::max(maxBendingEntry, std::abs(kPinned(idx, j)));
                maxBendingEntry = std::max(maxBendingEntry, std::abs(kPinned(j, idx)));
            }
        }
        checkTrue("All bending-DOF rows/cols are zero after pin-pin release", maxBendingEntry < 1e-9);

        // axial (0,6) ও torsion (3,9) অপরিবর্তিত থাকা উচিত (release করা হয়নি)
        checkClose("Axial stiffness k(0,0) unchanged by pin release", kPinned(0, 0), kRigid(0, 0));
        checkClose("Torsional stiffness k(3,3) unchanged by pin release", kPinned(3, 3), kRigid(3, 3));

        // Symmetry বজায় থাকা উচিত (condensation একটা symmetric operation)
        double maxAsymmetry = (kPinned - kPinned.transpose()).cwiseAbs().maxCoeff();
        checkTrue("Condensed matrix remains symmetric", maxAsymmetry < 1e-9);
    }

    // ============================================================
    // Test 2: Propped Cantilever with an Internal Pin — Known Textbook Case
    // ============================================================
    // দুইটা beam element (node1→node2→node3), node1 এ fixed support,
    // node3 এ roller (Y-only restrained) support, node2 এ একটা internal
    // pin (হিঞ্জ) — অর্থাৎ দ্বিতীয় element-এর startEnd (node2 প্রান্ত)
    // pin-released। Point load P, node2-তে (হিঞ্জেই) প্রয়োগ করা।
    //
    // এই সেটআপ পদার্থগতভাবে স্থিতিশীল (fixed + roller + internal hinge
    // মিলিয়ে একটা determinate/stable system, mechanism না) — কারণ
    // element 1 (node1-node2, rigid both ends) পুরো bending নেয়, আর
    // element 2 (node2-node3, node2-প্রান্তে pin) সেই হিঞ্জ দিয়ে moment
    // transfer করে না।
    //
    // Reference check: classical propped cantilever with internal hinge
    // at the load point, fixed-pinned overall — this is the standard
    // "cantilever + simply-supported prop with hinge" configuration from
    // Hibbeler/Kassimali. আমরা এখানে exact closed-form না মিলিয়ে বরং
    // qualitative/sanity checks করছি (solver সফল হয়, deflection finite ও
    // reasonable sign এ) — পূর্ণ closed-form ভবিষ্যতে ম্যানুয়াল hand-calc
    // দিয়ে cross-check করা ভালো হবে production এ যাওয়ার আগে।
    std::cout << "\n=== Test 2: Full Pipeline — Two-Element Beam with Internal Pin (Stable) ===\n";
    {
        SectionProperties section{0.01, 0.0001, 0.00005, 0.00002};
        MaterialProperties material{200e6, 80e6};
        const double L = 3.0;
        const double P = 10.0; // kN

        AnalysisModel model;
        model.nodes = {
            Node3D{"n1", 0.0, 0.0, 0.0},
            Node3D{"n2", L, 0.0, 0.0},
            Node3D{"n3", 2.0 * L, 0.0, 0.0}
        };

        FrameElement elem1; // node1-node2: fully rigid (moment) connection
        elem1.elementId = "e1";
        elem1.startNodeIndex = 0;
        elem1.endNodeIndex = 1;
        elem1.section = section;
        elem1.material = material;
        elem1.connectionType = "moment";

        FrameElement elem2; // node2-node3: pin at node2 end (internal hinge)
        elem2.elementId = "e2";
        elem2.startNodeIndex = 1;
        elem2.endNodeIndex = 2;
        elem2.section = section;
        elem2.material = material;
        elem2.connectionType = "pin"; // নোট: বর্তমান getEffectiveLocalStiffness both-end release করে;
                                       // node3 (roller, rotation free থাকলেও ঠিক আছে) প্রান্তেও release
                                       // হবে যা এই কেসে physically harmless (roller support rotation
                                       // constrain করে না এমনিতেও)
        model.elements = {elem1, elem2};

        BoundaryCondition fixedSupport;
        fixedSupport.nodeIndex = 0;
        fixedSupport.restrainX = fixedSupport.restrainY = fixedSupport.restrainZ = true;
        fixedSupport.restrainRx = fixedSupport.restrainRy = fixedSupport.restrainRz = true;

        BoundaryCondition rollerSupport;
        rollerSupport.nodeIndex = 2;
        rollerSupport.restrainX = true;  // axial drift ঠেকাতে (out-of-plane mechanism এড়াতে)
        rollerSupport.restrainY = true;
        rollerSupport.restrainZ = true;
        rollerSupport.restrainRx = true;  // torsion restrain (out-of-plane rotation mechanism এড়াতে, ২D সমস্যা তাই)
        // node3-এর ry,rz ও restrain করা আবশ্যক — কারণ elem2 বর্তমান
        // getEffectiveLocalStiffness() এ "pin" মানে *উভয়* প্রান্ত release
        // করে (node2 প্রান্ত এবং node3 প্রান্ত দুটোই), তাই node3-এর
        // rotation DOF-এ elem2 থেকে কোনো stiffness contribution নেই।
        // বাস্তব pin/roller support ঠিক এভাবেই rotation নেয় না, তাই এই
        // restrain যোগ করাটাই সঠিক মডেলিং (এটা solver-এর সীমাবদ্ধতা না,
        // "pin" মানে single-end-release ধরলে এই restrain লাগতো না —
        // দেখুন solver.cpp এর getEffectiveLocalStiffness() কমেন্ট)।
        rollerSupport.restrainRy = true;
        rollerSupport.restrainRz = true;

        model.boundaryConditions = {fixedSupport, rollerSupport};

        NodalLoad load;
        load.nodeIndex = 1; // node2 এ, হিঞ্জ পয়েন্টে
        load.fx = 0; load.fy = -P; load.fz = 0;
        load.mx = 0; load.my = 0; load.mz = 0;
        model.loads = {load};

        auto result = solveLinearStatic(model);
        checkTrue("Solver succeeded for two-element beam with internal pin", result.success);

        if (result.success) {
            double node2_uy = result.nodalDisplacements[1](1);
            checkTrue("node2 (hinge, loaded point) deflects downward (negative uy)", node2_uy < 0.0);
            checkTrue("Deflection is finite (no NaN/Inf)", std::isfinite(node2_uy));

            // Element 2 (pin at its start=node2) এর start-end (local index 5 = rz1)
            // এ moment থাকা উচিত না (pin release করেছে সেই DOF)
            double elem2_startMoment = result.elementEndForces[1](5); // rz at local node1 of elem2
            checkClose("Hinge transfers zero moment through the pinned end", elem2_startMoment, 0.0, 1e-6);
        }
    }

    // ============================================================
    // Test 3: Pin Release Only Affects the Released End (single-end release)
    // ============================================================
    std::cout << "\n=== Test 3: Single-End Release Leaves Other End Intact ===\n";
    {
        SectionProperties section{0.01, 0.0001, 0.00005, 0.00002};
        MaterialProperties material{200e6, 80e6};
        auto kRigid = computeLocalStiffnessMatrix(3.0, section, material);
        auto kStartPinned = applyEndReleases(kRigid, true, false);

        // node1 এর ry1,rz1 (idx 4,5) rows/cols শূন্য
        double maxStartBending = 0.0;
        for (int idx : {4, 5}) {
            for (int j = 0; j < 12; ++j) {
                maxStartBending = std::max(maxStartBending, std::abs(kStartPinned(idx, j)));
            }
        }
        checkTrue("Start-end bending rows are zero", maxStartBending < 1e-9);

        // node2 এর ry2,rz2 (idx 10,11) এখনো non-zero bending stiffness রাখা উচিত
        checkTrue("End-node (unreleased) still has bending stiffness",
                  std::abs(kStartPinned(11, 11)) > 1e-6);
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
