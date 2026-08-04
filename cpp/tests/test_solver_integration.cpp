#include "solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace civilos;

int testsPassed = 0;
int testsFailed = 0;

void checkClose(const std::string& name, double actual, double expected, double tolerance = 1e-4) {
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
    // Test: Simply Supported Beam, Center Point Load
    // ============================================================
    // ক্লাসিক্যাল সমাধান (Hibbeler, Timoshenko টেক্সটবুকে প্রতিষ্ঠিত):
    // একটা সিম্পলি-সাপোর্টেড বিম, span L, center এ point load P হলে
    // মধ্যবিন্দুতে max deflection δ = PL³/(48EI), এবং প্রতিটা support
    // এ reaction = P/2 (symmetry থেকে)।
    //
    // এই টেস্ট single-element cantilever টেস্টের (test_stiffness.cpp)
    // থেকে গুরুত্বপূর্ণভাবে ভিন্ন — এটা দুইটা element (তিনটা node:
    // left support, center load point, right support) ব্যবহার করে
    // global assembly পুরোপুরি পরীক্ষা করে, single-element static
    // condensation না।
    std::cout << "=== Integration Test: Simply Supported Beam (2 elements) ===\n";
    {
        const double L = 6.0;      // মোট span, মিটার (প্রতিটা element 3m)
        const double E = 200e6;    // kN/m² (steel)
        const double I = 8.0e-5;   // m⁴
        const double A = 0.005;    // m²
        const double G = 77e6;     // kN/m²
        const double J = 5e-7;     // m⁴
        const double P = 20.0;     // kN, center এ প্রয়োগ

        AnalysisModel model;

        // তিনটা node: 0 (left support), 1 (center), 2 (right support)
        // — সব X-অক্ষ বরাবর, Y=0, Z=0 (horizontal beam)
        model.nodes.push_back({"n0", 0.0, 0.0, 0.0});
        model.nodes.push_back({"n1", L / 2.0, 0.0, 0.0});
        model.nodes.push_back({"n2", L, 0.0, 0.0});

        SectionProperties section{A, I, I, J};
        MaterialProperties material{E, G};

        model.elements.push_back({"e0", 0, 1, section, material, "moment"});
        model.elements.push_back({"e1", 1, 2, section, material, "moment"});

        // Boundary conditions — simply supported মানে উভয় প্রান্তে
        // translation restrained (X,Y,Z সব দিকে, 3D stability এর জন্য
        // — যদিও লোড শুধু vertical plane এ, বাকি DOF ও restrain করা
        // দরকার নাহলে out-of-plane mechanism তৈরি হবে) কিন্তু rotation
        // মুক্ত (pin support, moment নেয় না)।
        BoundaryCondition leftSupport;
        leftSupport.nodeIndex = 0;
        leftSupport.restrainX = true;
        leftSupport.restrainY = true;
        leftSupport.restrainZ = true;
        leftSupport.restrainRx = true;  // torsion restrained (out-of-plane stability)
        leftSupport.restrainRy = false;
        leftSupport.restrainRz = false; // in-plane bending rotation মুক্ত — এটাই "pin" support কে বোঝায়
        model.boundaryConditions.push_back(leftSupport);

        BoundaryCondition rightSupport;
        rightSupport.nodeIndex = 2;
        rightSupport.restrainX = false; // একটা প্রান্ত axially মুক্ত থাকা উচিত (roller support), নাহলে thermal/axial লোডে over-constrained হয়ে যেত
        rightSupport.restrainY = true;
        rightSupport.restrainZ = true;
        rightSupport.restrainRx = true;
        rightSupport.restrainRy = false;
        rightSupport.restrainRz = false;
        model.boundaryConditions.push_back(rightSupport);

        // Center node এ শুধু vertical DOF ছাড়া বাকি সব মুক্ত থাকতে
        // পারতো, কিন্তু out-of-plane stability এর জন্য torsion/lateral
        // এখানেও restrain করা দরকার এই সরল 2D-in-3D টেস্ট কেসে (বাস্তব
        // building model এ slab/diaphragm এই কাজ করে, এই isolated
        // beam টেস্টে সেটা নেই তাই ম্যানুয়ালি restrain করা হচ্ছে)।
        BoundaryCondition centerLateral;
        centerLateral.nodeIndex = 1;
        centerLateral.restrainX = false;
        centerLateral.restrainY = false; // এখানেই load প্রয়োগ হবে, তাই Y মুক্ত রাখতে হবে
        centerLateral.restrainZ = true;  // out-of-plane (Z দিকে) নড়াচড়া বন্ধ
        centerLateral.restrainRx = true;
        centerLateral.restrainRy = false;
        centerLateral.restrainRz = false;
        model.boundaryConditions.push_back(centerLateral);

        NodalLoad centerLoad;
        centerLoad.nodeIndex = 1;
        centerLoad.fx = 0;
        centerLoad.fy = -P; // negative Y = নিচের দিকে (gravity direction)
        centerLoad.fz = 0;
        centerLoad.mx = 0;
        centerLoad.my = 0;
        centerLoad.mz = 0;
        model.loads.push_back(centerLoad);

        AnalysisResult result = solveLinearStatic(model);

        checkTrue("Solver succeeded", result.success);

        if (result.success) {
            const double centerDeflection = result.nodalDisplacements[1](1); // node 1, DOF index 1 = Y
            const double expectedDeflection = -(P * L * L * L) / (48.0 * E * I); // negative কারণ load negative-Y দিকে

            checkClose("Center deflection δ = PL³/(48EI)", centerDeflection, expectedDeflection);

            // Reaction force হিসাব: element end force থেকে left support
            // এ vertical reaction বের করা যায় (element 0 এর start-node
            // end force, local-y direction, কারণ element horizontal
            // তাই local-y = global-y এই কেসে)।
            const double leftReactionFromElement = result.elementEndForces[0](1); // element 0, local DOF 1 = v1 (shear at start)
            const double expectedReaction = P / 2.0;

            checkClose("Left support reaction ≈ P/2 (from element shear)",
                       std::abs(leftReactionFromElement), expectedReaction, 1e-3);

            // Phase 10n — result.reactionForces (penalty-method recovery)
            // কে দুইভাবে verify করা হচ্ছে: (1) উপরের element-shear-derived
            // reaction এর সাথে cross-check (দুইটা সম্পূর্ণ independent
            // পদ্ধতি একই সংখ্যা দিলে সেটা শক্তিশালী প্রমাণ যে নতুন
            // reactionForces সঠিক), (2) classical hand-calculation P/2
            // এর সাথে সরাসরি। boundaryConditions এ push হওয়া ক্রম অনুযায়ী
            // index 0 = leftSupport, index 1 = rightSupport, index 2 =
            // centerLateral (দেখুন উপরে push_back ক্রম)।
            checkTrue("reactionForces has one entry per boundary condition",
                      result.reactionForces.size() == model.boundaryConditions.size());

            if (result.reactionForces.size() == model.boundaryConditions.size()) {
                const double leftReactionY = result.reactionForces[0](1); // leftSupport, global Y
                const double rightReactionY = result.reactionForces[1](1); // rightSupport, global Y

                checkClose("reactionForces[left].Y = +P/2 (upward, penalty-method)",
                           leftReactionY, expectedReaction, 1e-3);
                checkClose("reactionForces[right].Y = +P/2 (upward, penalty-method)",
                           rightReactionY, expectedReaction, 1e-3);
                checkClose("reactionForces[left].Y matches element-shear-derived reaction",
                           leftReactionY, leftReactionFromElement, 1e-3);
                // Global equilibrium: সব vertical reaction এর যোগফল =
                // applied load এর বিপরীত দিকে ও সমান মান (ΣF=0 মানে
                // reaction = −applied — load টা −P (নিচের দিকে) ছিল,
                // তাই reaction যোগফল +P হওয়ার কথা upward দিকে, যেটা এই
                // sign convention এ +P; নিচের check সেটাই যাচাই করে)।
                checkClose("Sum of vertical reactions balances applied load (global equilibrium)",
                           leftReactionY + rightReactionY, P, 1e-3);
            }
        }
    }

    // ============================================================
    // Test: Unstable Structure Detection
    // ============================================================
    // কোনো boundary condition ছাড়া একটা structure সমাধান করার চেষ্টা
    // করলে solver স্পষ্টভাবে ব্যর্থ হওয়া উচিত, কোনো ভুল/অর্থহীন সংখ্যা
    // silently না দিয়ে — এটা একটা নিরাপত্তা বৈশিষ্ট্য, ভুল মডেল থেকে
    // বিপজ্জনক ডিজাইন সিদ্ধান্ত এড়াতে।
    std::cout << "\n=== Test: Unstable Structure (No Supports) Detection ===\n";
    {
        AnalysisModel model;
        model.nodes.push_back({"n0", 0.0, 0.0, 0.0});
        model.nodes.push_back({"n1", 3.0, 0.0, 0.0});

        SectionProperties section{0.005, 8e-5, 8e-5, 5e-7};
        MaterialProperties material{200e6, 77e6};
        model.elements.push_back({"e0", 0, 1, section, material, "moment"});
        // কোনো boundaryConditions যোগ করা হচ্ছে না — ইচ্ছাকৃতভাবে অস্থির

        AnalysisResult result = solveLinearStatic(model);
        checkTrue("Solver correctly reports failure for unrestrained structure", !result.success);
        std::cout << "    (error message: " << result.errorMessage << ")\n";
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
