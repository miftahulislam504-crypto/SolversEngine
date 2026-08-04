#include "stiffness.h"
#include "solver.h"
#include "types.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace civilos;

int testsPassed = 0;
int testsFailed = 0;

void checkClose(const std::string& name, double actual, double expected, double tolerance) {
    double relError = std::abs(expected) > 1e-12
        ? std::abs(actual - expected) / std::abs(expected)
        : std::abs(actual - expected);
    if (relError < tolerance) {
        std::cout << "  [PASS] " << name << " = " << actual << " (expected " << expected
                   << ", rel error = " << relError << ")\n";
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

    const double E = 200e6;   // kN/m2
    const double A = 0.01;    // m2
    const double Ixx = 0.0001; // m4 (strong axis, weak axis একই এই টেস্টে সরলতার জন্য)
    const double Iyy = 0.01;     // ইচ্ছাকৃতভাবে Ixx এর চেয়ে অনেক বড় — out-of-plane
                                  // (weak-axis, local z bending) buckling mode কে suppress
                                  // করতে, যাতে শুধু in-plane (strong-axis, Ixx-নির্ভর)
                                  // buckling mode-ই critical (সবচেয়ে ছোট |λ|) হয়। Iyy=Ixx
                                  // রাখলে একটা degenerate out-of-plane mode critical হয়ে
                                  // যায় (physically সঠিক কিন্তু এই টেস্টের উদ্দেশ্য — classic
                                  // in-plane Euler formula verify করা — এর সাথে অপ্রাসঙ্গিক)।
    const double J = 0.00005;
    const double G = E / (2 * 1.3);

    SectionProperties section{A, Ixx, Iyy, J};
    MaterialProperties material{E, G, 7.85}; // density এখানে ব্যবহৃত হচ্ছে না (buckling mass-independent), কিন্তু struct এ প্রয়োজন

    // ============================================================
    // Test 1: Pin-Pin Column — Euler Buckling Load (K=1)
    // ============================================================
    // Classic formula: Pcr = π² EI / (KL)², K=1 for pin-pin.
    // (Timoshenko, "Theory of Elastic Stability", বা যেকোনো standard
    // structural mechanics textbook — এটা structural engineering এর
    // সবচেয়ে সুপরিচিত সূত্র)
    //
    // মডেল: একটা vertical column, pin support উভয় প্রান্তে (rotation
    // free, translation restrained), unit reference axial load (P=1 kN
    // compression, ঋণাত্মক Y দিকে প্রয়োগ করা যাতে member এ compression
    // তৈরি হয়) প্রয়োগ করে critical load factor (λ) বের করা — তারপর
    // Pcr_actual = λ × P_applied।
    std::cout << "=== Test 1: Pin-Pin Column — Euler Buckling Load (10-element mesh) ===\n";
    {
        const double L = 4.0;
        const int numElements = 10;
        const double dL = L / numElements;
        const double appliedLoad = 1.0; // kN, reference load (compression)

        AnalysisModel model;
        for (int i = 0; i <= numElements; ++i) {
            model.nodes.push_back(Node3D{"n" + std::to_string(i), 0.0, dL * i, 0.0});
        }
        for (int i = 0; i < numElements; ++i) {
            FrameElement elem;
            elem.elementId = "e" + std::to_string(i);
            elem.startNodeIndex = i;
            elem.endNodeIndex = i + 1;
            elem.section = section;
            elem.material = material;
            elem.connectionType = "moment";
            model.elements.push_back(elem);
        }

        // Pin-pin: translation restrained at both ends, bending rotation
        // free (pin). এটা একটা VERTICAL column (local element x-axis =
        // global Y) — তাই torsion (element-এর নিজস্ব axis এর চারপাশে
        // rotation) আসলে global Ry, global Rx না (stiffness.cpp এর
        // computeTransformationMatrix docstring দেখুন, vertical member
        // এর local-axis mapping convention)। তাই torsion আটকাতে
        // restrainRy=true (Rx,Rz bending rotation free রাখা হচ্ছে,
        // pin support এর সংজ্ঞা অনুযায়ী)।
        BoundaryCondition bottom;
        bottom.nodeIndex = 0;
        bottom.restrainX = true; bottom.restrainY = true; bottom.restrainZ = true;
        bottom.restrainRx = false; bottom.restrainRy = true; bottom.restrainRz = false;

        BoundaryCondition top;
        top.nodeIndex = numElements;
        top.restrainX = true; top.restrainY = false; top.restrainZ = true;
        top.restrainRx = false; top.restrainRy = true; top.restrainRz = false;

        model.boundaryConditions = {bottom, top};

        // Compression load at the top node, pointing down (-Y) — compresses the column
        NodalLoad load;
        load.nodeIndex = numElements;
        load.fx = 0; load.fy = -appliedLoad; load.fz = 0;
        load.mx = 0; load.my = 0; load.mz = 0;
        model.loads = {load};

        auto result = solveLinearBuckling(model, 3);
        checkTrue("Buckling solve succeeded", result.success);

        if (result.success) {
            checkTrue("At least 1 mode computed", result.numModesComputed >= 1);

            // সবচেয়ে ছোট |λ| ধনাত্মক হওয়া উচিত (compression load এ
            // buckle করার জন্য positive multiplier লাগে, যেহেতু আমরা
            // ইতিমধ্যে compression প্রয়োগ করেছি — negative critical
            // factor মানে হতো বিপরীত দিকে (tension এ) buckling, যা এই
            // কলামের জন্য প্রযোজ্য না)
            checkTrue("First critical load factor is positive (compression buckling)",
                      result.criticalLoadFactors[0] > 0.0);

            const double K_pinpin = 1.0;
            const double expectedPcr = (M_PI * M_PI * E * Ixx) / std::pow(K_pinpin * L, 2);
            const double actualPcr = result.criticalLoadFactors[0] * appliedLoad;

            checkClose("Critical buckling load (kN) vs. Euler formula",
                       actualPcr, expectedPcr, 0.02); // 2% tolerance — FE mesh discretization
        }
    }

    // ============================================================
    // Test 2: Fixed-Free Cantilever Column — Euler Buckling Load (K=2)
    // ============================================================
    std::cout << "\n=== Test 2: Fixed-Free Cantilever Column — Euler Buckling Load (K=2) ===\n";
    {
        const double L = 3.0;
        const int numElements = 10;
        const double dL = L / numElements;
        const double appliedLoad = 1.0;

        AnalysisModel model;
        for (int i = 0; i <= numElements; ++i) {
            model.nodes.push_back(Node3D{"n" + std::to_string(i), 0.0, dL * i, 0.0});
        }
        for (int i = 0; i < numElements; ++i) {
            FrameElement elem;
            elem.elementId = "e" + std::to_string(i);
            elem.startNodeIndex = i;
            elem.endNodeIndex = i + 1;
            elem.section = section;
            elem.material = material;
            elem.connectionType = "moment";
            model.elements.push_back(elem);
        }

        BoundaryCondition fixedBase;
        fixedBase.nodeIndex = 0;
        fixedBase.restrainX = fixedBase.restrainY = fixedBase.restrainZ = true;
        fixedBase.restrainRx = fixedBase.restrainRy = fixedBase.restrainRz = true;
        model.boundaryConditions = {fixedBase};

        NodalLoad load;
        load.nodeIndex = numElements;
        load.fx = 0; load.fy = -appliedLoad; load.fz = 0;
        model.loads = {load};

        auto result = solveLinearBuckling(model, 3);
        checkTrue("Buckling solve succeeded", result.success);

        if (result.success) {
            const double K_fixedFree = 2.0;
            const double expectedPcr = (M_PI * M_PI * E * Ixx) / std::pow(K_fixedFree * L, 2);
            const double actualPcr = result.criticalLoadFactors[0] * appliedLoad;

            checkClose("Critical buckling load (kN) vs. Euler formula (K=2)",
                       actualPcr, expectedPcr, 0.02);
        }
    }

    // ============================================================
    // Test 3: Fixed-Fixed Column — Euler Buckling Load (K=0.5)
    // ============================================================
    std::cout << "\n=== Test 3: Fixed-Fixed Column — Euler Buckling Load (K=0.5) ===\n";
    {
        const double L = 4.0;
        const int numElements = 10;
        const double dL = L / numElements;
        const double appliedLoad = 1.0;

        AnalysisModel model;
        for (int i = 0; i <= numElements; ++i) {
            model.nodes.push_back(Node3D{"n" + std::to_string(i), 0.0, dL * i, 0.0});
        }
        for (int i = 0; i < numElements; ++i) {
            FrameElement elem;
            elem.elementId = "e" + std::to_string(i);
            elem.startNodeIndex = i;
            elem.endNodeIndex = i + 1;
            elem.section = section;
            elem.material = material;
            elem.connectionType = "moment";
            model.elements.push_back(elem);
        }

        BoundaryCondition bottom;
        bottom.nodeIndex = 0;
        bottom.restrainX = bottom.restrainY = bottom.restrainZ = true;
        bottom.restrainRx = bottom.restrainRy = bottom.restrainRz = true;

        // Top: translation restrained in X,Z (prevent sway — keep it a
        // pure fixed-fixed column, not a sway frame), free in Y (axial,
        // load applied there), rotation fully fixed
        BoundaryCondition top;
        top.nodeIndex = numElements;
        top.restrainX = true; top.restrainY = false; top.restrainZ = true;
        top.restrainRx = top.restrainRy = top.restrainRz = true;

        model.boundaryConditions = {bottom, top};

        NodalLoad load;
        load.nodeIndex = numElements;
        load.fx = 0; load.fy = -appliedLoad; load.fz = 0;
        model.loads = {load};

        auto result = solveLinearBuckling(model, 3);
        checkTrue("Buckling solve succeeded", result.success);

        if (result.success) {
            const double K_fixedFixed = 0.5;
            const double expectedPcr = (M_PI * M_PI * E * Ixx) / std::pow(K_fixedFixed * L, 2);
            const double actualPcr = result.criticalLoadFactors[0] * appliedLoad;

            checkClose("Critical buckling load (kN) vs. Euler formula (K=0.5)",
                       actualPcr, expectedPcr, 0.03); // মেশ discretization একটু বেশি এখানে (higher mode shape curvature)
        }
    }

    // ============================================================
    // Test 4: No Load Present → Graceful Failure (Not a Crash)
    // ============================================================
    std::cout << "\n=== Test 4: No Load Present — Graceful Failure ===\n";
    {
        AnalysisModel model;
        model.nodes = {
            Node3D{"n0", 0.0, 0.0, 0.0},
            Node3D{"n1", 0.0, 3.0, 0.0},
        };
        FrameElement elem;
        elem.elementId = "e0";
        elem.startNodeIndex = 0;
        elem.endNodeIndex = 1;
        elem.section = section;
        elem.material = material;
        elem.connectionType = "moment";
        model.elements = {elem};

        BoundaryCondition fixedBase;
        fixedBase.nodeIndex = 0;
        fixedBase.restrainX = fixedBase.restrainY = fixedBase.restrainZ = true;
        fixedBase.restrainRx = fixedBase.restrainRy = fixedBase.restrainRz = true;
        model.boundaryConditions = {fixedBase};
        // model.loads left empty deliberately

        auto result = solveLinearBuckling(model, 3);
        checkTrue("Solver fails gracefully (no crash) when no load is present", !result.success);
        if (!result.success) {
            checkTrue("Error message is non-empty", !result.errorMessage.empty());
        }
    }

    // ============================================================
    // Test 5: Geometric stiffness matrix is symmetric (mathematical sanity)
    // ============================================================
    std::cout << "\n=== Test 5: Geometric Stiffness Matrix Symmetry ===\n";
    {
        auto Kg = computeLocalGeometricStiffnessMatrix(3.0, -50.0); // compression
        double maxAsymmetry = (Kg - Kg.transpose()).cwiseAbs().maxCoeff();
        checkTrue("Geometric stiffness matrix is symmetric", maxAsymmetry < 1e-9);

        // axial(0,6) ও torsion(3,9) DOF এ কোনো geometric stiffness contribution
        // থাকা উচিত না (stiffness.h এ documented সরলীকরণ)
        checkTrue("No geometric stiffness on axial DOF",
                  std::abs(Kg(0, 0)) < 1e-12 && std::abs(Kg(6, 6)) < 1e-12);
        checkTrue("No geometric stiffness on torsion DOF",
                  std::abs(Kg(3, 3)) < 1e-12 && std::abs(Kg(9, 9)) < 1e-12);
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
