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

    const double E = 200e6;      // kN/m2 (steel-like)
    const double rho = 7.85;     // tonne/m3 (steel density, 7850 kg/m3 = 7.85 tonne/m3)
    const double A = 0.01;       // m2
    const double Ixx = 0.0001;   // m4 (strong axis)
    const double Iyy = 0.0001;   // m4 (symmetric section, weak axis = strong axis এখানে সরলতার জন্য)
    const double J = 0.00005;    // m4
    const double G = E / (2 * 1.3); // poisson ~0.3 হলে কাছাকাছি, exact মান গুরুত্বপূর্ণ না এই টেস্টে

    SectionProperties section{A, Ixx, Iyy, J};
    MaterialProperties material{E, G, rho};

    const double m = rho * A; // mass per unit length, tonne/m

    // ============================================================
    // Test 1: Fine-mesh cantilever beam — first natural frequency
    // ============================================================
    // Classic formula: ω1 = (β1*L)^2 * sqrt(EI / (m*L^4)), β1*L = 1.875104...
    // (Blevins, "Formulas for Natural Frequency and Mode Shape", Table 8-1)
    //
    // একটা একক 12-DOF element দিয়ে cantilever মডেল করলে coarse-mesh
    // discretization error আসবে (consistent mass FE approximation
    // সবসময় continuous-beam exact solution থেকে কিছুটা ভিন্ন হয়,
    // বিশেষত একক element এ) — তাই একাধিক element দিয়ে (mesh refinement)
    // মডেল করে সেই error কমানো হচ্ছে, যাতে ভাল convergence tolerance এ
    // (few %) মিলে।
    std::cout << "=== Test 1: Cantilever Beam — First Natural Frequency (10-element mesh) ===\n";
    {
        const double L = 4.0; // total length, m
        const int numElements = 10;
        const double dL = L / numElements;

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

        auto result = solveModalAnalysis(model, 3);
        checkTrue("Modal solve succeeded", result.success);

        if (result.success) {
            checkTrue("At least 1 mode computed", result.numModesComputed >= 1);

            const double betaL = 1.8751040687;
            const double expectedOmega1 = (betaL * betaL) * std::sqrt((E * Ixx) / (m * L * L * L * L));
            const double expectedFreqHz1 = expectedOmega1 / (2.0 * M_PI);

            checkClose("First natural frequency (Hz) vs. Blevins formula",
                       result.naturalFrequenciesHz[0], expectedFreqHz1, 0.02); // 2% tolerance — FE discretization error

            checkTrue("Frequencies are in ascending order",
                      result.naturalFrequenciesHz[0] <= result.naturalFrequenciesHz[1] &&
                      result.naturalFrequenciesHz[1] <= result.naturalFrequenciesHz[2]);

            checkTrue("All frequencies are positive and finite",
                      std::isfinite(result.naturalFrequenciesHz[0]) && result.naturalFrequenciesHz[0] > 0.0);
        }
    }

    // ============================================================
    // Test 2: Simply-supported beam — first natural frequency
    // ============================================================
    // Classic formula: ω1 = π^2 * sqrt(EI / (m*L^4))  (Blevins Table 8-1)
    std::cout << "\n=== Test 2: Simply-Supported Beam — First Natural Frequency (10-element mesh) ===\n";
    {
        const double L = 6.0;
        const int numElements = 10;
        const double dL = L / numElements;

        AnalysisModel model;
        for (int i = 0; i <= numElements; ++i) {
            model.nodes.push_back(Node3D{"n" + std::to_string(i), dL * i, 0.0, 0.0});
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

        // Pin-pin support: translation restrained (Y,Z), rotation free (X-axis
        // beam, so bending happens in the X-Y plane primarily — restrain Rx
        // (torsion) at both ends and X-translation at one end to prevent rigid
        // body modes, but leave Ry/Rz free for the classic simply-supported
        // boundary condition)
        BoundaryCondition left;
        left.nodeIndex = 0;
        left.restrainX = true; left.restrainY = true; left.restrainZ = true;
        left.restrainRx = true; left.restrainRy = false; left.restrainRz = false;

        BoundaryCondition right;
        right.nodeIndex = numElements;
        right.restrainX = false; right.restrainY = true; right.restrainZ = true;
        right.restrainRx = true; right.restrainRy = false; right.restrainRz = false;

        model.boundaryConditions = {left, right};

        auto result = solveModalAnalysis(model, 3);
        checkTrue("Modal solve succeeded", result.success);

        if (result.success) {
            const double expectedOmega1 = (M_PI * M_PI) * std::sqrt((E * Ixx) / (m * L * L * L * L));
            const double expectedFreqHz1 = expectedOmega1 / (2.0 * M_PI);

            checkClose("First natural frequency (Hz) vs. Blevins formula",
                       result.naturalFrequenciesHz[0], expectedFreqHz1, 0.02);
        }
    }

    // ============================================================
    // Test 3: Mode shapes are M-orthonormal (φᵢᵀMφⱼ = δᵢⱼ)
    // ============================================================
    // Modal analysis এর একটা মৌলিক গাণিতিক ধর্ম — ভিন্ন mode শুধু
    // orthogonal না, M-normalized ও (mass-normalization convention,
    // types.h এর ModalAnalysisResult docstring এ ব্যাখ্যা করা)। এটা
    // যাচাই করে confirm করা যায় solveModalAnalysis() সঠিকভাবে reduced
    // system থেকে full-DOF mode shape এ map করছে এবং eigen-solver এর
    // normalization convention আমাদের ধারণার সাথে মিলছে।
    std::cout << "\n=== Test 3: Mode Shape M-Orthonormality Check ===\n";
    {
        const double L = 3.0;
        AnalysisModel model;
        model.nodes = {
            Node3D{"n0", 0.0, 0.0, 0.0},
            Node3D{"n1", L, 0.0, 0.0},
            Node3D{"n2", 2 * L, 0.0, 0.0},
        };
        for (int i = 0; i < 2; ++i) {
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

        auto result = solveModalAnalysis(model, 4);
        checkTrue("Modal solve succeeded", result.success);

        if (result.success && result.numModesComputed >= 2) {
            // Reconstruct global mass matrix ও mode shape vector মিলিয়ে
            // φᵀMφ যাচাই করা (mode 0 এর সাথে নিজের — normalized হলে ≈1)
            Eigen::SparseMatrix<double> M = assembleGlobalMass(model);
            Eigen::MatrixXd Mdense = Eigen::MatrixXd(M);

            const int totalDOF = static_cast<int>(model.nodes.size()) * 6;

            auto flattenMode = [&](int modeIdx) {
                Eigen::VectorXd v(totalDOF);
                for (size_t i = 0; i < model.nodes.size(); ++i) {
                    v.segment(static_cast<int>(i) * 6, 6) = result.modeShapes[modeIdx][i];
                }
                return v;
            };

            Eigen::VectorXd phi0 = flattenMode(0);
            Eigen::VectorXd phi1 = flattenMode(1);

            double selfProduct0 = phi0.transpose() * Mdense * phi0;
            double crossProduct = phi0.transpose() * Mdense * phi1;

            checkClose("Mode 0 is M-normalized (φ0^T M φ0 ≈ 1)", selfProduct0, 1.0, 1e-6);
            checkTrue("Mode 0 and Mode 1 are M-orthogonal (φ0^T M φ1 ≈ 0)",
                      std::abs(crossProduct) < 1e-6);
        }
    }

    // ============================================================
    // Test 4: Zero/negative density is rejected gracefully (no crash, no NaN)
    // ============================================================
    std::cout << "\n=== Test 4: Zero Density Handled Without Crash ===\n";
    {
        MaterialProperties zeroMassMaterial{E, G, 0.0};
        AnalysisModel model;
        model.nodes = {
            Node3D{"n0", 0.0, 0.0, 0.0},
            Node3D{"n1", 3.0, 0.0, 0.0},
        };
        FrameElement elem;
        elem.elementId = "e0";
        elem.startNodeIndex = 0;
        elem.endNodeIndex = 1;
        elem.section = section;
        elem.material = zeroMassMaterial;
        elem.connectionType = "moment";
        model.elements = {elem};

        BoundaryCondition fixedBase;
        fixedBase.nodeIndex = 0;
        fixedBase.restrainX = fixedBase.restrainY = fixedBase.restrainZ = true;
        fixedBase.restrainRx = fixedBase.restrainRy = fixedBase.restrainRz = true;
        model.boundaryConditions = {fixedBase};

        auto result = solveModalAnalysis(model, 3);
        // zero mass → M matrix singular (positive-semidefinite না,
        // positive-definite না) → solver ব্যর্থ হওয়া উচিত, crash/NaN না
        checkTrue("Solver fails gracefully (no crash) for zero-density material", !result.success);
        if (!result.success) {
            checkTrue("Error message is non-empty", !result.errorMessage.empty());
        }
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
