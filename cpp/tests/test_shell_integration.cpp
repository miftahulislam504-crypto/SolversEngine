#include "shell.h"
#include "solver.h"
#include "types.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace civilos;

int testsPassed = 0;
int testsFailed = 0;

void checkTrue(const std::string& name, bool condition) {
    if (condition) {
        std::cout << "  [PASS] " << name << "\n";
        testsPassed++;
    } else {
        std::cout << "  [FAIL] " << name << "\n";
        testsFailed++;
    }
}

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

int main() {
    std::cout << std::setprecision(8);

    const double E = 200e6, nu = 0.2, t = 0.2;

    // ============================================================
    // Test 1: Cantilever Shear Wall — Full Pipeline via solveLinearStatic
    // ============================================================
    // একটা 4x4-element mesh দিয়ে একটা cantilever wall (base fixed,
    // lateral shear load top-এ) — যাচাই করে shellElements সঠিকভাবে
    // AnalysisModel এর মধ্য দিয়ে assembleGlobalStiffness() এ পৌঁছায়
    // ও solveLinearStatic() পূর্ণ pipeline (boundary condition, load
    // vector, sparse solve) দিয়ে সঠিক ফলাফল দেয়।
    std::cout << "=== Test 1: Cantilever Shear Wall — Full solveLinearStatic Pipeline ===\n";
    {
        const double a = 4.0, h = 4.0;
        const int n = 4;
        const double dLx = a / n, dLy = h / n;

        AnalysisModel model;
        auto nodeIndex = [&](int i, int j) { return j * (n + 1) + i; };
        for (int j = 0; j <= n; ++j)
            for (int i = 0; i <= n; ++i)
                model.nodes.push_back(Node3D{"n", i * dLx, j * dLy, 0.0});

        MaterialProperties material{E, E / (2 * (1 + nu)), 2.4, nu};
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                ShellElement shell;
                shell.elementId = "s" + std::to_string(j * n + i);
                shell.nodeIndices = {nodeIndex(i, j), nodeIndex(i + 1, j),
                                      nodeIndex(i + 1, j + 1), nodeIndex(i, j + 1)};
                shell.thickness = t;
                shell.material = material;
                model.shellElements.push_back(shell);
            }
        }

        for (int i = 0; i <= n; ++i) {
            BoundaryCondition bc;
            bc.nodeIndex = nodeIndex(i, 0);
            bc.restrainX = bc.restrainY = bc.restrainZ = true;
            bc.restrainRx = bc.restrainRy = bc.restrainRz = true;
            model.boundaryConditions.push_back(bc);
        }
        // Out-of-plane DOF restrained everywhere (pure in-plane wall test)
        for (size_t i = 0; i < model.nodes.size(); ++i) {
            BoundaryCondition bc;
            bc.nodeIndex = static_cast<int>(i);
            bc.restrainX = bc.restrainY = false;
            bc.restrainZ = true;
            bc.restrainRx = bc.restrainRy = true;
            bc.restrainRz = false;
            model.boundaryConditions.push_back(bc);
        }

        const double totalLoad = 10.0; // kN
        for (int i = 0; i <= n; ++i) {
            NodalLoad load;
            load.nodeIndex = nodeIndex(i, n);
            load.fx = totalLoad / (n + 1);
            load.fy = 0; load.fz = 0; load.mx = 0; load.my = 0; load.mz = 0;
            model.loads.push_back(load);
        }

        auto result = solveLinearStatic(model);
        checkTrue("Shell-only model solves successfully via solveLinearStatic", result.success);

        if (result.success) {
            const int topCenterNode = nodeIndex(n / 2, n);
            const int baseCenterNode = nodeIndex(n / 2, 0);

            checkTrue("Top node deflects laterally (ux > 0, load direction)",
                      result.nodalDisplacements[topCenterNode](0) > 0.0);
            checkTrue("Base node stays fixed (ux ≈ 0)",
                      std::abs(result.nodalDisplacements[baseCenterNode](0)) < 1e-9);
            checkTrue("Top displacement is finite (no NaN/Inf)",
                      std::isfinite(result.nodalDisplacements[topCenterNode](0)));

            // Deflection should increase monotonically with height (basic sanity
            // for a cantilever under lateral tip load)
            double midHeightDisp = result.nodalDisplacements[nodeIndex(n / 2, n / 2)](0);
            double topDisp = result.nodalDisplacements[topCenterNode](0);
            checkTrue("Top deflection exceeds mid-height deflection (monotonic cantilever behavior)",
                      topDisp > midHeightDisp);
        }
    }

    // ============================================================
    // Test 2: Mixed Frame + Shell Model — Both Element Types Coexist
    // ============================================================
    // একটা মডেল যেখানে frame element (column) ও shell element (wall
    // panel) একই node graph শেয়ার করে — যাচাই করে assembly উভয়
    // ধরনের element একসাথে সঠিকভাবে সমাধান করতে পারে (একটা shared
    // node এর stiffness contribution উভয় element থেকেই আসা উচিত)।
    std::cout << "\n=== Test 2: Mixed Frame + Shell Model ===\n";
    {
        SectionProperties section{0.01, 0.0001, 0.0001, 0.00005};
        MaterialProperties frameMaterial{E, E / (2 * 1.3), 7.85, 0.3};
        MaterialProperties shellMaterial{E, E / (2 * (1 + nu)), 2.4, nu};

        AnalysisModel model;
        model.nodes = {
            Node3D{"n0", 0.0, 0.0, 0.0},
            Node3D{"n1", 1.0, 0.0, 0.0},
            Node3D{"n2", 1.0, 1.0, 0.0},
            Node3D{"n3", 0.0, 1.0, 0.0},
            Node3D{"n4", 0.0, 2.0, 0.0}, // extra node for a frame element on top of the shell
        };

        ShellElement shell;
        shell.elementId = "wall1";
        shell.nodeIndices = {0, 1, 2, 3};
        shell.thickness = t;
        shell.material = shellMaterial;
        model.shellElements = {shell};

        FrameElement column;
        column.elementId = "col1";
        column.startNodeIndex = 3; // top-left of the shell
        column.endNodeIndex = 4;   // extends upward
        column.section = section;
        column.material = frameMaterial;
        column.connectionType = "moment";
        model.elements = {column};

        // Fix the base of the shell (nodes 0,1)
        for (int idx : {0, 1}) {
            BoundaryCondition bc;
            bc.nodeIndex = idx;
            bc.restrainX = bc.restrainY = bc.restrainZ = true;
            bc.restrainRx = bc.restrainRy = bc.restrainRz = true;
            model.boundaryConditions.push_back(bc);
        }
        // Restrain out-of-plane DOF everywhere (in-plane test only)
        for (int idx = 0; idx < 5; ++idx) {
            BoundaryCondition bc;
            bc.nodeIndex = idx;
            bc.restrainX = bc.restrainY = false;
            bc.restrainZ = true;
            bc.restrainRx = bc.restrainRy = true;
            bc.restrainRz = false;
            model.boundaryConditions.push_back(bc);
        }

        NodalLoad load;
        load.nodeIndex = 4; // tip of the column
        load.fx = 1.0; load.fy = 0; load.fz = 0;
        model.loads = {load};

        auto result = solveLinearStatic(model);
        checkTrue("Mixed frame+shell model solves successfully", result.success);

        if (result.success) {
            checkTrue("Column tip (node 4) deflects under load",
                      std::abs(result.nodalDisplacements[4](0)) > 1e-9);
            checkTrue("All displacements finite (no NaN/Inf anywhere)", [&]() {
                for (const auto& d : result.nodalDisplacements) {
                    for (int i = 0; i < 6; ++i) {
                        if (!std::isfinite(d(i))) return false;
                    }
                }
                return true;
            }());
        }
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
