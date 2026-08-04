#include "shell.h"
#include "types.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

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

    const double E = 200e6; // kN/m2
    const double nu = 0.2;
    const double t = 0.2;   // m
    const double L = 1.0;   // m

    // ============================================================
    // Test 1: Single Shell Element — Matrix Properties (Symmetry,
    // Rigid-Body Modes)
    // ============================================================
    std::cout << "=== Test 1: Single Shell Element — Matrix Symmetry & Rigid-Body Modes ===\n";
    {
        std::array<Node3D, 4> corners = {
            Node3D{"n0", 0.0, 0.0, 0.0},
            Node3D{"n1", L,   0.0, 0.0},
            Node3D{"n2", L,   L,   0.0},
            Node3D{"n3", 0.0, L,   0.0},
        };
        auto K = computeShellGlobalStiffness(corners, t, E, nu);

        double maxAsym = (K - K.transpose()).cwiseAbs().maxCoeff();
        checkTrue("24x24 shell stiffness matrix is symmetric", maxAsym < 1e-6);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 24, 24>> es(K);
        int negativeCount = 0;
        for (int i = 0; i < 24; ++i) {
            if (es.eigenvalues()(i) < -1e-6) negativeCount++;
        }
        checkTrue("No significantly negative eigenvalues (positive-semidefinite)", negativeCount == 0);

        int nearZeroCount = 0;
        for (int i = 0; i < 24; ++i) {
            if (std::abs(es.eigenvalues()(i)) < 1e-6) nearZeroCount++;
        }
        checkTrue("At least 6 near-zero eigenvalues (rigid-body modes present)", nearZeroCount >= 6);
    }

    // ============================================================
    // Test 2: Membrane — Uniaxial Tension vs. 1D Bar Approximation
    // ============================================================
    // একটা square panel এ symmetric nodal tension প্রয়োগ করে stretch
    // measure করা, 1D bar theory (delta = PL/EA) এর সাথে তুলনা। এটা
    // একটা approximation (discrete nodal load, uniform-stress না) —
    // তাই মাঝারি tolerance (~10%) ব্যবহার করা হয়েছে, exact match
    // প্রত্যাশিত না এই simplified loading এ।
    std::cout << "\n=== Test 2: Membrane — Uniaxial Tension vs. 1D Bar Approximation ===\n";
    {
        std::array<Node3D, 4> corners = {
            Node3D{"n0", 0.0, 0.0, 0.0},
            Node3D{"n1", L,   0.0, 0.0},
            Node3D{"n2", L,   L,   0.0},
            Node3D{"n3", 0.0, L,   0.0},
        };
        auto K = computeShellGlobalStiffness(corners, t, E, nu);

        std::vector<int> freeDOFs = {1 * 6 + 0, 2 * 6 + 0}; // ux1, ux2
        Eigen::Matrix2d Kreduced;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                Kreduced(i, j) = K(freeDOFs[i], freeDOFs[j]);

        const double P = 100.0; // kN
        Eigen::Vector2d F(P / 2.0, P / 2.0);
        Eigen::Vector2d u = Kreduced.colPivHouseholderQr().solve(F);

        const double A = t * L;
        const double expectedDelta = P * L / (E * A);

        checkClose("Uniaxial stretch (ux) vs. 1D bar theory (P*L/EA)", u(0), expectedDelta, 0.10);
        checkClose("Symmetric loading gives equal stretch at both nodes", u(1), u(0), 1e-9);
    }

    // ============================================================
    // Test 3: Plate Bending — Simply-Supported Square Plate vs. Timoshenko
    // ============================================================
    // Classic formula: w_center = alpha * q * a^4 / D, alpha=0.00406
    // for nu=0.3, simply-supported square plate (Timoshenko & Woinowsky-
    // Krieger, "Theory of Plates and Shells", Table).
    std::cout << "\n=== Test 3: Plate Bending — Simply-Supported Square Plate vs. Timoshenko ===\n";
    {
        const double nu_plate = 0.3;
        const double t_plate = 0.05; // thin plate: a/t = 80
        const double a = 4.0;        // m
        const double q = 1.0;        // kN/m2

        const int n = 8;
        const double dL = a / n;

        std::vector<Node3D> nodes;
        auto nodeIndex = [&](int i, int j) { return j * (n + 1) + i; };
        for (int j = 0; j <= n; ++j) {
            for (int i = 0; i <= n; ++i) {
                nodes.push_back(Node3D{"n", i * dL, j * dL, 0.0});
            }
        }

        const int totalDOF = static_cast<int>(nodes.size()) * 6;
        std::vector<Eigen::Triplet<double>> triplets;

        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                std::array<Node3D, 4> corners = {
                    nodes[nodeIndex(i, j)], nodes[nodeIndex(i + 1, j)],
                    nodes[nodeIndex(i + 1, j + 1)], nodes[nodeIndex(i, j + 1)],
                };
                auto Ke = computeShellGlobalStiffness(corners, t_plate, E, nu_plate);
                std::array<int, 4> gIdx = {nodeIndex(i, j), nodeIndex(i + 1, j),
                                            nodeIndex(i + 1, j + 1), nodeIndex(i, j + 1)};
                std::array<int, 24> dofIdx{};
                for (int node = 0; node < 4; ++node)
                    for (int d = 0; d < 6; ++d)
                        dofIdx[6 * node + d] = gIdx[node] * 6 + d;
                for (int r = 0; r < 24; ++r)
                    for (int c = 0; c < 24; ++c)
                        if (std::abs(Ke(r, c)) > 1e-15)
                            triplets.emplace_back(dofIdx[r], dofIdx[c], Ke(r, c));
            }
        }

        Eigen::SparseMatrix<double> K(totalDOF, totalDOF);
        K.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::VectorXd F = Eigen::VectorXd::Zero(totalDOF);
        const double elemArea = dL * dL;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                std::array<int, 4> gIdx = {nodeIndex(i, j), nodeIndex(i + 1, j),
                                            nodeIndex(i + 1, j + 1), nodeIndex(i, j + 1)};
                for (int node = 0; node < 4; ++node) {
                    F(gIdx[node] * 6 + 2) -= q * elemArea / 4.0;
                }
            }
        }

        const double PENALTY = 1e12;
        for (int j = 0; j <= n; ++j) {
            for (int i = 0; i <= n; ++i) {
                bool onEdge = (i == 0 || i == n || j == 0 || j == n);
                if (onEdge) {
                    int idx = nodeIndex(i, j) * 6 + 2;
                    K.coeffRef(idx, idx) += PENALTY;
                }
            }
        }
        for (size_t nIdx = 0; nIdx < nodes.size(); ++nIdx) {
            for (int d : {0, 1, 5}) { // in-plane translations + drilling restrained (pure-bending test simplification)
                int idx = static_cast<int>(nIdx) * 6 + d;
                K.coeffRef(idx, idx) += PENALTY;
            }
        }

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(K);
        checkTrue("Plate mesh stiffness matrix decomposition succeeded", solver.info() == Eigen::Success);

        if (solver.info() == Eigen::Success) {
            Eigen::VectorXd U = solver.solve(F);
            checkTrue("Plate mesh solve succeeded", solver.info() == Eigen::Success);

            const int centerIdx = nodeIndex(n / 2, n / 2);
            const double wCenter = U(centerIdx * 6 + 2);

            const double D = (E * t_plate * t_plate * t_plate) / (12.0 * (1.0 - nu_plate * nu_plate));
            const double alpha = 0.00406;
            const double expectedW = -alpha * q * std::pow(a, 4) / D;

            checkClose("Center deflection vs. Timoshenko simply-supported square plate formula",
                       wCenter, expectedW, 0.05); // ~2% observed, 5% tolerance for mesh-density safety margin
        }
    }

    // ============================================================
    // Test 4: Degenerate Geometry — Graceful Failure (Not a Crash)
    // ============================================================
    std::cout << "\n=== Test 4: Degenerate (Zero-Area) Geometry — Graceful Failure ===\n";
    {
        std::array<Node3D, 4> degenerateCorners = {
            Node3D{"n0", 0.0, 0.0, 0.0},
            Node3D{"n1", 1.0, 0.0, 0.0},
            Node3D{"n2", 2.0, 0.0, 0.0}, // collinear with n0,n1 — zero area
            Node3D{"n3", 3.0, 0.0, 0.0},
        };
        bool threw = false;
        try {
            auto K = computeShellGlobalStiffness(degenerateCorners, t, E, nu);
        } catch (const std::exception&) {
            threw = true;
        }
        checkTrue("Degenerate (collinear/zero-area) geometry throws instead of producing garbage",
                  threw);
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
