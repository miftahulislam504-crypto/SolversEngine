#include "solver.h"
#include "stiffness.h"
#include "shell.h"
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

#ifndef M_PI
// M_PI POSIX extension, সব compiler/standard library তে গ্যারান্টিড না
// (যেমন strict C++17 mode এ MSVC তে অনুপস্থিত থাকতে পারে) — নিজে
// define করে দেওয়া হচ্ছে যাতে এই ফাইল compiler-independent থাকে।
#define M_PI 3.14159265358979323846
#endif

namespace civilos {

namespace {

/**
 * Phase 10n — এই constant টা আগে applyBoundaryConditions() এর ভেতরে
 * local ছিল, এখন file-scope এ তোলা হয়েছে কারণ solveLinearStatic() এর
 * reaction-force recovery (নিচে দেখুন, AnalysisResult::reactionForces)
 * ও এটা প্রয়োজন — দুই জায়গায় আলাদা copy রাখলে ভবিষ্যতে একটা বদলালে
 * অন্যটা বদলাতে ভুলে গিয়ে reaction ভুল হয়ে যাওয়ার ঝুঁকি থাকত, তাই
 * single source of truth।
 */
const double PENALTY_FACTOR = 1e12;

/**
 * "moment" বা "pin" connectionType অনুযায়ী local stiffness matrix
 * modify করার ফাংশন। "pin" হলে উভয় প্রান্তের bending moment DOF
 * (ry, rz) static condensation দিয়ে release করা হয় (stiffness.h এর
 * applyEndReleases() দেখুন, সেখানে পদ্ধতির পূর্ণ ব্যাখ্যা আছে)।
 *
 * নোট: এই মুহূর্তে "pin" মানে উভয় প্রান্তে release (Phase 4a এর
 * Brace default আচরণ)। এক-প্রান্ত-pin/এক-প্রান্ত-moment (partial
 * release) সমর্থনের জন্য connectionType কে per-end string না রেখে
 * দুইটা আলাদা field (startRelease/endRelease) লাগবে — সেটা এখনো
 * frontend TypeScript টাইপে নেই (element.ts এ connectionType একটাই
 * string, per-element), তাই এখানে conservative সিদ্ধান্ত: "pin" =
 * both ends released। ভুল দিকে conservative না — বরং if anything,
 * both-end release করা বেশি flexible structure দেবে (moment বেশি
 * ধরবে না যা বাস্তবে নেই), এটাই নিরাপদ দিক।
 */
Eigen::Matrix<double, 12, 12> getEffectiveLocalStiffness(
    double length,
    const SectionProperties& section,
    const MaterialProperties& material,
    const std::string& connectionType
) {
    const Eigen::Matrix<double, 12, 12> kRigid = computeLocalStiffnessMatrix(length, section, material);

    if (connectionType == "pin") {
        return applyEndReleases(kRigid, /*releaseStart=*/true, /*releaseEnd=*/true);
    }

    return kRigid; // "moment" (default) — কোনো release না
}

} // anonymous namespace

Eigen::SparseMatrix<double> assembleGlobalStiffness(const AnalysisModel& model) {
    const int nodeCount = static_cast<int>(model.nodes.size());
    const int totalDOF = nodeCount * 6;

    // Triplet list ব্যবহার করা হচ্ছে sparse matrix build করার জন্য —
    // এটা Eigen এর সুপারিশকৃত পদ্ধতি বড় sparse matrix তৈরি করার সময়
    // (প্রতিটা element insertion সরাসরি sparse matrix এ করলে অনেক
    // বেশি ধীর হতো, কারণ sparse matrix এর internal storage প্রতিটা
    // insertion এ পুনর্বিন্যাস করতে হয়)।
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(model.elements.size() * 144 + model.shellElements.size() * 576); // frame: 12x12=144, shell: 24x24=576 (upper bound)

    for (const auto& element : model.elements) {
        const Node3D& startNode = model.nodes[element.startNodeIndex];
        const Node3D& endNode = model.nodes[element.endNodeIndex];

        const double dx = endNode.x - startNode.x;
        const double dy = endNode.y - startNode.y;
        const double dz = endNode.z - startNode.z;
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

        const auto kLocal = getEffectiveLocalStiffness(length, element.section, element.material, element.connectionType);
        const auto T = computeTransformationMatrix(startNode, endNode);
        const Eigen::Matrix<double, 12, 12> kGlobal = T.transpose() * kLocal * T;

        // Global DOF index mapping — element এর local 12 DOF কে global
        // system এর সঠিক position এ বসানো। প্রতিটা node এর 6 DOF
        // (node index * 6) থেকে শুরু হয়।
        std::array<int, 12> globalIndices{};
        for (int i = 0; i < 6; ++i) {
            globalIndices[i] = element.startNodeIndex * 6 + i;
            globalIndices[6 + i] = element.endNodeIndex * 6 + i;
        }

        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                if (std::abs(kGlobal(i, j)) > 1e-15) { // exact-zero entry বাদ দেওয়া, sparse efficiency বজায় রাখতে
                    triplets.emplace_back(globalIndices[i], globalIndices[j], kGlobal(i, j));
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Shell Element (Slab/Wall) Assembly — line element এর অনুরূপ যুক্তি,
    // কিন্তু 24x24 matrix ও 4-node connectivity (shell.h এর
    // computeShellGlobalStiffness ব্যবহার করে)।
    // ------------------------------------------------------------------
    for (const auto& shell : model.shellElements) {
        std::array<Node3D, 4> corners = {
            model.nodes[shell.nodeIndices[0]],
            model.nodes[shell.nodeIndices[1]],
            model.nodes[shell.nodeIndices[2]],
            model.nodes[shell.nodeIndices[3]],
        };

        const auto kGlobalShell = computeShellGlobalStiffness(
            corners, shell.thickness, shell.material.elasticModulus, shell.material.poissonsRatio);

        std::array<int, 24> globalIndices{};
        for (int node = 0; node < 4; ++node) {
            for (int d = 0; d < 6; ++d) {
                globalIndices[6 * node + d] = shell.nodeIndices[node] * 6 + d;
            }
        }

        for (int i = 0; i < 24; ++i) {
            for (int j = 0; j < 24; ++j) {
                if (std::abs(kGlobalShell(i, j)) > 1e-15) {
                    triplets.emplace_back(globalIndices[i], globalIndices[j], kGlobalShell(i, j));
                }
            }
        }
    }

    Eigen::SparseMatrix<double> globalStiffness(totalDOF, totalDOF);
    // setFromTriplets স্বয়ংক্রিয়ভাবে duplicate (একই position এ একাধিক
    // triplet) সমষ্টি করে — এটাই কাম্য আচরণ, কারণ একটা node একাধিক
    // element এর সংযোগস্থল হলে তাদের stiffness contribution যোগ হওয়া
    // উচিত (superposition principle, Direct Stiffness Method এর মূল ভিত্তি)।
    globalStiffness.setFromTriplets(triplets.begin(), triplets.end());

    return globalStiffness;
}

void applyBoundaryConditions(
    Eigen::SparseMatrix<double>& globalStiffness,
    const std::vector<BoundaryCondition>& boundaryConditions,
    int totalDOF
) {
    // Penalty factor — global stiffness matrix এর সবচেয়ে বড় diagonal
    // entry এর তুলনায় যথেষ্ট বড় হতে হবে যাতে সংশ্লিষ্ট DOF কার্যত
    // rigid হয়ে যায়, কিন্তু এত বড় না যে floating-point precision loss
    // হয় (numerical ill-conditioning)। 1e12 একটা প্রচলিত, প্রমাণিত মান
    // এই উদ্দেশ্যে (double precision এ ~15-16 significant digits থাকে,
    // তাই 1e12 factor এর সাথে normal-magnitude stiffness যোগ করলেও
    // যথেষ্ট precision অবশিষ্ট থাকে)। file-scope PENALTY_FACTOR constant
    // (উপরে দেখুন) ব্যবহৃত হচ্ছে, local না — solveLinearStatic() এর
    // reaction-force recovery ও এটাই ব্যবহার করে, single source of truth।

    (void)totalDOF; // বর্তমানে সরাসরি ব্যবহৃত হয় না (bound-check এর জন্য ভবিষ্যতে কাজে লাগতে পারে), কিন্তু signature এ রাখা হয়েছে explicit contract হিসেবে

    for (const auto& bc : boundaryConditions) {
        const int base = bc.nodeIndex * 6;
        if (bc.restrainX) globalStiffness.coeffRef(base + 0, base + 0) += PENALTY_FACTOR;
        if (bc.restrainY) globalStiffness.coeffRef(base + 1, base + 1) += PENALTY_FACTOR;
        if (bc.restrainZ) globalStiffness.coeffRef(base + 2, base + 2) += PENALTY_FACTOR;
        if (bc.restrainRx) globalStiffness.coeffRef(base + 3, base + 3) += PENALTY_FACTOR;
        if (bc.restrainRy) globalStiffness.coeffRef(base + 4, base + 4) += PENALTY_FACTOR;
        if (bc.restrainRz) globalStiffness.coeffRef(base + 5, base + 5) += PENALTY_FACTOR;
    }
}

Eigen::VectorXd assembleGlobalLoadVector(const AnalysisModel& model, int totalDOF) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(totalDOF);

    for (const auto& load : model.loads) {
        const int base = load.nodeIndex * 6;
        F(base + 0) += load.fx;
        F(base + 1) += load.fy;
        F(base + 2) += load.fz;
        F(base + 3) += load.mx;
        F(base + 4) += load.my;
        F(base + 5) += load.mz;
    }

    return F;
}

AnalysisResult solveLinearStatic(const AnalysisModel& model) {
    AnalysisResult result;
    result.success = false;

    if (model.nodes.empty()) {
        result.errorMessage = "Model has no nodes";
        return result;
    }
    if (model.elements.empty() && model.shellElements.empty()) {
        result.errorMessage = "Model has no elements (neither frame nor shell)";
        return result;
    }
    if (model.boundaryConditions.empty()) {
        result.errorMessage = "Model has no boundary conditions — structure is unrestrained (unstable)";
        return result;
    }

    const int totalDOF = static_cast<int>(model.nodes.size()) * 6;

    Eigen::SparseMatrix<double> K = assembleGlobalStiffness(model);
    applyBoundaryConditions(K, model.boundaryConditions, totalDOF);
    Eigen::VectorXd F = assembleGlobalLoadVector(model, totalDOF);

    // Sparse Cholesky (LDLT variant, যা indefinite matrix এও কাজ করে
    // যদিও আমাদের ক্ষেত্রে matrix positive-definite হওয়ার কথা একটা
    // সঠিকভাবে-restrained স্থিতিশীল structure এর জন্য)।
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K);

    if (solver.info() != Eigen::Success) {
        result.errorMessage = "Stiffness matrix decomposition failed — structure may be unstable "
                               "(insufficient supports, or a mechanism exists in the model)";
        return result;
    }

    Eigen::VectorXd U = solver.solve(F);

    if (solver.info() != Eigen::Success) {
        result.errorMessage = "Linear system solve failed after successful decomposition — "
                               "this is unexpected and may indicate numerical ill-conditioning";
        return result;
    }

    // প্রতিটা node এর জন্য 6-DOF displacement vector আলাদা করা
    result.nodalDisplacements.resize(model.nodes.size());
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        result.nodalDisplacements[i] = U.segment(static_cast<int>(i) * 6, 6);
    }

    // Phase 10n — Support reaction forces, penalty-method recovery।
    //
    // Global equilibrium: (K_structure + K_penalty) U = F_applied, তাই
    // K_structure U = F_applied − K_penalty U। ডান পাশের −K_penalty U
    // অংশটাই প্রকৃত reaction (ground/support structure কে যে বাড়তি
    // force দিচ্ছে, applied load এর বিপরীতে ভারসাম্য রাখতে) — তাই
    // reaction = −PENALTY_FACTOR × U(dof), পজিটিভ না।
    //
    // এই sign টা প্রথমে ভুল করে পজিটিভ রাখা হয়েছিল, hand-calculated
    // simply-supported-beam টেস্ট কেস (test_solver_integration.cpp,
    // center point load P, উভয় support এ reaction = P/2, upward/
    // positive-Y হওয়ার কথা যেহেতু load downward/negative-Y) দিয়ে
    // ধরা পড়ে — negative sign ছাড়া reaction load এর same direction এ
    // (downward) আসছিল, যেটা physically ভুল (একটা downward load কে
    // ধরে রাখতে support সবসময় upward push করে, বাস্তব সেতু/কলামে যেমন
    // হয়)। এখন negative sign সহ P/2 upward (positive) মেলে, ও
    // element-end-force থেকে independently derive করা reaction এর
    // সাথেও বিট-পারফেক্ট মেলে (দুইটা সম্পূর্ণ ভিন্ন পদ্ধতি একই সংখ্যা
    // দিচ্ছে, যা sign+magnitude উভয়ের জন্য শক্তিশালী প্রমাণ)।
    //
    // PENALTY_FACTOR এখানে applyBoundaryConditions() এর সাথে হুবহু
    // মিলিয়ে রাখা আবশ্যক — দুই জায়গায় আলাদা মান হলে reaction ভুল হবে,
    // তাই একই file-scope constant।
    result.reactionForces.resize(model.boundaryConditions.size());
    for (size_t i = 0; i < model.boundaryConditions.size(); ++i) {
        const auto& bc = model.boundaryConditions[i];
        const int base = bc.nodeIndex * 6;
        Eigen::VectorXd reaction = Eigen::VectorXd::Zero(6);
        if (bc.restrainX) reaction(0) = -PENALTY_FACTOR * U(base + 0);
        if (bc.restrainY) reaction(1) = -PENALTY_FACTOR * U(base + 1);
        if (bc.restrainZ) reaction(2) = -PENALTY_FACTOR * U(base + 2);
        if (bc.restrainRx) reaction(3) = -PENALTY_FACTOR * U(base + 3);
        if (bc.restrainRy) reaction(4) = -PENALTY_FACTOR * U(base + 4);
        if (bc.restrainRz) reaction(5) = -PENALTY_FACTOR * U(base + 5);
        result.reactionForces[i] = reaction;
    }

    // প্রতিটা element এর end force হিসাব — local coordinate এ
    // (f_local = k_local * T * u_global_element)
    result.elementEndForces.resize(model.elements.size());
    for (size_t e = 0; e < model.elements.size(); ++e) {
        const auto& element = model.elements[e];
        const Node3D& startNode = model.nodes[element.startNodeIndex];
        const Node3D& endNode = model.nodes[element.endNodeIndex];

        const double dx = endNode.x - startNode.x;
        const double dy = endNode.y - startNode.y;
        const double dz = endNode.z - startNode.z;
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

        const auto kLocal = getEffectiveLocalStiffness(length, element.section, element.material, element.connectionType);
        const auto T = computeTransformationMatrix(startNode, endNode);

        Eigen::VectorXd uElementGlobal(12);
        uElementGlobal.segment(0, 6) = U.segment(element.startNodeIndex * 6, 6);
        uElementGlobal.segment(6, 6) = U.segment(element.endNodeIndex * 6, 6);

        Eigen::VectorXd uElementLocal = T * uElementGlobal;
        result.elementEndForces[e] = kLocal * uElementLocal;
    }

    result.success = true;
    return result;
}

Eigen::SparseMatrix<double> assembleGlobalMass(const AnalysisModel& model) {
    const int nodeCount = static_cast<int>(model.nodes.size());
    const int totalDOF = nodeCount * 6;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(model.elements.size() * 144);

    for (const auto& element : model.elements) {
        const Node3D& startNode = model.nodes[element.startNodeIndex];
        const Node3D& endNode = model.nodes[element.endNodeIndex];

        const double dx = endNode.x - startNode.x;
        const double dy = endNode.y - startNode.y;
        const double dz = endNode.z - startNode.z;
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

        // নোট: mass matrix এ "pin" connectionType এর জন্য কোনো
        // condensation প্রয়োগ করা হয় না (stiffness.h এর
        // computeLocalMassMatrix() docstring এ ব্যাখ্যা করা কারণ) —
        // rigid ও pin-connected element এর mass matrix একই সূত্র
        // ব্যবহার করে গণনা করা হয়।
        const auto mLocal = computeLocalMassMatrix(length, element.section, element.material);
        const auto T = computeTransformationMatrix(startNode, endNode);
        const Eigen::Matrix<double, 12, 12> mGlobal = T.transpose() * mLocal * T;

        std::array<int, 12> globalIndices{};
        for (int i = 0; i < 6; ++i) {
            globalIndices[i] = element.startNodeIndex * 6 + i;
            globalIndices[6 + i] = element.endNodeIndex * 6 + i;
        }

        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                if (std::abs(mGlobal(i, j)) > 1e-15) {
                    triplets.emplace_back(globalIndices[i], globalIndices[j], mGlobal(i, j));
                }
            }
        }
    }

    Eigen::SparseMatrix<double> globalMass(totalDOF, totalDOF);
    globalMass.setFromTriplets(triplets.begin(), triplets.end());

    return globalMass;
}

ModalAnalysisResult solveModalAnalysis(const AnalysisModel& model, int numModes) {
    ModalAnalysisResult result;
    result.success = false;

    if (model.nodes.empty()) {
        result.errorMessage = "Model has no nodes";
        return result;
    }
    if (model.elements.empty()) {
        result.errorMessage = "Model has no elements";
        return result;
    }
    if (model.boundaryConditions.empty()) {
        result.errorMessage = "Model has no boundary conditions — structure is unrestrained (unstable)";
        return result;
    }

    const int totalDOF = static_cast<int>(model.nodes.size()) * 6;

    Eigen::SparseMatrix<double> Ksparse = assembleGlobalStiffness(model);
    Eigen::SparseMatrix<double> Msparse = assembleGlobalMass(model);

    // ------------------------------------------------------------------
    // Elimination method — restrained DOF সম্পূর্ণ বাদ দিয়ে reduced
    // system বানানো (উপরে solver.h এর কমেন্টে কারণ ব্যাখ্যা করা)।
    // ------------------------------------------------------------------
    std::vector<bool> isRestrained(totalDOF, false);
    for (const auto& bc : model.boundaryConditions) {
        const int base = bc.nodeIndex * 6;
        if (bc.restrainX) isRestrained[base + 0] = true;
        if (bc.restrainY) isRestrained[base + 1] = true;
        if (bc.restrainZ) isRestrained[base + 2] = true;
        if (bc.restrainRx) isRestrained[base + 3] = true;
        if (bc.restrainRy) isRestrained[base + 4] = true;
        if (bc.restrainRz) isRestrained[base + 5] = true;
    }

    // free (unrestrained) DOF এর global index → reduced-system index ম্যাপিং
    std::vector<int> freeDofGlobalIndex;
    freeDofGlobalIndex.reserve(totalDOF);
    for (int i = 0; i < totalDOF; ++i) {
        if (!isRestrained[i]) freeDofGlobalIndex.push_back(i);
    }

    const int reducedDOF = static_cast<int>(freeDofGlobalIndex.size());

    if (reducedDOF == 0) {
        result.errorMessage = "সব DOF restrained — কোনো dynamic degree of freedom অবশিষ্ট নেই "
                               "(পুরো structure সম্পূর্ণ fixed, modal analysis অর্থহীন)";
        return result;
    }

    // Dense matrix এ reduced system বসানো — GeneralizedSelfAdjointEigenSolver
    // dense matrix দাবি করে (sparse generalized eigensolver Eigen core এ
    // নেই, উপরে solver.h এ ব্যাখ্যা করা কারণে external dependency এড়ানো
    // হয়েছে)।
    Eigen::MatrixXd Kreduced(reducedDOF, reducedDOF);
    Eigen::MatrixXd Mreduced(reducedDOF, reducedDOF);

    Eigen::MatrixXd Kdense = Eigen::MatrixXd(Ksparse);
    Eigen::MatrixXd Mdense = Eigen::MatrixXd(Msparse);

    for (int i = 0; i < reducedDOF; ++i) {
        for (int j = 0; j < reducedDOF; ++j) {
            Kreduced(i, j) = Kdense(freeDofGlobalIndex[i], freeDofGlobalIndex[j]);
            Mreduced(i, j) = Mdense(freeDofGlobalIndex[i], freeDofGlobalIndex[j]);
        }
    }

    // ------------------------------------------------------------------
    // Generalized Symmetric Eigenvalue Problem: K φ = λ M φ, যেখানে
    // λ = ω² (angular frequency এর বর্গ)। M অবশ্যই positive-definite
    // হতে হবে এই solver এর জন্য (Cholesky-based ABx=λx পদ্ধতি ভেতরে
    // ব্যবহার করে) — একটা physically valid model এ (প্রতিটা element এর
    // material density > 0) এটা সবসময় সত্য।
    // ------------------------------------------------------------------
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(Kreduced, Mreduced);

    if (eigenSolver.info() != Eigen::Success) {
        result.errorMessage = "Generalized eigenvalue সমাধান ব্যর্থ হয়েছে — সম্ভাব্য কারণ: mass "
                               "matrix positive-definite না (কোনো element এর material density "
                               "শূন্য বা ঋণাত্মক আছে কিনা যাচাই করুন), অথবা numerical ill-conditioning";
        return result;
    }

    // Eigen এর GeneralizedSelfAdjointEigenSolver eigenvalue গুলো ছোট থেকে
    // বড় ক্রমে সাজানো ফেরত দেয় — এটাই আমাদের দরকার (সবচেয়ে নিচু
    // frequency mode গুলো সাধারণত সবচেয়ে গুরুত্বপূর্ণ, সেগুলো প্রথমে)।
    const Eigen::VectorXd& eigenvalues = eigenSolver.eigenvalues();  // λ = ω², rad²/s²
    const Eigen::MatrixXd& eigenvectors = eigenSolver.eigenvectors(); // প্রতিটা কলাম একটা mode shape (reduced DOF স্কেলে)

    const int availableModes = static_cast<int>(eigenvalues.size());
    const int modesToReturn = (numModes <= 0) ? availableModes : std::min(numModes, availableModes);

    result.naturalFrequenciesHz.reserve(modesToReturn);
    result.angularFrequenciesRadPerSec.reserve(modesToReturn);
    result.modeShapes.reserve(modesToReturn);

    for (int m = 0; m < modesToReturn; ++m) {
        double lambda = eigenvalues(m);

        // Numerical noise থেকে সামান্য ঋণাত্মক eigenvalue আসতে পারে
        // (rigid-body mode এর কাছাকাছি বা zero-frequency mode এ,
        // floating-point precision এর কারণে) — physically এটা শূন্য
        // ধরা উচিত, ঋণাত্মক frequency অর্থহীন। clamp করে 0 এ বসানো
        // হচ্ছে sqrt(negative) থেকে NaN এড়াতে।
        if (lambda < 0.0) lambda = 0.0;

        const double omega = std::sqrt(lambda); // rad/s
        const double freqHz = omega / (2.0 * M_PI);

        result.angularFrequenciesRadPerSec.push_back(omega);
        result.naturalFrequenciesHz.push_back(freqHz);

        // Reduced eigenvector কে পূর্ণ global DOF স্কেলে ফিরিয়ে আনা
        // (restrained DOF এ 0 বসিয়ে) — nodalDisplacements এর মতো একই
        // shape যাতে frontend একইভাবে consume করতে পারে।
        Eigen::VectorXd fullModeShape = Eigen::VectorXd::Zero(totalDOF);
        for (int i = 0; i < reducedDOF; ++i) {
            fullModeShape(freeDofGlobalIndex[i]) = eigenvectors(i, m);
        }

        // নোট: GeneralizedSelfAdjointEigenSolver থেকে পাওয়া eigenvector
        // ইতিমধ্যে M-normalized (φᵀMφ = 1) — Eigen library এর ডকুমেন্টেড
        // আচরণ generalized self-adjoint eigenproblem এর জন্য, তাই আলাদা
        // করে normalize করার প্রয়োজন নেই। এটা types.h এর
        // ModalAnalysisResult docstring এ উল্লেখিত convention এর সাথে
        // সঙ্গতিপূর্ণ।

        std::vector<Eigen::VectorXd> perNodeShape(model.nodes.size());
        for (size_t i = 0; i < model.nodes.size(); ++i) {
            perNodeShape[i] = fullModeShape.segment(static_cast<int>(i) * 6, 6);
        }
        result.modeShapes.push_back(std::move(perNodeShape));
    }

    result.numModesComputed = modesToReturn;
    result.success = true;
    return result;
}

BucklingAnalysisResult solveLinearBuckling(const AnalysisModel& model, int numModes) {
    BucklingAnalysisResult result;
    result.success = false;

    if (model.nodes.empty()) {
        result.errorMessage = "Model has no nodes";
        return result;
    }
    if (model.elements.empty()) {
        result.errorMessage = "Model has no elements";
        return result;
    }
    if (model.boundaryConditions.empty()) {
        result.errorMessage = "Model has no boundary conditions — structure is unrestrained (unstable)";
        return result;
    }
    if (model.loads.empty()) {
        result.errorMessage = "Model has no loads — Linear Buckling Analysis একটা নির্দিষ্ট load "
                               "pattern এর সাপেক্ষে axial force distribution দরকার (types.h এর "
                               "BucklingAnalysisResult docstring দেখুন); কোনো load ছাড়া সব axial "
                               "force শূন্য, buckling analysis অর্থহীন";
        return result;
    }

    // ------------------------------------------------------------------
    // ধাপ ১: দেওয়া load pattern দিয়ে Linear Static সমাধান — axial force
    // distribution বের করতে (উপরে solver.h এ ব্যাখ্যা করা কারণে)।
    // ------------------------------------------------------------------
    AnalysisResult staticResult = solveLinearStatic(model);
    if (!staticResult.success) {
        result.errorMessage = "Buckling Analysis এর জন্য প্রয়োজনীয় প্রাথমিক Linear Static সমাধান "
                               "ব্যর্থ হয়েছে: " + staticResult.errorMessage;
        return result;
    }

    const int totalDOF = static_cast<int>(model.nodes.size()) * 6;

    // ------------------------------------------------------------------
    // ধাপ ২: Global Geometric Stiffness (Kg) assembly — প্রতিটা element
    // এর axial force (staticResult.elementEndForces[e](0) = local index 0
    // = start-node axial, computeLocalStiffnessMatrix এর axial DOF
    // convention অনুযায়ী) ব্যবহার করে।
    // ------------------------------------------------------------------
    std::vector<Eigen::Triplet<double>> kgTriplets;
    kgTriplets.reserve(model.elements.size() * 144);

    for (size_t e = 0; e < model.elements.size(); ++e) {
        const auto& element = model.elements[e];
        const Node3D& startNode = model.nodes[element.startNodeIndex];
        const Node3D& endNode = model.nodes[element.endNodeIndex];

        const double dx = endNode.x - startNode.x;
        const double dy = endNode.y - startNode.y;
        const double dz = endNode.z - startNode.z;
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

        const double axialForce = staticResult.elementEndForces[e](0);

        const auto kgLocal = computeLocalGeometricStiffnessMatrix(length, axialForce);
        const auto T = computeTransformationMatrix(startNode, endNode);
        const Eigen::Matrix<double, 12, 12> kgGlobal = T.transpose() * kgLocal * T;

        std::array<int, 12> globalIndices{};
        for (int i = 0; i < 6; ++i) {
            globalIndices[i] = element.startNodeIndex * 6 + i;
            globalIndices[6 + i] = element.endNodeIndex * 6 + i;
        }

        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                if (std::abs(kgGlobal(i, j)) > 1e-15) {
                    kgTriplets.emplace_back(globalIndices[i], globalIndices[j], kgGlobal(i, j));
                }
            }
        }
    }

    Eigen::SparseMatrix<double> Kg(totalDOF, totalDOF);
    Kg.setFromTriplets(kgTriplets.begin(), kgTriplets.end());

    Eigen::SparseMatrix<double> Ksparse = assembleGlobalStiffness(model);

    // ------------------------------------------------------------------
    // ধাপ ৩: Elimination method দিয়ে boundary condition (Modal Analysis
    // এর মতোই কারণ — penalty method এখানেও eigenvalue দূষিত করবে)।
    // ------------------------------------------------------------------
    std::vector<bool> isRestrained(totalDOF, false);
    for (const auto& bc : model.boundaryConditions) {
        const int base = bc.nodeIndex * 6;
        if (bc.restrainX) isRestrained[base + 0] = true;
        if (bc.restrainY) isRestrained[base + 1] = true;
        if (bc.restrainZ) isRestrained[base + 2] = true;
        if (bc.restrainRx) isRestrained[base + 3] = true;
        if (bc.restrainRy) isRestrained[base + 4] = true;
        if (bc.restrainRz) isRestrained[base + 5] = true;
    }

    std::vector<int> freeDofGlobalIndex;
    freeDofGlobalIndex.reserve(totalDOF);
    for (int i = 0; i < totalDOF; ++i) {
        if (!isRestrained[i]) freeDofGlobalIndex.push_back(i);
    }

    const int reducedDOF = static_cast<int>(freeDofGlobalIndex.size());

    if (reducedDOF == 0) {
        result.errorMessage = "সব DOF restrained — কোনো buckling degree of freedom অবশিষ্ট নেই";
        return result;
    }

    Eigen::MatrixXd Kdense = Eigen::MatrixXd(Ksparse);
    Eigen::MatrixXd KgDense = Eigen::MatrixXd(Kg);

    Eigen::MatrixXd Kreduced(reducedDOF, reducedDOF);
    Eigen::MatrixXd KgReduced(reducedDOF, reducedDOF);

    for (int i = 0; i < reducedDOF; ++i) {
        for (int j = 0; j < reducedDOF; ++j) {
            Kreduced(i, j) = Kdense(freeDofGlobalIndex[i], freeDofGlobalIndex[j]);
            KgReduced(i, j) = KgDense(freeDofGlobalIndex[i], freeDofGlobalIndex[j]);
        }
    }

    // ------------------------------------------------------------------
    // ধাপ ৪: Generalized Eigenvalue Problem।
    //
    // মূল সমীকরণ: K φ = -λ Kg φ। সরাসরি এভাবে (Kreduced, -KgReduced)
    // solve করা যেত যদি -Kg positive-definite হতো — কিন্তু বাস্তবে
    // Kg সবসময় singular (rank-deficient): axial ও torsion DOF এ কোনো
    // geometric-stiffness contribution নেই (stiffness.h এর
    // computeLocalGeometricStiffnessMatrix() এর ইচ্ছাকৃত সরলীকরণ —
    // শুধু bending-coupled geometric effect ধরা হয়েছে)। ফলে -Kg এর
    // সেই dimension গুলোতে eigenvalue শূন্য, positive-definite হওয়া
    // গাণিতিকভাবেই অসম্ভব — এটা কোনো model-নির্দিষ্ট edge case না,
    // *প্রতিটা* মডেলেই ঘটবে এই sparsity pattern এর কারণে। তাই
    // Eigen::GeneralizedSelfAdjointEigenSolver(K, -Kg) directly ব্যবহার
    // করলে সবসময় ব্যর্থ হতো (Cholesky-based ভেতরের পদ্ধতি B কে strictly
    // positive-definite দাবি করে)।
    //
    // সমাধান: সমীকরণ উল্টে reformulate করা —
    //   K φ = -λ Kg φ
    //   ⟹ (-Kg) φ = (1/λ) K φ         [উভয় পাশে -1/λ দিয়ে গুণ, λ≠0 ধরে]
    //   ⟹ (-Kg) φ = μ K φ,  যেখানে μ = 1/λ
    //
    // এখন B matrix = Kreduced (properly-restrained structure এ সবসময়
    // strictly positive-definite — কোনো rigid-body mode অবশিষ্ট নেই
    // reduced/eliminated system এ) ব্যবহার করা যায়, A matrix = -KgReduced
    // (singular হলেও চলবে, শুধু B এর positive-definiteness দরকার এই
    // solver এ)। সমাধান শেষে μ থেকে λ = 1/μ তে ফিরিয়ে আনা হয় নিচে।
    //
    // এই reformulation numerically সম্পূর্ণ সমতুল্য (mathematically
    // identical eigenvector, eigenvalue শুধু reciprocal সম্পর্কে
    // ভিন্ন) — কোনো approximation বা accuracy loss নেই, শুধু matrix
    // pencil (A,B) এর ভূমিকা অদলবদল।
    // ------------------------------------------------------------------
    Eigen::MatrixXd negKgReduced = -KgReduced;

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(negKgReduced, Kreduced);

    if (eigenSolver.info() != Eigen::Success) {
        result.errorMessage = "Generalized eigenvalue সমাধান ব্যর্থ হয়েছে — সম্ভাব্য কারণ: "
                               "stiffness matrix (K) positive-definite না (structure unstable, "
                               "যদিও এটা আগেই solveLinearStatic ধাপে ধরা পড়ার কথা), অথবা numerical "
                               "ill-conditioning";
        return result;
    }

    const Eigen::VectorXd& muValues = eigenSolver.eigenvalues(); // μ = 1/λ
    const Eigen::MatrixXd& eigenvectors = eigenSolver.eigenvectors();

    // μ থেকে λ (critical load factor) এ ফিরিয়ে আনা: λ = 1/μ। μ≈0 হলে
    // (যা ঘটে সেই eigenvector গুলোর জন্য যেখানে -Kg এর null-space এর
    // দিকে align করে, অর্থাৎ axial/torsion-dominated non-physical
    // "buckling" mode — বাস্তবে buckling না, শুধু geometric stiffness
    // এর null-space artifact) λ অসীম হয়ে যায় (কোনো finite critical
    // load এ সেই mode ঘটবে না) — এই mode গুলো ফলাফল থেকে বাদ দেওয়া
    // হচ্ছে (physically অর্থহীন, |λ|=∞ কখনো governing case হতে পারে
    // না)।
    const double MU_ZERO_THRESHOLD = 1e-9;

    std::vector<std::pair<double, int>> validLambdaWithIndex; // (|λ|, original eigen-index)
    for (int i = 0; i < muValues.size(); ++i) {
        if (std::abs(muValues(i)) > MU_ZERO_THRESHOLD) {
            const double lambda = 1.0 / muValues(i);
            validLambdaWithIndex.emplace_back(std::abs(lambda), i);
        }
    }

    if (validLambdaWithIndex.empty()) {
        result.errorMessage = "কোনো finite critical load factor পাওয়া যায়নি — এই মডেলে buckling "
                               "প্রাসঙ্গিক geometric stiffness শূন্য হতে পারে (সব element এ axial "
                               "force প্রায় শূন্য, যেমন শুধু bending load থাকলে কোনো axial force "
                               "ছাড়া)";
        return result;
    }

    // সবচেয়ে ছোট |λ| (সবচেয়ে critical/nearest buckling mode) অনুযায়ী
    // sort করা — Modal Analysis এর মতোই যুক্তি, সবচেয়ে গুরুত্বপূর্ণ
    // mode প্রথমে।
    std::sort(validLambdaWithIndex.begin(), validLambdaWithIndex.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    const int availableModes = static_cast<int>(validLambdaWithIndex.size());
    const int modesToReturn = (numModes <= 0) ? availableModes : std::min(numModes, availableModes);

    result.criticalLoadFactors.reserve(modesToReturn);
    result.bucklingModeShapes.reserve(modesToReturn);

    for (int m = 0; m < modesToReturn; ++m) {
        const int eigenIdx = validLambdaWithIndex[m].second;
        const double lambda = 1.0 / muValues(eigenIdx);
        result.criticalLoadFactors.push_back(lambda);

        Eigen::VectorXd fullModeShape = Eigen::VectorXd::Zero(totalDOF);
        for (int i = 0; i < reducedDOF; ++i) {
            fullModeShape(freeDofGlobalIndex[i]) = eigenvectors(i, eigenIdx);
        }

        std::vector<Eigen::VectorXd> perNodeShape(model.nodes.size());
        for (size_t i = 0; i < model.nodes.size(); ++i) {
            perNodeShape[i] = fullModeShape.segment(static_cast<int>(i) * 6, 6);
        }
        result.bucklingModeShapes.push_back(std::move(perNodeShape));
    }

    result.numModesComputed = modesToReturn;
    result.success = true;
    return result;
}

PDeltaAnalysisResult solvePDelta(const AnalysisModel& model) {
    PDeltaAnalysisResult result;
    result.success = false;

    if (model.nodes.empty()) {
        result.errorMessage = "Model has no nodes";
        return result;
    }
    if (model.elements.empty()) {
        result.errorMessage = "Model has no elements";
        return result;
    }
    if (model.boundaryConditions.empty()) {
        result.errorMessage = "Model has no boundary conditions — structure is unrestrained (unstable)";
        return result;
    }
    if (model.loads.empty()) {
        result.errorMessage = "Model has no loads — P-Delta Analysis একটা নির্দিষ্ট load pattern এর "
                               "সাপেক্ষে axial force ও second-order effect হিসাব করে (types.h এর "
                               "PDeltaAnalysisResult docstring দেখুন); কোনো load ছাড়া P-Delta effect "
                               "শূন্য, Linear Static থেকে ফলাফল অভিন্ন হবে — এই ফাংশন ব্যবহারের "
                               "প্রয়োজন নেই সেক্ষেত্রে";
        return result;
    }

    // ------------------------------------------------------------------
    // ধাপ ১: First-order (Linear Static) সমাধান — axial force distribution
    // বের করতে (solveLinearBuckling এর মতোই একই প্রাথমিক ধাপ)।
    // ------------------------------------------------------------------
    AnalysisResult staticResult = solveLinearStatic(model);
    if (!staticResult.success) {
        result.errorMessage = "P-Delta Analysis এর জন্য প্রয়োজনীয় প্রাথমিক (first-order) Linear "
                               "Static সমাধান ব্যর্থ হয়েছে: " + staticResult.errorMessage;
        return result;
    }

    result.firstOrderAxialForces.reserve(model.elements.size());
    for (size_t e = 0; e < model.elements.size(); ++e) {
        result.firstOrderAxialForces.push_back(staticResult.elementEndForces[e](0));
    }

    const int totalDOF = static_cast<int>(model.nodes.size()) * 6;

    // ------------------------------------------------------------------
    // ধাপ ২: Global Geometric Stiffness (Kg) assembly — solveLinearBuckling
    // এর সাথে অভিন্ন যুক্তি (একই axial-force-থেকে-Kg পদ্ধতি)।
    // ------------------------------------------------------------------
    std::vector<Eigen::Triplet<double>> kgTriplets;
    kgTriplets.reserve(model.elements.size() * 144);

    for (size_t e = 0; e < model.elements.size(); ++e) {
        const auto& element = model.elements[e];
        const Node3D& startNode = model.nodes[element.startNodeIndex];
        const Node3D& endNode = model.nodes[element.endNodeIndex];

        const double dx = endNode.x - startNode.x;
        const double dy = endNode.y - startNode.y;
        const double dz = endNode.z - startNode.z;
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

        const double axialForce = staticResult.elementEndForces[e](0);

        const auto kgLocal = computeLocalGeometricStiffnessMatrix(length, axialForce);
        const auto T = computeTransformationMatrix(startNode, endNode);
        const Eigen::Matrix<double, 12, 12> kgGlobal = T.transpose() * kgLocal * T;

        std::array<int, 12> globalIndices{};
        for (int i = 0; i < 6; ++i) {
            globalIndices[i] = element.startNodeIndex * 6 + i;
            globalIndices[6 + i] = element.endNodeIndex * 6 + i;
        }

        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                if (std::abs(kgGlobal(i, j)) > 1e-15) {
                    kgTriplets.emplace_back(globalIndices[i], globalIndices[j], kgGlobal(i, j));
                }
            }
        }
    }

    Eigen::SparseMatrix<double> Kg(totalDOF, totalDOF);
    Kg.setFromTriplets(kgTriplets.begin(), kgTriplets.end());

    // ------------------------------------------------------------------
    // ধাপ ৩: Modified stiffness (K + Kg) দিয়ে পুনরায় সমাধান।
    //
    // সাইন নোট: Kg computeLocalGeometricStiffnessMatrix() এর আউটপুট,
    // যেখানে ইতিমধ্যে compression-positive input থেকে সঠিক দিকে (tension
    // stiffness বাড়ায়, compression কমায়) sign-adjusted করা হয়েছে
    // (stiffness.h এর ব্যাখ্যা দেখুন)। তাই এখানে সরাসরি K + Kg (যোগ,
    // বিয়োগ না) — Linear Buckling এর মতো আলাদা sign-flip আবশ্যক না,
    // কারণ এটা eigenvalue problem না, শুধু matrix-addition-based static
    // solve।
    //
    // Boundary condition: penalty method (solveLinearStatic এর মতোই,
    // eliminaton method না — solver.h docstring এ কারণ ব্যাখ্যা করা)।
    // ------------------------------------------------------------------
    Eigen::SparseMatrix<double> Kfirst = assembleGlobalStiffness(model);
    Eigen::SparseMatrix<double> Kmodified = Kfirst + Kg;
    applyBoundaryConditions(Kmodified, model.boundaryConditions, totalDOF);

    Eigen::VectorXd F = assembleGlobalLoadVector(model, totalDOF);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(Kmodified);

    if (solver.info() != Eigen::Success) {
        result.errorMessage = "P-Delta modified stiffness matrix (K + Kg) decomposition ব্যর্থ "
                               "হয়েছে — এটা সাধারণত মানে structure এই load এ geometric "
                               "instability এর কাছাকাছি বা তার বেশি চলে গেছে (critical buckling "
                               "load-এর কাছাকাছি বা অতিক্রম করেছে, types.h এর "
                               "BucklingAnalysisResult সংক্রান্ত দেখুন)। Linear Buckling Analysis "
                               "চালিয়ে critical load factor যাচাই করার পরামর্শ দেওয়া হচ্ছে।";
        return result;
    }

    Eigen::VectorXd Updelta = solver.solve(F);

    if (solver.info() != Eigen::Success) {
        result.errorMessage = "P-Delta linear system solve ব্যর্থ হয়েছে decomposition সফল হওয়ার "
                               "পরেও — এটা অপ্রত্যাশিত এবং numerical ill-conditioning এর ইঙ্গিত "
                               "দিতে পারে";
        return result;
    }

    result.nodalDisplacements.resize(model.nodes.size());
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        result.nodalDisplacements[i] = Updelta.segment(static_cast<int>(i) * 6, 6);
    }

    // Element end force হিসাব — Kg-modified stiffness না, বরং প্রতিটা
    // element এর মূল (material/geometric elastic) local stiffness
    // ব্যবহার করে (f_local = k_local * T * u_pdelta_global_element)।
    // এটা ইচ্ছাকৃত: element internal force সবসময় material stiffness
    // থেকে আসে (elastic constitutive relation), Kg কেবল global
    // equilibrium এ second-order geometric effect ধরার জন্য — সেটা
    // ইতিমধ্যে Updelta তে প্রতিফলিত হয়ে গেছে (বড় displacement)।
    result.elementEndForces.resize(model.elements.size());
    for (size_t e = 0; e < model.elements.size(); ++e) {
        const auto& element = model.elements[e];
        const Node3D& startNode = model.nodes[element.startNodeIndex];
        const Node3D& endNode = model.nodes[element.endNodeIndex];

        const double dx = endNode.x - startNode.x;
        const double dy = endNode.y - startNode.y;
        const double dz = endNode.z - startNode.z;
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

        const auto kLocal = getEffectiveLocalStiffness(length, element.section, element.material, element.connectionType);
        const auto T = computeTransformationMatrix(startNode, endNode);

        Eigen::VectorXd uElementGlobal(12);
        uElementGlobal.segment(0, 6) = Updelta.segment(element.startNodeIndex * 6, 6);
        uElementGlobal.segment(6, 6) = Updelta.segment(element.endNodeIndex * 6, 6);

        Eigen::VectorXd uElementLocal = T * uElementGlobal;
        result.elementEndForces[e] = kLocal * uElementLocal;
    }

    // Displacement amplification ratio — P-Delta displacement বনাম
    // first-order displacement এর সবচেয়ে বড় অনুপাত (P-Delta effect এর
    // overall মাত্রা caller কে দ্রুত বোঝাতে)।
    double maxRatio = 1.0;
    for (int i = 0; i < totalDOF; ++i) {
        const double firstOrderVal = staticResult.nodalDisplacements[i / 6](i % 6);
        const double pDeltaVal = Updelta(i);
        if (std::abs(firstOrderVal) > 1e-9) {
            const double ratio = std::abs(pDeltaVal / firstOrderVal);
            if (std::isfinite(ratio) && ratio > maxRatio) {
                maxRatio = ratio;
            }
        }
    }
    result.maxDisplacementAmplificationRatio = maxRatio;

    result.success = true;
    return result;
}

namespace {

/**
 * Design spectrum এর piecewise-linear interpolation — একটা period T এর
 * জন্য spectral acceleration Sa (g একক) বের করে। spectrum অবশ্যই
 * ছোট-থেকে-বড় periodSec ক্রমে সাজানো থাকতে হবে (caller এর দায়িত্ব,
 * response_spectrum.py এই ক্রম মেনে তৈরি করে)।
 *
 * T সীমার বাইরে গেলে (T < spectrum[0].periodSec অথবা T >
 * spectrum.back().periodSec) নিকটতম endpoint এর Sa ব্যবহার করা হয়
 * (flat extrapolation) — এটা একটা conservative-না-হলেও-reasonable
 * fallback, কারণ বাস্তবে খুব উঁচু mode (খুব ছোট T) এ Sa সাধারণত ছোট ও
 * প্রায়-constant থাকে, আর খুব বড় T এ ও তাই।
 */
double interpolateSpectralAcceleration(
    const std::vector<ResponseSpectrumPoint>& spectrum,
    double periodSec
) {
    if (spectrum.empty()) return 0.0;
    if (spectrum.size() == 1) return spectrum[0].spectralAccelerationG;

    if (periodSec <= spectrum.front().periodSec) return spectrum.front().spectralAccelerationG;
    if (periodSec >= spectrum.back().periodSec) return spectrum.back().spectralAccelerationG;

    for (size_t i = 0; i + 1 < spectrum.size(); ++i) {
        const auto& p0 = spectrum[i];
        const auto& p1 = spectrum[i + 1];
        if (periodSec >= p0.periodSec && periodSec <= p1.periodSec) {
            const double span = p1.periodSec - p0.periodSec;
            if (span < 1e-12) return p0.spectralAccelerationG; // duplicate period point, বিভাজন এড়াতে
            const double ratio = (periodSec - p0.periodSec) / span;
            return p0.spectralAccelerationG + ratio * (p1.spectralAccelerationG - p0.spectralAccelerationG);
        }
    }
    return spectrum.back().spectralAccelerationG; // অপৌঁছনীয়, কিন্তু নিরাপত্তার জন্য
}

/**
 * CQC (Complete Quadratic Combination) modal correlation coefficient —
 * Der Kiureghian (1981) সূত্র, সব commercial software (SAP2000, ETABS)
 * এ ব্যবহৃত standard formula:
 *
 *   ρᵢⱼ = 8ζ²(1+β)β^1.5 / [(1-β²)² + 4ζ²β(1+β)²]
 *
 * যেখানে β = ωⱼ/ωᵢ (frequency ratio), ζ = damping ratio (সব mode এর
 * জন্য একই ধরা হয়েছে এই সংস্করণে)। i=j হলে β=1, ρ=1 (স্বাভাবিক —
 * একটা mode এর নিজের সাথে correlation সম্পূর্ণ)।
 *
 * β→0 বা β→∞ (দূরের frequency) এ ρ→0, তখন CQC≈SRSS (Square Root of
 * Sum of Squares, শুধু diagonal term)। কাছাকাছি frequency (β≈1) এ
 * ρ significant থাকে, cross-term গুরুত্বপূর্ণ হয়ে ওঠে — এটাই CQC কে
 * SRSS এর চেয়ে বেশি নির্ভুল করে closely-spaced-mode structure এ
 * (যেমন torsionally-irregular building, বা symmetric building এ
 * repeated/near-repeated frequency)।
 */
double computeCQCCorrelationCoefficient(double omegaI, double omegaJ, double dampingRatio) {
    if (omegaI < 1e-9 || omegaJ < 1e-9) return (std::abs(omegaI - omegaJ) < 1e-9) ? 1.0 : 0.0;

    const double beta = omegaJ / omegaI;
    const double zeta = dampingRatio;

    const double numerator = 8.0 * zeta * zeta * (1.0 + beta) * std::pow(beta, 1.5);
    const double onemBeta2 = 1.0 - beta * beta;
    const double denominator = onemBeta2 * onemBeta2 + 4.0 * zeta * zeta * beta * (1.0 + beta) * (1.0 + beta);

    if (denominator < 1e-12) return 1.0; // β=1 এর কাছাকাছি সীমায়, numerically ρ→1
    return numerator / denominator;
}

} // anonymous namespace

ResponseSpectrumAnalysisResult solveResponseSpectrum(
    const AnalysisModel& model,
    const std::vector<ResponseSpectrumPoint>& spectrum,
    int directionDOF,
    double dampingRatio,
    int numModes
) {
    ResponseSpectrumAnalysisResult result;
    result.success = false;

    if (model.nodes.empty()) {
        result.errorMessage = "Model has no nodes";
        return result;
    }
    if (model.elements.empty()) {
        result.errorMessage = "Model has no elements";
        return result;
    }
    if (model.boundaryConditions.empty()) {
        result.errorMessage = "Model has no boundary conditions — structure is unrestrained (unstable)";
        return result;
    }
    if (spectrum.empty()) {
        result.errorMessage = "Response spectrum is empty — কমপক্ষে একটা (period, Sa) point দরকার";
        return result;
    }
    if (directionDOF < 0 || directionDOF > 2) {
        result.errorMessage = "directionDOF অবশ্যই 0 (X), 1 (Y), বা 2 (Z) হতে হবে — শুধু translational "
                               "DOF এ ground motion প্রয়োগযোগ্য";
        return result;
    }
    if (dampingRatio < 0.0 || dampingRatio >= 1.0) {
        result.errorMessage = "dampingRatio অবশ্যই [0, 1) সীমায় থাকতে হবে (সাধারণ মান 0.05 = 5%)";
        return result;
    }

    // ------------------------------------------------------------------
    // ধাপ ১: Modal Analysis — natural frequency ও mass-normalized mode
    // shape (φᵢᵀMφᵢ=1) বের করা। CQC-এ ব্যবহারযোগ্য যথেষ্ট mode পেতে,
    // এবং mass participation check এর জন্য, requested numModes ও মোট
    // available mode এর মধ্যে ছোটটা ব্যবহৃত হবে (solveModalAnalysis
    // নিজেই এই clamp করে)।
    // ------------------------------------------------------------------
    ModalAnalysisResult modalResult = solveModalAnalysis(model, numModes);
    if (!modalResult.success) {
        result.errorMessage = "Modal Analysis (RSA এর পূর্বশর্ত) ব্যর্থ হয়েছে: " + modalResult.errorMessage;
        return result;
    }

    const int totalDOF = static_cast<int>(model.nodes.size()) * 6;
    const int numComputedModes = modalResult.numModesComputed;

    // ------------------------------------------------------------------
    // ধাপ ২: Mass matrix — modal participation factor ও effective modal
    // mass বের করতে দরকার (Modal Analysis নিজে শুধু frequency/shape
    // ফেরত দেয়, mass matrix রাখে না — এখানে আবার assemble করা হচ্ছে)।
    // Total structure mass ও এখান থেকে বের করা হবে (mass participation
    // ratio এর denominator)।
    // ------------------------------------------------------------------
    Eigen::SparseMatrix<double> Msparse = assembleGlobalMass(model);
    Eigen::MatrixXd Mdense = Eigen::MatrixXd(Msparse);

    // Influence vector ι — directionDOF অনুযায়ী প্রতিটা node এর সেই
    // translational DOF এ 1, বাকি সব 0।
    Eigen::VectorXd influenceVector = Eigen::VectorXd::Zero(totalDOF);
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        influenceVector(static_cast<int>(i) * 6 + directionDOF) = 1.0;
    }

    // Total structure mass (directionDOF বরাবর) — ιᵀ M ι
    const double totalMass = influenceVector.transpose() * Mdense * influenceVector;
    if (totalMass < 1e-9) {
        result.errorMessage = "Total structure mass (directionDOF বরাবর) প্রায় শূন্য — material density "
                               "সঠিকভাবে সেট করা হয়েছে কিনা যাচাই করুন";
        return result;
    }

    // ------------------------------------------------------------------
    // ধাপ ৩: প্রতিটা mode এর Γᵢ, mᵢ*, Sa(Tᵢ), এবং peak modal response
    // (nodal displacement + element end force, উভয়ই local/global full
    // vector আকারে) বের করা — CQC combination এর আগে সব mode এর জন্য
    // আলাদাভাবে সংরক্ষণ করতে হবে।
    // ------------------------------------------------------------------
    std::vector<double> gammaFactors(numComputedModes);
    std::vector<double> effectiveMasses(numComputedModes);
    std::vector<double> spectralAccelerationsG(numComputedModes);
    std::vector<Eigen::VectorXd> modalPeakDisplacements(numComputedModes); // প্রতিটা: totalDOF-length vector
    std::vector<std::vector<Eigen::VectorXd>> modalPeakElementForces(numComputedModes); // [mode][element] = 12-DOF local force

    const double GRAVITY = 9.81; // m/s² — Sa(g একক) কে m/s² এ রূপান্তরের জন্য

    for (int m = 0; m < numComputedModes; ++m) {
        // fullModeShape পুনর্গঠন — modalResult.modeShapes[m] ইতিমধ্যে
        // per-node 6-DOF vector এর তালিকা (solveModalAnalysis এর
        // আউটপুট shape), এখানে একটা flat totalDOF vector এ concatenate
        // করা হচ্ছে matrix multiplication এর সুবিধার্থে।
        Eigen::VectorXd phi(totalDOF);
        for (size_t i = 0; i < model.nodes.size(); ++i) {
            phi.segment(static_cast<int>(i) * 6, 6) = modalResult.modeShapes[m][i];
        }

        // Γᵢ = φᵢᵀ M ι (φᵢᵀMφᵢ=1 mass-normalization থেকে সরলীকৃত, solver.h এ ব্যাখ্যা করা)
        const double gamma = phi.transpose() * Mdense * influenceVector;
        gammaFactors[m] = gamma;
        effectiveMasses[m] = gamma * gamma;

        const double omega = modalResult.angularFrequenciesRadPerSec[m];
        const double periodSec = (omega > 1e-9) ? (2.0 * M_PI / omega) : 0.0;
        const double saG = interpolateSpectralAcceleration(spectrum, periodSec);
        spectralAccelerationsG[m] = saG;

        // Dᵢ = Γᵢ · Sa / ωᵢ² (Sa কে m/s² এ রূপান্তর করে) — SDOF peak
        // generalized displacement। ω≈0 (rigid body mode, boundary
        // condition ঠিকমতো restrict না করলে তাত্ত্বিকভাবে ঘটতে পারে,
        // যদিও reduced-DOF elimination method এ সাধারণত ঘটে না) হলে
        // division-by-zero এড়াতে peak displacement 0 ধরা হচ্ছে।
        double modalPeakGeneralizedDisp = 0.0;
        if (omega > 1e-6) {
            const double saMetric = saG * GRAVITY; // m/s²
            modalPeakGeneralizedDisp = gamma * saMetric / (omega * omega);
        }

        modalPeakDisplacements[m] = modalPeakGeneralizedDisp * phi;

        // প্রতিটা element এর peak modal end force — mode shape কে একটা
        // "displacement" ধরে (Linear Static এর elementEndForces হিসাবের
        // মতোই পদ্ধতি: f_local = k_local · T · u_element), যেখানে
        // u_element = modalPeakDisplacements[m] থেকে element এর দুই
        // node এর 6-DOF অংশ।
        std::vector<Eigen::VectorXd> elementForcesThisMode(model.elements.size());
        for (size_t e = 0; e < model.elements.size(); ++e) {
            const auto& element = model.elements[e];
            const Node3D& startNode = model.nodes[element.startNodeIndex];
            const Node3D& endNode = model.nodes[element.endNodeIndex];

            const double dx = endNode.x - startNode.x;
            const double dy = endNode.y - startNode.y;
            const double dz = endNode.z - startNode.z;
            const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

            const auto kLocal = getEffectiveLocalStiffness(length, element.section, element.material, element.connectionType);
            const auto T = computeTransformationMatrix(startNode, endNode);

            Eigen::VectorXd uElementGlobal(12);
            uElementGlobal.segment(0, 6) = modalPeakDisplacements[m].segment(element.startNodeIndex * 6, 6);
            uElementGlobal.segment(6, 6) = modalPeakDisplacements[m].segment(element.endNodeIndex * 6, 6);

            Eigen::VectorXd uElementLocal = T * uElementGlobal;
            elementForcesThisMode[e] = kLocal * uElementLocal;
        }
        modalPeakElementForces[m] = std::move(elementForcesThisMode);
    }

    // ------------------------------------------------------------------
    // ধাপ ৪: Mass participation ratio — ব্যবহৃত সব mode এর effective
    // mass যোগফল / total structure mass। code-অনুযায়ী (BNBC/ASCE 7)
    // এটা সাধারণত ≥90% হওয়া আবশ্যক যথেষ্ট mode ব্যবহার করা হয়েছে কিনা
    // যাচাই করতে — caller (Python layer) এই মান দেখে numModes বাড়ানোর
    // পরামর্শ দিতে পারে যদি কম হয়।
    // ------------------------------------------------------------------
    double sumEffectiveMass = 0.0;
    for (double m : effectiveMasses) sumEffectiveMass += m;
    result.totalMassParticipationRatio = sumEffectiveMass / totalMass;

    // ------------------------------------------------------------------
    // ধাপ ৫: CQC Combination — nodal displacement (DOF-by-DOF)।
    // R = sqrt( ΣᵢΣⱼ ρᵢⱼ Rᵢ Rⱼ ) — প্রতিটা DOF-এ আলাদাভাবে প্রয়োগ করা।
    // ------------------------------------------------------------------
    result.nodalDisplacements.resize(model.nodes.size());
    for (int dof = 0; dof < totalDOF; ++dof) {
        double sumSquares = 0.0;
        for (int i = 0; i < numComputedModes; ++i) {
            for (int j = 0; j < numComputedModes; ++j) {
                const double rho = computeCQCCorrelationCoefficient(
                    modalResult.angularFrequenciesRadPerSec[i],
                    modalResult.angularFrequenciesRadPerSec[j],
                    dampingRatio
                );
                sumSquares += rho * modalPeakDisplacements[i](dof) * modalPeakDisplacements[j](dof);
            }
        }
        // sumSquares তাত্ত্বিকভাবে সবসময় ≥0 (CQC formula positive-
        // semidefinite quadratic form), কিন্তু floating-point noise এ
        // সামান্য ঋণাত্মক আসতে পারে খুব ছোট মানে — sqrt(negative) থেকে
        // NaN এড়াতে 0 এ clamp করা।
        const double combined = std::sqrt(std::max(0.0, sumSquares));
        const int nodeIdx = dof / 6;
        const int localDof = dof % 6;
        if (localDof == 0) result.nodalDisplacements[nodeIdx] = Eigen::VectorXd::Zero(6);
        result.nodalDisplacements[nodeIdx](localDof) = combined;
    }

    // ------------------------------------------------------------------
    // ধাপ ৬: CQC Combination — element end force (component-by-component,
    // types.h এর ResponseSpectrumAnalysisResult docstring এ ব্যাখ্যা করা
    // সীমাবদ্ধতা অনুযায়ী)।
    // ------------------------------------------------------------------
    result.elementEndForces.resize(model.elements.size());
    for (size_t e = 0; e < model.elements.size(); ++e) {
        Eigen::VectorXd combinedForce = Eigen::VectorXd::Zero(12);
        for (int component = 0; component < 12; ++component) {
            double sumSquares = 0.0;
            for (int i = 0; i < numComputedModes; ++i) {
                for (int j = 0; j < numComputedModes; ++j) {
                    const double rho = computeCQCCorrelationCoefficient(
                        modalResult.angularFrequenciesRadPerSec[i],
                        modalResult.angularFrequenciesRadPerSec[j],
                        dampingRatio
                    );
                    sumSquares += rho * modalPeakElementForces[i][e](component) * modalPeakElementForces[j][e](component);
                }
            }
            combinedForce(component) = std::sqrt(std::max(0.0, sumSquares));
        }
        result.elementEndForces[e] = combinedForce;
    }

    // ------------------------------------------------------------------
    // ধাপ ৭: Base shear — directionDOF বরাবর, base-supported node গুলোর
    // element end force থেকে যোগ করার বদলে, সহজ ও নির্ভরযোগ্য পদ্ধতি:
    // প্রতিটা mode এর modal base shear = Γᵢ · Sa(Tᵢ) · mᵢ* (SDOF-
    // equivalent base shear সূত্র, Chopra Eq. 13.2.5 এর প্রতিটা mode এ
    // প্রয়োগ), তারপর CQC দিয়ে combine। এটা element-force-summation এর
    // চেয়ে সরল ও সংখ্যাগতভাবে stable (element internal force sign
    // convention এর জটিলতা এড়ায়)।
    // ------------------------------------------------------------------
    std::vector<double> modalBaseShears(numComputedModes);
    for (int m = 0; m < numComputedModes; ++m) {
        const double saMetric = spectralAccelerationsG[m] * GRAVITY; // m/s²
        // Standard modal base shear সূত্র: Vᵢ = mᵢ* · Sa (Chopra Eq.
        // 13.2.5) — Γᵢ ফ্যাক্টর ইতিমধ্যে mᵢ*=Γᵢ² এর মধ্যে অন্তর্ভুক্ত।
        modalBaseShears[m] = effectiveMasses[m] * saMetric;
    }
    double baseShearSumSquares = 0.0;
    for (int i = 0; i < numComputedModes; ++i) {
        for (int j = 0; j < numComputedModes; ++j) {
            const double rho = computeCQCCorrelationCoefficient(
                modalResult.angularFrequenciesRadPerSec[i],
                modalResult.angularFrequenciesRadPerSec[j],
                dampingRatio
            );
            baseShearSumSquares += rho * modalBaseShears[i] * modalBaseShears[j];
        }
    }
    result.baseShear = std::sqrt(std::max(0.0, baseShearSumSquares));

    result.modalParticipationFactors = gammaFactors;
    result.effectiveModalMasses = effectiveMasses;
    result.modalSpectralAccelerations = spectralAccelerationsG;
    result.numModesComputed = numComputedModes;
    result.success = true;
    return result;
}

namespace {

/**
 * একটা element এর current hinge state (কোন প্রান্ত ইতিমধ্যে yielded)
 * অনুযায়ী effective local stiffness বের করে। hinge-yielded প্রান্তে
 * applyEndReleases() দিয়ে bending DOF release করা হয় — connectionType
 * ইতিমধ্যে "pin" হলে (Brace), সেটাও একসাথে ধরা হয় (pin AND yielded
 * উভয় শর্তে release প্রয়োজন, either one release করলেই যথেষ্ট কারণ
 * releaseStart/releaseEnd বুলিয়ান — OR logic)।
 */
Eigen::Matrix<double, 12, 12> getNonlinearEffectiveLocalStiffness(
    double length,
    const FrameElement& element,
    bool startYielded,
    bool endYielded
) {
    const Eigen::Matrix<double, 12, 12> kRigid = computeLocalStiffnessMatrix(length, element.section, element.material);

    const bool releaseStart = (element.connectionType == "pin") || startYielded;
    const bool releaseEnd = (element.connectionType == "pin") || endYielded;

    if (!releaseStart && !releaseEnd) return kRigid;
    return applyEndReleases(kRigid, releaseStart, releaseEnd);
}

} // anonymous namespace

NonlinearStaticAnalysisResult solveNonlinearStatic(
    const AnalysisModel& model,
    int numLoadSteps,
    int maxIterationsPerStep,
    double convergenceTolerance
) {
    NonlinearStaticAnalysisResult result;
    result.success = false;

    if (model.nodes.empty()) {
        result.errorMessage = "Model has no nodes";
        return result;
    }
    if (model.elements.empty()) {
        result.errorMessage = "Model has no elements";
        return result;
    }
    if (model.boundaryConditions.empty()) {
        result.errorMessage = "Model has no boundary conditions — structure is unrestrained (unstable)";
        return result;
    }
    if (model.loads.empty()) {
        result.errorMessage = "Model has no loads — Nonlinear Static Analysis একটা নির্দিষ্ট load "
                               "pattern এর সাপেক্ষে incremental-iterative সমাধান করে (types.h এর "
                               "NonlinearStaticAnalysisResult docstring দেখুন); কোনো load ছাড়া "
                               "yielding ঘটার সুযোগ নেই";
        return result;
    }
    if (numLoadSteps < 1) {
        result.errorMessage = "numLoadSteps অবশ্যই কমপক্ষে 1 হতে হবে";
        return result;
    }
    if (maxIterationsPerStep < 1) {
        result.errorMessage = "maxIterationsPerStep অবশ্যই কমপক্ষে 1 হতে হবে";
        return result;
    }
    if (convergenceTolerance <= 0.0) {
        result.errorMessage = "convergenceTolerance অবশ্যই positive হতে হবে";
        return result;
    }

    const int totalDOF = static_cast<int>(model.nodes.size()) * 6;
    const int numElements = static_cast<int>(model.elements.size());

    // প্রতিটা element এর দুই প্রান্তের hinge state — index [e][0]=start,
    // [e][1]=end। শুধু hingeAtStart/hingeAtEnd=true প্রান্তেই yield চেক
    // হবে, বাকিদের জন্য এই flag সবসময় false থেকে যাবে (ধাপ ঘ এর চেক এ
    // hingeAtStart/hingeAtEnd=false প্রান্ত স্কিপ করা হয়)।
    std::vector<std::array<bool, 2>> yielded(numElements, {false, false});

    const Eigen::VectorXd totalLoad = assembleGlobalLoadVector(model, totalDOF);
    const double loadIncrement = 1.0 / static_cast<double>(numLoadSteps);

    Eigen::VectorXd U = Eigen::VectorXd::Zero(totalDOF); // cumulative displacement, load step জুড়ে বজায় থাকে
    Eigen::VectorXd firstStepDisplacement; // amplification ratio হিসাবের জন্য প্রথম step এর পর সংরক্ষণ

    int totalIterations = 0;
    bool allStepsConverged = true;

    // ------------------------------------------------------------------
    // Load Step Loop — প্রতিটা step এ target load level বৃদ্ধি পায়,
    // cumulative displacement (U) এক step থেকে পরের step এ carry হয়
    // (progressive loading, প্রতিবার শূন্য থেকে শুরু হয় না)।
    // ------------------------------------------------------------------
    for (int step = 1; step <= numLoadSteps; ++step) {
        const Eigen::VectorXd targetLoad = totalLoad * (loadIncrement * static_cast<double>(step));

        bool stepConverged = false;

        // ------------------------------------------------------------------
        // Newton-Raphson Iteration Loop (এই load step এর জন্য) —
        // tangent stiffness rebuild + residual force reduce করা।
        // ------------------------------------------------------------------
        for (int iter = 0; iter < maxIterationsPerStep; ++iter) {
            ++totalIterations;

            // ধাপ ক: বর্তমান hinge state অনুযায়ী tangent stiffness (K_T) assemble
            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(model.elements.size() * 144);

            for (int e = 0; e < numElements; ++e) {
                const auto& element = model.elements[e];
                const Node3D& startNode = model.nodes[element.startNodeIndex];
                const Node3D& endNode = model.nodes[element.endNodeIndex];

                const double dx = endNode.x - startNode.x;
                const double dy = endNode.y - startNode.y;
                const double dz = endNode.z - startNode.z;
                const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

                const auto kLocal = getNonlinearEffectiveLocalStiffness(
                    length, element, yielded[e][0], yielded[e][1]
                );
                const auto T = computeTransformationMatrix(startNode, endNode);
                const Eigen::Matrix<double, 12, 12> kGlobal = T.transpose() * kLocal * T;

                std::array<int, 12> globalIndices{};
                for (int i = 0; i < 6; ++i) {
                    globalIndices[i] = element.startNodeIndex * 6 + i;
                    globalIndices[6 + i] = element.endNodeIndex * 6 + i;
                }
                for (int i = 0; i < 12; ++i) {
                    for (int j = 0; j < 12; ++j) {
                        if (std::abs(kGlobal(i, j)) > 1e-15) {
                            triplets.emplace_back(globalIndices[i], globalIndices[j], kGlobal(i, j));
                        }
                    }
                }
            }

            Eigen::SparseMatrix<double> Ktangent(totalDOF, totalDOF);
            Ktangent.setFromTriplets(triplets.begin(), triplets.end());
            applyBoundaryConditions(Ktangent, model.boundaryConditions, totalDOF);

            // ধাপ খ: Residual force R = F_target - F_internal(U), যেখানে
            // F_internal = K_T · U (tangent stiffness দিয়ে internal
            // resisting force approximate করা হচ্ছে — elastic-perfectly-
            // plastic hinge এর জন্য এই approximation যথেষ্ট নির্ভুল কারণ
            // released DOF তে stiffness contribution ইতিমধ্যে বাদ, তাই
            // internal force সেই DOF এ যোগ হয় না, বাস্তব yielded-moment-
            // capacity ধরে রাখার physical effect approximate হয়)।
            const Eigen::VectorXd internalForce = Ktangent * U;
            Eigen::VectorXd residual = targetLoad - internalForce;

            // Boundary condition এ residual force zero করা (penalty method
            // এর সাথে সামঞ্জস্যপূর্ণ — restrained DOF এ কোনো net force
            // থাকা উচিত না, reaction বাদে, যা এই residual equation এ
            // implicit ভাবে ধরা পড়ে penalty stiffness এর মাধ্যমে)।
            for (const auto& bc : model.boundaryConditions) {
                const int base = bc.nodeIndex * 6;
                if (bc.restrainX) residual(base + 0) = 0.0;
                if (bc.restrainY) residual(base + 1) = 0.0;
                if (bc.restrainZ) residual(base + 2) = 0.0;
                if (bc.restrainRx) residual(base + 3) = 0.0;
                if (bc.restrainRy) residual(base + 4) = 0.0;
                if (bc.restrainRz) residual(base + 5) = 0.0;
            }

            const double residualNorm = residual.norm();
            const double loadNorm = std::max(targetLoad.norm(), 1e-9);

            if (residualNorm / loadNorm < convergenceTolerance) {
                stepConverged = true;
                break;
            }

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(Ktangent);
            if (solver.info() != Eigen::Success) {
                result.errorMessage = "Nonlinear Static Analysis এর tangent stiffness matrix "
                                       "decomposition ব্যর্থ হয়েছে (load step " + std::to_string(step) +
                                       ") — এটা সাধারণত মানে যথেষ্ট hinge yield করে structure একটা "
                                       "mechanism (unstable) হয়ে গেছে এই load level এ। এটাই মূলত "
                                       "structural collapse এর গাণিতিক ইঙ্গিত — সেই load level এ "
                                       "structure এর capacity শেষ হয়ে গেছে।";
                return result;
            }

            const Eigen::VectorXd deltaU = solver.solve(residual);
            if (solver.info() != Eigen::Success) {
                result.errorMessage = "Nonlinear Static Analysis এর linear system solve ব্যর্থ হয়েছে "
                                       "decomposition সফল হওয়ার পরেও (load step " +
                                       std::to_string(step) + ") — numerical ill-conditioning এর ইঙ্গিত";
                return result;
            }

            U += deltaU;

            // ধাপ গ: নতুন U দিয়ে প্রতিটা hinge-assigned প্রান্তের moment
            // recompute করা, yield capacity ছাড়িয়ে গেছে কিনা চেক করা।
            // element এর নিজস্ব (rigid, hinge-release-ছাড়া) local
            // stiffness ব্যবহার করা হচ্ছে এই মোমেন্ট হিসাবে — কারণ moment
            // "capacity ছাড়িয়ে গেছে কিনা" প্রশ্নটা material capacity এর
            // সাপেক্ষে, released DOF এর কৃত্রিম zero-moment এর সাপেক্ষে না
            // (একবার release হয়ে গেলে সেই প্রান্তের moment ধরেই নেওয়া হয়
            // capacity তে "আটকে" আছে, নতুন করে চেক করার দরকার নেই — নিচে
            // yielded[e] true হলে স্কিপ করা হচ্ছে)।
            for (int e = 0; e < numElements; ++e) {
                const auto& element = model.elements[e];
                if (!element.hingeAtStart && !element.hingeAtEnd) continue; // কোনো hinge assign করা নেই এই element এ
                if (element.section.yieldMomentMzKNm <= 0.0) continue; // capacity নির্দিষ্ট করা নেই — চিরকাল elastic

                const Node3D& startNode = model.nodes[element.startNodeIndex];
                const Node3D& endNode = model.nodes[element.endNodeIndex];
                const double dx = endNode.x - startNode.x;
                const double dy = endNode.y - startNode.y;
                const double dz = endNode.z - startNode.z;
                const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

                const auto kRigidLocal = computeLocalStiffnessMatrix(length, element.section, element.material);
                const auto T = computeTransformationMatrix(startNode, endNode);

                Eigen::VectorXd uElementGlobal(12);
                uElementGlobal.segment(0, 6) = U.segment(element.startNodeIndex * 6, 6);
                uElementGlobal.segment(6, 6) = U.segment(element.endNodeIndex * 6, 6);
                const Eigen::VectorXd uElementLocal = T * uElementGlobal;
                const Eigen::VectorXd fElementLocal = kRigidLocal * uElementLocal;

                // local index 5 = start-node Mz, local index 11 = end-node Mz
                // (stiffness.h এর DOF ordering: [u,v,w,rx,ry,rz] প্রতি node)
                if (element.hingeAtStart && !yielded[e][0] && std::abs(fElementLocal(5)) > element.section.yieldMomentMzKNm) {
                    yielded[e][0] = true;
                }
                if (element.hingeAtEnd && !yielded[e][1] && std::abs(fElementLocal(11)) > element.section.yieldMomentMzKNm) {
                    yielded[e][1] = true;
                }
            }
        }

        if (!stepConverged) {
            allStepsConverged = false;
            // non-convergence হলেও এগিয়ে যাওয়া হচ্ছে (crash না করে) —
            // এই load step এর সেরা প্রাপ্ত U ব্যবহার করে পরবর্তী step এ
            // এগোনো, চূড়ান্ত ফলাফলে converged=false দিয়ে caller কে
            // জানানো (types.h এর NonlinearStaticAnalysisResult docstring
            // অনুযায়ী)।
        }

        if (step == 1) {
            firstStepDisplacement = U;
        }
    }

    // ------------------------------------------------------------------
    // চূড়ান্ত ফলাফল সংকলন — displacement, element end force, hinge state
    // ------------------------------------------------------------------
    result.nodalDisplacements.resize(model.nodes.size());
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        result.nodalDisplacements[i] = U.segment(static_cast<int>(i) * 6, 6);
    }

    result.elementEndForces.resize(numElements);
    for (int e = 0; e < numElements; ++e) {
        const auto& element = model.elements[e];
        const Node3D& startNode = model.nodes[element.startNodeIndex];
        const Node3D& endNode = model.nodes[element.endNodeIndex];
        const double dx = endNode.x - startNode.x;
        const double dy = endNode.y - startNode.y;
        const double dz = endNode.z - startNode.z;
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

        // চূড়ান্ত element end force — dual পদ্ধতি ব্যবহার করা হচ্ছে:
        //   ১. Non-yielded DOF: hinge-aware effective (released যেখানে
        //      প্রযোজ্য) stiffness দিয়ে হিসাব — Linear Static এর মতোই।
        //   ২. Yielded moment DOF (Mz, index 5/11): effective stiffness
        //      সেই DOF এ 0 (released বলে), তাই সরাসরি সেখান থেকে force
        //      পড়লে ভুলভাবে 0 আসবে (released মানে "structure locally ঐ
        //      DOF এ কোনো additional moment resist করে না", কিন্তু hinge
        //      নিজে তখনও তার yield capacity তে moment ধরে রাখে —
        //      elastic-perfectly-plastic সংজ্ঞা অনুযায়ী)। তাই yielded
        //      প্রান্তে moment = ±yieldCapacity (rigid-stiffness হিসাবের
        //      sign অনুসরণ করে) সরাসরি বসানো হচ্ছে।
        const auto kEffective = getNonlinearEffectiveLocalStiffness(length, element, yielded[e][0], yielded[e][1]);
        const auto T = computeTransformationMatrix(startNode, endNode);

        Eigen::VectorXd uElementGlobal(12);
        uElementGlobal.segment(0, 6) = U.segment(element.startNodeIndex * 6, 6);
        uElementGlobal.segment(6, 6) = U.segment(element.endNodeIndex * 6, 6);
        const Eigen::VectorXd uElementLocal = T * uElementGlobal;
        Eigen::VectorXd finalForce = kEffective * uElementLocal;

        if (yielded[e][0] || yielded[e][1]) {
            // Sign নির্ধারণের জন্য rigid (release-ছাড়া) stiffness দিয়ে
            // moment recompute করা — এই sign yield মুহূর্তে যে দিকে
            // moment demand ছিল তা ধরে রাখে (elastic-perfectly-plastic
            // এ moment magnitude capacity তে আটকায়, কিন্তু sign/direction
            // অপরিবর্তিত থাকে)।
            const auto kRigidLocal = computeLocalStiffnessMatrix(length, element.section, element.material);
            const Eigen::VectorXd rigidForce = kRigidLocal * uElementLocal;

            if (yielded[e][0]) {
                const double sign = (rigidForce(5) >= 0.0) ? 1.0 : -1.0;
                finalForce(5) = sign * element.section.yieldMomentMzKNm;
            }
            if (yielded[e][1]) {
                const double sign = (rigidForce(11) >= 0.0) ? 1.0 : -1.0;
                finalForce(11) = sign * element.section.yieldMomentMzKNm;
            }
        }

        result.elementEndForces[e] = finalForce;

        // hinge state রিপোর্ট করা (শুধু hinge-assigned প্রান্ত)
        if (element.hingeAtStart) {
            PlasticHingeState h;
            h.elementIndex = e;
            h.isAtStartNode = true;
            h.yielded = yielded[e][0];
            h.finalMomentKNm = result.elementEndForces[e](5);
            result.hingeStates.push_back(h);
        }
        if (element.hingeAtEnd) {
            PlasticHingeState h;
            h.elementIndex = e;
            h.isAtStartNode = false;
            h.yielded = yielded[e][1];
            h.finalMomentKNm = result.elementEndForces[e](11);
            result.hingeStates.push_back(h);
        }
    }

    // Displacement amplification ratio — চূড়ান্ত (সম্পূর্ণ load) বনাম
    // প্রথম load step (near-elastic, যেহেতু yielding সাধারণত পরবর্তী
    // step এ শুরু হয়) displacement এর সবচেয়ে বড় অনুপাত — nonlinear
    // effect এর overall মাত্রার একটা দ্রুত ইঙ্গিত (P-Delta এর
    // maxDisplacementAmplificationRatio এর সমতুল্য ধারণা)।
    double maxRatio = 1.0;
    if (firstStepDisplacement.size() == totalDOF) {
        for (int i = 0; i < totalDOF; ++i) {
            const double firstVal = firstStepDisplacement(i) * static_cast<double>(numLoadSteps); // extrapolate to full load (elastic assumption)
            const double finalVal = U(i);
            if (std::abs(firstVal) > 1e-9) {
                const double ratio = std::abs(finalVal / firstVal);
                if (std::isfinite(ratio) && ratio > maxRatio) {
                    maxRatio = ratio;
                }
            }
        }
    }
    result.maxDisplacementAmplificationRatio = maxRatio;

    result.totalLoadSteps = numLoadSteps;
    result.totalNewtonIterations = totalIterations;
    result.converged = allStepsConverged;
    result.success = true;
    return result;
}

PushoverAnalysisResult solvePushover(
    const AnalysisModel& model,
    int controlNodeIndex,
    int controlDOF,
    double targetControlDisplacementM,
    double loadStepIncrement,
    int maxPushSteps,
    int maxIterationsPerStep,
    double convergenceTolerance
) {
    PushoverAnalysisResult result;
    result.success = false;

    if (model.nodes.empty()) {
        result.errorMessage = "Model has no nodes";
        return result;
    }
    if (model.elements.empty()) {
        result.errorMessage = "Model has no elements";
        return result;
    }
    if (model.boundaryConditions.empty()) {
        result.errorMessage = "Model has no boundary conditions — structure is unrestrained (unstable)";
        return result;
    }
    if (model.loads.empty()) {
        result.errorMessage = "Model has no loads — Pushover একটা fixed lateral load pattern push করে "
                               "(solver.h এর solvePushover() docstring দেখুন); কোনো load ছাড়া push "
                               "করার কিছু নেই";
        return result;
    }
    if (controlNodeIndex < 0 || controlNodeIndex >= static_cast<int>(model.nodes.size())) {
        result.errorMessage = "controlNodeIndex model এর node সংখ্যার সীমার বাইরে";
        return result;
    }
    if (controlDOF < 0 || controlDOF > 2) {
        result.errorMessage = "controlDOF অবশ্যই 0 (ux), 1 (uy), বা 2 (uz) হতে হবে";
        return result;
    }
    if (targetControlDisplacementM <= 0.0) {
        result.errorMessage = "targetControlDisplacementM অবশ্যই positive হতে হবে (push direction "
                               "নির্বিশেষে magnitude হিসেবে ব্যবহৃত হয়)";
        return result;
    }
    if (loadStepIncrement <= 0.0 || loadStepIncrement > 1.0) {
        result.errorMessage = "loadStepIncrement অবশ্যই (0, 1] সীমায় থাকতে হবে (full load pattern এর fraction)";
        return result;
    }
    if (maxPushSteps < 1) {
        result.errorMessage = "maxPushSteps অবশ্যই কমপক্ষে 1 হতে হবে";
        return result;
    }
    if (maxIterationsPerStep < 1) {
        result.errorMessage = "maxIterationsPerStep অবশ্যই কমপক্ষে 1 হতে হবে";
        return result;
    }
    if (convergenceTolerance <= 0.0) {
        result.errorMessage = "convergenceTolerance অবশ্যই positive হতে হবে";
        return result;
    }

    const int totalDOF = static_cast<int>(model.nodes.size()) * 6;
    const int numElements = static_cast<int>(model.elements.size());
    const int controlGlobalDOF = controlNodeIndex * 6 + controlDOF;

    std::vector<std::array<bool, 2>> yielded(numElements, {false, false});

    const Eigen::VectorXd fullLoad = assembleGlobalLoadVector(model, totalDOF);
    const double fullLoadNorm = fullLoad.norm();
    if (fullLoadNorm < 1e-9) {
        result.errorMessage = "প্রয়োগকৃত load এর magnitude প্রায় শূন্য — push করার কিছু নেই";
        return result;
    }

    Eigen::VectorXd U = Eigen::VectorXd::Zero(totalDOF);
    result.capacityCurve.push_back(PushoverCurvePoint{0.0, 0.0, 0}); // push শুরুর বিন্দু (origin)

    int totalIterations = 0;
    double currentLoadFactor = 0.0; // 0 থেকে 1+ পর্যন্ত (1.0 = পূর্ণ model.loads pattern একবার প্রয়োগ)

    // ------------------------------------------------------------------
    // Push Step Loop — প্রতিটা step এ load factor loadStepIncrement
    // দিয়ে বাড়ে, target control displacement এ পৌঁছানো বা tangent
    // stiffness singular (collapse) না হওয়া পর্যন্ত।
    // ------------------------------------------------------------------
    for (int step = 1; step <= maxPushSteps; ++step) {
        const double nextLoadFactor = currentLoadFactor + loadStepIncrement;
        const Eigen::VectorXd targetLoad = fullLoad * nextLoadFactor;

        bool stepConverged = false;
        Eigen::VectorXd Ustep = U; // এই step এর trial displacement, শুধু convergence এ commit হবে

        for (int iter = 0; iter < maxIterationsPerStep; ++iter) {
            ++totalIterations;

            std::vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(model.elements.size() * 144);

            for (int e = 0; e < numElements; ++e) {
                const auto& element = model.elements[e];
                const Node3D& startNode = model.nodes[element.startNodeIndex];
                const Node3D& endNode = model.nodes[element.endNodeIndex];

                const double dx = endNode.x - startNode.x;
                const double dy = endNode.y - startNode.y;
                const double dz = endNode.z - startNode.z;
                const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

                const auto kLocal = getNonlinearEffectiveLocalStiffness(
                    length, element, yielded[e][0], yielded[e][1]
                );
                const auto T = computeTransformationMatrix(startNode, endNode);
                const Eigen::Matrix<double, 12, 12> kGlobal = T.transpose() * kLocal * T;

                std::array<int, 12> globalIndices{};
                for (int i = 0; i < 6; ++i) {
                    globalIndices[i] = element.startNodeIndex * 6 + i;
                    globalIndices[6 + i] = element.endNodeIndex * 6 + i;
                }
                for (int i = 0; i < 12; ++i) {
                    for (int j = 0; j < 12; ++j) {
                        if (std::abs(kGlobal(i, j)) > 1e-15) {
                            triplets.emplace_back(globalIndices[i], globalIndices[j], kGlobal(i, j));
                        }
                    }
                }
            }

            Eigen::SparseMatrix<double> Ktangent(totalDOF, totalDOF);
            Ktangent.setFromTriplets(triplets.begin(), triplets.end());
            applyBoundaryConditions(Ktangent, model.boundaryConditions, totalDOF);

            const Eigen::VectorXd internalForce = Ktangent * Ustep;
            Eigen::VectorXd residual = targetLoad - internalForce;

            for (const auto& bc : model.boundaryConditions) {
                const int base = bc.nodeIndex * 6;
                if (bc.restrainX) residual(base + 0) = 0.0;
                if (bc.restrainY) residual(base + 1) = 0.0;
                if (bc.restrainZ) residual(base + 2) = 0.0;
                if (bc.restrainRx) residual(base + 3) = 0.0;
                if (bc.restrainRy) residual(base + 4) = 0.0;
                if (bc.restrainRz) residual(base + 5) = 0.0;
            }

            const double residualNorm = residual.norm();
            const double loadNorm = std::max(targetLoad.norm(), 1e-9);

            if (residualNorm / loadNorm < convergenceTolerance) {
                stepConverged = true;
                break;
            }

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(Ktangent);
            if (solver.info() != Eigen::Success) {
                // Tangent stiffness singular — structure এই load level এ
                // আর কোনো additional load নিতে পারছে না (mechanism/
                // collapse)। এটা crash/error না, বরং pushover এর একটা
                // স্বাভাবিক, অর্থপূর্ণ থামার শর্ত — capacity curve এর
                // শেষ বিন্দু হলো structure এর ultimate capacity।
                result.structureCollapsed = true;
                break;
            }

            const Eigen::VectorXd deltaU = solver.solve(residual);
            if (solver.info() != Eigen::Success) {
                result.structureCollapsed = true;
                break;
            }

            Ustep += deltaU;

            for (int e = 0; e < numElements; ++e) {
                const auto& element = model.elements[e];
                if (!element.hingeAtStart && !element.hingeAtEnd) continue;
                if (element.section.yieldMomentMzKNm <= 0.0) continue;

                const Node3D& startNode = model.nodes[element.startNodeIndex];
                const Node3D& endNode = model.nodes[element.endNodeIndex];
                const double dx = endNode.x - startNode.x;
                const double dy = endNode.y - startNode.y;
                const double dz = endNode.z - startNode.z;
                const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

                const auto kRigidLocal = computeLocalStiffnessMatrix(length, element.section, element.material);
                const auto T = computeTransformationMatrix(startNode, endNode);

                Eigen::VectorXd uElementGlobal(12);
                uElementGlobal.segment(0, 6) = Ustep.segment(element.startNodeIndex * 6, 6);
                uElementGlobal.segment(6, 6) = Ustep.segment(element.endNodeIndex * 6, 6);
                const Eigen::VectorXd uElementLocal = T * uElementGlobal;
                const Eigen::VectorXd fElementLocal = kRigidLocal * uElementLocal;

                if (element.hingeAtStart && !yielded[e][0] && std::abs(fElementLocal(5)) > element.section.yieldMomentMzKNm) {
                    yielded[e][0] = true;
                }
                if (element.hingeAtEnd && !yielded[e][1] && std::abs(fElementLocal(11)) > element.section.yieldMomentMzKNm) {
                    yielded[e][1] = true;
                }
            }
        }

        if (result.structureCollapsed) {
            // এই step commit হয় না (Ustep বাতিল) — collapse-এর ঠিক আগের
            // শেষ successful step ই চূড়ান্ত ফলাফল হিসেবে থেকে যাবে।
            break;
        }

        if (!stepConverged) {
            // এই নির্দিষ্ট step এ maxIterationsPerStep এর মধ্যে converge
            // করেনি (কিন্তু matrix singular হয়নি) — এটাও practically
            // একটা near-collapse অবস্থা, push বন্ধ করা হচ্ছে এখানে,
            // এই step বাতিল করে আগের commit করা অবস্থা রাখা হচ্ছে।
            result.structureCollapsed = true;
            break;
        }

        // Step সফল — commit করা এবং capacity curve এ point যোগ করা
        U = Ustep;
        currentLoadFactor = nextLoadFactor;

        const double baseShearKN = fullLoadNorm * currentLoadFactor; // equilibrium: applied lateral load sum = base reaction sum
        const double controlDisp = U(controlGlobalDOF);

        int numYielded = 0;
        for (const auto& pair : yielded) {
            if (pair[0]) ++numYielded;
            if (pair[1]) ++numYielded;
        }

        result.capacityCurve.push_back(PushoverCurvePoint{baseShearKN, controlDisp, numYielded});
        result.totalPushSteps = step;

        if (std::abs(controlDisp) >= targetControlDisplacementM) {
            result.reachedTargetDisplacement = true;
            break;
        }
    }

    // ------------------------------------------------------------------
    // চূড়ান্ত ফলাফল সংকলন — শেষ successfully-committed U দিয়ে
    // ------------------------------------------------------------------
    result.finalNodalDisplacements.resize(model.nodes.size());
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        result.finalNodalDisplacements[i] = U.segment(static_cast<int>(i) * 6, 6);
    }

    result.finalElementEndForces.resize(numElements);
    for (int e = 0; e < numElements; ++e) {
        const auto& element = model.elements[e];
        const Node3D& startNode = model.nodes[element.startNodeIndex];
        const Node3D& endNode = model.nodes[element.endNodeIndex];
        const double dx = endNode.x - startNode.x;
        const double dy = endNode.y - startNode.y;
        const double dz = endNode.z - startNode.z;
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

        const auto kEffective = getNonlinearEffectiveLocalStiffness(length, element, yielded[e][0], yielded[e][1]);
        const auto T = computeTransformationMatrix(startNode, endNode);

        Eigen::VectorXd uElementGlobal(12);
        uElementGlobal.segment(0, 6) = U.segment(element.startNodeIndex * 6, 6);
        uElementGlobal.segment(6, 6) = U.segment(element.endNodeIndex * 6, 6);
        const Eigen::VectorXd uElementLocal = T * uElementGlobal;
        Eigen::VectorXd finalForce = kEffective * uElementLocal;

        if (yielded[e][0] || yielded[e][1]) {
            const auto kRigidLocal = computeLocalStiffnessMatrix(length, element.section, element.material);
            const Eigen::VectorXd rigidForce = kRigidLocal * uElementLocal;

            if (yielded[e][0]) {
                const double sign = (rigidForce(5) >= 0.0) ? 1.0 : -1.0;
                finalForce(5) = sign * element.section.yieldMomentMzKNm;
            }
            if (yielded[e][1]) {
                const double sign = (rigidForce(11) >= 0.0) ? 1.0 : -1.0;
                finalForce(11) = sign * element.section.yieldMomentMzKNm;
            }
        }

        result.finalElementEndForces[e] = finalForce;

        if (element.hingeAtStart) {
            PlasticHingeState h;
            h.elementIndex = e;
            h.isAtStartNode = true;
            h.yielded = yielded[e][0];
            h.finalMomentKNm = result.finalElementEndForces[e](5);
            result.finalHingeStates.push_back(h);
        }
        if (element.hingeAtEnd) {
            PlasticHingeState h;
            h.elementIndex = e;
            h.isAtStartNode = false;
            h.yielded = yielded[e][1];
            h.finalMomentKNm = result.finalElementEndForces[e](11);
            result.finalHingeStates.push_back(h);
        }
    }

    result.totalNewtonIterations = totalIterations;
    result.success = true; // collapse/non-convergence-এ থামলেও success=true — সেটা একটা বৈধ pushover ফলাফল (ultimate capacity পাওয়া গেছে), caller structureCollapsed/reachedTargetDisplacement দেখে বুঝবে
    return result;
}

} // namespace civilos
