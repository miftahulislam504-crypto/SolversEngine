#include "shell.h"
#include <cmath>
#include <stdexcept>

namespace civilos {

namespace {

// ============================================================
// সাধারণ (Common) helper — 2D shape function ও derivative, membrane
// ও plate bending উভয় অংশেই ব্যবহৃত হয় (isoparametric mapping)।
// ============================================================

/**
 * Bilinear শেপ ফাংশন (4-node quad), natural coordinate (ξ, η) এ,
 * প্রতিটা -1..1 রেঞ্জে। N[i](ξ,η) — node i এর শেপ ফাংশন মান।
 * কনভেনশন: node 0=(-1,-1), 1=(1,-1), 2=(1,1), 3=(-1,1) (standard Q4
 * counter-clockwise numbering, types.h এর ShellElement এ node ordering
 * এর সাথে সঙ্গতিপূর্ণ)।
 */
std::array<double, 4> shapeFunctions(double xi, double eta) {
    return {
        0.25 * (1 - xi) * (1 - eta),
        0.25 * (1 + xi) * (1 - eta),
        0.25 * (1 + xi) * (1 + eta),
        0.25 * (1 - xi) * (1 + eta),
    };
}

/** শেপ ফাংশনের ξ ও η এর সাপেক্ষে derivative — natural coordinate এ। */
std::array<std::array<double, 2>, 4> shapeFunctionDerivatives(double xi, double eta) {
    // dN[i] = {dN_i/dξ, dN_i/dη}
    return {{
        {-0.25 * (1 - eta), -0.25 * (1 - xi)},
        { 0.25 * (1 - eta), -0.25 * (1 + xi)},
        { 0.25 * (1 + eta),  0.25 * (1 + xi)},
        {-0.25 * (1 + eta),  0.25 * (1 - xi)},
    }};
}

/**
 * Jacobian matrix (2x2) — natural (ξ,η) থেকে local physical (x,y)
 * coordinate এ mapping এর derivative, এবং তার determinant (element
 * area scaling factor, integration এ dx·dy = det(J)·dξ·dη)।
 */
struct JacobianResult {
    Eigen::Matrix2d J;
    double detJ;
    Eigen::Matrix2d Jinv;
};

JacobianResult computeJacobian(
    const std::array<Eigen::Vector3d, 4>& corners,
    double xi, double eta
) {
    auto dN = shapeFunctionDerivatives(xi, eta);

    Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
    for (int i = 0; i < 4; ++i) {
        J(0, 0) += dN[i][0] * corners[i].x(); // dx/dξ
        J(0, 1) += dN[i][0] * corners[i].y(); // dy/dξ
        J(1, 0) += dN[i][1] * corners[i].x(); // dx/dη
        J(1, 1) += dN[i][1] * corners[i].y(); // dy/dη
    }

    JacobianResult result;
    result.J = J;
    result.detJ = J.determinant();

    if (std::abs(result.detJ) < 1e-12) {
        throw std::invalid_argument(
            "Shell element Jacobian determinant প্রায় শূন্য — element degenerate (zero area, বা "
            "corner গুলো collinear/duplicate) হতে পারে, mesh geometry যাচাই করুন");
    }

    result.Jinv = J.inverse();
    return result;
}

// 2x2 Gauss quadrature points ও weights (membrane integration এর জন্য) — standard,
// প্রতিটা weight=1.0 (2-point Gauss rule প্রতি dimension এ, ±1/√3 এ)
constexpr double GAUSS_2PT = 0.5773502691896258; // 1/√3

const std::array<std::array<double, 2>, 4> GAUSS_2X2_POINTS = {{
    {-GAUSS_2PT, -GAUSS_2PT},
    { GAUSS_2PT, -GAUSS_2PT},
    { GAUSS_2PT,  GAUSS_2PT},
    {-GAUSS_2PT,  GAUSS_2PT},
}};

// ============================================================
// Membrane (Plane Stress) Stiffness — Bilinear Q4
// ============================================================

/**
 * Plane-stress constitutive matrix (3x3) — isotropic material,
 * standard সূত্র (Cook et al., Ch. 6): D = E/(1-ν²) * [[1,ν,0],[ν,1,0],[0,0,(1-ν)/2]]
 */
Eigen::Matrix3d planeStressConstitutiveMatrix(double E, double nu) {
    Eigen::Matrix3d D = Eigen::Matrix3d::Zero();
    const double factor = E / (1.0 - nu * nu);
    D(0, 0) = factor;
    D(0, 1) = factor * nu;
    D(1, 0) = factor * nu;
    D(1, 1) = factor;
    D(2, 2) = factor * (1.0 - nu) / 2.0;
    return D;
}

/**
 * Membrane stiffness (8x8, শুধু ux,uy DOF প্রতি node — drilling rz
 * এখানে না, সেটা আলাদাভাবে যোগ হয় পরে) — 2x2 Gauss quadrature দিয়ে
 * numerical integration: K = ∫∫ B^T D B t dξdη (t=thickness, B=strain-
 * displacement matrix, ব্যাখ্যা নিচে লুপের ভিতরে)।
 */
Eigen::Matrix<double, 8, 8> computeMembraneStiffness(
    const std::array<Eigen::Vector3d, 4>& corners,
    double thickness,
    double E,
    double nu
) {
    Eigen::Matrix<double, 8, 8> K = Eigen::Matrix<double, 8, 8>::Zero();
    Eigen::Matrix3d D = planeStressConstitutiveMatrix(E, nu);

    for (const auto& gp : GAUSS_2X2_POINTS) {
        const double xi = gp[0];
        const double eta = gp[1];

        auto jac = computeJacobian(corners, xi, eta);
        auto dN = shapeFunctionDerivatives(xi, eta);

        // dN/dx, dN/dy (physical derivative) = Jinv * dN/dξdη
        std::array<double, 4> dNdx{}, dNdy{};
        for (int i = 0; i < 4; ++i) {
            Eigen::Vector2d dNi_natural(dN[i][0], dN[i][1]);
            Eigen::Vector2d dNi_physical = jac.Jinv * dNi_natural;
            dNdx[i] = dNi_physical(0);
            dNdy[i] = dNi_physical(1);
        }

        // Strain-displacement matrix B (3x8) — প্রতিটা node এর (ux,uy)
        // থেকে strain (εxx, εyy, γxy) এ ম্যাপিং। Standard plane-stress
        // B matrix গঠন (Cook et al. Ch. 6, eq. 6.2-x):
        Eigen::Matrix<double, 3, 8> B = Eigen::Matrix<double, 3, 8>::Zero();
        for (int i = 0; i < 4; ++i) {
            B(0, 2 * i)     = dNdx[i];  // εxx = ∂u/∂x
            B(1, 2 * i + 1) = dNdy[i];  // εyy = ∂v/∂y
            B(2, 2 * i)     = dNdy[i];  // γxy = ∂u/∂y + ∂v/∂x
            B(2, 2 * i + 1) = dNdx[i];
        }

        // Gauss weight (2x2 rule, প্রতিটা point এ weight=1.0×1.0=1.0)
        const double weight = 1.0;
        K += B.transpose() * D * B * thickness * jac.detJ * weight;
    }

    return K;
}

// ============================================================
// Drilling DOF (rz) Stabilization — Allman-type penalty stiffness
// ============================================================

/**
 * প্রতিটা node এর rz (in-plane/drilling rotation) DOF এ একটা ছোট
 * penalty stiffness বসানো — bilinear membrane element প্রকৃতপক্ষে এই
 * DOF ব্যবহার করে না (কোনো natural stiffness নেই), কিন্তু আমাদের
 * global system এ প্রতিটা node এ ৬-DOF দরকার (frame element এর সাথে
 * সঙ্গতিপূর্ণ, এবং সব-coplanar shell mesh এ rz DOF সম্পূর্ণ unconstrained
 * থাকলে stiffness matrix singular হয়ে যাবে সেই DOF এ)।
 *
 * পদ্ধতি: একটা ছোট, physically-negligible কিন্তু numerically-non-zero
 * stiffness যোগ করা প্রতিটা node এর rz এ (diagonal penalty, off-
 * diagonal coupling ছাড়া, simplification হিসেবে — পূর্ণ Allman
 * drilling formulation এ rz প্রকৃত membrane strain energy এর সাথে
 * coupled থাকে, কিন্তু সেটা উল্লেখযোগ্যভাবে বেশি জটিল)। Magnitude
 * membrane stiffness এর গড় মানের একটা ছোট ভগ্নাংশ (~1e-3 থেকে 1e-4
 * factor) হিসেবে নেওয়া হয়েছে — যথেষ্ট বড় যাতে matrix non-singular
 * থাকে, যথেষ্ট ছোট যাতে প্রকৃত membrane/bending behavior এ উল্লেখযোগ্য
 * প্রভাব না ফেলে।
 *
 * এই সরলীকরণের ফলে: rz DOF এ প্রয়োগ করা সরাসরি moment বাস্তবসম্মতভাবে
 * ধরা পড়বে না (drilling moment এর জন্য কোনো real stiffness নেই) —
 * এটা একটা known সীমাবদ্ধতা, shell.h এর computeShellLocalStiffness
 * docstring এ নোট করা প্রয়োজন যদি এই ফাংশন থেকে আলাদা করে ডাকা হয়।
 */
double drillingPenaltyStiffness(double membraneStiffnessMagnitude) {
    constexpr double PENALTY_FACTOR = 1e-3;
    return membraneStiffnessMagnitude * PENALTY_FACTOR;
}

// ============================================================
// Plate Bending Stiffness — Mindlin-Reissner, Selective Reduced Integration
// ============================================================

/**
 * Bending constitutive matrix (3x3) — isotropic thin/moderately-thick
 * plate, standard সূত্র (Cook et al. Ch. 13, eq. 13.2-x):
 *   Db = (E*t³ / (12*(1-ν²))) * [[1,ν,0],[ν,1,0],[0,0,(1-ν)/2]]
 * plane-stress constitutive matrix (planeStressConstitutiveMatrix)
 * এর গঠন একই, শুধু t³/12 factor (bending rigidity) দিয়ে গুণ করা।
 */
Eigen::Matrix3d bendingConstitutiveMatrix(double E, double nu, double thickness) {
    const double factor = (E * thickness * thickness * thickness) / (12.0 * (1.0 - nu * nu));
    Eigen::Matrix3d Db = Eigen::Matrix3d::Zero();
    Db(0, 0) = factor;
    Db(0, 1) = factor * nu;
    Db(1, 0) = factor * nu;
    Db(1, 1) = factor;
    Db(2, 2) = factor * (1.0 - nu) / 2.0;
    return Db;
}

/**
 * Transverse shear constitutive matrix (2x2) — Mindlin plate theory,
 * standard সূত্র: Ds = k * G * t * I₂, যেখানে k = shear correction
 * factor (5/6, thin homogeneous plate এর জন্য প্রচলিত মান — Reissner
 * এর মূল derivation থেকে), G = shear modulus = E/(2(1+ν)) (isotropic
 * material)।
 */
Eigen::Matrix2d shearConstitutiveMatrix(double E, double nu, double thickness) {
    constexpr double SHEAR_CORRECTION_FACTOR = 5.0 / 6.0;
    const double G = E / (2.0 * (1.0 + nu));
    const double factor = SHEAR_CORRECTION_FACTOR * G * thickness;
    Eigen::Matrix2d Ds = Eigen::Matrix2d::Zero();
    Ds(0, 0) = factor;
    Ds(1, 1) = factor;
    return Ds;
}

/**
 * Mindlin plate bending element stiffness (12x12 — 3 DOF [w,θx,θy]
 * প্রতি node, 4 node)।
 *
 * DOF sign convention (এই ফাংশনের নিজস্ব local convention, shell.h এর
 * ShellElement docstring এ ভবিষ্যতে reference করার জন্য এখানে
 * স্পষ্টভাবে লেখা হলো):
 *   w = local z (normal) দিকে deflection
 *   θx = local x-axis এর চারপাশে rotation (positive: right-hand rule,
 *        θx বাড়লে +y দিকে w কমে — অর্থাৎ κy = -∂θx/∂y একটা প্রচলিত
 *        Mindlin-plate sign convention)
 *   θy = local y-axis এর চারপাশে rotation (κx = ∂θy/∂x)
 *
 * Bending strain (curvature): κx=∂θy/∂x, κy=-∂θx/∂y, κxy=∂θy/∂y-∂θx/∂x
 * Shear strain: γxz=∂w/∂x-θy... প্রকৃতপক্ষে convention ভেদে চিহ্ন
 * পরিবর্তিত হয়, তাই নিচে সরাসরি Cook et al. Ch. 13 এর eq. 13.3
 * অনুসরণ করা হয়েছে হুবহু।
 *
 * Selective Reduced Integration: bending অংশ 2×2 full Gauss quadrature
 * দিয়ে integrate করা হয়, কিন্তু shear অংশ 1×1 reduced Gauss (কেন্দ্র
 * বিন্দু) দিয়ে — এটা "shear locking" (thin plate limit এ bilinear
 * element এর একটা সুপরিচিত numerical pathology, যেখানে full
 * integration কৃত্রিমভাবে অতিরিক্ত-stiff ফলাফল দেয়) এড়ানোর একটা
 * প্রতিষ্ঠিত, ব্যাপকভাবে ব্যবহৃত কৌশল (Hughes, "The Finite Element
 * Method", Ch. 5; Zienkiewicz & Taylor)।
 */
Eigen::Matrix<double, 12, 12> computePlateBendingStiffness(
    const std::array<Eigen::Vector3d, 4>& corners,
    double thickness,
    double E,
    double nu
) {
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();

    Eigen::Matrix3d Db = bendingConstitutiveMatrix(E, nu, thickness);
    Eigen::Matrix2d Ds = shearConstitutiveMatrix(E, nu, thickness);

    // --- Bending part: full 2x2 Gauss integration ---
    for (const auto& gp : GAUSS_2X2_POINTS) {
        const double xi = gp[0];
        const double eta = gp[1];

        auto jac = computeJacobian(corners, xi, eta);
        auto dN = shapeFunctionDerivatives(xi, eta);

        std::array<double, 4> dNdx{}, dNdy{};
        for (int i = 0; i < 4; ++i) {
            Eigen::Vector2d dNi_natural(dN[i][0], dN[i][1]);
            Eigen::Vector2d dNi_physical = jac.Jinv * dNi_natural;
            dNdx[i] = dNi_physical(0);
            dNdy[i] = dNi_physical(1);
        }

        // Bending B matrix (3x12) — [w,θx,θy] প্রতি node থেকে
        // curvature [κx,κy,κxy] এ ম্যাপিং। w কলামে শূন্য (bending
        // curvature সরাসরি w থেকে আসে না Mindlin theory তে, শুধু
        // rotation থেকে — এটাই Kirchhoff থেকে মূল পার্থক্য)।
        Eigen::Matrix<double, 3, 12> Bb = Eigen::Matrix<double, 3, 12>::Zero();
        for (int i = 0; i < 4; ++i) {
            const int col = 3 * i;
            // κx = ∂θy/∂x
            Bb(0, col + 2) = dNdx[i];
            // κy = -∂θx/∂y
            Bb(1, col + 1) = -dNdy[i];
            // κxy = ∂θy/∂y - ∂θx/∂x
            Bb(2, col + 1) = -dNdx[i];
            Bb(2, col + 2) = dNdy[i];
        }

        const double weight = 1.0;
        K += Bb.transpose() * Db * Bb * jac.detJ * weight;
    }

    // --- Shear part: reduced 1x1 Gauss integration (center point, ξ=η=0) ---
    {
        const double xi = 0.0;
        const double eta = 0.0;

        auto jac = computeJacobian(corners, xi, eta);
        auto N = shapeFunctions(xi, eta);
        auto dN = shapeFunctionDerivatives(xi, eta);

        std::array<double, 4> dNdx{}, dNdy{};
        for (int i = 0; i < 4; ++i) {
            Eigen::Vector2d dNi_natural(dN[i][0], dN[i][1]);
            Eigen::Vector2d dNi_physical = jac.Jinv * dNi_natural;
            dNdx[i] = dNi_physical(0);
            dNdy[i] = dNi_physical(1);
        }

        // Shear B matrix (2x12) — [w,θx,θy] থেকে shear strain [γxz,γyz] এ
        Eigen::Matrix<double, 2, 12> Bs = Eigen::Matrix<double, 2, 12>::Zero();
        for (int i = 0; i < 4; ++i) {
            const int col = 3 * i;
            // γxz = ∂w/∂x - θy
            Bs(0, col + 0) = dNdx[i];
            Bs(0, col + 2) = -N[i];
            // γyz = ∂w/∂y + θx
            Bs(1, col + 0) = dNdy[i];
            Bs(1, col + 1) = N[i];
        }

        // 1-point reduced integration এ weight = পুরো element area
        // (2x2 natural domain এর পূর্ণ Jacobian-weighted area, একক
        // Gauss point এ পুরো ∫∫dξdη=4 ধরে — center point এ detJ ব্যবহার
        // করে সমগ্র element area আনুমানিক করা, যা reduced integration
        // এর standard practice)
        const double weight = 4.0; // ∫∫ dξ dη over [-1,1]x[-1,1] = 4, single-point rule এ পুরো weight একটাই point এ
        K += Bs.transpose() * Ds * Bs * jac.detJ * weight;
    }

    return K;
}

} // anonymous namespace

Eigen::Matrix<double, 24, 24> computeShellLocalStiffness(
    const std::array<Eigen::Vector3d, 4>& localCorners,
    double thickness,
    double elasticModulus,
    double poissonsRatio
) {
    if (thickness <= 0.0) {
        throw std::invalid_argument("Shell thickness must be positive (got " + std::to_string(thickness) + ")");
    }

    Eigen::Matrix<double, 24, 24> K = Eigen::Matrix<double, 24, 24>::Zero();

    // ---------------------------------------------------------------
    // ধাপ ১: Membrane stiffness (ux, uy) — local DOF index মধ্যে
    // প্রতিটা node এর [ux,uy,uz,rx,ry,rz] থেকে (ux,uy) নেওয়া, অর্থাৎ
    // node i এর জন্য global index 6i+0 (ux), 6i+1 (uy)।
    // ---------------------------------------------------------------
    Eigen::Matrix<double, 8, 8> Kmembrane =
        computeMembraneStiffness(localCorners, thickness, elasticModulus, poissonsRatio);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // node i এর (ux,uy) → node j এর (ux,uy) coupling
            K(6 * i + 0, 6 * j + 0) += Kmembrane(2 * i + 0, 2 * j + 0);
            K(6 * i + 0, 6 * j + 1) += Kmembrane(2 * i + 0, 2 * j + 1);
            K(6 * i + 1, 6 * j + 0) += Kmembrane(2 * i + 1, 2 * j + 0);
            K(6 * i + 1, 6 * j + 1) += Kmembrane(2 * i + 1, 2 * j + 1);
        }
    }

    // ---------------------------------------------------------------
    // ধাপ ২: Drilling DOF (rz) penalty stabilization — diagonal only
    // ---------------------------------------------------------------
    double avgMembraneDiag = 0.0;
    for (int i = 0; i < 8; ++i) avgMembraneDiag += Kmembrane(i, i);
    avgMembraneDiag /= 8.0;

    const double drillingK = drillingPenaltyStiffness(avgMembraneDiag);
    for (int i = 0; i < 4; ++i) {
        K(6 * i + 5, 6 * i + 5) += drillingK; // rz = local index 5 প্রতি node
    }

    // ---------------------------------------------------------------
    // ধাপ ৩: Plate bending stiffness (uz, rx, ry) — Mindlin-Reissner
    // formulation, selective reduced integration (shear locking এড়াতে)
    // ---------------------------------------------------------------
    Eigen::Matrix<double, 12, 12> Kbending =
        computePlateBendingStiffness(localCorners, thickness, elasticModulus, poissonsRatio);

    // Bending DOF per node: uz (index 2), rx (index 3), ry (index 4)
    // — এই ৩টা DOF Kbending এর local ordering এ [w,θx,θy] প্রতি node,
    // তাই node i এর জন্য Kbending index 3i+0(w), 3i+1(θx), 3i+2(θy)
    // → global local shell index 6i+2(uz), 6i+3(rx), 6i+4(ry)
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    const int shellRowOffset = (a == 0) ? 2 : (a == 1 ? 3 : 4);
                    const int shellColOffset = (b == 0) ? 2 : (b == 1 ? 3 : 4);
                    K(6 * i + shellRowOffset, 6 * j + shellColOffset) += Kbending(3 * i + a, 3 * j + b);
                }
            }
        }
    }

    return K;
}

ShellLocalGeometry computeShellLocalGeometry(
    const std::array<Node3D, 4>& globalCorners
) {
    Eigen::Vector3d p0(globalCorners[0].x, globalCorners[0].y, globalCorners[0].z);
    Eigen::Vector3d p1(globalCorners[1].x, globalCorners[1].y, globalCorners[1].z);
    Eigen::Vector3d p2(globalCorners[2].x, globalCorners[2].y, globalCorners[2].z);
    Eigen::Vector3d p3(globalCorners[3].x, globalCorners[3].y, globalCorners[3].z);

    // Local x = corner0 → corner1 দিকের unit vector
    Eigen::Vector3d localX = (p1 - p0).normalized();

    // Local z (normal) = diagonal cross product দিয়ে (p2-p0) × (p3-p1)
    // — এটা 4 corner গড়ে averaged normal দেয় (non-planar/সামান্য
    // warped quad এর জন্যও reasonable), pure (p1-p0)×(p3-p0) এর চেয়ে
    // বেশি robust।
    Eigen::Vector3d diagonal1 = p2 - p0;
    Eigen::Vector3d diagonal2 = p3 - p1;
    Eigen::Vector3d localZ = diagonal1.cross(diagonal2).normalized();

    // Local y = localZ × localX (right-handed orthogonal system সম্পূর্ণ করতে)
    Eigen::Vector3d localY = localZ.cross(localX).normalized();
    // localX পুনরায় ঠিক করা (localY, localZ এর সাথে সঠিকভাবে orthogonal
    // থাকতে, যদি মূল localX সামান্য non-planar approximation এর কারণে
    // perfectly orthogonal না থাকে)
    localX = localY.cross(localZ).normalized();

    if (!std::isfinite(localX.sum()) || !std::isfinite(localY.sum()) || !std::isfinite(localZ.sum())) {
        throw std::invalid_argument(
            "Shell element local coordinate system গণনা ব্যর্থ হয়েছে — degenerate geometry "
            "(zero-area element বা collinear corners) সন্দেহ করা হচ্ছে");
    }

    ShellLocalGeometry result;

    // প্রতিটা corner কে local coordinate এ project করা (origin = corner0)
    for (int i = 0; i < 4; ++i) {
        Eigen::Vector3d p(globalCorners[i].x, globalCorners[i].y, globalCorners[i].z);
        Eigen::Vector3d relative = p - p0;
        double localXCoord = relative.dot(localX);
        double localYCoord = relative.dot(localY);
        double localZCoord = relative.dot(localZ); // flat element ধরে নিলে এটা ≈0 হওয়ার কথা

        result.localCorners[i] = Eigen::Vector3d(localXCoord, localYCoord, localZCoord);
    }

    // 24x24 transformation matrix — প্রতিটা node এ 6x6 rotation block
    // (৩টা translation + ৩টা rotation DOF, একই rotation matrix উভয়ের
    // জন্য প্রযোজ্য — standard rigid-rotation transformation, frame
    // element এর computeTransformationMatrix() এর মতোই যুক্তি)
    Eigen::Matrix3d R;
    R.row(0) = localX;
    R.row(1) = localY;
    R.row(2) = localZ;

    Eigen::Matrix<double, 24, 24> T = Eigen::Matrix<double, 24, 24>::Zero();
    for (int node = 0; node < 4; ++node) {
        T.block<3, 3>(6 * node, 6 * node) = R;
        T.block<3, 3>(6 * node + 3, 6 * node + 3) = R;
    }
    result.transformationMatrix = T;

    return result;
}

Eigen::Matrix<double, 24, 24> computeShellGlobalStiffness(
    const std::array<Node3D, 4>& globalCorners,
    double thickness,
    double elasticModulus,
    double poissonsRatio
) {
    auto geometry = computeShellLocalGeometry(globalCorners);
    auto kLocal = computeShellLocalStiffness(geometry.localCorners, thickness, elasticModulus, poissonsRatio);
    const auto& T = geometry.transformationMatrix;
    return T.transpose() * kLocal * T;
}

} // namespace civilos
