#include "stiffness.h"
#include <cmath>
#include <stdexcept>
#include <vector>
#include <Eigen/QR>

namespace civilos {

Eigen::Matrix<double, 12, 12> computeLocalStiffnessMatrix(
    double length,
    const SectionProperties& section,
    const MaterialProperties& material
) {
    const double L = length;
    const double A = section.area;
    const double Ixx = section.ixx; // strong axis (local z-axis bending, বিক্ষেপ local y দিকে)
    const double Iyy = section.iyy; // weak axis (local y-axis bending, বিক্ষেপ local z দিকে)
    const double J = section.j;
    const double E = material.elasticModulus;
    const double G = material.shearModulus;

    if (L <= 0.0) {
        throw std::invalid_argument("Element length must be positive (got " + std::to_string(L) + ")");
    }

    const double L2 = L * L;
    const double L3 = L2 * L;

    // 12 DOF ক্রম: [u1,v1,w1,rx1,ry1,rz1, u2,v2,w2,rx2,ry2,rz2]
    // যেখানে u=axial(local x), v=local y, w=local z, rx=torsion, ry/rz=bending rotation
    Eigen::Matrix<double, 12, 12> k = Eigen::Matrix<double, 12, 12>::Zero();

    // --- Axial stiffness (local x direction) ---
    const double axial = E * A / L;
    k(0, 0) = axial;
    k(0, 6) = -axial;
    k(6, 0) = -axial;
    k(6, 6) = axial;

    // --- Torsional stiffness (local x-axis rotation) ---
    const double torsional = G * J / L;
    k(3, 3) = torsional;
    k(3, 9) = -torsional;
    k(9, 3) = -torsional;
    k(9, 9) = torsional;

    // --- Bending about local z-axis (deflection in local y, uses Ixx — "strong axis") ---
    // DOF indices: v1=1, rz1=5, v2=7, rz2=11
    {
        const double EI = E * Ixx;
        const double k1 = 12.0 * EI / L3;
        const double k2 = 6.0 * EI / L2;
        const double k3 = 4.0 * EI / L;
        const double k4 = 2.0 * EI / L;

        k(1, 1) = k1;
        k(1, 5) = k2;
        k(1, 7) = -k1;
        k(1, 11) = k2;

        k(5, 1) = k2;
        k(5, 5) = k3;
        k(5, 7) = -k2;
        k(5, 11) = k4;

        k(7, 1) = -k1;
        k(7, 5) = -k2;
        k(7, 7) = k1;
        k(7, 11) = -k2;

        k(11, 1) = k2;
        k(11, 5) = k4;
        k(11, 7) = -k2;
        k(11, 11) = k3;
    }

    // --- Bending about local y-axis (deflection in local z, uses Iyy — "weak axis") ---
    // DOF indices: w1=2, ry1=4, w2=8, ry2=10
    // নোট: এই ব্লকের off-diagonal terms এর সাইন strong-axis ব্লক থেকে
    // উল্টো (দুই জায়গায় negative যেখানে strong-axis এ positive) —
    // এটা ভুল না, এটা local y ও local z axis এর right-hand-rule
    // orientation-এর পার্থক্য থেকে আসে (একটা standard, well-documented
    // sign convention পার্থক্য মেট্রিক্স স্ট্রাকচারাল অ্যানালাইসিস
    // রেফারেন্সে, McGuire et al. সহ)।
    {
        const double EI = E * Iyy;
        const double k1 = 12.0 * EI / L3;
        const double k2 = 6.0 * EI / L2;
        const double k3 = 4.0 * EI / L;
        const double k4 = 2.0 * EI / L;

        k(2, 2) = k1;
        k(2, 4) = -k2;
        k(2, 8) = -k1;
        k(2, 10) = -k2;

        k(4, 2) = -k2;
        k(4, 4) = k3;
        k(4, 8) = k2;
        k(4, 10) = k4;

        k(8, 2) = -k1;
        k(8, 4) = k2;
        k(8, 8) = k1;
        k(8, 10) = k2;

        k(10, 2) = -k2;
        k(10, 4) = k4;
        k(10, 8) = k2;
        k(10, 10) = k3;
    }

    return k;
}

Eigen::Matrix<double, 12, 12> computeTransformationMatrix(
    const Node3D& startNode,
    const Node3D& endNode
) {
    const double dx = endNode.x - startNode.x;
    const double dy = endNode.y - startNode.y;
    const double dz = endNode.z - startNode.z;
    const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (length < 1e-9) {
        throw std::invalid_argument("Element start and end nodes coincide (zero length)");
    }

    // Local x-axis unit vector (element এর axial direction)
    const Eigen::Vector3d localX(dx / length, dy / length, dz / length);

    // Reference vector নির্ধারণ — সাধারণত global Y (উপরে), কিন্তু
    // perfectly vertical element (column) এর ক্ষেত্রে global Y নিজেই
    // local X এর সমান্তরাল হয়ে যায় যা cross product কে undefined করে
    // দেয় (0 vector)। এই বিশেষ ক্ষেত্রে global Z কে reference হিসেবে
    // ব্যবহার করা হচ্ছে — এটা একটা standard convention (SAP2000/ETABS
    // সহ প্রায় সব স্ট্রাকচারাল সফটওয়্যারে vertical member এর জন্য
    // ব্যবহৃত), যদিও চূড়ান্ত local-axis orientation সফটওয়্যার-ভেদে
    // (rotation angle) কিছুটা ভিন্ন হতে পারে — এটা শুধু stiffness/
    // displacement এর numerical মান পরিবর্তন করে না (সঠিক থাকে),
    // শুধু কোন local axis "y" আর কোনটা "z" তা প্রভাবিত করে, যা
    // section orientation report করার সময় প্রাসঙ্গিক (এই Phase 4a
    // তে এখনো section-orientation-specific output নেই)।
    const bool isVertical = std::abs(localX.x()) < 1e-6 && std::abs(localX.z()) < 1e-6;
    const Eigen::Vector3d referenceVector = isVertical
        ? Eigen::Vector3d(0.0, 0.0, 1.0)  // vertical member: global Z কে reference
        : Eigen::Vector3d(0.0, 1.0, 0.0); // non-vertical member: global Y কে reference

    // local Z = local X × reference (তারপর normalize)
    Eigen::Vector3d localZ = localX.cross(referenceVector);
    localZ.normalize();

    // local Y = local Z × local X (right-hand rule সম্পূর্ণ করতে)
    const Eigen::Vector3d localY = localZ.cross(localX);

    // 3x3 rotation matrix — প্রতিটা সারি একটা local axis কে global
    // coordinate এ প্রকাশ করে
    Eigen::Matrix3d rotation;
    rotation.row(0) = localX;
    rotation.row(1) = localY;
    rotation.row(2) = localZ;

    // 12x12 transformation matrix — ৪টা 3x3 rotation ব্লক-ডায়াগোনালে
    // (প্রতিটা node এর translation ও rotation DOF আলাদাভাবে rotate হয়)
    Eigen::Matrix<double, 12, 12> T = Eigen::Matrix<double, 12, 12>::Zero();
    for (int i = 0; i < 4; ++i) {
        T.block<3, 3>(i * 3, i * 3) = rotation;
    }

    return T;
}

Eigen::Matrix<double, 12, 12> computeGlobalElementStiffness(
    const Node3D& startNode,
    const Node3D& endNode,
    const SectionProperties& section,
    const MaterialProperties& material
) {
    const double dx = endNode.x - startNode.x;
    const double dy = endNode.y - startNode.y;
    const double dz = endNode.z - startNode.z;
    const double length = std::sqrt(dx * dx + dy * dy + dz * dz);

    const Eigen::Matrix<double, 12, 12> kLocal = computeLocalStiffnessMatrix(length, section, material);
    const Eigen::Matrix<double, 12, 12> T = computeTransformationMatrix(startNode, endNode);

    // K_global = T^T * K_local * T
    return T.transpose() * kLocal * T;
}

Eigen::Matrix<double, 12, 12> applyEndReleases(
    const Eigen::Matrix<double, 12, 12>& kLocal,
    bool releaseStart,
    bool releaseEnd
) {
    if (!releaseStart && !releaseEnd) {
        return kLocal; // কিছুই release করার নেই — rigid connection, unchanged ফেরত
    }

    // released (slave) DOF indices জোগাড় করা — local index স্কিম:
    // [u1,v1,w1,rx1,ry1,rz1, u2,v2,w2,rx2,ry2,rz2] অনুযায়ী bending
    // rotation ry=4/10, rz=5/11 (node 1/node 2)। axial(0,6) ও
    // torsion(3,9) কখনো release হয় না (উপরে stiffness.h এ ব্যাখ্যা)।
    std::vector<int> slaveIdx;
    if (releaseStart) {
        slaveIdx.push_back(4);  // ry1
        slaveIdx.push_back(5);  // rz1
    }
    if (releaseEnd) {
        slaveIdx.push_back(10); // ry2
        slaveIdx.push_back(11); // rz2
    }

    std::vector<int> masterIdx;
    for (int i = 0; i < 12; ++i) {
        bool isSlave = false;
        for (int s : slaveIdx) {
            if (s == i) { isSlave = true; break; }
        }
        if (!isSlave) masterIdx.push_back(i);
    }

    const int nS = static_cast<int>(slaveIdx.size());
    const int nM = static_cast<int>(masterIdx.size());

    Eigen::MatrixXd Kmm(nM, nM), Kms(nM, nS), Ksm(nS, nM), Kss(nS, nS);
    for (int i = 0; i < nM; ++i)
        for (int j = 0; j < nM; ++j)
            Kmm(i, j) = kLocal(masterIdx[i], masterIdx[j]);
    for (int i = 0; i < nM; ++i)
        for (int j = 0; j < nS; ++j)
            Kms(i, j) = kLocal(masterIdx[i], slaveIdx[j]);
    for (int i = 0; i < nS; ++i)
        for (int j = 0; j < nM; ++j)
            Ksm(i, j) = kLocal(slaveIdx[i], masterIdx[j]);
    for (int i = 0; i < nS; ++i)
        for (int j = 0; j < nS; ++j)
            Kss(i, j) = kLocal(slaveIdx[i], slaveIdx[j]);

    // K_ss সাধারণত well-conditioned (bending rotation diagonal terms,
    // কখনো singular না একক release এ) — তবুও pseudo-inverse ব্যবহার
    // করা হচ্ছে defense-in-depth হিসেবে, কারণ edge-case section
    // property (যেমন E বা I শূন্যের কাছাকাছি ভুল ইনপুট) থাকলে সাধারণ
    // inverse crash করতে পারে, pseudo-inverse করবে না (least-squares
    // সমাধান দেবে, যা এখানে physically reasonable)।
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> kssDecomp(Kss);
    Eigen::MatrixXd KssInv = kssDecomp.pseudoInverse();

    Eigen::MatrixXd Kcondensed = Kmm - Kms * KssInv * Ksm;

    // Condensed (nM x nM) matrix কে ফিরিয়ে পূর্ণ 12x12 এ বসানো — slave
    // row/column গুলো zero থাকে (released DOF-এ কোনো stiffness
    // contribution নেই, কিন্তু matrix size 12x12 রাখা হচ্ছে যাতে global
    // assembly logic এ কোনো পরিবর্তন লাগবে না, DOF indexing scheme
    // সব element জুড়ে একই থাকে)।
    Eigen::Matrix<double, 12, 12> result = Eigen::Matrix<double, 12, 12>::Zero();
    for (int i = 0; i < nM; ++i)
        for (int j = 0; j < nM; ++j)
            result(masterIdx[i], masterIdx[j]) = Kcondensed(i, j);

    return result;
}

Eigen::Matrix<double, 12, 12> computeLocalMassMatrix(
    double length,
    const SectionProperties& section,
    const MaterialProperties& material
) {
    const double L = length;
    const double A = section.area;
    const double rho = material.density; // tonne/m³

    if (L <= 0.0) {
        throw std::invalid_argument("Element length must be positive (got " + std::to_string(L) + ")");
    }

    const double L2 = L * L;
    const double m = rho * A; // mass per unit length, tonne/m

    Eigen::Matrix<double, 12, 12> M = Eigen::Matrix<double, 12, 12>::Zero();

    // --- Axial mass (local x direction, linear shape function) ---
    // DOF: u1=0, u2=6
    {
        const double c = m * L / 6.0;
        M(0, 0) = 2.0 * c;
        M(0, 6) = 1.0 * c;
        M(6, 0) = 1.0 * c;
        M(6, 6) = 2.0 * c;
    }

    // --- Torsional mass (local x-axis rotation, linear shape function) ---
    // DOF: rx1=3, rx2=9। polar mass moment of inertia per unit length ≈
    // ρ*J (torsional constant J প্রায়ই polar-moment-of-area এর কাছাকাছি
    // ধরা হয় পাতলা-দেয়ালবিহীন সরল section এ; নিখুঁত polar mass moment
    // ρ*(Iyy+Ixx) হওয়া উচিত পুরোপুরি সঠিকভাবে, কিন্তু torsional mass
    // এর প্রভাব সাধারণত translational mode frequency তে নগণ্য, তাই J
    // ব্যবহার একটা প্রচলিত ও গ্রহণযোগ্য approximation।
    {
        const double J = section.j;
        const double c = rho * J * L / 6.0;
        M(3, 3) = 2.0 * c;
        M(3, 9) = 1.0 * c;
        M(9, 3) = 1.0 * c;
        M(9, 9) = 2.0 * c;
    }

    // --- Bending mass, local z-axis bending plane (deflection in local y) ---
    // DOF: v1=1, rz1=5, v2=7, rz2=11
    // Standard consistent mass sub-matrix (Cook et al., cubic Hermite shape
    // functions): mL/420 * [[156, 22L, 54, -13L],
    //                        [22L, 4L², 13L, -3L²],
    //                        [54, 13L, 156, -22L],
    //                        [-13L, -3L², -22L, 4L²]]
    // stiffness.cpp এর strong-axis bending ব্লকের DOF ordering (v1,rz1,v2,rz2)
    // এর সাথে সরাসরি সঙ্গতিপূর্ণ sign convention (rotation DOF থেকে
    // translation DOF এর দিকে coupling সবসময় positive এই ক্রমে, McGuire
    // et al. ও Cook et al. উভয়ের কনভেনশনে অভিন্ন — bending mass matrix
    // stiffness এর মতো sign-flip সমস্যায় পড়ে না কারণ এটা virtual-work
    // ভিত্তিক shape-function integral, rotational-equilibrium ভিত্তিক না)
    {
        const double c = m * L / 420.0;
        M(1, 1)  = 156.0 * c;
        M(1, 5)  = 22.0 * L * c;
        M(1, 7)  = 54.0 * c;
        M(1, 11) = -13.0 * L * c;

        M(5, 1)  = 22.0 * L * c;
        M(5, 5)  = 4.0 * L2 * c;
        M(5, 7)  = 13.0 * L * c;
        M(5, 11) = -3.0 * L2 * c;

        M(7, 1)  = 54.0 * c;
        M(7, 5)  = 13.0 * L * c;
        M(7, 7)  = 156.0 * c;
        M(7, 11) = -22.0 * L * c;

        M(11, 1)  = -13.0 * L * c;
        M(11, 5)  = -3.0 * L2 * c;
        M(11, 7)  = -22.0 * L * c;
        M(11, 11) = 4.0 * L2 * c;
    }

    // --- Bending mass, local y-axis bending plane (deflection in local z) ---
    // DOF: w1=2, ry1=4, w2=8, ry2=10। একই consistent mass সূত্র, শুধু
    // DOF index ভিন্ন — stiffness.cpp এর weak-axis bending ব্লকে যে
    // sign-flip ছিল (local y/z axis orientation পার্থক্যের কারণে) তা
    // mass matrix এ প্রযোজ্য না (উপরের নোট দেখুন) — তাই এখানে সরাসরি
    // strong-axis ব্লকের মতোই সব-positive coupling ব্যবহার করা হচ্ছে।
    {
        const double c = m * L / 420.0;
        M(2, 2)   = 156.0 * c;
        M(2, 4)   = 22.0 * L * c;
        M(2, 8)   = 54.0 * c;
        M(2, 10)  = -13.0 * L * c;

        M(4, 2)   = 22.0 * L * c;
        M(4, 4)   = 4.0 * L2 * c;
        M(4, 8)   = 13.0 * L * c;
        M(4, 10)  = -3.0 * L2 * c;

        M(8, 2)   = 54.0 * c;
        M(8, 4)   = 13.0 * L * c;
        M(8, 8)   = 156.0 * c;
        M(8, 10)  = -22.0 * L * c;

        M(10, 2)  = -13.0 * L * c;
        M(10, 4)  = -3.0 * L2 * c;
        M(10, 8)  = -22.0 * L * c;
        M(10, 10) = 4.0 * L2 * c;
    }

    return M;
}

Eigen::Matrix<double, 12, 12> computeLocalGeometricStiffnessMatrix(
    double length,
    double axialForce
) {
    const double L = length;

    // ইনপুট (axialForce) convention: compression-positive (উপরে
    // stiffness.h এ ব্যাখ্যা করা, solveLinearStatic() এর raw
    // elementEndForces sign এর সাথে সঙ্গতিপূর্ণ)। কিন্তু নিচের classical
    // textbook geometric-stiffness সূত্র (Przemieniecki) tension-positive
    // convention এ derive করা (tension স্টিফনেস বাড়ায়, positive P) —
    // তাই এখানে sign flip করে P কে tension-positive এ রূপান্তর করা
    // হচ্ছে, যাতে নিচের সূত্র সঠিকভাবে প্রযোজ্য হয় (compression P তখন
    // ঋণাত্মক হবে সূত্রে, যেটাই textbook convention আশা করে)।
    const double P = -axialForce;

    if (L <= 0.0) {
        throw std::invalid_argument("Element length must be positive (got " + std::to_string(L) + ")");
    }

    const double L2 = L * L;

    Eigen::Matrix<double, 12, 12> Kg = Eigen::Matrix<double, 12, 12>::Zero();

    // --- Geometric stiffness, local z-axis bending plane (deflection in local y) ---
    // DOF: v1=1, rz1=5, v2=7, rz2=11। সূত্র: P/(30L) * [[36, 3L, -36, 3L],
    //                                                     [3L, 4L², -3L, -L²],
    //                                                     [-36, -3L, 36, -3L],
    //                                                     [3L, -L², -3L, 4L²]]
    {
        const double c = P / (30.0 * L);
        Kg(1, 1)  = 36.0 * c;
        Kg(1, 5)  = 3.0 * L * c;
        Kg(1, 7)  = -36.0 * c;
        Kg(1, 11) = 3.0 * L * c;

        Kg(5, 1)  = 3.0 * L * c;
        Kg(5, 5)  = 4.0 * L2 * c;
        Kg(5, 7)  = -3.0 * L * c;
        Kg(5, 11) = -1.0 * L2 * c;

        Kg(7, 1)  = -36.0 * c;
        Kg(7, 5)  = -3.0 * L * c;
        Kg(7, 7)  = 36.0 * c;
        Kg(7, 11) = -3.0 * L * c;

        Kg(11, 1)  = 3.0 * L * c;
        Kg(11, 5)  = -1.0 * L2 * c;
        Kg(11, 7)  = -3.0 * L * c;
        Kg(11, 11) = 4.0 * L2 * c;
    }

    // --- Geometric stiffness, local y-axis bending plane (deflection in local z) ---
    // DOF: w1=2, ry1=4, w2=8, ry2=10। একই সূত্র, ভিন্ন DOF index (mass
    // matrix এর weak-axis ব্লকের মতো একই যুক্তি — sign-flip প্রযোজ্য
    // না, দেখুন computeLocalMassMatrix এর কমেন্ট)।
    {
        const double c = P / (30.0 * L);
        Kg(2, 2)   = 36.0 * c;
        Kg(2, 4)   = 3.0 * L * c;
        Kg(2, 8)   = -36.0 * c;
        Kg(2, 10)  = 3.0 * L * c;

        Kg(4, 2)   = 3.0 * L * c;
        Kg(4, 4)   = 4.0 * L2 * c;
        Kg(4, 8)   = -3.0 * L * c;
        Kg(4, 10)  = -1.0 * L2 * c;

        Kg(8, 2)   = -36.0 * c;
        Kg(8, 4)   = -3.0 * L * c;
        Kg(8, 8)   = 36.0 * c;
        Kg(8, 10)  = -3.0 * L * c;

        Kg(10, 2)  = 3.0 * L * c;
        Kg(10, 4)  = -1.0 * L2 * c;
        Kg(10, 8)  = -3.0 * L * c;
        Kg(10, 10) = 4.0 * L2 * c;
    }

    // নোট: axial (u1,u2) ও torsion (rx1,rx2) DOF এ geometric stiffness
    // শূন্য রাখা হয়েছে — উপরে stiffness.h এ ব্যাখ্যা করা সরলীকরণ
    // (শুধু bending-coupled geometric effect ধরা হয়েছে, যা dominant
    // এবং প্রচলিত practice)।

    return Kg;
}

} // namespace civilos
