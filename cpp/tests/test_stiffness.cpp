#include "stiffness.h"
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
    // Test 1: Matrix Symmetry
    // ============================================================
    // যেকোনো বৈধ stiffness matrix অবশ্যই symmetric হতে হবে (Maxwell-
    // Betti reciprocal theorem থেকে আসা একটা মৌলিক ভৌত সত্য — energy
    // conservation-এর সরাসরি ফলাফল)। এটা fail করলে বুঝতে হবে matrix
    // এ ভুল আছে, তা মান যতই "যুক্তিসঙ্গত" দেখাক না কেন।
    std::cout << "=== Test 1: Local Stiffness Matrix Symmetry ===\n";
    {
        SectionProperties section{0.09, 0.000675, 0.0003, 0.0005}; // 300x300mm rectangular এর কাছাকাছি
        MaterialProperties material{25e6, 10.4e6}; // E=25000 MPa concrete, G≈E/2.4 এর কাছাকাছি
        auto k = computeLocalStiffnessMatrix(4.0, section, material);

        double maxAsymmetry = (k - k.transpose()).cwiseAbs().maxCoeff();
        checkTrue("Matrix is symmetric (max asymmetry < 1e-9)", maxAsymmetry < 1e-9);
    }

    // ============================================================
    // Test 2: Cantilever Beam — Known Textbook Solution
    // ============================================================
    // একটা cantilever beam (এক প্রান্তে fixed, অন্য প্রান্তে free),
    // free end এ একটা point load P প্রয়োগ করলে ক্লাসিক্যাল বিম থিওরি
    // অনুযায়ী: tip deflection δ = PL³/(3EI), tip rotation θ = PL²/(2EI)
    // এই সূত্র যেকোনো Strength of Materials টেক্সটবুকে (Hibbeler,
    // Timoshenko) পাওয়া যায় — সলভার এই একই মান দেয় কিনা যাচাই করা
    // হচ্ছে stiffness matrix সরাসরি ব্যবহার করে (static condensation
    // দিয়ে node 1 fixed ধরে node 2 এর free-DOF সমাধান করে)।
    std::cout << "\n=== Test 2: Cantilever Beam vs. Textbook Solution ===\n";
    {
        const double L = 3.0;      // m
        const double E = 200e6;    // kN/m² (steel, 200000 MPa)
        const double I = 8.0e-5;   // m⁴ (একটা মাঝারি steel W-shape এর কাছাকাছি Ixx)
        const double A = 0.005;    // m²
        const double G = 77e6;     // kN/m²
        const double J = 5e-7;     // m⁴

        SectionProperties section{A, I, I, J}; // Ixx=Iyy ধরা হয়েছে এই টেস্টে সরলতার জন্য
        MaterialProperties material{E, G};

        auto k = computeLocalStiffnessMatrix(L, section, material);

        // Node 1 (index 0-5) fixed, Node 2 (index 6-11) free।
        // Static condensation: প্রথম 6 DOF eliminate করে বাকি 6x6
        // sub-matrix বের করা হচ্ছে (এই সাব-ম্যাট্রিক্সই আসলে fixed-end
        // cantilever এর reduced stiffness — DOF 6-11 এর জন্য)।
        Eigen::Matrix<double, 6, 6> kReduced = k.block<6, 6>(6, 6);

        // Free end এ শুধু local-y দিকে (strong axis bending) একটা
        // point load P প্রয়োগ — DOF ইনডেক্স: v2=1 (reduced system এ,
        // যা মূল system এর DOF 7 এর সাথে সঙ্গতিপূর্ণ)
        const double P = 10.0; // kN
        Eigen::VectorXd F = Eigen::VectorXd::Zero(6);
        F(1) = P; // local-y force at node 2

        Eigen::VectorXd displacement = kReduced.colPivHouseholderQr().solve(F);

        const double tipDeflection = displacement(1); // v2
        const double tipRotation = displacement(5);   // rz2

        const double expectedDeflection = (P * L * L * L) / (3.0 * E * I);
        const double expectedRotation = (P * L * L) / (2.0 * E * I);

        checkClose("Tip deflection δ = PL³/(3EI)", tipDeflection, expectedDeflection);
        checkClose("Tip rotation θ = PL²/(2EI)", tipRotation, expectedRotation);
    }

    // ============================================================
    // Test 3: Rigid Body Motion — একটা মৌলিক invariant
    // ============================================================
    // যদি পুরো element কে কোনো বিকৃতি ছাড়াই সরানো হয় (rigid body
    // translation — উভয় node একই দিকে একই পরিমাণ সরে), তাহলে internal
    // force শূন্য হওয়া উচিত (কোনো strain নেই, তাই কোনো stress নেই)।
    // এটা F = K * u থেকে যাচাই করা হচ্ছে rigid translation vector দিয়ে।
    std::cout << "\n=== Test 3: Rigid Body Translation Produces Zero Force ===\n";
    {
        SectionProperties section{0.09, 0.000675, 0.0003, 0.0005};
        MaterialProperties material{25e6, 10.4e6};
        auto k = computeLocalStiffnessMatrix(4.0, section, material);

        // উভয় node কে local-x দিকে 0.01m সরানো (rigid translation)
        Eigen::VectorXd rigidTranslation(12);
        rigidTranslation << 0.01, 0, 0, 0, 0, 0,  0.01, 0, 0, 0, 0, 0;

        Eigen::VectorXd force = k * rigidTranslation;
        double maxForce = force.cwiseAbs().maxCoeff();

        checkTrue("Rigid translation → zero internal force (max < 1e-9)", maxForce < 1e-9);
    }

    // ============================================================
    // Test 4: Transformation Matrix Orthogonality
    // ============================================================
    // একটা valid rotation matrix সবসময় orthogonal (T^T * T = Identity)
    // — এটা না হলে length/angle preserve হয় না, যা geometrically অসম্ভব।
    std::cout << "\n=== Test 4: Transformation Matrix Orthogonality ===\n";
    {
        Node3D n1{"n1", 0, 0, 0};
        Node3D n2{"n2", 3, 4, 0}; // একটা diagonal element (3-4-5 triangle, length=5)
        auto T = computeTransformationMatrix(n1, n2);

        Eigen::Matrix<double, 12, 12> shouldBeIdentity = T.transpose() * T;
        Eigen::Matrix<double, 12, 12> identity = Eigen::Matrix<double, 12, 12>::Identity();
        double maxDeviation = (shouldBeIdentity - identity).cwiseAbs().maxCoeff();

        checkTrue("T^T * T = Identity (max deviation < 1e-9)", maxDeviation < 1e-9);
    }

    // ============================================================
    // Test 5: Vertical Column — Special Case Handling
    // ============================================================
    // একটা perfectly vertical member (column) এ reference vector
    // বিশেষভাবে হ্যান্ডল করা হয় (কমেন্টে ব্যাখ্যা করা আছে stiffness.cpp
    // এ) — এই টেস্ট নিশ্চিত করে যে সেই বিশেষ কেসেও transformation
    // matrix বৈধ (orthogonal, non-degenerate) থাকে, crash বা NaN না দেয়।
    std::cout << "\n=== Test 5: Vertical Column Edge Case ===\n";
    {
        Node3D n1{"n1", 0, 0, 0};
        Node3D n2{"n2", 0, 5, 0}; // সম্পূর্ণ vertical, Y-দিকে
        auto T = computeTransformationMatrix(n1, n2);

        bool hasNaN = T.hasNaN();
        checkTrue("Vertical column transformation has no NaN", !hasNaN);

        Eigen::Matrix<double, 12, 12> shouldBeIdentity = T.transpose() * T;
        Eigen::Matrix<double, 12, 12> identity = Eigen::Matrix<double, 12, 12>::Identity();
        double maxDeviation = (shouldBeIdentity - identity).cwiseAbs().maxCoeff();
        checkTrue("Vertical column T is still orthogonal", maxDeviation < 1e-9);
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
