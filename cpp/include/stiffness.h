#pragma once

#include "types.h"
#include <Eigen/Dense>

namespace civilos {

/**
 * Local Stiffness Matrix — একটা 3D frame element এর 12x12 stiffness
 * matrix, element-এর নিজস্ব local coordinate system এ (local x = axial
 * direction, node1→node2)।
 *
 * সূত্রের উৎস: এই সূত্রগুলো standard structural analysis textbook
 * (McGuire, Gallagher & Ziemian, "Matrix Structural Analysis", অথবা
 * Przemieniecki, "Theory of Matrix Structural Analysis") এ প্রতিষ্ঠিত
 * classical 3D frame element stiffness matrix — কোনো নতুন derivation
 * না, well-known engineering mechanics।
 *
 * ধরে নেওয়া হয়েছে (Phase 4a এর সীমাবদ্ধতা):
 *   - Prismatic element (দৈর্ঘ্য জুড়ে uniform cross-section)
 *   - Linear elastic material (Hooke's law, ছোট deformation)
 *   - Euler-Bernoulli beam theory (shear deformation ignore করা হয়েছে
 *     — Timoshenko beam এর shear correction factor প্রয়োগ করা হয়নি;
 *     এটা slender member এ (span/depth > ~10, যা বেশিরভাগ RC/steel
 *     beam এর জন্য সত্য) নির্ভুল, কিন্তু deep beam এ error বাড়াতে পারে)
 *   - "moment" connectionType — rigid (fully moment-connected) জোড়া
 *     ধরে নেওয়া হয়েছে
 *
 * ⚠️ "pin" connectionType (Brace-এর ডিফল্ট, Phase 2c) এই Phase 4a তে
 * এখনো প্রয়োগ করা হয়নি — অর্থাৎ এই মুহূর্তে সব element (Brace সহ)
 * rigid/moment-connected হিসেবে সলভ হবে, connectionType field দেখা
 * সত্ত্বেও। এটা একটা সৎভাবে স্বীকার করা সীমাবদ্ধতা: pin release এর
 * জন্য static condensation (matrix partition করে rotational DOF
 * eliminate করা) প্রয়োজন, যেটা এই ফাইলে এখনো implement করা হয়নি।
 * এই মুহূর্তে Brace ব্যবহার করলে ফলাফল conservative না হয়ে ভুল দিকে
 * যেতে পারে (moment capacity থাকার মতো stiff দেখাবে, বাস্তবে যা নেই)
 * — তাই Brace-সহ মডেলের ফলাফল এই Phase এ ব্যবহার করা উচিত না,
 * পরবর্তী ধাপে pin release যোগ না হওয়া পর্যন্ত।
 */
Eigen::Matrix<double, 12, 12> computeLocalStiffnessMatrix(
    double length,
    const SectionProperties& section,
    const MaterialProperties& material
);

/**
 * Transformation Matrix — local coordinate থেকে global coordinate এ
 * রূপান্তরের জন্য 12x12 matrix (৪টা 3x3 rotation matrix ব্লক-ডায়াগোনালে)।
 *
 * এই ফাংশনের জন্য element এর local x-axis (node1→node2 direction) এবং
 * একটা "reference vector" লাগে যা local y-axis এর দিক নির্ধারণ করে
 * (কারণ শুধু axial direction দিয়ে element এর ঘূর্ণন/orientation পুরো
 * নির্ধারিত হয় না — একটা সরলরেখার চারপাশে অসীম সংখ্যক orientation
 * সম্ভব)। Beam/Column এর জন্য এই reference vector সাধারণত global
 * Y-axis (উপরে) থেকে derive করা হয়, perfectly vertical column ছাড়া
 * (সেক্ষেত্রে বিশেষ ব্যবস্থা লাগে, নিচে দেখুন)।
 */
Eigen::Matrix<double, 12, 12> computeTransformationMatrix(
    const Node3D& startNode,
    const Node3D& endNode
);

/**
 * Global Stiffness Matrix (element level) — local matrix কে
 * transformation matrix দিয়ে rotate করে global coordinate এ আনে।
 * K_global = T^T * K_local * T
 */
Eigen::Matrix<double, 12, 12> computeGlobalElementStiffness(
    const Node3D& startNode,
    const Node3D& endNode,
    const SectionProperties& section,
    const MaterialProperties& material
);

/**
 * End Release (Static Condensation) — "pin" connectionType এর জন্য
 * bending moment DOF (ry, rz — উভয় প্রান্তে) release করে। Torsion (rx)
 * ও axial (u) DOF release করা হয় না — বাস্তব Brace/truss member এ
 * axial force ও কিছুটা torsional continuity থাকতে পারে, শুধু bending
 * moment transfer করে না (pin/hinge connection এর সংজ্ঞা অনুযায়ী)।
 *
 * পদ্ধতি: Static Condensation (Guyan reduction এর বিশেষ ক্ষেত্র)।
 * released DOF গুলোকে "slave" ধরে, বাকি DOF ("master") এর সাপেক্ষে
 * elimination করা হয়:
 *
 *   K_condensed = K_mm - K_ms * K_ss^{-1} * K_sm
 *
 * যেখানে subscript m=master (retained), s=slave (released)। এটা
 * পদার্থগতভাবে বোঝায়: released DOF-এ কোনো external moment নেই
 * (pin এ moment=0 by definition), তাই সেই DOF এর equilibrium equation
 * থেকে slave displacement/rotation কে master-এর function হিসেবে
 * সমাধান করে বাকি matrix-এ substitute করা হয় — released DOF সম্পূর্ণ
 * বাদ না দিয়ে (matrix size 12x12-ই থাকে, released row/column
 * effectively zero stiffness contribute করে বাকি structure-এ, কিন্তু
 * global assembly-তে size পরিবর্তন করতে হয় না)।
 *
 * releaseStart/releaseEnd = true হলে সেই প্রান্তের ry ও rz (bending
 * rotation DOF, local indices 4,5 প্রথম node-এ ও 10,11 দ্বিতীয়
 * node-এ) release করা হয়। একটা প্রান্ত release (releaseStart=true,
 * releaseEnd=false) সমর্থিত — যেমন one moment-connected/one-pinned
 * bracing, যদিও Phase 4a এর Brace default হলো উভয় প্রান্ত pin।
 *
 * K_ss singular হলে (উভয় প্রান্তে release করা torsion-free element,
 * বিরল edge case) pseudo-inverse (Eigen::CompleteOrthogonalDecomposition)
 * ব্যবহার করা হয় সাধারণ inverse এর বদলে, যাতে crash না করে।
 */
Eigen::Matrix<double, 12, 12> applyEndReleases(
    const Eigen::Matrix<double, 12, 12>& kLocal,
    bool releaseStart,
    bool releaseEnd
);

/**
 * Local Consistent Mass Matrix — একটা 3D frame element এর 12x12 mass
 * matrix, local coordinate system এ (K matrix এর মতো একই DOF ordering:
 * [u1,v1,w1,rx1,ry1,rz1, u2,v2,w2,rx2,ry2,rz2])।
 *
 * "Consistent" মানে কী: mass কে node-এ lump (শুধু translational DOF এ
 * mass/2 করে ভাগ করে বসানো) না করে, সেই একই shape function ব্যবহার
 * করা হয়েছে যা stiffness matrix derive করতে ব্যবহৃত হয়েছিল (axial এর
 * জন্য linear shape function, bending এর জন্য cubic Hermite shape
 * function)। এটা lumped mass এর চেয়ে বেশি accurate, বিশেষত উচ্চতর
 * mode (2nd, 3rd mode) এর frequency তে — এটাই SAP2000/ETABS এর
 * ডিফল্ট পদ্ধতি।
 *
 * সূত্রের উৎস: Cook, Malkus, Plesha & Witt, "Concepts and Applications
 * of Finite Element Analysis" (4th ed.), Chapter 11 — অথবা সমতুল্যভাবে
 * Przemieniecki "Theory of Matrix Structural Analysis" — classical
 * consistent mass matrix for 3D Euler-Bernoulli frame element। কোনো
 * নতুন derivation না, well-established textbook সূত্র।
 *
 * ধরে নেওয়া হয়েছে (Phase 4 এর এই ধাপের সীমাবদ্ধতা):
 *   - Rotary inertia (rotational DOF এর নিজস্ব mass moment of inertia)
 *     ignore করা হয়েছে — শুধু translational mass distribution বিবেচনা
 *     করা হয়েছে। বেশিরভাগ slender frame member এ (span/depth > ~10)
 *     rotary inertia এর প্রভাব transverse (translational) inertia এর
 *     তুলনায় নগণ্য, তাই এটা একটা সাধারণ ও গ্রহণযোগ্য সরলীকরণ (অনেক
 *     commercial software এও এটা একটা default/optional সেটিং)। তবে
 *     এটা একটা approximation — খুব খাটো, গভীর member এ সামান্য error
 *     আনতে পারে।
 *   - Torsional mass distribution axial এর মতোই linear shape function
 *     দিয়ে ধরা হয়েছে (এটাও standard practice)।
 *   - "pin" connectionType এর mass matrix-এ কোনো condensation প্রয়োগ
 *     করা হয় না — কারণ mass matrix-এ moment DOF এর condensation করলে
 *     eigenvalue problem এর সামঞ্জস্য নষ্ট হতে পারে (K ও M একই
 *     transformation এর সাপেক্ষে condense করতে হবে একসাথে, নাহলে
 *     generalized eigenproblem ভুল ফলাফল দেবে)। এই মুহূর্তে (Modal
 *     Analysis প্রথম সংস্করণ) pin-connected element এর mass matrix
 *     rigid element এর মতোই গণনা করা হয় — এটা conservative না হয়ে
 *     সামান্য ভুল দিকেও যেতে পারে pin-heavy মডেলে, নিচে warning এ
 *     জানানো থাকবে (analysis_orchestration.py)।
 */
Eigen::Matrix<double, 12, 12> computeLocalMassMatrix(
    double length,
    const SectionProperties& section,
    const MaterialProperties& material
);

/**
 * Local Geometric Stiffness Matrix — একটা 3D frame element এর 12x12
 * geometric (stress) stiffness matrix, local coordinate এ (K ও M
 * matrix এর মতো একই DOF ordering)।
 *
 * উদ্দেশ্য: axial force (P) থাকলে member এর transverse bending
 * stiffness পরিবর্তিত হয় — compression (P ঋণাত্মক, এই ফাংশনের
 * convention এ নিচে দেখুন) effective stiffness কমায় (member আরও
 * সহজে buckle/deflect করে), tension বাড়ায়। এই প্রভাব stiffness.cpp
 * এর computeLocalStiffnessMatrix() এ ধরা পড়ে না (সেটা pure material/
 * geometric elastic stiffness, force-independent) — এই আলাদা matrix
 * সেই gap পূরণ করে।
 *
 * সাইন কনভেনশন (numerically যাচাই করা, test_buckling_analysis.cpp এর
 * ডিবাগ সেশনে): এই ফাংশনের axialForce প্যারামিটার agrees করে
 * solveLinearStatic() এর AnalysisResult::elementEndForces[e](0) (local
 * index 0 = start-node axial force) এর convention এর সাথে —
 * **positive elementEndForces[e](0) মানে compression, ঋণাত্মক মানে
 * tension** (এটা solveLinearStatic() এর f_local = k_local·u_local
 * থেকে আসা raw sign, যা conventional engineering "tension-positive"
 * থেকে উল্টো — কারণ এখানে end-1 এর force component member এর নিজের
 * উপর ক্রিয়াশীল internal force না, বরং node-এর উপর element কর্তৃক
 * প্রযুক্ত reaction force হিসেবে বের করা হয়েছে, যার sign compression
 * এ positive আসে এই particular DOF-ordering ও force-extraction
 * পদ্ধতিতে)। এই ফাংশন সেই raw convention *সরাসরি* গ্রহণ করে (কোনো
 * sign-flip নিজে করে না) — caller (solver.cpp এর solveLinearBuckling)
 * দায়িত্বশীল সঠিক sign পাস করার জন্য, নিচে সেখানে বিস্তারিত কমেন্ট
 * আছে।
 *
 * এই কনভেনশন compression এ Kg কে সঠিক দিকে (stiffness কমানোর দিকে)
 * কাজ করায় solver.cpp এর buckling eigenvalue সমীকরণে (K φ = -λ·Kg φ)।
 *
 * সূত্রের উৎস: Przemieniecki, "Theory of Matrix Structural Analysis",
 * Chapter 10 — classical consistent geometric stiffness matrix for a
 * 3D Euler-Bernoulli frame element (McGuire, Gallagher & Ziemian এও
 * সমতুল্য derivation)। কোনো নতুন derivation না।
 *
 * সরলীকরণ (এই ধাপের সীমাবদ্ধতা): শুধু axial force এর প্রভাব ধরা
 * হয়েছে bending geometric stiffness এ (সবচেয়ে গুরুত্বপূর্ণ ও প্রচলিত
 * সরলীকরণ — সাধারণ column/frame buckling এ dominant term)। Shear
 * force ও torsion এর geometric stiffness contribution (যা সাধারণত
 * অনেক ছোট, বিশেষ ক্ষেত্র ছাড়া) এই সংস্করণে অন্তর্ভুক্ত করা হয়নি।
 */
Eigen::Matrix<double, 12, 12> computeLocalGeometricStiffnessMatrix(
    double length,
    double axialForce
);

} // namespace civilos
