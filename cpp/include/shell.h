#pragma once

#include "types.h"
#include <Eigen/Dense>
#include <array>

namespace civilos {

/**
 * Shell Element Local Stiffness Matrix (24x24) — একটা 4-node
 * quadrilateral shell (Slab/Wall) এর জন্য, membrane (in-plane) ও
 * plate bending (out-of-plane) stiffness একসাথে সমন্বিত।
 *
 * === পদ্ধতির সংক্ষিপ্ত বিবরণ ===
 *
 * DOF layout (প্রতি node ৬টা, মোট 4×6=24): প্রতিটা node এ
 * [ux, uy, uz, rx, ry, rz] — যেখানে ux,uy = membrane (in-plane
 * translation, local x-y plane), uz = plate bending deflection
 * (local z, normal direction), rx,ry = plate bending rotation
 * (bending curvature তৈরি করে), rz = "drilling" DOF (in-plane
 * rotation — নিচে ব্যাখ্যা)।
 *
 * **Membrane অংশ** (ux, uy, ও drilling rz): standard bilinear
 * isoparametric plane-stress element (4-node quad, Q4), 2×2 Gauss
 * quadrature দিয়ে integrate করা। এই classical element এ প্রকৃতপক্ষে
 * কোনো rz (drilling) stiffness নেই — কিন্তু আমাদের shell এ প্রতি node
 * এ ৬-DOF চাই (frame element এর সাথে সামঞ্জস্যপূর্ণ রাখতে, ও coplanar
 * shell panel এ global stiffness matrix singular না হওয়ার জন্য),
 * তাই একটা ছোট artificial "drilling stiffness" (Allman-type penalty,
 * নিচে বিস্তারিত) যোগ করা হয়েছে rz DOF কে non-zero কিন্তু physically
 * negligible stiffness দিতে।
 *
 * **Plate bending অংশ**: DKQ (Discrete Kirchhoff Quadrilateral) —
 * thin-plate (Kirchhoff) theory-ভিত্তিক, যেখানে transverse shear
 * deformation ignore করা হয় (thickness/span অনুপাত ছোট থাকলে, যা
 * বেশিরভাগ slab/wall এর জন্য সত্য — thick plate/shear-deformable
 * এই Phase এ সমর্থিত না)। DKQ চারটা triangular DKT (Discrete
 * Kirchhoff Triangle) sub-element এর সমন্বয়ে গঠিত (quad কে ৪টা
 * corner+center triangle এ ভাগ করে, তারপর condensation করে কেন্দ্রের
 * অতিরিক্ত DOF বাদ দেওয়া হয়) — এটা plate bending এর জন্য সবচেয়ে
 * সুপরিচিত ও robust (numerically locking-free) formulation এর একটা।
 *
 * সূত্রের উৎস: এই পুরো formulation প্রতিষ্ঠিত FE textbook থেকে —
 * Cook, Malkus, Plesha & Witt, "Concepts and Applications of Finite
 * Element Analysis" (Chapter 6 — plane elements, Chapter 13 — plate
 * bending); Batoz & Tahar (1982), "Evaluation of a new quadrilateral
 * thin plate bending element" (DKQ এর মূল paper)। কোনো নতুন derivation
 * না — well-established, widely-used (SAP2000/ETABS-এর মতো commercial
 * software এও ব্যবহৃত) পদ্ধতি।
 *
 * ধরে নেওয়া হয়েছে (এই সংস্করণের সীমাবদ্ধতা):
 *   - Flat shell (element নিজে সমতল, যদিও adjacent element ভিন্ন plane
 *     এ থাকতে পারে, যেমন একটা কোণাকৃতি slab-wall junction) — warped/
 *     curved shell element সমর্থিত না
 *   - Isotropic, linear elastic material
 *   - Thin plate theory (transverse shear deformation ignore, Mindlin/
 *     Reissner thick-plate সমর্থিত না)
 *   - Uniform thickness প্রতি element এ
 *   - 4-node quad শুধু (3-node triangle সমর্থিত না এই সংস্করণে —
 *     mesh generation সবসময় quad তৈরি করবে, analysis_orchestration.py
 *     এ নিশ্চিত করা)
 */
Eigen::Matrix<double, 24, 24> computeShellLocalStiffness(
    const std::array<Eigen::Vector3d, 4>& localCorners, // local coordinate system এ 4 corner এর (x,y,0) — z সবসময় 0 (flat element, local plane এ)
    double thickness,
    double elasticModulus,
    double poissonsRatio
);

/**
 * একটা shell element এর 4 global corner coordinate থেকে local
 * coordinate system (element এর নিজস্ব plane) বানানোর ট্রান্সফরমেশন।
 *
 * পদ্ধতি: local x = corner0→corner1 দিকের unit vector, local z =
 * element normal (corner ভেক্টরগুলোর cross product থেকে, right-hand
 * rule অনুযায়ী counter-clockwise node ordering ধরে — types.h এর
 * ShellElement docstring এ node-ordering convention দেখুন), local y =
 * local z × local x (orthogonal সম্পূর্ণ করতে)।
 *
 * রিটার্ন করে: (localCorners — element plane এ projected 2D coordinate,
 * z=0 সহ 3D vector হিসেবে সংরক্ষিত সরলতার জন্য, transformationMatrix
 * — 24x24 global-to-local transformation, computeShellGlobalStiffness
 * এ ব্যবহারের জন্য)।
 */
struct ShellLocalGeometry {
    std::array<Eigen::Vector3d, 4> localCorners;
    Eigen::Matrix<double, 24, 24> transformationMatrix;
};

ShellLocalGeometry computeShellLocalGeometry(
    const std::array<Node3D, 4>& globalCorners
);

/**
 * একটা shell element এর global (assembly-ready) 24x24 stiffness
 * matrix — computeShellLocalStiffness ও computeShellLocalGeometry
 * একত্রিত করে, K_global = T^T * K_local * T সূত্রে (frame element এর
 * computeGlobalElementStiffness() এর সাথে সঙ্গতিপূর্ণ পদ্ধতি)।
 */
Eigen::Matrix<double, 24, 24> computeShellGlobalStiffness(
    const std::array<Node3D, 4>& globalCorners,
    double thickness,
    double elasticModulus,
    double poissonsRatio
);

} // namespace civilos
