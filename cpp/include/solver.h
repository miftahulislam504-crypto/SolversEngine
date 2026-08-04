#pragma once

#include "types.h"
#include <Eigen/Sparse>

namespace civilos {

/**
 * Global Stiffness Matrix Assembly — প্রতিটা frame element এর 12x12
 * এবং প্রতিটা shell element এর 24x24 global stiffness কে পুরো
 * structure এর NxN global matrix এ সঠিক position এ যোগ করে (N = 6 *
 * node সংখ্যা)।
 *
 * Sparse matrix ব্যবহার করা হচ্ছে dense এর বদলে — বাস্তব building
 * model এ hundreds/thousands node থাকতে পারে, তাই NxN dense matrix
 * (N=6000 হলে 36 million entries) মেমরিতে ও গণনায় অপ্রয়োজনীয়ভাবে
 * ভারী হবে, যেখানে বাস্তবে বেশিরভাগ entry শূন্য (একটা node সাধারণত
 * মাত্র কয়েকটা element এর সাথে সংযুক্ত, সব node এর সাথে না)।
 *
 * ⚠️ Shell element (Slab/Wall) এই মুহূর্তে শুধু এই ফাংশনে (তাই শুধু
 * solveLinearStatic এ) সমর্থিত — Modal/Buckling/P-Delta এখনো
 * shellElements ignore করে (mass matrix ও geometric stiffness এখনো
 * shell এর জন্য implement করা হয়নি, shell.h এর ভবিষ্যত সম্প্রসারণ)।
 */
Eigen::SparseMatrix<double> assembleGlobalStiffness(const AnalysisModel& model);

/**
 * Boundary Condition প্রয়োগ — Penalty method ব্যবহার করা হচ্ছে
 * (restrained DOF এর diagonal এ একটা অতি-বড় সংখ্যা যোগ করা, matrix
 * থেকে row/column সম্পূর্ণ মুছে ফেলার বদলে)। কারণ: elimination method
 * (row/column বাদ দেওয়া) matrix resize এবং index remapping দাবি করে
 * যা sparse matrix এ ব্যয়বহুল ও bug-prone; penalty method এ matrix
 * size অপরিবর্তিত থাকে, শুধু constrained DOF এ stiffness এত বেশি
 * বসানো হয় যে সেই DOF এর displacement কার্যত শূন্যের সমান হয়ে যায়
 * (যথেষ্ট বড় penalty factor সহ, numerically stable থাকে)।
 */
void applyBoundaryConditions(
    Eigen::SparseMatrix<double>& globalStiffness,
    const std::vector<BoundaryCondition>& boundaryConditions,
    int totalDOF
);

/**
 * Global Load Vector তৈরি — NodalLoad গুলোকে সঠিক position এ বসিয়ে
 * একটা N-length vector বানায়।
 */
Eigen::VectorXd assembleGlobalLoadVector(const AnalysisModel& model, int totalDOF);

/**
 * সম্পূর্ণ Linear Static Analysis — assembly, boundary condition,
 * সমাধান (K*u=F), এবং element end force এর হিসাব একসাথে করে।
 *
 * Solver পদ্ধতি: Sparse Cholesky decomposition (Eigen::SimplicialLDLT)
 * — একটা symmetric positive-definite matrix (যেকোনো বৈধ, স্থিতিশীল
 * structure এর global stiffness matrix, boundary condition প্রয়োগের
 * পর) এর জন্য এটা দ্রুততম ও সবচেয়ে নির্ভরযোগ্য পদ্ধতির একটা, generic
 * LU decomposition এর চেয়ে।
 *
 * যদি matrix singular হয় (কাঠামো অস্থির — যেমন কোনো support নেই,
 * বা একটা mechanism তৈরি হয়ে গেছে ভুল modeling এর কারণে), সমাধান
 * ব্যর্থ হবে এবং AnalysisResult.success = false হবে, কোনো ভুল সংখ্যা
 * silently না দিয়ে।
 */
AnalysisResult solveLinearStatic(const AnalysisModel& model);

/**
 * Global Mass Matrix Assembly — assembleGlobalStiffness এর সমান্তরাল
 * পদ্ধতি, কিন্তু element এর consistent mass matrix (stiffness.h এর
 * computeLocalMassMatrix) ব্যবহার করে। DOF indexing scheme, sparse
 * triplet-assembly পদ্ধতি — সবকিছু stiffness assembly এর অনুরূপ, যাতে
 * K ও M একই global DOF ordering এ থাকে (generalized eigenproblem
 * K φ = ω² M φ এর জন্য এটা আবশ্যক — দুই matrix এর সারি/কলাম
 * mismatch হলে ফলাফল সম্পূর্ণ ভুল হবে)।
 */
Eigen::SparseMatrix<double> assembleGlobalMass(const AnalysisModel& model);

/**
 * Modal Analysis — Natural Frequency এবং Mode Shape বের করার জন্য
 * generalized eigenvalue problem সমাধান করে: K φ = ω² M φ
 *
 * Boundary condition handling: Linear Static এর penalty method এখানে
 * ব্যবহার করা যাবে না — কারণ penalty method কৃত্রিমভাবে বিশাল
 * stiffness যোগ করে (কিন্তু সেই DOF এর mass তো আসল/ছোট মানই থেকে
 * যায়), ফলে ω² = huge_K/small_M একটা spurious অতি-উচ্চ-frequency
 * mode তৈরি করবে যা আসল mode গুলোর সাথে eigen-solver কে বিভ্রান্ত
 * করবে বা numerically unstable করে দেবে। তাই এখানে elimination
 * method ব্যবহার করা হয়েছে — restrained DOF সম্পূর্ণ matrix থেকে বাদ
 * দিয়ে (row ও column দুটোই মুছে) একটা ছোট reduced system বানানো হয়,
 * তারপর সেই reduced K,M এ eigen solve করা হয়।
 *
 * সমাধান পদ্ধতি: Dense generalized symmetric eigensolver
 * (Eigen::GeneralizedSelfAdjointEigenSolver)। Sparse eigensolver
 * (যেমন Spectra library) এখানে ব্যবহার করা হয়নি কারণ এটা একটা
 * external dependency যোগ করত (এই প্রজেক্টের network/build
 * সীমাবদ্ধতায় নতুন dependency যোগ করা এড়ানো হয়েছে), এবং Phase 4a এর
 * টার্গেট মডেল সাইজ (কয়েকশো element পর্যন্ত, README দেখুন) এর জন্য
 * dense solve যথেষ্ট দ্রুত। বড় মডেলে (হাজার হাজার DOF) dense solve
 * ধীর হয়ে যাবে — সেটা একটা future optimization (sparse iterative
 * eigensolver, শুধু প্রথম N mode দরকার হলে যথেষ্ট, পুরো spectrum না)।
 *
 * numModes: কতগুলো সবচেয়ে-নিচু-frequency mode ফেরত দিতে হবে (0 বা
 * ঋণাত্মক দিলে সব mode — reduced system এর DOF সংখ্যা যতগুলো)।
 * সাধারণত building analysis এ প্রথম ৬-১২টা mode-ই যথেষ্ট (মোট mass
 * এর সিংহভাগ ধরার জন্য, modal participation factor অনুযায়ী), তাই
 * ডিফল্ট আচরণ প্রথম কয়েকটা mode ফেরত দেওয়া, পুরো spectrum না।
 */
ModalAnalysisResult solveModalAnalysis(const AnalysisModel& model, int numModes = 12);

/**
 * Linear (Eigenvalue) Buckling Analysis — critical load factor ও
 * buckling mode shape বের করার জন্য।
 *
 * পদ্ধতি:
 *   ১. দেওয়া load pattern (model.loads) দিয়ে প্রথমে Linear Static
 *      সমাধান চালানো হয় — প্রতিটা element এর axial force বের করতে।
 *      এটা একটা প্রয়োজনীয় পূর্বধাপ (types.h এর BucklingAnalysisResult
 *      docstring এ ব্যাখ্যা করা কারণে — buckling load-নির্ভর, তাই কোন
 *      load pattern তা জানা লাগবেই)।
 *   ২. সেই axial force দিয়ে প্রতিটা element এর geometric stiffness
 *      matrix (stiffness.h এর computeLocalGeometricStiffnessMatrix)
 *      বানিয়ে global geometric stiffness (Kg) assemble করা হয়।
 *   ৩. তাত্ত্বিক সমীকরণ: K φ = -λ Kg φ। কিন্তু এটা সরাসরি
 *      Eigen::GeneralizedSelfAdjointEigenSolver(K, -Kg) আকারে solve
 *      করা যায় না — কারণ Kg *সবসময়* singular (axial ও torsion DOF এ
 *      কোনো geometric stiffness নেই, computeLocalGeometricStiffnessMatrix
 *      এর ইচ্ছাকৃত সরলীকরণ), আর এই solver B matrix (এখানে -Kg) কে
 *      strictly positive-definite দাবি করে। তাই সমীকরণ উল্টে
 *      reformulate করা হয়েছে: (-Kg) φ = μ K φ (μ = 1/λ), যেখানে K
 *      (properly-restrained structure এ) সবসময় positive-definite।
 *      সমাধান শেষে λ = 1/μ তে ফিরিয়ে আনা হয়, এবং μ≈0 (অর্থাৎ λ→∞,
 *      non-physical/irrelevant mode) বাদ দেওয়া হয়। Modal Analysis এর
 *      মতোই elimination-method boundary condition ব্যবহার করা হয়েছে
 *      (penalty method না — একই কারণে, spurious eigenvalue এড়াতে)।
 *      বিস্তারিত reformulation logic ও derivation solver.cpp এ।
 *
 * নোট: এই ফাংশন Linear Static কে internally আরেকবার কল করে (ধাপ ১),
 * তাই caller কে আলাদাভাবে আগে থেকে Linear Static চালাতে হয় না।
 *
 * numModes: কতগুলো সবচেয়ে-ছোট-|λ|-এর mode ফেরত দিতে হবে (ডিফল্ট 6 —
 * buckling analysis এ সাধারণত শুধু প্রথম কয়েকটা critical mode
 * প্রাসঙ্গিক হয়)।
 *
 * ⚠️ একটা মডেলিং সতর্কতা (numerical verification এর সময় ধরা পড়েছে,
 * test_buckling_analysis.cpp এ documented): যদি একটা member এর strong-
 * axis (Ixx) ও weak-axis (Iyy) moment of inertia কাছাকাছি মানের হয়
 * (প্রায় square/circular section), এবং out-of-plane দিকে কোনো
 * restraint না থাকে, তাহলে সবচেয়ে critical (সবচেয়ে ছোট |λ|) mode
 * out-of-plane (weak-axis) buckling হতে পারে, in-plane (strong-axis,
 * সাধারণত যেটা প্রকৌশলী আশা করেন) না। এটা কোনো bug না — physically
 * সঠিক ফলাফল (সত্যিকারের সবচেয়ে দুর্বল দিকই buckle করবে প্রথমে) —
 * কিন্তু ব্যবহারকারীর প্রত্যাশার সাথে না মিললে বিভ্রান্তিকর হতে পারে,
 * তাই বেশ কয়েকটা mode (numModes>1) ফেরত দেখে কোনটা in-plane vs
 * out-of-plane তা mode shape থেকে যাচাই করা ভালো অভ্যাস।
 */
BucklingAnalysisResult solveLinearBuckling(const AnalysisModel& model, int numModes = 6);

/**
 * P-Delta (Second-Order/Geometric Nonlinear Static) Analysis —
 * single-iteration পদ্ধতি (types.h এর PDeltaAnalysisResult docstring
 * এ তত্ত্ব ব্যাখ্যা করা)।
 *
 * পদ্ধতি:
 *   ১. Linear Static সমাধান — প্রথম-ক্রম axial force বের করতে
 *      (solveLinearBuckling এর মতোই একই প্রাথমিক ধাপ)।
 *   ২. প্রতিটা element এর axial force দিয়ে geometric stiffness (Kg)
 *      বানিয়ে global Kg assemble করা।
 *   ৩. Modified system (K + Kg) দিয়ে পুনরায় সমাধান — একই load vector
 *      (F) ব্যবহার করে, যেহেতু load pattern পরিবর্তন হয়নি, শুধু
 *      stiffness।
 *
 * Boundary condition: এখানে Linear Static এর মতোই penalty method
 * ব্যবহার করা হয়েছে (Modal/Buckling এর elimination method না) —
 * কারণ P-Delta একটা static solve (K+Kg)U=F, eigenvalue problem না,
 * তাই penalty method এর numerical issue (যা eigenvalue solve এ
 * সমস্যা করত) এখানে প্রযোজ্য না। এটা solveLinearStatic() এর সাথে
 * consistent রাখে (একই boundary-condition পদ্ধতি ব্যবহার করলে
 * first-order ও second-order ফলাফল সরাসরি তুলনাযোগ্য থাকে)।
 *
 * ⚠️ গুরুত্বপূর্ণ সতর্কতা (numerical testing এ ধরা পড়েছে,
 * test_pdelta_analysis.cpp এ documented): "(K + Kg) সবসময় positive-
 * definite না হলে ফাংশন পরিষ্কারভাবে ব্যর্থ হবে" — এই ধারণা আংশিক
 * ভুল। যেহেতু penalty method ব্যবহার করা হয়েছে (উপরে ব্যাখ্যা করা),
 * boundary-condition DOF এ কৃত্রিম বিশাল stiffness (1e12) যোগ হয়,
 * যা matrix decomposition কে numerically "সফল" রাখতে পারে এমনকি
 * load critical buckling load ছাড়িয়ে গেলেও (কারণ penalty term
 * diagonal-dominant থাকে, প্রকৃত structural ill-conditioning সত্ত্বেও)।
 * সরাসরি পরীক্ষায় দেখা গেছে: load fraction 99% Pcr এ amplification
 * ratio ঠিকভাবে বড় হয় (~100x, তাত্ত্বিক 1/(1-0.99)=100 এর কাছাকাছি),
 * কিন্তু load 150% Pcr (critical load পার হয়ে যাওয়া) এ ratio
 * বিভ্রান্তিকরভাবে *ছোট* হয়ে যায় (~2x) — কারণ 1/(1-x) সূত্র x>1 এ
 * negative হয়ে "wrap around" করে। অর্থাৎ:
 *   - load critical load এর কাছাকাছি কিন্তু না-পেরোলে: solve সফল হয়,
 *     amplification ratio (result.maxDisplacementAmplificationRatio)
 *     বড় হতে থাকে — এটাই নির্ভরযোগ্য সতর্কতা সংকেত।
 *   - load critical load সম্পূর্ণ পেরিয়ে গেলে: solve *প্রায়ই* এখনো
 *     সফল হবে (penalty method এর কারণে), কিন্তু ফলাফল physically
 *     অর্থহীন এবং amplification ratio ছোট দেখাতে পারে (বিভ্রান্তিকর,
 *     সতর্কতা signal হিসেবে নির্ভরযোগ্য না)।
 * তাই caller এর উচিত সবসময় আলাদাভাবে solveLinearBuckling() চালিয়ে
 * critical load factor (λ) যাচাই করা এবং প্রয়োগকৃত load যদি λ এর
 * তুলনায় খুব কাছাকাছি বা বেশি হয় (যেমন load/Pcr > 0.9), P-Delta
 * ফলাফলে সতর্কতার সাথে নির্ভর করা — শুধু solvePDelta() এর success/
 * failure এর উপর নির্ভর না করে।
 */
PDeltaAnalysisResult solvePDelta(const AnalysisModel& model);

/**
 * Response Spectrum Analysis (RSA) — একটা design response spectrum
 * (period vs spectral acceleration, ইউজার/BNBC-প্রদত্ত tabulated curve)
 * দিয়ে ভবনের peak seismic response (displacement, element force, base
 * shear) বের করে, পূর্ণ time history জানা ছাড়াই।
 *
 * পদ্ধতি (standard modal response spectrum method, SAP2000/ETABS এর
 * ভিত্তি — Chopra, "Dynamics of Structures", Chapter 13):
 *   ১. Modal Analysis চালানো (internally solveModalAnalysis() কল করে)
 *      — natural frequency/period ও mass-normalized mode shape (φᵢ,
 *      φᵢᵀMφᵢ=1) বের করা।
 *   ২. প্রতিটা mode i এর জন্য modal participation factor বের করা:
 *
 *        Γᵢ = φᵢᵀ M ι / (φᵢᵀ M φᵢ) = φᵢᵀ M ι   [যেহেতু φᵢᵀMφᵢ=1]
 *
 *      যেখানে ι (influence vector) হলো একটা vector যাতে ground-motion-
 *      direction এর translational DOF এ 1, বাকি সব DOF এ 0 (যেমন X-
 *      direction ground motion হলে প্রতিটা node এর ux DOF এ 1, বাকি
 *      uy,uz,rx,ry,rz এ 0)। Γᵢ measure করে mode i কতটা "excited" হয়
 *      সেই direction এর ground motion দ্বারা।
 *   ৩. Effective modal mass: mᵢ* = Γᵢ² (φᵢᵀMφᵢ) = Γᵢ² [যেহেতু φᵢᵀMφᵢ=1
 *      mass-normalization থেকে] — কিন্তু এটা conceptually "কতটা mass
 *      এই mode participate করছে" বোঝায়, সব mode এর effective mass
 *      যোগফল ভাগ মোট mass ≈ 1.0 হওয়া উচিত যথেষ্ট mode নিলে (mass
 *      participation check, নিচে দেখুন)।
 *   ৪. প্রতিটা mode এর peak modal displacement (generalized coordinate):
 *        Dᵢ = Γᵢ · Sa(Tᵢ) / ωᵢ²   [Sa(Tᵢ) design spectrum থেকে
 *        interpolate করা, g একক থেকে m/s² এ রূপান্তরিত]
 *      এবং peak physical response (mode i এর একার অবদান):
 *        uᵢ = Dᵢ · φᵢ   (nodal displacement vector)
 *   ৫. CQC (Complete Quadratic Combination) দিয়ে সব mode এর peak
 *      response একত্র করা:
 *        R = sqrt( ΣᵢΣⱼ ρᵢⱼ Rᵢ Rⱼ )
 *      যেখানে Rᵢ = mode i এর peak response (যেকোনো quantity —
 *      displacement DOF, element force ইত্যাদি — একই সূত্র প্রযোজ্য),
 *      এবং ρᵢⱼ হলো modal correlation coefficient (Der Kiureghian, 1981
 *      formula, damping ratio ζ এর function) — কাছাকাছি frequency
 *      mode এর মধ্যে (β=ωⱼ/ωᵢ≈1) correlation বেশি (ρ→1), দূরের
 *      frequency mode এর মধ্যে কম (ρ→0, তখন CQC≈SRSS)। এই সূত্র নিচে
 *      computeCQCCorrelationCoefficient() এ implement করা।
 *
 * ⚠️ গুরুত্বপূর্ণ সীমাবদ্ধতা: element end force CQC combination সরাসরি
 * প্রতিটা mode এর element end force বের করে (Linear Static এর মতো,
 * mode shape কে "displacement" ধরে stiffness দিয়ে force বের করা) তারপর
 * সেই per-mode force গুলো CQC করা হয় — এটা standard practice, কিন্তু
 * এর অর্থ elementEndForces এর প্রতিটা component (axial, shear, moment)
 * *আলাদাভাবে* CQC combine হয়, একসাথে না। ফলে combined element end
 * force vector নিজে কোনো একটা বাস্তব equilibrium state represent করে
 * না (axial আর moment আলাদা mode থেকে আসতে পারে) — এটা design envelope
 * হিসেবে ব্যবহারযোগ্য (প্রতিটা component এর independent worst-case),
 * কিন্তু একটা একক coherent deformed shape না। এটা RSA method এর সহজাত
 * সীমাবদ্ধতা (সব commercial software এও একই আচরণ), বাগ না।
 *
 * dampingRatio: সব mode এর জন্য একই damping ratio ধরা হয়েছে (সাধারণত
 * concrete structure এ 0.05 = 5%, ASCE 7/BNBC এর common default) —
 * mode-specific ভিন্ন damping (Rayleigh damping ইত্যাদি) এই সংস্করণে
 * সমর্থিত না।
 *
 * directionDOF: 0=X, 1=Y, 2=Z (translational DOF index, node-local 6-
 * DOF ordering এর প্রথম তিনটার একটা) — কোন দিকে ground motion প্রযোজ্য
 * তা নির্দিষ্ট করে। BNBC 2020 সাধারণত X ও Z (দুই horizontal direction)
 * উভয় দিকে আলাদা RSA চালানো দাবি করে (vertical, Y, সাধারণত আলাদাভাবে
 * প্রয়োজন না হলে বাদ দেওয়া হয় — এই সংস্করণে caller প্রয়োজনে Y ও পাস
 * করতে পারে, কোনো hardcoded restriction নেই)।
 */
ResponseSpectrumAnalysisResult solveResponseSpectrum(
    const AnalysisModel& model,
    const std::vector<ResponseSpectrumPoint>& spectrum,
    int directionDOF,
    double dampingRatio = 0.05,
    int numModes = 12
);

/**
 * Nonlinear Static Analysis (Concentrated Plastic Hinge পদ্ধতি) —
 * material yielding (moment hinge formation) ধরে load-displacement
 * সম্পর্ক বের করে, Linear Static এর মতো একটা single matrix solve না
 * করে বরং incremental-iterative পদ্ধতিতে।
 *
 * অ্যালগরিদম (Load-Control Newton-Raphson, standard nonlinear FE
 * পদ্ধতি — Bathe, "Finite Element Procedures", Chapter 6 এর সাধারণ
 * কাঠামো অনুসরণ করে):
 *
 *   ১. সম্পূর্ণ applied load কে numLoadSteps সমান increment এ ভাগ করা
 *      (উদাহরণ: numLoadSteps=10 হলে, প্রতি step এ 10% load যোগ হয়)।
 *      এটা কেন দরকার: yielding একটা load-history-নির্ভর, non-smooth
 *      ঘটনা (একটা hinge হঠাৎ "on" হয়ে যায় যখন moment capacity
 *      ছাড়িয়ে যায়) — পুরো load একবারে প্রয়োগ করলে Newton-Raphson
 *      iteration ভুল hinge sequence "ধরে ফেলতে" পারে (দুইটা hinge
 *      একসাথে ভুল ক্রমে yield করলে ভিন্ন ফলাফল আসতে পারে বাস্তব
 *      progressive yielding থেকে)। Incremental loading এই sequence
 *      সঠিকভাবে ট্র্যাক করে।
 *
 *   ২. প্রতিটা load step এ, Newton-Raphson iteration:
 *      ক. বর্তমান hinge state (কোন কোন প্রান্ত ইতিমধ্যে yielded)
 *         অনুযায়ী tangent stiffness matrix (K_T) assemble করা —
 *         yielded প্রান্তে applyEndReleases() দিয়ে সেই প্রান্তের
 *         bending DOF release করা হয় (existing pin-release
 *         infrastructure পুনর্ব্যবহার, নিচে সীমাবদ্ধতায় বিস্তারিত)।
 *      খ. K_T · ΔU = R (residual force) সমাধান করে displacement
 *         increment বের করা, cumulative displacement আপডেট করা।
 *      গ. নতুন displacement দিয়ে প্রতিটা element এর প্রান্তের moment
 *         recompute করা, কোনো non-yielded hinge-assigned প্রান্ত তার
 *         yield capacity অতিক্রম করেছে কিনা চেক করা — করলে সেই
 *         প্রান্তকে yielded চিহ্নিত করে পরবর্তী iteration এ tangent
 *         stiffness rebuild করা (ধাপ ক এ ফিরে)।
 *      ঘ. residual force norm একটা tolerance এর নিচে না আসা পর্যন্ত
 *         (খ)-(গ) পুনরাবৃত্তি, অথবা maxIterationsPerStep এ পৌঁছালে
 *         non-convergence হিসেবে চিহ্নিত (কিন্তু crash না করে
 *         সর্বোত্তম প্রাপ্ত ফলাফল সহ এগিয়ে যাওয়া, caller কে
 *         result.converged=false দিয়ে জানানো)।
 *
 *   ৩. সব load step সম্পন্ন হলে, চূড়ান্ত displacement/force/hinge
 *      state ফেরত দেওয়া।
 *
 * Elastic-Perfectly-Plastic hinge model: একটা hinge-assigned প্রান্তের
 * moment |M| তার yield capacity My ছাড়িয়ে গেলে, সেই প্রান্তে bending
 * DOF release করা হয় (applyEndReleases() দিয়ে, উভয় axis — major ও
 * minor — একসাথে, নিচে সীমাবদ্ধতায় ব্যাখ্যা করা)। এর ফলে সেই প্রান্তের
 * moment "capacity তে আটকে" থাকে (perfectly plastic — কোনো strain
 * hardening/softening নেই) ও বাকি structure যেন সেই প্রান্তে একটা
 * internal pin আছে এমনভাবে আচরণ করে।
 *
 * ⚠️ গুরুত্বপূর্ণ সীমাবদ্ধতা (types.h এর NonlinearStaticAnalysisResult
 * docstring এও উল্লেখ করা):
 *   - শুধু major-axis moment (Mz, section.yieldMomentMzKNm) এর
 *     বিপরীতে yield চেক করা হয় — কিন্তু hinge trigger হলে
 *     applyEndReleases() উভয় axis (Mz ও My) release করে, কারণ
 *     existing pin-release infrastructure single-axis release
 *     সমর্থন করে না। এটা একটা conservative approximation (My তে
 *     capacity থাকলেও তা "হারানো" ধরা হচ্ছে) — ভবিষ্যতে single-axis
 *     release যোগ হলে এই সীমাবদ্ধতা দূর হবে। section.yieldMomentMyKNm
 *     field রাখা আছে ভবিষ্যতের জন্য, এই সংস্করণে সক্রিয়ভাবে ব্যবহৃত
 *     হয় না।
 *   - Axial-moment interaction (P-M diagram) নেই — column hinge এর
 *     capacity axial force নির্বিশেষে constant ধরা হয়।
 *   - Unloading path সমর্থিত না (monotonic loading assumed) — একবার
 *     yield হলে hinge সেই load step থেকে বাকি সব step এ released
 *     থাকে, load কমলেও (বর্তমান load-control পদ্ধতিতে load সবসময়
 *     বাড়ে, কমার সুযোগ নেই, তাই এটা practical সীমাবদ্ধতা না এই
 *     সংস্করণে, কিন্তু ভবিষ্যতে displacement-control বা cyclic
 *     loading যোগ হলে গুরুত্বপূর্ণ হয়ে উঠবে)।
 *   - Shear/axial hinge নেই, শুধু moment hinge।
 */
NonlinearStaticAnalysisResult solveNonlinearStatic(
    const AnalysisModel& model,
    int numLoadSteps = 10,
    int maxIterationsPerStep = 30,
    double convergenceTolerance = 1e-4
);

/**
 * Pushover Analysis — একটা fixed lateral load pattern (model.loads এ
 * caller-প্রদত্ত) ধীরে ধীরে scale-up করে push করা হয় (single load
 * pattern, একাধিক load case না) যতক্ষণ না control node তার target
 * displacement এ পৌঁছায় অথবা structure collapse করে (tangent
 * stiffness singular), যেটা আগে ঘটে। প্রতিটা successful push step এ
 * base shear ও control-node displacement capacityCurve এ যোগ হয় —
 * এটাই pushover এর signature output (classic base-shear-vs-roof-
 * displacement curve, FEMA 356/ASCE 41 এর performance-based design এর
 * ভিত্তি)।
 *
 * তাত্ত্বিক ভিত্তি ও hinge model solveNonlinearStatic() এর সাথে অভিন্ন
 * (একই docstring এর ব্যাখ্যা প্রযোজ্য — Concentrated Plastic Hinge,
 * elastic-perfectly-plastic, Load-Control Newton-Raphson, same
 * সীমাবদ্ধতা: শুধু Mz hinge, hinge trigger হলে উভয় axis release, কোনো
 * P-M interaction, কোনো unloading path)।
 *
 * solveNonlinearStatic() থেকে মূল পার্থক্য: numLoadSteps দিয়ে fixed-
 * count বিভাজনের বদলে push adaptively চলতে থাকে (displacement-target
 * বা collapse-এ থামা পর্যন্ত), এবং প্রতিটা step এ capacity curve point
 * capture হয়। ভবিষ্যতে gravity-load-first-then-lateral-push (বাস্তব
 * pushover practice) যোগ করা যেতে পারে — এই সংস্করণে শুধু single
 * lateral pattern push করা হয় (সরলীকরণ, নিচে সীমাবদ্ধতায় বলা আছে)।
 *
 * controlNodeIndex/controlDOF: capacity curve এর displacement অক্ষ
 * কোন node/DOF থেকে পড়া হবে (সাধারণত roof-level node, push
 * direction এর translational DOF: 0=ux, 1=uy, 2=uz)।
 *
 * loadStepIncrement: প্রতিটা push step এ load pattern কতটুকু বাড়বে
 * (fraction of full model.loads pattern, যেমন 0.02 মানে প্রতি step এ
 * ২%)। ছোট মান মানে বেশি নিখুঁত capacity curve কিন্তু বেশি step/সময়।
 *
 * maxPushSteps: safety limit — targetControlDisplacementM এ কখনো না
 * পৌঁছালে (বা loadStepIncrement এত ছোট যে বহু step লাগবে) অসীম loop
 * এড়াতে।
 *
 * ⚠️ সীমাবদ্ধতা: gravity load সমর্থিত না আলাদাভাবে — শুধু single
 * lateral load pattern push করা হয় model.loads থেকে। বাস্তব pushover
 * practice এ সাধারণত gravity load (dead+partial live) আগে থেকে ধ্রুবক
 * রাখা হয়, তারপর lateral load push করা হয় — এই সংস্করণে caller চাইলে
 * gravity+lateral উভয়ই model.loads এ একসাথে দিয়ে সেই combined pattern
 * push করতে পারে (approximation, gravity ও lateral আলাদা রাখা যায় না
 * এই সংস্করণে)।
 */
PushoverAnalysisResult solvePushover(
    const AnalysisModel& model,
    int controlNodeIndex,
    int controlDOF,
    double targetControlDisplacementM,
    double loadStepIncrement = 0.02,
    int maxPushSteps = 200,
    int maxIterationsPerStep = 30,
    double convergenceTolerance = 1e-4
);

} // namespace civilos
