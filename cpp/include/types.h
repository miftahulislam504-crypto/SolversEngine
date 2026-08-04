#pragma once

/**
 * FE Solver Core Types
 * =====================
 * এই ফাইলের struct গুলো frontend-এর TypeScript টাইপের (src/lib/types/
 * element.ts, material.ts, section.ts) সাথে সরাসরি সামঞ্জস্যপূর্ণ,
 * কারণ FastAPI layer JSON payload থেকে এই struct-এ deserialize করবে।
 *
 * একক (Units) — সততার সাথে স্পষ্ট করা জরুরি, কারণ ভুল একক দিয়ে পুরো
 * সলভার নীরবে ভুল ফলাফল দেবে:
 *   - দৈর্ঘ্য/কোঅর্ডিনেট: মিটার (m) — frontend viewport-এর সাথে সামঞ্জস্যপূর্ণ
 *   - Force: kilonewton (kN)
 *   - Section properties (Area, I, J): মিটার⁴/মিটার² এককে রূপান্তরিত
 *     করে এখানে আসে (frontend mm এককে হিসাব করে, solver-এ পাঠানোর
 *     আগে m এককে convert করতে হবে — এই conversion দায়িত্ব FastAPI
 *     layer-এর, C++ কোড ধরে নেয় ইনপুট ইতিমধ্যে SI (m) এককে আছে)
 *   - Elastic Modulus: kilonewton per square meter (kN/m²) — material.ts
 *     এ MPa তে হিসাব হয়, MPa × 1000 = kN/m² (রূপান্তর নিচে ব্যাখ্যা করা)
 */

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <array>

namespace civilos {

struct Node3D {
    std::string nodeId;
    double x, y, z; // meters
};

/**
 * Section geometric properties — frontend এর SectionProperties
 * (src/lib/types/section.ts, computeSectionProperties এর রিটার্ন
 * টাইপ) এর সাথে সরাসরি মিল, কিন্তু এককে ভিন্ন (frontend mm/mm⁴,
 * এখানে m/m⁴ — FastAPI layer conversion করবে)।
 */
struct SectionProperties {
    double area;      // m²
    double ixx;       // m⁴ — strong axis moment of inertia
    double iyy;       // m⁴ — weak axis moment of inertia
    double j;          // m⁴ — torsional constant

    // Concentrated Plastic Hinge (Nonlinear Static Analysis) এর জন্য —
    // ঐচ্ছিক yield moment capacity, kN·m একক। 0.0 (ডিফল্ট) মানে "এই
    // section এর জন্য hinge capacity নির্দিষ্ট করা হয়নি" — solveNonlinearStatic()
    // তখন এই section ব্যবহারকারী সব element কে elastic (কখনো yield না
    // করা) ধরবে। মান design-code-specific হিসাব (যেমন steel: Zx·Fy,
    // concrete: ACI 318 nominal moment capacity Mn) এর দায়িত্ব caller
    // এর (frontend/Design Engine স্তরে) — solver নিজে fy বা rebar
    // থেকে recompute করে না, শুধু একটা প্রদত্ত capacity ব্যবহার করে
    // plastic hinge behavior simulate করে।
    //
    // strong-axis (Mz, ixx এর সাথে সংশ্লিষ্ট) ও weak-axis (My, iyy এর
    // সাথে সংশ্লিষ্ট) হিসেবে আলাদা রাখা হয়েছে কারণ বাস্তব section এ এই
    // দুই মান ভিন্ন (rectangular/I-section এ significant পার্থক্য)।
    double yieldMomentMzKNm = 0.0;
    double yieldMomentMyKNm = 0.0;
};

struct MaterialProperties {
    double elasticModulus; // E, kN/m²
    double shearModulus;   // G, kN/m²
    double density;        // ρ, mass density — ton/m³ (কেন ton/m³: force এককে kN
                            // ব্যবহার করা হচ্ছে (ইউরোপীয়/SI ইঞ্জিনিয়ারিং কনভেনশন),
                            // তাই সামঞ্জস্যপূর্ণ mass একক হলো "tonne" (1000 kg) —
                            // কারণ F=ma, kN = tonne·(m/s²), সরাসরি kg ব্যবহার করলে
                            // মিশ্র একক (kN বনাম kg) সমীকরণে ভুল আনবে। frontend এর
                            // unitWeight (kN/m³, বাংলাদেশি/ACI কনভেনশনে ব্যবহৃত)
                            // থেকে density = unitWeight / g (g=9.81 m/s²) সূত্রে
                            // রূপান্তর FastAPI layer এ হয় (model_conversion.py)।
    double poissonsRatio;  // ν — শুধু ShellElement এ ব্যবহৃত (plane-stress ও plate
                            // bending constitutive matrix এ দরকার); FrameElement
                            // এর material এ এই field সাধারণত অব্যবহৃত থাকে (E,G
                            // থেকে বেশিরভাগ frame calculation চলে) কিন্তু struct
                            // reuse করার জন্য এখানে রাখা হয়েছে, আলাদা material
                            // struct দুইবার সংজ্ঞায়িত না করে।
};

/**
 * একটা 3D Frame Element (Beam/Column/Brace/Pile) — frontend এর
 * LineElement (src/lib/types/element.ts) এর সাথে সামঞ্জস্যপূর্ণ।
 * connectionType বহন করা হচ্ছে যাতে ভবিষ্যতে pin-ended (truss-এর
 * মতো) element এর জন্য partial stiffness release প্রয়োগ করা যায়
 * (Phase 4a তে এখনো এটা প্রয়োগ করা হয়নি — নিচে বিস্তারিত কারণ)।
 */
struct FrameElement {
    std::string elementId;
    int startNodeIndex; // globalNodes ভেক্টরে index, node ID string না — পারফরম্যান্সের জন্য (repeated string lookup এড়াতে)
    int endNodeIndex;
    SectionProperties section;
    MaterialProperties material;
    std::string connectionType; // "moment" অথবা "pin"

    // Concentrated Plastic Hinge assignment — কোন প্রান্তে hinge বসবে।
    // ডিফল্ট false মানে সেই প্রান্তে কোনো hinge নেই (সবসময় elastic,
    // Linear Static এর মতোই আচরণ করবে সেই প্রান্তে)। উভয় false হলে,
    // অথবা section.yieldMomentMzKNm/yieldMomentMyKNm উভয়ই 0.0 হলে,
    // সম্পূর্ণ element effectively elastic (solveNonlinearStatic()
    // তখন এই element এর জন্য কোনো iteration overhead ছাড়াই Linear
    // Static এর মতো আচরণ করবে)।
    bool hingeAtStart = false;
    bool hingeAtEnd = false;
};

/**
 * একটা 4-node Quadrilateral Shell Element (Slab/Wall) — frontend এর
 * AreaElement (src/lib/types/element.ts) এর সাথে সামঞ্জস্যপূর্ণ, mesh
 * discretization এর পর একটা quad panel।
 *
 * "Shell" মানে membrane (in-plane, axial+shear — plane stress) ও
 * plate bending (out-of-plane, flexure) দুটোই একসাথে ধরা — Slab
 * (মূলত out-of-plane load, gravity) এবং Wall (মূলত in-plane load,
 * lateral shear) উভয়ই একই element formulation দিয়ে মডেল করা যায়,
 * যেহেতু বাস্তবে দুটো effect আলাদা করে বলা কঠিন (একটা shear wall এ
 * out-of-plane bending-ও থাকতে পারে সামান্য, ইত্যাদি) — একটা পূর্ণ
 * shell formulation universal ও সবচেয়ে সম্পূর্ণ সমাধান।
 *
 * Node ordering: counter-clockwise (viewed from local +z, normal
 * direction) — nodeIndices[0..3] যথাক্রমে corner 1,2,3,4। এই ক্রম
 * ভুল হলে element এর normal direction উল্টে যাবে এবং area ঋণাত্মক
 * (বা shape function জ্যামিতিকভাবে ভুল) হয়ে যাবে — mesh generation
 * কোডে (analysis_orchestration.py) এই ক্রম বজায় রাখা একটা critical
 * invariant।
 *
 * thickness: slab/wall এর পুরুত্ব (m) — membrane stiffness thickness
 * এর সমানুপাতিক (t), plate bending stiffness t³ এর সমানুপাতিক (thin
 * plate theory অনুযায়ী, thickness³/12 factor সহ স্ট্যান্ডার্ড
 * plate-bending stiffness সূত্রে)।
 */
struct ShellElement {
    std::string elementId;
    std::array<int, 4> nodeIndices; // counter-clockwise, উপরে ব্যাখ্যা করা
    double thickness; // m
    MaterialProperties material; // elasticModulus, poissonsRatio এখানে ব্যবহৃত (shearModulus/density ঐচ্ছিক প্রাসঙ্গিকতা)
};

/**
 * Boundary Condition — কোন node-এর কোন DOF restrained (এখন পর্যন্ত
 * শুধু support condition, elastic spring support Phase 4 এর পরের
 * ধাপে আসবে যখন Spring/Damper element যোগ হবে)।
 */
struct BoundaryCondition {
    int nodeIndex;
    bool restrainX, restrainY, restrainZ;
    bool restrainRx, restrainRy, restrainRz;
};

/**
 * Nodal Load — একটা node-এ সরাসরি প্রযুক্ত force/moment। Uniform
 * line load কে equivalent nodal load এ রূপান্তর করার দায়িত্ব একটা
 * আলাদা ফাংশনের (loadConversion.h/cpp), এই struct শুধু চূড়ান্ত,
 * ইতিমধ্যে-nodal-এ-রূপান্তরিত load ধরে রাখে।
 */
struct NodalLoad {
    int nodeIndex;
    double fx, fy, fz;    // kN
    double mx, my, mz;    // kN·m
};

struct AnalysisModel {
    std::vector<Node3D> nodes;
    std::vector<FrameElement> elements;
    std::vector<ShellElement> shellElements; // Slab/Wall, mesh-discretized quad panel
    std::vector<BoundaryCondition> boundaryConditions;
    std::vector<NodalLoad> loads;
};

struct AnalysisResult {
    bool success;
    std::string errorMessage; // success=false হলে কারণ (যেমন "unstable structure — singular stiffness matrix")
    std::vector<Eigen::VectorXd> nodalDisplacements; // প্রতিটা node এর জন্য ৬-DOF displacement vector [ux,uy,uz,rx,ry,rz]
    std::vector<Eigen::VectorXd> elementEndForces;   // প্রতিটা element এর জন্য ১২-DOF end force vector (local coordinate এ)

    /**
     * Phase 10n — Support reaction forces।
     *
     * শুধু boundaryConditions এ যেসব node আছে সেগুলোর জন্য পপুলেট হয়
     * (unsupported node এ entry থাকে না — reactionForces.size() ==
     * boundaryConditions.size(), positionally একই ক্রমে)। প্রতিটা entry
     * একটা ৬-DOF vector [Fx,Fy,Fz,Mx,My,Mz] (kN, kN·m), global
     * coordinate এ (local axis transform এর দরকার নেই, কারণ এটা সরাসরি
     * global boundary condition DOF এর সাথে সামঞ্জস্যপূর্ণ)।
     *
     * গণনা পদ্ধতি: আমাদের boundary condition penalty method দিয়ে
     * apply হয় (দেখুন applyBoundaryConditions — restrained DOF এ একটা
     * বিশাল PENALTY_FACTOR স্টিফনেস যোগ করা হয়, row/column elimination
     * না)। Global equilibrium (K_structure + K_penalty) U = F দিয়ে
     * reaction = −PENALTY_FACTOR × U[dof] (ঋণাত্মক sign আবশ্যক —
     * K_structure U = F − K_penalty U, ডান পাশের −K_penalty U অংশটাই
     * প্রকৃত reaction)। hand-calculated simply-supported-beam টেস্ট
     * কেসে (test_solver_integration.cpp) verify করা হয়েছে: downward
     * center load এ উভয় support এ reaction upward/positive আসে ঋণাত্মক
     * sign সহ, যা physically সঠিক (positive sign দিলে ভুলভাবে downward
     * reaction আসত, যা ধরা পড়েছিল এই hand-verification এই ধাপে)।
     */
    std::vector<Eigen::VectorXd> reactionForces;
};

/**
 * Modal Analysis এর ফলাফল — natural frequency ও mode shape।
 *
 * mode shape গুলো global DOF (N = 6 × node সংখ্যা) স্কেলে সংরক্ষিত,
 * restrained DOF সহ (সেই position গুলোতে সবসময় 0, যেহেতু boundary
 * condition অনুযায়ী সেই DOF এ কোনো displacement সম্ভব না) — এটা
 * nodalDisplacements এর সাথে সামঞ্জস্যপূর্ণ রাখা হয়েছে (একই indexing
 * scheme, node-ভিত্তিক ৬-DOF ভেক্টরের তালিকা), যাতে frontend এই দুই
 * ধরনের ফলাফল একই ভাবে render করতে পারে।
 *
 * Mode shape গুলো mass-normalized (φᵀMφ = 1 প্রতিটা mode-এর জন্য) —
 * এটা structural dynamics-এ একটা standard normalization convention
 * (SAP2000/ETABS সহ প্রায় সব commercial software এটা ব্যবহার করে),
 * যা পরবর্তী ধাপে (Response Spectrum Analysis, modal combination)
 * প্রয়োজন হবে।
 */
struct ModalAnalysisResult {
    bool success;
    std::string errorMessage;
    std::vector<double> naturalFrequenciesHz;         // প্রতিটা mode-এর natural frequency, Hz (ছোট থেকে বড় ক্রমে)
    std::vector<double> angularFrequenciesRadPerSec;  // ω = 2πf, rad/s — মূল eigenvalue (ω²) থেকে সরাসরি, রাউন্ডিং এড়াতে আলাদা রাখা হয়েছে
    std::vector<std::vector<Eigen::VectorXd>> modeShapes; // modeShapes[i] = i-তম mode এর জন্য প্রতিটা node এর ৬-DOF shape vector
    int numModesComputed = 0;
};

/**
 * Linear (Eigenvalue) Buckling Analysis এর ফলাফল — critical load factor
 * ও buckling mode shape।
 *
 * তত্ত্ব সংক্ষেপে: একটা structure এ axial force (compression/tension)
 * থাকলে তার effective bending stiffness পরিবর্তিত হয় — compression
 * stiffness কমায় (P-Delta effect এর মূল কারণ), tension বাড়ায়। এই
 * প্রভাব ধরার জন্য geometric stiffness matrix (Kg) লাগে, যা element
 * এর axial force এর একটা রৈখিক function। Linear buckling সমাধান করে:
 *
 *   (K + λ·Kg) φ = 0   ⟹   K φ = -λ·Kg φ
 *
 * যেখানে λ হলো critical load factor — অর্থাৎ applied load কে λ গুণ
 * করলে structure buckle করবে (buckling load = λ × applied load
 * pattern)। এটা একটা generalized eigenvalue problem, Modal Analysis
 * এর মতোই গাণিতিক কাঠামো (K φ = ω² M φ) কিন্তু M এর বদলে -Kg।
 *
 * গুরুত্বপূর্ণ পার্থক্য Modal Analysis থেকে: এখানে একটা নির্দিষ্ট load
 * pattern দরকার (কোন load এর প্রভাবে buckling ঘটছে তা নির্দিষ্ট করতে
 * হয়) — তাই এই analysis একটা load case-নির্ভর, Modal Analysis এর মতো
 * load-independent (শুধু mass/stiffness distribution-নির্ভর) না।
 * প্রথমে সেই load pattern দিয়ে একটা Linear Static সমাধান চালিয়ে
 * প্রতিটা element এর axial force বের করা হয়, তারপর সেই axial force
 * দিয়ে Kg বানানো হয়।
 *
 * critical load factor এর ব্যাখ্যা: এখানে আমাদের axial force convention
 * হলো compression-positive (solveLinearStatic() এর raw
 * elementEndForces sign, stiffness.h এর computeLocalGeometricStiffnessMatrix()
 * docstring এ বিস্তারিত ব্যাখ্যা)। তাই:
 *   - λ > 0: প্রয়োগকৃত load pattern এ compression member এর জন্য একটা
 *     বৈধ (finite, positive) critical load multiplier আছে — λ × applied
 *     load এ structure buckle করবে। λ যত বড়, নিরাপত্তা margin তত বেশি।
 *   - λ < 1 কিন্তু > 0: প্রয়োগকৃত load ইতিমধ্যে critical buckling load
 *     ছাড়িয়ে গেছে (তাত্ত্বিকভাবে, linear buckling এর assumption এর
 *     মধ্যে)
 *   - λ < 0: বিপরীত দিকের load pattern এ (প্রয়োগকৃত load উল্টো করলে)
 *     buckle করবে — বর্তমান load direction এ প্রাসঙ্গিক না
 */
struct BucklingAnalysisResult {
    bool success;
    std::string errorMessage;
    std::vector<double> criticalLoadFactors; // λ, প্রতিটা mode এর জন্য (absolute value অনুযায়ী ছোট থেকে বড় ক্রমে)
    std::vector<std::vector<Eigen::VectorXd>> bucklingModeShapes; // প্রতিটা mode এর জন্য প্রতিটা node এর ৬-DOF shape vector
    int numModesComputed = 0;
};

/**
 * P-Delta (Second-Order/Geometric Nonlinear Static) Analysis এর ফলাফল।
 *
 * তত্ত্ব সংক্ষেপে: Linear Static Analysis first-order — অর্থাৎ axial
 * force এর কারণে bending stiffness পরিবর্তন (P-Delta effect) ধরে না,
 * শুধু undeformed geometry তে equilibrium ধরে। P-Delta Analysis সেই
 * প্রভাব অন্তর্ভুক্ত করে — একটা compression member এ deflection হলে
 * সেই deflection axial load দিয়ে গুণিত হয়ে অতিরিক্ত moment তৈরি করে
 * ("P times Delta"), যা আবার deflection বাড়ায় — একটা second-order
 * effect যা lateral-load-heavy বা slender column-heavy structure এ
 * (উঁচু বিল্ডিং, slender frame) গুরুত্বপূর্ণ।
 *
 * পদ্ধতি (এই সংস্করণ — single-iteration, ETABS/SAP2000 এর ডিফল্ট
 * পদ্ধতির অনুরূপ, solver.h এর solvePDelta() docstring এ বিস্তারিত):
 *   ১. Linear Static সমাধান করে প্রতিটা element এর axial force বের করা
 *   ২. সেই axial force দিয়ে geometric stiffness (Kg) বানানো (Linear
 *      Buckling Analysis এর মতোই, computeLocalGeometricStiffnessMatrix)
 *   ৩. Modified stiffness (K + Kg) দিয়ে পুনরায় সমাধান — এটাই P-Delta-
 *      modified displacement/force দেয়
 *
 * এটা একটা approximation (single-iteration) — পুরোপুরি converged
 * nonlinear সমাধান না (যা একাধিক iteration দাবি করত, প্রতি iteration এ
 * নতুন axial force দিয়ে Kg আপডেট করে, পরিবর্তন যথেষ্ট ছোট না হওয়া
 * পর্যন্ত)। তবে single-iteration পদ্ধতি বেশিরভাগ ব্যবহারিক ক্ষেত্রে
 * (moderate P-Delta effect, যা বেশিরভাগ building code দাবি করে) যথেষ্ট
 * নির্ভুল এবং এটাই ETABS/SAP2000 এর ডিফল্ট আচরণ।
 */
struct PDeltaAnalysisResult {
    bool success;
    std::string errorMessage;
    std::vector<Eigen::VectorXd> nodalDisplacements; // P-Delta-modified displacement, nodalDisplacements এর মতোই shape (AnalysisResult এর সাথে সঙ্গতিপূর্ণ)
    std::vector<Eigen::VectorXd> elementEndForces;   // P-Delta-modified element end force
    std::vector<double> firstOrderAxialForces;       // প্রতিটা element এর প্রাথমিক (first-order, ধাপ ১ থেকে) axial force — caller কে amplification factor বুঝতে সাহায্য করার জন্য
    double maxDisplacementAmplificationRatio = 1.0;  // সবচেয়ে বড় |P-Delta displacement| / |first-order displacement| অনুপাত, যেকোনো DOF জুড়ে (P-Delta effect এর overall মাত্রা বোঝার একটা দ্রুত সূচক)
};

/**
 * একটা Design Response Spectrum point — period (T, সেকেন্ড) বনাম
 * spectral acceleration (Sa, g একক)। solveResponseSpectrum() এই
 * point-তালিকার মধ্যে piecewise-linear interpolation করে যেকোনো mode
 * এর period এর Sa বের করে। BNBC/ASCE/Eurocode-নির্দিষ্ট নয় — spectrum
 * shape generation Python layer এ (app/response_spectrum.py), C++ শুধু
 * generic tabulated point নেয়।
 */
struct ResponseSpectrumPoint {
    double periodSec = 0.0;
    double spectralAccelerationG = 0.0;
};

/**
 * Response Spectrum Analysis এর ফলাফল — প্রতিটা mode কে independent
 * SDOF oscillator ধরে design spectrum থেকে peak response (Sa(Tᵢ)) বের
 * করা হয়, modal participation factor ও effective modal mass দিয়ে
 * scale করা হয়, তারপর CQC (Complete Quadratic Combination) দিয়ে
 * একত্র করে সামগ্রিক peak response (envelope) বানানো হয়।
 *
 * ফলাফল সবসময় ≥0 (magnitude, peak) — Linear Static এর মতো signed না,
 * কারণ প্রতিটা mode এর peak আলাদা সময়ে ঘটে, sign/direction তথ্য নেই।
 * এটা response spectrum method এর একটা সহজাত boundary, বাগ না।
 */
struct ResponseSpectrumAnalysisResult {
    bool success = false;
    std::string errorMessage;
    std::vector<Eigen::VectorXd> nodalDisplacements;
    std::vector<Eigen::VectorXd> elementEndForces;
    double baseShear = 0.0;
    double totalMassParticipationRatio = 0.0;
    std::vector<double> modalParticipationFactors;
    std::vector<double> effectiveModalMasses;
    std::vector<double> modalSpectralAccelerations;
    int numModesComputed = 0;
};

/**
 * একটা প্লাস্টিক হিঞ্জের চূড়ান্ত (converged) অবস্থা — Nonlinear Static
 * Analysis এর ফলাফলের অংশ, caller (frontend) কে দেখাতে দেয় কোন কোন
 * hinge yield করেছে (performance-based design/damage assessment এর
 * প্রাথমিক ভিত্তি — FEMA 356/ASCE 41 এর hinge status রিপোর্টিং এর
 * সাথে ধারণাগতভাবে সাদৃশ্যপূর্ণ, যদিও এই সংস্করণে শুধু "yielded/not"
 * বাইনারি স্ট্যাটাস, Immediate Occupancy/Life Safety/Collapse
 * Prevention এর মতো performance-level threshold এখনো নেই)।
 */
struct PlasticHingeState {
    int elementIndex = 0;
    bool isAtStartNode = true; // true হলে element এর start node প্রান্তের hinge, false হলে end node
    bool yielded = false;
    double finalMomentKNm = 0.0;   // converged moment এই hinge এ (yielded হলে ≈ yieldMoment, চিহ্ন সহ)
    double plasticRotationRad = 0.0; // yielded হলে accumulated plastic rotation, না হলে 0
};

/**
 * Nonlinear Static Analysis (Concentrated Plastic Hinge পদ্ধতি) এর
 * ফলাফল।
 *
 * তত্ত্ব সংক্ষেপে: Linear Static এ ধরে নেওয়া হয় stiffness matrix (K)
 * সম্পূর্ণ load range জুড়ে constant থাকে (material সবসময় elastic)।
 * বাস্তবে কোনো element এর প্রান্তের moment তার yield capacity
 * (section.yieldMomentMzKNm/yieldMomentMyKNm, hingeAtStart/hingeAtEnd
 * দিয়ে চিহ্নিত প্রান্তে) ছাড়িয়ে গেলে, সেই প্রান্তে একটা "plastic hinge"
 * তৈরি হয় — moment সেই capacity তে "আটকে" যায় (elastic-perfectly-
 * plastic idealization, কোনো strain hardening না), এবং সেই প্রান্তের
 * rotational stiffness কার্যত শূন্য হয়ে যায় (একটা internal pin এর
 * মতো আচরণ করে, কিন্তু moment=0 না — moment=yield capacity তে ধরে
 * রাখে)।
 *
 * এই non-constant, load-history-নির্ভর stiffness এর কারণে, সমাধান একটা
 * single matrix solve দিয়ে হয় না (Linear Static/P-Delta এর মতো) —
 * বরং load কে ছোট ছোট increment এ ভাগ করে (loadSteps), প্রতিটা
 * increment এ Newton-Raphson iteration চালিয়ে equilibrium (residual
 * force ≈ 0) না পাওয়া পর্যন্ত tangent stiffness rebuild করা হয় (নতুন
 * hinge তৈরি হলে stiffness matrix বদলে যায়)। solveNonlinearStatic()
 * এর docstring এ (solver.h) সম্পূর্ণ algorithm ব্যাখ্যা করা।
 *
 * ⚠️ গুরুত্বপূর্ণ সীমাবদ্ধতা:
 *   - শুধু major-axis moment (Mz) ভিত্তিক hinge সমর্থিত এই সংস্করণে —
 *     axial-moment interaction (P-M interaction diagram, column এর
 *     জন্য গুরুত্বপূর্ণ) এখনো নেই। প্রতিটা hinge শুধু তার নিজস্ব moment
 *     capacity এর বিপরীতে independently চেক হয়, axial force এর প্রভাব
 *     capacity তে ধরা হয় না।
 *   - Unloading/cyclic behavior সমর্থিত না (monotonic loading assumed)
 *     — একবার yield হলে, load কমলেও hinge elastic এ ফিরে আসে না এই
 *     সরলীকৃত সংস্করণে (প্রকৃত elastic-perfectly-plastic unloading
 *     path একটা ভবিষ্যৎ উন্নতি)।
 *   - Shear hinge/axial hinge নেই, শুধু moment hinge।
 */
struct NonlinearStaticAnalysisResult {
    bool success = false;
    std::string errorMessage;
    std::vector<Eigen::VectorXd> nodalDisplacements; // চূড়ান্ত (সম্পূর্ণ load প্রয়োগের পর) converged displacement
    std::vector<Eigen::VectorXd> elementEndForces;   // চূড়ান্ত converged element end force
    std::vector<PlasticHingeState> hingeStates;       // শুধু hinge-assigned প্রান্তগুলোর চূড়ান্ত অবস্থা
    int totalLoadSteps = 0;
    int totalNewtonIterations = 0; // সব load step মিলিয়ে মোট iteration সংখ্যা (convergence performance এর একটা ইঙ্গিত)
    bool converged = false;        // সব load step এ convergence tolerance এর মধ্যে পৌঁছেছে কিনা
    double maxDisplacementAmplificationRatio = 1.0; // nonlinear বনাম প্রাথমিক (প্রথম load step, near-elastic) displacement এর অনুপাত
};

/**
 * Pushover capacity curve এর একটা point — base shear (kN) বনাম control
 * node এর displacement (m), একটা নির্দিষ্ট push step এ। পুরো curve
 * (সব step এর point একসাথে) হলো pushover এর signature deliverable —
 * এটা দিয়েই structure এর overall lateral capacity, ductility, ও
 * performance point (FEMA 356/ASCE 41 এর সাথে ধারণাগতভাবে সাদৃশ্যপূর্ণ,
 * যদিও এই সংস্করণে performance-level threshold নিজে নেই) বোঝা যায়।
 */
struct PushoverCurvePoint {
    double baseShearKN = 0.0;
    double controlDisplacementM = 0.0;
    int numHingesYielded = 0; // এই point এ কতগুলো hinge yielded অবস্থায় আছে — capacity curve এর সাথে "কোন step এ কোন hinge প্রথম yield করলো" বোঝার একটা সহজ সূচক
};

/**
 * Pushover Analysis এর ফলাফল — Nonlinear Static Analysis এরই একটা
 * বিশেষ প্রয়োগ (একই Concentrated Plastic Hinge পদ্ধতি, একই Load-
 * Control Newton-Raphson iteration, solver.h এর solveNonlinearStatic()
 * docstring এ ব্যাখ্যা করা তাত্ত্বিক ভিত্তি প্রযোজ্য), কিন্তু দুটো
 * গুরুত্বপূর্ণ পার্থক্য সহ:
 *
 *   ১. Load pattern সবসময় একটা single, fixed lateral load shape হয়
 *      (সাধারণত storey-wise distributed lateral force, BNBC/ASCE 41
 *      এর triangular বা uniform pattern) — model.loads এ যা আছে তাই
 *      ব্যবহার করা হয় (caller এর দায়িত্ব সঠিক lateral pattern সেট করা,
 *      solveNonlinearStatic() এর মতোই)। প্রতিটা step এ এই পুরো pattern
 *      একটা scale factor দিয়ে বাড়ানো হয় — gravity load এই push
 *      pattern এর অংশ না (বাস্তব pushover এ gravity আগে থেকেই প্রয়োগ
 *      করা থাকে, constant রাখা হয়, তারপর lateral push শুরু হয় — এই
 *      সংস্করণে gravity load আলাদাভাবে সমর্থিত না, শুধু single lateral
 *      pattern push করা হয়, এটা একটা সরলীকরণ)।
 *   ২. numLoadSteps দিয়ে fixed load বিভাজনের বদলে, push a) একটা
 *      targetControlDisplacementM এ পৌঁছানো পর্যন্ত, বা b) structure
 *      আর কোনো additional load নিতে না পারা (tangent stiffness
 *      singular — collapse/mechanism) পর্যন্ত চলে, যেটা আগে ঘটে।
 *      প্রতিটা successful step এ capacityCurve এ একটা point যোগ হয়।
 *
 * controlNodeIndex ও controlDOF (0=ux, 1=uy, 2=uz) caller নির্দিষ্ট
 * করে — সাধারণত roof-level এর কোনো node, lateral push direction এর
 * DOF এ (roof displacement vs base shear, classic pushover curve)।
 */
struct PushoverAnalysisResult {
    bool success = false;
    std::string errorMessage;
    std::vector<PushoverCurvePoint> capacityCurve; // push শুরু (0,0) থেকে চূড়ান্ত অবস্থা পর্যন্ত, ক্রমানুসারে
    std::vector<Eigen::VectorXd> finalNodalDisplacements; // চূড়ান্ত push step এর displacement
    std::vector<Eigen::VectorXd> finalElementEndForces;
    std::vector<PlasticHingeState> finalHingeStates;
    bool reachedTargetDisplacement = false; // targetControlDisplacementM এ পৌঁছেছে কিনা (নাকি collapse/non-convergence এ থেমেছে)
    bool structureCollapsed = false;        // tangent stiffness singular হয়ে push থামলে true (mechanism/collapse এর গাণিতিক ইঙ্গিত)
    int totalPushSteps = 0;
    int totalNewtonIterations = 0;
};

} // namespace civilos
