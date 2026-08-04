# CivilOS Structural — Compute Microservice

Phase 4a: C++ FE solver (Linear Static + Modal + Linear Buckling +
P-Delta Analysis, Line Element **ও** Shell Element/Slab/Wall) এখন
পুরোপুরি ইন্টিগ্রেটেড। FastAPI দিয়ে বানানো, Cloud Run/Railway-তে
deploy হবে।

> **যাচাই নোট (সততার সাথে):** সব C++ কোড এই sandbox-এ সরাসরি কম্পাইল
> ও চালিয়ে যাচাই করা হয়েছে (g++ ম্যানুয়ালি, এবং CMake দিয়ে যা
> production Dockerfile ব্যবহার করবে) — cantilever beam, simply-supported
> beam, এবং একটা multi-element portal frame, প্রতিটার ফলাফল classical
> textbook সূত্রের সাথে মিলিয়ে (সাত সংখ্যা পর্যন্ত নির্ভুল)। পুরো
> Python↔C++ bridge (pybind11) এবং সম্পূর্ণ HTTP API (FastAPI সার্ভার
> চালিয়ে curl দিয়ে) টেস্ট করা হয়েছে। **কিন্তু Dockerfile নিজে এই
> sandbox-এ build করে দেখা যায়নি** (Docker এখানে নেই, Phase 0-এর মতোই
> সীমাবদ্ধতা) — প্রতিটা ধাপ native environment-এ ম্যানুয়ালি replicate
> করে verify করা হয়েছে (apt install cmake/eigen, venv setup, cmake
> build), কিন্তু Docker multi-stage COPY-এর প্রকৃত আচরণ যাচাই হয়নি।
> প্রথম Railway deploy-ই এই Dockerfile-এর প্রথম real build — সেটা fail
> করলে build log পাঠালে ঠিক করে দেব।

---

## Phase 4a-তে কী বসেছে

### C++ Solver (`cpp/`)
- **`include/types.h`** — Node, FrameElement, **ShellElement** (4-node quad, Slab/Wall), Material (density+poissonsRatio সহ), Section, Load, এবং ModalAnalysisResult/BucklingAnalysisResult/PDeltaAnalysisResult এর core struct
- **`include/stiffness.h` / `src/stiffness.cpp`** — 3D frame element এর 12x12 local stiffness matrix, **consistent mass matrix**, ও **geometric stiffness matrix** (সবগুলো McGuire/Przemieniecki/Cook et al. এর প্রতিষ্ঠিত সূত্র, নতুন derivation না), coordinate transformation matrix, pin-release static condensation
- **`include/shell.h` / `src/shell.cpp`** — 4-node quadrilateral shell element (24x24 local stiffness), **membrane** (bilinear Q4, plane-stress, Allman-type drilling DOF stabilization) ও **plate bending** (Mindlin-Reissner, selective reduced integration — shear locking এড়াতে) একসাথে সমন্বিত (Slab ও Wall উভয়ের জন্য একই formulation)
- **`include/solver.h` / `src/solver.cpp`** — Global sparse assembly (stiffness — frame ও shell উভয়ই, mass, geometric stiffness — শুধু frame এই মুহূর্তে), penalty-method boundary condition (Linear Static/P-Delta এ) / elimination method (Modal/Buckling এ), Sparse Cholesky (SimplicialLDLT) সমাধান, Generalized Symmetric Eigenvalue সমাধান
- **`src/bindings.cpp`** — pybind11 bridge, চারটা analysis function expose করা: `solve_linear_static` (frame+shell), `solve_modal_analysis`, `solve_linear_buckling`, `solve_pdelta` (এই তিনটা এখনো frame-only)
- **`tests/`** — আটটা test suite (stiffness matrix, full integration, pin-release, modal, buckling, P-Delta, shell element, shell integration), CMake/CTest দিয়ে চলে

### Python Integration (`app/`)
- **`model_conversion.py`** — frontend এর mm/MPa একক থেকে solver এর m/kN এককে রূপান্তর, এবং unitWeight (kN/m³) থেকে density (tonne/m³) রূপান্তর, self-verifying assertion সহ
- **`mesh_generation.py`** — Slab/Wall এর polygon vertices কে quad mesh এ রূপান্তর: **ear-clipping triangulation** (general simple polygon, concave সহ) + **triangle-to-quad conversion** (centroid+edge-midpoint দিয়ে প্রতিটা triangle কে ৩টা quad এ ভাগ)
- **`analysis_orchestration.py`** — frontend এর element/material/section/loadCase JSON থেকে node graph বানানো (coincident endpoint merge করে, mid-span point load এর জন্য element split করে, Slab/Wall মেশ করে), solver input তৈরি, ফলাফল ব্যাখ্যাযোগ্য warning সহ ফেরত দেওয়া
- **`main.py`** — `/jobs/analysis` এখন চারটা analysis_type synchronously সমাধান করে: "linear-static" (frame+shell), "modal", "buckling", "pdelta" (এই তিনটা এখনো frame-only)

### যাচাইকৃত numerical ফলাফল (এই sandbox-এ)
- Cantilever beam: tip deflection = PL³/3EI, rotation = PL²/2EI — উভয়ই exact মিল
- Simply-supported beam (2-element): center deflection = PL³/48EI, support reaction = P/2 — উভয়ই exact মিল
- Rigid-body motion → zero force (mathematical invariant) — pass
- Transformation matrix orthogonality (T^T·T = I) — pass, vertical-column edge case সহ
- Portal frame (mm/MPa raw input দিয়ে, full HTTP pipeline দিয়ে) — numerically সংগত ফলাফল, base moment = force × arm length ভেরিফাই করা হয়েছে
- Mid-span point load (element split): simply-supported beam, center ও off-center load — PL³/48EI ও Pa²b²/3EIL সূত্রের সাথে ~1e-9 relative error এ মিল
- **Modal Analysis**: cantilever beam (10-element mesh) প্রথম natural frequency Blevins এর সূত্র (ω₁=(1.875)²√(EI/mL⁴)) এর সাথে 0.34% এর মধ্যে মিলেছে; simply-supported beam প্রথম frequency (ω₁=π²√(EI/mL⁴)) সাথে মিলেছে ~7e-6 relative error এ; mode shape M-orthonormality numerically যাচাই করা হয়েছে
- **Linear Buckling**: pin-pin, fixed-free (K=2), ও fixed-fixed (K=0.5) column — Euler formula (Pcr=π²EI/(KL)²) এর সাথে 1e-5 থেকে 1e-3 relative error এ মিল
- **P-Delta**: cantilever column moment amplification factor 1/(1-P/Pcr) সূত্রের সাথে ~5.4% error এ মিল (approximate formula vs FE, প্রত্যাশিত); zero-axial-force sanity check exact মিলেছে (0% error)
- **Shell Element (Membrane)**: uniaxial tension vs. 1D bar theory (PL/EA) — ~4% এর মধ্যে (discrete-nodal-load approximation, প্রত্যাশিত)
- **Shell Element (Plate Bending)**: simply-supported square plate (8×8 mesh), center deflection vs. Timoshenko এর classical formula (α·q·a⁴/D) — ~1.9% error এ মিল
- **Mesh Generation**: rectangle, concave L-shape, vertical wall, triangle, pentagon — সব ক্ষেত্রে মোট mesh area মূল polygon area এর সাথে exact মিলেছে (area-conservation sanity check)

---

## Modal Analysis — নতুন সংযোজন

Natural frequency (Hz) ও mode shape বের করার জন্য generalized eigenvalue
problem সমাধান করা হয় (K φ = ω² M φ)।

- **Mass matrix**: Consistent mass matrix (SAP2000/ETABS এর ডিফল্ট
  পদ্ধতি, lumped mass না) — Cook et al. এর প্রতিষ্ঠিত সূত্র। Rotary
  inertia ignore করা হয়েছে (slender member এ গ্রহণযোগ্য approximation)।
- **Boundary condition**: Linear Static এর penalty method এখানে ব্যবহার
  করা হয়নি (কারণ penalty stiffness + আসল mass মিলিয়ে spurious
  অতি-উচ্চ-frequency mode তৈরি করত) — এখানে elimination method (restrained
  DOF সম্পূর্ণ বাদ) ব্যবহৃত।
- **Eigensolver**: Dense `Eigen::GeneralizedSelfAdjointEigenSolver`
  ব্যবহার করা হয়েছে (sparse eigensolver library — Spectra — এড়ানো
  হয়েছে external dependency কমাতে)। ছোট/মাঝারি মডেলে (কয়েকশো element)
  যথেষ্ট দ্রুত; হাজার+ DOF মডেলে ভবিষ্যতে sparse iterative eigensolver
  প্রয়োজন হতে পারে।
- **API**: `POST /jobs/analysis` এ `analysis_type: "modal"` পাঠান,
  ঐচ্ছিক `num_modes` (ডিফল্ট 12)। Result shape: `{analysisType,
  numModesComputed, modes: [{naturalFrequencyHz,
  angularFrequencyRadPerSec, modeShape}, ...], ...}`।
- **⚠️ জানা সীমাবদ্ধতা**: material এর `unitWeight` (frontend এর কিছু
  মান, যেমন steel=78.5 kN/m³) নিজেই g≈10 ধরে রাউন্ড করা মান — নিখুঁত
  g=9.80665 দিয়ে বিপরীত রূপান্তর করলে density-তে ~2% পার্থক্য আসে (steel
  এ 7.85 এর বদলে ~8.00 tonne/m³)। Frequency তে এর প্রভাব ~1% এর কম
  (frequency ∝ 1/√mass) — প্রকৌশল উদ্দেশ্যে গ্রহণযোগ্য, বিস্তারিত
  `model_conversion.py` এর `convert_unit_weight_to_density()` docstring এ।
- **⚠️ pin-connected element এর mass matrix**: static condensation
  (stiffness এ যা করা হয়) mass matrix এ প্রয়োগ করা হয়নি — rigid ও
  pin-connected element এর mass matrix একই সূত্রে গণনা হয়। pin-heavy
  মডেলে (অনেক Brace) এটা frequency তে সামান্য প্রভাব ফেলতে পারে।

---

## Linear Buckling Analysis — নতুন সংযোজন

Critical load factor (λ) ও buckling mode shape বের করার জন্য
geometric-stiffness-based eigenvalue problem সমাধান করা হয়।

- **পদ্ধতি**: প্রথমে দেওয়া load pattern দিয়ে Linear Static চালিয়ে
  প্রতিটা element এর axial force বের করা হয়, তারপর সেই axial force
  দিয়ে geometric stiffness (Kg) বানানো হয় (Przemieniecki সূত্র,
  bending-coupled effect শুধু — shear/torsion geometric stiffness
  অন্তর্ভুক্ত না)। তাত্ত্বিক সমীকরণ `K φ = -λ Kg φ` সরাসরি সমাধান করা
  যায় না (Kg ইনহেরেন্টলি singular, axial/torsion DOF এ geometric
  stiffness শূন্য বলে) — তাই `(-Kg) φ = μ K φ` আকারে reformulate করে
  সমাধান করা হয় (μ=1/λ), যেখানে K সবসময় positive-definite।
- **API**: `analysis_type: "buckling"`, ঐচ্ছিক `num_modes` (ডিফল্ট 6)।
  `loadCases` অখালি হতে হবে (buckling load-নির্ভর, কোন load pattern
  তা জানা আবশ্যক)।
- **⚠️ মডেলিং সতর্কতা**: যদি একটা member এর Ixx ও Iyy কাছাকাছি মানের
  হয় (প্রায় square section) এবং out-of-plane এ কোনো restraint না
  থাকে, সবচেয়ে critical mode out-of-plane (weak-axis) হতে পারে,
  in-plane (যা সাধারণত প্রত্যাশিত) না। এটা physically সঠিক, bug না —
  কিন্তু একাধিক mode দেখে যাচাই করা ভালো অভ্যাস।
- **Sign convention নোট**: `elementEndForces[e](0)` (axial force)
  compression-positive, conventional tension-positive না — এটা
  numerically যাচাই করে code comment এ স্পষ্ট করা হয়েছে (আগে একটা
  sign-related bug ধরা পড়েছিল এই আবিষ্কারের ফলে)।

---

## P-Delta Analysis — নতুন সংযোজন

Second-order (geometric nonlinear) static analysis — axial force এর
কারণে bending stiffness পরিবর্তনের প্রভাব ধরে।

- **পদ্ধতি**: single-iteration (ETABS/SAP2000-এর ডিফল্ট পদ্ধতি) —
  প্রথমে Linear Static থেকে axial force বের করে, geometric stiffness
  (Kg, buckling-এর কোডই পুনর্ব্যবহার) বানিয়ে, `(K+Kg)U=F` দিয়ে পুনরায়
  সমাধান। পুরোপুরি converged nonlinear সমাধান না (একাধিক iteration
  দাবি করত), কিন্তু বেশিরভাগ ব্যবহারিক ক্ষেত্রে যথেষ্ট নির্ভুল।
- **Boundary condition**: penalty method (Modal/Buckling এর
  elimination method না) — কারণ এটা একটা static solve, eigenvalue
  problem না।
- **API**: `analysis_type: "pdelta"`। Result এ `nodalDisplacements`,
  `elementEndForces` এর পাশাপাশি `firstOrderAxialForces` ও
  `maxDisplacementAmplificationRatio` থাকে।
- **⚠️ গুরুত্বপূর্ণ সীমাবদ্ধতা (numerical testing এ ধরা পড়েছে)**:
  penalty method এর কারণে, load critical buckling load ছাড়িয়ে গেলেও
  solve প্রায়ই এখনো "সফল" হয় (কৃত্রিম penalty stiffness matrix কে
  numerically stable রাখে)। Near-critical load (~99% Pcr) এ
  amplification ratio সঠিকভাবে বড় হয় (~100x), কিন্তু load critical
  load সম্পূর্ণ ছাড়িয়ে গেলে (>100% Pcr) ratio বিভ্রান্তিকরভাবে *ছোট*
  দেখাতে পারে (formula 1/(1-x) এর sign-wrap এর কারণে, x>1 এ)। তাই
  `main.py` amplification ratio > 3.0 হলে একটা warning যোগ করে, এবং
  caller কে সবসময় আলাদাভাবে `analysis_type: "buckling"` চালিয়ে
  critical load factor cross-check করার পরামর্শ দেওয়া হচ্ছে — শুধু
  solve success/failure এর উপর নির্ভর না করে।

---

## Shell Element (Slab/Wall) — নতুন সংযোজন

Slab, Wall, Shear Wall, Core Wall — এই ৪টা Area Element এখন Linear
Static Analysis এ solve করা যায় (Modal/Buckling/P-Delta এখনো
frame-only, shell এর mass matrix/geometric stiffness এখনো implement
করা হয়নি)।

- **Element formulation**: 4-node quadrilateral shell — **membrane**
  (bilinear Q4, plane-stress, in-plane force/deformation ধরে) ও
  **plate bending** (Mindlin-Reissner, out-of-plane bending/shear ধরে)
  একসাথে সমন্বিত একটা পূর্ণ shell (Slab মূলত out-of-plane load handle
  করে, Wall মূলত in-plane — একই formulation উভয়ের জন্য ব্যবহারযোগ্য)।
  Drilling DOF (in-plane rotation, rz) একটা ছোট penalty stiffness দিয়ে
  stabilize করা (bilinear membrane এ প্রাকৃতিক drilling stiffness
  নেই, কিন্তু ৬-DOF/node রাখতে ও coplanar mesh এ singular matrix এড়াতে
  প্রয়োজন)। Plate bending এ shear locking এড়াতে selective reduced
  integration ব্যবহার করা হয়েছে (bending: 2×2 full Gauss, shear: 1×1
  reduced Gauss)।
- **Mesh generation**: Slab/Wall এর polygon vertices (frontend থেকে,
  যেকোনো সংখ্যক vertex, concave সহ) কে quad mesh এ রূপান্তর করা হয়
  ear-clipping triangulation + triangle-to-quad conversion দিয়ে —
  `mesh_generation.py`। এটা একটা সর্বজনীন (general polygon) সমাধান,
  শুধু rectangle-এ সীমাবদ্ধ না।
- **API**: `analysis_type: "linear-static"` — কোনো নতুন analysis_type
  লাগেনি, Slab/Wall element গুলো স্বয়ংক্রিয়ভাবে shell element এ mesh
  হয়ে যায় ও frame element এর সাথে একই solve pipeline এ যোগ হয়। Mixed
  frame+shell model (যেমন shear wall + সংযুক্ত beam) সমর্থিত।
- **⚠️ গুরুত্বপূর্ণ সীমাবদ্ধতা**:
  - Mesh density uniform/adaptive না — শুধু polygon boundary অনুযায়ী
    triangulate করে প্রতিটা triangle কে ৩-quad এ ভাগ করা হয়, কোনো
    target-element-size-ভিত্তিক refinement নেই। বড় slab এ কম, uneven
    element হবে।
  - Element internal force (moment/shear per unit width, stress
    recovery) এই ধাপে রিপোর্ট করা হয় না — শুধু nodal displacement।
    এটা একটা future ধাপ (shell.cpp এ B-matrix ইতিমধ্যে আছে, কিন্তু
    Gauss-point stress থেকে nodal-force output বানানো এখনো বাকি)।
  - Shell এ সরাসরি load application (area/point load slab এ) এখনো
    সমর্থিত না — শুধু geometry/stiffness contribution ধরা হয়, shell
    কে target করা কোনো Point Load থাকলে সেটা warning সহ বাদ দেওয়া হয়।
  - Modal/Buckling/P-Delta এখনো shell ignore করে (frame-only)।
  - Hole/opening-সহ polygon (একটা slab এর মধ্যে void) সমর্থিত না।

---

## গুরুত্বপূর্ণ, honest সীমাবদ্ধতা (Phase 4a)

এই ফিচারগুলো **এখনো নেই** — ব্যবহার করার আগে জানা দরকার:

1. **শুধু Line Element** (Beam/Column/Brace/Pile) সলভ হয়। Slab/Wall/
   Shear Wall/Core Wall/Footing এই মুহূর্তে **skip** হয় (FE mesh প্রয়োজন,
   এখনো implement করা হয়নি) — response এ warning আসে কতগুলো element
   skip হয়েছে।

2. **✅ "pin" connectionType এখন প্রয়োগ করা হয়** (static condensation
   দিয়ে, `stiffness.h`/`applyEndReleases()`)। `connectionType: "pin"`
   হলে element-এর **উভয় প্রান্তের** bending moment DOF (ry, rz) release
   হয় — axial ও torsion stiffness অক্ষত থাকে। ⚠️ **সীমাবদ্ধতা: এখন শুধু
   both-end release সমর্থিত, single-end (একদিকে moment, অন্যদিকে pin)
   না** — কারণ frontend-এর `connectionType` field element-এ একটাই
   string (per-end আলাদা field নেই)। যদি আপনার মডেলে one-end-pin দরকার
   হয় (যেমন internal hinge শুধু এক প্রান্তে), সেটা এখনো সঠিকভাবে মডেল
   করা যাবে না এই Phase-এ — উভয় প্রান্তই release হয়ে যাবে।
   Numerical verification: `test_pin_release.cpp`-এ ৩টা টেস্ট (bending
   DOF শূন্য হওয়া, symmetric matrix, ও একটা two-element internal-hinge
   পূর্ণ pipeline টেস্ট যেখানে হিঞ্জ পয়েন্টে zero-moment condition
   যাচাই করা হয়েছে) — সব pass করেছে।

3. **Mid-span Point Load সঠিকভাবে হ্যান্ডল হয় না** — element split করে
   intermediate node বসানো এখনো নেই। `positionRatio` 0/1 থেকে দূরে হলে
   load কে নিকটতম প্রান্তে snap করা হয় এবং একটা 🔴 warning আসে element
   নাম-সহ। **আপাতত শুধু element-এর ঠিক শুরু বা শেষে (positionRatio 0
   বা 1) load প্রয়োগ করুন নির্ভরযোগ্য ফলাফলের জন্য।**

4. **Support condition auto-detected, manually-defined না** — Y≈0 তে
   থাকা সব node স্বয়ংক্রিয়ভাবে fully-fixed support ধরা হয়। এটা একটা
   preliminary heuristic; যদি আপনার মডেলের base অন্য elevation-এ থাকে
   বা কোনো base-level column আসলে transfer beam-এর মাধ্যমে সাপোর্টেড,
   এই heuristic ভুল ফলাফল দিতে পারে। প্রতিটা সফল সমাধানেও এই fact
   response-এ ℹ️ note হিসেবে থাকে।

5. **Uniform Line/Area Load convert হয় না** — শুধু Point Load
   ব্যবহারযোগ্য এই Phase-এ।

6. **শুধু "linear-static"** — বাকি ১৮টা analysis type (Modal, P-Delta,
   Response Spectrum, Nonlinear, Pushover, ইত্যাদি) 501 Not Implemented
   দেয়, কোনো fake ফলাফল না।

7. **Synchronous solving** — বড় মডেলে (হাজার+ element) request timeout
   হতে পারে। Asynchronous job queue এখনো নেই।

---

## লোকাল টেস্ট (deploy করার আগে যাচাই করতে চাইলে)

### C++ অংশ শুধু
```bash
# System dependencies (একবারের কাজ)
sudo apt-get install -y build-essential cmake libeigen3-dev
pip install pybind11

cd cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
```
সব টেস্ট pass করা উচিত (`test_stiffness`, `test_solver_integration`)।

### পূর্ণাঙ্গ FastAPI সার্ভিস (C++ + Python একসাথে)
```bash
# উপরের C++ build ধাপ আগে করুন, তারপর:
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

cd app
uvicorn main:app --reload --port 8080
```

```bash
curl http://localhost:8080/health

curl -X POST http://localhost:8080/jobs/analysis \
  -H "Content-Type: application/json" \
  -d '{
    "project_id": "test",
    "analysis_type": "linear-static",
    "model_payload": {
      "elements": [{"elementId": "e0", "category": "column", "label": "C1", "materialId": "m1", "sectionId": "s1", "startPoint": {"x":0,"y":0,"z":0}, "endPoint": {"x":0,"y":3,"z":0}, "connectionType": "moment"}],
      "materials": [{"materialId": "m1", "type": "concrete", "fc": 28, "poissonsRatio": 0.2}],
      "sections": [{"sectionId": "s1", "properties": {"area": 90000, "ixx": 675000000, "iyy": 300000000, "j": 500000000}}],
      "loadCases": [{"loadCaseId": "l1", "elementId": "e0", "applicationType": "point", "forceX": 5, "forceY": 0, "forceZ": 0, "positionRatio": 1.0}]
    }
  }'
```
ফেরত আসা `job_id` দিয়ে `GET /jobs/{job_id}` করলে displacement/force/warnings দেখা যাবে।

---

## GitHub-এ পুশ করা ও Deploy

```bash
git add -A
git commit -m "Phase 4a: C++ FE solver — Linear Static Analysis (Line Elements)"
git push
```

Railway-তে যদি ইতিমধ্যে কানেক্ট করা থাকে, push করলেই নতুন build শুরু
হবে — এবার build অনেক বেশি সময় নেবে আগের তুলনায় (C++ কম্পাইল হওয়ার
কারণে, সাধারণত ২-৫ মিনিট)। Railway build log-এ `cpp-builder` stage-এর
progress দেখা যাবে।

**যদি build fail করে:** build log-এর শেষ কয়েক লাইন পাঠান — `cmake`/
`g++` error, নাকি `civilos_solver import check` (Dockerfile-এর শেষে
থাকা sanity check) fail করেছে, তা থেকে সমস্যার জায়গা নির্দিষ্ট করা
যাবে।

---

## API এন্ডপয়েন্ট (Phase 4a)

| Method | Path | কাজ |
|---|---|---|
| GET | `/health` | লাইভনেস/রেডিনেস চেক |
| POST | `/jobs/analysis` | Linear Static Analysis সাবমিট ও সমাধান (synchronous) |
| GET | `/jobs/{job_id}` | ফলাফল (displacement, element forces, warnings) ফেরত দেয় |

`model_payload` এর প্রত্যাশিত shape `app/main.py`-এর `AnalysisJobRequest`
docstring এ, এবং সীমাবদ্ধতা `app/analysis_orchestration.py`-এর মডিউল
docstring এ বিস্তারিত।
