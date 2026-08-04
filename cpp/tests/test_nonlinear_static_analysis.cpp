#include "stiffness.h"
#include "solver.h"
#include "types.h"
#include <iostream>
#include <iomanip>
#include <cmath>

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

/**
 * একটা সাধারণ portal frame (দুই কলাম + এক beam) model বানায় — উভয়
 * কলামের base fixed, উপরে একটা beam দুই কলামকে জোড়া দেয়, top-left
 * node এ lateral (global Z) force।
 *
 * একক cantilever (base fixed, top free) এর বদলে portal frame ব্যবহার
 * করার কারণ: single cantilever statically-determinate — base hinge
 * release হলেই পুরো element একটা free-rotating mechanism হয়ে যায়
 * (কোনো redundancy নেই, তাই hinge yield করলে solve singular/unstable
 * হয়ে পড়ে)। Portal frame statically-indeterminate (একাধিক load path),
 * তাই একটা hinge yield করলেও বাকি structure (অন্য কলাম, beam) load
 * carry করা চালিয়ে যায় — শুধু structure softer হয়ে যায় (progressive
 * yielding এর বাস্তব আচরণ), সম্পূর্ণ mechanism হয় না। এটাই Nonlinear
 * Static Analysis এর বাস্তব ব্যবহারের ক্ষেত্র (একটা redundant বহুতল
 * ফ্রেমে progressive hinge formation), তাই test structure ও সেই
 * ধরনের হওয়া বাস্তবসম্মত।
 *
 * Node layout: n0=(0,0,0) base-left [fixed], n1=(0,L,0) top-left,
 * n2=(W,L,0) top-right, n3=(W,0,0) base-right [fixed]।
 * hingeAtStart=true হলে left column এর base (n0 প্রান্ত) এ hinge বসে।
 */
AnalysisModel buildPortalFrameModel(
    double L, double W, const SectionProperties& section, const MaterialProperties& material,
    bool hingeAtBaseOfLeftColumn, double lateralForceKN
) {
    AnalysisModel model;
    model.nodes = {
        Node3D{"n0", 0.0, 0.0, 0.0},
        Node3D{"n1", 0.0, L, 0.0},
        Node3D{"n2", W, L, 0.0},
        Node3D{"n3", W, 0.0, 0.0},
    };

    FrameElement leftCol;
    leftCol.elementId = "leftCol";
    leftCol.startNodeIndex = 0;
    leftCol.endNodeIndex = 1;
    leftCol.section = section;
    leftCol.material = material;
    leftCol.connectionType = "moment";
    leftCol.hingeAtStart = hingeAtBaseOfLeftColumn;
    leftCol.hingeAtEnd = false;

    FrameElement beam;
    beam.elementId = "beam";
    beam.startNodeIndex = 1;
    beam.endNodeIndex = 2;
    beam.section = section;
    beam.material = material;
    beam.connectionType = "moment";

    FrameElement rightCol;
    rightCol.elementId = "rightCol";
    rightCol.startNodeIndex = 3;
    rightCol.endNodeIndex = 2;
    rightCol.section = section;
    rightCol.material = material;
    rightCol.connectionType = "moment";

    model.elements = {leftCol, beam, rightCol};

    BoundaryCondition fixedLeft, fixedRight;
    fixedLeft.nodeIndex = 0;
    fixedLeft.restrainX = fixedLeft.restrainY = fixedLeft.restrainZ = true;
    fixedLeft.restrainRx = fixedLeft.restrainRy = fixedLeft.restrainRz = true;
    fixedRight.nodeIndex = 3;
    fixedRight.restrainX = fixedRight.restrainY = fixedRight.restrainZ = true;
    fixedRight.restrainRx = fixedRight.restrainRy = fixedRight.restrainRz = true;
    model.boundaryConditions = {fixedLeft, fixedRight};

    NodalLoad load;
    load.nodeIndex = 1; // top-left node
    load.fz = lateralForceKN; // global Z — কলাম local axis convention অনুযায়ী Mz bending তৈরি করে (উপরের ব্যাখ্যা)
    load.fx = load.fy = load.mx = load.my = load.mz = 0.0;
    model.loads = {load};

    return model;
}

int main() {
    std::cout << std::setprecision(8);

    const double E = 200e6;      // kN/m2
    const double rho = 7.85;     // tonne/m3
    const double A = 0.01;       // m2
    const double Ixx = 0.0001;   // m4
    const double Iyy = 0.0001;   // m4
    const double J = 0.00005;    // m4
    const double G = E / (2 * 1.3);

    SectionProperties section{A, Ixx, Iyy, J};
    MaterialProperties material{E, G, rho};

    const double L = 3.0;
    const double W = 4.0;
    const double forceKN = 10.0;

    // ============================================================
    // Test 1: No hinge assigned (or capacity=0, i.e. always elastic)
    // — Nonlinear Static result must match Linear Static exactly
    // ============================================================
    std::cout << "=== Test 1: No-Hinge Case Matches Linear Static Exactly ===\n";
    {
        auto model = buildPortalFrameModel(L, W, section, material, /*hingeAtBaseOfLeftColumn=*/false, forceKN);
        auto linear = solveLinearStatic(model);
        auto nonlinear = solveNonlinearStatic(model, /*numLoadSteps=*/5, 30, 1e-6);

        checkTrue("Linear Static solve succeeded", linear.success);
        checkTrue("Nonlinear Static solve succeeded (no hinge)", nonlinear.success);

        if (linear.success && nonlinear.success) {
            checkClose("Top-left-node uz: Nonlinear vs Linear Static (no hinge)",
                       nonlinear.nodalDisplacements[1](2), linear.nodalDisplacements[1](2), 1e-4);
            checkClose("Left-column base moment (Mz): Nonlinear vs Linear Static (no hinge)",
                       nonlinear.elementEndForces[0](5), linear.elementEndForces[0](5), 1e-4);
            checkTrue("Converged (no hinge, purely elastic)", nonlinear.converged);
            checkTrue("No hinges reported (hingeAtStart=false)", nonlinear.hingeStates.empty());
        }
    }

    // ============================================================
    // Test 2: Hinge assigned but capacity far exceeds demand — still
    // matches elastic Linear Static, hinge reported as not-yielded
    // ============================================================
    std::cout << "\n=== Test 2: Hinge Assigned, Capacity Never Exceeded ===\n";
    {
        SectionProperties strongSection = section;
        strongSection.yieldMomentMzKNm = 1000.0; // অনেক বড়, কখনো yield করবে না এই load এ

        auto model = buildPortalFrameModel(L, W, strongSection, material, /*hingeAtBaseOfLeftColumn=*/true, forceKN);
        auto linear = solveLinearStatic(model);
        auto nonlinear = solveNonlinearStatic(model, 5, 30, 1e-6);

        checkTrue("Linear Static solve succeeded", linear.success);
        checkTrue("Nonlinear Static solve succeeded (high-capacity hinge)", nonlinear.success);

        if (linear.success && nonlinear.success) {
            checkClose("Top-left-node uz: Nonlinear vs Linear Static (hinge never yields)",
                       nonlinear.nodalDisplacements[1](2), linear.nodalDisplacements[1](2), 1e-4);
            checkTrue("Converged", nonlinear.converged);
            checkTrue("Exactly one hinge state reported", nonlinear.hingeStates.size() == 1);
            if (!nonlinear.hingeStates.empty()) {
                checkTrue("Hinge NOT yielded (capacity far exceeds demand)", !nonlinear.hingeStates[0].yielded);
            }
        }
    }

    // ============================================================
    // Test 3: Hinge triggers — Linear Static (elastic) base moment
    // exceeds an intentionally-low yield capacity
    // ============================================================
    // Portal frame এ (statically indeterminate) hand-calculable simple
    // statics নেই cantilever এর মতো, তাই এখানে Linear Static সমাধানকেই
    // elastic-demand reference হিসেবে ব্যবহার করা হচ্ছে — yield capacity
    // ইচ্ছাকৃতভাবে সেই elastic demand এর অনেক নিচে সেট করে (নিশ্চিত
    // yield নিশ্চিত করতে)।
    std::cout << "\n=== Test 3: Hinge Triggers When Elastic Demand Exceeds Capacity ===\n";
    {
        auto elasticRefModel = buildPortalFrameModel(L, W, section, material, false, forceKN);
        auto elasticRef = solveLinearStatic(elasticRefModel);
        checkTrue("Elastic reference (Linear Static) solve succeeded", elasticRef.success);
        const double elasticBaseMoment = elasticRef.success ? std::abs(elasticRef.elementEndForces[0](5)) : 0.0;

        SectionProperties weakSection = section;
        weakSection.yieldMomentMzKNm = 0.5 * elasticBaseMoment; // নিশ্চিতভাবে yield করবে

        auto model = buildPortalFrameModel(L, W, weakSection, material, /*hingeAtBaseOfLeftColumn=*/true, forceKN);
        auto nonlinear = solveNonlinearStatic(model, /*numLoadSteps=*/20, 50, 1e-6);

        checkTrue("Nonlinear Static solve succeeded (weak hinge)", nonlinear.success);

        if (elasticRef.success && nonlinear.success) {
            checkTrue("Converged despite yielding", nonlinear.converged);
            checkTrue("Exactly one hinge state reported", nonlinear.hingeStates.size() == 1);
            if (!nonlinear.hingeStates.empty()) {
                checkTrue("Hinge DID yield (elastic demand > capacity)", nonlinear.hingeStates[0].yielded);
                checkClose("Final base moment ≈ yield capacity (elastic-perfectly-plastic)",
                           std::abs(nonlinear.hingeStates[0].finalMomentKNm), weakSection.yieldMomentMzKNm, 0.02);
            }

            // Yielded structure is softer → nonlinear top displacement > linear elastic prediction
            checkTrue("Nonlinear top-left-node uz > Linear Static uz (softening from yielding)",
                      std::abs(nonlinear.nodalDisplacements[1](2)) > std::abs(elasticRef.nodalDisplacements[1](2)));

            checkTrue("Displacement amplification ratio > 1.0 (nonlinear softening detected)",
                      nonlinear.maxDisplacementAmplificationRatio > 1.0);
        }
    }

    // ============================================================
    // Test 4: Result is reasonably insensitive to numLoadSteps once
    // yielding has fully developed (converged physics shouldn't
    // depend heavily on step count for this simple monotonic case)
    // ============================================================
    std::cout << "\n=== Test 4: Converged Result Is Stable Across Different Load-Step Counts ===\n";
    {
        auto elasticRefModel = buildPortalFrameModel(L, W, section, material, false, forceKN);
        auto elasticRef = solveLinearStatic(elasticRefModel);
        const double elasticBaseMoment = elasticRef.success ? std::abs(elasticRef.elementEndForces[0](5)) : 0.0;

        SectionProperties weakSection = section;
        weakSection.yieldMomentMzKNm = 0.5 * elasticBaseMoment;

        auto modelCoarse = buildPortalFrameModel(L, W, weakSection, material, true, forceKN);
        auto modelFine = buildPortalFrameModel(L, W, weakSection, material, true, forceKN);

        auto coarse = solveNonlinearStatic(modelCoarse, /*numLoadSteps=*/5, 50, 1e-6);
        auto fine = solveNonlinearStatic(modelFine, /*numLoadSteps=*/40, 50, 1e-6);

        checkTrue("Coarse-step solve succeeded", coarse.success);
        checkTrue("Fine-step solve succeeded", fine.success);

        if (coarse.success && fine.success) {
            checkClose("Top-left-node uz: coarse (5 steps) vs fine (40 steps) load stepping",
                       coarse.nodalDisplacements[1](2), fine.nodalDisplacements[1](2), 0.02);
            checkTrue("Both report the hinge as yielded",
                      coarse.hingeStates[0].yielded && fine.hingeStates[0].yielded);
        }
    }


    // ============================================================
    // Test 5: Invalid inputs are rejected gracefully
    // ============================================================
    std::cout << "\n=== Test 5: Invalid Inputs Rejected Without Crash ===\n";
    {
        auto model = buildPortalFrameModel(L, W, section, material, true, forceKN);

        auto noSteps = solveNonlinearStatic(model, /*numLoadSteps=*/0, 30, 1e-6);
        checkTrue("numLoadSteps=0 rejected", !noSteps.success);

        auto noIters = solveNonlinearStatic(model, 10, /*maxIterationsPerStep=*/0, 1e-6);
        checkTrue("maxIterationsPerStep=0 rejected", !noIters.success);

        auto badTol = solveNonlinearStatic(model, 10, 30, /*convergenceTolerance=*/-1.0);
        checkTrue("Negative convergenceTolerance rejected", !badTol.success);

        AnalysisModel noLoadModel = model;
        noLoadModel.loads.clear();
        auto noLoad = solveNonlinearStatic(noLoadModel, 10, 30, 1e-6);
        checkTrue("Model with no loads rejected", !noLoad.success);

        AnalysisModel noBcModel = model;
        noBcModel.boundaryConditions.clear();
        auto noBc = solveNonlinearStatic(noBcModel, 10, 30, 1e-6);
        checkTrue("Model with no boundary conditions rejected", !noBc.success);
    }

    // ============================================================
    // Test 6: Near-zero hinge capacity — hinge yields almost
    // immediately, but the (redundant) portal frame remains stable
    // ============================================================
    std::cout << "\n=== Test 6: Near-Zero Capacity Yields Immediately, Structure Remains Stable (Redundancy) ===\n";
    {
        SectionProperties tinySection = section;
        tinySection.yieldMomentMzKNm = 1e-6; // কার্যত শূন্য capacity — প্রথম load step এই yield করা উচিত

        auto model = buildPortalFrameModel(L, W, tinySection, material, true, forceKN);
        auto result = solveNonlinearStatic(model, 10, 30, 1e-6);

        // Portal frame redundant (একাধিক load path) — বাম কলামের base
        // hinge release হয়ে গেলেও beam ও ডান কলাম load carry করা চালিয়ে
        // যায়, তাই solve সফল হওয়া উচিত (mechanism না, single-cantilever
        // test এর বিপরীতে যেখানে এটা সত্যিকারের mechanism হয়ে যেত)।
        checkTrue("Solve succeeds (redundant structure remains stable after hinge yields)", result.success);
        if (result.success) {
            checkTrue("Hinge reported as yielded (near-zero capacity)",
                      !result.hingeStates.empty() && result.hingeStates[0].yielded);
            checkClose("Final base moment ≈ near-zero capacity",
                       std::abs(result.hingeStates[0].finalMomentKNm), tinySection.yieldMomentMzKNm, 0.5);
        }
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
