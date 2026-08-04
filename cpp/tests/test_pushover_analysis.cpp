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
 * Portal frame model (test_nonlinear_static_analysis.cpp এর মতোই —
 * statically indeterminate, redundant, hinge release করলেও mechanism
 * হয় না)। lateralForceKN হলো base pattern magnitude — solvePushover()
 * নিজেই এটাকে বিভিন্ন load factor দিয়ে scale করবে ভিতরে ভিতরে, তাই
 * এখানে একটা মাঝারি reference magnitude দিলেই যথেষ্ট (push শুরু হবে
 * load factor≈0 থেকে, বাড়তে থাকবে)।
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
    load.fz = lateralForceKN; // global Z — কলাম local axis convention অনুযায়ী Mz bending তৈরি করে
    load.fx = load.fy = load.mx = load.my = load.mz = 0.0;
    model.loads = {load};

    return model;
}

int main() {
    std::cout << std::setprecision(8);

    const double E = 200e6;
    const double rho = 7.85;
    const double A = 0.01;
    const double Ixx = 0.0001;
    const double Iyy = 0.0001;
    const double J = 0.00005;
    const double G = E / (2 * 1.3);

    SectionProperties section{A, Ixx, Iyy, J};
    MaterialProperties material{E, G, rho};

    const double L = 3.0;
    const double W = 4.0;
    const double referenceForceKN = 10.0;

    // ============================================================
    // Test 1: No hinge — capacity curve stays perfectly linear
    // (elastic stiffness constant) up to the target displacement
    // ============================================================
    std::cout << "=== Test 1: Elastic (No-Hinge) Pushover — Linear Capacity Curve ===\n";
    {
        auto model = buildPortalFrameModel(L, W, section, material, /*hingeAtBaseOfLeftColumn=*/false, referenceForceKN);

        // ছোট target displacement — নিশ্চিতভাবে elastic range এর মধ্যে
        auto result = solvePushover(model, /*controlNodeIndex=*/1, /*controlDOF=*/2, /*targetDisp=*/0.002,
                                     /*loadStepIncrement=*/0.05, 200, 30, 1e-6);

        checkTrue("Pushover solve succeeded", result.success);
        checkTrue("Target displacement reached (elastic case)", result.reachedTargetDisplacement);
        checkTrue("No collapse (elastic case)", !result.structureCollapsed);
        checkTrue("Capacity curve has multiple points", result.capacityCurve.size() > 2);

        if (result.success && result.capacityCurve.size() > 2) {
            // Elastic হলে base shear / control displacement অনুপাত (secant
            // stiffness) সব point এ প্রায় constant থাকা উচিত (linear
            // elastic behavior, কোনো hinge yield হয়নি)।
            const double k0 = result.capacityCurve[1].baseShearKN / result.capacityCurve[1].controlDisplacementM;
            const double kLast = result.capacityCurve.back().baseShearKN / result.capacityCurve.back().controlDisplacementM;
            checkClose("Secant stiffness constant across elastic pushover (point 1 vs last)", kLast, k0, 0.01);
            checkTrue("Zero hinges yielded throughout (no hinge assigned)",
                      result.finalHingeStates.empty());
        }
    }

    // ============================================================
    // Test 2: Weak hinge — capacity curve softens (kink) after yield
    // ============================================================
    std::cout << "\n=== Test 2: Weak-Hinge Pushover — Capacity Curve Softens After Yield ===\n";
    {
        // প্রথমে elastic reference base moment বের করা (Linear Static দিয়ে)
        auto elasticRefModel = buildPortalFrameModel(L, W, section, material, false, referenceForceKN);
        auto elasticRef = solveLinearStatic(elasticRefModel);
        checkTrue("Elastic reference solve succeeded", elasticRef.success);
        const double elasticBaseMoment = elasticRef.success ? std::abs(elasticRef.elementEndForces[0](5)) : 0.0;

        SectionProperties weakSection = section;
        weakSection.yieldMomentMzKNm = 0.3 * elasticBaseMoment; // অবশ্যই push এর মধ্যেই yield করবে

        auto model = buildPortalFrameModel(L, W, weakSection, material, /*hingeAtBaseOfLeftColumn=*/true, referenceForceKN);

        // target displacement যথেষ্ট বড় যাতে hinge yield করার পরেও push চলতে থাকে
        auto result = solvePushover(model, 1, 2, /*targetDisp=*/0.05, 0.02, 300, 30, 1e-6);

        checkTrue("Pushover solve succeeded (weak hinge)", result.success);
        checkTrue("Capacity curve has enough points to detect softening", result.capacityCurve.size() > 10);

        if (result.success && result.capacityCurve.size() > 10) {
            // প্রাথমিক (elastic) secant stiffness বনাম শেষের দিকের secant
            // stiffness — hinge yield করার পর structure softer হয়ে যাওয়া
            // উচিত, তাই stiffness কমা উচিত।
            const auto& earlyPoint = result.capacityCurve[2];
            const auto& latePoint = result.capacityCurve.back();
            const double earlyStiffness = earlyPoint.baseShearKN / earlyPoint.controlDisplacementM;
            const double lateStiffness = latePoint.baseShearKN / latePoint.controlDisplacementM;

            checkTrue("Secant stiffness decreases after yielding (structure softens)",
                      lateStiffness < earlyStiffness);

            checkTrue("At least one hinge yielded by the end of the push",
                      latePoint.numHingesYielded > 0);

            checkTrue("Exactly one hinge assigned, and it yielded",
                      result.finalHingeStates.size() == 1 && result.finalHingeStates[0].yielded);

            checkClose("Final hinge moment clamps to yield capacity",
                       std::abs(result.finalHingeStates[0].finalMomentKNm), weakSection.yieldMomentMzKNm, 0.02);
        }
    }

    // ============================================================
    // Test 3: numHingesYielded count only increases (monotonic) as
    // the push progresses — a basic sanity/consistency check
    // ============================================================
    std::cout << "\n=== Test 3: Yielded-Hinge Count Is Monotonically Non-Decreasing Along the Push ===\n";
    {
        auto elasticRefModel = buildPortalFrameModel(L, W, section, material, false, referenceForceKN);
        auto elasticRef = solveLinearStatic(elasticRefModel);
        const double elasticBaseMoment = elasticRef.success ? std::abs(elasticRef.elementEndForces[0](5)) : 0.0;

        SectionProperties weakSection = section;
        weakSection.yieldMomentMzKNm = 0.3 * elasticBaseMoment;

        auto model = buildPortalFrameModel(L, W, weakSection, material, true, referenceForceKN);
        auto result = solvePushover(model, 1, 2, 0.05, 0.02, 300, 30, 1e-6);

        checkTrue("Pushover solve succeeded", result.success);
        if (result.success) {
            bool monotonic = true;
            for (size_t i = 1; i < result.capacityCurve.size(); ++i) {
                if (result.capacityCurve[i].numHingesYielded < result.capacityCurve[i - 1].numHingesYielded) {
                    monotonic = false;
                    break;
                }
            }
            checkTrue("numHingesYielded never decreases along the capacity curve", monotonic);
        }
    }

    // ============================================================
    // Test 4: Near-zero-capacity hinge on a redundant frame — the
    // push still proceeds (structure doesn't collapse immediately,
    // thanks to redundancy) and the capacity curve reflects an
    // early, low-stiffness plateau
    // ============================================================
    std::cout << "\n=== Test 4: Near-Zero Capacity — Frame Redundancy Keeps Push Going ===\n";
    {
        SectionProperties tinySection = section;
        tinySection.yieldMomentMzKNm = 1e-6;

        auto model = buildPortalFrameModel(L, W, tinySection, material, true, referenceForceKN);
        auto result = solvePushover(model, 1, 2, 0.02, 0.02, 300, 30, 1e-6);

        checkTrue("Pushover solve succeeded (near-zero capacity)", result.success);
        if (result.success) {
            checkTrue("Capacity curve has more than just the origin point", result.capacityCurve.size() > 1);
            checkTrue("At least one hinge yielded almost immediately",
                      !result.finalHingeStates.empty() && result.finalHingeStates[0].yielded);
        }
    }

    // ============================================================
    // Test 5: Invalid inputs are rejected gracefully
    // ============================================================
    std::cout << "\n=== Test 5: Invalid Inputs Rejected Without Crash ===\n";
    {
        auto model = buildPortalFrameModel(L, W, section, material, true, referenceForceKN);

        auto badControlNode = solvePushover(model, /*controlNodeIndex=*/99, 2, 0.01, 0.02, 100, 30, 1e-6);
        checkTrue("Out-of-range controlNodeIndex rejected", !badControlNode.success);

        auto badControlDOF = solvePushover(model, 1, /*controlDOF=*/7, 0.01, 0.02, 100, 30, 1e-6);
        checkTrue("Invalid controlDOF rejected", !badControlDOF.success);

        auto badTarget = solvePushover(model, 1, 2, /*targetDisp=*/-0.01, 0.02, 100, 30, 1e-6);
        checkTrue("Negative targetControlDisplacementM rejected", !badTarget.success);

        auto badIncrement = solvePushover(model, 1, 2, 0.01, /*loadStepIncrement=*/1.5, 100, 30, 1e-6);
        checkTrue("loadStepIncrement > 1.0 rejected", !badIncrement.success);

        auto badMaxSteps = solvePushover(model, 1, 2, 0.01, 0.02, /*maxPushSteps=*/0, 30, 1e-6);
        checkTrue("maxPushSteps=0 rejected", !badMaxSteps.success);

        AnalysisModel noLoadModel = model;
        noLoadModel.loads.clear();
        auto noLoad = solvePushover(noLoadModel, 1, 2, 0.01, 0.02, 100, 30, 1e-6);
        checkTrue("Model with no loads rejected", !noLoad.success);
    }

    // ============================================================
    // Test 6: maxPushSteps acts as a safety limit — a target that
    // can never be reached within the step budget still returns a
    // valid (partial) result, not an infinite loop or crash
    // ============================================================
    std::cout << "\n=== Test 6: maxPushSteps Safety Limit Stops the Push Gracefully ===\n";
    {
        auto model = buildPortalFrameModel(L, W, section, material, false, referenceForceKN);
        // অত্যন্ত বড় target displacement + অল্প কিছু step বাজেট — target
        // এ পৌঁছানো অসম্ভব এই বাজেটে।
        auto result = solvePushover(model, 1, 2, /*targetDisp=*/100.0, 0.05, /*maxPushSteps=*/5, 30, 1e-6);

        checkTrue("Solve completes without hanging/crashing", result.success);
        checkTrue("Target NOT reached within the tiny step budget", !result.reachedTargetDisplacement);
        checkTrue("totalPushSteps respects the maxPushSteps budget", result.totalPushSteps <= 5);
    }

    std::cout << "\n========================================\n";
    std::cout << "Results: " << testsPassed << " passed, " << testsFailed << " failed\n";
    std::cout << "========================================\n";

    return testsFailed > 0 ? 1 : 0;
}
