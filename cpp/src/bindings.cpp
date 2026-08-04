/**
 * Python Bindings (pybind11)
 * ===========================
 * এই ফাইল C++ solver কে একটা Python module হিসেবে expose করে
 * (`civilos_solver` নামে), যাতে বিদ্যমান FastAPI layer (app/main.py)
 * থেকে সরাসরি import করে কল করা যায় — কোনো subprocess/IPC overhead
 * ছাড়া, একই process এর ভিতরে।
 *
 * এখন দুইটা analysis type expose করা হয়েছে: solve_linear_static এবং
 * solve_modal_analysis (Modal — natural frequency ও mode shape)।
 *
 * ডিজাইন সিদ্ধান্ত: C++ struct এর বদলে Python dict/list ব্যবহার করা
 * হচ্ছে ইনপুট/আউটপুট হিসেবে (pybind11 এর নিজস্ব class-binding এর বদলে)
 * — কারণ FastAPI ইতিমধ্যে JSON কে Python dict এ parse করে দেয়
 * (Pydantic models এর মাধ্যমে), তাই dict/list ভিত্তিক interface টা
 * FastAPI layer এর সাথে সবচেয়ে কম friction এ কাজ করে, এবং কোনো
 * অতিরিক্ত serialization layer লাগে না।
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "solver.h"

namespace py = pybind11;
using namespace civilos;

namespace {

AnalysisModel parseModelFromDict(const py::dict& input) {
    AnalysisModel model;

    for (auto item : input["nodes"].cast<py::list>()) {
        py::dict nodeDict = item.cast<py::dict>();
        Node3D node;
        node.nodeId = nodeDict["nodeId"].cast<std::string>();
        node.x = nodeDict["x"].cast<double>();
        node.y = nodeDict["y"].cast<double>();
        node.z = nodeDict["z"].cast<double>();
        model.nodes.push_back(node);
    }

    for (auto item : input["elements"].cast<py::list>()) {
        py::dict elemDict = item.cast<py::dict>();
        FrameElement element;
        element.elementId = elemDict["elementId"].cast<std::string>();
        element.startNodeIndex = elemDict["startNodeIndex"].cast<int>();
        element.endNodeIndex = elemDict["endNodeIndex"].cast<int>();
        element.connectionType = elemDict["connectionType"].cast<std::string>();

        py::dict sectionDict = elemDict["section"].cast<py::dict>();
        element.section.area = sectionDict["area"].cast<double>();
        element.section.ixx = sectionDict["ixx"].cast<double>();
        element.section.iyy = sectionDict["iyy"].cast<double>();
        element.section.j = sectionDict["j"].cast<double>();
        // yieldMomentMzKNm/yieldMomentMyKNm ঐচ্ছিক (Nonlinear Static
        // Analysis এর জন্য, types.h এর SectionProperties docstring
        // দেখুন) — না থাকলে 0.0 (অর্থাৎ "hinge capacity নির্দিষ্ট করা
        // নেই, সবসময় elastic")। পুরনো solver_input dict (এই field
        // ছাড়া) এখনো ভাঙবে না।
        element.section.yieldMomentMzKNm = sectionDict.contains("yieldMomentMzKNm")
            ? sectionDict["yieldMomentMzKNm"].cast<double>()
            : 0.0;
        element.section.yieldMomentMyKNm = sectionDict.contains("yieldMomentMyKNm")
            ? sectionDict["yieldMomentMyKNm"].cast<double>()
            : 0.0;

        // hingeAtStart/hingeAtEnd ঐচ্ছিক — না থাকলে false (কোনো hinge
        // নেই, elastic আচরণ)।
        element.hingeAtStart = elemDict.contains("hingeAtStart")
            ? elemDict["hingeAtStart"].cast<bool>()
            : false;
        element.hingeAtEnd = elemDict.contains("hingeAtEnd")
            ? elemDict["hingeAtEnd"].cast<bool>()
            : false;

        py::dict materialDict = elemDict["material"].cast<py::dict>();
        element.material.elasticModulus = materialDict["elasticModulus"].cast<double>();
        element.material.shearModulus = materialDict["shearModulus"].cast<double>();
        // density ঐচ্ছিক — Linear Static Analysis এ ব্যবহৃতই হয় না (mass
        // matrix লাগে না), তাই পুরনো solver_input dict (density-বিহীন)
        // এখনো ভাঙবে না। Modal Analysis এ density না থাকলে 0.0 ধরা হবে,
        // যা assembleGlobalMass কে একটা সব-শূন্য (singular) mass matrix
        // দেবে — solveModalAnalysis সেটা gracefully ব্যর্থ করবে (দেখুন
        // test_modal_analysis.cpp এর "Zero Density" টেস্ট), silently
        // ভুল ফলাফল না দিয়ে।
        element.material.density = materialDict.contains("density")
            ? materialDict["density"].cast<double>()
            : 0.0;

        model.elements.push_back(element);
    }

    // shellElements ঐচ্ছিক — Linear Static ছাড়া বাকি analysis type
    // (Modal/Buckling/P-Delta) এখনো shell সমর্থন করে না, এবং পুরনো
    // solver_input dict (shellElements key ছাড়া) এখনো ভাঙবে না।
    if (input.contains("shellElements")) {
        for (auto item : input["shellElements"].cast<py::list>()) {
            py::dict shellDict = item.cast<py::dict>();
            ShellElement shell;
            shell.elementId = shellDict["elementId"].cast<std::string>();

            py::list nodeIndicesList = shellDict["nodeIndices"].cast<py::list>();
            if (nodeIndicesList.size() != 4) {
                throw std::invalid_argument(
                    "Shell element '" + shell.elementId + "' এর nodeIndices অবশ্যই ঠিক ৪টা হতে "
                    "হবে (4-node quad), পাওয়া গেছে " + std::to_string(nodeIndicesList.size()) + "টা");
            }
            for (int i = 0; i < 4; ++i) {
                shell.nodeIndices[i] = nodeIndicesList[i].cast<int>();
            }

            shell.thickness = shellDict["thickness"].cast<double>();

            py::dict materialDict = shellDict["material"].cast<py::dict>();
            shell.material.elasticModulus = materialDict["elasticModulus"].cast<double>();
            shell.material.shearModulus = materialDict.contains("shearModulus")
                ? materialDict["shearModulus"].cast<double>()
                : 0.0; // shell এ সরাসরি ব্যবহৃত হয় না (poissonsRatio থেকে G derive করা হয় shell.cpp এ)
            shell.material.density = materialDict.contains("density")
                ? materialDict["density"].cast<double>()
                : 0.0;
            shell.material.poissonsRatio = materialDict["poissonsRatio"].cast<double>();

            model.shellElements.push_back(shell);
        }
    }

    for (auto item : input["boundaryConditions"].cast<py::list>()) {
        py::dict bcDict = item.cast<py::dict>();
        BoundaryCondition bc;
        bc.nodeIndex = bcDict["nodeIndex"].cast<int>();
        bc.restrainX = bcDict["restrainX"].cast<bool>();
        bc.restrainY = bcDict["restrainY"].cast<bool>();
        bc.restrainZ = bcDict["restrainZ"].cast<bool>();
        bc.restrainRx = bcDict["restrainRx"].cast<bool>();
        bc.restrainRy = bcDict["restrainRy"].cast<bool>();
        bc.restrainRz = bcDict["restrainRz"].cast<bool>();
        model.boundaryConditions.push_back(bc);
    }

    for (auto item : input["loads"].cast<py::list>()) {
        py::dict loadDict = item.cast<py::dict>();
        NodalLoad load;
        load.nodeIndex = loadDict["nodeIndex"].cast<int>();
        load.fx = loadDict["fx"].cast<double>();
        load.fy = loadDict["fy"].cast<double>();
        load.fz = loadDict["fz"].cast<double>();
        load.mx = loadDict["mx"].cast<double>();
        load.my = loadDict["my"].cast<double>();
        load.mz = loadDict["mz"].cast<double>();
        model.loads.push_back(load);
    }

    return model;
}

py::dict resultToDict(const AnalysisResult& result, const std::vector<BoundaryCondition>& boundaryConditions) {
    py::dict output;
    output["success"] = result.success;
    output["errorMessage"] = result.errorMessage;

    py::list displacements;
    for (const auto& d : result.nodalDisplacements) {
        py::dict dofs;
        dofs["ux"] = d(0);
        dofs["uy"] = d(1);
        dofs["uz"] = d(2);
        dofs["rx"] = d(3);
        dofs["ry"] = d(4);
        dofs["rz"] = d(5);
        displacements.append(dofs);
    }
    output["nodalDisplacements"] = displacements;

    py::list endForces;
    for (const auto& f : result.elementEndForces) {
        py::dict forces;
        // Local 12-DOF end force: [Fx1,Fy1,Fz1,Mx1,My1,Mz1, Fx2,Fy2,Fz2,Mx2,My2,Mz2]
        forces["startAxial"] = f(0);
        forces["startShearY"] = f(1);
        forces["startShearZ"] = f(2);
        forces["startTorsion"] = f(3);
        forces["startMomentY"] = f(4);
        forces["startMomentZ"] = f(5);
        forces["endAxial"] = f(6);
        forces["endShearY"] = f(7);
        forces["endShearZ"] = f(8);
        forces["endTorsion"] = f(9);
        forces["endMomentY"] = f(10);
        forces["endMomentZ"] = f(11);
        endForces.append(forces);
    }
    output["elementEndForces"] = endForces;

    // Phase 10n — Support reaction forces। nodeIndex ও DOF vector
    // দুটোই দেওয়া হচ্ছে যাতে frontend সরাসরি জানতে পারে কোন node এর
    // reaction (boundaryConditions[i].nodeIndex, result.reactionForces[i]
    // positionally একই ক্রমে, দেখুন types.h এর doc-comment)। global
    // coordinate এ, kN/kN·m একক।
    py::list reactions;
    for (size_t i = 0; i < result.reactionForces.size() && i < boundaryConditions.size(); ++i) {
        py::dict r;
        r["nodeIndex"] = boundaryConditions[i].nodeIndex;
        const auto& f = result.reactionForces[i];
        r["fx"] = f(0);
        r["fy"] = f(1);
        r["fz"] = f(2);
        r["mx"] = f(3);
        r["my"] = f(4);
        r["mz"] = f(5);
        reactions.append(r);
    }
    output["reactionForces"] = reactions;

    return output;
}

/**
 * ModalAnalysisResult → Python dict। Output shape:
 * {
 *   "success": bool,
 *   "errorMessage": str,
 *   "numModesComputed": int,
 *   "modes": [
 *     {
 *       "naturalFrequencyHz": float,
 *       "angularFrequencyRadPerSec": float,
 *       "modeShape": [{"ux":..,"uy":..,"uz":..,"rx":..,"ry":..,"rz":..}, ...]  // প্রতিটা node এর জন্য
 *     }, ...
 *   ]
 * }
 *
 * নোট: nodalDisplacements এর মতো flat list না রেখে "modes" এর ভেতরে
 * per-mode grouping রাখা হয়েছে — কারণ Modal Analysis এর ফলাফল
 * স্বভাবতই per-mode (প্রতিটা mode এর নিজস্ব frequency ও shape), আর
 * frontend কে (ভবিষ্যতে) mode-select করে একটার পর একটা shape animate/
 * visualize করতে হবে বলে এই grouping সরাসরি UI logic এর সাথে মিলবে।
 */
py::dict modalResultToDict(const ModalAnalysisResult& result) {
    py::dict output;
    output["success"] = result.success;
    output["errorMessage"] = result.errorMessage;
    output["numModesComputed"] = result.numModesComputed;

    py::list modes;
    for (int m = 0; m < result.numModesComputed; ++m) {
        py::dict modeDict;
        modeDict["naturalFrequencyHz"] = result.naturalFrequenciesHz[m];
        modeDict["angularFrequencyRadPerSec"] = result.angularFrequenciesRadPerSec[m];

        py::list shapeList;
        for (const auto& nodeShape : result.modeShapes[m]) {
            py::dict dofs;
            dofs["ux"] = nodeShape(0);
            dofs["uy"] = nodeShape(1);
            dofs["uz"] = nodeShape(2);
            dofs["rx"] = nodeShape(3);
            dofs["ry"] = nodeShape(4);
            dofs["rz"] = nodeShape(5);
            shapeList.append(dofs);
        }
        modeDict["modeShape"] = shapeList;

        modes.append(modeDict);
    }
    output["modes"] = modes;

    return output;
}

/**
 * BucklingAnalysisResult → Python dict। Output shape:
 * {
 *   "success": bool,
 *   "errorMessage": str,
 *   "numModesComputed": int,
 *   "modes": [
 *     {
 *       "criticalLoadFactor": float,
 *       "bucklingModeShape": [{"ux":..,"uy":..,"uz":..,"rx":..,"ry":..,"rz":..}, ...]
 *     }, ...
 *   ]
 * }
 *
 * modalResultToDict এর মতোই per-mode grouping — একই যুক্তি (frontend
 * mode-select করে shape animate/visualize করবে)।
 */
py::dict bucklingResultToDict(const BucklingAnalysisResult& result) {
    py::dict output;
    output["success"] = result.success;
    output["errorMessage"] = result.errorMessage;
    output["numModesComputed"] = result.numModesComputed;

    py::list modes;
    for (int m = 0; m < result.numModesComputed; ++m) {
        py::dict modeDict;
        modeDict["criticalLoadFactor"] = result.criticalLoadFactors[m];

        py::list shapeList;
        for (const auto& nodeShape : result.bucklingModeShapes[m]) {
            py::dict dofs;
            dofs["ux"] = nodeShape(0);
            dofs["uy"] = nodeShape(1);
            dofs["uz"] = nodeShape(2);
            dofs["rx"] = nodeShape(3);
            dofs["ry"] = nodeShape(4);
            dofs["rz"] = nodeShape(5);
            shapeList.append(dofs);
        }
        modeDict["bucklingModeShape"] = shapeList;

        modes.append(modeDict);
    }
    output["modes"] = modes;

    return output;
}

/**
 * PDeltaAnalysisResult → Python dict। Output shape:
 * {
 *   "success": bool,
 *   "errorMessage": str,
 *   "nodalDisplacements": [{"ux":..,"uy":..,...}, ...],   // resultToDict এর মতোই shape
 *   "elementEndForces": [{"startAxial":..,...}, ...],      // resultToDict এর মতোই shape
 *   "firstOrderAxialForces": [float, ...],                  // প্রতিটা element এর first-order axial force
 *   "maxDisplacementAmplificationRatio": float
 * }
 *
 * nodalDisplacements ও elementEndForces resultToDict() এর সাথে অভিন্ন
 * shape ব্যবহার করে (dict key নাম মিল রাখা হয়েছে ইচ্ছাকৃতভাবে) — যাতে
 * frontend একই parsing logic (runAnalysis.ts এর ParsedAnalysisResult)
 * P-Delta ফলাফলেও পুনর্ব্যবহার করতে পারে।
 */
py::dict pdeltaResultToDict(const PDeltaAnalysisResult& result) {
    py::dict output;
    output["success"] = result.success;
    output["errorMessage"] = result.errorMessage;

    py::list displacements;
    for (const auto& d : result.nodalDisplacements) {
        py::dict dofs;
        dofs["ux"] = d(0);
        dofs["uy"] = d(1);
        dofs["uz"] = d(2);
        dofs["rx"] = d(3);
        dofs["ry"] = d(4);
        dofs["rz"] = d(5);
        displacements.append(dofs);
    }
    output["nodalDisplacements"] = displacements;

    py::list endForces;
    for (const auto& f : result.elementEndForces) {
        py::dict forces;
        forces["startAxial"] = f(0);
        forces["startShearY"] = f(1);
        forces["startShearZ"] = f(2);
        forces["startTorsion"] = f(3);
        forces["startMomentY"] = f(4);
        forces["startMomentZ"] = f(5);
        forces["endAxial"] = f(6);
        forces["endShearY"] = f(7);
        forces["endShearZ"] = f(8);
        forces["endTorsion"] = f(9);
        forces["endMomentY"] = f(10);
        forces["endMomentZ"] = f(11);
        endForces.append(forces);
    }
    output["elementEndForces"] = endForces;

    py::list firstOrderAxial;
    for (double v : result.firstOrderAxialForces) {
        firstOrderAxial.append(v);
    }
    output["firstOrderAxialForces"] = firstOrderAxial;

    output["maxDisplacementAmplificationRatio"] = result.maxDisplacementAmplificationRatio;

    return output;
}

/**
 * ResponseSpectrumAnalysisResult → Python dict। Output shape:
 * {
 *   "success": bool,
 *   "errorMessage": str,
 *   "nodalDisplacements": [{"ux":..,"uy":..,...}, ...],   // resultToDict এর মতোই shape, কিন্তু magnitude (≥0)
 *   "elementEndForces": [{"startAxial":..,...}, ...],      // resultToDict এর মতোই shape, কিন্তু magnitude (≥0)
 *   "baseShear": float,
 *   "totalMassParticipationRatio": float,
 *   "numModesComputed": int,
 *   "modalDetails": [
 *     {"participationFactor": float, "effectiveMass": float, "spectralAccelerationG": float}, ...
 *   ]
 * }
 *
 * nodalDisplacements ও elementEndForces resultToDict() এর সাথে অভিন্ন key
 * shape ব্যবহার করে (frontend এর একই parsing logic পুনর্ব্যবহারের জন্য),
 * কিন্তু মান সবসময় ≥0 (CQC peak magnitude convention, types.h এর
 * ResponseSpectrumAnalysisResult docstring এ ব্যাখ্যা করা)। modalDetails
 * frontend কে প্রতিটা mode এর অবদান (participation) দেখাতে দেয় — যেমন
 * mass participation ratio কম হলে ব্যবহারকারীকে numModes বাড়ানোর
 * পরামর্শ দেওয়ার UI বানানো যাবে।
 */
py::dict responseSpectrumResultToDict(const ResponseSpectrumAnalysisResult& result) {
    py::dict output;
    output["success"] = result.success;
    output["errorMessage"] = result.errorMessage;

    py::list displacements;
    for (const auto& d : result.nodalDisplacements) {
        py::dict dofs;
        dofs["ux"] = d(0);
        dofs["uy"] = d(1);
        dofs["uz"] = d(2);
        dofs["rx"] = d(3);
        dofs["ry"] = d(4);
        dofs["rz"] = d(5);
        displacements.append(dofs);
    }
    output["nodalDisplacements"] = displacements;

    py::list endForces;
    for (const auto& f : result.elementEndForces) {
        py::dict forces;
        forces["startAxial"] = f(0);
        forces["startShearY"] = f(1);
        forces["startShearZ"] = f(2);
        forces["startTorsion"] = f(3);
        forces["startMomentY"] = f(4);
        forces["startMomentZ"] = f(5);
        forces["endAxial"] = f(6);
        forces["endShearY"] = f(7);
        forces["endShearZ"] = f(8);
        forces["endTorsion"] = f(9);
        forces["endMomentY"] = f(10);
        forces["endMomentZ"] = f(11);
        endForces.append(forces);
    }
    output["elementEndForces"] = endForces;

    output["baseShear"] = result.baseShear;
    output["totalMassParticipationRatio"] = result.totalMassParticipationRatio;
    output["numModesComputed"] = result.numModesComputed;

    py::list modalDetails;
    for (int m = 0; m < result.numModesComputed; ++m) {
        py::dict detail;
        detail["participationFactor"] = result.modalParticipationFactors[m];
        detail["effectiveMass"] = result.effectiveModalMasses[m];
        detail["spectralAccelerationG"] = result.modalSpectralAccelerations[m];
        modalDetails.append(detail);
    }
    output["modalDetails"] = modalDetails;

    return output;
}

/**
 * NonlinearStaticAnalysisResult → Python dict। Output shape:
 * {
 *   "success": bool,
 *   "errorMessage": str,
 *   "nodalDisplacements": [{"ux":..,...}, ...],  // resultToDict এর মতোই shape
 *   "elementEndForces": [{"startAxial":..,...}, ...],  // resultToDict এর মতোই shape
 *   "hingeStates": [
 *     {"elementIndex": int, "isAtStartNode": bool, "yielded": bool,
 *      "finalMomentKNm": float, "plasticRotationRad": float}, ...
 *   ],
 *   "totalLoadSteps": int,
 *   "totalNewtonIterations": int,
 *   "converged": bool,
 *   "maxDisplacementAmplificationRatio": float
 * }
 */
py::dict nonlinearStaticResultToDict(const NonlinearStaticAnalysisResult& result) {
    py::dict output;
    output["success"] = result.success;
    output["errorMessage"] = result.errorMessage;

    py::list displacements;
    for (const auto& d : result.nodalDisplacements) {
        py::dict dofs;
        dofs["ux"] = d(0);
        dofs["uy"] = d(1);
        dofs["uz"] = d(2);
        dofs["rx"] = d(3);
        dofs["ry"] = d(4);
        dofs["rz"] = d(5);
        displacements.append(dofs);
    }
    output["nodalDisplacements"] = displacements;

    py::list endForces;
    for (const auto& f : result.elementEndForces) {
        py::dict forces;
        forces["startAxial"] = f(0);
        forces["startShearY"] = f(1);
        forces["startShearZ"] = f(2);
        forces["startTorsion"] = f(3);
        forces["startMomentY"] = f(4);
        forces["startMomentZ"] = f(5);
        forces["endAxial"] = f(6);
        forces["endShearY"] = f(7);
        forces["endShearZ"] = f(8);
        forces["endTorsion"] = f(9);
        forces["endMomentY"] = f(10);
        forces["endMomentZ"] = f(11);
        endForces.append(forces);
    }
    output["elementEndForces"] = endForces;

    py::list hingeStates;
    for (const auto& h : result.hingeStates) {
        py::dict hingeDict;
        hingeDict["elementIndex"] = h.elementIndex;
        hingeDict["isAtStartNode"] = h.isAtStartNode;
        hingeDict["yielded"] = h.yielded;
        hingeDict["finalMomentKNm"] = h.finalMomentKNm;
        hingeDict["plasticRotationRad"] = h.plasticRotationRad;
        hingeStates.append(hingeDict);
    }
    output["hingeStates"] = hingeStates;

    output["totalLoadSteps"] = result.totalLoadSteps;
    output["totalNewtonIterations"] = result.totalNewtonIterations;
    output["converged"] = result.converged;
    output["maxDisplacementAmplificationRatio"] = result.maxDisplacementAmplificationRatio;

    return output;
}

/**
 * PushoverAnalysisResult → Python dict। Output shape:
 * {
 *   "success": bool,
 *   "errorMessage": str,
 *   "capacityCurve": [
 *     {"baseShearKN": float, "controlDisplacementM": float, "numHingesYielded": int}, ...
 *   ],
 *   "finalNodalDisplacements": [{"ux":..,...}, ...],
 *   "finalElementEndForces": [{"startAxial":..,...}, ...],
 *   "finalHingeStates": [{"elementIndex":.., "isAtStartNode":.., "yielded":.., "finalMomentKNm":.., "plasticRotationRad":..}, ...],
 *   "reachedTargetDisplacement": bool,
 *   "structureCollapsed": bool,
 *   "totalPushSteps": int,
 *   "totalNewtonIterations": int
 * }
 */
py::dict pushoverResultToDict(const PushoverAnalysisResult& result) {
    py::dict output;
    output["success"] = result.success;
    output["errorMessage"] = result.errorMessage;

    py::list capacityCurve;
    for (const auto& point : result.capacityCurve) {
        py::dict pointDict;
        pointDict["baseShearKN"] = point.baseShearKN;
        pointDict["controlDisplacementM"] = point.controlDisplacementM;
        pointDict["numHingesYielded"] = point.numHingesYielded;
        capacityCurve.append(pointDict);
    }
    output["capacityCurve"] = capacityCurve;

    py::list displacements;
    for (const auto& d : result.finalNodalDisplacements) {
        py::dict dofs;
        dofs["ux"] = d(0);
        dofs["uy"] = d(1);
        dofs["uz"] = d(2);
        dofs["rx"] = d(3);
        dofs["ry"] = d(4);
        dofs["rz"] = d(5);
        displacements.append(dofs);
    }
    output["finalNodalDisplacements"] = displacements;

    py::list endForces;
    for (const auto& f : result.finalElementEndForces) {
        py::dict forces;
        forces["startAxial"] = f(0);
        forces["startShearY"] = f(1);
        forces["startShearZ"] = f(2);
        forces["startTorsion"] = f(3);
        forces["startMomentY"] = f(4);
        forces["startMomentZ"] = f(5);
        forces["endAxial"] = f(6);
        forces["endShearY"] = f(7);
        forces["endShearZ"] = f(8);
        forces["endTorsion"] = f(9);
        forces["endMomentY"] = f(10);
        forces["endMomentZ"] = f(11);
        endForces.append(forces);
    }
    output["finalElementEndForces"] = endForces;

    py::list hingeStates;
    for (const auto& h : result.finalHingeStates) {
        py::dict hingeDict;
        hingeDict["elementIndex"] = h.elementIndex;
        hingeDict["isAtStartNode"] = h.isAtStartNode;
        hingeDict["yielded"] = h.yielded;
        hingeDict["finalMomentKNm"] = h.finalMomentKNm;
        hingeDict["plasticRotationRad"] = h.plasticRotationRad;
        hingeStates.append(hingeDict);
    }
    output["finalHingeStates"] = hingeStates;

    output["reachedTargetDisplacement"] = result.reachedTargetDisplacement;
    output["structureCollapsed"] = result.structureCollapsed;
    output["totalPushSteps"] = result.totalPushSteps;
    output["totalNewtonIterations"] = result.totalNewtonIterations;

    return output;
}

} // anonymous namespace

/**
 * Python থেকে কল করার জন্য প্রধান entry point।
 * ইনপুট dict shape (FastAPI Pydantic model থেকে .dict() করে আসবে):
 * {
 *   "nodes": [{"nodeId": str, "x": float, "y": float, "z": float}, ...],
 *   "elements": [{"elementId": str, "startNodeIndex": int, "endNodeIndex": int,
 *                 "connectionType": str, "section": {...}, "material": {...}}, ...],
 *   "boundaryConditions": [{"nodeIndex": int, "restrainX": bool, ...}, ...],
 *   "loads": [{"nodeIndex": int, "fx": float, ...}, ...]
 * }
 */
py::dict solveLinearStaticPy(py::dict input) {
    try {
        AnalysisModel model = parseModelFromDict(input);
        AnalysisResult result = solveLinearStatic(model);
        return resultToDict(result, model.boundaryConditions);
    } catch (const std::exception& e) {
        py::dict errorOutput;
        errorOutput["success"] = false;
        errorOutput["errorMessage"] = std::string("C++ solver exception: ") + e.what();
        errorOutput["nodalDisplacements"] = py::list();
        errorOutput["elementEndForces"] = py::list();
        errorOutput["reactionForces"] = py::list();
        return errorOutput;
    }
}

/**
 * Python থেকে কল করার জন্য Modal Analysis entry point। ইনপুট dict shape
 * solveLinearStaticPy এর মতোই (একই AnalysisModel), শুধু প্রতিটা element
 * এর material dict এ "density" key থাকা আবশ্যক (না থাকলে 0.0 ধরা হবে,
 * যা solver কে gracefully ব্যর্থ করবে — parseModelFromDict এর কমেন্ট
 * দেখুন)।
 *
 * numModes প্যারামিটার: কতগুলো mode ফেরত দিতে হবে (ডিফল্ট 12 —
 * solver.h এর solveModalAnalysis() docstring এ ব্যাখ্যা করা কারণে)।
 */
py::dict solveModalAnalysisPy(py::dict input, int numModes) {
    try {
        AnalysisModel model = parseModelFromDict(input);
        ModalAnalysisResult result = solveModalAnalysis(model, numModes);
        return modalResultToDict(result);
    } catch (const std::exception& e) {
        py::dict errorOutput;
        errorOutput["success"] = false;
        errorOutput["errorMessage"] = std::string("C++ solver exception: ") + e.what();
        errorOutput["numModesComputed"] = 0;
        errorOutput["modes"] = py::list();
        return errorOutput;
    }
}

/**
 * Python থেকে কল করার জন্য Linear Buckling Analysis entry point। ইনপুট
 * dict shape solveLinearStaticPy এর মতোই, এবং model.loads অবশ্যই
 * অখালি হতে হবে (solveLinearBuckling() এর docstring দেখুন — কোন load
 * pattern এর সাপেক্ষে buckling হচ্ছে তা জানা আবশ্যক)। density এখানে
 * ব্যবহৃত হয় না (buckling mass-independent), তাই material dict এ
 * density না থাকলেও সমস্যা নেই।
 *
 * numModes প্যারামিটার: কতগুলো mode ফেরত দিতে হবে (ডিফল্ট 6)।
 */
py::dict solveLinearBucklingPy(py::dict input, int numModes) {
    try {
        AnalysisModel model = parseModelFromDict(input);
        BucklingAnalysisResult result = solveLinearBuckling(model, numModes);
        return bucklingResultToDict(result);
    } catch (const std::exception& e) {
        py::dict errorOutput;
        errorOutput["success"] = false;
        errorOutput["errorMessage"] = std::string("C++ solver exception: ") + e.what();
        errorOutput["numModesComputed"] = 0;
        errorOutput["modes"] = py::list();
        return errorOutput;
    }
}

/**
 * Python থেকে কল করার জন্য P-Delta Analysis entry point। ইনপুট dict
 * shape solveLinearStaticPy এর মতোই, এবং model.loads অবশ্যই অখালি
 * হতে হবে (solvePDelta() এর docstring দেখুন)। density এখানে ব্যবহৃত
 * হয় না (P-Delta mass-independent), তাই material dict এ density না
 * থাকলেও সমস্যা নেই। কোনো numModes প্যারামিটার নেই — P-Delta একটা
 * static (non-eigenvalue) solve, single ফলাফল দেয়।
 */
py::dict solvePDeltaPy(py::dict input) {
    try {
        AnalysisModel model = parseModelFromDict(input);
        PDeltaAnalysisResult result = solvePDelta(model);
        return pdeltaResultToDict(result);
    } catch (const std::exception& e) {
        py::dict errorOutput;
        errorOutput["success"] = false;
        errorOutput["errorMessage"] = std::string("C++ solver exception: ") + e.what();
        errorOutput["nodalDisplacements"] = py::list();
        errorOutput["elementEndForces"] = py::list();
        errorOutput["firstOrderAxialForces"] = py::list();
        errorOutput["maxDisplacementAmplificationRatio"] = 1.0;
        return errorOutput;
    }
}

/**
 * Python থেকে কল করার জন্য Response Spectrum Analysis entry point। ইনপুট
 * dict shape solveLinearStaticPy এর মতোই, এবং প্রতিটা element এর material
 * dict এ "density" থাকা আবশ্যক (Modal Analysis এর মতোই কারণে — mass
 * matrix দরকার, density না থাকলে 0.0 ধরা হবে ও solver gracefully ব্যর্থ
 * হবে)।
 *
 * অতিরিক্ত প্যারামিটার:
 *   spectrum: [{"periodSec": float, "spectralAccelerationG": float}, ...]
 *             — ছোট-থেকে-বড় periodSec ক্রমে সাজানো তালিকা (BNBC 2020
 *             design spectrum, app/response_spectrum.py এ তৈরি হয়)
 *   direction_dof: 0 (X), 1 (Y), বা 2 (Z) — কোন দিকে ground motion
 *   damping_ratio: ডিফল্ট 0.05 (5%, concrete structure এর সাধারণ মান)
 *   num_modes: ডিফল্ট 12
 */
py::dict solveResponseSpectrumPy(py::dict input, py::list spectrum, int directionDof, double dampingRatio, int numModes) {
    try {
        AnalysisModel model = parseModelFromDict(input);

        std::vector<ResponseSpectrumPoint> spectrumPoints;
        for (auto item : spectrum) {
            py::dict pointDict = item.cast<py::dict>();
            ResponseSpectrumPoint point;
            point.periodSec = pointDict["periodSec"].cast<double>();
            point.spectralAccelerationG = pointDict["spectralAccelerationG"].cast<double>();
            spectrumPoints.push_back(point);
        }

        ResponseSpectrumAnalysisResult result = solveResponseSpectrum(
            model, spectrumPoints, directionDof, dampingRatio, numModes
        );
        return responseSpectrumResultToDict(result);
    } catch (const std::exception& e) {
        py::dict errorOutput;
        errorOutput["success"] = false;
        errorOutput["errorMessage"] = std::string("C++ solver exception: ") + e.what();
        errorOutput["nodalDisplacements"] = py::list();
        errorOutput["elementEndForces"] = py::list();
        errorOutput["baseShear"] = 0.0;
        errorOutput["totalMassParticipationRatio"] = 0.0;
        errorOutput["numModesComputed"] = 0;
        errorOutput["modalDetails"] = py::list();
        return errorOutput;
    }
}

/**
 * Python থেকে কল করার জন্য Nonlinear Static Analysis entry point। ইনপুট
 * dict shape solveLinearStaticPy এর মতোই, কিন্তু element এর section
 * dict এ ঐচ্ছিক "yieldMomentMzKNm"/"yieldMomentMyKNm" (default 0.0) ও
 * element dict এ ঐচ্ছিক "hingeAtStart"/"hingeAtEnd" (default false)
 * থাকতে পারে (parseModelFromDict দেখুন)।
 */
py::dict solveNonlinearStaticPy(py::dict input, int numLoadSteps, int maxIterationsPerStep, double convergenceTolerance) {
    try {
        AnalysisModel model = parseModelFromDict(input);
        NonlinearStaticAnalysisResult result = solveNonlinearStatic(
            model, numLoadSteps, maxIterationsPerStep, convergenceTolerance
        );
        return nonlinearStaticResultToDict(result);
    } catch (const std::exception& e) {
        py::dict errorOutput;
        errorOutput["success"] = false;
        errorOutput["errorMessage"] = std::string("C++ solver exception: ") + e.what();
        errorOutput["nodalDisplacements"] = py::list();
        errorOutput["elementEndForces"] = py::list();
        errorOutput["hingeStates"] = py::list();
        errorOutput["totalLoadSteps"] = 0;
        errorOutput["totalNewtonIterations"] = 0;
        errorOutput["converged"] = false;
        errorOutput["maxDisplacementAmplificationRatio"] = 1.0;
        return errorOutput;
    }
}

/**
 * Python থেকে কল করার জন্য Pushover Analysis entry point। ইনপুট dict
 * shape solveNonlinearStaticPy এর মতোই (element এ ঐচ্ছিক
 * hingeAtStart/hingeAtEnd, section এ ঐচ্ছিক yieldMomentMzKNm), plus
 * push-নির্দিষ্ট parameter।
 */
py::dict solvePushoverPy(
    py::dict input, int controlNodeIndex, int controlDof, double targetControlDisplacementM,
    double loadStepIncrement, int maxPushSteps, int maxIterationsPerStep, double convergenceTolerance
) {
    try {
        AnalysisModel model = parseModelFromDict(input);
        PushoverAnalysisResult result = solvePushover(
            model, controlNodeIndex, controlDof, targetControlDisplacementM,
            loadStepIncrement, maxPushSteps, maxIterationsPerStep, convergenceTolerance
        );
        return pushoverResultToDict(result);
    } catch (const std::exception& e) {
        py::dict errorOutput;
        errorOutput["success"] = false;
        errorOutput["errorMessage"] = std::string("C++ solver exception: ") + e.what();
        errorOutput["capacityCurve"] = py::list();
        errorOutput["finalNodalDisplacements"] = py::list();
        errorOutput["finalElementEndForces"] = py::list();
        errorOutput["finalHingeStates"] = py::list();
        errorOutput["reachedTargetDisplacement"] = false;
        errorOutput["structureCollapsed"] = false;
        errorOutput["totalPushSteps"] = 0;
        errorOutput["totalNewtonIterations"] = 0;
        return errorOutput;
    }
}

PYBIND11_MODULE(civilos_solver, m) {
    m.doc() = "CivilOS Structural — C++ FE Solver (Phase 4a: Linear Static + Modal + Linear Buckling + "
               "P-Delta + Response Spectrum + Nonlinear Static + Pushover Analysis)";
    m.def("solve_linear_static", &solveLinearStaticPy,
          "Solve a linear static structural analysis problem given nodes, elements, "
          "boundary conditions, and loads. Returns nodal displacements and element end forces.");
    m.def("solve_modal_analysis", &solveModalAnalysisPy,
          py::arg("input"), py::arg("num_modes") = 12,
          "Solve a modal (eigenvalue) analysis problem given nodes, elements, boundary "
          "conditions, and per-element material density. Returns natural frequencies (Hz) "
          "and mass-normalized mode shapes for the requested number of lowest modes.");
    m.def("solve_linear_buckling", &solveLinearBucklingPy,
          py::arg("input"), py::arg("num_modes") = 6,
          "Solve a linear (eigenvalue) buckling analysis problem given nodes, elements, "
          "boundary conditions, and a load pattern. Returns critical load factors and "
          "buckling mode shapes for the requested number of lowest-|factor| modes.");
    m.def("solve_pdelta", &solvePDeltaPy,
          py::arg("input"),
          "Solve a P-Delta (second-order geometric nonlinear static) analysis problem "
          "given nodes, elements, boundary conditions, and a load pattern. Uses a "
          "single-iteration approach (K+Kg)U=F. Returns P-Delta-modified nodal "
          "displacements, element end forces, first-order axial forces, and a "
          "displacement amplification ratio.");
    m.def("solve_response_spectrum", &solveResponseSpectrumPy,
          py::arg("input"), py::arg("spectrum"), py::arg("direction_dof"),
          py::arg("damping_ratio") = 0.05, py::arg("num_modes") = 12,
          "Solve a Response Spectrum Analysis (RSA) problem given nodes, elements, "
          "boundary conditions, per-element material density, and a tabulated design "
          "spectrum (list of {periodSec, spectralAccelerationG} points). Combines modal "
          "peak responses via CQC (Complete Quadratic Combination). direction_dof: 0=X, "
          "1=Y, 2=Z. Returns peak (magnitude, always >=0) nodal displacements, element "
          "end forces, base shear, mass participation ratio, and per-mode details.");
    m.def("solve_nonlinear_static", &solveNonlinearStaticPy,
          py::arg("input"), py::arg("num_load_steps") = 10, py::arg("max_iterations_per_step") = 30,
          py::arg("convergence_tolerance") = 1e-4,
          "Solve a Nonlinear Static Analysis problem (Concentrated Plastic Hinge method) "
          "given nodes, elements (optionally with hingeAtStart/hingeAtEnd and "
          "section.yieldMomentMzKNm), boundary conditions, and a load pattern. Uses "
          "incremental Load-Control Newton-Raphson iteration. Returns final converged "
          "nodal displacements, element end forces, per-hinge yield state, and "
          "convergence diagnostics.");
    m.def("solve_pushover", &solvePushoverPy,
          py::arg("input"), py::arg("control_node_index"), py::arg("control_dof"),
          py::arg("target_control_displacement_m"), py::arg("load_step_increment") = 0.02,
          py::arg("max_push_steps") = 200, py::arg("max_iterations_per_step") = 30,
          py::arg("convergence_tolerance") = 1e-4,
          "Solve a Pushover Analysis problem — pushes a fixed lateral load pattern "
          "(model loads) up to a target displacement at control_node_index/control_dof "
          "(0=ux, 1=uy, 2=uz), or until the structure collapses (singular tangent "
          "stiffness), whichever comes first. Uses the same Concentrated Plastic Hinge "
          "method as Nonlinear Static. Returns the base-shear-vs-control-displacement "
          "capacity curve, final converged state, and collapse/target-reached flags.");
}
