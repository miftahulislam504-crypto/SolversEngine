"""
CivilOS Structural — Compute Microservice
==========================================

এই সার্ভিস Cloud Run/Railway-তে চলে এবং ভারী স্ট্রাকচারাল ক্যালকুলেশনের
জন্য দায়ী। Next.js ফ্রন্টএন্ড এই সার্ভিসে JSON payload পাঠায়, একটা
job_id ফেরত পায়, এবং সেই job_id দিয়ে status/result poll করে।

Phase 4a (এই সংস্করণ): C++ FE solver (cpp/ ডিরেক্টরি, pybind11 দিয়ে
`civilos_solver` module হিসেবে expose করা) এখন সরাসরি কল হয় — Linear
Static Analysis এবং Modal Analysis, শুধু Line Element (Beam/Column/
Brace/Pile) এর জন্য। Slab/Wall (Area Element) এখনো solve হয় না (FE
mesh প্রয়োজন, Phase 4 এর পরের ধাপ)। বিস্তারিত সীমাবদ্ধতা
app/analysis_orchestration.py এর docstring এ।

সলভার synchronous ভাবে কল হয় (job submit করার সাথেই সমাধান শেষ হয়ে
যায়, Phase 0 এর মতো placeholder echo না, কিন্তু asynchronous queue-ও
না) — ছোট/মাঝারি মডেলের (কয়েকশো element পর্যন্ত) জন্য এটা যথেষ্ট
দ্রুত (sparse Cholesky সমাধান সাধারণত মিলিসেকেন্ডে শেষ হয়)। অনেক বড়
মডেলে (হাজার+ element) request timeout এড়াতে asynchronous job queue
(Cloud Tasks/Pub-Sub) প্রয়োজন হবে, যা এখনো implement করা হয়নি।
"""

import os
import sys
import time
import uuid
from datetime import datetime, timezone
from enum import Enum
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from app.analysis_orchestration import ModelParsingError, build_solver_model
from app.response_spectrum import build_bnbc_2020_spectrum

# civilos_solver (কম্পাইল করা C++ pybind11 module) খুঁজে পাওয়ার জন্য
# দুটো সম্ভাব্য জায়গা:
#   1. production (Dockerfile দিয়ে বিল্ড): .so ফাইল সরাসরি Python এর
#      site-packages এ কপি করা থাকে (Dockerfile দেখুন), তাই সাধারণ
#      `import civilos_solver` ই যথেষ্ট, কোনো path manipulation লাগে না।
#   2. local development/sandbox testing: module শুধু cpp/build/ এ
#      কম্পাইল হয়ে থাকে (এখনো কোথাও install করা হয়নি) — এই fallback
#      path সেই ক্ষেত্রে কাজ করে, যাতে site-packages এ কপি না করেও
#      লোকালি টেস্ট করা যায়।
# এই দুটো path explicit ভাবে লেখা হয়েছে (implicit sys.path নির্ভরতার
# বদলে) যাতে import ব্যর্থ হলে ঠিক কোথায় খুঁজেছে তা ডিবাগ করা সহজ হয়।
_LOCAL_BUILD_PATH = os.path.join(os.path.dirname(__file__), "..", "cpp", "build")
if os.path.isdir(_LOCAL_BUILD_PATH) and _LOCAL_BUILD_PATH not in sys.path:
    sys.path.insert(0, _LOCAL_BUILD_PATH)

try:
    import civilos_solver
except ImportError as import_error:
    # সলভার module কম্পাইল না হলে (dev environment এ cpp/build/ না
    # থাকলে, বা Docker build এ কোনো সমস্যা হলে) পুরো FastAPI app crash
    # করার বদলে একটা স্পষ্ট, অ্যাকশনযোগ্য error state এ থাকা ভালো —
    # /health endpoint তখনও কাজ করবে (ডিবাগ করা সহজ হয়), কিন্তু
    # /jobs/analysis একটা informative 503 দেবে।
    civilos_solver = None
    _SOLVER_IMPORT_ERROR = str(import_error)
else:
    _SOLVER_IMPORT_ERROR = None

app = FastAPI(
    title="CivilOS Structural Compute Service",
    version="0.1.0",
    description="Phase 0 scaffold — job submission/polling only, no real solver yet.",
)

# CORS: শুধু Vercel-এ deploy হওয়া Next.js app থেকে কল অ্যালাউ করা হবে।
# Phase 0-তে dev-এর সুবিধার জন্য চওড়া রাখা হলো; production-এ এটা
# নির্দিষ্ট origin-এ কমিয়ে আনা জরুরি (নিচে ENV var দিয়ে করা যায়)।
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # TODO(Phase 0 hardening): নির্দিষ্ট Vercel domain বসাতে হবে
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


class JobStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class AnalysisJobRequest(BaseModel):
    project_id: str = Field(..., description="Hub project ID")
    analysis_type: str = Field(
        ..., description='Phase 4a তে "linear-static", "modal", "buckling", "pdelta", '
        '"response-spectrum", "nonlinear-static" ও "pushover" সমর্থিত। অন্য মান (time-history '
        'ইত্যাদি) accept হবে কিন্তু 501 Not Implemented দেবে।'
    )
    model_payload: dict[str, Any] = Field(
        default_factory=dict,
        description="Expected keys: 'elements' (StructuralElement[]), 'materials' "
        "(StructuralMaterial[]), 'sections' (StructuralSection[], প্রতিটাতে precomputed "
        "'properties' field সহ — frontend এর computeSectionProperties() এর আউটপুট), "
        "'loadCases' (LoadCase[])। সব frontend এর src/lib/types/ এর টাইপ শেপ অনুসরণ করে। "
        "analysis_type='buckling'/'pdelta'/'nonlinear-static' এর জন্য loadCases অখালি থাকা "
        "আবশ্যক (কোন load pattern এর সাপেক্ষে গণনা হচ্ছে তা নির্দিষ্ট করতে হয়)।",
    )
    num_modes: int = Field(
        default=12,
        description="analysis_type='modal'/'buckling'/'response-spectrum' হলে কতগুলো mode "
        "ফেরত দিতে হবে ('modal'/'response-spectrum' এর ডিফল্ট 12, 'buckling' এ ব্যবহারিকভাবে "
        "6 যথেষ্ট)। 'linear-static'/'pdelta'/'nonlinear-static' এর জন্য এই field ignore করা হয়।",
    )
    seismic_zone: str | None = Field(
        default=None,
        description="analysis_type='response-spectrum' এর জন্য আবশ্যক — BNBC 2020 সিসমিক জোন "
        "('1'..'4')। response_spectrum.py এর build_bnbc_2020_spectrum() এ পাস হয়।",
    )
    site_class: str | None = Field(
        default=None,
        description="analysis_type='response-spectrum' এর জন্য আবশ্যক — BNBC 2020 site class "
        "('SA'..'SE')।",
    )
    direction_dof: int = Field(
        default=0,
        description="analysis_type='response-spectrum' এর জন্য — ground motion direction "
        "(0=X, 1=Y, 2=Z)। ডিফল্ট 0 (X-direction)।",
    )
    damping_ratio: float = Field(
        default=0.05,
        description="analysis_type='response-spectrum' এর জন্য — modal damping ratio, "
        "ডিফল্ট 0.05 (5%, concrete structure এর সাধারণ মান)।",
    )
    num_load_steps: int = Field(
        default=10,
        description="analysis_type='nonlinear-static' এর জন্য — সম্পূর্ণ load কে কতগুলো "
        "সমান increment এ ভাগ করা হবে (Newton-Raphson এর জন্য, solver.h এর "
        "solveNonlinearStatic() docstring দেখুন)। বেশি step মানে বেশি নিখুঁত hinge-yielding "
        "sequence, কিন্তু বেশি solve time।",
    )
    max_iterations_per_step: int = Field(
        default=30,
        description="analysis_type='nonlinear-static' এর জন্য — প্রতিটা load step এ সর্বোচ্চ "
        "Newton-Raphson iteration সংখ্যা, এর মধ্যে convergence না হলে non-convergence "
        "হিসেবে চিহ্নিত (crash না করে)।",
    )
    convergence_tolerance: float = Field(
        default=1e-4,
        description="analysis_type='nonlinear-static'/'pushover' এর জন্য — residual force norm / "
        "load norm এর অনুপাত এই মানের নিচে এলে convergence ধরা হয়।",
    )
    control_point_x: float | None = Field(
        default=None,
        description="analysis_type='pushover' এর জন্য আবশ্যক (control_point_y/z সহ) — capacity "
        "curve এর displacement অক্ষ কোন geometric point (x,y,z মিটার) থেকে পড়া হবে (সাধারণত "
        "roof-level node)। solver_input এর 'nodes' array এ coordinate ম্যাচ করে node index "
        "resolve করা হয় (frontend থেকে raw solver node index পাঠানো সম্ভব না, কারণ node merging/"
        "split ordering internal — build_solver_model() এর NodeGraph দেখুন)।",
    )
    control_point_y: float | None = Field(default=None, description="দেখুন control_point_x।")
    control_point_z: float | None = Field(default=None, description="দেখুন control_point_x।")
    control_dof: int = Field(
        default=2,
        description="analysis_type='pushover' এর জন্য — push direction এর translational DOF "
        "(0=X, 1=Y, 2=Z)। ডিফল্ট 2 (Z), কারণ frontend এর vertical-column local-axis convention "
        "অনুযায়ী lateral push সাধারণত Z-direction force দিয়ে হয় (solver.h এর buildCantileverModel "
        "টেস্ট কমেন্ট দেখুন, cpp repo)।",
    )
    target_control_displacement_m: float | None = Field(
        default=None,
        description="analysis_type='pushover' এর জন্য আবশ্যক — control node এর টার্গেট displacement "
        "(মিটার, magnitude), এখানে পৌঁছালে push বন্ধ হবে (অথবা structure collapse করলে আগেই থামবে)।",
    )
    load_step_increment: float = Field(
        default=0.02,
        description="analysis_type='pushover' এর জন্য — প্রতিটা push step এ load pattern কতটুকু "
        "বাড়বে (0-1 এর মধ্যে fraction, ডিফল্ট 0.02 = প্রতি step এ ২%)।",
    )
    max_push_steps: int = Field(
        default=200,
        description="analysis_type='pushover' এর জন্য — সর্বোচ্চ push step সংখ্যা (safety limit, "
        "target displacement এ কখনো না পৌঁছালে অসীম loop এড়াতে)।",
    )
    support_overrides: list[dict[str, Any]] | None = Field(
        default=None,
        description="Phase 5 — hardcoded 'base = fixed' (Y≈0) heuristic override করার ঐচ্ছিক "
        "input। প্রতিটা entry: {'x': float, 'y': float, 'z': float, 'supportType': "
        "'fixed'|'pinned'|'free'|'custom', 'restrainX'/'restrainY'/'restrainZ'/'restrainRx'/"
        "'restrainRy'/'restrainRz': bool (শুধু supportType='custom' হলে আবশ্যক)}। (x,y,z) কে "
        "solver model এর কোনো node coordinate এর সাথে ম্যাচ করে resolve করা হয় (control_point_x "
        "এর মতো একই 3-decimal rounding পদ্ধতিতে — build_solver_model() এর node_list এ, "
        "raw solver node index পাঠানো সম্ভব না কারণ node merging/split internal)। কোনো coordinate "
        "কোনো node এর সাথে না মিললে সেই entry নীরবে বাদ যায় না — একটা warning এ জানানো হয়, request "
        "ব্যর্থ হয় না (partial override এর ক্ষেত্রে বাকি override গুলো তবুও প্রয়োগ হয়)। None বা খালি "
        "list দিলে পুরনো Y≈0-based heuristic অপরিবর্তিতভাবে চলবে (backward compatible)। একই "
        "coordinate একাধিকবার override হলে শেষেরটা জিতবে।",
    )


class AnalysisJobResponse(BaseModel):
    job_id: str
    status: JobStatus
    submitted_at: str


class JobStatusResponse(BaseModel):
    job_id: str
    status: JobStatus
    submitted_at: str
    completed_at: str | None = None
    result: dict[str, Any] | None = None
    error: str | None = None


# Phase 0 এর জন্য in-memory store যথেষ্ট। Phase 4-এ এটা Firestore বা
# Cloud Tasks/Pub-Sub ভিত্তিক প্রকৃত job queue-তে সরানো হবে, কারণ
# Cloud Run instance রিস্টার্ট হলে in-memory ডেটা হারিয়ে যায়।
_job_store: dict[str, JobStatusResponse] = {}


@app.get("/health")
def health_check() -> dict[str, str]:
    """Cloud Run readiness/liveness probe এর জন্য।"""
    return {"status": "ok", "service": "civilos-structural-solver"}


@app.post("/jobs/analysis", response_model=AnalysisJobResponse)
def submit_analysis_job(request: AnalysisJobRequest) -> AnalysisJobResponse:
    """
    একটা বিশ্লেষণ (analysis) জব সাবমিট করে এবং synchronously সমাধান করে
    (docstring এ ফাইলের শুরুতে ব্যাখ্যা করা — বড় মডেলে এটা পরিবর্তন
    হবে asynchronous queue এ)।

    Phase 4a: "linear-static", "modal", "buckling", "pdelta",
    "response-spectrum", "nonlinear-static" ও "pushover" প্রকৃতভাবে
    সমর্থিত। বাকি analysis type (Master Plan এর Section 6 — Time
    History ইত্যাদি) 501 Not Implemented দেয়, কোনো fake ফলাফল না দিয়ে।
    """
    job_id = str(uuid.uuid4())
    now = datetime.now(timezone.utc).isoformat()

    if not request.project_id:
        raise HTTPException(status_code=400, detail="project_id is required")

    if request.analysis_type not in ("linear-static", "modal", "buckling", "pdelta", "response-spectrum", "nonlinear-static", "pushover"):
        raise HTTPException(
            status_code=501,
            detail=f"Analysis type '{request.analysis_type}' এখনো সমর্থিত না। Phase 4a তে "
                   f"'linear-static', 'modal', 'buckling', 'pdelta', 'response-spectrum', "
                   f"'nonlinear-static' ও 'pushover' পাওয়া যায়। বাকি সব analysis type "
                   f"(time-history ইত্যাদি) পরবর্তী ধাপে যোগ হবে।",
        )

    if civilos_solver is None:
        raise HTTPException(
            status_code=503,
            detail=f"C++ solver module লোড করা যায়নি (import error: {_SOLVER_IMPORT_ERROR})। "
                   f"এটা একটা deployment সমস্যা — cpp/build/ ঠিকভাবে হয়েছে কিনা Dockerfile "
                   f"যাচাই করুন।",
        )

    if request.analysis_type == "response-spectrum":
        if not request.seismic_zone or not request.site_class:
            raise HTTPException(
                status_code=422,
                detail="analysis_type='response-spectrum' এর জন্য 'seismic_zone' ('1'-'4') ও "
                       "'site_class' ('SA'-'SE') উভয়ই আবশ্যক।",
            )
        if request.direction_dof not in (0, 1, 2):
            raise HTTPException(
                status_code=422,
                detail="'direction_dof' অবশ্যই 0 (X), 1 (Y), বা 2 (Z) হতে হবে।",
            )

    if request.analysis_type == "pushover":
        if request.control_point_x is None or request.control_point_y is None or request.control_point_z is None:
            raise HTTPException(
                status_code=422,
                detail="analysis_type='pushover' এর জন্য 'control_point_x', 'control_point_y' ও "
                       "'control_point_z' — তিনটাই আবশ্যক — capacity curve এর displacement অক্ষ "
                       "কোন node থেকে পড়া হবে তা নির্দিষ্ট করে (সাধারণত push করা roof node এর "
                       "coordinate)।",
            )
        if request.target_control_displacement_m is None or request.target_control_displacement_m <= 0:
            raise HTTPException(
                status_code=422,
                detail="analysis_type='pushover' এর জন্য 'target_control_displacement_m' আবশ্যক ও "
                       "positive হতে হবে।",
            )
        if request.control_dof not in (0, 1, 2):
            raise HTTPException(
                status_code=422,
                detail="'control_dof' অবশ্যই 0 (X), 1 (Y), বা 2 (Z) হতে হবে।",
            )

    elements = request.model_payload.get("elements", [])
    materials = request.model_payload.get("materials", [])
    sections = request.model_payload.get("sections", [])
    load_cases = request.model_payload.get("loadCases", [])

    started_at = time.monotonic()

    try:
        solver_input, warnings, sub_element_registry = build_solver_model(
            elements=elements,
            materials=materials,
            sections=sections,
            load_cases=load_cases,
            supported_node_ids=set(),  # Phase 4a তে explicit support marking নেই, docstring দেখুন analysis_orchestration.py এ
            support_overrides=request.support_overrides,  # Phase 5 — ঐচ্ছিক, None হলে heuristic অপরিবর্তিত থাকে
        )
    except ModelParsingError as parsing_error:
        raise HTTPException(
            status_code=422,
            detail=f"মডেল পার্স করতে ব্যর্থ: {parsing_error}",
        )

    if request.analysis_type == "modal":
        solver_output = civilos_solver.solve_modal_analysis(solver_input, request.num_modes)
    elif request.analysis_type == "buckling":
        solver_output = civilos_solver.solve_linear_buckling(solver_input, request.num_modes)
    elif request.analysis_type == "pdelta":
        solver_output = civilos_solver.solve_pdelta(solver_input)
    elif request.analysis_type == "response-spectrum":
        try:
            spectrum = build_bnbc_2020_spectrum(request.seismic_zone, request.site_class)
        except ValueError as spectrum_error:
            raise HTTPException(status_code=422, detail=str(spectrum_error))
        solver_output = civilos_solver.solve_response_spectrum(
            solver_input, spectrum, request.direction_dof, request.damping_ratio, request.num_modes
        )
    elif request.analysis_type == "nonlinear-static":
        solver_output = civilos_solver.solve_nonlinear_static(
            solver_input, request.num_load_steps, request.max_iterations_per_step, request.convergence_tolerance
        )
    elif request.analysis_type == "pushover":
        # control_point (x,y,z) কে solver_input["nodes"] এ coordinate ম্যাচ
        # করে node index এ resolve করা — analysis_orchestration.py এর
        # NodeGraph.index_of() এর মতোই 3-decimal rounding দিয়ে (floating-
        # point coordinate তুলনার সময় সামঞ্জস্যপূর্ণ থাকতে)।
        control_node_idx = None
        target_key = (
            round(request.control_point_x, 3),
            round(request.control_point_y, 3),
            round(request.control_point_z, 3),
        )
        for i, node in enumerate(solver_input["nodes"]):
            node_key = (round(node["x"], 3), round(node["y"], 3), round(node["z"], 3))
            if node_key == target_key:
                control_node_idx = i
                break
        if control_node_idx is None:
            raise HTTPException(
                status_code=422,
                detail=f"control_point ({request.control_point_x}, {request.control_point_y}, "
                       f"{request.control_point_z}) মডেলের কোনো node এর সাথে মেলেনি। এই coordinate "
                       f"এ কোনো element endpoint আছে কিনা যাচাই করুন।",
            )
        solver_output = civilos_solver.solve_pushover(
            solver_input, control_node_idx, request.control_dof,
            request.target_control_displacement_m, request.load_step_increment,
            request.max_push_steps, request.max_iterations_per_step, request.convergence_tolerance
        )
    else:
        solver_output = civilos_solver.solve_linear_static(solver_input)

    elapsed_seconds = time.monotonic() - started_at

    if not solver_output["success"]:
        job_record = JobStatusResponse(
            job_id=job_id,
            status=JobStatus.FAILED,
            submitted_at=now,
            completed_at=datetime.now(timezone.utc).isoformat(),
            error=solver_output["errorMessage"],
        )
        _job_store[job_id] = job_record
        return AnalysisJobResponse(job_id=job_id, status=JobStatus.FAILED, submitted_at=now)

    if request.analysis_type == "modal":
        result_payload = _build_modal_result_payload(solver_input, solver_output, warnings, elapsed_seconds)
    elif request.analysis_type == "buckling":
        result_payload = _build_buckling_result_payload(solver_input, solver_output, warnings, elapsed_seconds)
    elif request.analysis_type == "pdelta":
        result_payload = _build_pdelta_result_payload(solver_input, solver_output, warnings, elapsed_seconds)
    elif request.analysis_type == "response-spectrum":
        result_payload = _build_response_spectrum_result_payload(solver_input, solver_output, warnings, elapsed_seconds)
    elif request.analysis_type == "nonlinear-static":
        result_payload = _build_nonlinear_static_result_payload(
            solver_input, solver_output, warnings, sub_element_registry, elapsed_seconds
        )
    elif request.analysis_type == "pushover":
        result_payload = _build_pushover_result_payload(
            solver_input, solver_output, warnings, sub_element_registry, elapsed_seconds
        )
    else:
        result_payload = _build_linear_static_result_payload(
            solver_input, solver_output, warnings, sub_element_registry, elapsed_seconds
        )

    job_record = JobStatusResponse(
        job_id=job_id,
        status=JobStatus.COMPLETED,
        submitted_at=now,
        completed_at=datetime.now(timezone.utc).isoformat(),
        result=result_payload,
    )
    _job_store[job_id] = job_record

    return AnalysisJobResponse(job_id=job_id, status=JobStatus.COMPLETED, submitted_at=now)


def _build_linear_static_result_payload(
    solver_input: dict[str, Any],
    solver_output: dict[str, Any],
    warnings: list[str],
    sub_element_registry: list[dict[str, Any]],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """একটা সফল Linear Static solver_output কে frontend-বোধগম্য result dict এ রূপান্তর করে।"""
    # C++ solver শুধু positional index জানে (elementEndForces[i] মানে
    # solver_input["elements"][i])। sub_element_registry দিয়ে প্রতিটা
    # entry কে তার আসল (frontend-এর জানা) elementId ও sub-element
    # position (subStartRatio/subEndRatio) সহ ফেরত দেওয়া হচ্ছে — মূলত
    # mid-span-load-split element গুলোর ক্ষেত্রে গুরুত্বপূর্ণ, যেখানে
    # একটা original elementId একাধিক solver sub-element এ ভাগ হয়ে
    # গেছে। Split না হলে প্রতিটা elementId এর জন্য একটাই entry থাকবে
    # (subStartRatio=0.0, subEndRatio=1.0) — আগের (positional-only)
    # আচরণের সাথে backward-compatible, শুধু elementId এখন explicit।
    element_end_forces_with_id = []
    for registry_entry in sub_element_registry:
        idx = registry_entry["solverElementIndex"]
        forces = dict(solver_output["elementEndForces"][idx])
        forces["elementId"] = registry_entry["originalElementId"]
        forces["subStartRatio"] = registry_entry["subStartRatio"]
        forces["subEndRatio"] = registry_entry["subEndRatio"]
        element_end_forces_with_id.append(forces)

    return {
        "analysisType": "linear-static",
        "nodalDisplacements": solver_output["nodalDisplacements"],
        "elementEndForces": element_end_forces_with_id,
        # Phase 10n — support reaction forces (nodeIndex + 6-DOF vector,
        # global coordinate)। C++ resultToDict() থেকে সরাসরি pass-through,
        # কোনো transformation দরকার নেই এখানে (positional split-element
        # concern reaction এ প্রযোজ্য না — reaction শুধু boundaryCondition
        # node এ, sub-element split এর সাথে সম্পর্কহীন)।
        "reactionForces": solver_output["reactionForces"],
        # nodes: node_list[i] এর coordinate সরাসরি nodalDisplacements[i] এর
        # সাথে ম্যাপ হয় (Phase 8a) — frontend আগে শুধু positional index
        # পেত, কোন displacement কোন story/grid point তা জানার কোনো উপায়
        # ছিল না। এখন প্রতিটা analysis result এ full node_list ফেরত যায়।
        "nodes": solver_input["nodes"],
        "nodeCount": len(solver_input["nodes"]),
        "elementCount": len(solver_input["elements"]) + len(solver_input.get("shellElements", [])),
        "solveTimeSeconds": round(elapsed_seconds, 4),
        "warnings": warnings,
    }


def _build_modal_result_payload(
    solver_input: dict[str, Any],
    solver_output: dict[str, Any],
    warnings: list[str],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """একটা সফল Modal Analysis solver_output কে frontend-বোধগম্য result dict এ রূপান্তর করে।"""
    return {
        "analysisType": "modal",
        "numModesComputed": solver_output["numModesComputed"],
        "modes": solver_output["modes"],  # প্রতিটা: {naturalFrequencyHz, angularFrequencyRadPerSec, modeShape}
        "nodes": solver_input["nodes"],  # modeShape[i] এর coordinate — Phase 8a
        "nodeCount": len(solver_input["nodes"]),
        "elementCount": len(solver_input["elements"]) + len(solver_input.get("shellElements", [])),
        "solveTimeSeconds": round(elapsed_seconds, 4),
        "warnings": warnings,
    }


def _build_buckling_result_payload(
    solver_input: dict[str, Any],
    solver_output: dict[str, Any],
    warnings: list[str],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """একটা সফল Linear Buckling Analysis solver_output কে frontend-বোধগম্য result dict এ রূপান্তর করে।"""
    buckling_warnings = list(warnings)
    # সাহায্যকারী warning: সবচেয়ে critical mode positive না হলে (অর্থাৎ
    # বর্তমান load direction এ buckling প্রাসঙ্গিক না হলে) ব্যবহারকারীকে
    # জানানো — types.h এর BucklingAnalysisResult docstring এ ব্যাখ্যা করা
    # sign convention অনুযায়ী।
    modes = solver_output.get("modes", [])
    if modes and modes[0]["criticalLoadFactor"] < 0:
        buckling_warnings.append(
            "⚠️ সবচেয়ে critical mode এর load factor ঋণাত্মক — অর্থাৎ প্রয়োগকৃত load pattern এর "
            "বিপরীত দিকে (load উল্টো করলে) buckling ঘটবে, বর্তমান দিকে না। বর্তমান load direction "
            "এ কার্যত এই মডেলে buckling risk কম হতে পারে (অথবা load pattern পুনর্বিবেচনা করুন)।"
        )

    return {
        "analysisType": "buckling",
        "numModesComputed": solver_output["numModesComputed"],
        "modes": modes,  # প্রতিটা: {criticalLoadFactor, bucklingModeShape}
        "nodes": solver_input["nodes"],  # bucklingModeShape[i] এর coordinate — Phase 8a
        "nodeCount": len(solver_input["nodes"]),
        "elementCount": len(solver_input["elements"]) + len(solver_input.get("shellElements", [])),
        "solveTimeSeconds": round(elapsed_seconds, 4),
        "warnings": buckling_warnings,
    }


def _build_pdelta_result_payload(
    solver_input: dict[str, Any],
    solver_output: dict[str, Any],
    warnings: list[str],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """একটা সফল P-Delta Analysis solver_output কে frontend-বোধগম্য result dict এ রূপান্তর করে।"""
    pdelta_warnings = list(warnings)
    # সাহায্যকারী warning: amplification ratio অস্বাভাবিক বড় হলে (structure
    # buckling এর কাছাকাছি হতে পারে) ব্যবহারকারীকে জানানো — solver.h এর
    # solvePDelta() docstring এ ব্যাখ্যা করা penalty-method সীমাবদ্ধতা
    # অনুযায়ী (matrix decomposition নিজে থেকে সবসময় এই অবস্থা ধরতে
    # পারে না, তাই এই heuristic threshold টা একটা proactive সতর্কতা)।
    amplification = solver_output.get("maxDisplacementAmplificationRatio", 1.0)
    if amplification > 3.0:
        pdelta_warnings.append(
            f"⚠️ Displacement amplification ratio ({amplification:.2f}x) অস্বাভাবিক বেশি — এটা "
            f"ইঙ্গিত দিতে পারে যে structure এই load এ critical buckling load এর কাছাকাছি চলে "
            f"গেছে। Linear Buckling Analysis (analysis_type='buckling') চালিয়ে critical load "
            f"factor যাচাই করার পরামর্শ দেওয়া হচ্ছে — শুধু এই ফলাফলের সফলতার উপর নির্ভর করবেন না।"
        )

    return {
        "analysisType": "pdelta",
        "nodalDisplacements": solver_output["nodalDisplacements"],
        "elementEndForces": solver_output["elementEndForces"],
        "firstOrderAxialForces": solver_output["firstOrderAxialForces"],
        "maxDisplacementAmplificationRatio": amplification,
        "nodes": solver_input["nodes"],  # Phase 8a
        "nodeCount": len(solver_input["nodes"]),
        "elementCount": len(solver_input["elements"]) + len(solver_input.get("shellElements", [])),
        "solveTimeSeconds": round(elapsed_seconds, 4),
        "warnings": pdelta_warnings,
    }


def _build_response_spectrum_result_payload(
    solver_input: dict[str, Any],
    solver_output: dict[str, Any],
    warnings: list[str],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """একটা সফল Response Spectrum Analysis solver_output কে frontend-বোধগম্য result dict এ রূপান্তর করে।"""
    rsa_warnings = list(warnings)
    # সাহায্যকারী warning: mass participation ratio কম হলে (BNBC/ASCE 7
    # সাধারণত ≥90% দাবি করে) ব্যবহারকারীকে numModes বাড়ানোর পরামর্শ —
    # solver.h এর solveResponseSpectrum() docstring এ ব্যাখ্যা করা কারণ
    # অনুযায়ী।
    participation = solver_output.get("totalMassParticipationRatio", 0.0)
    if participation < 0.90:
        rsa_warnings.append(
            f"⚠️ Mass participation ratio ({participation * 100:.1f}%) 90% এর কম — BNBC 2020 ও "
            f"অধিকাংশ design code অনুযায়ী এটা যথেষ্ট না। num_modes বাড়িয়ে (বর্তমান: "
            f"{solver_output.get('numModesComputed', 0)}) আবার চালানোর পরামর্শ দেওয়া হচ্ছে, যাতে "
            f"আরও বেশি mode এর অবদান captured হয়।"
        )

    # nodalDisplacements এখানে CQC (Complete Quadratic Combination) দিয়ে
    # combine করা — প্রতিটা DOF-এ ফলাফল sqrt(non-negative quadratic form),
    # তাই সবসময় ≥0 (কোনো sign/direction তথ্য নেই, solver.cpp এর
    # computeCQCCorrelationCoefficient() এর পরের ধাপ দেখুন)। Story Drift
    # (Phase 8c) দুইটা story-level displacement সরাসরি বিয়োগ করে বের করা
    # হয়, কিন্তু magnitude-only মান বিয়োগ করলে ভুল sign/মান আসতে পারে
    # (দুই story একই দিকে সমান দূরত্ব move করলেও CQC তাদের individual
    # sign হারিয়ে ফেলেছে)। displacementIsMagnitudeOnly=True flag দিয়ে
    # এটা explicit করা হলো, যাতে ভবিষ্যতে কোনো consumer (frontend বা
    # অন্য কোনো ক্লায়েন্ট) না জেনে ভুল subtraction না করে।
    return {
        "analysisType": "response-spectrum",
        "nodalDisplacements": solver_output["nodalDisplacements"],
        "displacementIsMagnitudeOnly": True,
        "elementEndForces": solver_output["elementEndForces"],
        "baseShear": solver_output["baseShear"],
        "totalMassParticipationRatio": participation,
        "numModesComputed": solver_output["numModesComputed"],
        "modalDetails": solver_output["modalDetails"],  # প্রতিটা: {participationFactor, effectiveMass, spectralAccelerationG}
        "nodes": solver_input["nodes"],  # Phase 8a
        "nodeCount": len(solver_input["nodes"]),
        "elementCount": len(solver_input["elements"]) + len(solver_input.get("shellElements", [])),
        "solveTimeSeconds": round(elapsed_seconds, 4),
        "warnings": rsa_warnings,
    }


def _build_nonlinear_static_result_payload(
    solver_input: dict[str, Any],
    solver_output: dict[str, Any],
    warnings: list[str],
    sub_element_registry: list[dict[str, Any]],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """একটা সফল Nonlinear Static Analysis solver_output কে frontend-বোধগম্য result dict এ রূপান্তর করে।"""
    nonlinear_warnings = list(warnings)

    if not solver_output.get("converged", False):
        nonlinear_warnings.append(
            f"⚠️ সব load step convergence tolerance এর মধ্যে পৌঁছায়নি "
            f"({solver_output.get('totalNewtonIterations', 0)}টা মোট Newton-Raphson iteration এর পরেও)। "
            f"এর মানে হতে পারে structure এই load level এ প্রায়-mechanism (অত্যধিক hinge yield করেছে) "
            f"অথবা num_load_steps/max_iterations_per_step বাড়ানো প্রয়োজন। ফলাফল সতর্কতার সাথে ব্যবহার করুন।"
        )

    # Phase 10r — elementEndForces ও hingeStates দুটোই backend থেকে শুধু
    # positional/solver-internal elementIndex দিয়ে আসে (frontend elementId
    # না) — sub_element_registry দিয়ে সেই gap পূরণ করা হচ্ছে, ঠিক
    # _build_linear_static_result_payload() এ elementEndForces এর জন্য
    # যেভাবে করা হয়েছিল সেই একই প্যাটার্নে (দেখুন সেই function এর
    # doc-comment বিস্তারিত ব্যাখ্যার জন্য)। এটা linear-static এর
    # আগে থেকেই gap ছিল, Phase 10r এর জন্য এখন ঠিক করা হলো কারণ hinge
    # state visualization এ elementId ছাড়া কোনো element highlight করা
    # যাবে না।
    element_end_forces_with_id = []
    for registry_entry in sub_element_registry:
        idx = registry_entry["solverElementIndex"]
        forces = dict(solver_output["elementEndForces"][idx])
        forces["elementId"] = registry_entry["originalElementId"]
        forces["subStartRatio"] = registry_entry["subStartRatio"]
        forces["subEndRatio"] = registry_entry["subEndRatio"]
        element_end_forces_with_id.append(forces)

    # hingeStates এ elementIndex সরাসরি sub_element_registry এর
    # solverElementIndex এর সাথে মেলে (একই solver-internal frame-element
    # array, hinge C++ লুপে যে e ব্যবহার হয় সেটাই সেই index)।
    registry_by_solver_index = {r["solverElementIndex"]: r for r in sub_element_registry}
    hinge_states_with_id = []
    for h in solver_output.get("hingeStates", []):
        registry_entry = registry_by_solver_index.get(h["elementIndex"])
        if registry_entry is None:
            continue  # defensive — সাধারণত hওয়ার কথা না, প্রতিটা solver element registry তে থাকা উচিত
        enriched = dict(h)
        enriched["elementId"] = registry_entry["originalElementId"]
        hinge_states_with_id.append(enriched)

    yielded_hinges = [h for h in hinge_states_with_id if h.get("yielded")]
    if yielded_hinges:
        nonlinear_warnings.append(
            f"ℹ️ {len(yielded_hinges)}টা hinge yield করেছে প্রয়োগকৃত load এর অধীনে "
            f"({len(hinge_states_with_id)}টার মধ্যে assign করা hinge)। "
            f"hingeStates এ প্রতিটার বিস্তারিত অবস্থা দেখুন।"
        )

    return {
        "analysisType": "nonlinear-static",
        "nodalDisplacements": solver_output["nodalDisplacements"],
        "displacementIsMagnitudeOnly": False,  # signed, সরাসরি drift subtraction-safe (RSA এর মতো CQC না)
        "elementEndForces": element_end_forces_with_id,
        "hingeStates": hinge_states_with_id,  # প্রতিটা: {elementId, elementIndex, isAtStartNode, yielded, finalMomentKNm, plasticRotationRad}
        "totalLoadSteps": solver_output["totalLoadSteps"],
        "totalNewtonIterations": solver_output["totalNewtonIterations"],
        "converged": solver_output["converged"],
        "maxDisplacementAmplificationRatio": solver_output["maxDisplacementAmplificationRatio"],
        "nodes": solver_input["nodes"],  # Phase 8a
        "nodeCount": len(solver_input["nodes"]),
        "elementCount": len(solver_input["elements"]) + len(solver_input.get("shellElements", [])),
        "solveTimeSeconds": round(elapsed_seconds, 4),
        "warnings": nonlinear_warnings,
    }


def _build_pushover_result_payload(
    solver_input: dict[str, Any],
    solver_output: dict[str, Any],
    warnings: list[str],
    sub_element_registry: list[dict[str, Any]],
    elapsed_seconds: float,
) -> dict[str, Any]:
    """একটা সফল Pushover Analysis solver_output কে frontend-বোধগম্য result dict এ রূপান্তর করে।"""
    pushover_warnings = list(warnings)

    if solver_output.get("structureCollapsed", False):
        pushover_warnings.append(
            "🔴 Push এর সময় structure এর tangent stiffness singular হয়ে গেছে (mechanism/collapse) — "
            "target displacement এ পৌঁছানোর আগেই push থেমেছে। capacityCurve এর শেষ বিন্দু হলো "
            "structure এর ultimate lateral capacity (এই hinge configuration এর অধীনে)।"
        )
    elif not solver_output.get("reachedTargetDisplacement", False):
        pushover_warnings.append(
            f"⚠️ target_control_displacement_m এ পৌঁছানো যায়নি max_push_steps "
            f"({solver_output.get('totalPushSteps', 0)}টা step ব্যবহৃত) এর মধ্যে। max_push_steps বা "
            f"load_step_increment বাড়ানোর পরামর্শ দেওয়া হচ্ছে।"
        )

    # Phase 10r — finalElementEndForces ও finalHingeStates দুটোই backend
    # থেকে শুধু positional/solver-internal elementIndex দিয়ে আসে, ঠিক
    # nonlinear-static এর মতোই একই gap — একই sub_element_registry
    # enrichment প্যাটার্ন এখানেও প্রয়োগ করা হলো।
    final_element_end_forces_with_id = []
    for registry_entry in sub_element_registry:
        idx = registry_entry["solverElementIndex"]
        forces = dict(solver_output["finalElementEndForces"][idx])
        forces["elementId"] = registry_entry["originalElementId"]
        forces["subStartRatio"] = registry_entry["subStartRatio"]
        forces["subEndRatio"] = registry_entry["subEndRatio"]
        final_element_end_forces_with_id.append(forces)

    registry_by_solver_index = {r["solverElementIndex"]: r for r in sub_element_registry}
    final_hinge_states_with_id = []
    for h in solver_output.get("finalHingeStates", []):
        registry_entry = registry_by_solver_index.get(h["elementIndex"])
        if registry_entry is None:
            continue
        enriched = dict(h)
        enriched["elementId"] = registry_entry["originalElementId"]
        final_hinge_states_with_id.append(enriched)

    yielded_hinges = [h for h in final_hinge_states_with_id if h.get("yielded")]
    if yielded_hinges:
        pushover_warnings.append(
            f"ℹ️ চূড়ান্ত push অবস্থায় {len(yielded_hinges)}টা hinge yield করেছে "
            f"({len(final_hinge_states_with_id)}টার মধ্যে assign করা hinge)।"
        )

    return {
        "analysisType": "pushover",
        "capacityCurve": solver_output["capacityCurve"],  # প্রতিটা: {baseShearKN, controlDisplacementM, numHingesYielded}
        "finalNodalDisplacements": solver_output["finalNodalDisplacements"],
        "finalElementEndForces": final_element_end_forces_with_id,
        "finalHingeStates": final_hinge_states_with_id,  # প্রতিটা: {elementId, elementIndex, isAtStartNode, yielded, finalMomentKNm, plasticRotationRad}
        "reachedTargetDisplacement": solver_output["reachedTargetDisplacement"],
        "structureCollapsed": solver_output["structureCollapsed"],
        "totalPushSteps": solver_output["totalPushSteps"],
        "totalNewtonIterations": solver_output["totalNewtonIterations"],
        "displacementIsMagnitudeOnly": False,  # finalNodalDisplacements signed (nonlinear static এর মতো)
        "nodes": solver_input["nodes"],  # Phase 8a
        "nodeCount": len(solver_input["nodes"]),
        "elementCount": len(solver_input["elements"]) + len(solver_input.get("shellElements", [])),
        "solveTimeSeconds": round(elapsed_seconds, 4),
        "warnings": pushover_warnings,
    }



@app.get("/jobs/{job_id}", response_model=JobStatusResponse)
def get_job_status(job_id: str) -> JobStatusResponse:
    """একটা জবের বর্তমান স্ট্যাটাস ও ফলাফল (থাকলে) ফেরত দেয়।"""
    job = _job_store.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return job
