"""
CivilOS Structural — Compute Microservice
==========================================

এই সার্ভিস Cloud Run-এ চলবে এবং ভারী স্ট্রাকচারাল ক্যালকুলেশনের জন্য দায়ী।
Next.js ফ্রন্টএন্ড এই সার্ভিসে JSON payload পাঠায়, একটা job_id ফেরত পায়,
এবং সেই job_id দিয়ে status/result poll করে।

Phase 0-তে এই ফাইলে শুধু scaffolding আছে:
  - health check endpoint
  - job submission endpoint (in-memory job store, echo-back validation)
  - job status/result endpoint

Phase 4 (Analysis Engine)-এ আসল সলভার লজিক এখানে বসবে — তখন
C++ core (Eigen + OpenBLAS/MKL) কে Python থেকে pybind11 বা subprocess
দিয়ে কল করা হবে, অথবা ভারী অংশ পুরোপুরি C++ binary হিসেবে থাকবে এবং
এই FastAPI layer শুধু orchestration করবে।
"""

import uuid
from datetime import datetime, timezone
from enum import Enum
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

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
    allow_credentials=True,
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
        ..., description='e.g. "linear-static", "modal", "response-spectrum"'
    )
    model_payload: dict[str, Any] = Field(
        default_factory=dict,
        description="Structural model data (elements, materials, sections, loads). "
        "Phase 1-3 এ এই শেপ চূড়ান্ত হবে।",
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
    একটা বিশ্লেষণ (analysis) জব সাবমিট করে।

    Phase 0: শুধু validate করে, job_id দেয়, এবং সাথে সাথেই একটা
    placeholder "completed" স্ট্যাটাসে সেট করে দেয় (echo-back), যাতে
    Next.js ↔ Cloud Run যোগাযোগ end-to-end টেস্ট করা যায়।

    Phase 4: এখানে আসল async job dispatch হবে (Cloud Tasks queue এ পাঠানো,
    worker যেটা C++ solver কল করবে)।
    """
    job_id = str(uuid.uuid4())
    now = datetime.now(timezone.utc).isoformat()

    if not request.project_id:
        raise HTTPException(status_code=400, detail="project_id is required")

    # --- Phase 0 placeholder: এখানে কোনো প্রকৃত সলভার কল নেই ---
    placeholder_result = {
        "note": "Phase 0 scaffold — no real analysis performed yet.",
        "echoed_analysis_type": request.analysis_type,
        "element_count_received": len(request.model_payload.get("elements", [])),
    }

    job_record = JobStatusResponse(
        job_id=job_id,
        status=JobStatus.COMPLETED,
        submitted_at=now,
        completed_at=now,
        result=placeholder_result,
    )
    _job_store[job_id] = job_record

    return AnalysisJobResponse(job_id=job_id, status=JobStatus.COMPLETED, submitted_at=now)


@app.get("/jobs/{job_id}", response_model=JobStatusResponse)
def get_job_status(job_id: str) -> JobStatusResponse:
    """একটা জবের বর্তমান স্ট্যাটাস ও ফলাফল (থাকলে) ফেরত দেয়।"""
    job = _job_store.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return job
