# CivilOS Structural — Compute Microservice

Phase 0 scaffold। FastAPI দিয়ে বানানো, Cloud Run-এ deploy হবে। এখনো কোনো
আসল সলভার লজিক নেই — শুধু job submission/polling এর কাঠামো, যাতে
Next.js ↔ Cloud Run যোগাযোগ প্রথমেই টেস্ট করা যায়। Phase 4-এ এখানে
C++ (Eigen/OpenBLAS/MKL) সলভার যোগ হবে।

> **যাচাই নোট:** `app/main.py`-এর FastAPI কোড লোকালি uvicorn দিয়ে চালিয়ে
> `/health`, `/jobs/analysis`, `/jobs/{job_id}` — তিনটা endpoint-ই টেস্ট
> করা হয়েছে এবং ঠিকভাবে কাজ করছে। কিন্তু `Dockerfile`-টা এই sandbox-এ
> build করে দেখা যায়নি (Docker এখানে নেই)। ধাপ ৪-এর
> `gcloud run deploy --source .` কমান্ডটাই আসলে প্রথম real build —
> সেটা Cloud Build ব্যবহার করে, তাই লোকাল Docker ছাড়াই চলবে। যদি সেই
> কমান্ড fail করে, আমাকে error log পাঠালে আমি ঠিক করে দেব।

---

## ১. লোকাল টেস্ট (ঐচ্ছিক, কিন্তু deploy করার আগে করে দেখা ভালো)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8080
```

তারপর আরেকটা টার্মিনালে:

```bash
curl http://localhost:8080/health
```

`{"status":"ok","service":"civilos-structural-solver"}` — এটা দেখলে ঠিক আছে।

---

## ২. GitHub-এ রিপো বানানো

```bash
cd civilos-structural-solver
git init
git add .
git commit -m "Phase 0: solver microservice scaffold"
git branch -M main
git remote add origin https://github.com/<আপনার-ইউজারনেম>/civilos-structural-solver.git
git push -u origin main
```

---

## ৩. GCP প্রজেক্ট প্রস্তুত করা (এক-বারের কাজ)

এই ধাপগুলা আপনার নিজের টার্মিনালে করতে হবে, যেখানে `gcloud` CLI ইনস্টল করা
আছে ([ইনস্টল লিংক](https://cloud.google.com/sdk/docs/install))। আমি (Claude)
এই sandbox থেকে আপনার GCP অ্যাকাউন্টে লগইন করতে পারি না — এই কমান্ডগুলা
নিজে চালাতে হবে।

```bash
# লগইন
gcloud auth login

# প্রজেক্ট সিলেক্ট করুন (না থাকলে তৈরি করুন)
gcloud config set project <আপনার-GCP-PROJECT-ID>

# দরকারি API চালু করুন
gcloud services enable run.googleapis.com artifactregistry.googleapis.com

# Artifact Registry-তে একটা Docker repo বানান (একবারই লাগবে)
gcloud artifacts repositories create civilos \
  --repository-format=docker \
  --location=asia-southeast1 \
  --description="CivilOS microservices"
```

> **অঞ্চল (region) নোট:** `asia-southeast1` (সিঙ্গাপুর) বাংলাদেশ থেকে
> সবচেয়ে কাছের Cloud Run region। চাইলে অন্য region ব্যবহার করতে পারেন,
> কিন্তু নিচের সব জায়গায় একই region বসাতে হবে।

---

## ৪. প্রথমবার ম্যানুয়াল deploy (GitHub Actions সেটআপ করার আগে টেস্ট করতে)

```bash
gcloud run deploy civilos-structural-solver \
  --source . \
  --region asia-southeast1 \
  --platform managed \
  --allow-unauthenticated \
  --memory 512Mi \
  --cpu 1
```

এটা শেষে একটা URL দেবে, যেমন:
`https://civilos-structural-solver-xxxxx-as.a.run.app`

এই URL টা কপি করে রাখুন — এটাই Vercel-এর env var-এ বসবে (নিচে দেখুন)।

টেস্ট করুন:
```bash
curl https://civilos-structural-solver-xxxxx-as.a.run.app/health
```

---

## ৫. স্বয়ংক্রিয় deploy (GitHub Actions) সেটআপ

প্রতিবার `main` branch-এ push করলেই যাতে অটোমেটিক Cloud Run-এ deploy হয়,
তার জন্য একটা Service Account বানাতে হবে:

```bash
gcloud iam service-accounts create civilos-deployer \
  --display-name="CivilOS CI/CD Deployer"

# প্রয়োজনীয় permission দিন
gcloud projects add-iam-policy-binding <আপনার-GCP-PROJECT-ID> \
  --member="serviceAccount:civilos-deployer@<আপনার-GCP-PROJECT-ID>.iam.gserviceaccount.com" \
  --role="roles/run.admin"

gcloud projects add-iam-policy-binding <আপনার-GCP-PROJECT-ID> \
  --member="serviceAccount:civilos-deployer@<আপনার-GCP-PROJECT-ID>.iam.gserviceaccount.com" \
  --role="roles/artifactregistry.writer"

gcloud projects add-iam-policy-binding <আপনার-GCP-PROJECT-ID> \
  --member="serviceAccount:civilos-deployer@<আপনার-GCP-PROJECT-ID>.iam.gserviceaccount.com" \
  --role="roles/iam.serviceAccountUser"

# একটা key ফাইল বানান
gcloud iam service-accounts keys create key.json \
  --iam-account=civilos-deployer@<আপনার-GCP-PROJECT-ID>.iam.gserviceaccount.com
```

তারপর GitHub repo-তে যান → **Settings → Secrets and variables → Actions →
New repository secret**, এবং এই তিনটা secret যোগ করুন:

| Secret নাম | মান |
|---|---|
| `GCP_PROJECT_ID` | আপনার GCP প্রজেক্ট আইডি |
| `GCP_SA_KEY` | `key.json` ফাইলের **পুরো কন্টেন্ট** (raw JSON, base64 না) |
| `GCP_REGION` | `asia-southeast1` (বা আপনি যেটা বেছেছেন) |

⚠️ **`key.json` ফাইলটা কখনো git-এ কমিট করবেন না।** Secret যোগ করার পর
লোকাল থেকে ফাইলটা ডিলিট করে দিন: `rm key.json`

এরপর থেকে `main`-এ যেকোনো push স্বয়ংক্রিয়ভাবে নতুন revision deploy করবে
(`.github/workflows/deploy.yml` দেখুন)।

---

## ৬. যাচাই করা

```bash
# Health check
curl <YOUR_CLOUD_RUN_URL>/health

# একটা টেস্ট জব সাবমিট করা
curl -X POST <YOUR_CLOUD_RUN_URL>/jobs/analysis \
  -H "Content-Type: application/json" \
  -d '{"project_id": "test-123", "analysis_type": "linear-static", "model_payload": {"elements": []}}'

# ফেরত আসা job_id দিয়ে status চেক
curl <YOUR_CLOUD_RUN_URL>/jobs/<JOB_ID>
```

---

## API এন্ডপয়েন্ট (Phase 0)

| Method | Path | কাজ |
|---|---|---|
| GET | `/health` | লাইভনেস/রেডিনেস চেক |
| POST | `/jobs/analysis` | একটা analysis job সাবমিট করে, job_id ফেরত দেয় |
| GET | `/jobs/{job_id}` | job-এর status ও ফলাফল ফেরত দেয় |

Phase 0-এ `/jobs/analysis` সাথে সাথেই একটা placeholder ফলাফল দিয়ে
`completed` স্ট্যাটাসে চলে যায় — আসল সলভার নেই। এটা শুধু end-to-end
যোগাযোগ যাচাই করার জন্য।
