# CivilOS Structural — Compute Microservice

Phase 0 scaffold। FastAPI দিয়ে বানানো, **Railway**-তে deploy হবে। এখনো কোনো
আসল সলভার লজিক নেই — শুধু job submission/polling এর কাঠামো, যাতে
Next.js ↔ Railway যোগাযোগ প্রথমেই টেস্ট করা যায়। Phase 4-এ এখানে
C++ (Eigen/OpenBLAS/MKL) সলভার যোগ হবে।

> **কেন Cloud Run না, Railway?** Google Cloud Run কারিগরিভাবে ঠিক
> প্ল্যাটফর্ম হলেও, GCP প্রজেক্ট চালু করতে বিলিং একাউন্ট (কার্ড) বাধ্যতামূলক —
> এমনকি ফ্রি টায়ারের জন্যও। Railway কার্ড ছাড়াই GitHub দিয়ে সাইন আপ করে
> ফ্রি টায়ারে (মাসে $৫ ক্রেডিট) deploy করা যায়, এবং পুরো প্রসেসটা ব্রাউজার
> dashboard থেকেই হয় — `gcloud` CLI বা কোনো টার্মিনাল লাগে না। কোডে কোনো
> Google-নির্দিষ্ট জিনিস ছিল না (plain FastAPI + Dockerfile), তাই এই সরানো
> কোনো লজিক পরিবর্তন করেনি।

> **যাচাই নোট:** `app/main.py`-এর FastAPI কোড লোকালি uvicorn দিয়ে চালিয়ে
> `/health`, `/jobs/analysis`, `/jobs/{job_id}` — তিনটা endpoint-ই টেস্ট
> করা হয়েছে এবং ঠিকভাবে কাজ করছে। `Dockerfile`-টা এই sandbox-এ build করে
> দেখা যায়নি (Docker এখানে নেই), কিন্তু Railway নিজেই push করার পর build
> করে — সেটাই প্রথম real build। যদি সেটা fail করে, Railway dashboard-এর
> "Deployments" ট্যাব থেকে build log কপি করে আমাকে পাঠালে ঠিক করে দেব।

---

## ১. লোকাল টেস্ট (ঐচ্ছিক, কিন্তু deploy করার আগে করে দেখা ভালো)

যদি লোকালি টেস্ট করার সুযোগ থাকে (যেমন অন্য কারো কম্পিউটারে):

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

আপনার যেহেতু লোকাল টার্মিনাল নেই, এই ধাপ স্কিপ করে সরাসরি ধাপ ২-এ যেতে পারেন —
Railway নিজেই build করে দেখাবে কাজ করছে কিনা।

---

## ২. GitHub-এ রিপো বানানো / আপডেট করা

এই ফোল্ডারের কনটেন্ট আপনার GitHub রিপোতে commit ও push করুন (আপনার
স্বাভাবিক workflow অনুযায়ী — zip থেকে GitHub-এ আপলোড করে commit)।

⚠️ যদি আগে GCP-এর জন্য বানানো রিপোতেই কাজ করেন, নিশ্চিত করুন পুরনো
`.github/workflow/` (singular, "s" ছাড়া) ফোল্ডারটা ডিলিট করে দিয়েছেন —
নতুন সঠিক ফোল্ডার হলো `.github/workflows/` (plural)। GitHub Actions শুধু
plural নামের ফোল্ডার থেকেই workflow পড়ে; singular নামটা থাকলে GitHub সেটা
চুপচাপ ignore করে, কোনো error ছাড়াই — তাই এই ভুলটা ধরা কঠিন। এই zip-এ
সঠিক ফোল্ডার নামেই workflow ফাইল দেওয়া আছে।

---

## ৩. Railway প্রজেক্ট বানানো (এক-বারের কাজ, কার্ড লাগবে না)

সব ধাপ ব্রাউজারে, [railway.app](https://railway.app)-এ:

1. **railway.app**-এ যান → **"Login"** → **"Login with GitHub"** দিয়ে সাইন
   ইন করুন (কোনো কার্ড চাইবে না)।
2. Dashboard-এ **"New Project"** → **"Deploy from GitHub repo"** সিলেক্ট
   করুন।
3. যে রিপোতে এই সলভার কোড push করেছেন সেটা সিলেক্ট করুন। যদি রিপো লিস্টে
   না দেখায়, **"Configure GitHub App"** থেকে Railway-কে সেই রিপো অ্যাক্সেস
   দিন।
4. Railway স্বয়ংক্রিয়ভাবে `Dockerfile` খুঁজে পাবে এবং build শুরু করবে।
   প্রথম deploy শেষ হতে ১-৩ মিনিট লাগতে পারে।
5. Deploy শেষ হলে, প্রজেক্টের সেই service-এ ক্লিক করুন → **"Settings"**
   ট্যাব → **"Networking"** সেকশন → **"Generate Domain"** বাটনে ক্লিক করুন।
   এটা একটা পাবলিক URL দেবে, যেমন:
   `https://civilos-structural-solver-production.up.railway.app`

এই URL টা কপি করে রাখুন — এটাই Vercel-এর env var-এ বসবে (নিচে দেখুন)।

টেস্ট করুন (ব্রাউজারে সরাসরি URL-এর শেষে `/health` জুড়ে দিয়ে খুললেই হবে,
অথবা curl থাকলে):
```bash
curl https://civilos-structural-solver-production.up.railway.app/health
```

> **ফ্রি টায়ার সম্পর্কে:** Railway-র ফ্রি "Trial" প্ল্যানে সাইন আপের সময়
> কিছু একবারের ক্রেডিট দেওয়া হয় (কার্ড ছাড়া)। এই Phase 0 সার্ভিসটা খুবই
> হালকা (512Mi মেমরি, প্রায় idle) বলে এই ক্রেডিটে অনেকদিন চলার কথা। ক্রেডিট
> শেষ হয়ে গেলে Railway dashboard-এ notification দেখাবে — তখন সিদ্ধান্ত
> নেওয়া যাবে পরবর্তী পদক্ষেপ কী হবে।

---

## ৪. স্বয়ংক্রিয় deploy (GitHub Actions) সেটআপ — ঐচ্ছিক

**লক্ষ্য করুন:** ধাপ ৩-এই Railway আপনার GitHub রিপোর সাথে connect হয়ে
গেছে, তাই **এমনিতেই প্রতিবার `main`-এ push করলে Railway স্বয়ংক্রিয়ভাবে
নতুন deploy করবে** — এই ধাপ ছাড়াই। এই ধাপটা শুধু তখনই দরকার যদি আপনি
GitHub Actions ট্যাবেও deploy-এর status/log দেখতে চান (যেমন আগে GCP
Actions workflow-এ ছিল)।

যদি এটা চান:

1. Railway dashboard-এ, উপরে ডানদিকে আপনার প্রোফাইল আইকনে ক্লিক করুন →
   **"Account Settings"** → **"Tokens"** ট্যাব → **"Create Token"**।
2. টোকেনটা কপি করুন (এটা শুধু একবারই দেখানো হবে)।
3. GitHub repo-তে যান → **Settings → Secrets and variables → Actions →
   New repository secret**:

   | Secret নাম | মান |
   |---|---|
   | `RAILWAY_TOKEN` | ধাপ ২-এ কপি করা টোকেন |

4. এরপর থেকে `main`-এ যেকোনো push GitHub Actions ট্যাবেও দেখানো একটা
   deploy চালাবে (`.github/workflows/deploy.yml` দেখুন)।

⚠️ যদি এই ধাপ স্কিপ করেন, কোনো সমস্যা নেই — Railway-র নিজস্ব GitHub
ইন্টিগ্রেশনই যথেষ্ট deploy-এর জন্য। এটা শুধু বাড়তি visibility দেয়।

---

## ৫. যাচাই করা

```bash
# Health check (ব্রাউজারেও খোলা যাবে)
curl <YOUR_RAILWAY_URL>/health

# একটা টেস্ট জব সাবমিট করা
curl -X POST <YOUR_RAILWAY_URL>/jobs/analysis \
  -H "Content-Type: application/json" \
  -d '{"project_id": "test-123", "analysis_type": "linear-static", "model_payload": {"elements": []}}'

# ফেরত আসা job_id দিয়ে status চেক
curl <YOUR_RAILWAY_URL>/jobs/<JOB_ID>
```

`curl` না থাকলে [Postman](https://www.postman.com/) বা মোবাইলের কোনো
HTTP client অ্যাপ দিয়েও এই তিনটা রিকোয়েস্ট পাঠানো যায়।

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

---

## পরবর্তী ধাপ: Vercel-এর সাথে যুক্ত করা

Next.js অ্যাপ থেকে এই সার্ভিস কল করতে, Vercel প্রজেক্টে একটা environment
variable যোগ করুন:

1. Vercel dashboard → আপনার প্রজেক্ট → **Settings → Environment Variables**
2. যোগ করুন: `SOLVER_SERVICE_URL` = আপনার Railway URL
   (যেমন `https://civilos-structural-solver-production.up.railway.app`)
3. Next.js কোডে (API route বা server action-এ) এই env var ব্যবহার করে
   `fetch(`${process.env.SOLVER_SERVICE_URL}/jobs/analysis`, ...)` কল
   করা যাবে।
