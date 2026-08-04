# CivilOS Structural — Compute Microservice
# Multi-stage build: প্রথম stage এ C++ solver কম্পাইল হয় (build
# toolchain সহ), দ্বিতীয় (final) stage এ শুধু compiled .so ফাইল ও
# Python runtime থাকে — এভাবে final image size ছোট থাকে (cmake/g++/
# Eigen headers ইত্যাদি production এ বহন করার দরকার নেই)।

# ============================================================
# Stage 1: C++ Build
# ============================================================
FROM python:3.12-slim AS cpp-builder

WORKDIR /build

# C++ build toolchain + Eigen (header-only, apt দিয়ে সহজে পাওয়া যায়)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

# pybind11 — CMakeLists.txt এটা `python3 -m pybind11 --cmakedir` দিয়ে
# খুঁজে পায়, তাই build stage এও pip install করা প্রয়োজন।
RUN pip install --no-cache-dir pybind11

COPY cpp ./cpp

RUN cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && \
    cmake --build cpp/build --target civilos_solver -j$(nproc)

# ============================================================
# Stage 2: Final Runtime Image
# ============================================================
FROM python:3.12-slim

WORKDIR /service

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Compiled Python extension module (.so) কে সরাসরি Python এর
# site-packages এ কপি করা হচ্ছে — এভাবে app/main.py তে কোনো path
# manipulation ছাড়াই `import civilos_solver` কাজ করবে (main.py এর
# fallback sys.path insertion cpp/build/ খুঁজবে, কিন্তু production এ
# এই কপি থাকায় সেই fallback এর দরকারই পড়বে না, এটাই primary path)।
COPY --from=cpp-builder /build/cpp/build/civilos_solver*.so /usr/local/lib/python3.12/site-packages/

COPY app ./app

# Render/Cloud Run/Railway PORT env var দিয়ে পোর্ট বলে দেয়; ডিফল্ট
# 8080 রাখা হলো লোকাল টেস্টের জন্য।
ENV PORT=8080
EXPOSE 8080

# Container startup এই সময় একটা quick sanity check চালানো হচ্ছে —
# যদি civilos_solver import ব্যর্থ হয় (কোনো build সমস্যার কারণে),
# container সাথে সাথে স্পষ্ট error দিয়ে exit করবে, uvicorn চালু হয়ে
# পরে প্রতিটা request এ 503 দেওয়ার বদলে — এটা deployment সমস্যা দ্রুত
# ধরতে সাহায্য করে (container log এই এরর সরাসরি দেখা যাবে)।
#
# গুরুত্বপূর্ণ: কোনো sys.path manipulation ছাড়াই import করা হচ্ছে
# ইচ্ছাকৃতভাবে — .so ফাইল site-packages এ কপি করা হয়েছে উপরে, যা
# ইতিমধ্যে Python এর default import path এ আছে, তাই এই check ঠিক
# সেই একই কন্ডিশনে চলে যেভাবে production এ uvicorn চলার সময় import
# ঘটবে। কোনো path hint দিলে এই check false-positive পাস করতে পারত
# even যদি site-packages কপি আসলে ব্যর্থ হয়ে থাকে।
RUN python3 -c "import civilos_solver; print('civilos_solver import check: OK')"

CMD exec uvicorn app.main:app --host 0.0.0.0 --port ${PORT}
