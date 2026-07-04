# CivilOS Structural — Compute Microservice
# Railway এর জন্য optimized: PORT env var respect করে, slim base image ব্যবহার করে।
#
# Phase 4-এ যখন C++ solver (Eigen/OpenBLAS) যোগ হবে, তখন এই Dockerfile-এ
# একটা build stage যোগ করতে হবে যেটা C++ কম্পাইল করবে এবং binary টা
# এই ইমেজে কপি করবে। এখন Phase 0-তে শুধু Python layer।

FROM python:3.12-slim

WORKDIR /service

# System dependencies — Phase 4-এ C++ বিল্ড টুলচেইন (build-essential, cmake,
# libeigen3-dev) এখানে যোগ হবে।
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

# Railway PORT env var দিয়ে পোর্ট বলে দেয়; ডিফল্ট 8080 রাখা হলো লোকাল টেস্টের জন্য।
ENV PORT=8080
EXPOSE 8080

CMD exec uvicorn app.main:app --host 0.0.0.0 --port ${PORT}
