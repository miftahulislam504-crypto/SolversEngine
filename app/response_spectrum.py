"""
BNBC 2020 Design Response Spectrum Generator
==============================================

Response Spectrum Analysis (RSA, cpp/src/solver.cpp এর
solveResponseSpectrum()) এর ইনপুট হলো একটা তালিকা — period (T) বনাম
spectral acceleration (Sa, g একক) — যা এখানে BNBC 2020 (Part 6,
Chapter 2, Section 2.5.7.1-2) এর design spectrum সূত্র অনুযায়ী তৈরি করা
হয়।

এই ফাইল ইচ্ছাকৃতভাবে C++ solver থেকে আলাদা রাখা হয়েছে (solver.h এর
solveResponseSpectrum() docstring এও ব্যাখ্যা করা) — solver শুধু generic
tabulated (T, Sa) point নেয়, কোনো নির্দিষ্ট design code জানে না। এতে
ভবিষ্যতে ASCE 7 বা Eurocode 8 এর spectrum ও একই solver দিয়ে চালানো
যাবে, শুধু এই স্তরে একটা নতুন generator ফাংশন যোগ করে।

BNBC 2020 spectrum shape (তিনটা অংশ):
  ১. T < T0 (short-period rising branch):
       Sa(T) = S_DS * [0.4 + 0.6*(T/T0)]
  ২. T0 <= T <= TS (constant-acceleration plateau):
       Sa(T) = S_DS
  ৩. TS < T <= TL (constant-velocity, 1/T decay):
       Sa(T) = S_D1 / T
  (T > TL, constant-displacement/1/T² region, BNBC 2020 এ সাধারণত
  বাংলাদেশের জন্য প্রাসঙ্গিক building period range এর বাইরে বলে এই
  সংস্করণে অন্তর্ভুক্ত করা হয়নি — TL কে flat-continue করা হয়েছে, যা
  conservative)

যেখানে:
  S_DS = (2/3) * S_MS,  S_MS = Fa * Ss
  S_D1 = (2/3) * S_M1,  S_M1 = Fv * S1
  T0 = 0.2 * (S_D1 / S_DS)
  TS = S_D1 / S_DS

Ss, S1 (short-period ও 1-second period mapped spectral acceleration,
g একক) BNBC 2020 Table 6.2.13 থেকে Bangladesh-এর সিসমিক জোন অনুযায়ী
আসে। Fa, Fv (site coefficient) BNBC 2020 Table 6.2.16/6.2.17 থেকে site
class ও Ss/S1 এর সাপেক্ষে আসে।

সততার সাথে সীমাবদ্ধতা: BNBC 2020-এর পূর্ণাঙ্গ Table 6.2.16/6.2.17 এ Fa/Fv
আসলে Ss/S1 এর val অনুযায়ী নানা ধাপে ভিন্ন (non-linear table lookup,
প্রতিটা site class এ ৫টা ধাপ)। এখানে ব্যবহৃত মান টেবিলের কাছাকাছি single-
representative মান (moderate Ss/S1 range এর জন্য প্রযোজ্য মান) — খুব
উঁচু বা খুব নিচু Ss/S1 এ প্রকৃত BNBC টেবিল অনুযায়ী ভিন্ন Fa/Fv হতে পারে।
frontend/src/lib/loads/seismicLoad.ts এর getSiteAmplificationFactor() ও
একই ধরনের সরলীকরণ ব্যবহার করে (সেখানে ELF এর জন্য, এখানে RSA এর জন্য) —
এই দুই মডিউল ইচ্ছাকৃতভাবে একে অপরের থেকে independent রাখা হয়েছে (ভিন্ন
ভাষা/প্রসেস, কোনো shared import সম্ভব না এই architecture এ), তাই কোনো
future BNBC টেবিল আপডেট উভয় জায়গায় আলাদাভাবে করতে হবে — এটা একটা known
duplication, single-source-of-truth না।
"""

from __future__ import annotations

SeismicZone = str  # "1" | "2" | "3" | "4" — frontend এর SeismicZone টাইপের সাথে সঙ্গতিপূর্ণ
SiteClass = str    # "SA" | "SB" | "SC" | "SD" | "SE"


# BNBC 2020 Table 6.2.13 এর কাছাকাছি — zone অনুযায়ী Ss (short-period)
# ও S1 (1-second period) mapped spectral acceleration, g একক। এই মান
# frontend এর seismicLoad.ts এর Zone Coefficient Z (0.12/0.2/0.28/0.36)
# এর সাথে সামঞ্জস্যপূর্ণভাবে scale করা (Ss ≈ 2.5*Z এর কাছাকাছি অনুপাত,
# BNBC এর সাধারণ mapping)।
_ZONE_SS_S1: dict[SeismicZone, tuple[float, float]] = {
    "1": (0.30, 0.12),
    "2": (0.50, 0.20),
    "3": (0.70, 0.28),
    "4": (0.90, 0.36),
}

# BNBC 2020 Table 6.2.16 (Fa) ও Table 6.2.17 (Fv) এর একটা representative
# (moderate Ss/S1 range) subset — সীমাবদ্ধতা উপরে docstring এ ব্যাখ্যা করা।
_SITE_FA: dict[SiteClass, float] = {
    "SA": 0.8,
    "SB": 1.0,
    "SC": 1.2,
    "SD": 1.6,
    "SE": 2.5,
}
_SITE_FV: dict[SiteClass, float] = {
    "SA": 0.8,
    "SB": 1.0,
    "SC": 1.6,
    "SD": 2.4,
    "SE": 3.5,
}


def get_zone_ss_s1(zone: SeismicZone) -> tuple[float, float]:
    """একটা সিসমিক জোনের (Ss, S1) — mapped spectral acceleration (g একক) ফেরত দেয়।"""
    if zone not in _ZONE_SS_S1:
        raise ValueError(f"Unknown seismic zone: {zone!r} — must be one of {sorted(_ZONE_SS_S1)}")
    return _ZONE_SS_S1[zone]


def build_bnbc_2020_spectrum(
    seismic_zone: SeismicZone,
    site_class: SiteClass,
    num_points: int = 60,
    max_period_sec: float = 4.0,
) -> list[dict[str, float]]:
    """
    BNBC 2020 design response spectrum (T, Sa) point তালিকা তৈরি করে —
    civilos_solver.solve_response_spectrum() এর `spectrum` argument এ
    সরাসরি পাস করার জন্য প্রস্তুত shape এ
    ([{"periodSec": ..., "spectralAccelerationG": ...}, ...])।

    num_points: T=0 থেকে max_period_sec পর্যন্ত কতগুলো sample point —
    বেশি point মানে solver এর piecewise-linear interpolation বেশি
    নিখুঁতভাবে আসল curve অনুসরণ করবে (curve নিজেই piecewise-linear/
    piecewise-1/T, তাই খুব বেশি point দরকার নেই — 60 বেশিরভাগ building
    period range (T0 এর নিচে থেকে TL পর্যন্ত) এ যথেষ্ট রেজোলিউশন দেয়)।

    max_period_sec: সাধারণত বাংলাদেশের building period range (নিচু-
    থেকে-মাঝারি-উঁচু ভবন) 4 সেকেন্ডের মধ্যেই থাকে — এর বেশি লম্বা period
    এর মোড থাকলে (খুব নমনীয়/উঁচু ভবন), caller এই argument বাড়িয়ে দিতে
    পারে।
    """
    if num_points < 2:
        raise ValueError("num_points must be at least 2")
    if max_period_sec <= 0:
        raise ValueError("max_period_sec must be positive")

    ss, s1 = get_zone_ss_s1(seismic_zone)

    if site_class not in _SITE_FA:
        raise ValueError(f"Unknown site class: {site_class!r} — must be one of {sorted(_SITE_FA)}")
    fa = _SITE_FA[site_class]
    fv = _SITE_FV[site_class]

    sms = fa * ss
    sm1 = fv * s1
    sds = (2.0 / 3.0) * sms
    sd1 = (2.0 / 3.0) * sm1

    if sds < 1e-9:
        # শূন্য বা প্রায়-শূন্য seismic hazard (তাত্ত্বিকভাবে সম্ভব যদি
        # কেউ ভুলবশত অসামঞ্জস্যপূর্ণ ইনপুট দেয়) — একটা flat zero spectrum
        # ফেরত দেওয়া হচ্ছে, division-by-zero (T0 হিসাবে) এড়াতে।
        return [
            {"periodSec": 0.0, "spectralAccelerationG": 0.0},
            {"periodSec": max_period_sec, "spectralAccelerationG": 0.0},
        ]

    t0 = 0.2 * (sd1 / sds)
    ts = sd1 / sds

    def sa_at(t: float) -> float:
        if t < 1e-9:
            return sds * 0.4  # T=0 এ formula অনুযায়ী 0.4*S_DS (T0 branch এর T=0 limit)
        if t < t0:
            return sds * (0.4 + 0.6 * (t / t0))
        if t <= ts:
            return sds
        # T > ts: 1/T decay, TL এর পরেও flat-continue করা হচ্ছে (docstring
        # এ ব্যাখ্যা করা conservative সরলীকরণ)
        return sd1 / t

    points: list[dict[str, float]] = []
    for i in range(num_points):
        t = (i / (num_points - 1)) * max_period_sec
        points.append({"periodSec": t, "spectralAccelerationG": sa_at(t)})

    return points
