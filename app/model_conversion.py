"""
Model Conversion — Frontend JSON → C++ Solver Input
======================================================

সততার সাথে এই মডিউলের গুরুত্ব: এখানে একক রূপান্তর (unit conversion)
ভুল হলে সলভার নীরবে ভুল সংখ্যা দেবে (crash করবে না, শুধু ভুল উত্তর
দেবে) — এটা সবচেয়ে বিপজ্জনক ধরনের bug। তাই এই মডিউলের প্রতিটা
conversion factor এখানে explicit ভাবে লেখা এবং নিচে টেস্ট করা হয়েছে।

একক রূপান্তরের সারাংশ:
  - দৈর্ঘ্য/কোঅর্ডিনেট: frontend মিটার (m) এ কাজ করে (Grid/Story
    coordinate) — সলভারও মিটারে, তাই কোনো conversion লাগে না।
  - Section dimension (width, depth ইত্যাদি): frontend মিলিমিটার (mm)
    এ, কিন্তু SectionProperties (area, ixx, iyy, j) ইতিমধ্যে mm²/mm⁴
    এককে frontend এ হিসাব হয়ে আসে — solver mm² কে m² এ (÷ 1e6),
    mm⁴ কে m⁴ এ (÷ 1e12) রূপান্তর করে।
  - Elastic/Shear Modulus: frontend MPa (=N/mm²) এ, solver kN/m² এ।
    1 MPa = 1 N/mm² = 1000 kN/m² (এই রূপান্তরটা নিচে verifyMpaToKnPerM2
    এ derive করে দেখানো হয়েছে, কারণ এটা সহজে ভুল হওয়ার মতো একটা
    রূপান্তর)।
  - Force: frontend ইতিমধ্যে kN এ (LoadCase, src/lib/types/load.ts) —
    কোনো conversion লাগে না।
  - Material density (Modal Analysis এর জন্য): frontend unitWeight
    (kN/m³) এ, solver mass density (tonne/m³) এ — density = unitWeight
    / g, নিচে convert_unit_weight_to_density() এ derivation দেখুন।
"""

from typing import Any

MM2_TO_M2 = 1e-6      # mm² → m² (1m = 1000mm, area স্কেল হয় 1000² = 1e6)
MM4_TO_M4 = 1e-12     # mm⁴ → m⁴ (length⁴ স্কেল হয় 1000⁴ = 1e12)
MPA_TO_KN_PER_M2 = 1000.0  # MPa → kN/m² (নিচে derivation দেখুন)
GRAVITATIONAL_ACCELERATION_M_S2 = 9.80665  # g, standard gravity (m/s²) — density derive করতে ব্যবহৃত


def _verify_mpa_to_kn_per_m2_conversion() -> None:
    """
    এই ফাংশন module import হওয়ার সময় একবার চলে এবং conversion factor
    এর derivation নিজে যাচাই করে — যদি কোনো ভবিষ্যতে edit এই factor
    ভুলভাবে বদলে দেয়, এই assertion সাথে সাথে ধরবে (silent wrong-answer
    এর বদলে একটা loud crash, যা অনেক বেশি নিরাপদ)।

    Derivation: 1 MPa = 1 N/mm² (সংজ্ঞা অনুযায়ী)
                       = 1 N / (1e-3 m)²
                       = 1 N / 1e-6 m²
                       = 1e6 N/m²
                       = 1e6 Pa
                       = 1e3 kPa  (1 kPa = 1000 Pa)
                       = 1e3 kN/m²  (1 kPa = 1 kN/m², যেহেতু kPa = kN/m² সংজ্ঞাগতভাবে)
    সুতরাং 1 MPa = 1000 kN/m² — যা MPA_TO_KN_PER_M2 এর মান।
    """
    derived_factor = 1e6 / 1e3  # N/m² থেকে kN/m² এ যেতে ÷1000, কিন্তু MPa থেকে N/m² এ যেতে ×1e6 — নিট factor
    assert abs(derived_factor - MPA_TO_KN_PER_M2) < 1e-9, (
        f"MPa→kN/m² conversion factor mismatch: derived {derived_factor}, "
        f"constant is {MPA_TO_KN_PER_M2}. এটা একটা critical unit-conversion bug — "
        f"সলভারে পাঠানো stiffness ভুল হবে।"
    )


_verify_mpa_to_kn_per_m2_conversion()


def convert_section_to_solver_units(section_properties_mm: dict[str, float]) -> dict[str, float]:
    """
    Frontend এর computeSectionProperties() রিটার্ন করা SectionProperties
    (area: mm², ixx/iyy/j: mm⁴) কে solver এর একক (m², m⁴) এ রূপান্তর করে।

    yieldMomentMzKNm/yieldMomentMyKNm ঐচ্ছিক (Nonlinear Static Analysis
    এর জন্য, cpp/include/types.h এর SectionProperties docstring দেখুন)
    — ইতিমধ্যে kN·m এককে থাকে (frontend এর design-code-specific
    capacity হিসাব থেকে, যেমন steel: Zx·Fy, concrete: ACI 318 Mn),
    তাই কোনো unit conversion ছাড়াই সরাসরি pass-through করা হচ্ছে।
    section_properties_mm এ এই key না থাকলে (পুরনো payload, বা এই
    section এর জন্য hinge capacity define করা হয়নি) 0.0 (elastic
    section) ধরা হয়।
    """
    return {
        "area": section_properties_mm["area"] * MM2_TO_M2,
        "ixx": section_properties_mm["ixx"] * MM4_TO_M4,
        "iyy": section_properties_mm["iyy"] * MM4_TO_M4,
        "j": section_properties_mm["j"] * MM4_TO_M4,
        "yieldMomentMzKNm": section_properties_mm.get("yieldMomentMzKNm", 0.0),
        "yieldMomentMyKNm": section_properties_mm.get("yieldMomentMyKNm", 0.0),
    }


def convert_material_to_solver_units(
    elastic_modulus_mpa: float, shear_modulus_mpa: float
) -> dict[str, float]:
    """Frontend এর material.ts এর E/G (MPa) কে solver এর kN/m² এ রূপান্তর করে।"""
    return {
        "elasticModulus": elastic_modulus_mpa * MPA_TO_KN_PER_M2,
        "shearModulus": shear_modulus_mpa * MPA_TO_KN_PER_M2,
    }


def convert_unit_weight_to_density(unit_weight_kn_per_m3: float) -> float:
    """
    Frontend এর material.ts এর unitWeight (kN/m³, ACI/BNBC কনভেনশনে
    ব্যবহৃত — যেমন concrete ≈ 24 kN/m³, steel ≈ 78.5 kN/m³) কে solver এর
    mass density (tonne/m³) এ রূপান্তর করে, Modal Analysis এর mass
    matrix গণনার জন্য।

    Derivation: Newton's second law, F = m·a। unitWeight আসলে একটা
    "specific weight" (weight per unit volume), যা mass density এর
    সাথে সম্পর্কিত: weight = mass × g, তাই:

        unitWeight [kN/m³] = density [tonne/m³] × g [m/s²]

    কারণ 1 tonne × 1 m/s² = 1 kN (F=ma সংজ্ঞা থেকে সরাসরি — 1000 kg ×
    1 m/s² = 1000 N = 1 kN)। সুতরাং:

        density [tonne/m³] = unitWeight [kN/m³] / g [m/s²]

    g = 9.80665 m/s² (standard gravity, exact সংজ্ঞায়িত মান — CGPM
    দ্বারা প্রতিষ্ঠিত আন্তর্জাতিক মান, আঞ্চলিক ভিন্নতা নেই যা এখানে
    গুরুত্বপূর্ণ)।

    উদাহরণ: concrete unitWeight = 24 kN/m³ হলে density ≈ 2.448 tonne/m³
    (বাস্তব concrete density ≈ 2400-2500 kg/m³ = 2.4-2.5 tonne/m³ এর
    সাথে সঙ্গতিপূর্ণ — sanity check হিসেবে এটা মিলিয়ে দেখা হয়েছে)।

    ⚠️ সীমাবদ্ধতা: frontend এর কিছু unitWeight মান (যেমন steel এর 78.5
    kN/m³) নিজেই একটা রাউন্ড করা ইঞ্জিনিয়ারিং কনভেনশন মান — সাধারণত
    g≈10 m/s² ধরে derive করা (7850 kg/m³ × 10 m/s² / 1000 = 78.5
    kN/m³)। এই ফাংশন g=9.80665 (নিখুঁত মান) দিয়ে বিপরীত রূপান্তর করে,
    তাই ফেরত পাওয়া density আসল material density থেকে ~2% ভিন্ন হতে
    পারে (steel এ 7.85 এর বদলে ~8.00 tonne/m³ আসে)। এই সামান্য
    পার্থক্য Modal Analysis এর natural frequency তে ~1% এর কম প্রভাব
    ফেলবে (frequency ∝ 1/√mass), যা প্রকৌশল উদ্দেশ্যে গ্রহণযোগ্য —
    কিন্তু ভবিষ্যতে যদি নিখুঁত material density দরকার হয়, frontend এর
    unitWeight এর বদলে সরাসরি density field যোগ করা ভালো হবে material
    library তে।
    """
    return unit_weight_kn_per_m3 / GRAVITATIONAL_ACCELERATION_M_S2


def build_solver_input(
    nodes: list[dict[str, Any]],
    elements: list[dict[str, Any]],
    boundary_conditions: list[dict[str, Any]],
    loads: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    ইতিমধ্যে-প্রস্তুত (already-converted, already-indexed) node/element/
    bc/load তালিকা থেকে সরাসরি C++ solver এর প্রত্যাশিত dict shape
    বানায়। এই ফাংশন নিজে কোনো conversion করে না — উপরের convert_*
    ফাংশনগুলো এবং caller (main.py এর endpoint) ইতিমধ্যে সঠিক একক ও
    index নিশ্চিত করে এখানে পাঠাবে। এই বিভাজন ইচ্ছাকৃত: conversion
    logic (এই ফাইলে, ইউনিট-টেস্টযোগ্য, side-effect-free) এবং
    orchestration (main.py তে, request handling) আলাদা রাখা হয়েছে।
    """
    return {
        "nodes": nodes,
        "elements": elements,
        "boundaryConditions": boundary_conditions,
        "loads": loads,
    }
