"""
test_response_spectrum.py — response_spectrum.py এর জন্য unit test।
চালানো: python3 test_response_spectrum.py
"""

import os
import sys

# app/ থেকে সরাসরি চালালে কাজ করে (app/ নিজেই sys.path এ থাকে), কিন্তু
# pytest root থেকে চালালে "app" একটা package হিসেবে treat হয় ও app/
# ভিতরের bare import ("response_spectrum" মডিউল top-level এ খোঁজে)
# ব্যর্থ হয়। এই লাইন app/ কে নিজেই path এ যোগ করে (idempotent), যাতে
# দুই সময়েই (standalone ও pytest) bare import কাজ করে।
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from response_spectrum import build_bnbc_2020_spectrum, get_zone_ss_s1

testsPassed = 0
testsFailed = 0


def check(name: str, condition: bool) -> None:
    global testsPassed, testsFailed
    if condition:
        print(f"  [PASS] {name}")
        testsPassed += 1
    else:
        print(f"  [FAIL] {name}")
        testsFailed += 1
        # pytest collection এর সময় module-level কোড import হিসেবে চলে —
        # raised AssertionError ছাড়া pytest বুঝতে পারবে না কিছু fail
        # করেছে (check()/check_close() শুধু print/counter, exception না)।
        # standalone চালানোর output/exit code অপরিবর্তিত থাকে (নিচের
        # __main__ guard দেখুন)।
        raise AssertionError(f"check failed: {name}")


def check_close(name: str, actual: float, expected: float, tol: float) -> None:
    rel_error = abs(actual - expected) / abs(expected) if abs(expected) > 1e-12 else abs(actual - expected)
    if rel_error < tol:
        print(f"  [PASS] {name} = {actual:.6f} (expected {expected:.6f}, rel error = {rel_error:.6f})")
        global testsPassed
        testsPassed += 1
    else:
        print(f"  [FAIL] {name} = {actual:.6f} (expected {expected:.6f}, rel error = {rel_error:.6f})")
        global testsFailed
        testsFailed += 1
        raise AssertionError(f"check_close failed: {name} (actual={actual}, expected={expected})")


print("=== Test 1: Spectrum Shape — Rise, Plateau, Decay ===")
spectrum = build_bnbc_2020_spectrum("3", "SC", num_points=100, max_period_sec=3.0)
values = [p["spectralAccelerationG"] for p in spectrum]
periods = [p["periodSec"] for p in spectrum]

peak_idx = values.index(max(values))
check("Spectrum has a positive peak", values[peak_idx] > 0)
check("Values before peak are non-decreasing (rising branch)",
      all(values[i] <= values[i + 1] + 1e-9 for i in range(peak_idx)))
check("Values after peak (past plateau) eventually decay (last < peak)",
      values[-1] < values[peak_idx])
check("First point (T=0) is positive but less than peak (short-period rising branch)",
      0 < values[0] < values[peak_idx])

print("\n=== Test 2: Higher Seismic Zone Gives Higher Spectrum Everywhere ===")
spectrum_zone1 = build_bnbc_2020_spectrum("1", "SC", num_points=50, max_period_sec=2.0)
spectrum_zone4 = build_bnbc_2020_spectrum("4", "SC", num_points=50, max_period_sec=2.0)
all_higher = all(
    z4["spectralAccelerationG"] >= z1["spectralAccelerationG"]
    for z1, z4 in zip(spectrum_zone1, spectrum_zone4)
)
check("Zone 4 (highest hazard) spectrum >= Zone 1 (lowest) at every period", all_higher)

print("\n=== Test 3: Softer Site Class Amplifies the Spectrum ===")
spectrum_rock = build_bnbc_2020_spectrum("3", "SA", num_points=50, max_period_sec=2.0)
spectrum_soft = build_bnbc_2020_spectrum("3", "SE", num_points=50, max_period_sec=2.0)
peak_rock = max(p["spectralAccelerationG"] for p in spectrum_rock)
peak_soft = max(p["spectralAccelerationG"] for p in spectrum_soft)
check("Soft soil (SE) peak Sa > rock (SA) peak Sa (site amplification)", peak_soft > peak_rock)

print("\n=== Test 4: Output Points Are Sorted by Period (solver requirement) ===")
spectrum_check = build_bnbc_2020_spectrum("2", "SD", num_points=30, max_period_sec=2.5)
periods_check = [p["periodSec"] for p in spectrum_check]
check("Periods are strictly ascending (solver.h interpolation precondition)",
      all(periods_check[i] < periods_check[i + 1] for i in range(len(periods_check) - 1)))

print("\n=== Test 5: Invalid Inputs Raise Clear Errors ===")
try:
    build_bnbc_2020_spectrum("99", "SC")
    check("Invalid zone raises ValueError", False)
except ValueError:
    check("Invalid zone raises ValueError", True)

try:
    build_bnbc_2020_spectrum("3", "SZ")
    check("Invalid site class raises ValueError", False)
except ValueError:
    check("Invalid site class raises ValueError", True)

try:
    build_bnbc_2020_spectrum("3", "SC", num_points=1)
    check("num_points=1 raises ValueError", False)
except ValueError:
    check("num_points=1 raises ValueError", True)

print("\n=== Test 6: get_zone_ss_s1 Returns Consistent Values Used Internally ===")
ss, s1 = get_zone_ss_s1("4")
check("Zone 4 has highest Ss among all zones", ss == max(get_zone_ss_s1(z)[0] for z in ("1", "2", "3", "4")))
check("Zone 4 has highest S1 among all zones", s1 == max(get_zone_ss_s1(z)[1] for z in ("1", "2", "3", "4")))

print("\n========================================")
print(f"Results: {testsPassed} passed, {testsFailed} failed")
print("========================================")

# __main__ guard: pytest এই ফাইল import করলে sys.exit() রেজ হয়ে
# INTERNALERROR দিত (pytest নিজের exit flow না এমন SystemExit ক্র্যাশ
# হিসেবে ধরে) — standalone চালানোর (python3 test_response_spectrum.py)
# সময়ই শুধু এই কল সক্রিয়।
if __name__ == "__main__":
    import sys
    sys.exit(1 if testsFailed > 0 else 0)
