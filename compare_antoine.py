import math

# Test both sets of coefficients at 273 K
T = 273.0
P_amb = 2e5  # Pa

print("="*70)
print("Antoine Coefficient Comparison at T = 273 K, P_amb = 2 bar")
print("="*70)

# Set 1: From euler_explicit.py
A1 = 6.67956
B1 = 1002.711
C1 = 25.215
log10_P_kPa_1 = A1 - B1/(T - C1)
P_sat_1 = 10**log10_P_kPa_1 * 1000  # kPa to Pa
superheat_1 = P_sat_1 - P_amb

print("\nSet 1 (euler_explicit.py):")
print(f"  A = {A1}, B = {B1}, C = {C1}")
print(f"  P_sat(273 K) = {P_sat_1/1e5:.3f} bar")
print(f"  Superheat = {superheat_1/1e5:.3f} bar")

# Set 2: From test_rpe.c
A2 = 9.96268
B2 = 1617.9
C2 = 6.65
log10_P_kPa_2 = A2 - B2/(T - C2)
P_sat_2 = 10**log10_P_kPa_2 * 1000  # kPa to Pa
superheat_2 = P_sat_2 - P_amb

print("\nSet 2 (test_rpe.c - NIST):")
print(f"  A = {A2}, B = {B2}, C = {C2}")
print(f"  P_sat(273 K) = {P_sat_2/1e5:.3f} bar")
print(f"  Superheat = {superheat_2/1e5:.3f} bar")

print("\n" + "="*70)
print("Temperature sweep comparison:")
print("="*70)
print(f"{'T (K)':<10} {'Set 1 (bar)':<15} {'Set 2 (bar)':<15} {'Ratio':<10}")
print("-"*70)

for T in [273, 283, 293, 303, 313, 323]:
    log10_P1 = A1 - B1/(T - C1)
    P1 = 10**log10_P1  # bar
    
    log10_P2 = A2 - B2/(T - C2)
    P2 = 10**log10_P2  # bar
    
    ratio = P2 / P1
    print(f"{T:<10} {P1:<15.3f} {P2:<15.3f} {ratio:<10.1f}")

print("\nConclusion:")
print("Set 1 gives MUCH LOWER P_sat values (correct for NH3)")
print("Set 2 gives MUCH HIGHER P_sat values (possibly for a different substance or unit system)")

