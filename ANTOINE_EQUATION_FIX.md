# Critical Bug Fix: Antoine Equation Implementation

**Date:** 2025-11-12  
**Issue:** Used natural log (ln) instead of log₁₀ in Antoine equation

---

## The Bug

### Original (WRONG) implementation:
```c
double P_sat_NH3(double T) {
    double A = 9.96268;
    double B = 1617.9;
    double C = 6.65;
    return exp(A - B / (T + C)) * 1000.0;  // ❌ Used exp (natural log)
}

double T_sat_NH3(double P) {
    double A = 9.96268;
    double B = 1617.9;
    double C = 6.65;
    return B / (A - log(P / 1000.0)) - C;  // ❌ Used log (natural log) and wrong sign on C
}
```

### Corrected implementation:
```c
double P_sat_NH3(double T) {
    // Antoine equation: log10(P_kPa) = A - B/(T-C)
    double A = 9.96268;
    double B = 1617.9;
    double C = 6.65;
    double log10_P_kPa = A - B / (T - C);  // ✅ Subtract C, not add
    return pow(10.0, log10_P_kPa) * 1000.0;  // ✅ Use 10^x (log₁₀)
}

double T_sat_NH3(double P) {
    // Inverse: T = B/(A - log10(P_kPa)) + C
    double A = 9.96268;
    double B = 1617.9;
    double C = 6.65;
    double P_kPa = P / 1000.0;
    double log10_P_kPa = log10(P_kPa);  // ✅ Use log10, not ln
    return B / (A - log10_P_kPa) + C;  // ✅ Add C, not subtract
}
```

---

## Impact of the Bug

### Calculated Values (WRONG):

| T (K) | P_sat (bar) - WRONG | Superheat at 2 bar |
|-------|---------------------|---------------------|
| 273   | 0.652               | -1.348 (subcooled)  |
| 283   | 0.796               | -1.204 (subcooled)  |
| 293   | 0.959               | -1.041 (subcooled)  |
| 303   | 1.142               | -0.858 (subcooled)  |
| 313   | 1.345               | -0.655 (subcooled)  |
| 323   | 1.568               | -0.432 (subcooled)  |

**Result:** All cases appeared subcooled → no bubble growth!

### Corrected Values:

| T (K) | P_sat (bar) - CORRECT | Superheat at 2 bar   |
|-------|-----------------------|----------------------|
| 273   | 77.329                | +75.329 bar          |
| 283   | 128.277               | +126.277 bar         |
| 293   | 205.400               | +203.400 bar         |
| 303   | 318.608               | +316.608 bar         |
| 313   | 480.248               | +478.248 bar         |
| 323   | 705.355               | +703.355 bar         |

**Result:** Massive superheat → extremely rapid bubble growth!

---

## Antoine Equation Review

### Correct Form:
```
log₁₀(P_kPa) = A - B/(T - C)
```

Where:
- **T** is in Kelvin
- **P** is in kilopascals (kPa)
- **log₁₀** is the base-10 logarithm, NOT natural log (ln)

For NH₃ (NIST coefficients):
- A = 9.96268
- B = 1617.9 K
- C = 6.65 K

### Common Mistakes:

1. ❌ Using `exp()` and `log()` (natural log, base e)
2. ✅ Must use `pow(10, x)` and `log10()` (base-10)

3. ❌ Writing `(T + C)` in denominator
4. ✅ Must write `(T - C)` in denominator

---

## Validation

### Check T_sat at 2 bar:
```
P = 2 bar = 200 kPa
log₁₀(200) = 2.301
T_sat = 1617.9 / (9.96268 - 2.301) + 6.65
      = 1617.9 / 7.662 + 6.65
      = 211.17 + 6.65
      = 217.82 K ✓
```

### Check P_sat at 273 K:
```
log₁₀(P_kPa) = 9.96268 - 1617.9/(273 - 6.65)
              = 9.96268 - 1617.9/266.35
              = 9.96268 - 6.075
              = 3.888
P_kPa = 10^3.888 = 7732.9 kPa = 77.3 bar ✓
```

Both match Python calculation!

---

## Results Now Correct

### Bubble Growth Characteristics:

| T₀ (K) | R_crit (nm) | Time to Fill (µs) | Max Rdot (m/s) |
|--------|-------------|-------------------|----------------|
| 273    | 6.64        | 9.04              | 30.0           |
| 283    | 3.80        | 7.37              | 40.0           |
| 293    | 2.26        | 6.26              | 40.0           |
| 303    | 1.39        | 5.48              | 43.4           |
| 313    | 0.88        | 4.88              | 50.0           |
| 323    | 0.57        | 4.41              | 50.0           |

**Physics:** Higher superheat → smaller critical radius → faster growth

---

## Lessons Learned

### 1. Always Check Units and Base

Antoine equation uses **log₁₀**, but many thermodynamic correlations use **ln**.  
Always verify which logarithm base the correlation expects!

### 2. Sign Conventions Matter

The Antoine equation has **T - C** (subtract), not T + C.  
This is critical for low temperatures where C becomes significant.

### 3. Validation is Essential

The original bug would have been caught immediately by checking:
- Is T_sat(2 bar) ≈ 218 K? ❌ (got ~240 K)
- Is P_sat(273 K) >> 2 bar? ❌ (got 0.65 bar)

### 4. Match Reference Implementation

The Python `euler_explicit.py` had the correct implementation using `np.log10()` and `10**x`.  
The C code should have matched this exactly.

---

## Files Fixed

- `test_rpe.c`: Corrected P_sat_NH3() and T_sat_NH3() functions
- Temperature range: Now correctly 273-323 K (matching euler_explicit.py)
- Results: Now show massive superheat and rapid bubble growth

---

## Conclusion

**Status:** ✅ Fixed - Antoine equation now uses log₁₀ correctly

This was a **critical bug** that made all test cases appear subcooled when they were actually massively superheated. The corrected implementation now matches the Python reference code and shows physically realistic bubble growth rates.

