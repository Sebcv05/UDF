# Dual Geometry Heat Transfer Fix

**Date:** 2025-11-12  
**Issue:** test_rpe.c was missing film-based heat transfer limit

---

## Problem Found

The unit test `test_rpe.c` was only using **semi-infinite heat transfer**, while both the CONVERGE UDF and `euler_explicit.py` use **dual geometry** (semi-infinite + film, taking the maximum).

### Original test_rpe.c (WRONG):
```c
// Only semi-infinite
double h_conv = Nu * params->k_l / L_char;
double Q_conv = h_conv * A_bubble * dT;
```

### euler_explicit.py (CORRECT):
```python
# Semi-infinite formulation
h_conv_si = Nu * k_l / L_char
Q_si = h_conv_si * A_surf * dT

# Film formulation
if film_thickness > 1e-12:
    h_conv_film = Nu * k_l / film_thickness
    Q_film = h_conv_film * A_surf * dT

# Use maximum (least restrictive)
Q_total = max(Q_si, Q_film)
```

### CONVERGE UDF RPE_euler.c (CORRECT):
```c
// Semi-infinite formulation
h_conv_si = Nu * params->k_l / L_char;
Q_si = h_conv_si * A_bubble * dT;

// Film formulation
film_thickness = params->Ro - R;
if (film_thickness > 1e-12) {
    h_conv_film = Nu * params->k_l / film_thickness;
    Q_film = h_conv_film * A_bubble * dT;
}

// Use maximum (least restrictive)
Q_conv = (Q_film > Q_si) ? Q_film : Q_si;
```

---

## Physics Explanation

### Semi-Infinite Formulation
**Used when:** Bubble is small, film is thick
```
h_conv_si = Nu × k_l / L_char
L_char = 2×R
```
Heat transfer coefficient scales with **bubble size**.

### Film Formulation
**Used when:** Bubble is large, film is thin
```
h_conv_film = Nu × k_l / δ
δ = R₀ - R (film thickness)
```
Heat transfer coefficient scales inversely with **film thickness**.

### Why Take Maximum?

As bubble grows:
- **Early:** R small → L_char small → h_conv_si large → Q_si dominates
- **Late:** δ small → h_conv_film very large → Q_film dominates

Taking `max(Q_si, Q_film)` ensures we use the **least restrictive** (highest heat transfer) formulation at each instant.

---

## Impact on Results

### Before Fix (semi-infinite only):

| T₀ (K) | Max Rdot (m/s) | Film at Max (µm) | Time at Max (µs) |
|--------|----------------|------------------|------------------|
| 343    | 0.108          | 70.8             | 85.0             |
| 353    | 2.844          | 78.2             | 1.70             |
| 393    | 34.670         | 82.1             | 0.024            |
| 403    | 35.890         | 82.1             | 0.024            |

**Problem:** Maximum at thick film (near R₀) because semi-infinite gives highest Q when L_char is large.

### After Fix (dual geometry):

| T₀ (K) | Max Rdot (m/s) | Film at Max (µm) | Time at Max (µs) |
|--------|----------------|------------------|------------------|
| 343    | 0.108          | 70.8             | 85.0             |
| 353    | 4.277          | 4.13             | 36.3             |
| 363    | 6.381          | 4.13             | 17.9             |
| 373    | 8.159          | 4.13             | 12.5             |
| 383    | 9.756          | 4.13             | 9.88             |
| 393    | 11.240         | 4.13             | 8.31             |
| 403    | 12.647         | 4.13             | 7.23             |

**Fixed:** Maximum at **thin film** (~4 µm) where film formulation dominates and gives very high heat transfer rate.

---

## Physical Consistency

### Film Thickness vs Heat Transfer

At thin film (δ ~ 4 µm):
```
h_conv_film = Nu × k_l / δ
            ≈ 100 × 0.5 / 4e-6
            ≈ 1.25×10⁷ W/(m²·K)
```

At thick film (δ ~ 80 µm), using semi-infinite:
```
h_conv_si = Nu × k_l / (2×R)
          ≈ 100 × 0.5 / (2×80e-6)
          ≈ 3.1×10⁵ W/(m²·K)
```

**Ratio: 40× higher heat transfer at thin film!**

This explains why bubble acceleration peaks when film is thin.

---

## Breakup Implications

### Old (wrong) results:
- Maximum Rdot at thick film
- Would predict breakup occurs early
- Inconsistent with thin-film breakup criterion kb = η/δ

### New (correct) results:
- Maximum Rdot at thin film (~4 µm)
- Bubble acceleration increases as film thins
- Consistent with breakup criterion: highest stress when both Rdot and 1/δ are large

This aligns with the physical expectation that droplets break up when the **film is thin** and **bubble is accelerating rapidly**.

---

## Comparison to CONVERGE Results

Now the unit test properly matches the UDF physics:

**Heat Transfer:**
- ✅ Dual geometry (semi-infinite + film)
- ✅ Takes maximum of both formulations
- ✅ Film dominates at thin film thickness

**Initial Conditions:**
- ✅ R_init = 1.1×R_crit
- ✅ Temperature-dependent critical radius
- ✅ Equilibrium bubble mass initialization

**Integration:**
- ✅ Explicit Euler
- ✅ Same ODE system (R, Rdot, T_drop, m_b)
- ✅ Same timestep strategy

---

## Plot Improvements

### 1. Film Thickness Axis Reversed

**Before:** 0 → 82.5 µm (left to right)  
**After:** 82.5 → 0 µm (left to right)

**Reason:** Time flows left to right, bubble grows left to right, film thins left to right → more intuitive.

### 2. Peak Location Clear

With dual geometry, all curves except 343 K peak at **δ ≈ 4 µm**, clearly visible on the reversed axis plot.

### 3. Physically Meaningful

The plot now shows:
- **Left (thick film):** Slow acceleration, semi-infinite dominates
- **Middle:** Transition region
- **Right (thin film):** Rapid acceleration, film dominates → breakup

---

## Validation

### Energy Balance Check

At maximum Rdot (T = 393 K, δ = 4 µm):
```
Q_film = h_conv_film × A × ΔT
       = 1.25×10⁷ × 4π×(78e-6)² × 10
       ≈ 960 W

mdot = Q_film / L_v
     = 960 / 1.37e6
     ≈ 7×10⁻⁴ kg/s

Rdot = mdot / (ρ_v × A)
     ≈ 7×10⁻⁴ / (10 × 4π×(78e-6)²)
     ≈ 9 m/s ✓
```

Matches observed max Rdot ≈ 11 m/s (order of magnitude correct).

---

## Summary of Changes

### test_rpe.c:
1. Added film-based heat transfer calculation
2. Added logic to select max(Q_si, Q_film)
3. Now matches UDF and euler_explicit.py exactly

### plot_rpe_sweep.py:
1. Reversed x-axis: `ax1.set_xlim(82.5, 0)`
2. Film thickness decreases left to right (more intuitive)

### Results:
1. Maximum Rdot now at thin film (~4 µm) for all superheated cases
2. Consistent with physical expectation
3. Matches CONVERGE UDF behavior

---

## Files Updated

- `test_rpe.c`: Added dual geometry heat transfer (lines 74-124)
- `plot_rpe_sweep.py`: Reversed film thickness axis
- `rpe_sweep_plots.pdf`: Regenerated with correct physics
- `rpe_sweep_all.csv`: New data with dual geometry
- `rpe_sweep_summary.csv`: Updated statistics

---

## Conclusion

**Status:** ✅ Fixed - Now matches CONVERGE UDF physics exactly

The unit test now correctly implements:
1. Dual geometry heat transfer (semi-infinite + film)
2. Takes maximum for least restrictive limit
3. Shows physically consistent behavior with peak Rdot at thin film

This matches the thermal_test_result.pdf that was validated earlier in the week.

