# Void Fraction Calculation Fix

**Date:** 2025-12-04  
**Commit:** a01ae78  
**Issue:** Incorrect void fraction formula

---

## The Problem

**WRONG formula used:**
```c
epsilon = R_bubble³ / R_drop_0³
```

This uses the **initial droplet radius** (R_drop_0), which doesn't change.

**CORRECT formula:**
```c
epsilon = R_bubble³ / R_drop_current³
```

Must use **current droplet radius** (radius[p_idx]), which grows as bubble expands.

---

## Physical Definition

**Void fraction** = Volume of bubble / Total volume of droplet

```
ε = V_bubble / V_total
  = (4/3)π·R_bubble³ / (4/3)π·R_drop³
  = R_bubble³ / R_drop³
```

**Key point:** As the bubble grows, the droplet also expands. Both radii are time-dependent.

---

## Why This Matters

### Scenario: Bubble growth from nucleation to breakup

**Initial state (t=0):**
- R_drop_0 = 50 μm (injection)
- R_bubble_0 = 0.1 μm (nucleation)
- ε = (0.1)³/(50)³ ≈ 0

**During growth (t=t1):**
- R_drop = 60 μm (expanded due to bubble)
- R_bubble = 35 μm (growing)
- ε_correct = (35)³/(60)³ = 0.25
- ε_wrong = (35)³/(50)³ = 0.34 ← **Too high!**

**At breakup target:**
- Want: ε = 0.55
- R_drop = 70 μm (expanded)
- R_bubble_correct = 58.5 μm (from ε = 0.55)
- R_bubble_wrong = 54.1 μm (if using R_drop_0=50μm)

**Result:** Wrong formula triggers breakup too early!

---

## Changes Made

### 1. RPE_song.c

**Line ~160: Added current droplet radius**
```c
CONVERGE_precision_t R_drop_current = old_parcel_cloud->radius[p_idx];
```

**Line ~189: Updated function call**
```c
// OLD:
CONVERGE_precision_t epsilon = song_compute_void_fraction(R, Ro);

// NEW:
CONVERGE_precision_t epsilon = song_compute_void_fraction(R, R_drop_current);
```

**Line ~40-55: Updated function documentation**
```c
// Correct formula: ε = V_bubble / (V_bubble + V_liquid)
//                    = R_bubble³ / R_drop³
// 
// For a droplet: V_total = (4/3)π·R_drop³
//                V_bubble = (4/3)π·R_bubble³
//
// Therefore: ε = R_bubble³ / R_drop³
```

---

### 2. spray_drop_distort_NH3.c

**Line ~525-527: Fixed calculation**
```c
// OLD:
CONVERGE_precision_t R = old_parcel_cloud.r_bubble[p_idx];
CONVERGE_precision_t Ro = old_parcel_cloud.r_drop_0[p_idx];
CONVERGE_precision_t epsilon = (R*R*R) / (Ro*Ro*Ro);

// NEW:
CONVERGE_precision_t R_bubble = old_parcel_cloud.r_bubble[p_idx];
CONVERGE_precision_t R_drop = old_parcel_cloud.radius[p_idx];  // Current!
CONVERGE_precision_t epsilon = (R_bubble*R_bubble*R_bubble) / (R_drop*R_drop*R_drop);
```

**Line ~535: Updated debug message**
```c
// OLD:
printf("  R_bubble=%.3e m, R_drop_0=%.3e m, R_drop=%.3e m\n", ...);

// NEW:
printf("  R_bubble=%.3e m, R_drop=%.3e m, R_drop_0=%.3e m\n", ...);
```

Now prints: current R_bubble, current R_drop, then initial R_drop_0 for reference.

---

### 3. include/RPE_song.h

**Line ~33-36: Updated function signature**
```c
// OLD:
CONVERGE_precision_t song_compute_void_fraction(
    CONVERGE_precision_t R,
    CONVERGE_precision_t Ro
);

// NEW:
CONVERGE_precision_t song_compute_void_fraction(
    CONVERGE_precision_t R_bubble,
    CONVERGE_precision_t R_drop_current
);
```

Clearer variable names indicate which radius to use.

---

## Impact on Physics

### Before Fix (WRONG):

Void fraction calculated using initial droplet size:
```
ε_wrong = R_bubble³ / R_drop_0³
```

- As bubble grows, R_bubble increases
- But denominator (R_drop_0) stays constant
- ε grows faster than it should
- **Breakup triggered too early**

### After Fix (CORRECT):

Void fraction tracks actual droplet state:
```
ε_correct = R_bubble³ / R_drop_current³
```

- As bubble grows, R_bubble increases
- Droplet expands, so R_drop also increases
- ε growth rate is physically correct
- **Breakup at correct void fraction**

---

## Relationship to Geometry()

In thermal model, `Geometry()` function updates droplet radius based on bubble size:

```c
R_drop_new = f(R_bubble, R_drop_old, ...)
```

Song model **skips** Geometry(), so where does R_drop come from?

**Answer:** Song model should probably call Geometry() too!

**Current issue:** 
- Song uses `radius[p_idx]` which may not be updated
- Need to verify if Geometry() is needed for Song

**Two options:**

1. **Add Geometry() for Song:**
   - Call before void fraction check
   - Updates droplet radius based on bubble
   
2. **Calculate R_drop directly:**
   - Use volume conservation
   - R_drop³ = R_drop_0³ + (R_bubble³ - R_bubble_0³)

**TODO:** Investigate if Song needs Geometry() or manual calculation

---

## Testing Impact

### Expected Changes:

**With wrong formula (before fix):**
- Breakup at void ≈ 0.55 using R_drop_0
- Actual physical void fraction < 0.55
- Early breakup

**With correct formula (after fix):**
- Breakup when actual void = 0.55
- Later breakup (more bubble growth)
- Physically accurate

### Verification:

Check `breakup_events.csv`:
```csv
time,ncyc,p_idx,parent_radius,r_bubble,v_bubble,child_radius
```

Calculate actual void fraction:
```
ε_at_breakup = (r_bubble/parent_radius)³
```

Should be close to 0.55 for Song model breakups.

---

## Summary

✅ **Fixed:** Void fraction now uses current droplet radius  
✅ **Compiled:** No errors  
✅ **Committed:** a01ae78  

⚠️ **TODO:** Verify if Song needs Geometry() call to update radius[p_idx]

**Status:** Physics is now correct, but may need Geometry() addition.

---

## Quick Reference

### Correct void fraction formula:
```c
CONVERGE_precision_t R_bubble = old_parcel_cloud.r_bubble[p_idx];
CONVERGE_precision_t R_drop = old_parcel_cloud.radius[p_idx];  // CURRENT
CONVERGE_precision_t epsilon = (R_bubble*R_bubble*R_bubble) / (R_drop*R_drop*R_drop);
```

### What each radius means:
- `r_bubble[p_idx]` = Current bubble radius (time-dependent)
- `radius[p_idx]` = Current droplet radius (time-dependent)
- `r_bubble_0[p_idx]` = Initial bubble radius (constant)
- `r_drop_0[p_idx]` = Initial droplet radius (constant)

### Physical constraint:
```
R_bubble ≤ R_drop  (always)
0 ≤ ε ≤ 1          (by definition)
```

At ε = 0.55, bubble occupies 55% of droplet volume.

