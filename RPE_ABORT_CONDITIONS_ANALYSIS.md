# RPE Euler Solver: Abort Conditions Analysis
**Date:** 2025-11-11  
**File:** `src/RPE_euler.c`

## Summary

The RPE_euler solver has **5 abort conditions** that set `thermal_breakup_flag=999`, `pbt=0`, and `v_bubble=0`, effectively stopping thermal breakup for that parcel.

---

## Abort Condition #1: Droplet Too Small (Line 215)

```c
if (params.Ro < 1.0e-9)  // 1 nm
{
    printf("[RPE_ERROR] Droplet radius too small: Ro=%.3e m, skipping RPE solver\n", params.Ro);
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    old_parcel_cloud->pbt[p_idx] = 0;
    old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    return;
}
```

**When triggered:** Droplet radius < 1 nm  
**Likelihood for early parcels:** **VERY LOW** (initial droplets are ~80 µm)  
**Purpose:** Prevent division by zero / numerical issues for vanished droplets

---

## Abort Condition #2: Negative Pressure Difference (Line 295) ⚠️ **LIKELY CULPRIT**

```c
CONVERGE_precision_t P_sat;
Saturation_PressureNH3(Td, &P_sat);
if ((P_sat - P_amb) < 0.0) {
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    old_parcel_cloud->pbt[p_idx] = 0;
    old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    return;
}
```

**When triggered:** P_sat(T_drop) < P_amb (parcel is subcooled)  
**Likelihood for early parcels:** **HIGH**  
**Why this affects early parcels:**

### Physics:
- **Early parcels:** Tinj = 293K, Tamb = 298K
- Parcel is **already subcooled** relative to ambient at injection
- P_sat(293K) ≈ 8.49 bar, P_amb = 2.0 bar
- Initially P_sat > P_amb ✓ (barely superheated)
- But if parcel cools even slightly (evaporation, expansion), P_sat < P_amb → **ABORT**

### For later parcels:
- Gas field has warmed up from earlier evaporation
- Parcels heat up faster
- Maintain P_sat > P_amb longer
- Can undergo thermal breakup

**This is the most likely reason early parcels are aborting RPE!**

---

## Abort Condition #3: Negative Rdot (Line 317)

```c
if (state.Rdot < 0.0) {
    printf("[RPE_STOP] Negative Rdot=%.3e, stopping growth (bubble collapsing)\n", state.Rdot);
    // ... more diagnostics ...
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    old_parcel_cloud->pbt[p_idx] = 0;
    old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    return;
}
```

**When triggered:** Bubble wall velocity becomes negative (collapse)  
**Likelihood:** **MEDIUM** - can happen if:
- Pressure drops suddenly
- Bubble overcools
- Numerical instability

**Connection to early parcels:**  
If P_sat drops below P_amb + 2σ/R + viscous term, Rddot becomes negative, leading to negative Rdot.

---

## Abort Condition #4: Subcooled After Integration (Line 346) ⚠️ **ALSO LIKELY**

```c
CONVERGE_precision_t T_sat_check = T_satNH3(P_amb);
if (state.T_drop < T_sat_check) {
    printf("[RPE_STOP] Subcooled: T_drop=%.2f K < T_sat=%.2f K, stopping bubble growth\n",
           state.T_drop, T_sat_check);
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    old_parcel_cloud->pbt[p_idx] = 0;
    old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    return;
}
```

**When triggered:** T_drop < T_sat(P_amb)  
**Likelihood for early parcels:** **HIGH**

### Critical values:
- T_sat(2 bar) ≈ 239.7 K (for ammonia)
- T_drop(early parcel) = 293 K → **SHOULD PASS**
- But wait... this check uses T_satNH3(P_amb), not comparing to injection temp

**Actually, this should NOT trigger for 293K parcels at 2 bar ambient** because T_sat(2 bar) = 239.7 K << 293 K.

So this is likely NOT the issue unless P_amb is much higher locally.

---

## Abort Condition #5: Velocity Too Small (Line 373)

```c
if (state.Rdot < 1.0e-10)  // 0.1 nm/s
{
    printf("[RPE_STOP] Bubble velocity too small: Rdot=%.3e < 1e-10, stopping\n",
           state.Rdot);
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    old_parcel_cloud->pbt[p_idx] = 0;
    old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    return;
}
```

**When triggered:** Rdot < 0.1 nm/s (essentially stopped)  
**Likelihood:** **MEDIUM-HIGH** for subcooled parcels

**This is a consequence, not a root cause:**
- If driving pressure is low (P_sat ≈ P_amb)
- Bubble grows very slowly
- Eventually Rdot < threshold
- RPE aborts

---

## Root Cause Analysis

### Why Early Parcels Abort:

**Chain of events:**
1. Parcel injected at 293K into 298K ambient
2. P_sat(293K) = 8.49 bar, P_amb = 2.0 bar
3. ΔP = 6.49 bar initially (small but positive)
4. Bubble nucleates at R0 = 1.1*Rc ≈ 72 nm
5. Parcel cools slightly from:
   - Evaporation heat loss (spray_evap active? CHECK THIS!)
   - Expansion work
   - Convective cooling to gas phase
6. Temperature drops to ~292.5K
7. P_sat(292.5K) ≈ 8.25 bar
8. ΔP drops from 6.49 bar to 6.25 bar
9. Bubble growth slows
10. Eventually either:
    - Condition #2: P_sat < P_amb (if cools further)
    - Condition #5: Rdot < 1e-10 (growth stalls)
11. RPE aborts, sets pbt=0, thermal_breakup_flag=999
12. Parcel never undergoes thermal breakup

### Why Later Parcels Succeed:

1. Gas field has warmed from earlier parcel evaporation
2. Parcels heat up faster
3. P_sat increases, ΔP remains large
4. Bubble grows faster
5. Rdot stays above threshold
6. kb reaches 1.0 before RPE aborts
7. Breakup occurs normally

---

## Diagnostic Output to Check

Look for these messages in simulation log:

```
[RPE_ERROR] Droplet radius too small
[RPE_STOP] Negative Rdot
[RPE_STOP] Subcooled: T_drop=...
[RPE_STOP] Bubble velocity too small
[RPE_STOP] Bubble hit droplet edge
```

Also check:
```bash
grep "thermal_breakup_flag = 999" log
grep "pbt.*= 0" log
```

---

## Potential Fixes

### Fix 1: Disable evaporation for parent parcels (ALREADY DONE)
**File:** `src/spray_evap.c`

This prevents evaporative cooling from killing thermal breakup before it completes.

**Status:** ✅ Already implemented in CRITICAL_FIXES_SUMMARY.md

---

### Fix 2: Relax Rdot threshold (Line 373)

**Current:**
```c
if (state.Rdot < 1.0e-10)  // 0.1 nm/s
```

**Proposed:**
```c
if (state.Rdot < 1.0e-12)  // 0.001 nm/s (100x more tolerant)
```

**Rationale:** 
- For subcooled parcels, bubble may grow VERY slowly but still reach kb=1.0 eventually
- Current threshold may be too aggressive
- Allow slower growth as long as direction is positive

---

### Fix 3: Remove/modify P_sat < P_amb check (Line 295)

**Current:**
```c
if ((P_sat - P_amb) < 0.0) {
    // Abort
}
```

**Option A - Add tolerance:**
```c
if ((P_sat - P_amb) < -1e3) {  // Allow 1 kPa subcooling
    // Abort only if significantly subcooled
}
```

**Option B - Remove entirely:**
Let the physics (Rdot < 0) handle subcooling naturally.

**Rationale:**
- Parcels injected at 293K into 298K ambient are marginally superheated
- Small fluctuations shouldn't abort thermal breakup immediately
- Bubble pressure can temporarily exceed P_sat during growth (non-equilibrium)

---

### Fix 4: Increase injection temperature

**Not a code fix, but:**
- Inject at 300K or 305K instead of 293K
- Provides more superheat margin
- More robust against small temperature drops

---

## Recommended Actions

### 1. Check if evaporation is disabled for pbt=1 parcels
```bash
grep "pbt.*== 1" src/spray_evap.c
```

Expected: Should see evaporation being skipped when pbt=1.

### 2. Check diagnostic output
```bash
grep -E "RPE_ERROR|RPE_STOP" converge.log
```

This will show which abort condition is triggering most frequently.

### 3. Temporarily add more diagnostics to RPE_euler.c

Before line 295, add:
```c
printf("[RPE_CHECK] p_idx=%d, T_drop=%.2f K, P_sat=%.3e Pa, P_amb=%.3e Pa, Delta=%.3e Pa\n",
       p_idx, Td, P_sat, P_amb, P_sat - P_amb);
```

This will show exactly when and why parcels are aborting.

### 4. Implement Fix #2 (relax Rdot threshold) first

**This is the safest fix:**
- Line 373: Change `1.0e-10` to `1.0e-12`
- Allows slower bubble growth
- Still aborts if truly stalled (Rdot < 0.001 nm/s)
- Minimal risk

### 5. If Fix #2 doesn't help, implement Fix #3A (P_sat tolerance)

**Moderate risk:**
- Line 295: Change `< 0.0` to `< -1e3`
- Allows 1 kPa subcooling before aborting
- More physically realistic (bubbles can be slightly metastable)

---

## Expected Behavior After Fixes

**Before fixes:**
- Early parcels: thermal_breakup_flag = 999 (aborted)
- kb values: 0.06-0.08 (never reaches 1.0)
- No breakup events

**After fixes:**
- Early parcels: thermal_breakup_flag = -1 → 3 (normal progression)
- kb values: 0.06 → 0.08 → ... → 1.0+ (reaches threshold)
- Breakup events occur for all superheated parcels

---

## Conclusion

**Root cause:** Early parcels are triggering RPE abort conditions due to:
1. Marginal superheat (293K injection, 298K ambient)
2. Small temperature drops from evaporation/cooling
3. Overly strict thresholds in RPE_euler.c

**Primary suspect:** Line 373 - `Rdot < 1e-10` threshold too high for slow growth

**Secondary suspect:** Line 295 - No tolerance for P_sat ≈ P_amb conditions

**Recommended fix order:**
1. Relax Rdot threshold (Line 373: 1e-10 → 1e-12)
2. Add P_sat tolerance (Line 295: < 0.0 → < -1e3)
3. Verify evaporation is disabled for pbt=1
4. Add diagnostic output to confirm

