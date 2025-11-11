# Persistent Parcel Fixes - 2025-11-11

## Problem Summary

**Issue:** 10-20 parcels persisting to 0.5 ms without undergoing thermal breakup, despite normal breakup occurring at ~15-20 µs for most parcels.

**Root Cause:** Three abort conditions in `spray_drop_distort_NH3.c` were improperly disabling thermal breakup, causing parcels to persist indefinitely with `pbt=0` but still visible in domain.

---

## Fixes Implemented

### Fix #1: Temperature > 300K Check (Line 328)

**Problem:** Parcel was "removed" but still persisted with `pbt=0`.

**Before:**
```c
if(Td > 300.0){
   printf("\n\n before P_sat, Td = %f, ...");
   printf("\nRemoving Parcel\n");
   old_parcel_cloud.num_drop[p_idx] = 0.0;
   old_parcel_cloud.radius[p_idx] = 0.0;
   // ... sets many fields to 0 ...
   old_parcel_cloud.pbt[p_idx] = 0;
   old_parcel_cloud.thermal_breakup_flag[p_idx] = 0;  // Should be 999
   continue;
}
```

**After:**
```c
if(Td > 300.0){
   static int td_high_count = 0;
   if (td_high_count < 10) {
      printf("[THERMAL_ABORT] p_idx=%li, Td=%.2f K > 300K, removing parcel (lifetime=%.3e s, radius=%.3e m)\n",
             p_idx, Td, old_parcel_cloud.lifetime[p_idx], old_parcel_cloud.radius[p_idx]);
      td_high_count++;
   }
   // ... sets fields to 0 ...
   old_parcel_cloud.pbt[p_idx] = 0;
   old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;  // FIXED: Was 0
   continue;
}
```

**Changes:**
- Added diagnostic output with rate limiting
- Changed `thermal_breakup_flag` from 0 to 999 for consistency
- Cleaner, more informative message

---

### Fix #2: Invalid P_sat Check (Line 366) ⚠️ **CRITICAL**

**Problem:** Did NOT set `pbt=0` or `thermal_breakup_flag=999`, so parcel thought it was still in thermal breakup but never processed.

**Before:**
```c
if (P_sat < 1)
{
   printf("\nParcel Temperature outside of range for P_sat Correlation, continuing...");
   old_parcel_cloud.v_bubble[p_idx] = 0.0;
   // ... sets some fields but NOT pbt or thermal_breakup_flag ...
   continue;
}
```

**After:**
```c
if (P_sat < 1)
{
   static int psat_low_count = 0;
   if (psat_low_count < 10) {
      printf("[THERMAL_ABORT] p_idx=%li, P_sat=%.3e Pa < 1 Pa, Td=%.2f K outside Antoine range (lifetime=%.3e s)\n",
             p_idx, P_sat, Td, old_parcel_cloud.lifetime[p_idx]);
      psat_low_count++;
   }
   old_parcel_cloud.v_bubble[p_idx] = 0.0;
   // ... other field resets ...
   old_parcel_cloud.pbt[p_idx] = 0;  // CRITICAL FIX: Was missing
   old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;  // CRITICAL FIX: Was missing
   continue;
}
```

**Changes:**
- **CRITICAL:** Added `pbt[p_idx] = 0` (was missing!)
- **CRITICAL:** Added `thermal_breakup_flag[p_idx] = 999` (was missing!)
- Added diagnostic output with rate limiting
- More informative message

**This was likely catching the majority of the 10-20 stuck parcels.**

---

### Fix #3: Negative Bubble Radius Check (Line 388)

**Problem:** Logic was correct (sets `pbt=0`), but lacked diagnostic output.

**Before:**
```c
if(old_parcel_cloud.r_bubble[p_idx]<0.0)
{
   old_parcel_cloud.r_bubble[p_idx]=0.0;
   old_parcel_cloud.thermal_breakup_flag[p_idx]=999;
   old_parcel_cloud.pbt[p_idx]=0;
}
```

**After:**
```c
if(old_parcel_cloud.r_bubble[p_idx]<0.0)
{
   static int rbubble_neg_count = 0;
   if (rbubble_neg_count < 10) {
      printf("[THERMAL_ABORT] p_idx=%li, r_bubble<0 (subcooled), disabling thermal breakup (lifetime=%.3e s, Td=%.2f K)\n",
             p_idx, old_parcel_cloud.lifetime[p_idx], Td);
      rbubble_neg_count++;
   }
   old_parcel_cloud.r_bubble[p_idx]=0.0;
   old_parcel_cloud.thermal_breakup_flag[p_idx]=999;
   old_parcel_cloud.pbt[p_idx]=0;
}
```

**Changes:**
- Added diagnostic output with rate limiting
- Logic unchanged (was already correct)

---

### Fix #4: Diagnostic for Stuck Parcels (Line 425)

**Purpose:** Identify any remaining stuck parcels that slip through.

**Added:**
```c
// DIAGNOSTIC: Identify stuck parcels that persisted >100 µs without breaking up
if (old_parcel_cloud.lifetime[p_idx] > 1.0e-4 && 
    old_parcel_cloud.pbt[p_idx] == 0 && 
    old_parcel_cloud.thermal_breakup_flag[p_idx] != 3 && 
    old_parcel_cloud.is_child[p_idx] == 0 &&
    old_parcel_cloud.radius[p_idx] > 5.0e-5) {
   static int stuck_parcel_count = 0;
   if (stuck_parcel_count < 20) {
      printf("[STUCK_PARCEL] p_idx=%li, lifetime=%.3e s, radius=%.3e m, pbt=%d, tbf=%d, Td=%.2f K, r_bubble=%.3e m\n",
             p_idx, old_parcel_cloud.lifetime[p_idx], old_parcel_cloud.radius[p_idx],
             old_parcel_cloud.pbt[p_idx], old_parcel_cloud.thermal_breakup_flag[p_idx], 
             Td, old_parcel_cloud.r_bubble[p_idx]);
      stuck_parcel_count++;
   }
}
```

**Triggers when:**
- Lifetime > 100 µs (0.0001 s)
- pbt = 0 (thermal breakup disabled)
- thermal_breakup_flag ≠ 3 (not a successful kb-triggered breakup)
- is_child = 0 (parent parcel, not child from breakup)
- radius > 50 µm (still a large droplet)

**Output limited to first 20 detections** to avoid log spam.

---

## Expected Results After Fixes

### Diagnostic Output You Should See:

```bash
# If parcels hit the abort conditions, you'll see:
[THERMAL_ABORT] p_idx=42, Td=301.23 K > 300K, removing parcel (lifetime=3.456e-05 s, radius=8.123e-05 m)
[THERMAL_ABORT] p_idx=17, P_sat=8.234e-01 Pa < 1 Pa, Td=250.12 K outside Antoine range (lifetime=2.345e-05 s)
[THERMAL_ABORT] p_idx=99, r_bubble<0 (subcooled), disabling thermal breakup (lifetime=1.234e-05 s, Td=292.45 K)

# If any parcels still get stuck despite fixes:
[STUCK_PARCEL] p_idx=5, lifetime=5.123e-04 s, radius=8.234e-05 m, pbt=0, tbf=999, Td=293.15 K, r_bubble=0.000e+00 m
```

### In Tecplot at 0.5 ms:

**Before fixes:**
- 10-20 large parcels (radius ~80 µm) with pbt=0, thermal_breakup_flag=999 or 0
- Persisting indefinitely

**After fixes:**
- No stuck parcels (or very few if there's another edge case)
- All parent parcels either:
  - Broken up (thermal_breakup_flag=3)
  - Properly removed (num_drop=0, radius=0)
  - Legitimately disabled (with clear diagnostic message)

---

## Rate Limiting Strategy

All diagnostic messages use **static counters** to limit output to first 10-20 occurrences:

```c
static int td_high_count = 0;
if (td_high_count < 10) {
   printf("[THERMAL_ABORT] ...");
   td_high_count++;
}
```

**Why:** Prevents log file explosion if many parcels trigger same condition.

**Result:** You'll see each abort type up to 10 times, giving you good statistics without overwhelming output.

---

## Files Modified

**File:** `/home/apollo19/Desktop/Dan_B/UDF/src/spray_drop_distort_NH3.c`

**Lines changed:**
- 328-347: Fix #1 (Td > 300K)
- 366-380: Fix #2 (P_sat < 1) - **CRITICAL**
- 388-412: Fix #3 (r_bubble < 0)
- 425-442: Fix #4 (Stuck parcel diagnostic)

**Total additions:** ~40 lines of diagnostic/fix code

---

## Testing Plan

### 1. Compile
```bash
cd /path/to/case/directory
./upc2.sh
```

### 2. Run Simulation
```bash
./run.sh
```

### 3. Check for Abort Messages
```bash
grep "THERMAL_ABORT" log | head -20
grep "STUCK_PARCEL" log
```

### 4. Count Abort Types
```bash
grep "\[THERMAL_ABORT\]" log | awk '{print $3}' | sort | uniq -c
```

Expected output:
```
  5 Td=...    (Td > 300K)
 10 P_sat=... (Invalid P_sat) <- Likely the main culprit
  2 r_bubble<0 (Subcooled)
```

### 5. Verify No Stuck Parcels
```bash
grep "STUCK_PARCEL" log
```

Expected: **No output** (or very few if there's an edge case we missed)

### 6. Check Tecplot at 0.5 ms
- Load post file at t=0.5 ms
- Filter parcels by: radius > 50 µm AND is_child = 0
- **Should see zero or very few remaining**

---

## If Issues Persist

### If you still see [STUCK_PARCEL] messages:

1. **Check the diagnostic output** - it will tell you:
   - Which p_idx
   - lifetime, radius, Td
   - pbt and thermal_breakup_flag values
   - r_bubble value

2. **Look for patterns:**
   - All stuck parcels have same thermal_breakup_flag value?
   - Similar temperatures?
   - Similar lifetimes?

3. **Add more diagnostics** at the specific line where that thermal_breakup_flag is set.

### If [THERMAL_ABORT] messages are overwhelming:

Increase the rate limit from 10 to 50 or 100:
```c
if (td_high_count < 100) {  // Was 10
```

---

## Summary

**Primary fix:** Added `pbt=0` and `thermal_breakup_flag=999` to the P_sat < 1 check (Fix #2).

**This was the missing piece** that allowed parcels to persist with thermal breakup "enabled" but never actually processing.

**All three abort conditions now properly:**
1. Set `pbt=0` (disable thermal breakup)
2. Set `thermal_breakup_flag=999` (mark as aborted)
3. Print diagnostic message (for debugging)

**Expected impact:** Should eliminate or drastically reduce the 10-20 persistent parcels you observed at 0.5 ms.

---

## Git Status

- Branch: `v3.1.12`
- Modified: `src/spray_drop_distort_NH3.c`
- Status: Changes not yet committed
- Ready for: `./upc2.sh` compilation

---

## Next Steps

1. Compile with `./upc2.sh`
2. Run simulation with `./run.sh`
3. Check log for `[THERMAL_ABORT]` and `[STUCK_PARCEL]` messages
4. Visualize in Tecplot at 0.5 ms
5. Report results

If stuck parcels still appear, the diagnostic will tell us exactly where to look next!

