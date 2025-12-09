# Stuck Superheated Parcel Fix - December 9, 2024

## Problem
Parcels observed persisting long into the simulation (0.1+ ms) without breaking up despite being superheated:
- `is_child == 0` (parent parcels)
- `T_drop = 280-290 K` 
- `P_amb = 1 bar` → `T_sat = 254 K`
- **Superheat:** ΔT = 26-36 K (should trigger thermal breakup)

These parcels had incorrect flag combinations preventing entry into the thermal breakup routine.

## Root Cause
For a parcel to enter thermal breakup, it needs:
1. `is_child == 0` (must be parent)
2. `thermal_breakup_flag < 0` (must be in active breakup mode, typically -1)
3. `pbt == 1` (pre-breakup tag must be enabled)

Stuck parcels had either:
- `pbt == 0` (pre-breakup disabled)
- `thermal_breakup_flag >= 0` (breakup disabled)

## Existing Related Logic

### Clause 1: Disable breakup for subcooled parcels (lines 415-428)
```c
if (P_sat_new < P_amb) {
    // Not superheated - disable thermal breakup for this parcel
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Not superheated");
    continue;
}
```
**Purpose:** Catches parcels that have **cooled below** saturation and disables thermal breakup.

### Clause 2: Force large stuck parents to children (lines 884-901)
```c
if (old_parcel_cloud.is_child[p_idx] == 0 && 
    old_parcel_cloud.radius[p_idx] > 70e-6) {
    // Convert to child
    old_parcel_cloud.is_child[p_idx] = 1;
    old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;
}
```
**Purpose:** Catches large parents that never entered thermal breakup and converts to children.

## New Solution: Re-enable Thermal Breakup (INVERSE of Clause 1)

### File Modified: `/home/apollo19/Desktop/Dan_B/UDF/src/spray_drop_distort_NH3.c`

**Location:** Lines 430-457 (immediately after "not superheated" check)

**Change Added:**
```c
// INVERSE CHECK: Re-enable thermal breakup for stuck superheated parcels
// Catch parcels that ARE superheated but have wrong flags preventing thermal breakup
if (old_parcel_cloud.is_child[p_idx] == 0 && 
    P_sat_new >= P_amb &&  // Parcel IS superheated
    old_parcel_cloud.lifetime[p_idx] > 1.0e-5 && // Has existed for > 10 μs
    (old_parcel_cloud.pbt[p_idx] == 0 || 
     old_parcel_cloud.thermal_breakup_flag[p_idx] >= 0)) {
    // This parcel should be in thermal breakup but isn't - reset flags to enable it
    static int stuck_superheat_count = 0;
    if (stuck_superheat_count < 10) {
        printf("[STUCK_SUPERHEATED_FIX] p_idx=%li, superheated but not in thermal breakup - resetting flags\n", p_idx);
        printf("                         P_sat=%.3e > P_amb=%.3e Pa, T_drop=%.2f K, lifetime=%.3e s\n",
               P_sat_new, P_amb, old_parcel_cloud.temp[p_idx], old_parcel_cloud.lifetime[p_idx]);
        printf("                         OLD: pbt=%d, tbt=%d, tbf=%d, R=%.3e m\n",
               old_parcel_cloud.pbt[p_idx], old_parcel_cloud.tbt[p_idx],
               old_parcel_cloud.thermal_breakup_flag[p_idx], old_parcel_cloud.radius[p_idx]);
        stuck_superheat_count++;
    }
    // Reset flags to enable thermal breakup entry
    old_parcel_cloud.pbt[p_idx] = 1;  // Enable pre-breakup tag
    old_parcel_cloud.tbt[p_idx] = 0;  // Reset thermal breakup tag
    old_parcel_cloud.thermal_breakup_flag[p_idx] = -1;  // Enable active breakup mode
    if (stuck_superheat_count <= 10) {
        printf("                         NEW: pbt=%d, tbt=%d, tbf=%d (ready for thermal breakup)\n",
               old_parcel_cloud.pbt[p_idx], old_parcel_cloud.tbt[p_idx],
               old_parcel_cloud.thermal_breakup_flag[p_idx]);
    }
}
```

## Logic Breakdown

### Conditions for Flag Reset:
1. **`is_child == 0`** - Must be a parent parcel (not already broken up)
2. **`P_sat >= P_amb`** - Must be superheated (opposite of subcooled check)
3. **`lifetime > 1.0e-5 s`** - Has existed for > 10 μs (avoid newly injected parcels)
4. **`pbt == 0 OR thermal_breakup_flag >= 0`** - Has wrong flags preventing thermal breakup

### Actions Taken:
1. **`pbt = 1`** - Enable pre-breakup tag (required for thermal breakup entry)
2. **`tbt = 0`** - Reset thermal breakup tag (indicates not yet broken)
3. **`thermal_breakup_flag = -1`** - Set to active breakup mode (enables processing)

### Why 10 μs Lifetime Threshold?
- Newly injected parcels start with `pbt=1, thermal_breakup_flag=-1` (correct flags)
- Only older parcels that have had flags changed incorrectly need this fix
- 10 μs is ~10-100 CFD timesteps (enough to identify truly stuck parcels)
- Prevents unnecessary flag resets on fresh parcels

## Expected Behavior After Fix

### During Simulation:
1. First 10 stuck superheated parcels will print `[STUCK_SUPERHEATED_FIX]` messages
2. Messages show old vs new flag values
3. Parcels will **enter thermal breakup** on the next iteration (once flags are reset)

### Log Output Example:
```
[STUCK_SUPERHEATED_FIX] p_idx=12345, superheated but not in thermal breakup - resetting flags
                         P_sat=1.90e+06 > P_amb=1.00e+05 Pa, T_drop=285.00 K, lifetime=1.234e-04 s
                         OLD: pbt=0, tbt=0, tbf=0, R=5.000e-05 m
                         NEW: pbt=1, tbt=0, tbf=-1 (ready for thermal breakup)
```

### Physical Impact:
- Stuck superheated parcels will now undergo thermal breakup
- Bubble growth and breakup will proceed normally
- Prevents large superheated parcels from persisting unrealistically
- Should reduce the number of large parcels at late times (> 0.1 ms)

## Comparison with Existing Fixes

| Fix | Condition | Action | Purpose |
|-----|-----------|--------|---------|
| **Subcooled (line 415)** | P_sat < P_amb | Set `is_child=1`, disable breakup | Stop breakup for cooled parcels |
| **Stuck large (line 884)** | R > 70 μm, never entered | Set `is_child=1`, disable breakup | Convert stuck parents |
| **NEW: Stuck superheat** | P_sat ≥ P_amb, wrong flags | Reset `pbt=1, tbf=-1` | **Enable breakup** for superheated |

**Key difference:** This fix **enables** thermal breakup rather than disabling it.

## Testing Checklist

After compiling:
1. ✅ Check for `[STUCK_SUPERHEATED_FIX]` messages in log
2. ✅ Verify flags are reset: `OLD: pbt=0/tbf≥0` → `NEW: pbt=1, tbf=-1`
3. ✅ Confirm parcels enter thermal breakup after flag reset
4. ✅ Monitor parcel counts at late times (should decrease)
5. ✅ Check Tecplot: fewer large parcels with `is_child=0` at t > 0.1 ms

## Compilation

From case directory:
```bash
cd /home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev/
./upc2.sh
```

## Notes

- **Non-invasive:** Only resets flags, doesn't change physics
- **Targeted:** Only affects stuck superheated parent parcels
- **Safe:** 10 μs lifetime check prevents affecting fresh parcels
- **Logged:** First 10 instances provide diagnostic information
- **Complementary:** Works alongside existing stuck parcel fixes

## Summary

✅ **Added:** Inverse check that re-enables thermal breakup for stuck superheated parcels  
✅ **Location:** Lines 430-457 in spray_drop_distort_NH3.c  
✅ **Action:** Resets `pbt=1, tbt=0, thermal_breakup_flag=-1`  
✅ **Trigger:** Superheated parent with wrong flags and lifetime > 10 μs  
✅ **Impact:** Enables thermal breakup instead of converting to child

**Date:** 2024-12-09  
**Author:** Based on user observation of stuck parcels in Tecplot
