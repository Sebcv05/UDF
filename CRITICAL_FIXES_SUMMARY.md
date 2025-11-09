# Critical Fixes for Thermal Breakup Issues

## Date: 2025-11-09

## Problem Identified

1. **No breakup events after 3.5e-5 s** (expected around 2e-5 s)
2. **Suspected cause**: Evaporation cooling droplets too fast
3. **Concern**: RPE solver allowing T_drop < T_sat (unphysical)

## Solutions Implemented

### Fix 1: Temperature Check in RPE_euler.c

**Added subcooling check:**
```c
// After Euler step, check if droplet cooled below saturation
CONVERGE_precision_t T_sat_check = T_satNH3(P_amb);
if (state.T_drop < T_sat_check) {
    // Droplet subcooled - stop bubble growth
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    old_parcel_cloud->pbt[p_idx] = 0;
    old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    return;
}
```

**Why this helps:**
- Prevents unphysical bubble growth in subcooled droplets
- If evaporation cools droplet below T_sat, bubble growth stops
- More physically consistent

### Fix 2: Disable Evaporation During Thermal Breakup (spray_evap.c)

**Added at start of parcel loop (line ~670):**
```c
// Skip evaporation for parcels undergoing thermal breakup
if(parcel_cloud.pbt[i_pc] == 1)
{
    // Zero out evaporation rate for all species
    for(CONVERGE_index_t isp = 0; isp < num_parcel_species; isp++)
    {
        parcel_cloud.drdt[i_pc * num_parcel_species + isp] = 0.0;
        parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = 0.0;
    }
    
    printf("[EVAP_SKIP] Evaporation disabled for thermal breakup parcels\n");
    continue;  // Skip to next parcel
}
```

**Why this helps:**
- **Thermal breakup is MUCH faster than evaporation** (~microseconds vs milliseconds)
- Evaporation was competing with bubble growth, cooling droplets too fast
- Now thermal breakup takes precedence
- Droplet stays hot long enough for bubble to grow and trigger breakup

## Expected Behavior After Fixes

### What you should see:

1. **`[EVAP_SKIP]`** message when first parcel enters thermal breakup
2. **`[RPE_DEBUG]`** message when tracking starts
3. **`[BREAKUP]`** messages when breakup events occur

### Timeline:
```
t = 0         : Injection
t = 1-5 µs    : Parcels heat up, enter thermal breakup (pbt=1)
              → [EVAP_SKIP] message appears
              → [RPE_DEBUG] starts tracking oldest parcel
t = 10-20 µs  : Bubble grows (RPE solver)
              → No evaporation cooling
              → Temperature stays high
t = 20-30 µs  : kb > threshold
              → [BREAKUP] messages appear
```

## Debugging Output

### Check for these messages:

```bash
# Run simulation and filter output:
./run.sh 2>&1 | grep -E "\[RPE_DEBUG\]|\[BREAKUP\]|\[EVAP_SKIP\]"

# Should see:
[EVAP_SKIP] Evaporation disabled for parcels in thermal breakup (pbt=1)
[RPE_DEBUG] Opened log file: rpe_parcel_debug.csv
[RPE_DEBUG] Tracking parcel with lifetime = 1.234567e-05 s
[BREAKUP] Cycle 150, Rank 0, Time 2.345678e-05 s: 1 breakups (Total: 1)
[BREAKUP] Cycle 151, Rank 0, Time 2.456789e-05 s: 2 breakups (Total: 3)
```

### If still no breakups:

1. **Check RPE debug file:**
   ```bash
   python plot_rpe_debug.py
   ```
   - Is R growing?
   - Is T_drop staying high or dropping?
   - Is Rdot > 0?

2. **Check temperatures:**
   - Is T_drop > T_sat throughout?
   - If not, evap is still too strong (may need to adjust)

3. **Check DGRE/breakup criterion:**
   - Is kb approaching threshold?
   - What is kb_threshold value?

## Rollback if Needed

```bash
git checkout e507d35  # Before these fixes
# Or
git checkout Inertial_RPE  # Original code
```

## Files Modified

- `src/RPE_euler.c` - Added T_drop < T_sat check
- `src/spray_evap.c` - Disabled evap for pbt=1 parcels

## Commit Hash

7080848

## Next Steps

1. Recompile: `./upc2.sh`
2. Run simulation: `./run.sh`
3. Monitor output for [BREAKUP] messages
4. Plot debug data: `python plot_rpe_debug.py`
5. Report results

## Physical Justification

**Why disable evaporation during thermal breakup?**

- Thermal breakup timescale: **~10-50 µs** (rapid)
- Evaporation timescale: **~1-10 ms** (slow)
- Bubble growth velocity: **10-30 m/s** (explosive)
- Evaporation cooling rate: **10³-10⁴ K/s** (gradual)

**During the brief thermal breakup event:**
- Bubble pressure drives rapid growth
- Heat transfer to bubble interface dominates
- External evaporation is negligible by comparison
- Setting drdt=0 for pbt=1 is physically justified

**After breakup:**
- Child droplets have pbt=0
- Evaporation resumes normally
- No impact on post-breakup behavior

