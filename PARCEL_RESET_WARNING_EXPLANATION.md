# PARCEL_RESET_WARNING Explanation

## What the Warning Means

```
[PARCEL_RESET_WARNING] Invalid r_drop_0=8.250e-05 m, using current radius 3.611e-05 m (p_idx=1, reason: Not superheated)
```

This warning appears when the `reset_parcel_to_child()` function is called and detects an inconsistency with the `r_drop_0` variable.

## Root Cause

### What is `r_drop_0`?

`r_drop_0` is meant to store the **initial droplet radius at injection** for each parcel. However, for **child parcels** (created by breakup), this variable gets set differently:

**From `parcel_prop.c` line 286:**
```c
// In case of KH-RT breakup
// R_D_0 - should not be used for child parcels
parcel_cloud.r_drop_0[passed_child_parcel_idx] = parcel_cloud.radius[passed_parent_parcel_idx];
```

So for child parcels, `r_drop_0` stores the **parent's radius at breakup**, NOT the child's initial radius.

### The Problem Sequence

1. **Parent parcel injected:**
   - `r_drop_0 = 82.5 μm` (injection radius)
   - `radius = 82.5 μm` (current radius)

2. **KH-RT breakup occurs (somewhere in simulation):**
   - Creates child parcels
   - Child's `r_drop_0 = 82.5 μm` (parent radius, NOT child radius)
   - Child's `radius = 10-40 μm` (actual child size, smaller!)

3. **Parcel evaporates and cools:**
   - Current `radius = 26-40 μm` (shrunk from evaporation)
   - Temperature drops below saturation → "Not superheated"
   - `reset_parcel_to_child()` is called

4. **Reset function tries to restore "injection" radius:**
   - Wants to set `radius = r_drop_0 = 82.5 μm`
   - But current radius is only 26-40 μm
   - **Fails safety check:** `R_new > R_old * 2.0` (82.5 > 2 × 36)

5. **Warning triggered:**
   - Safety check prevents expanding parcel artificially
   - Uses current radius instead of invalid `r_drop_0`

## Why This Happens

The `reset_parcel_to_child()` function is designed for **parent parcels** that:
1. Were injected at some radius (e.g., 82.5 μm)
2. Grew a bubble but never broke up
3. Cooled below saturation → need to be reset

But it's being called on **child parcels** (already broken up) that:
1. Were created at smaller radius (e.g., 35 μm)
2. Have `r_drop_0` pointing to parent's radius (82.5 μm)
3. This creates a mismatch

## Is This a Problem?

**Short answer: No, the safety check prevents any real issue.**

### What the Safety Check Does (line 36-45 in parcel_reset.c):
```c
// Safety check: if injection radius is invalid, use current radius
if (R_new <= 0.0 || R_new > R_old * 2.0) {
    printf("[PARCEL_RESET_WARNING] ...\n");
    R_new = R_old;  // Keep current radius
    Nd_new = Nd_old;  // Keep current num_drop
}
```

**Result:**
- Parcel radius stays at current value (36 μm, not 82.5 μm) ✅
- Number of drops stays constant (mass conserved) ✅
- Bubble is still zeroed out ✅
- Flags still set correctly (`is_child=1`, `pbt=0`, `tbf=999`) ✅

## Why "Not superheated"?

The warning shows `reason: Not superheated`, which means:

1. These are parcels that were in or approaching thermal breakup
2. They cooled below saturation temperature during evaporation
3. The subcooled check (line 415-428 in `spray_drop_distort_NH3.c`) caught them:
   ```c
   if (P_sat_new < P_amb) {
       // Not superheated - disable thermal breakup
       reset_parcel_to_child(&old_parcel_cloud, p_idx, "Not superheated");
   }
   ```

This is **expected behavior** for parcels that:
- Evaporate rapidly in a cool environment
- Lose enough heat to drop below saturation
- Should exit the thermal breakup routine

## Physical Interpretation

**These parcels:**
1. Were created as children from KH-RT breakup (small, ~30-40 μm)
2. Cooled during evaporation (T_drop < T_sat(P_amb))
3. Are being correctly disabled from thermal breakup
4. Safety check prevents artificially expanding them back to parent size

**This is correct physics!** Small droplets evaporate faster and cool below saturation.

## Summary

| Item | Status | Explanation |
|------|--------|-------------|
| **Warning severity** | ⚠️ INFO | Not an error, just diagnostic |
| **Mass conservation** | ✅ OK | Safety check maintains current mass |
| **Parcel behavior** | ✅ OK | Correctly disabled from thermal breakup |
| **Flags** | ✅ OK | Properly set to child mode |
| **Physics** | ✅ OK | Small droplets cooling is expected |

## Why You're Seeing Many of These

At `ncyc=2300` (t ≈ 157 μs):
- Many child parcels have been created from breakup
- These are evaporating and cooling
- As they drop below saturation, they trigger the subcooled check
- Each one generates this warning (first 5 only, then silenced)

**This is normal for a spray simulation** where:
- Ambient temperature < droplet boiling point
- Small droplets evaporate and cool rapidly
- Thermal breakup is only relevant for hot, superheated parcels

## Should You Worry?

**No.** The warning is informational only and the safety check ensures correct behavior.

If you want to **silence the warnings**, you could:
1. Reduce the warning count from 5 to 1 in `parcel_reset.c` line 38
2. Or remove the printf entirely (but keep the safety check logic)

## Alternative: Fix `r_drop_0` for Child Parcels

If you want to eliminate the root cause, modify `parcel_prop.c` line 314:
```c
// Change from:
parcel_cloud.r_drop_0[passed_child_parcel_idx] = parcel_cloud.radius[passed_child_parcel_idx];

// To:
parcel_cloud.r_drop_0[passed_child_parcel_idx] = 0.0;  // Mark as child (no injection radius)
```

Then in `parcel_reset.c`, add check:
```c
if (R_new <= 0.0 || R_new > R_old * 2.0) {
    R_new = R_old;  // Already does this
    Nd_new = Nd_old;
}
```

But this isn't necessary since the current safety check already handles it correctly.

---

**Conclusion:** These warnings are **expected and harmless**. The safety check ensures parcels behave correctly even when `r_drop_0` is mismatched for child parcels.
