# Tiny Parcel Skip Fix - December 5, 2024

## Problem
```
Breakup.c: Error: Denominator in rad_denom calculation is close to zero. Aborting.
radius = 9.483395e-08
rububble = 9.009225e-08
```

Tiny droplets (radius ~ 95 nm, close to bubble radius) were causing division by zero in the Breakup routine when calculating `rad_denom = radius - r_bubble`.

## Root Cause
Parcels that have nearly completely evaporated (radius < 1 μm) were still being processed through the breakup routine, causing numerical issues when:
- Film thickness → 0
- Geometric calculations fail (nearly spherical bubble in tiny droplet)
- Division by zero in `rad_denom`

## Solution

Added **early exit check** at the start of the parcel loop in `spray_drop_distort_NH3.c`.

### File Modified: `/home/apollo19/Desktop/Dan_B/UDF/src/spray_drop_distort_NH3.c`

**Location:** Line 186 (start of parcel loop)

**Change:**
```c
for (int p_idx = 0; p_idx < num_parcels_in_cloud; p_idx++)
{
   // Skip parcels that are too small (avoid numerical issues in breakup routines)
   if (old_parcel_cloud.radius[p_idx] < 1.0e-6) {
      static int skip_count = 0;
      if (skip_count < 10) {
         printf("[DISTORT_SKIP] Parcel too small: p_idx=%d, radius=%.3e m < 1e-6 m (skipping)\n",
                p_idx, old_parcel_cloud.radius[p_idx]);
         skip_count++;
      }
      continue;  // Skip this parcel
   }
   
   // ... rest of loop
}
```

## Why 1e-6 m (1 μm) Threshold?

1. **Physical reasoning:**
   - Initial droplet radius ~ 50 μm
   - When radius < 1 μm, droplet has lost 99.9998% of its volume
   - Breakup is not relevant for such small droplets
   - These parcels are essentially fully evaporated

2. **Numerical stability:**
   - Bubble radius initialization ~ 20-30 nm
   - At R_droplet ~ 100 nm, bubble occupies most of the droplet
   - Film thickness ~ nm scale → geometric calculations unreliable
   - Prevents division by zero in `rad_denom = R_drop - R_bubble`

3. **Computational efficiency:**
   - Skip unnecessary calculations for nearly-evaporated parcels
   - Avoids costly error handling and aborts

## Expected Behavior

1. **During simulation:**
   - First 10 small parcels will print `[DISTORT_SKIP]` messages
   - After that, skips silently to avoid log spam
   - No MPI aborts from tiny parcels

2. **Physical impact:**
   - Negligible: parcels < 1 μm have essentially no mass
   - These parcels would evaporate completely in the next timestep anyway
   - No breakup physics are relevant at this scale

## Testing

Compiled successfully. The error message from your run shows:
- `radius = 9.483395e-08 m` (94.8 nm)
- This is **well below** the 1e-6 m threshold
- Will now be skipped with a diagnostic message

## Alternative Thresholds

If you find 1 μm is too aggressive, you can adjust:
- `1.0e-6` = 1 μm (current, very safe)
- `5.0e-7` = 500 nm (more aggressive)
- `1.0e-7` = 100 nm (very aggressive, only skip extreme cases)

The 1 μm threshold is recommended because:
- Safe: well above bubble radius scale
- Conservative: only skips truly tiny parcels
- Prevents all observed errors (smallest error was at 95 nm)

## Summary

✅ **Added:** Early exit for parcels with radius < 1 μm  
✅ **Location:** Top of parcel loop in spray_drop_distort_NH3.c  
✅ **Impact:** Prevents division by zero in breakup routines  
✅ **Performance:** Negligible (skips only nearly-evaporated parcels)  
✅ **Compiled:** Successfully
