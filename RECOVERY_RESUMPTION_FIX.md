# Critical Fix: Recovery Resumption Issue

**Date:** 2025-11-17  
**Status:** Fixed - ready for testing

---

## Problem Identified

After implementing time-based recovery, parcels that experienced collapse were NOT resuming bubble growth after the recovery period. They appeared "stuck" with small bubble sizes.

### Root Cause

The issue was a **logic conflict between recovery and the "very small velocity" check**:

1. During recovery, we set `v_bubble = 0.0` to keep the bubble stable
2. On the next timestep, RPE_euler_solver is called with `state.Rdot = 0.0`
3. The code checks: `if (state.Rdot < 1.0e-10)` → TRUE!
4. It then sets `pbt = 0` and `thermal_breakup_flag = 999`, **permanently disabling RPE**
5. Parcel never attempts bubble growth again

### Why This Happened

The sequence of events:

```
Timestep N:   Collapse detected → Recovery applied → v_bubble = 0.0 → flag = 888
Timestep N+1: flag = 888 → spray_drop_distort sets flag = -999 → break sub-loop
Timestep N+2: RPE enters with Rdot = 0.0 → "too small" check triggers → pbt = 0, flag = 999
              ^^^ PARCEL PERMANENTLY DISABLED ^^^
```

The parcel could never restart because it was killed by the small velocity check before it could even check if it should attempt recovery again.

---

## Solution Implemented

### 1. Set Small Positive Velocity During Recovery

Instead of `v_bubble = 0.0`, use `v_bubble = 1.0e-9` m/s:

**In recovery wait period:**
```c
old_parcel_cloud->v_bubble[p_idx] = 1.0e-9;  // Small positive velocity
```

**When applying recovery:**
```c
old_parcel_cloud->v_bubble[p_idx] = 1.0e-9;  // Small positive velocity
```

This velocity is:
- Large enough to pass the `< 1.0e-10` check
- Small enough not to cause significant bubble growth during recovery wait
- Negligible compared to normal bubble growth rates (1-20 m/s)

### 2. Skip "Very Small Velocity" Check During Recovery

Modified the check to skip parcels that are in recovery mode:

```c
// Check for very small velocity (but skip if in recovery mode)
if (state.Rdot < 1.0e-10 && old_parcel_cloud->recovery_time[p_idx] == 0.0) {
    // Kill parcel only if NOT in recovery
    old_parcel_cloud->v_bubble[p_idx] = 0.0;
    old_parcel_cloud->pbt[p_idx] = 0;
    old_parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    return;
}
```

If `recovery_time > 0`, the parcel is either:
- Currently in recovery wait period, OR
- Has recovered but not yet reset the timer (will reset when Rdot > 0)

Either way, we should NOT kill it for having small velocity.

---

## Files Modified

**`src/RPE_euler.c`:**
- Line ~352: Changed `v_bubble = 0.0` → `1.0e-9` (recovery wait)
- Line ~392: Changed `v_bubble = 0.0` → `1.0e-9` (recovery apply)
- Line ~463: Added `&& recovery_time[p_idx] == 0.0` to small velocity check

---

## Expected Behavior After Fix

### Timestep N: Collapse and Recovery
```
[RPE_COLLAPSE] Negative Rdot=-7.470e-03, attempting recovery
[RECOVERY #1] R_bubble: 3.812e-05 -> 1.191e-07 m
              v_bubble set to 1.0e-9 m/s
              recovery_time = 1.234567e-05 s
```

### Timestep N+1 to N+X: Recovery Wait Period (< 50 μs)
```
[RPE_RECOVERY_WAIT] Still in recovery period: 1.234e-05 / 5.000e-05 s elapsed
                    Rdot=-2.345e-03 (negative but waiting)
                    v_bubble maintained at 1.0e-9 m/s
                    Passes "too small" check because recovery_time > 0
```

### Timestep N+Y: After Recovery Period (> 50 μs)
If still collapsing:
```
[RPE_COLLAPSE] Negative Rdot=-1.234e-03, attempting recovery
[RECOVERY #2] R_bubble: 1.191e-07 -> 1.191e-07 m
              v_bubble set to 1.0e-9 m/s again
```

If conditions improved (Rdot > 0):
```
[RPE_RECOVERY_SUCCESS] Bubble recovered! Rdot=1.234e-02 m/s
                       Resetting recovery timer
                       Normal bubble growth resumes
```

---

## Why This Works

1. **Small positive velocity prevents false "too small" trigger**
   - 1.0e-9 m/s passes the `< 1.0e-10` threshold
   - But causes negligible growth: 1.0e-9 m/s × 5e-5 s = 5e-14 m (0.05 nanometers!)

2. **Explicit recovery mode check**
   - Even if velocity is truly small, don't kill parcel if it's in recovery
   - Recovery timer acts as a "protected" flag

3. **Normal operation unaffected**
   - Parcels not in recovery (recovery_time = 0) still get killed if Rdot < 1e-10
   - This prevents parcels from lingering with no growth when they shouldn't

---

## Testing Checklist

- [ ] Compile successfully
- [ ] See `[RPE_COLLAPSE]` messages for collapsing parcels
- [ ] See `[RPE_RECOVERY_WAIT]` messages during 50 μs wait period
- [ ] See multiple `[RECOVERY #2]`, `[RECOVERY #3]`, etc. if still collapsing
- [ ] See `[RPE_RECOVERY_SUCCESS]` when parcel recovers
- [ ] Verify parcel continues bubble growth after recovery
- [ ] Confirm normal parcels (99%) unaffected

---

## Notes

- The 1.0e-9 m/s velocity is 10× larger than the kill threshold (1.0e-10)
- Over 50 μs recovery period: 1e-9 × 5e-5 = 5e-14 m growth (negligible)
- Over 1000 recovery periods: 5e-11 m = 0.05 microns (still negligible)
- Normal bubble growth: 1-20 m/s → 50-1000 microns in 50 μs

---

## Comparison: Before vs After

### Before (Broken):
```
Cycle 82:  Collapse → Recovery → v_bubble = 0
Cycle 83:  Enter RPE → Rdot = 0 → "too small" → KILL (pbt=0, flag=999)
Cycle 84+: Never enters RPE again → Parcel stuck at small size
```

### After (Fixed):
```
Cycle 82:  Collapse → Recovery → v_bubble = 1e-9, recovery_time = t
Cycle 83:  Enter RPE → Rdot = 1e-9 → passes "too small" check
           Check recovery_time → still in wait period → maintain
Cycle 132: Enter RPE → recovery period elapsed → attempt recovery #2
           OR bubble recovered → reset timer → normal growth
```

---

## Conclusion

This was a critical bug that prevented the recovery system from working at all. The parcel was being killed immediately after recovery by an unrelated safety check. The fix ensures that parcels in recovery mode are protected from premature termination while still allowing normal safety checks for other parcels.
