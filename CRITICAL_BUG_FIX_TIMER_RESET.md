# CRITICAL BUG FIX: Recovery Timer Premature Reset

**Date:** 2025-11-17  
**Status:** Fixed

---

## The Problem

Looking at the log (line 473+):

```
Cycle 82, Time=3.239393e-05 s: [RPE_COLLAPSE] ... [RECOVERY #1] logged at 3.239393e-05 s
Cycle 83, Time=3.239393e-05 s: [RPE_RECOVERY_SUCCESS] after 6.143e-08 s (61 nanoseconds!)
```

**The recovery timer was being reset immediately** - only 61 nanoseconds after recovery, instead of waiting 20 microseconds!

---

## Root Cause

The recovery success check was:
```c
if (recovery_time > 0 && Rdot > 0) {
    // Reset timer immediately!
}
```

This check ran **every time RPE was called**, including:
- Within the same global timestep
- Within sub-cycling loops  
- Before the 20 μs recovery period could elapse

### What Was Happening:

1. **Cycle 82, Sub-step 500:** Collapse detected → Recovery applied → `recovery_time = 3.239393e-05`
2. **Cycle 82, Sub-step 501:** (61 ns later) RPE called again → Rdot is now positive → Timer RESET!
3. **Result:** Recovery period never enforced, timer reset instantly

---

## The Fix

Added a check to ensure the recovery period has elapsed:

```c
// RESET recovery timer if bubble is growing successfully after recovery
// BUT only if the recovery period has elapsed
if (old_parcel_cloud->recovery_time[p_idx] > 0.0 && state.Rdot > 0.0) {
    CONVERGE_precision_t current_time = CONVERGE_simulation_time_sec();
    CONVERGE_precision_t recovery_start_time = old_parcel_cloud->recovery_time[p_idx];
    CONVERGE_precision_t time_since_recovery = current_time - recovery_start_time;
    
    const CONVERGE_precision_t RECOVERY_PERIOD = 2.0e-5;  // 20 microseconds
    
    // Only declare success if recovery period has elapsed
    if (time_since_recovery >= RECOVERY_PERIOD) {
        printf("[RPE_RECOVERY_SUCCESS] ...");
        old_parcel_cloud->recovery_time[p_idx] = 0.0;  // NOW it's safe to reset
        old_parcel_cloud->thermal_breakup_flag[p_idx] = -1;
    }
}
```

---

## Expected Behavior After Fix

### Scenario 1: Recovery Succeeds Quickly (Conditions Improve)

```
Time=3.239393e-05 s: [RPE_COLLAPSE] ... [RECOVERY #1]
Time=3.239393e-05 s to 5.239393e-05 s: Rdot positive but timer NOT reset (< 20 μs)
Time=5.239393e-05 s: [RPE_RECOVERY_SUCCESS] 20 μs elapsed, Rdot still positive → Reset timer
```

### Scenario 2: Parcel Collapses Again (Conditions Still Bad)

```
Time=3.239393e-05 s: [RPE_COLLAPSE] ... [RECOVERY #1]
Time=3.239393e-05 s to 5.239393e-05 s: [RPE_RECOVERY_WAIT] ... waiting
Time=5.239393e-05 s: Rdot still negative → [RPE_COLLAPSE] ... [RECOVERY #2]
Time=5.239393e-05 s to 7.239393e-05 s: [RPE_RECOVERY_WAIT] ... waiting
Time=7.239393e-05 s: [RPE_COLLAPSE] ... [RECOVERY #3]
... continues until conditions improve or max attempts
```

### Scenario 3: Partial Recovery (Oscillating Conditions)

```
Time=3.239393e-05 s: [RPE_COLLAPSE] ... [RECOVERY #1]
Time=3.5e-05 s to 5.0e-05 s: Rdot positive (recovery working) but timer NOT reset
Time=5.1e-05 s: Rdot goes negative again → [RPE_COLLAPSE] ... [RECOVERY #2]
Time=7.2e-05 s: Rdot positive for > 20 μs → [RPE_RECOVERY_SUCCESS]
```

---

## Why This Matters

Without this fix:
- Recovery period was NEVER enforced
- Timer reset after 61 nanoseconds instead of 20 microseconds
- Parcels that needed multiple recovery attempts only got one
- System couldn't handle parcels in difficult flow conditions

With this fix:
- Recovery period properly enforced (20 μs)
- Parcels get multiple attempts if needed
- Success only declared after sustained positive growth
- Better handling of pressure waves and adverse conditions

---

## Files Modified

**`src/RPE_euler.c`:**
- Line ~435-450: Added recovery period check before resetting timer
- Now requires `time_since_recovery >= RECOVERY_PERIOD` for success

---

## Testing

The next run should show:

1. **Recovery attempts spaced 20 μs apart:**
   ```
   [RECOVERY #1] at time T
   [RECOVERY #2] at time T + 2e-5 s
   [RECOVERY #3] at time T + 4e-5 s
   ```

2. **Recovery success only after 20 μs of positive growth:**
   ```
   [RPE_RECOVERY_SUCCESS] ... 2.123e-05 s after recovery (> 20 μs)
   ```

3. **No more premature resets:**
   - Should NOT see "6.143e-08 s after recovery"
   - All success messages should be > 2e-5 s

---

## Impact

This was a **critical bug** that completely disabled the recovery wait period. The system was working as a "single attempt recovery" instead of "time-based recovery with wait period."

With this fix, parcels experiencing prolonged adverse conditions will get multiple recovery attempts spaced 20 μs apart, giving the flow time to improve between attempts.

---

## Related Issues

This bug is why you were seeing:
- Only one recovery message per parcel
- No "REPEAT COLLAPSE" messages
- Recovery appearing to work but parcels still being rogue

The recovery WAS working (resetting the bubble), but the timer was being reset immediately, so the system thought the parcel had fully recovered when it might still be in adverse flow conditions.
