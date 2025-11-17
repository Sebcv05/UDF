# Code Reverted - Recovery Skip Architecture Removed

**Date:** 2025-11-17  
**Status:** Reverted due to crash

---

## What Was Attempted

Implemented a "skip parcel during recovery" architecture where parcels in recovery mode would be completely bypassed in the main loop using `continue` statement.

## Why It Was Reverted

The simulation crashed with:
- SIGNAL 11 (Segmentation fault) on rank 42
- SIGNAL 9 (Killed) on all other ranks
- Error message: "recovering .... because transport equations did not converge"

While the root cause appears to be a CONVERGE solver convergence failure (not directly our code), the crash did not occur before this change, so it was reverted as a precaution.

---

## Changes Reverted in spray_drop_distort_NH3.c

### Removed: Recovery Skip Check (lines ~461-488)
```c
// This was REMOVED:
if (old_parcel_cloud.recovery_time[p_idx] > 0.0 && old_parcel_cloud.pbt[p_idx] == 1) {
   CONVERGE_precision_t current_time = CONVERGE_simulation_time_sec();
   CONVERGE_precision_t time_since_recovery = current_time - old_parcel_cloud.recovery_time[p_idx];
   const CONVERGE_precision_t RECOVERY_PERIOD = 2.0e-5;
   
   if (time_since_recovery < RECOVERY_PERIOD) {
      printf("[RECOVERY_SKIP] ...");
      continue;  // This may have caused issues
   } else {
      printf("[RECOVERY_ELAPSED] ...");
   }
}
```

### Removed: RPE Entry Logging (lines ~464-477)
```c
// This was REMOVED:
if (old_parcel_cloud.recovery_time[p_idx] > 0.0) {
    printf("[RPE_ENTRY_IN_RECOVERY] ...");
}
```

### Removed: Recovery Break Logging (lines ~537-542)
```c
// This was REMOVED:
printf("[RECOVERY_BREAK] p_idx=%li, flag 888->-999, ...");
```

---

## Changes KEPT in RPE_euler.c

The following diagnostic additions were KEPT because they're just printf statements and shouldn't cause crashes:

1. `[RPE_KILL_IN_RECOVERY]` - Logs when parcels in recovery get killed by safety checks
2. `[RPE_SUBCOOL_IN_RECOVERY]` - Logs when subcooled check is skipped
3. `[RPE_SMALL_VEL_IN_RECOVERY]` - Logs when small velocity check is skipped

These are useful for understanding what happens to parcels during recovery and don't modify control flow.

---

## Current State

The code is back to the state before the "skip architecture" was implemented:
- Parcels in recovery continue through RPE with protection checks
- Recovery wait logic in RPE_euler.c sets flag=888 and returns early
- Diagnostic logging remains to help debug recovery issues

---

## Potential Cause of Crash (Speculation)

The `continue` statement may have interacted poorly with:
1. Other control flow logic in the parcel loop
2. Memory management or state updates expected for every parcel
3. Some CONVERGE internal bookkeeping

Without deeper investigation, the exact cause is unclear, but the crash correlation with the change necessitated the revert.

---

## Next Steps

To fix the recovery issue without the skip architecture:
1. Keep the diagnostic logging to understand what's happening
2. Run simulation to see what kills recovered parcels
3. Add targeted protections based on diagnostic output
4. Consider alternative approaches that don't use `continue`

---

## Lessons Learned

- The `continue` approach, while elegant, may have unintended consequences in complex loops
- CONVERGE's internal state management may expect certain operations for all parcels
- Changes that modify control flow need more careful testing
- Reverting quickly when issues arise is the right approach

---

## Files Modified (Revert)

**`src/spray_drop_distort_NH3.c`:**
- Removed recovery skip check with `continue`
- Removed RPE entry logging
- Removed recovery break logging
- Back to original structure

**`src/RPE_euler.c`:**
- NO CHANGES (diagnostic logging kept)

The code should now run without crashing.
