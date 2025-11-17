# Recovery Period Check After RPE_euler

**Date:** 2025-11-17  
**Status:** Implemented - Alternative Approach

---

## Problem

After recovery is applied, parcels continue through the sub-cycle loop but shouldn't actually do anything until the 20 μs recovery period elapses. The previous approach (using `continue` at the top level) caused crashes.

---

## Solution

Add a check **immediately after RPE_euler_solver** returns to break out of the sub-cycle loop if the recovery period hasn't elapsed yet.

### Location

**`spray_drop_distort_NH3.c`** - Line ~523 (right after RPE_euler_solver call)

### Code Added

```c
// CHECK: If parcel is in recovery mode and period hasn't elapsed, break from sub-cycle
if (old_parcel_cloud.recovery_time[p_idx] > 0.0) {
   CONVERGE_precision_t current_time = CONVERGE_simulation_time_sec();
   CONVERGE_precision_t time_since_recovery = current_time - old_parcel_cloud.recovery_time[p_idx];
   const CONVERGE_precision_t RECOVERY_PERIOD = 2.0e-5;  // 20 microseconds
   
   if (time_since_recovery < RECOVERY_PERIOD) {
      // Still in recovery period - skip rest of sub-cycle for this parcel
      printf("[RECOVERY_WAIT_SUBCYCLE] p_idx=%li, time_since_recovery=%.3e/%.3e s, breaking sub-cycle\n",
             p_idx, time_since_recovery, RECOVERY_PERIOD);
      break;  // Exit sub-cycle loop, move to next parcel
   } else {
      // Recovery period has elapsed - check if recovered or needs retry in RPE_euler
      printf("[RECOVERY_PERIOD_ELAPSED] p_idx=%li, time_since_recovery=%.3e s, continuing sub-cycle\n",
             p_idx, time_since_recovery);
   }
}
```

---

## How It Works

### Cycle When Recovery Is Applied (flag=888)

1. RPE_euler detects collapse (Rdot < 0)
2. Applies recovery (R_bubble → 1.1*Rc, R_drop → r_drop_0)
3. Sets `recovery_time = current_time`
4. Sets `thermal_breakup_flag = 888`
5. Returns from RPE_euler
6. **Our new check:** Sees recovery_time > 0 and time < 20 μs → break
7. Lower check sees flag=888 → sets flag=-999, break (redundant now)
8. Sub-cycle ends, moves to next parcel

### Subsequent Cycles (0 < time < 20 μs)

1. Parcel enters RPE block (flag=-999, pbt=1)
2. RPE_euler is called
3. RPE_euler sees `recovery_time > 0` and checks wait period
4. If still waiting, RPE_euler sets flag=888 and returns early
5. **Our new check:** Sees recovery_time > 0 and time < 20 μs → break
6. Sub-cycle ends immediately, no processing

**Result:** Parcel does minimal work during wait period

### After 20 μs Elapses

1. Parcel enters RPE block (flag=-999, pbt=1)
2. RPE_euler is called
3. RPE_euler sees recovery period elapsed, checks Rdot:
   - **If Rdot > 0:** Success! Resets recovery_time=0, flag=-1, continues
   - **If Rdot < 0:** Still collapsing, applies recovery again (recovery #2)
4. **Our new check:**
   - If Rdot > 0: recovery_time now 0, check passes, continues normally
   - If Rdot < 0: recovery_time reset to current_time, time < 20 μs, break

---

## Why This Works Better Than Previous Approach

### Previous Approach (Caused Crash)
- Check at **top level** before entering RPE block
- Used `continue` to skip entire parcel processing
- May have interfered with CONVERGE bookkeeping

### New Approach (Safer)
- Check **inside sub-cycle loop** after RPE_euler returns
- Uses `break` to exit sub-cycle (already used elsewhere in code)
- Parcel still goes through RPE_euler (which handles recovery logic)
- Minimal changes to control flow
- No interference with outer loop structure

---

## Expected Log Output

### During Recovery Wait:
```
Cycle 82:  [RPE_COLLAPSE] ... [RECOVERY #1]
Cycle 83:  [RECOVERY_WAIT_SUBCYCLE] 1μs/20μs, breaking sub-cycle
Cycle 84:  [RECOVERY_WAIT_SUBCYCLE] 2μs/20μs, breaking sub-cycle
...
Cycle 132: [RECOVERY_WAIT_SUBCYCLE] 19μs/20μs, breaking sub-cycle
```

### When Period Elapses:
```
Cycle 133: [RECOVERY_PERIOD_ELAPSED] 21μs elapsed, continuing sub-cycle
```

**Followed by one of:**
- `[RPE_RECOVERY_SUCCESS]` - If Rdot > 0 (recovered)
- `[REPEAT COLLAPSE] ... [RECOVERY #2]` - If Rdot < 0 (needs another attempt)

---

## Advantages

1. **Safer:** Uses `break` in sub-cycle loop (already used elsewhere)
2. **Cleaner:** Check is right where it matters (after RPE_euler)
3. **Self-contained:** Doesn't affect outer parcel loop structure
4. **Compatible:** Works with existing flag=888 logic
5. **Minimal:** Small addition, no major restructuring

---

## Testing

Compile and run. Look for:

1. **`[RECOVERY_WAIT_SUBCYCLE]` messages** - Should appear every cycle during wait
2. **`[RECOVERY_PERIOD_ELAPSED]` message** - Should appear once when 20 μs passes
3. **`[RPE_RECOVERY_SUCCESS]` or `[REPEAT COLLAPSE]`** - Shows outcome after period elapses
4. **No crashes** - This approach doesn't change outer loop control flow

---

## Files Modified

**`src/spray_drop_distort_NH3.c`:**
- Line ~523-551: Added recovery period check after RPE_euler_solver call
- Uses `break` to exit sub-cycle loop if period hasn't elapsed
- Logs diagnostic messages for wait and elapsed states

---

## Comparison to Failed Approach

| Aspect | Failed Approach | New Approach |
|--------|----------------|--------------|
| Location | Before RPE block | After RPE_euler call |
| Control | `continue` (skip parcel) | `break` (exit sub-cycle) |
| RPE Entry | Never enters | Enters but exits quickly |
| Safety | Changed outer loop | Uses existing pattern |
| Result | Crashed | Should work |

The key difference is that the parcel still goes through RPE_euler (which has its own wait logic), but we prevent any further processing in the sub-cycle loop until the period elapses.
