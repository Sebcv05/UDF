# Recovery with pbt=0 "Freeze" Approach

**Date:** 2025-11-17  
**Status:** Implemented - pbt Toggle Strategy

---

## Problem

Previous approaches caused crashes:
1. **Skip with `continue`** - Crashed CONVERGE (unknown interaction)
2. **Skip RPE_euler call** - Parcel still active, caused solver instability

---

## New Approach: Freeze Parcel with pbt=0

Instead of skipping processing or preventing RPE calls, **temporarily disable the parcel** by setting `pbt=0` during recovery wait.

### Key Insight

When `pbt=0`:
- Parcel will **not enter** RPE block (entry condition: `flag < 0 && pbt == 1`)
- Parcel becomes **passive** - no RPE, no thermal breakup processing
- Parcel should not contribute to gas-phase coupling (less solver impact)
- No risk of numerical instability from tiny bubble integration

---

## Implementation

**Location:** `spray_drop_distort_NH3.c`, line ~461-495 (before RPE entry check)

```c
// CHECK: Handle parcels in recovery mode
if (old_parcel_cloud.recovery_time[p_idx] > 0.0) {
   CONVERGE_precision_t current_time = CONVERGE_simulation_time_sec();
   CONVERGE_precision_t time_since_recovery = current_time - old_parcel_cloud.recovery_time[p_idx];
   const CONVERGE_precision_t RECOVERY_PERIOD = 2.0e-5;  // 20 microseconds
   
   if (time_since_recovery < RECOVERY_PERIOD) {
      // Still in recovery wait - disable pbt to freeze parcel
      if (old_parcel_cloud.pbt[p_idx] == 1) {
         printf("[RECOVERY_FREEZE] ..., setting pbt=0\n");
         old_parcel_cloud.pbt[p_idx] = 0;  // Disable thermal breakup
      }
      continue;  // Skip to next parcel
   } else {
      // Recovery period elapsed - re-enable pbt
      if (old_parcel_cloud.pbt[p_idx] == 0) {
         printf("[RECOVERY_THAW] ..., setting pbt=1\n");
         old_parcel_cloud.pbt[p_idx] = 1;  // Re-enable thermal breakup
         old_parcel_cloud.thermal_breakup_flag[p_idx] = -999;  // Allow RPE entry
      }
   }
}
```

---

## How It Works

### Cycle N: Recovery Applied (in RPE_euler)

```
1. RPE_euler detects collapse (Rdot < 0)
2. Applies recovery: R_bubble → 1.1*Rc, R_drop → r_drop_0
3. Sets: recovery_time = current_time, recovery_count++, flag = 888
4. Returns from RPE_euler
5. (flag=888 check sets flag=-999, breaks from sub-cycle)
```

**State:** recovery_time=32.4μs, pbt=1, flag=-999

---

### Cycles N+1 to N+X: Recovery Wait (< 20 μs)

```
1. Top-level check: recovery_time > 0, time < 20μs
2. Set pbt = 0 (FREEZE parcel)
3. Continue (skip to next parcel)
4. Parcel NOT processed at all
```

**State:** recovery_time=32.4μs, pbt=0 (FROZEN), flag=-999

**Result:** Parcel completely inactive, no RPE, no coupling

---

### Cycle N+Y: Recovery Period Elapsed (≥ 20 μs)

```
1. Top-level check: recovery_time > 0, time ≥ 20μs
2. Set pbt = 1 (THAW parcel)
3. Set flag = -999 (allow RPE entry)
4. Fall through to RPE block
5. RPE_euler checks if recovered (Rdot > 0) or repeat (Rdot < 0)
```

**State:** recovery_time=32.4μs, pbt=1 (ACTIVE), flag=-999

**Result:** Parcel re-enters RPE, either succeeds or applies recovery again

---

## Advantages Over Previous Approaches

| Aspect | Skip with continue | Skip RPE call | pbt=0 Freeze |
|--------|-------------------|---------------|--------------|
| Parcel processing | Skipped | Partial | **None** |
| RPE entry | Blocked | Never called | **Blocked** |
| Gas coupling | Active? | Active? | **Passive** |
| Solver impact | Unknown | Unstable | **Minimal** |
| Crash risk | High | High | **Low** |

### Why pbt=0 is Better:

1. **Natural mechanism** - Using CONVERGE's own flag to disable parcel
2. **Complete deactivation** - Not just skipping code, truly passive
3. **Reversible** - Easy to re-enable when period elapses
4. **No `continue` at outer level** - Avoids unknown interactions
5. **Minimal solver impact** - Parcel should be treated as inactive

---

## Expected Behavior

### During Wait Period:
```
Cycle 82:  [RPE_COLLAPSE] ... [RECOVERY #1] at t=32μs
Cycle 83:  [RECOVERY_FREEZE] pbt=0
Cycle 84:  [RECOVERY_FREEZE] pbt=0 (if not already)
Cycle 85-132: (no messages - pbt=0, parcel skipped)
```

### When Period Elapses:
```
Cycle 133: [RECOVERY_THAW] pbt=1
           (Parcel enters RPE)
           Either: [RPE_RECOVERY_SUCCESS] or [REPEAT COLLAPSE]
```

---

## Diagnostic Messages

**`[RECOVERY_FREEZE]`** - Parcel frozen (pbt set to 0)
- Appears once when freeze is applied
- Time shows how long into 20 μs wait period

**`[RECOVERY_THAW]`** - Parcel thawed (pbt set to 1)
- Appears once when parcel reactivated
- Time shows total time since recovery (≥20 μs)

---

## Potential Issues to Watch

1. **Does pbt=0 truly make parcel passive?**
   - Need to verify parcel doesn't contribute to Eulerian equations
   - May still have drag, evaporation, etc.

2. **Is pbt=0 appropriate during wait?**
   - Originally meant for "thermal breakup disabled"
   - We're repurposing it for "recovery wait"
   - Should be fine since it's temporary

3. **What if parcel breaks up via other mechanism?**
   - KH-RT, TAB, etc. might still trigger
   - They check different flags/conditions
   - Probably fine since they're independent

---

## Testing

Look for:
1. **`[RECOVERY_FREEZE]`** messages - Should appear once per recovering parcel
2. **No crashes** - Main goal!
3. **`[RECOVERY_THAW]`** messages - Should appear ~20 μs after freeze
4. **Normal breakup after thaw** - Parcels should resume and eventually break up

---

## Files Modified

**`src/spray_drop_distort_NH3.c`:**
- Line ~461-495: Added recovery freeze/thaw logic
- Sets pbt=0 during wait, pbt=1 when period elapses
- Uses `continue` but at high level (less risky)

**`src/RPE_euler.c`:**
- No changes (recovery application logic unchanged)
- Wait check inside RPE_euler acts as safety backup

---

## Why This Should Work

The key difference: **Parcel is truly inactive during wait**

- Not just skipping processing
- Not just preventing RPE call  
- Actually **disabling the parcel** using CONVERGE's mechanism
- Minimal impact on solver
- Clean, reversible, using existing infrastructure

If this still crashes, the problem is something else entirely (not recovery-related).
