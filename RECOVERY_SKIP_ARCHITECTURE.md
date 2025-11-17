# Recovery Architecture: Skip Sub-Cycling During Wait Period

**Date:** 2025-11-17  
**Status:** Implemented - Major Restructure

---

## Problem with Previous Approach

The old approach allowed parcels in recovery to continue through the sub-cycling loop:
1. Recovery applied → bubble set to 1.1*Rc → v_bubble set to small value
2. Parcel continues sub-cycling with tiny bubble
3. Safety checks could kill it (negative P_sat-P_amb, subcooled, small velocity, etc.)
4. Needed complex "protection" logic to skip all safety checks
5. Still risky - any new safety check could break recovery

**Result:** Fragile system prone to breaking with code changes

---

## New Architecture

### High-Level Flow

```
For each parcel:
  ├─ Check if in recovery mode (recovery_time > 0)
  │  ├─ If time_since_recovery < 20 μs:
  │  │  └─ Skip this parcel entirely (continue to next)
  │  └─ If time_since_recovery >= 20 μs:
  │     └─ Allow parcel to proceed to RPE
  │
  └─ Normal processing (RPE, breakup checks, etc.)
```

### Implementation Location

**`spray_drop_distort_NH3.c`** - Line ~461-488:

```c
// CHECK: Is this parcel in recovery mode waiting for the recovery period to elapse?
if (old_parcel_cloud.recovery_time[p_idx] > 0.0 && old_parcel_cloud.pbt[p_idx] == 1) {
   CONVERGE_precision_t current_time = CONVERGE_simulation_time_sec();
   CONVERGE_precision_t time_since_recovery = current_time - old_parcel_cloud.recovery_time[p_idx];
   const CONVERGE_precision_t RECOVERY_PERIOD = 2.0e-5;  // 20 microseconds
   
   if (time_since_recovery < RECOVERY_PERIOD) {
      // Still in recovery period - skip this parcel entirely
      printf("[RECOVERY_SKIP] ...");
      continue;  // Skip to next parcel
   } else {
      // Recovery period has elapsed - allow parcel to proceed
      printf("[RECOVERY_ELAPSED] ...");
      // Don't reset recovery_time yet - let RPE_euler handle success/failure
   }
}

// Now proceed with RPE block (only if recovery period elapsed or not in recovery)
if (thermal_breakup_flag < 0 && pbt == 1) {
   // RPE processing...
}
```

---

## Key Advantages

### 1. **No Safety Check Conflicts**
- Parcel never enters RPE during wait period
- No need to protect against individual safety checks
- No risk of new safety checks breaking recovery

### 2. **Cleaner Code**
- Single check at top level
- No complex protection logic scattered throughout RPE_euler
- Easy to understand: "waiting" vs "active"

### 3. **Physically Correct**
- Parcel is "dormant" during recovery wait
- No sub-cycling with unstable tiny bubble
- Bubble size, velocity, etc. frozen until recovery period elapses

### 4. **Robust to Code Changes**
- Adding new safety checks won't break recovery
- Changing RPE logic won't affect recovery
- Self-contained at high level

---

## Recovery Lifecycle

### Timestep N: Collapse Detected

**In RPE_euler.c:**
```c
if (Rdot < 0) {
   // Apply recovery: R_bubble → 1.1*Rc, R_drop → r_drop_0
   recovery_time = current_time;  // Log time
   recovery_count++;
   thermal_breakup_flag = 888;   // Signal recovery
   return;                       // Exit RPE
}
```

**In spray_drop_distort.c:**
```c
if (flag == 888) {
   flag = -999;  // Reset for next timestep
   break;        // Exit sub-loop
}
```

**Result:** Parcel exits immediately, no further processing this timestep

---

### Timesteps N+1 to N+X: Recovery Wait (< 20 μs)

**In spray_drop_distort.c (BEFORE RPE block):**
```c
if (recovery_time > 0 && time_since_recovery < 20 μs) {
   printf("[RECOVERY_SKIP] ...");
   continue;  // Skip to next parcel
}
```

**Result:** Parcel not processed at all, completely bypassed

---

### Timestep N+Y: Recovery Period Elapsed (>= 20 μs)

**In spray_drop_distort.c:**
```c
if (recovery_time > 0 && time_since_recovery >= 20 μs) {
   printf("[RECOVERY_ELAPSED] ...");
   // Allow to proceed (don't continue)
}
```

**Parcel enters RPE block, one of two outcomes:**

#### Outcome A: Still Collapsing (Rdot < 0)
```c
// In RPE_euler.c
if (Rdot < 0 && time_since_recovery >= RECOVERY_PERIOD) {
   printf("[REPEAT COLLAPSE] attempting recovery #2");
   // Apply recovery again
   recovery_time = current_time;  // Update time
   recovery_count++;              // Increment counter
   return;
}
```
**Result:** Another 20 μs wait begins

#### Outcome B: Recovered Successfully (Rdot > 0)
```c
// In RPE_euler.c
if (recovery_time > 0 && Rdot > 0 && time_since_recovery >= RECOVERY_PERIOD) {
   printf("[RPE_RECOVERY_SUCCESS] ...");
   recovery_time = 0.0;  // Reset timer
   flag = -1;            // Normal operation
}
// Continue with normal bubble growth
```
**Result:** Recovery complete, normal RPE resumes

---

## Diagnostic Messages

### During Wait Period:
```
[RECOVERY_SKIP] p_idx=X, time_since_recovery=5.234e-06/2.000e-05 s, skipping to next parcel
```
- Appears every timestep during 20 μs wait
- Shows parcel is being bypassed

### When Period Elapses:
```
[RECOVERY_ELAPSED] p_idx=X, time_since_recovery=2.123e-05 s, allowing RPE entry
```
- Appears once when parcel allowed back into RPE
- Followed by either REPEAT COLLAPSE or RECOVERY_SUCCESS

---

## Comparison: Old vs New

### Old Approach (Fragile):
```
Collapse → Recovery applied
Next timestep:
  ├─ Enter RPE with tiny bubble
  ├─ Check 1: Subcooled? → Protected (complex logic)
  ├─ Check 2: Small velocity? → Protected (complex logic)  
  ├─ Check 3: Negative P? → Protected (complex logic)
  ├─ Check 4: ???  → Could break recovery!
  └─ Update state with tiny bubble
```

### New Approach (Robust):
```
Collapse → Recovery applied
Next 20 μs:
  └─ Skip parcel entirely (no RPE entry, no checks, nothing)

After 20 μs:
  └─ Enter RPE normally, check if recovered or needs another attempt
```

---

## What This Fixes

1. ✅ **Parcels can't be killed during recovery** - They're not processed at all
2. ✅ **No protection logic needed** - Skip happens before any checks
3. ✅ **Clean separation** - "Waiting" vs "Active" states
4. ✅ **Future-proof** - Adding new checks won't break recovery
5. ✅ **Physically correct** - Dormant period while conditions stabilize

---

## Files Modified

**`src/spray_drop_distort_NH3.c`:**
- Line ~461-488: Added recovery skip check BEFORE RPE block
- Uses `continue` to skip to next parcel if in recovery wait

**`src/RPE_euler.c`:**
- Recovery wait check (line ~364-375) now acts as redundant safety check
- Main recovery logic unchanged

---

## Testing

Look for these patterns in the log:

### Successful Recovery:
```
Cycle 82:  [RPE_COLLAPSE] ... [RECOVERY #1] at t=32 μs
Cycle 83:  [RECOVERY_SKIP] 1μs/20μs elapsed
Cycle 84:  [RECOVERY_SKIP] 2μs/20μs elapsed
...
Cycle 132: [RECOVERY_SKIP] 19μs/20μs elapsed
Cycle 133: [RECOVERY_ELAPSED] 21μs elapsed
           [RPE_RECOVERY_SUCCESS] Rdot > 0
```

### Multiple Recovery Attempts:
```
Cycle 82:  [RPE_COLLAPSE] ... [RECOVERY #1]
Cycle 83-132: [RECOVERY_SKIP] ...
Cycle 133: [RECOVERY_ELAPSED] 21μs
           [REPEAT COLLAPSE] ... [RECOVERY #2]
Cycle 134-183: [RECOVERY_SKIP] ...
Cycle 184: [RECOVERY_ELAPSED] 21μs
           [RPE_RECOVERY_SUCCESS]
```

### What Should NOT Appear:
- ❌ `[RPE_KILL_IN_RECOVERY]` - Parcels can't be killed if they're skipped
- ❌ `[RPE_IN_RECOVERY]` during wait - Parcel shouldn't enter RPE
- ❌ Any safety check triggers for parcels in recovery wait

---

## Benefits

This is a **much better architecture** because:
1. Simple and robust
2. Self-documenting code
3. Physically meaningful (dormant period)
4. Impossible to break with future changes
5. No complex interdependencies

The parcel literally does nothing during the recovery wait period - it's completely bypassed in the main loop. Only after 20 μs does it re-enter the normal processing flow.
