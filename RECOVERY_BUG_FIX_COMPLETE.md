# RECOVERY STRUCTURE - COMPLETE ANALYSIS AND BUG FIX

**Date:** 2025-11-17  
**Status:** BUG IDENTIFIED AND FIXED

---

## THE BUG

**Symptom:** Parcels apply recovery but never resume RPE after the 20 μs wait period.

**Root Cause:** The flag 888 → -999 conversion was being skipped, preventing re-entry.

### What Was Happening:

```
Cycle 82 (t=32.4 μs): Recovery applied
  ├─ RPE_euler: Rdot < 0 → set flag=888, recovery_time=32.4 μs
  ├─ Back in spray_drop_distort: My check breaks due to time < 20 μs
  └─ Flag 888 → -999 conversion (line 563) SKIPPED!

Cycle 83 (t=32.5 μs): Try to enter RPE
  ├─ Entry check: flag < 0 && pbt == 1
  ├─ Flag = 888 (NOT < 0!)
  └─ BLOCKED from entering RPE!

Cycles 84-412: Forever blocked
  └─ Flag stays at 888, never enters RPE again
```

### The Fix:

Added flag reset BEFORE breaking:
```c
if (time_since_recovery < RECOVERY_PERIOD) {
   // CRITICAL: Reset flag to -999 so parcel can re-enter RPE next cycle
   if (old_parcel_cloud.thermal_breakup_flag[p_idx] == 888) {
      old_parcel_cloud.thermal_breakup_flag[p_idx] = -999;
   }
   break;
}
```

---

## COMPLETE RECOVERY STRUCTURE

### 1. Entry Condition (spray_drop_distort_NH3.c, line 464)

```c
if (thermal_breakup_flag < 0 && pbt == 1)
```

**Parcels can enter if:**
- `flag = -1` (normal bubble growth)
- `flag = -999` (post-recovery, waiting or checking)

**Parcels CANNOT enter if:**
- `flag = 0` (not initialized)
- `flag = 888` (just recovered - TEMPORARY)
- `flag = 999` (aborted)
- `flag = 3` (broken up)
- `flag >=  5` (various breakup conditions)
- `pbt = 0` (RPE disabled)

---

### 2. RPE_euler Processing

#### A. Negative Rdot Detected (Collapse)

**Location:** RPE_euler.c, line ~349-432

```c
if (state.Rdot < 0) {
    // Check if in recovery wait period
    if (recovery_time > 0 && time_since_recovery < RECOVERY_PERIOD) {
        // STILL WAITING - don't apply new recovery
        v_bubble = 1.0e-9;  // Keep stable
        flag = 888;         // Signal wait
        return;
    }
    
    // Check if recovery period just elapsed
    if (recovery_time > 0 && time_since_recovery >= RECOVERY_PERIOD) {
        // REPEAT COLLAPSE - apply recovery again
        printf("[REPEAT COLLAPSE] ...");
        // Apply recovery (R_bubble → 1.1*Rc, etc.)
        recovery_time = current_time;  // Reset timer
        recovery_count++;
        flag = 888;
        return;
    }
    
    // FIRST COLLAPSE - apply recovery
    printf("[RPE_COLLAPSE] ...");
    printf("[RECOVERY #1] ...");
    // Apply recovery
    recovery_time = current_time;  // Start timer
    recovery_count = 1;
    flag = 888;
    return;
}
```

#### B. Positive Rdot (Growing)

**Location:** RPE_euler.c, line ~443-458

```c
// Check if this is recovery success
if (recovery_time > 0 && state.Rdot > 0 && time_since_recovery >= RECOVERY_PERIOD) {
    printf("[RPE_RECOVERY_SUCCESS] Rdot=%.3e, %.3e s after recovery\n", ...);
    recovery_time = 0.0;  // Reset timer
    flag = -1;            // Normal operation
}
// Continue with normal bubble growth
```

---

### 3. Post-RPE Checks (spray_drop_distort_NH3.c)

#### A. Recovery Period Check (line ~524-548) - **MY NEW CODE**

```c
if (recovery_time > 0) {
    if (time_since_recovery < RECOVERY_PERIOD) {
        // STILL WAITING
        printf("[RECOVERY_WAIT_SUBCYCLE] ...");
        if (flag == 888) flag = -999;  // CRITICAL: Allow re-entry next cycle
        break;
    } else {
        // PERIOD ELAPSED
        printf("[RECOVERY_PERIOD_ELAPSED] ...");
        // Continue (RPE_euler handles success/repeat)
    }
}
```

#### B. Small Velocity Check (line ~550-556)

```c
if (v_bubble < 1.0e-10) {
    pbt = 0;
    flag = 999;
    break;
}
```

#### C. Flag 888 Check (line ~558-565) - **REDUNDANT NOW**

```c
if (flag == 888) {
    flag = -999;
    break;
}
```

**Note:** This is now redundant because my check above already handles it.

---

## RECOVERY LIFECYCLE - CORRECT FLOW

### Cycle N: First Collapse

```
1. Parcel enters RPE (flag=-1, pbt=1)
2. RPE_euler: Rdot < 0, recovery_time=0 → FIRST COLLAPSE
3. Apply recovery: R_bubble → 1.1*Rc, R_drop → r_drop_0
4. Set: recovery_time = current_time, recovery_count = 1, flag = 888
5. Return from RPE_euler
6. My check: recovery_time > 0, time=0 < 20μs
   - Set flag = -999 (CRITICAL!)
   - Break
7. End of cycle
```

**State:** flag=-999, recovery_time=32.4μs, recovery_count=1, pbt=1

---

### Cycles N+1 to N+X: Recovery Wait (< 20 μs)

```
1. Parcel enters RPE (flag=-999 < 0, pbt=1) ✓
2. RPE_euler: Rdot < 0, recovery_time > 0, time < 20μs → STILL WAITING
3. Set: v_bubble = 1.0e-9, flag = 888
4. Return from RPE_euler
5. My check: recovery_time > 0, time < 20μs
   - Set flag = -999 (CRITICAL!)
   - Break
6. End of cycle
```

**State:** flag=-999, recovery_time=32.4μs (unchanged), recovery_count=1, pbt=1

**This repeats every cycle for 20 μs**

---

### Cycle N+Y: Recovery Period Elapsed (≥ 20 μs)

#### Scenario A: Bubble Recovered (Rdot > 0)

```
1. Parcel enters RPE (flag=-999 < 0, pbt=1) ✓
2. RPE_euler: Rdot > 0, recovery_time > 0, time ≥ 20μs → SUCCESS!
3. Set: recovery_time = 0, flag = -1
4. Continue with normal bubble growth
5. Return from RPE_euler
6. My check: recovery_time = 0 → Skip (not in recovery)
7. Continue sub-cycling normally
```

**State:** flag=-1, recovery_time=0, recovery_count=1, pbt=1

**Recovery complete, normal operation resumes**

#### Scenario B: Still Collapsing (Rdot < 0)

```
1. Parcel enters RPE (flag=-999 < 0, pbt=1) ✓
2. RPE_euler: Rdot < 0, recovery_time > 0, time ≥ 20μs → REPEAT COLLAPSE
3. Apply recovery again: R_bubble → 1.1*Rc, R_drop → r_drop_0
4. Set: recovery_time = current_time (NEW), recovery_count = 2, flag = 888
5. Return from RPE_euler
6. My check: recovery_time > 0, time=0 < 20μs (reset!)
   - Set flag = -999
   - Break
7. End of cycle
```

**State:** flag=-999, recovery_time=52.4μs (NEW), recovery_count=2, pbt=1

**Another 20 μs wait begins**

---

## KEY INSIGHTS

### 1. Flag Management is Critical

| Flag | Meaning | Can Enter RPE? |
|------|---------|----------------|
| -1   | Normal growth | ✓ Yes (flag < 0) |
| -999 | Post-recovery wait | ✓ Yes (flag < 0) |
| 0    | Not initialized | ✗ No |
| 888  | Just recovered (temp) | ✗ No (flag > 0!) |
| 999  | Aborted | ✗ No |
| 3+   | Breakup states | ✗ No |

**The 888 flag is temporary and MUST be converted to -999 before the next cycle!**

### 2. Two Places Set Flag to 888

1. **RPE_euler.c** (line ~373, ~399, ~431)
   - First collapse
   - Recovery wait
   - Repeat collapse

2. Both need flag → -999 conversion before next cycle

### 3. My Fix Ensures Conversion

Before my fix:
- flag=888 → -999 only happened at line 563
- My break at line 538 skipped this
- Flag stayed at 888 → blocked from re-entry

After my fix:
- Explicitly convert flag=888 → -999 before breaking
- Parcel can re-enter next cycle

---

## DIAGNOSTIC MESSAGES

### Normal Flow:

```
[RPE_COLLAPSE] ... [RECOVERY #1] at t=32μs
[RECOVERY_WAIT_SUBCYCLE] 0μs/20μs
[RECOVERY_WAIT_SUBCYCLE] 1μs/20μs
[RECOVERY_WAIT_SUBCYCLE] 2μs/20μs
...
[RECOVERY_WAIT_SUBCYCLE] 19μs/20μs
[RECOVERY_PERIOD_ELAPSED] 21μs elapsed
[RPE_RECOVERY_SUCCESS] Rdot > 0
```

### Multiple Attempts:

```
[RPE_COLLAPSE] ... [RECOVERY #1] at t=32μs
[RECOVERY_WAIT_SUBCYCLE] ... (20μs)
[RECOVERY_PERIOD_ELAPSED] 21μs
[REPEAT COLLAPSE] ... [RECOVERY #2]
[RECOVERY_WAIT_SUBCYCLE] ... (20μs)
[RECOVERY_PERIOD_ELAPSED] 21μs
[RPE_RECOVERY_SUCCESS] Rdot > 0
```

---

## FILES MODIFIED

**`src/spray_drop_distort_NH3.c`:**
- Line ~524-548: Added recovery period check after RPE_euler
- Line ~539: CRITICAL FIX - Reset flag 888 → -999 before breaking

**`src/RPE_euler.c`:**
- Line ~364-375: Recovery wait logic (sets flag=888)
- Line ~377-432: Recovery application logic
- Line ~443-458: Recovery success check

---

## WHY IT FAILED BEFORE

The parcel was stuck in an infinite loop:
1. Apply recovery → flag=888
2. Break before flag→-999 conversion
3. Try to enter next cycle: flag=888 > 0 → BLOCKED!
4. Never enters RPE again
5. Never checks if period elapsed
6. Never attempts recovery success

The fix: **Always convert flag 888 → -999 before breaking**, ensuring the parcel can re-enter on the next cycle.

---

## TESTING

After this fix, you should see:
1. ✓ `[RECOVERY_WAIT_SUBCYCLE]` messages for ~20 μs
2. ✓ `[RECOVERY_PERIOD_ELAPSED]` message when period elapses
3. ✓ `[RPE_RECOVERY_SUCCESS]` or `[REPEAT COLLAPSE]` messages
4. ✓ Parcels breaking up eventually (not stuck)

If you don't see these, there's another issue.
