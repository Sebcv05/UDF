# Logging Cleanup and Recovery Tracking

**Date:** 2025-11-17  
**Status:** Implemented

---

## Changes Made

### 1. Disabled Verbose Breakup Logging

**Purpose:** Clean up log file to focus on recovery behavior

**Files Modified:**
- `src/spray_drop_distort_NH3.c` - Commented out detailed `[BREAKUP]` message
- `src/Breakup.c` - Commented out:
  - `BREAKUP_COMPLETE` messages
  - `[RR_SAMPLE]` messages (5 print statements)

These messages were cluttering the log and making it difficult to trace recovery behavior.

---

### 2. Added Recovery Status Tracking

**Purpose:** Show when parcels are in recovery mode and when RPE resumes

**Added to `src/RPE_euler.c`:**

#### A. Recovery Mode Status (at start of RPE call)
```c
[RPE_IN_RECOVERY] p_idx=X, recovery_count=Y, time_since_recovery=Z.ZZZe-XX s, R=X.XXXe-XX m
```
- Prints up to 20 times to avoid spam
- Shows which parcels are currently in recovery
- Displays time elapsed since last recovery
- Shows current bubble radius

#### B. Recovery Period Elapsed Message
```c
[RECOVERY PERIOD ELAPSED] X.XXXe-XX s since last recovery, attempting recovery #Y
```
- Prints when recovery period (50 μs) has passed
- Shows this is a repeat recovery attempt
- Indicates which attempt number this is

---

## Expected Log Output

### First Collapse
```
[RPE_COLLAPSE] Negative Rdot=-7.470e-03, attempting recovery (bubble collapsing)
               Time: 3.239393e-05 s, Last recovery: 0.000000e+00 s (3.239e-05 s ago)
               T_drop=293.06 K, T_sat(P_amb)=274.75 K, P_sat(T_drop)=8.468e+05 Pa, P_amb=4.518e+05 Pa
               R=3.812e-05 m, Ro=8.318e-05 m, dRdt=1.264e-03 m/s, dRdotdt=-8.882e+06 m/s²
               [RECOVERY #1] R_bubble: 3.812e-05 -> 1.191e-07 m (1.1*Rc=1.083e-07)
               [RECOVERY #1] R_drop: 8.318e-05 -> 8.250e-05 m (r_drop_0)
               [RECOVERY #1] num_drop: 9.263e-01 -> 3.260e-01 (mass conserved)
               [RECOVERY #1] Recovery time logged: 3.239393e-05 s, wait 5.000e-05 s before next
```

### During Recovery Period
```
[RPE_IN_RECOVERY] p_idx=0, recovery_count=1, time_since_recovery=1.234e-05 s, R=1.191e-07 m
[RPE_RECOVERY_WAIT] Still in recovery period: 1.234e-05 / 5.000e-05 s elapsed
                    Rdot=-2.345e-03 (negative but waiting), R=1.191e-07 m
```

### After Recovery Period (If Still Collapsing)
```
[RPE_IN_RECOVERY] p_idx=0, recovery_count=1, time_since_recovery=5.234e-05 s, R=1.191e-07 m
[RPE_COLLAPSE] Negative Rdot=-1.234e-03, attempting recovery (bubble collapsing)
               Time: 8.473782e-05 s, Last recovery: 3.239393e-05 s (5.234e-05 s ago)
               [RECOVERY PERIOD ELAPSED] 5.234e-05 s since last recovery, attempting recovery #2
               T_drop=293.05 K, T_sat(P_amb)=274.75 K, P_sat(T_drop)=8.467e+05 Pa, P_amb=4.518e+05 Pa
               [RECOVERY #2] R_bubble: 1.191e-07 -> 1.191e-07 m (1.1*Rc=1.083e-07)
               [RECOVERY #2] Recovery time logged: 8.473782e-05 s, wait 5.000e-05 s before next
```

### Successful Recovery
```
[RPE_IN_RECOVERY] p_idx=0, recovery_count=2, time_since_recovery=1.234e-05 s, R=1.191e-07 m
[RPE_RECOVERY_SUCCESS] Bubble recovered! Rdot=9.908e+00 m/s, 1.234e-05 s after recovery
                       Resetting recovery timer (was 8.473782e-05 s)
(recovery_time reset to 0, normal growth resumes)
```

---

## What to Look For in Logs

1. **First collapse:** See `[RPE_COLLAPSE]` with recovery #1
2. **Wait period:** See `[RPE_IN_RECOVERY]` and `[RPE_RECOVERY_WAIT]` messages
3. **Elapsed period:** See `[RECOVERY PERIOD ELAPSED]` with increasing attempt numbers
4. **Success:** See `[RPE_RECOVERY_SUCCESS]` when Rdot becomes positive

If you see:
- Only recovery #1 with no #2, #3, etc. → recovery period not elapsing properly
- No `[RPE_IN_RECOVERY]` messages → parcel not entering RPE after recovery
- `[RPE_RECOVERY_SUCCESS]` immediately → recovery is working!

---

## Benefits

1. **Cleaner logs:** No breakup/sampling spam
2. **Clear recovery tracking:** Can see exactly when parcels are in recovery mode
3. **Easy debugging:** Can trace full recovery sequence
4. **Attempt counting:** See how many recovery attempts each parcel needs

---

## Notes

- `[RPE_IN_RECOVERY]` prints limited to 20 messages to prevent spam
- Breakup messages still written to `breakup_debug.csv` file
- Logging can be re-enabled by uncommenting the printf statements
- All original functionality preserved, just output suppressed
