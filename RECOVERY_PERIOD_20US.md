# Recovery Period Reduced + Additional Diagnostics

**Date:** 2025-11-17  
**Status:** Implemented for testing

---

## Changes Made

### 1. Recovery Period Reduced
Changed from **50 microseconds to 20 microseconds**:
```c
const CONVERGE_precision_t RECOVERY_PERIOD = 2.0e-5;  // 20 μs (was 5e-5)
```

### 2. Enhanced Recovery Success Message
Added more detail to see parcel state after recovery:
```c
[RPE_RECOVERY_SUCCESS] Bubble recovered! Rdot=X.XXXe+XX m/s, X.XXXe-XX s after recovery
                       R_bubble=X.XXXe-XX m, R_drop=X.XXXe-XX m, recovery_count was X
                       Resetting recovery timer (was X.XXXe-XX s), flag=X
```

### 3. Flag Reset on Recovery Success
Added explicit flag reset to ensure normal operation resumes:
```c
old_parcel_cloud->thermal_breakup_flag[p_idx] = -1;  // Reset to normal
```

### 4. Enhanced Repeat Collapse Message
More alarming message if parcel keeps collapsing:
```c
[REPEAT COLLAPSE] X.XXXe-XX s since last recovery (period=2.000e-05 s), attempting recovery #X
[REPEAT COLLAPSE] Parcel keeps collapsing! Check conditions.
```

---

## Analysis of Previous Log Output

Looking at your earlier log:
```
Cycle 82:  [RPE_COLLAPSE] ... [RECOVERY #1] ... time logged: 3.239393e-05 s
Cycle 83:  [RPE_RECOVERY_SUCCESS] Rdot=9.908e+00 m/s, 6.143e-08 s after recovery
```

**Key Observation:** The recovery "succeeded" after only **61 nanoseconds** (6.143e-08 s)!

This suggests:
1. Recovery was applied: R_bubble reset to 1.19e-07 m (119 nm)
2. On the very next timestep, Rdot was already positive (9.9 m/s)
3. Recovery timer was reset to 0
4. Parcel resumed normal growth

**But you still see rogue parcels in Tecplot?**

---

## Possible Issues

### Issue 1: Parcel Collapses Again Later
- Recovery succeeds temporarily
- Parcel grows for a while
- Hits another pressure wave and collapses again
- Because `recovery_time` was reset to 0, it counts as a "new" first recovery
- This could repeat indefinitely

**What to look for:** Multiple `[RPE_COLLAPSE]` messages with recovery #1 for the same parcel at different times

### Issue 2: Recovered Bubble Too Small
- R_critical ≈ 1.08e-07 m (108 nm)
- Recovery sets R_bubble = 1.1 * R_c ≈ 1.19e-07 m (119 nm)
- This is **extremely small** compared to typical bubble sizes (100-200 μm)
- Maybe the bubble can't grow fast enough from this tiny size?

**What to look for:** Does the bubble actually grow after recovery success?

### Issue 3: Flag Not Reset Properly (NOW FIXED)
- Previously, flag might have been stuck at -999
- Now explicitly reset to -1 on recovery success

### Issue 4: Parcel Never Attempts Breakup
- Parcel recovers and grows
- But never reaches the breakup criterion?
- Stays as a large parent with small bubble?

**What to look for:** Check if recovered parcels ever reach breakup

---

## What to Check in Next Run

### 1. Does recovery actually happen multiple times for the SAME parcel?
Look for:
```
p_idx=X: [RPE_COLLAPSE] ... [RECOVERY #1] ... time=T1
p_idx=X: [REPEAT COLLAPSE] ... [RECOVERY #2] ... time=T2
p_idx=X: [REPEAT COLLAPSE] ... [RECOVERY #3] ... time=T3
```

If you DON'T see this, the problem is NOT the recovery period - it's something else.

### 2. What happens after recovery success?
Look for:
```
[RPE_RECOVERY_SUCCESS] ... R_bubble=1.191e-07 m, R_drop=8.250e-05 m
(20-50 μs later)
[BREAKUP] Parcel X: ... R_bubble=??? m
```

Does the bubble grow? Does it break up?

### 3. Check recovery_count in final rogue parcels
In Tecplot, visualize the `recovery_count` field:
- If count = 0: Parcel never collapsed
- If count = 1: Parcel collapsed once, recovered, then grew normally
- If count > 1: Parcel collapsed multiple times

---

## Hypothesis

Based on your log showing immediate recovery success, I suspect:

**The recovery IS working** - the parcel recovers and starts growing with positive Rdot.

**BUT** one of these might be happening:
1. Bubble grows very slowly from tiny size (119 nm) and takes too long to reach breakup criterion
2. Parcel hits another pressure wave later and collapses AGAIN (new recovery #1)
3. Something else prevents the parcel from breaking up despite having a growing bubble

---

## Testing Strategy

Run the simulation with these changes and check:

1. **Count total recovery attempts:**
   ```bash
   grep "RECOVERY #" log | wc -l
   ```

2. **Count unique parcels that recovered:**
   ```bash
   grep "RECOVERY #1" log | wc -l
   ```

3. **Check if any parcel has multiple recoveries:**
   ```bash
   grep "REPEAT COLLAPSE" log
   ```

4. **Check recovery success rate:**
   ```bash
   grep "RPE_RECOVERY_SUCCESS" log | wc -l
   ```

If you see:
- Many RECOVERY #1 but no REPEAT COLLAPSE → Each parcel collapses once, recovers, then something else happens
- REPEAT COLLAPSE messages → Parcels are collapsing multiple times (recovery period issue OR environment issue)
- No REPEAT COLLAPSE and few recoveries → Recovery working fine, rogue parcels have different cause

---

## Notes

- Recovery period now 20 μs (was 50 μs) - should see multiple attempts faster if needed
- Flag explicitly reset on success - should prevent any "stuck flag" issues
- Enhanced logging - should clearly show if parcels are collapsing repeatedly
- Recovery counter tracks history even after success - can check in Tecplot

---

## Next Steps if Issue Persists

If you still see rogue parcels after this:

1. **Check if they ever attempted recovery** - look at recovery_count field
2. **Check if recovery succeeded** - look for RPE_RECOVERY_SUCCESS messages
3. **Track a single rogue parcel through the log** - see its full history
4. **Visualize R_bubble vs time** for a recovered parcel - see if it actually grows

The problem might not be the recovery system - it might be that recovered parcels can't grow fast enough or hit other issues later.
