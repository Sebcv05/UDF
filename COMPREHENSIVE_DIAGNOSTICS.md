# Comprehensive Recovery Diagnostics

**Date:** 2025-11-17  
**Status:** Diagnostics added for debugging

---

## Purpose

Track exactly what happens to parcels in recovery mode to find why they stop attempting bubble growth after the initial recovery.

---

## New Diagnostic Messages

### 1. When pbt Gets Set to 0 in RPE (Parcel Killed)

**`[RPE_KILL_IN_RECOVERY]`** - Logs when a parcel in recovery mode gets killed by a safety check

#### Location: RPE_euler.c

**a) Droplet too small:**
```c
[RPE_KILL_IN_RECOVERY] p_idx=X, Reason: Ro too small (X.XXe-XX m), recovery_time=X.XXe-XX s, recovery_count=X
```

**b) Negative pressure difference:**
```c
[RPE_KILL_IN_RECOVERY] p_idx=X, Reason: Negative P_sat-P_amb (X.XXe-XX Pa), recovery_time=X.XXe-XX s, recovery_count=X
```

### 2. When Recovery Checks Are Skipped (Parcel Protected)

**`[RPE_SUBCOOL_IN_RECOVERY]`** - Parcel would be subcooled but is protected:
```c
[RPE_SUBCOOL_IN_RECOVERY] p_idx=X, T_drop=XXX.XX K < T_sat=XXX.XX K, but protected (recovery_time=X.XXe-XX s)
```

**`[RPE_SMALL_VEL_IN_RECOVERY]`** - Parcel has small velocity but is protected:
```c
[RPE_SMALL_VEL_IN_RECOVERY] p_idx=X, Rdot=X.XXe-XX < 1e-10, but protected (recovery_time=X.XXe-XX s)
```

### 3. Recovery Flow Control

**`[RECOVERY_BREAK]`** - When flag=888 causes break from sub-loop:
```c
[RECOVERY_BREAK] p_idx=X, flag 888->-999, breaking sub-loop, recovery_time=X.XXe-XX s, R_bubble=X.XXe-XX m
```

**`[RPE_ENTRY_IN_RECOVERY]`** - Parcel in recovery enters RPE block:
```c
[RPE_ENTRY_IN_RECOVERY] p_idx=X, flag=X, pbt=X, time_since_recovery=X.XXe-XX s, R_bubble=X.XXe-XX m
```

---

## What To Look For in Logs

### Successful Recovery Sequence:
```
Cycle N:   [RPE_COLLAPSE] ... [RECOVERY #1]
           [RECOVERY_BREAK] flag 888->-999
Cycle N+1: [RPE_ENTRY_IN_RECOVERY] flag=-999, pbt=1
           [RPE_IN_RECOVERY] bubble growing
Cycle N+2: [RPE_ENTRY_IN_RECOVERY] flag=-999, pbt=1
           [RPE_IN_RECOVERY] bubble growing
...
Cycle N+X: [RPE_RECOVERY_SUCCESS] time_since_recovery > 2e-5 s
```

### Failed Recovery - Parcel Killed:
```
Cycle N:   [RPE_COLLAPSE] ... [RECOVERY #1]
           [RECOVERY_BREAK] flag 888->-999
Cycle N+1: [RPE_ENTRY_IN_RECOVERY] flag=-999, pbt=1
           [RPE_KILL_IN_RECOVERY] Reason: ... ← PROBLEM HERE
```

### Failed Recovery - Parcel Blocked:
```
Cycle N:   [RPE_COLLAPSE] ... [RECOVERY #1]
           [RECOVERY_BREAK] flag 888->-999
Cycle N+1: [RPE_BLOCKED_IN_RECOVERY] flag=XXX, pbt=X ← pbt changed or flag changed!
```

### Protected Parcel (Working as Intended):
```
Cycle N:   [RPE_ENTRY_IN_RECOVERY] ...
           [RPE_SUBCOOL_IN_RECOVERY] but protected
           [RPE_IN_RECOVERY] bubble still growing (protection working)
```

---

## Debug Strategy

Run the simulation and grep for these patterns:

### 1. Find if parcels are being killed:
```bash
grep "RPE_KILL_IN_RECOVERY" log
```
- If you see this, a safety check is killing recovered parcels
- Check which reason appears

### 2. Find if parcels can't enter RPE:
```bash
grep "RPE_BLOCKED_IN_RECOVERY" log
```
- If you see this, either `flag` or `pbt` changed unexpectedly
- The message shows what they changed to

### 3. Check if protection is working:
```bash
grep "SUBCOOL_IN_RECOVERY\|SMALL_VEL_IN_RECOVERY" log
```
- If you see these, protection is activating correctly
- Parcels should continue growing despite adverse conditions

### 4. Track full recovery lifecycle:
```bash
grep "p_idx=1.*RECOVERY\|p_idx=1.*RPE" log
```
- Replace "1" with specific parcel index
- Shows full sequence for one parcel

---

## Expected Output Patterns

### If Recovery Is Working:
- Many `[RPE_ENTRY_IN_RECOVERY]` messages
- Many `[RPE_IN_RECOVERY]` messages
- Some `[RPE_RECOVERY_SUCCESS]` messages after 20+ μs
- Possibly some `[RPE_SUBCOOL_IN_RECOVERY]` (protection activating)
- **NO** `[RPE_KILL_IN_RECOVERY]` messages
- **NO** `[RPE_BLOCKED_IN_RECOVERY]` messages

### If Something Is Breaking Recovery:
- `[RPE_KILL_IN_RECOVERY]` → A safety check is killing parcels
- `[RPE_BLOCKED_IN_RECOVERY]` → Flag or pbt changed unexpectedly
- Few or no `[RPE_ENTRY_IN_RECOVERY]` after first cycle
- No `[RPE_RECOVERY_SUCCESS]` messages

---

## Files Modified

**`src/RPE_euler.c`:**
- Added `[RPE_KILL_IN_RECOVERY]` to 2 locations where pbt=0 for small Ro or negative P_sat-P_amb
- Added `[RPE_SUBCOOL_IN_RECOVERY]` to log when subcooled check is skipped
- Added `[RPE_SMALL_VEL_IN_RECOVERY]` to log when small velocity check is skipped

**`src/spray_drop_distort_NH3.c`:**
- Added `[RECOVERY_BREAK]` when flag 888→-999 transition happens
- Added `[RPE_ENTRY_IN_RECOVERY]` when parcel in recovery enters RPE block
- Note: `[RPE_BLOCKED_IN_RECOVERY]` was removed (would require structural changes)

---

## Diagnostic Limits

To prevent log spam, each diagnostic is limited:
- `[RPE_KILL_IN_RECOVERY]`: No limit (we need to see all kills)
- `[RPE_SUBCOOL_IN_RECOVERY]`: 5 messages
- `[RPE_SMALL_VEL_IN_RECOVERY]`: 5 messages
- `[RECOVERY_BREAK]`: 10 messages
- `[RPE_ENTRY_IN_RECOVERY]`: 20 messages

These can be adjusted if needed.

---

## Next Steps

1. Compile with these diagnostics
2. Run simulation
3. Search log for the new diagnostic messages
4. Identify where/why recovered parcels stop entering RPE
5. Fix the identified issue
6. Remove or reduce diagnostics once fixed

---

## Notes

- These diagnostics only apply to parcels with `recovery_time > 0`
- Normal parcels (not in recovery) won't generate these messages
- The diagnostics show the full lifecycle from recovery → success or failure
- Once the issue is identified, these printf statements can be commented out
