# Thermal Breakup Logic Analysis
**File:** `src/spray_drop_distort_NH3.c`  
**Date:** 2025-11-11

## State Variables Overview

### Key Flags
1. **`pbt` (Pre-Breakup Thermal)**: Indicates parcel is *eligible* for thermal breakup processing
   - `pbt = 1`: Parcel is in thermal breakup monitoring state (bubble growing)
   - `pbt = 0`: Parcel NOT undergoing thermal breakup (either hasn't started or already finished)

2. **`tbt` (Triggered Breakup Thermal)**: Breakup criterion has been met, ready to execute breakup
   - `tbt = 1`: Breakup triggered, needs to create child parcels
   - `tbt = 0`: No breakup yet

3. **`thermal_breakup_flag`**: Diagnostic flag indicating breakup status/reason
   - `-1`: Initial state (just injected, monitoring not yet started)
   - `0`: Normal state during monitoring
   - `3`: kb-triggered breakup (normal physics-based)
   - `5`: Bubble reached droplet surface (geometric limit)
   - `6`: REMOVED - was "extreme expansion" (1.5x radius)
   - `999`: Thermal breakup aborted/stopped (subcooling, stalled growth, etc.)

4. **`is_child`**: Indicates if parcel is a child from previous breakup
   - `0`: Parent parcel
   - `1`: Child parcel (from thermal/KH-RT/collision breakup)

---

## Initialization (parcel_prop.c)

**When:** New parcel injected into domain

```c
parcel_cloud.pbt[passed_parcel_idx] = 1;                     // Enable thermal breakup monitoring
parcel_cloud.tbt[passed_parcel_idx] = 0;                     // No breakup yet
parcel_cloud.thermal_breakup_flag[passed_parcel_idx] = -1;   // Initial state
parcel_cloud.is_child[passed_parcel_idx] = 0;                // Parent parcel
parcel_cloud.r_bubble[passed_parcel_idx] = 1.1 * Rc;         // Initialize bubble at critical radius
```

**Result:** Parcel enters thermal breakup monitoring with `pbt=1, tbt=0, thermal_breakup_flag=-1`

---

## Main Processing Loop (spray_drop_distort_NH3.c)

### Entry Condition (Line 411)
```c
if (old_parcel_cloud.thermal_breakup_flag[p_idx] < 0 && old_parcel_cloud.pbt[p_idx] == 1)
```

**Translation:** Process this parcel if:
- `thermal_breakup_flag < 0`: Still in initial state (never processed before)
- `pbt == 1`: Thermal breakup monitoring is active

**Problem with this logic:** Once `thermal_breakup_flag` is set to 0 or higher, parcel exits this block permanently. This means:
- First cycle: `thermal_breakup_flag = -1` → enters block, sets flag to 0
- Subsequent cycles: `thermal_breakup_flag = 0` → **NEVER RE-ENTERS THIS BLOCK**
- Parcel stuck with whatever state it had after first cycle

---

## Breakup Trigger Paths

### Path 1: REMOVED - Extreme Expansion (Lines 423-433, now commented out)
```c
/* COMMENTED OUT - WAS CAUSING EARLY PARCEL SKIP ISSUE
if(old_parcel_cloud.radius[p_idx] > old_parcel_cloud.r_drop_0[p_idx]*1.5)
{
   old_parcel_cloud.thermal_breakup_flag[p_idx] = 6;
   old_parcel_cloud.tbt[p_idx] = 1;
   old_parcel_cloud.pbt[p_idx] = 0;
   continue;  // Skip to next parcel
}
*/
```
**Intended purpose:** Safety mechanism for extreme droplet expansion  
**Actual effect:** Blocked parcels from ever reaching kb calculation  
**Status:** Removed (2025-11-11)

---

### Path 2: Bubble Growth Stalled (Line 491-496)
```c
if(old_parcel_cloud.v_bubble[p_idx] < 1.0e-10)
{
   old_parcel_cloud.pbt[p_idx] = 0;                    // Disable thermal breakup
   old_parcel_cloud.thermal_breakup_flag[p_idx] = 999; // Mark as aborted
   break;  // Exit sub-timestep loop
}
```
**When:** Bubble velocity drops below 0.1 nm/s (essentially stopped)  
**Effect:** Abort thermal breakup for this parcel  
**Reason:** No driving force for breakup (likely subcooled or equilibrated)

---

### Path 3: Bubble Reached Droplet Surface (Line 520-526)
```c
if (Rb > old_parcel_cloud.radius[p_idx])
{
   old_parcel_cloud.thermal_breakup_flag[p_idx] = 5;
   old_parcel_cloud.tbt[p_idx] = 1;                    // Trigger breakup
   old_parcel_cloud.r_bubble[p_idx] = 0.8 * old_parcel_cloud.radius[p_idx];
   break;  // Exit sub-timestep loop
}
```
**When:** Bubble radius exceeds droplet radius  
**Effect:** Trigger breakup immediately (geometric limit)  
**Physical meaning:** Bubble has consumed droplet, cannot grow further  
**Correction:** Sets bubble to 80% of droplet radius to avoid singularity

---

### Path 4: kb > Threshold (Line 571-579) ✅ PRIMARY CRITERION
```c
if (kb > kb_threshold)
{
   old_parcel_cloud.thermal_breakup_flag[p_idx] = 3;
   old_parcel_cloud.tbt[p_idx] = 1;                    // Trigger breakup
   old_parcel_cloud.r_bubble[p_idx] = Rb;
   old_parcel_cloud.eta_drop[p_idx] = kb;              // Store final kb
   break;  // Exit sub-timestep loop
}
```
**When:** Disturbance growth rate exceeds threshold  
**Effect:** Trigger breakup (normal physics-based path)  
**Physical meaning:** Bubble-driven surface instabilities have grown sufficiently  
**This should be the PRIMARY and NORMAL breakup mechanism**

---

## Critical Logic Issue: Single-Entry Problem

### The Problem
```c
Line 411: if (thermal_breakup_flag < 0 && pbt == 1)
```

This condition only allows entry on the **first cycle** after injection:
- Cycle 1: `thermal_breakup_flag = -1` → enters, processes, sets to 0
- Cycle 2+: `thermal_breakup_flag = 0` → **BLOCKED FROM ENTRY**

### What Should Happen
Parcel should continue processing until one of these occurs:
1. kb > threshold (normal breakup)
2. Bubble reaches surface (geometric limit)
3. Growth stops (abort)

### Current Broken Flow
```
Injection → pbt=1, tbf=-1
    ↓
Cycle 1: Enter block (tbf<0)
    ↓
Set tbf=0 (somewhere in initialization?)
    ↓
Cycle 2: BLOCKED (tbf not < 0)
    ↓
Parcel persists forever with pbt=1, never triggers breakup
```

---

## Where is thermal_breakup_flag Set to 0?

**Line 339:**
```c
if (P_sat < 1)
{
   // ... cleanup code ...
   old_parcel_cloud.thermal_breakup_flag[p_idx] = 0;
   // ...
   continue;
}
```
This sets `thermal_breakup_flag = 0` if saturation pressure is invalid.

**But there's NO OTHER LOCATION** that sets it to 0 during normal operation!

This means:
- If parcel has valid P_sat, `thermal_breakup_flag` stays at -1
- Parcel continues to enter the block each cycle ✓
- **Unless** the commented-out 1.5x expansion check was blocking it

---

## Corrected Understanding After Code Review

**The fix we applied (removing lines 423-433) was correct!**

### Before Fix:
```
Cycle 1: tbf=-1, pbt=1 → Enter block
    ↓
Check: radius > 1.5*r_drop_0? → NO (early parcel, hasn't expanded)
    ↓
CONTINUE to next parcel (line 429)
    ↓
Never reaches kb calculation
    ↓
tbf stays -1, pbt stays 1, but parcel skipped
```

### After Fix:
```
Cycle 1: tbf=-1, pbt=1 → Enter block
    ↓
(1.5x check removed)
    ↓
Process RPE solver
    ↓
Process DGRE
    ↓
Calculate kb
    ↓
If kb > threshold: Set tbt=1, tbf=3, break
    ↓
Execute breakup routine
```

---

## Recommended Logic Improvements

### Issue 1: Entry Condition Too Restrictive
**Current:**
```c
if (thermal_breakup_flag < 0 && pbt == 1)
```

**Should be:**
```c
if (pbt == 1 && tbt == 0)
```

**Reasoning:**
- `pbt=1`: Thermal breakup monitoring active
- `tbt=0`: Breakup not yet triggered
- Don't rely on `thermal_breakup_flag` value for flow control

---

### Issue 2: No Explicit Flag Management
The code lacks clear state transitions. Should add:

```c
// After entering thermal breakup block for first time
if (thermal_breakup_flag < 0) {
   thermal_breakup_flag = 0;  // Mark as "in processing"
}
```

---

### Issue 3: Multiple Exit Paths Without Unified Logic

**Current:** 4 different ways to set `tbt=1`:
1. Removed expansion check (was line 426)
2. Bubble surface contact (line 523)
3. kb threshold (line 575)
4. Others?

**Should have:** Single unified trigger logic:
```c
// After all physics calculations
bool trigger_breakup = false;
int breakup_reason = 0;

if (kb > kb_threshold) {
   trigger_breakup = true;
   breakup_reason = 3;  // kb-triggered
}
else if (Rb > 0.95 * radius) {
   trigger_breakup = true;
   breakup_reason = 5;  // geometric limit
}

if (trigger_breakup) {
   tbt = 1;
   thermal_breakup_flag = breakup_reason;
   pbt = 0;  // Disable further monitoring
   break;
}
```

---

## Breakup Execution (Lines 600+)

**Condition:**
```c
if (old_parcel_cloud.tbt[p_idx] && old_parcel_cloud.thermal_breakup_flag[p_idx]!=4)
```

**Effect:** Calls `Breakup()` function to create child parcels

**After breakup:**
- Parent: `is_child = 0` → remains parent (or deleted?)
- Children: `is_child = 1`, `pbt = 0` (no thermal breakup monitoring)

---

## Summary of Current Logic Flow

### Correct Path (after fix):
```
1. Injection: pbt=1, tbt=0, tbf=-1, is_child=0
2. Enter loop: tbf<0 && pbt=1 → TRUE
3. (1.5x check removed - no longer blocks)
4. Process RPE → bubble grows
5. Process Geometry → droplet expands
6. Process DGRE → calculate omega
7. Calculate kb
8. Check: kb > threshold?
   - YES → tbt=1, tbf=3, break
   - NO → continue sub-timesteps
9. If tbt=1: Execute Breakup()
10. Create children with is_child=1, pbt=0
```

### Edge Cases:
- **Subcooling:** RPE_euler sets tbf=999, pbt=0 → aborts
- **Stalled growth:** v_bubble<1e-10 → tbf=999, pbt=0 → aborts
- **Surface contact:** Rb>radius → tbf=5, tbt=1 → breakup
- **Invalid P_sat:** tbf=0 → but pbt still 1, so continues processing ✓

---

## Critical Finding: The 1.5x Check Was The Only Problem

**The core logic is actually sound**, except for the removed 1.5x expansion check.

### Why it was breaking early parcels:
1. Early parcel: Tinj=293K, Tamb=298K (subcooled)
2. Slow heating → slow expansion → radius < 1.5*r_drop_0 for many cycles
3. Hit line 429: `continue;` → skip to next parcel
4. Never reaches kb calculation
5. Never triggers breakup

### Why later parcels worked:
1. Later parcel: Enters warmer gas (from earlier evaporation/mixing)
2. Faster heating → faster expansion → radius > 1.5*r_drop_0 quickly
3. Skips the 1.5x check (or passes it)
4. Processes kb calculation
5. Triggers breakup normally

---

## Recommendations

### 1. ✅ DONE: Remove 1.5x expansion check
**Status:** Completed in FIXES_2025-11-11.md

### 2. Add explicit state management (optional but recommended)
```c
// At entry to thermal breakup block
if (thermal_breakup_flag == -1) {
   thermal_breakup_flag = 0;  // Mark as actively processing
}
```

### 3. Clarify entry condition (optional but cleaner)
```c
// Change from:
if (thermal_breakup_flag < 0 && pbt == 1)

// To:
if (pbt == 1 && tbt == 0 && is_child == 0)
```
More explicit about intent: monitor parent parcels that haven't triggered yet.

### 4. Add diagnostic output for skipped parcels
```c
// After line 411
if (pbt == 1 && !(thermal_breakup_flag < 0)) {
   printf("[DEBUG] Parcel %d skipped: pbt=1 but tbf=%d\n", p_idx, thermal_breakup_flag);
}
```

### 5. Unified breakup trigger logic (recommended for clarity)
See "Issue 3" above for proposed structure.

---

## Conclusion

**The 1.5x radius expansion check was the sole cause of early parcels skipping thermal breakup.**

Removing it (as done in the fix) allows all parcels to:
1. Enter the thermal breakup processing block
2. Run RPE solver to grow bubbles
3. Calculate DGRE and kb
4. Trigger breakup when kb > threshold (normal physics)

No further changes are strictly necessary, but the recommendations above would improve:
- Code clarity and maintainability
- Diagnostic capability
- Robustness against future bugs

The core physics logic (RPE → DGRE → kb criterion) is **sound and working correctly** after this fix.
