# Flag System Analysis and Simplification Proposal

**Date:** December 9, 2024  
**Issue:** Flag logic in `spray_drop_distort_NH3.c` has become overly complicated

---

## Current Flag System

### 1. **is_child** (int)
**Purpose:** Distinguish parent parcels from child parcels

**Values:**
- `0` = Parent parcel (injected, not yet broken up)
- `1` = Child parcel (created from breakup)

**Set in:**
- `parcel_prop.c` (injection): `0`
- `parcel_prop.c` (child creation): `1`
- `spray_drop_distort_NH3.c` (stuck parcel fixes): `1`
- `parcel_reset.c` (reset function): `1`

**Checked:** 8+ locations in `spray_drop_distort_NH3.c`

---

### 2. **thermal_breakup_flag** (int)
**Purpose:** Track thermal breakup state and why it was disabled

**Values (9 different states!):**
- `-1` = Active thermal breakup (normal operation)
- `-999` = Recovery attempted (temporary)
- `0` = ??? (appears in comparisons but unclear meaning)
- `3` = Bubble reached threshold, ready for breakup
- `4` = Child from thermal breakup (disabled)
- `5` = ??? (set in code)
- `6` = ??? (set in code)
- `888` = Recovery flag (bubble collapsed, trying to recover)
- `999` = Thermal breakup permanently disabled

**Set in:**
- `parcel_prop.c` (injection): `-1`
- `parcel_prop.c` (child creation): `4`
- `spray_drop_distort_NH3.c` (breakup): `3`, `5`, `6`, `888`, `-999`, `999`
- `parcel_reset.c` (reset): `999`
- `RPE_euler.c` (various conditions): `888`, others

**Checked:** 10+ locations

**PROBLEM:** Too many magic numbers, unclear semantics, overlapping purposes

---

### 3. **pbt** (Pre-Breakup Tag) (int)
**Purpose:** Enable/disable entry into thermal breakup routine

**Values:**
- `0` = Disabled (parcel cannot enter thermal breakup)
- `1` = Enabled (parcel can enter thermal breakup)

**Set in:**
- `parcel_prop.c` (injection): `1`
- `parcel_prop.c` (child creation): `0`
- `spray_drop_distort_NH3.c` (stuck parcel fix): `1`, `0`
- `parcel_reset.c` (reset): `0`

**Checked:** Multiple locations

**OVERLAP:** Redundant with `thermal_breakup_flag` - both control thermal breakup eligibility

---

### 4. **tbt** (Thermal Breakup Tag) (int)
**Purpose:** Mark parcels that have completed thermal breakup

**Values:**
- `0` = Not yet broken via thermal mechanism
- `1` = Has undergone thermal breakup

**Set in:**
- `parcel_prop.c` (injection): Not set
- `parcel_prop.c` (child creation): `0`
- `spray_drop_distort_NH3.c` (breakup): `1`, `0`
- `parcel_reset.c` (reset): Not set

**Checked:** Few locations

**OVERLAP:** Partially redundant with `thermal_breakup_flag = 3`

---

### 5. **film_flag** (int) - HIJACKED
**Purpose (original):** Film boiling model flag (CONVERGE built-in)

**Purpose (hijacked):** Store `is_child` for .h5 output (diagnostic)

**Values:**
- `1` = Child parcel (for visualization)

**Set in:**
- `spray_drop_distort_NH3.c`: `1`
- `parcel_prop.c` (child creation): `1`

**PROBLEM:** Using CONVERGE built-in variable for unintended purpose

---

### 6. **recovery_time** (double)
**Purpose:** Track time spent in bubble collapse recovery

**Values:**
- `0.0` = Not in recovery
- `> 0.0` = Recovery period remaining (seconds)

**Set in:**
- `RPE_euler.c` (recovery logic)

**Checked:** `RPE_euler.c` (recovery conditions)

---

### 7. **recovery_count** (int)
**Purpose:** Count number of recovery attempts

**Values:**
- `0` = No recovery attempts
- `> 0` = Number of times recovery attempted

**Set in:**
- `RPE_euler.c` (recovery logic)

**Checked:** `RPE_euler.c` (recovery diagnostics)

---

## Problems with Current System

### 1. **Too Many Flags with Overlapping Purpose**
- `is_child`, `pbt`, `tbt`, `thermal_breakup_flag` all relate to breakup state
- Confusing which flag to check in which context
- Easy to miss updating one flag when state changes

### 2. **Magic Numbers**
`thermal_breakup_flag` has 9 different values with unclear meanings:
- What's the difference between `4` and `999`?
- Why both `-999` and `888` for recovery?
- What do `5` and `6` represent?

### 3. **Redundant Checks**
Entry to thermal breakup requires:
```c
is_child == 0 && thermal_breakup_flag < 0 && pbt == 1
```
Three conditions checking essentially the same thing!

### 4. **Inconsistent Updates**
When disabling thermal breakup, must remember to set:
- `is_child = 1`
- `thermal_breakup_flag = 999`
- `pbt = 0`
- `tbt = 0` (sometimes)
- `film_flag = 1` (for diagnostics)

Miss one → bugs!

### 5. **Variable Hijacking**
Using `film_flag` for `is_child` output is a hack that will break if film model is enabled.

---

## Simplified Proposal: 2-Flag System

### Design Principle
**One flag for state, one flag for phase**

---

### **Flag 1: parcel_state** (int) - WHO AM I?

**Purpose:** Fundamental parcel identity

**Values:**
- `0` = PARENT (injected, not broken up)
- `1` = CHILD_KH_RT (created from KH-RT breakup)
- `2` = CHILD_THERMAL (created from thermal breakup)

**Replaces:** `is_child`, `tbt`, `film_flag` (hijack)

**Set once:** At creation (injection or breakup)

**Never changes:** Immutable after creation

**Benefits:**
- Clear identity: know exactly how parcel was created
- No confusion about breakup history
- Can distinguish thermal vs KH-RT children (useful for diagnostics)

---

### **Flag 2: breakup_phase** (int) - WHAT AM I DOING?

**Purpose:** Current thermal breakup phase (state machine)

**Values:**
```c
// Thermal breakup state machine
#define BREAKUP_DISABLED       0  // Cannot do thermal breakup
#define BREAKUP_ELIGIBLE       1  // Can enter thermal breakup
#define BREAKUP_ACTIVE         2  // Currently growing bubble
#define BREAKUP_RECOVERY       3  // Bubble collapsed, trying to recover
#define BREAKUP_READY          4  // Bubble ready to fragment
#define BREAKUP_COMPLETE       5  // Breakup executed, done
```

**Replaces:** `thermal_breakup_flag`, `pbt`, recovery state

**Benefits:**
- Clear state machine semantics
- Named constants instead of magic numbers
- Single source of truth for thermal breakup state
- Recovery is just another phase (no special flags)

---

## Proposed Variable Mapping

| Old Flags | New Flag | New Value |
|-----------|----------|-----------|
| `is_child=0` | `parcel_state` | `PARENT (0)` |
| `is_child=1` (KH-RT) | `parcel_state` | `CHILD_KH_RT (1)` |
| `is_child=1` (thermal) | `parcel_state` | `CHILD_THERMAL (2)` |
| `thermal_breakup_flag=-1, pbt=1` | `breakup_phase` | `BREAKUP_ELIGIBLE (1)` |
| `thermal_breakup_flag=-1, pbt=1` (in loop) | `breakup_phase` | `BREAKUP_ACTIVE (2)` |
| `thermal_breakup_flag=888` | `breakup_phase` | `BREAKUP_RECOVERY (3)` |
| `thermal_breakup_flag=3` | `breakup_phase` | `BREAKUP_READY (4)` |
| `thermal_breakup_flag=999, pbt=0` | `breakup_phase` | `BREAKUP_DISABLED (0)` |
| `thermal_breakup_flag=4` (child) | `breakup_phase` | `BREAKUP_COMPLETE (5)` |
| `tbt=1` | `parcel_state` | `PARENT` + `breakup_phase=COMPLETE` |
| `film_flag=1` | `parcel_state` | `CHILD_*` (anything != PARENT) |

---

## Code Simplification Examples

### Current Code (Complex)
```c
// Entry check - 3 flags!
if (old_parcel_cloud.is_child[p_idx] == 0 &&
    old_parcel_cloud.thermal_breakup_flag[p_idx] < 0 && 
    old_parcel_cloud.pbt[p_idx] == 1) {
    // Enter thermal breakup
}

// Disable thermal breakup - must set 5 things!
old_parcel_cloud.is_child[p_idx] = 1;
old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;
old_parcel_cloud.pbt[p_idx] = 0;
old_parcel_cloud.tbt[p_idx] = 0;
old_parcel_cloud.film_flag[p_idx] = 1;

// Check if in recovery - magic number!
if (old_parcel_cloud.thermal_breakup_flag[p_idx] == 888) {
    // Handle recovery
}
```

### Simplified Code (Proposed)
```c
// Entry check - 1 condition!
if (old_parcel_cloud.breakup_phase[p_idx] == BREAKUP_ELIGIBLE) {
    // Enter thermal breakup
    old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_ACTIVE;
}

// Disable thermal breakup - single assignment!
old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_DISABLED;
// parcel_state unchanged (stays PARENT)

// Check if in recovery - named constant!
if (old_parcel_cloud.breakup_phase[p_idx] == BREAKUP_RECOVERY) {
    // Handle recovery
}

// Check if child - clear semantics!
if (old_parcel_cloud.parcel_state[p_idx] != PARENT) {
    // This is a child parcel
}
```

---

## Implementation Header (globals.h)

```c
// Parcel state enumeration (immutable identity)
typedef enum {
    PARENT = 0,           // Injected parcel, not yet broken
    CHILD_KH_RT = 1,      // Created from KH-RT breakup
    CHILD_THERMAL = 2     // Created from thermal breakup
} ParcelState;

// Thermal breakup phase (state machine)
typedef enum {
    BREAKUP_DISABLED = 0,   // Cannot do thermal breakup
    BREAKUP_ELIGIBLE = 1,   // Can enter thermal breakup (ready)
    BREAKUP_ACTIVE = 2,     // Growing bubble (in sub-timestep loop)
    BREAKUP_RECOVERY = 3,   // Bubble collapsed, attempting recovery
    BREAKUP_READY = 4,      // Bubble reached threshold, ready to fragment
    BREAKUP_COMPLETE = 5    // Breakup executed, done
} BreakupPhase;
```

---

## Migration Strategy

### Phase 1: Add New Flags (Parallel)
1. Add `parcel_state` and `breakup_phase` to `load_spray_env.c`
2. Initialize both alongside old flags
3. Update both systems in parallel (redundant)

### Phase 2: Migrate Code Sections
1. Start with `parcel_prop.c` (initialization)
2. Then `spray_drop_distort_NH3.c` (main logic)
3. Then `RPE_euler.c`, `RPE_song.c` (solvers)
4. Update `Breakup.c`, `parcel_reset.c` (utilities)

### Phase 3: Remove Old Flags
1. Verify all code uses new flags
2. Remove old flag checks
3. Keep old flags in `load_spray_env.c` for .h5 compatibility (can output but not use)

---

## Recovery Tracking (Optional 3rd Flag)

If recovery needs more detail:

### **Flag 3: recovery_state** (int) - RECOVERY DETAILS

**Values:**
```c
#define RECOVERY_NONE      0  // Not in recovery
#define RECOVERY_ACTIVE    1  // Currently in recovery period
#define RECOVERY_FAILED    2  // Recovery exceeded max attempts
```

**Benefits:**
- Separate concern (not part of main state machine)
- Only checked in `RPE_euler.c`
- Could replace `recovery_time` + `recovery_count` combo

---

## Comparison: Old vs New

### Number of Flags
- **Old:** 7 flags (is_child, pbt, tbt, thermal_breakup_flag, film_flag, recovery_time, recovery_count)
- **New:** 2-3 flags (parcel_state, breakup_phase, [optional: recovery_state])

### Number of States
- **Old:** 9 values for thermal_breakup_flag alone
- **New:** 3 parcel states, 6 breakup phases (all named, clear purpose)

### Entry Condition
- **Old:** `is_child==0 && thermal_breakup_flag<0 && pbt==1` (3 checks)
- **New:** `breakup_phase==BREAKUP_ELIGIBLE` (1 check)

### Disable Thermal Breakup
- **Old:** Set 5 different flags
- **New:** Set 1 flag: `breakup_phase = BREAKUP_DISABLED`

### Code Clarity
- **Old:** Magic numbers everywhere (888, 999, -999, 4, 3...)
- **New:** Named constants (BREAKUP_RECOVERY, BREAKUP_READY, etc.)

---

## Benefits Summary

✅ **Fewer flags:** 7 → 2-3  
✅ **Clear semantics:** Named enums instead of magic numbers  
✅ **No redundancy:** Single source of truth for each concept  
✅ **Easier to understand:** State machine is explicit  
✅ **Fewer bugs:** Can't forget to update multiple flags  
✅ **Better diagnostics:** `parcel_state` distinguishes thermal vs KH-RT children  
✅ **Maintainable:** New developers can understand quickly  

---

## State Machine Diagram (Proposed)

```
PARENT (parcel_state)
   │
   │ (Superheated, meets conditions)
   ↓
breakup_phase: ELIGIBLE → ACTIVE → READY → COMPLETE
                   ↓          ↑
                   ↓          │ (Recovery attempt)
                   └─ RECOVERY ─┘
                   │
                   │ (Subcooled, too small, etc.)
                   ↓
                DISABLED

After breakup:
   parcel_state → CHILD_THERMAL
   breakup_phase → DISABLED (children can't break again)
```

---

## Questions for Consideration

1. **Do we need to distinguish CHILD_THERMAL vs CHILD_KH_RT?**
   - Pro: Better diagnostics, can track which mechanism dominates
   - Con: Adds complexity if not needed

2. **Should recovery be part of breakup_phase or separate?**
   - Current proposal: Part of phase (BREAKUP_RECOVERY)
   - Alternative: Separate `recovery_state` flag (cleaner separation)

3. **What about Song model?**
   - Same system works: Song uses `breakup_phase` identically
   - Song vs thermal is controlled by `use_song_rpe` flag, not parcel state

4. **Backward compatibility for .h5 output?**
   - Keep old flags in `load_spray_env.c` but populate from new flags
   - Or: Add new output variables and deprecate old ones

---

## Recommendation

**Implement 2-flag system:**
- `parcel_state` (identity: PARENT, CHILD_KH_RT, CHILD_THERMAL)
- `breakup_phase` (thermal breakup state machine: 6 phases)

**Migration:**
- Parallel implementation first (keep old flags working)
- Migrate code section by section
- Test at each step
- Remove old flags last

**Timeline:**
- Phase 1 (add new flags): 1 hour
- Phase 2 (migrate code): 4-6 hours (systematic replacement)
- Phase 3 (remove old): 1 hour

**Total:** ~1 day of careful refactoring

---

**Would this simplification make the code more maintainable?**
