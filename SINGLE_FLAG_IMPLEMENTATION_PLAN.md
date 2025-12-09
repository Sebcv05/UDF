# Single-Flag Implementation: Complete Change List

**Branch Created:** `Song_v1` (backup of current state)  
**Date:** December 9, 2024

---

## Summary

**Changes Required:**
- **37** locations reference `is_child`
- **53** locations reference `thermal_breakup_flag`
- **30** locations reference `pbt`
- **26** locations reference `tbt`
- **Total:** ~146 changes across codebase

---

## Files to Modify

### 1. **include/env.h** (or globals.h)
**Purpose:** Define new enum

**Add:**
```c
// Breakup phase state machine (replaces is_child, pbt, tbt, thermal_breakup_flag)
typedef enum {
    BREAKUP_DISABLED  = 0,  // Parent, not eligible (subcooled, too small, etc.)
    BREAKUP_ELIGIBLE  = 1,  // Parent, superheated, ready to enter
    BREAKUP_ACTIVE    = 2,  // Parent, growing bubble (in sub-timestep loop)
    BREAKUP_RECOVERY  = 3,  // Parent, recovering from bubble collapse
    BREAKUP_READY     = 4,  // Parent, bubble at threshold, ready to fragment
    BREAKUP_COMPLETE  = 5   // CHILD (result of breakup, any mechanism)
} BreakupPhase;
```

**Changes:** 1 addition

---

### 2. **src/load_spray_env.c**
**Purpose:** Register new variable, remove old ones

**Line ~57:** Register new variable
```c
// NEW: Single breakup phase variable
CONVERGE_variable_register("breakup_phase", CONVERGE_INT, 
                           DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);

// REMOVE these:
// CONVERGE_variable_register("thermal_breakup_flag", CONVERGE_INT, ...);
// CONVERGE_variable_register("tbt", CONVERGE_INT, ...);
// CONVERGE_variable_register("is_child", CONVERGE_INT, ...);
// CONVERGE_variable_register("pbt", CONVERGE_INT, ...);
```

**Line ~127:** Get field ID
```c
// NEW
BREAKUP_PHASE = CONVERGE_lagrangian_field_id("breakup_phase");

// REMOVE these:
// THERMAL_BREAKUP_FLAG = CONVERGE_lagrangian_field_id("thermal_breakup_flag");
// TBT = CONVERGE_lagrangian_field_id("tbt");
// IS_CHILD = CONVERGE_lagrangian_field_id("is_child");
// PBT = CONVERGE_lagrangian_field_id("pbt");
```

**Line ~476:** Load into parcel_cloud structure
```c
// NEW
parcel_cloud_loc->breakup_phase = (int *)CONVERGE_cloud_get_field_data(c, BREAKUP_PHASE);

// REMOVE these:
// parcel_cloud_loc->tbt = (int *)CONVERGE_cloud_get_field_data(c, TBT);
// parcel_cloud_loc->is_child = (int *)CONVERGE_cloud_get_field_data(c, IS_CHILD);
// parcel_cloud_loc->pbt = (int *)CONVERGE_cloud_get_field_data(c, PBT);
// parcel_cloud_loc->thermal_breakup_flag = (int *)CONVERGE_cloud_get_field_data(c, THERMAL_BREAKUP_FLAG);
```

**Changes:** ~12 lines (add 3, remove 9)

---

### 3. **include/parcel_cloud.h** (or wherever ParcelCloud struct is defined)
**Purpose:** Update structure definition

```c
struct ParcelCloud {
    // ... existing fields ...
    
    // NEW: Single breakup phase
    int *breakup_phase;
    
    // REMOVE these:
    // int *is_child;
    // int *pbt;
    // int *tbt;
    // int *thermal_breakup_flag;
    
    // ... rest of fields ...
};
```

**Changes:** ~4 lines (add 1, remove 4)

---

### 4. **src/parcel_prop.c** (Parcel Injection/Creation)
**Purpose:** Initialize breakup_phase for new parcels

#### At injection (line ~187)
**Current:**
```c
parcel_cloud.is_child[passed_parcel_idx] = 0;
parcel_cloud.pbt[passed_parcel_idx] = 1;
parcel_cloud.thermal_breakup_flag[passed_parcel_idx] = -1;
```

**New:**
```c
parcel_cloud.breakup_phase[passed_parcel_idx] = BREAKUP_ELIGIBLE;
```

#### At child creation from KH-RT (line ~286)
**Current:**
```c
parcel_cloud.is_child[child_idx] = 1;
parcel_cloud.pbt[child_idx] = 0;
parcel_cloud.tbt[child_idx] = 0;
parcel_cloud.thermal_breakup_flag[child_idx] = 4;
```

**New:**
```c
parcel_cloud.breakup_phase[child_idx] = BREAKUP_COMPLETE;
```

#### At child creation from thermal breakup (line ~301)
**Current:**
```c
parcel_cloud.pbt[child_idx] = 0;
parcel_cloud.tbt[child_idx] = 0;
parcel_cloud.thermal_breakup_flag[child_idx] = 4;
```

**New:**
```c
parcel_cloud.breakup_phase[child_idx] = BREAKUP_COMPLETE;
```

**Changes:** ~15 lines affected

---

### 5. **src/spray_drop_distort_NH3.c** (Main Breakup Logic)
**Purpose:** Replace all flag checks and assignments

#### Line 329: Tracked parcel check
**Current:**
```c
if (!tracking_initialized && old_parcel_cloud.is_child[p_idx] == 0 && 
    old_parcel_cloud.pbt[p_idx] == 1)
```

**New:**
```c
if (!tracking_initialized && 
    old_parcel_cloud.breakup_phase[p_idx] >= BREAKUP_ELIGIBLE &&
    old_parcel_cloud.breakup_phase[p_idx] < BREAKUP_COMPLETE)
```

#### Line 369: Temperature too high reset
**Current:**
```c
reset_parcel_to_child(&old_parcel_cloud, p_idx, "Temperature too high");
```

**New:**
```c
old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Line 402: P_sat too low reset
**Current:**
```c
reset_parcel_to_child(&old_parcel_cloud, p_idx, "Pre-check: P_sat < P_amb");
```

**New:**
```c
old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Line 426: Not superheated
**Current:**
```c
reset_parcel_to_child(&old_parcel_cloud, p_idx, "Not superheated");
```

**New:**
```c
old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Line 430-457: Stuck superheat fix (YOUR RECENT ADDITION)
**Current:**
```c
if (old_parcel_cloud.is_child[p_idx] == 0 && 
    P_sat_new >= P_amb &&
    old_parcel_cloud.lifetime[p_idx] > 1.0e-5 && 
    (old_parcel_cloud.pbt[p_idx] == 0 || 
     old_parcel_cloud.thermal_breakup_flag[p_idx] >= 0)) {
    old_parcel_cloud.pbt[p_idx] = 1;
    old_parcel_cloud.tbt[p_idx] = 0;
    old_parcel_cloud.thermal_breakup_flag[p_idx] = -1;
}
```

**New:**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == BREAKUP_DISABLED && 
    P_sat_new >= P_amb &&
    old_parcel_cloud.lifetime[p_idx] > 1.0e-5) {
    old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_ELIGIBLE;
}
```

#### Line 447-463: Stuck parcel fix (large parent)
**Current:**
```c
if (old_parcel_cloud.is_child[p_idx] == 0 && 
    old_parcel_cloud.radius[p_idx] > 70e-6 &&
    old_parcel_cloud.thermal_breakup_flag[p_idx] >= 0) {
    old_parcel_cloud.is_child[p_idx] = 1;
    old_parcel_cloud.film_flag[p_idx] = 1;
    old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;
}
```

**New:**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == BREAKUP_DISABLED && 
    old_parcel_cloud.radius[p_idx] > 70e-6) {
    old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_COMPLETE;
}
```

#### Line 469-481: Child reentry diagnostic
**Current:**
```c
if (old_parcel_cloud.is_child[p_idx] == 1 && 
    old_parcel_cloud.thermal_breakup_flag[p_idx] < 0 && 
    old_parcel_cloud.pbt[p_idx] == 1)
```

**New:**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == BREAKUP_COMPLETE)
    // This shouldn't happen - log it
```

#### Line 484-486: ENTRY TO THERMAL BREAKUP (CRITICAL)
**Current:**
```c
if (old_parcel_cloud.is_child[p_idx] == 0 &&
    old_parcel_cloud.thermal_breakup_flag[p_idx] < 0 && 
    old_parcel_cloud.pbt[p_idx] == 1)
```

**New:**
```c
if (old_parcel_cloud.breakup_phase[p_idx] >= BREAKUP_ELIGIBLE &&
    old_parcel_cloud.breakup_phase[p_idx] < BREAKUP_COMPLETE)
```

#### Line 657: Recovery flag check
**Current:**
```c
if (old_parcel_cloud.thermal_breakup_flag[p_idx] == 888)
```

**New:**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == BREAKUP_RECOVERY)
```

#### Line 681: Post-RPE v_bubble small
**Current:**
```c
reset_parcel_to_child(&old_parcel_cloud, p_idx, "Post-RPE: v_bubble too small");
```

**New:**
```c
old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Line 857: Breakup threshold reached
**Current:**
```c
old_parcel_cloud.thermal_breakup_flag[p_idx] = 3;
old_parcel_cloud.tbt[p_idx] = 1;
```

**New:**
```c
old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_READY;
```

#### Line 884-901: Stuck no thermal
**Current:**
```c
if (old_parcel_cloud.is_child[p_idx] == 0 && 
    old_parcel_cloud.radius[p_idx] > 70e-6) {
    old_parcel_cloud.is_child[p_idx] = 1;
    old_parcel_cloud.film_flag[p_idx] = 1;
    old_parcel_cloud.thermal_breakup_flag[p_idx] = 999;
    old_parcel_cloud.pbt[p_idx] = 0;
    old_parcel_cloud.tbt[p_idx] = 0;
}
```

**New:**
```c
if (old_parcel_cloud.breakup_phase[p_idx] < BREAKUP_COMPLETE && 
    old_parcel_cloud.radius[p_idx] > 70e-6) {
    old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_COMPLETE;
}
```

#### Line 906: Check tbt and thermal_breakup_flag
**Current:**
```c
if (old_parcel_cloud.tbt[p_idx] && old_parcel_cloud.thermal_breakup_flag[p_idx] != 4)
```

**New:**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == BREAKUP_READY)
```

**Changes:** ~30 locations in this file

---

### 6. **src/RPE_euler.c** (Thermal Breakup Solver)
**Purpose:** Replace flag checks and recovery logic

#### Line 254: Droplet too small
**Current:**
```c
reset_parcel_to_child(old_parcel_cloud, p_idx, "Droplet too small");
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Line 335: P_sat < P_amb (subcooled)
**Current:**
```c
reset_parcel_to_child(old_parcel_cloud, p_idx, "P_sat < P_amb (subcooled)");
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Line 371: Recovered parcel
**Current:**
```c
reset_parcel_to_child(old_parcel_cloud, p_idx, "Recovered parcel in RPE");
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Line 389: Bubble collapse
**Current:**
```c
reset_parcel_to_child(old_parcel_cloud, p_idx, "Bubble collapse (Rdot < 0)");
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_RECOVERY;
// Or BREAKUP_COMPLETE if not attempting recovery
```

#### Line 422: Subcooled (T < T_sat)
**Current:**
```c
reset_parcel_to_child(old_parcel_cloud, p_idx, "Subcooled (T < T_sat)");
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Line 456: Rdot too small
**Current:**
```c
reset_parcel_to_child(old_parcel_cloud, p_idx, "Rdot too small");
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

#### Recovery flag assignments (~line 888 equivalent)
**Current:**
```c
old_parcel_cloud->thermal_breakup_flag[p_idx] = 888;
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_RECOVERY;
```

**Changes:** ~15 locations in this file

---

### 7. **src/RPE_song.c** (Song Solver)
**Purpose:** Replace flag checks

#### Line ~9 comment: Update description
**Current:**
```c
// Void fraction = 0.55 triggers reset_parcel_to_child, after which DGRE takes over.
```

**New:**
```c
// Void fraction = 0.55 triggers breakup_phase = BREAKUP_COMPLETE.
```

#### Void fraction threshold reached
**Current:**
```c
reset_parcel_to_child(&old_parcel_cloud, p_idx, "Song: void fraction reached");
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_READY;
```

**Changes:** ~5 locations in this file

---

### 8. **src/Vb.c** (Bubble velocity calculation)
**Purpose:** Replace flag check

#### Line 22: P_sat check
**Current:**
```c
reset_parcel_to_child(old_parcel_cloud, p_idx, "Vb calc: P_sat < P_amb");
```

**New:**
```c
old_parcel_cloud->breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

**Changes:** 1-2 locations

---

### 9. **src/parcel_reset.c** (Reset utility function)
**Purpose:** Update or remove this function

#### Option A: Update function
**Current:**
```c
void reset_parcel_to_child(...) {
    parcel_cloud->is_child[p_idx] = 1;
    parcel_cloud->thermal_breakup_flag[p_idx] = 999;
    parcel_cloud->pbt[p_idx] = 0;
    // ...
}
```

**New:**
```c
void reset_parcel_to_child(...) {
    parcel_cloud->breakup_phase[p_idx] = BREAKUP_COMPLETE;
    // Rest stays same (r_bubble = 0, etc.)
}
```

#### Option B: Remove function entirely
Replace all calls with direct assignment:
```c
old_parcel_cloud.breakup_phase[p_idx] = BREAKUP_COMPLETE;
```

**Changes:** Update function or remove entirely (~10 call sites)

---

### 10. **src/Breakup.c** (Breakup execution)
**Purpose:** Set breakup_phase for children

When creating children, set:
```c
child_parcel_cloud.breakup_phase[child_idx] = BREAKUP_COMPLETE;
```

And parent either:
- Becomes BREAKUP_COMPLETE (if converted to first child)
- Or is deleted

**Changes:** ~5 locations

---

### 11. **include/parcel_reset.h**
**Purpose:** Update function signature or remove

```c
// Update signature
void reset_parcel_to_child(struct ParcelCloud* parcel_cloud, 
                           CONVERGE_index_t p_idx, 
                           const char* reason);

// Or remove if function is eliminated
```

**Changes:** 1 file

---

## Diagnostic/Output Files

### 12. **Any printf/logging statements**
Update diagnostics to use breakup_phase:

**Current:**
```c
printf("is_child=%d, pbt=%d, tbt=%d, tbf=%d", 
       is_child, pbt, tbt, thermal_breakup_flag);
```

**New:**
```c
printf("breakup_phase=%d", breakup_phase);
```

**Changes:** ~10-20 locations

---

## Summary Table

| File | Approx. Changes | Type |
|------|----------------|------|
| `include/env.h` | 1 | Add enum |
| `src/load_spray_env.c` | 12 | Register/load variable |
| `include/parcel_cloud.h` | 4 | Update struct |
| `src/parcel_prop.c` | 15 | Initialization |
| `src/spray_drop_distort_NH3.c` | 30 | Main logic |
| `src/RPE_euler.c` | 15 | Solver logic |
| `src/RPE_song.c` | 5 | Solver logic |
| `src/Vb.c` | 2 | Utility |
| `src/parcel_reset.c` | 10 | Utility function |
| `src/Breakup.c` | 5 | Child creation |
| `include/parcel_reset.h` | 1 | Header |
| Diagnostics/logging | 20 | Printf statements |
| **TOTAL** | **~120** | **changes** |

---

## Migration Strategy

### Phase 1: Add New System (Parallel)
1. Add `breakup_phase` to env.h, load_spray_env.c, parcel_cloud.h
2. Initialize `breakup_phase` alongside old flags in parcel_prop.c
3. Keep both systems working (redundant but safe)
4. Test: Compile and run, verify no crashes

### Phase 2: Migrate Code Section-by-Section
1. **Start with parcel_prop.c** (initialization only)
2. **Then spray_drop_distort_NH3.c** (main logic, critical)
3. **Then RPE_euler.c and RPE_song.c** (solvers)
4. **Then utilities** (Breakup.c, Vb.c, parcel_reset.c)
5. Test after each section

### Phase 3: Remove Old System
1. Comment out old flag registrations in load_spray_env.c
2. Comment out old struct members in parcel_cloud.h
3. Test: Should still compile and run
4. Remove commented code if successful

### Phase 4: Cleanup
1. Remove parcel_reset.c function (if not needed)
2. Clean up diagnostic messages
3. Update documentation

---

## Testing Checklist

After each phase:
- [ ] Code compiles without errors
- [ ] Code links without errors
- [ ] Simulation runs without crashes
- [ ] Parcels enter thermal breakup (check for `[THERMAL_ENTRY]` messages)
- [ ] Parcels break up correctly (children created)
- [ ] Children don't re-enter thermal breakup
- [ ] Stuck parcel fixes still work
- [ ] Recovery logic still works (if applicable)
- [ ] Output files (.h5) contain breakup_phase variable

---

## Rollback Plan

If implementation fails:

```bash
cd /home/apollo19/Desktop/Dan_B/UDF
git checkout Song_v1  # Return to backup branch
```

All changes will be reverted to pre-simplification state.

---

## Estimated Time

- Phase 1: 1 hour (add new system)
- Phase 2: 4-6 hours (migrate code, section by section)
- Phase 3: 1 hour (remove old system)
- Phase 4: 1 hour (cleanup)

**Total:** 7-9 hours for complete implementation and testing

---

## Risk Assessment

**Low Risk:**
- Adding new variable (doesn't break anything)
- Parallel implementation (both systems work)

**Medium Risk:**
- Entry logic changes (line 484 in spray_drop_distort_NH3.c)
- Recovery logic changes (may need careful testing)

**High Risk:**
- Removing old system before verification complete

**Mitigation:** Follow phased approach, test thoroughly at each step

---

## Next Steps

1. ✅ Backup branch created: `Song_v1`
2. Confirm you want to proceed
3. Start Phase 1: Add new variable system
4. Test compilation
5. Begin Phase 2 migration (section by section)

Ready to proceed?
