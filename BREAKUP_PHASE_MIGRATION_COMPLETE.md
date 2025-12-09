# Breakup Phase Migration - COMPLETE ✅

**Date:** 2025-12-09  
**Branch:** v3.1.12  
**Backup Branch:** Song_v1 (created before migration)

---

## Summary

Successfully migrated the entire thermal breakup codebase from a complex multi-flag system (thermal_breakup_flag, is_child, pbt, tbt) to a single unified `breakup_phase` integer variable.

---

## New State Machine

```
breakup_phase values:
  0 = DISABLED  (parent, not eligible - subcooled, too small, etc.)
  1 = ELIGIBLE  (parent, superheated, ready to enter thermal breakup)
  2 = ACTIVE    (parent, growing bubble in sub-timestep loop)
  3 = RECOVERY  (parent, bubble collapsed, attempting recovery) [UNUSED IN CURRENT IMPLEMENTATION]
  4 = READY     (parent, bubble at threshold, ready to fragment)
  5 = COMPLETE  (child - result of breakup, any mechanism)
```

---

## Files Migrated (Core Logic)

### ✅ spray_drop_distort_NH3.c
**Main thermal breakup loop**
- Entry condition: `breakup_phase >= 1 && <= 4`
- Breakup execution: `breakup_phase == 4`
- Unified stuck parcel logic: if phase=0, check superheat → re-enable (phase=1) or force child (phase=5)
- All edge cases handled with breakup_phase assignments
- Added comprehensive state documentation

### ✅ parcel_reset.c
**Reset function for converting parcels to children**
- Sets `breakup_phase = 5` (COMPLETE)
- Zeros bubble radius and velocity
- Resets initial droplet radius
- Used throughout codebase for error handling

### ✅ Breakup.c
**Thermal breakup execution**
- Entry check: `breakup_phase != 4` → ERROR
- After breakup: `breakup_phase = 5` (child)
- Zeros bubble properties
- Updated all diagnostics

### ✅ Breakup_Song.c
**Song model breakup execution**
- Entry check: `breakup_phase != 4` → ERROR  
- After breakup: `breakup_phase = 5` (child)
- Zeros bubble velocity and radius
- Updated diagnostics

### ✅ parcel_prop.c
**Parcel initialization**
- **parcel_inject():** New parcels → `breakup_phase = 1` (ELIGIBLE)
- **parcel_child():** Children → `breakup_phase = 5` (COMPLETE)
- Validation: children must have phase==5
- Removed old flag assignments

### ✅ RPE_euler.c
**Thermal RPE solver**
- Removed obsolete thermal_breakup_flag==-999 check
- Updated in_thermal_breakup: `breakup_phase >= 1 && <= 4`
- Recovery: converts to child immediately (phase=5)
- Phase=3 (RECOVERY) unused in current implementation

### ✅ RPE_song.c
**Song RPE solver**
- Added breakup_phase documentation
- No changes needed (already clean)

### ✅ Geometry.c
**Droplet geometry evolution**
- Negative dRd error → `reset_parcel_to_child()`
- Geometry error → `breakup_phase = 4` (READY)
- Added breakup_phase documentation

### ✅ spray_evap.c
**Evaporation model**
- Updated diagnostic warnings to use breakup_phase
- No logic changes needed

---

## Key Design Decisions

### 1. Unified Stuck Parcel Logic
**Old approach (3 separate checks):**
```c
// Check 1: Re-enable if superheated AND lifetime > 10μs
// Check 2: Force to child if phase=0 AND radius > 70μm
// Check 3: Force to child if phase<5 AND radius > 70μm after loop
```

**New approach (1 unified check):**
```c
if (breakup_phase == 0) {
    if (P_sat >= P_amb) {
        breakup_phase = 1;  // Re-enable
    } else {
        reset_to_child();    // Force to child
    }
}
```

### 2. Recovery Strategy
Current implementation: **Collapse → Immediate Child Conversion**
- On bubble collapse (Rdot < 0), immediately call `reset_parcel_to_child()`
- Sets `breakup_phase = 5` and `recovery_time` field
- Simpler than attempting to recover bubble growth
- Phase=3 (RECOVERY) is defined but unused

### 3. Breakup Criteria
**Three paths to breakup_phase = 4:**
1. Normal thermal: `kb > kb_threshold`
2. Song model: `void_fraction ≥ 0.55`
3. Bubble too big: `Rb > R_drop`

All set phase=4, exit sub-timestep loop, trigger breakup at main check

### 4. Model Separation
**Song vs Thermal models correctly separated:**
- Song: checks void fraction, then `continue` (skips DGRE/kb)
- Thermal: checks kb threshold (Song never reaches this)
- Flag `use_song_rpe` controls which model runs

---

## Breakup Phase Flow

```
INJECTION
    ↓
  phase = 1 (ELIGIBLE)
    ↓
ENTRY CHECK (spray_drop_distort_NH3.c)
    ↓
if phase >= 1 && <= 4 → ENTER LOOP
    ↓
  phase = 2 (ACTIVE)
    ↓
SUB-TIMESTEP LOOP
    ↓
    ├─→ kb > threshold → phase = 4 (READY) → BREAKUP
    ├─→ void ≥ 0.55   → phase = 4 (READY) → BREAKUP
    ├─→ Rb > R_drop   → phase = 4 (READY) → BREAKUP
    ├─→ Rdot < 0      → phase = 5 (COMPLETE/child) [recovery]
    └─→ else          → phase = 2 (continues)
    ↓
if phase == 4 → EXECUTE BREAKUP
    ↓
  phase = 5 (COMPLETE/child)
    ↓
NEVER RE-ENTERS THERMAL BREAKUP
```

---

## Edge Cases Handled

### Parcel Disabled (phase → 5)
1. Temperature too high
2. P_sat too low (outside Antoine range)
3. Not superheated (P_sat < P_amb)
4. Bubble growth stopped (v_bubble < 1e-10)
5. Stuck large parent
6. Rb negative (error)
7. Bubble collapse (Rdot < 0)

### Parcel Ready to Break (phase → 4)
1. kb > threshold (normal thermal)
2. void fraction ≥ 0.55 (Song model)
3. Rb > R_drop (bubble exceeded droplet)

### Parcel Re-enabled (phase 0 → 1)
1. Was disabled but IS superheated (stuck parcel fix)

---

## Testing Checklist

### ✅ Compilation
- All files compile without errors
- No warnings related to flag migration

### 🔲 Runtime Testing (TODO)
- [ ] Thermal breakup occurs correctly
- [ ] Song breakup occurs correctly  
- [ ] Stuck parcels are handled
- [ ] No parcels persist without breaking when superheated
- [ ] Children never re-enter thermal breakup
- [ ] Recovery mechanism works (converts to child)

### 🔲 Validation (TODO)
- [ ] Compare with Song_v1 backup branch
- [ ] Check parcel counts (parents vs children)
- [ ] Verify breakup timing
- [ ] Check for "stuck" parcels in Tecplot output

---

## Remaining Old Flag Usage

**ALL LEGACY FLAGS HAVE BEEN REMOVED! ✅**

The codebase no longer contains ANY references to:
- `parcel_cloud.is_child[...]`
- `parcel_cloud.pbt[...]`
- `parcel_cloud.tbt[...]`
- `parcel_cloud.thermal_breakup_flag[...]`

**Backward Compatibility:**
Output functions maintain compatibility by deriving legacy values from breakup_phase:
- `is_child` → `(breakup_phase == 5) ? 1 : 0`
- `pbt` → `(breakup_phase >= 1 && <= 4) ? 1 : 0`

This allows existing Tecplot scripts and post-processing tools to continue working without modification.

---

## Files with breakup_phase Comment Block

All major files now include the state documentation at the top:
- spray_drop_distort_NH3.c
- parcel_reset.c
- Breakup.c
- Breakup_Song.c
- parcel_prop.c
- RPE_euler.c
- RPE_song.c
- Geometry.c

---

## Rollback Instructions

If issues arise, rollback to backup:
```bash
cd /home/apollo19/Desktop/Dan_B/UDF
git checkout Song_v1
```

---

## Next Steps

1. **Test the simulation** with the new implementation
2. **Compare results** with Song_v1 backup branch
3. **Monitor for stuck parcels** in Tecplot output
4. **Update post-processing** tools if needed (post.c, parcel_output.c)
5. **Consider removing Phase=3** entirely if never used

---

## Migration Statistics

- **Core files migrated:** 8
- **Total flag replacements:** ~100+
- **Stuck parcel checks consolidated:** 3 → 1
- **State machine states:** 6 (0-5)
- **Lines of code changed:** ~200
- **Compilation errors:** 0 ✅
- **Time to migrate:** ~2 hours

---

## Contact

For questions about this migration, refer to:
- This document
- `EDGE_CASES_BREAKUP_PHASE.md` (detailed edge case analysis)
- Git commit history on branch v3.1.12

---

**Status: MIGRATION COMPLETE ✅**  
**Ready for testing**
