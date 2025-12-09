# Edge Cases and breakup_phase Assignments

## Summary of All Edge Cases in spray_drop_distort_NH3.c

### 1. **Temperature Too High (Line 379)**
```c
if (Td > temp_drop_0 + 2.0) {
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Temperature too high");
    old_parcel_cloud.num_drop[p_idx] = 0.0;  // Zero out for removal
    old_parcel_cloud.radius[p_idx] = 0.0;
}
```
**Result:** `breakup_phase = 5` (COMPLETE/child)
**Reason:** Droplet too hot - unphysical, remove from thermal breakup

---

### 2. **P_sat Too Low (Line 412)**
```c
if (P_sat < 1.0) {  // Below Antoine equation range
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Pre-check: P_sat < P_amb");
}
```
**Result:** `breakup_phase = 5` (COMPLETE/child)
**Reason:** Temperature outside Antoine range - disable thermal breakup

---

### 3. **Not Superheated (Line 436)**
```c
if (P_sat_new < P_amb) {  // Subcooled
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Not superheated");
}
```
**Result:** `breakup_phase = 5` (COMPLETE/child)
**Reason:** Parcel cooled below saturation - no longer superheated

---

### 4. **Stuck Superheated - Re-enable (Line 444)**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == 0 &&  // Disabled
    P_sat_new >= P_amb &&  // But IS superheated
    lifetime > 1.0e-5) {
    old_parcel_cloud.breakup_phase[p_idx] = 1;  // Re-enable
}
```
**Result:** `breakup_phase = 1` (ELIGIBLE)
**Reason:** Parcel was incorrectly disabled but is actually superheated - fix it

---

### 5. **Stuck Large Parent (Line 481)**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == 0 &&  // Disabled
    old_parcel_cloud.radius[p_idx] > 70e-6) {      // Large
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Stuck parcel fix");
}
```
**Result:** `breakup_phase = 5` (COMPLETE/child)
**Reason:** Large parent never entered thermal breakup - force to child so KH-RT/evaporation can work

---

### 6. **Child Trying to Re-enter (Line 499)**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == 5) {  // Is child
    // Log warning - this shouldn't happen
}
```
**Result:** No change (stays at 5)
**Reason:** Diagnostic - children should never try to re-enter thermal breakup

---

### 7. **Song Model: Void Fraction Reached (Line 628)**
```c
if (epsilon >= 0.55) {  // Song breakup criterion
    old_parcel_cloud.breakup_phase[p_idx] = 4;  // READY
    break;
}
```
**Result:** `breakup_phase = 4` (READY)
**Reason:** Song model void fraction threshold reached - ready to break

---

### 8. **Song Model: v_bubble Too Small (Line 663)**
```c
if (old_parcel_cloud.v_bubble[p_idx] < 1.0e-10) {  // Song edge case
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Song: v_bubble too small");
    break;
}
```
**Result:** `breakup_phase = 5` (COMPLETE/child)
**Reason:** Song model - bubble stopped growing

---

### 9. **Post-RPE: v_bubble Too Small (Line 678)**
```c
if (old_parcel_cloud.v_bubble[p_idx] < 1.0e-10) {  // Thermal model
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Post-RPE: v_bubble too small");
    break;
}
```
**Result:** `breakup_phase = 5` (COMPLETE/child)
**Reason:** Thermal model - bubble growth stopped

---

### 10. **In Recovery Mode (Line 683)**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == 3) {  // RECOVERY
    // Skip rest of timestep, try again next time
    break;
}
```
**Result:** `breakup_phase = 3` (RECOVERY - no change)
**Reason:** Bubble collapsed, attempting recovery - stay in recovery state

---

### 11. **Rb Negative (Line 706)** ⚠️ ERROR CASE
```c
if (Rb < 0.0) {
    Rb = 0.0;
    old_parcel_cloud.breakup_phase[p_idx] = 5;  // Mark as child
    printf("Rb negative after RPE solver\n");
    CONVERGE_mpi_abort();
}
```
**Result:** `breakup_phase = 5` (COMPLETE/child) then ABORT
**Reason:** Error condition - negative bubble radius is unphysical

---

### 12. **Rb > R_drop (Bubble Exceeds Droplet) (Line 714)**
```c
if (Rb > old_parcel_cloud.radius[p_idx]) {
    old_parcel_cloud.breakup_phase[p_idx] = 4;  // READY to break
    old_parcel_cloud.r_bubble[p_idx] = 0.8 * old_parcel_cloud.radius[p_idx];
    break;
}
```
**Result:** `breakup_phase = 4` (READY)
**Reason:** Bubble exceeded droplet radius - immediate breakup trigger
**Action:** Caps bubble at 80% of droplet, exits loop → triggers breakup at line 925

---

### 13. **kb > Threshold (Line 881)**
```c
if (kb > kb_threshold) {  // Breakup criterion
    old_parcel_cloud.breakup_phase[p_idx] = 4;  // READY
    old_parcel_cloud.r_bubble[p_idx] = Rb;
    old_parcel_cloud.eta_drop[p_idx] = kb;
    break;
}
```
**Result:** `breakup_phase = 4` (READY)
**Reason:** Normal thermal breakup - kb criterion met
**Action:** Exits sub-timestep loop → triggers breakup at line 925

---

### 14. **Stuck After Loop: Never Entered Thermal (Line 906)**
```c
if (old_parcel_cloud.breakup_phase[p_idx] < 5 &&  // Is parent
    old_parcel_cloud.radius[p_idx] > 70e-6) {      // Large
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Stuck: never entered thermal");
}
```
**Result:** `breakup_phase = 5` (COMPLETE/child)
**Reason:** Large parent exited loop without entering thermal breakup - force to child

---

### 15. **Breakup Execution (Line 925)**
```c
if (old_parcel_cloud.breakup_phase[p_idx] == 4) {  // READY to break
    // Execute Breakup() or Breakup_Song()
}
```
**Result:** After breakup, parent → children all have `breakup_phase = 5` (set in Breakup.c)
**Reason:** Parcel fragmentation

---

## State Machine Flow

```
INJECTION → 1 (ELIGIBLE)
             ↓
         CONDITIONS CHECK
             ↓
    ┌────────┴────────┐
    │                 │
  FAIL              PASS
    │                 │
    ↓                 ↓
 0 (DISABLED)    2 (ACTIVE)
    │                 │
    │         Sub-timestep Loop
    │                 │
    │           ┌─────┴─────┐
    │           │           │
    │       COLLAPSE    kb > threshold
    │           │           │
    │           ↓           ↓
    │      3 (RECOVERY) 4 (READY)
    │           │           │
    │           ↓           ↓
    │      (try next    BREAKUP
    │       timestep)   Execution
    │                       │
    ↓                       ↓
 5 (COMPLETE/CHILD) ← 5 (COMPLETE/CHILD)
```

---

## Critical Paths to breakup_phase = 4 (READY)

1. **Normal thermal:** kb > threshold (line 881)
2. **Song model:** epsilon ≥ 0.55 (line 628)
3. **Bubble too big:** Rb > R_drop (line 714)

All three exit the sub-timestep loop with phase = 4, then execute breakup at line 925.

---

## Critical Paths to breakup_phase = 5 (CHILD)

### Direct Assignment:
1. Temperature too high (line 379)
2. P_sat too low (line 412)
3. Not superheated (line 436)
4. Stuck large parent (line 481, 906)
5. v_bubble too small (line 663, 678)
6. Rb negative error (line 706)

### Via Breakup Execution:
7. After successful breakup (line 925 → Breakup.c)

---

## Recovery Path (breakup_phase = 3)

Set in **RPE_euler.c** when bubble collapses (Rdot < 0):
- Bubble shrinking detected
- `breakup_phase = 3` (RECOVERY)
- Next timestep: check at line 683
- Stays in recovery until bubble starts growing again
- Then returns to phase 2 (ACTIVE) or 1 (ELIGIBLE)

---

## Validation Checklist

### ✅ Parcels that should break (phase → 4):
- [x] kb > threshold (normal thermal)
- [x] Song void fraction ≥ 0.55
- [x] Bubble exceeds droplet radius

### ✅ Parcels that should be disabled (phase → 5):
- [x] Subcooled (not superheated)
- [x] Temperature out of range
- [x] Bubble growth stopped
- [x] Stuck large parents

### ✅ Parcels that should recover (phase → 3):
- [x] Bubble collapse detected (in RPE_euler.c)
- [x] Recovery attempt logic

### ✅ Parcels that should re-enable (phase 0 → 1):
- [x] Stuck superheated fix

---

## Questions/Concerns

### 1. Line 714: Is capping bubble at 0.8*R_drop correct?
**Current behavior:** When Rb > R_drop, set phase=4, cap bubble, trigger breakup
**Is this right?** Seems like a safety check to prevent unphysical state

### 2. Should phase=0 (DISABLED) ever transition back to phase=1?
**Current answer:** YES - line 444 does this for stuck superheated parcels
**Is this correct?** This catches parcels that were incorrectly disabled

### 3. Recovery state details
**What sets phase=3?** RPE_euler.c when Rdot < 0 (not migrated yet)
**What exits recovery?** Needs to be checked in RPE_euler.c migration

---

## Files Still Using Old Flags (Not Migrated)

1. **parcel_prop.c** - initialization
2. **RPE_euler.c** - recovery logic, sets phase=3
3. **RPE_song.c** - Song solver
4. **Breakup.c** - sets children to phase=5
5. **Vb.c** - bubble velocity checks

These will be migrated in next phase.
