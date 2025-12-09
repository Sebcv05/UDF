# Film Flag Hijack - UPDATED ✅

**Date:** 2025-12-09  
**File:** parcel_prop.c  
**Status:** Now uses breakup_phase

---

## What is the "Hijack"?

The `film_flag` field is a built-in CONVERGE field intended for tracking wall films. However, **we don't use the wall film model** in this simulation.

Instead, we **hijack** this field to store **breakup_phase** so it can be:
1. Visualized in Tecplot
2. Exported to .h5 files
3. Post-processed for analysis

This is safe because the film model is not active.

---

## Updated Implementation

### parcel_inject() - Line 188
```c
// DIAGNOSTIC: Hijack film_flag to store breakup_phase (safe if not using film model)
parcel_cloud.film_flag[passed_parcel_idx] = parcel_cloud.breakup_phase[passed_parcel_idx];
```

**New parcels:** `breakup_phase = 1` (ELIGIBLE) → `film_flag = 1`

### parcel_child() - Line 312
```c
// DIAGNOSTIC: Hijack film_flag to store breakup_phase (safe if not using film model)
parcel_cloud.film_flag[passed_child_parcel_idx] = parcel_cloud.breakup_phase[passed_child_parcel_idx];
```

**Child parcels:** `breakup_phase = 5` (COMPLETE) → `film_flag = 5`

---

## OLD Implementation (Before Migration)

### Before:
```c
// parcel_inject()
parcel_cloud.film_flag[passed_parcel_idx] = 0;  // 0 = parent

// parcel_child()
parcel_cloud.film_flag[passed_child_parcel_idx] = 1;  // 1 = child
```

This only gave binary information (parent/child), not the full breakup state.

---

## What You'll See in .h5 Files

**FILM_FLAG values now mean:**

| FILM_FLAG | breakup_phase | Meaning |
|-----------|---------------|---------|
| 0 | DISABLED | Parent, not eligible for breakup (subcooled/small) |
| 1 | ELIGIBLE | Parent, ready for thermal breakup |
| 2 | ACTIVE | Parent, bubble actively growing |
| 3 | RECOVERY | Parent, bubble collapsed (unused) |
| 4 | READY | Parent, ready to fragment |
| 5 | COMPLETE | Child parcel (post-breakup) |

---

## Analysis Results (Post-Fix)

After rerunning with the updated code, you should see a distribution like:

```
FILM_FLAG = 0: XX%  (disabled parents)
FILM_FLAG = 1: XX%  (eligible parents)
FILM_FLAG = 2: XX%  (active parents - growing bubble)
FILM_FLAG = 4: XX%  (ready to break)
FILM_FLAG = 5: XX%  (children)
```

**Before the fix**, you would only see:
- FILM_FLAG = 0: ~78-88% (all parents lumped together)
- FILM_FLAG = 1: ~12-22% (all children)

**After the fix**, you get full breakup state visibility!

---

## Why This Matters

### Debugging Stuck Parcels
You can now identify in Tecplot:
- **FILM_FLAG = 0**: Parcels that are disabled (may be stuck if they should be superheated)
- **FILM_FLAG = 1**: Parcels waiting to enter thermal breakup
- **FILM_FLAG = 2**: Parcels actively growing bubbles
- **FILM_FLAG = 4**: Parcels about to break (should disappear quickly)
- **FILM_FLAG = 5**: Children (should not re-enter breakup)

### Visualization
Color parcels by FILM_FLAG in Tecplot to see:
- Which parcels are stuck in DISABLED state
- How many are actively in thermal breakup
- The spatial distribution of breakup states

---

## Other Hijacked Fields

**film_thickness** is also hijacked:
```c
parcel_cloud.film_thickness[idx] = parcel_cloud.r_bubble[idx];
```

This allows visualizing bubble radius in Tecplot as "FILM_THICKNESS" field.

---

## Verification

After rerunning the simulation, check:
```bash
python3 analyze_film_flag.py
```

You should now see values 0-5, not just 0-1!

---

## Files Modified

1. **parcel_prop.c** - Updated hijack to use breakup_phase
2. **analyze_film_flag.py** - Script to analyze .h5 outputs

---

## Backward Compatibility

Old .h5 files (before this fix) will show:
- FILM_FLAG = 0 or 1 only (binary parent/child)

New .h5 files (after this fix) will show:
- FILM_FLAG = 0, 1, 2, 3, 4, or 5 (full breakup state)

---

**Status: HIJACK UPDATED ✅**  
**Rerun simulation to see new breakup_phase values in FILM_FLAG field**
