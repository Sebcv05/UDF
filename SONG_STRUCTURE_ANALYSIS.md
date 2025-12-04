# Song Model Structure Analysis

## Current Problem

**Song model is currently embedded in the thermal breakup loop:**
- It shares the same sub-timestep loop
- After Song RPE, code still calls Geometry(), DGRE(), BreakupCriterion()
- Breakup is determined by kb > threshold, NOT void fraction = 0.55

**This is WRONG because:**
1. Song model already determines breakup (void = 0.55)
2. DGRE/kb criteria are redundant
3. Song calls `reset_parcel_to_child()` but then code continues through DGRE/kb
4. Two competing breakup criteria (void fraction vs kb)

---

## Two Implementation Options

### Option 1: Separate Song UDF Function
**Create:** `src/spray_drop_distort_song.c`

**Pros:**
- Clean separation of physics
- No interference between models
- Easier to maintain/debug
- Clear code flow for each model
- Can have different sub-timestep strategies

**Cons:**
- Code duplication (some common functions)
- Need to register new UDF
- More files to manage

**Structure:**
```c
CONVERGE_UDF(drop_distort_song, ...) {
    // For each parcel in thermal breakup:
    //   1. RPE_song_solver() in sub-timestep loop
    //   2. If void >= 0.55:
    //        - Call Breakup() directly
    //        - NO DGRE, NO kb criterion
    //   3. Update parcel states
}
```

---

### Option 2: Fork Inside spray_drop_distort_NH3.c
**Keep same file, but fork logic early**

**Pros:**
- Single UDF registration
- Shared initialization/cleanup
- Less code duplication
- Easier to switch models

**Cons:**
- More complex control flow
- Risk of accidentally sharing wrong logic
- Harder to follow for Song-specific behavior

**Structure:**
```c
static void spray_distort_cell_NH3(...) {
    // Common setup
    
    if (use_song_rpe) {
        // === SONG PATH ===
        // Sub-timestep loop:
        //   - RPE_song_solver()
        //   - Check void fraction
        //   - If void >= 0.55:
        //       * Call Breakup() directly
        //       * Skip DGRE/kb entirely
        
    } else {
        // === THERMAL PATH (current code) ===
        // Sub-timestep loop:
        //   - RPE_euler_solver()
        //   - Geometry()
        //   - DGRE()
        //   - BreakupCriterion()
        //   - If kb > threshold: Breakup()
    }
}
```

---

## Recommended Approach: **Option 2 (Fork Inside)**

**Rationale:**
1. Minimal disruption to existing code
2. Easier to compare thermal vs Song in same file
3. Can share common functions (Breakup, property lookups)
4. Single UDF means simpler build system

**Implementation plan:**

### Step 1: Identify Fork Point
Fork after line 450 (where thermal_breakup_flag check happens)

### Step 2: Song-Specific Loop
```c
if (use_song_rpe) {
    // SONG MODEL: Simple void fraction criterion
    
    // Sub-timestep loop
    for (int sub = 0; sub < n_sub; ++sub) {
        // Grow bubble
        RPE_song_solver(...);
        
        // Check void fraction
        CONVERGE_precision_t epsilon = pow(R_bubble/R_drop_0, 3);
        
        if (epsilon >= 0.55) {
            // Breakup immediately - NO DGRE, NO kb
            old_parcel_cloud.thermal_breakup_flag[p_idx] = 3;
            old_parcel_cloud.tbt[p_idx] = 1;
            
            // Call Breakup() directly
            Breakup(&old_parcel_cloud, p_idx, ...);
            break;
        }
    }
    continue; // Skip thermal path
}

// THERMAL MODEL: Current complex logic
// ... existing code ...
```

### Step 3: Key Differences

| Aspect | Thermal | Song |
|--------|---------|------|
| RPE solver | `RPE_euler_solver()` | `RPE_song_solver()` |
| Geometry() | ✅ Yes | ❌ No |
| DGRE() | ✅ Yes | ❌ No |
| BreakupCriterion() | ✅ Yes (kb) | ❌ No |
| Breakup trigger | kb > threshold | void >= 0.55 |
| Temperature | Evolves | Constant |

---

## What Needs to Change

### 1. Remove `reset_parcel_to_child()` from RPE_song.c
Currently Song solver calls this, but breakup should happen in main loop

### 2. Add Void Fraction Check in Main Loop
After RPE_song_solver(), check epsilon directly

### 3. Call Breakup() Directly
Don't go through DGRE/kb path

### 4. Skip Geometry()?
Song model is isothermal - does droplet expand?
- If yes: Need Geometry()
- If no: Skip it

**Decision needed:** Does Song model need Geometry() call?

---

## Code Changes Required

### File: `src/RPE_song.c`
**REMOVE these lines (~182):**
```c
reset_parcel_to_child(old_parcel_cloud, p_idx, "Void=0.55 (Song)");
return;
```

**REPLACE with:**
```c
// Don't reset here - let main loop handle breakup
// Just return, main loop will check void fraction
return;
```

### File: `src/spray_drop_distort_NH3.c`
**ADD after line 516 (after RPE solver call):**
```c
if (use_song_rpe) {
    // Song model: check void fraction for breakup
    CONVERGE_precision_t R = old_parcel_cloud.r_bubble[p_idx];
    CONVERGE_precision_t Ro = old_parcel_cloud.r_drop_0[p_idx];
    CONVERGE_precision_t epsilon = (R*R*R) / (Ro*Ro*Ro);
    
    if (epsilon >= 0.55) {
        printf("[SONG_BREAKUP] Void=%.4f >= 0.55, triggering breakup\n", epsilon);
        old_parcel_cloud.thermal_breakup_flag[p_idx] = 3;
        old_parcel_cloud.tbt[p_idx] = 1;
        // Skip DGRE/kb, jump to breakup
        break;  // Exit sub-timestep loop
    }
    
    // Continue sub-timestep loop (void < 0.55)
    continue;  // Skip Geometry/DGRE/kb for Song
}

// THERMAL MODEL continues with existing code
// (Geometry, DGRE, BreakupCriterion, kb check)
```

---

## Testing Plan

### Test 1: Verify Song Skips DGRE
```bash
grep "DGRE" converge.log  # Should NOT appear with Song model
```

### Test 2: Verify Void Fraction Breakup
```bash
grep "SONG_BREAKUP" converge.log  # Should show epsilon values ~0.55
```

### Test 3: Compare Breakup Statistics
- Run thermal: Count thermal breakups
- Run Song: Count Song breakups
- Should be similar numbers (validating Song works)

### Test 4: Verify No kb Criterion
```bash
grep "kb" converge.log  # Should NOT appear with Song model
```

---

## Questions to Resolve

1. **Does Song need Geometry()?**
   - Isothermal → no thermal expansion
   - But bubble growth changes droplet volume
   - Probably YES, need Geometry()

2. **What about DGRE after breakup?**
   - Once parcel breaks up, it becomes children
   - Children follow normal KH-RT or other models
   - DGRE not relevant after breakup

3. **Should Song have its own breakup flag?**
   - Currently uses thermal_breakup_flag = 3
   - Could create song_breakup_flag for clarity
   - Or keep same flag (simpler)

---

## Summary

**RECOMMENDED:** Implement Option 2 (fork inside spray_drop_distort_NH3.c)

**Changes needed:**
1. Remove `reset_parcel_to_child()` from RPE_song.c
2. Add void fraction check after RPE_song_solver() call
3. Skip Geometry/DGRE/kb path for Song
4. Jump directly to Breakup() when void >= 0.55

**Benefits:**
- Song model operates independently of DGRE/kb
- Clean separation of physics
- Easier to validate
- Minimal code disruption

**This makes Song model work as intended: pure void fraction criterion.**
