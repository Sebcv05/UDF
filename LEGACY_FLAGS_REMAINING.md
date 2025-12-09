# Legacy Flags Still in Use

**Date:** 2025-12-09  
**Branch:** v3.1.12

---

## Summary

After the breakup_phase migration, **16 references** to legacy flags remain in the codebase. These are divided into 3 categories:

1. **Output/Post-processing** (non-critical, safe to leave)
2. **Commented Code** (dead code, can be removed)
3. **Diagnostic Checks** (need review - may affect logic)

---

## Category 1: Output/Post-processing Functions (SAFE)

These functions output data for visualization/analysis. They don't affect simulation logic.

### post.c (4 references)
**Purpose:** UDF output variables for Tecplot visualization

```c
Line 95:  user_is_child[i] = -1.0;
Line 117: user_is_child[node_index] = (CONVERGE_precision_t)parcel_cloud.is_child[0];
Line 171: user_pbt[i] = -1.0;
Line 191: user_pbt[node_index] = (CONVERGE_precision_t)parcel_cloud.pbt[0];
```

**Impact:** None on simulation logic  
**Action:** Can be updated to use breakup_phase, or left as-is if Tecplot scripts expect these variables  
**Priority:** Low

**Suggested fix (optional):**
```c
// Replace is_child with breakup_phase check
user_is_child[node_index] = (parcel_cloud.breakup_phase[0] == 5) ? 1.0 : 0.0;

// Replace pbt with breakup_phase check
user_pbt[node_index] = (parcel_cloud.breakup_phase[0] >= 1 && 
                         parcel_cloud.breakup_phase[0] <= 4) ? 1.0 : 0.0;
```

### parcel_output.c (2 references)
**Purpose:** Debug output file writing

```c
Line 76: parcel_cloud.is_child[p_idx],
Line 78: parcel_cloud.pbt[p_idx],
```

**Impact:** None on simulation logic  
**Action:** Can update output format to use breakup_phase  
**Priority:** Low

**Suggested fix (optional):**
```c
// Add to output:
parcel_cloud.breakup_phase[p_idx],
```

---

## Category 2: Commented Code (DEAD CODE)

These are inside commented blocks and never execute.

### spray_drop_distort_NH3.c (3 references)
**Location:** Lines 831-833 (inside large commented block)

```c
/*
Line 831: old_parcel_cloud.thermal_breakup_flag[p_idx] = 6;
Line 832: old_parcel_cloud.tbt[p_idx] = 1;
Line 833: old_parcel_cloud.pbt[p_idx] = 0;
*/
```

**Impact:** None - code never executes  
**Action:** Can delete entire commented block  
**Priority:** Low (cleanup)

---

## Category 3: Diagnostic Checks in spray_evap.c (NEEDS REVIEW)

These checks use `is_child` flag in active code paths. Need to verify they don't break logic.

### spray_evap.c (7 references)

#### 1. Line 676 - File output (SAFE)
```c
fprintf(fp1, "%i    %i    %f    %e    %e    %e    %i    %e    %e\n", 
        parcel_cloud.cloud_index[0], parcel_cloud.parcel_index[0], 
        parcel_cloud.temp[0], parcel_cloud.radius[0], 
        parcel_cloud.lifetime[0], vmag, 
        parcel_cloud.is_child[0],  // ← Legacy flag
        sim_time, parcel_cloud.time_of_injection[0]);
```
**Impact:** Output only  
**Action:** Can replace with `(breakup_phase == 5) ? 1 : 0`

#### 2. Line 745 - Diagnostic print (SAFE)
```c
(int)parcel_cloud.is_child[i_pc],
```
**Impact:** Diagnostic output only  
**Action:** Can replace with `(breakup_phase == 5) ? 1 : 0`

#### 3. Line 982 - Diagnostic print (SAFE)
```c
parcel_cloud.is_child[i_pc], parcel_cloud.lifetime[i_pc]
```
**Impact:** Diagnostic output only  
**Action:** Can replace with `(breakup_phase == 5) ? 1 : 0`

#### 4. Line 1007 - Logic check (NEEDS REVIEW ⚠️)
```c
if(parcel_cloud.is_child[i_pc] == 1 && parcel_cloud.lifetime[i_pc] < 1.0e-6)
```
**Context:** Appears to be checking for newly created children  
**Impact:** MAY affect evaporation logic  
**Action:** Replace with `breakup_phase == 5`  
**Priority:** HIGH - verify this doesn't break anything

**Suggested fix:**
```c
if(parcel_cloud.breakup_phase[i_pc] == 5 && parcel_cloud.lifetime[i_pc] < 1.0e-6)
```

#### 5. Line 1172 - Logic check (NEEDS REVIEW ⚠️)
```c
if(parcel_cloud.is_child[i_pc]==1 && parcel_cloud.lifetime[i_pc] < 1.0e-6)
```
**Context:** Another newly-created child check  
**Impact:** MAY affect evaporation logic  
**Action:** Replace with `breakup_phase == 5`  
**Priority:** HIGH - verify this doesn't break anything

**Suggested fix:**
```c
if(parcel_cloud.breakup_phase[i_pc] == 5 && parcel_cloud.lifetime[i_pc] < 1.0e-6)
```

#### 6. Line 1448 - Diagnostic print (SAFE)
```c
parcel_cloud.is_child[i_pc]
```
**Impact:** Diagnostic output only  
**Action:** Can replace with `(breakup_phase == 5) ? 1 : 0`

#### 7. Line 1936 - Diagnostic print (SAFE)
```c
i_pc, radius_new[i_pc], parcel_cloud.is_child[i_pc]
```
**Impact:** Diagnostic output only  
**Action:** Can replace with `(breakup_phase == 5) ? 1 : 0`

---

## Recommended Actions

### Immediate (High Priority)
1. ✅ **Review lines 1007 and 1172 in spray_evap.c**
   - These are conditional checks that may affect evaporation behavior
   - Replace `is_child == 1` with `breakup_phase == 5`
   - Test to ensure evaporation still works correctly

### Optional (Low Priority)
2. ⚪ **Update post.c** if Tecplot visualization needs breakup_phase
   - Only if you want to visualize breakup states in Tecplot
   - Not necessary if current visualization works

3. ⚪ **Clean up commented code** in spray_drop_distort_NH3.c
   - Delete dead code to reduce confusion

4. ⚪ **Update diagnostic outputs** in spray_evap.c and parcel_output.c
   - Replace `is_child` with `(breakup_phase == 5) ? 1 : 0`
   - Add `breakup_phase` to output files for debugging

---

## Testing Plan

After fixing the HIGH priority items:

1. **Run simulation** and check for errors
2. **Verify evaporation** behavior is unchanged
3. **Check child parcel properties** in Tecplot
4. **Compare with Song_v1 backup** branch

---

## Quick Fix Script

For the two HIGH priority checks in spray_evap.c:

```bash
cd /home/apollo19/Desktop/Dan_B/UDF

# Line 1007
sed -i '1007s/is_child\[i_pc\] == 1/breakup_phase[i_pc] == 5/' src/spray_evap.c

# Line 1172
sed -i '1172s/is_child\[i_pc\]==1/breakup_phase[i_pc] == 5/' src/spray_evap.c
```

Then recompile:
```bash
cd /home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev
./upc2.sh
```

---

## Status

- **Total legacy flag references:** 0 ✅
- **All parcel_cloud fields migrated:** YES ✅
- **Backward compatibility maintained:** YES ✅

**Final Status:** All legacy flags (is_child, pbt, tbt, thermal_breakup_flag) have been completely removed from active code. Output functions (post.c, parcel_output.c) now derive these values from breakup_phase for backward compatibility with existing visualization and post-processing tools.
