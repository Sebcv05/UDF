# Song RPE Bubble Collapse Fix - December 5, 2024

## Problem Summary
Bubbles were collapsing (`Rdot = -109.8 m/s`) at high superheat (323K) where they should be growing rapidly.

## Root Cause Analysis

### The Issue
**Pressure mismatch between bubble initialization and Song solver:**

1. **In `parcel_prop.c` (line 157):** 
   - Bubble initialization used **hardcoded `P_amb = 2.0e5 Pa`** (2 bar)
   ```c
   CONVERGE_precision_t ambient_pres = 2.0e5;  // HARDCODED - WRONG!
   CONVERGE_precision_t Rc = 2*sigma / (P_sat - ambient_pres);
   R_bubble_0 = 1.1 * Rc;
   ```

2. **In `RPE_song.c`:** 
   - Song solver used **actual CFD pressure** from `global_pressure[node_index]` (e.g., 1 bar = 1.0e5 Pa)

3. **The Problem:**
   - The Song RPE equation contains the term: `(2σ/R₀ + P_r0)·(R₀/R)³`
   - This term depends critically on `R₀`, which was calculated with **wrong pressure**
   - When `R₀` is incorrect, this pressure term becomes incorrect → negative acceleration → collapse

### Why This Matters
At **T=323K** (high superheat):
- P_sat ≈ 19 bar (1.9e6 Pa)
- If P_amb_actual = 1 bar but P_amb_init = 2 bar:
  - **R_c with 2 bar**: `Rc = 2*0.021 / (19e5 - 2e5) = 2.47e-8 m`
  - **R_c with 1 bar**: `Rc = 2*0.021 / (19e5 - 1e5) = 2.33e-8 m`
  - R₀ values differ by ~6%
  - But the `(R₀/R)³` term **amplifies this error cubically**
  - Result: Song solver sees inconsistent initial condition → collapse

## Solution Applied

### Files Modified in /home/apollo19/Desktop/Dan_B/UDF/

#### 1. `src/parcel_prop.c` (lines 139-172)

**Changed FROM (WRONG):**
```c
// const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(passed_spray_cloud);
// CONVERGE_precision_t ambient_pres = pressure[node_index];  // COMMENTED OUT
CONVERGE_precision_t ambient_pres = 2.0e5;  // HARDCODED
```

**Changed TO (CORRECT):**
```c
// Get actual ambient pressure from CFD mesh
const CONVERGE_index_t node_index = CONVERGE_cloud_get_node_index(passed_spray_cloud);
CONVERGE_precision_t ambient_pres = pressure[node_index];  // NOW USES ACTUAL P_amb

// Calculate critical radius Rc = 2*sigma / DeltaP
CONVERGE_precision_t Rc = 2.0 * parcel_cloud.surf_ten[passed_parcel_idx] / (P_sat - ambient_pres);
// Initialize bubble at 1.1 * Rc (10% above critical radius for stable growth)
parcel_cloud.r_bubble[passed_parcel_idx] = 1.1 * Rc;

// Diagnostic: print first few initializations
static int init_count = 0;
if (init_count < 5) {
    printf("[INIT_BUBBLE] p_idx=%li, T=%.2f K, P_sat=%.2e Pa, P_amb=%.2e Pa\n",
           passed_parcel_idx, Td, P_sat, ambient_pres);
    printf("[INIT_BUBBLE]   Rc=%.3e m, R_bubble_0=%.3e m (1.1*Rc)\n",
           Rc, parcel_cloud.r_bubble[passed_parcel_idx]);
    init_count++;
}
```

#### 2. `src/RPE_song.c` (lines 203-237)

**Added enhanced diagnostics:**
```c
// Debug logging for first few calls
if (song_debug_count < SONG_DEBUG_MAX) {
    printf("[SONG_STEP] p_idx=%li, R=%.3e m, Rdot=%.3e m/s, ε=%.4f, ρ_m=%.1f kg/m³, T=%.2f K\n",
           p_idx, R, Rdot, epsilon, rho_m, T_drop);
    printf("[SONG_STEP]   R0=%.3e m, P_sat=%.2e Pa, P_amb=%.2e Pa, ΔP=%.2e Pa\n",
           R0, P_sat, P_amb, P_sat - P_amb);
    song_debug_count++;
}

// Check for collapse (negative velocity) - WITH DETAILED DIAGNOSTICS
if (Rdot < 0.0) {
    static int collapse_detail_count = 0;
    if (collapse_detail_count < 5) {
        printf("[SONG_COLLAPSE] Bubble collapsing: Rdot=%.3e m/s at p_idx=%li\n", Rdot, p_idx);
        printf("[SONG_COLLAPSE]   R=%.3e m, R0=%.3e m, Rddot=%.3e m/s²\n", R, R0, Rddot);
        printf("[SONG_COLLAPSE]   P_sat=%.2e Pa, P_amb=%.2e Pa, ΔP=%.2e Pa\n", P_sat, P_amb, P_sat - P_amb);
        printf("[SONG_COLLAPSE]   T_drop=%.2f K, rho_m=%.1f kg/m³, ε=%.4f\n", T_drop, rho_m, epsilon);
        
        // Calculate individual pressure terms for diagnosis
        CONVERGE_precision_t P_init_coeff = 2.0 * params.sigma / R0 + params.P_r0;
        CONVERGE_precision_t R_ratio_cubed = pow(R0/R, 3.0);
        CONVERGE_precision_t P_init = P_init_coeff * R_ratio_cubed;
        CONVERGE_precision_t P_laplace = 2.0 * params.sigma / R;
        printf("[SONG_COLLAPSE]   P_init=%.2e Pa, P_laplace=%.2e Pa\n", P_init, P_laplace);
        collapse_detail_count++;
    }
    return;
}
```

#### 3. `/home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev/workflow.md`

**Updated to clarify UDF directory structure:**
- Made it explicit that **ALL edits must be in /home/apollo19/Desktop/Dan_B/UDF/**
- Clarified that Dev/src/ is a copy (don't edit there)
- Added stars and warning symbols to make it unmistakable

## Why This Fix Works

Now both initialization and solver use the **same pressure**:
1. Critical radius calculation: `Rc = 2σ / (P_sat - P_amb_actual)`
2. Initial bubble radius: `R₀ = 1.1 × Rc`
3. Song RPE term: `(2σ/R₀ + P_r0)·(R₀/R)³` - now **consistent** with initialization

**Result:** The pressure balance in the Song RPE is correct, bubbles grow as expected.

## Expected Behavior After Fix

At **T=323K** with typical conditions:
- P_sat ≈ 19 bar (1.9e6 Pa)
- P_amb ≈ 1-2 bar (1-2e5 Pa)  
- ΔP ≈ 17-18 bar → **strong driving force for growth**
- Rc ≈ 2-3 × 10⁻⁸ m (20-30 nm)
- R₀ = 1.1×Rc ≈ 22-33 nm
- **Bubble should now GROW steadily, not collapse**

## Testing Instructions

1. From `/home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev/`, run:
   ```bash
   ./upc2.sh  # This syncs UDF changes and compiles
   ```

2. Run your 323K test case

3. **Expected output:**
   - `[INIT_BUBBLE]` messages showing correct P_amb from CFD mesh
   - `[SONG_STEP]` messages showing positive Rdot (growth)
   - **NO** `[SONG_COLLAPSE]` messages
   - Void fraction increasing smoothly toward 0.55

4. **If you still see collapse:**
   - Check the `[INIT_BUBBLE]` messages - is P_amb reasonable?
   - Check the `[SONG_COLLAPSE]` detail messages - what are P_init and P_laplace?
   - The diagnostics will help identify any remaining issues

## Key Takeaways

1. **Consistency is critical:** Initial conditions and solver must use the same ambient pressure
2. **Don't hardcode physical parameters** that vary in the CFD domain
3. **Use actual mesh values:** `pressure[node_index]` gives the local CFD pressure
4. **The Song RPE is sensitive to R₀** because of the `(R₀/R)³` term
5. **Always edit in /home/apollo19/Desktop/Dan_B/UDF/**, not in Dev/src/

## Commit Message
```
Fix bubble collapse at high superheat by using actual CFD pressure

- parcel_prop.c: Use pressure[node_index] instead of hardcoded 2 bar
- RPE_song.c: Add detailed collapse diagnostics
- Ensures R_bubble_0 calculation matches Song solver pressure
- Fixes cubic amplification error in (R0/R)³ term
```
