# Song RPE CONVERGE UDF Integration Plan

## Overview
This document outlines how to integrate the standalone Song RPE solver (`test_song_rpe.c`) into the CONVERGE UDF framework.

---

## Key Differences: Standalone vs CONVERGE

### Standalone Test (`test_song_rpe.c`)
```c
ResultArrays* solve_song_rpe(const Condition* cond) {
    // Read from test condition struct
    params.rho_l = RHO_L;  // Hardcoded constant
    params.T0 = cond->T0;  // From test condition
    
    // Initialize state
    state.R = R_bubble_0;   // Calculated from R_c
    state.Rdot = 0.0;       // Start from rest
    
    // Time loop
    while (t < TMAX) {
        // Compute acceleration
        // Euler step
        // Store to CSV array
    }
}
```

### CONVERGE UDF (to be created: `RPE_song.c`)
```c
void RPE_song_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
) {
    // Read from parcel_cloud[p_idx]
    params.rho_l = old_parcel_cloud->density[p_idx];
    params.T_drop = old_parcel_cloud->temp[p_idx];
    
    // Read current state
    state.R = old_parcel_cloud->r_bubble[p_idx];
    state.Rdot = old_parcel_cloud->v_bubble[p_idx];
    state.R0 = old_parcel_cloud->r_bubble_0[p_idx];
    state.Ro = old_parcel_cloud->r_drop_0[p_idx];
    
    // Single timestep integration (dt_sub provided by CONVERGE)
    // Compute acceleration
    // Euler step
    
    // Write back to parcel_cloud[p_idx]
    old_parcel_cloud->r_bubble[p_idx] = state.R;
    old_parcel_cloud->v_bubble[p_idx] = state.Rdot;
    // Temperature unchanged (isothermal)
}
```

---

## Files to Create/Modify

### 1. Create: `src/RPE_song.c`
**Purpose:** Song RPE solver for CONVERGE

**Key Functions:**
```c
// Main solver (called once per parcel per substep)
void RPE_song_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
);

// Helper functions (can copy from test_song_rpe.c):
CONVERGE_precision_t compute_void_fraction(
    CONVERGE_precision_t R, 
    CONVERGE_precision_t Ro
);

CONVERGE_precision_t compute_mixture_density(
    CONVERGE_precision_t epsilon,
    CONVERGE_precision_t rho_v,
    CONVERGE_precision_t rho_l
);

CONVERGE_precision_t compute_song_acceleration(
    CONVERGE_precision_t R,
    CONVERGE_precision_t Rdot,
    CONVERGE_precision_t R0,
    CONVERGE_precision_t rho_m,
    CONVERGE_precision_t P_sat,
    const SongParams* params
);
```

**Data Flow:**
```
INPUT (from parcel_cloud[p_idx]):
  - radius[p_idx]           → Current droplet radius
  - r_bubble[p_idx]         → Current bubble radius (R)
  - v_bubble[p_idx]         → Bubble wall velocity (Rdot)
  - r_bubble_0[p_idx]       → Initial bubble radius (R0)
  - r_drop_0[p_idx]         → Initial droplet radius (Ro)
  - temp[p_idx]             → Droplet temperature (isothermal)
  - density[p_idx]          → Liquid density (ρ_l)
  - viscosity[p_idx]        → Liquid viscosity (μ_l)
  - surf_ten[p_idx]         → Surface tension (σ)

OUTPUT (to parcel_cloud[p_idx]):
  - r_bubble[p_idx]         → Updated bubble radius
  - v_bubble[p_idx]         → Updated bubble velocity
  - temp[p_idx]             → UNCHANGED (isothermal)
```

### 2. Create: `include/RPE_song.h`
**Purpose:** Header file for Song RPE

```c
#ifndef RPE_SONG_H
#define RPE_SONG_H

#include "lagrangian/env.h"
#include <CONVERGE/udf.h>

// Song-specific parameters
typedef struct {
    CONVERGE_precision_t rho_l;       // Liquid density
    CONVERGE_precision_t mu_l;        // Dynamic viscosity
    CONVERGE_precision_t sigma;       // Surface tension
    CONVERGE_precision_t R_spec;      // Specific gas constant (488.2 for NH3)
    CONVERGE_precision_t P_amb;       // Ambient pressure
    CONVERGE_precision_t P_r0;        // Residual pressure (1e6 Pa)
    CONVERGE_precision_t kappa;       // Surface viscosity (0.0)
} SongParams;

// Function declarations
void RPE_song_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
);

CONVERGE_precision_t compute_void_fraction(
    CONVERGE_precision_t R, 
    CONVERGE_precision_t Ro
);

CONVERGE_precision_t compute_mixture_density(
    CONVERGE_precision_t epsilon,
    CONVERGE_precision_t rho_v,
    CONVERGE_precision_t rho_l
);

CONVERGE_precision_t compute_song_acceleration(
    CONVERGE_precision_t R,
    CONVERGE_precision_t Rdot,
    CONVERGE_precision_t R0,
    CONVERGE_precision_t rho_m,
    CONVERGE_precision_t P_sat,
    const SongParams* params
);

#endif // RPE_SONG_H
```

### 3. Modify: `src/spray_drop_distort_NH3.c`
**Purpose:** Add call to Song solver with switching logic

**Location:** Around line 505 where `RPE_euler_solver` is called

**Current code:**
```c
// Line ~505
RPE_euler_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                 hvap_table, cp_table, num_parcel_species);
```

**Modified code:**
```c
#include <RPE_song.h>  // Add to top of file

// Around line 505, replace with:
if (USE_SONG_RPE) {
    // Song isothermal model
    RPE_song_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                    hvap_table, cp_table, num_parcel_species);
} else {
    // Current thermal model
    RPE_euler_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                     hvap_table, cp_table, num_parcel_species);
}
```

### 4. Modify: `src/read_input.c` (or wherever user inputs are read)
**Purpose:** Add flag to switch between models

Add to user_inputs.in:
```
USE_SONG_RPE = 0    # 0=thermal RPE, 1=Song isothermal RPE
```

Read in code:
```c
// In read_input.c or globals.c
int USE_SONG_RPE = 0;  // Global variable

// In input reading function:
USE_SONG_RPE = CONVERGE_input_int("USE_SONG_RPE", 0);
printf("Using %s RPE model\n", USE_SONG_RPE ? "Song isothermal" : "thermal");
```

### 5. Modify: `include/globals.h`
**Purpose:** Declare global flag

```c
// Add to globals.h
extern int USE_SONG_RPE;
```

---

## Implementation Steps

### Step 1: Create Header File
```bash
cd /home/apollo19/Desktop/Dan_B/UDF
# Create include/RPE_song.h (see structure above)
```

### Step 2: Create Source File
**Strategy:** Copy from `test_song_rpe.c` and adapt

**Mapping:**
| Standalone Function | UDF Function | Changes |
|---------------------|--------------|---------|
| `solve_song_rpe()` | `RPE_song_solver()` | Remove time loop, read from parcel, single timestep |
| `compute_song_acceleration()` | Same | Keep as-is (just change types) |
| `compute_void_fraction()` | Same | Keep as-is |
| `compute_mixture_density()` | Same | Keep as-is |
| `Psat_from_Antoine()` | Use existing `Saturation_PressureNH3()` | Already in UDF |

**Key Changes:**
```c
// STANDALONE (test_song_rpe.c):
double compute_song_acceleration(
    const BubbleState* state,      // ← struct with R, Rdot, R0
    const SongParams* params,
    double rho_m,
    double P_sat
)

// UDF (RPE_song.c):
CONVERGE_precision_t compute_song_acceleration(
    CONVERGE_precision_t R,         // ← individual parameters
    CONVERGE_precision_t Rdot,
    CONVERGE_precision_t R0,
    CONVERGE_precision_t rho_m,
    CONVERGE_precision_t P_sat,
    const SongParams* params        // ← still have params struct
)
```

### Step 3: Main Solver Structure
```c
void RPE_song_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    CONVERGE_table_t** hvap_table,
    CONVERGE_table_t** cp_table,
    CONVERGE_size_t num_parcel_species
) {
    // 1. INITIALIZE PARAMETERS from parcel
    SongParams params;
    params.rho_l = old_parcel_cloud->density[p_idx];
    params.mu_l = old_parcel_cloud->viscosity[p_idx];
    params.sigma = old_parcel_cloud->surf_ten[p_idx];
    params.R_spec = 488.2;  // NH3 specific gas constant
    params.P_amb = P_amb;
    params.P_r0 = 1.0e6;    // Residual pressure
    params.kappa = 0.0;     // Surface viscosity
    
    // 2. READ STATE from parcel
    CONVERGE_precision_t R = old_parcel_cloud->r_bubble[p_idx];
    CONVERGE_precision_t Rdot = old_parcel_cloud->v_bubble[p_idx];
    CONVERGE_precision_t R0 = old_parcel_cloud->r_bubble_0[p_idx];
    CONVERGE_precision_t Ro = old_parcel_cloud->r_drop_0[p_idx];
    CONVERGE_precision_t T_drop = old_parcel_cloud->temp[p_idx];
    
    // 3. CALCULATE P_sat (isothermal - T_drop constant)
    CONVERGE_precision_t P_sat;
    Saturation_PressureNH3(T_drop, &P_sat);
    
    // 4. CHECK SUPERHEAT (same as thermal model)
    if (P_sat <= P_amb) {
        reset_parcel_to_child(old_parcel_cloud, p_idx, "Not superheated (Song)");
        return;
    }
    
    // 5. COMPUTE VOID FRACTION
    CONVERGE_precision_t epsilon = compute_void_fraction(R, Ro);
    
    // 6. CHECK TERMINATION (void = 0.99)
    if (epsilon >= 0.99) {
        reset_parcel_to_child(old_parcel_cloud, p_idx, "Void fraction = 0.99 (Song)");
        return;
    }
    
    // 7. COMPUTE VAPOR DENSITY (ideal gas)
    CONVERGE_precision_t rho_v = P_sat / (params.R_spec * T_drop);
    
    // 8. COMPUTE MIXTURE DENSITY
    CONVERGE_precision_t rho_m = compute_mixture_density(epsilon, rho_v, params.rho_l);
    
    // 9. COMPUTE ACCELERATION
    CONVERGE_precision_t Rddot = compute_song_acceleration(
        R, Rdot, R0, rho_m, P_sat, &params
    );
    
    // 10. EXPLICIT EULER STEP (single timestep)
    Rdot += Rddot * dt_sub;
    R += Rdot * dt_sub;
    
    // 11. SAFETY CHECKS
    if (R < 1e-12) R = 1e-12;
    if (Rdot < 0.0) {
        reset_parcel_to_child(old_parcel_cloud, p_idx, "Bubble collapse (Song)");
        return;
    }
    if (R > 0.999 * Ro) {
        R = 0.999 * Ro;
        Rdot = 0.0;
    }
    
    // 12. UPDATE PARCEL (R and Rdot only, T unchanged)
    old_parcel_cloud->r_bubble[p_idx] = R;
    old_parcel_cloud->v_bubble[p_idx] = Rdot;
    // Temperature NOT updated (isothermal model)
}
```

### Step 4: Integration with spray_drop_distort_NH3.c
**Location:** Line 505

```c
// Before (current):
RPE_euler_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                 hvap_table, cp_table, num_parcel_species);

// After (with switch):
if (USE_SONG_RPE) {
    RPE_song_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                    hvap_table, cp_table, num_parcel_species);
} else {
    RPE_euler_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                     hvap_table, cp_table, num_parcel_species);
}
```

### Step 5: Add to CMakeLists.txt
Ensure RPE_song.c is compiled:
```cmake
# Should be automatically picked up if in src/ directory
# Verify with: ls src/RPE_song.c
```

---

## Testing Strategy

### Phase 1: Compilation Test
```bash
cd /home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev
./upc2.sh  # Compile with new Song solver
# Check for errors in compile_log
```

### Phase 2: Simple Test Case
**Setup:**
- Single parcel injection
- Known superheat (e.g., T=293K, P=2bar)
- USE_SONG_RPE = 1
- Short simulation time

**Expected:**
- Bubble grows monotonically
- No crashes
- Parcel eventually reaches void=0.99 or converts to child

### Phase 3: Comparison Test
**Setup:**
- Run same case with USE_SONG_RPE = 0 (thermal)
- Run same case with USE_SONG_RPE = 1 (Song)
- Compare outputs

**Expected differences:**
- Song: Faster growth (no thermal limiting)
- Song: Temperature constant
- Song: Stops at void=0.99
- Thermal: Slower growth (thermal limiting)
- Thermal: Temperature decreases
- Thermal: Multiple stopping criteria

### Phase 4: Full Spray Test
- Multiple parcels
- Realistic injection conditions
- Monitor for stability

---

## Debugging Tips

### 1. Add Diagnostic Logging
```c
// In RPE_song_solver, add:
static int call_count = 0;
if (call_count < 10) {
    printf("[SONG_RPE] p_idx=%li, R=%.3e, Rdot=%.3e, epsilon=%.4f, T=%.2f K\n",
           p_idx, R, Rdot, epsilon, T_drop);
    call_count++;
}
```

### 2. Check for NaN/Inf
```c
if (isnan(R) || isinf(R) || isnan(Rdot) || isinf(Rdot)) {
    printf("[SONG_ERROR] NaN/Inf detected! p_idx=%li\n", p_idx);
    reset_parcel_to_child(old_parcel_cloud, p_idx, "NaN/Inf");
    return;
}
```

### 3. Monitor Void Fraction
```c
if (epsilon > 0.95) {
    printf("[SONG_NEAR_TARGET] p_idx=%li, epsilon=%.4f, R=%.3e μm\n",
           p_idx, epsilon, R*1e6);
}
```

---

## Key Points

### ✅ Advantages of Song Model
1. **Simpler:** Only 2 ODEs vs 4
2. **Faster:** No thermal calculations
3. **Physically motivated:** Based on Song et al. experiments
4. **Clear stopping:** void = 0.99

### ⚠️ Considerations
1. **Isothermal assumption:** May not be valid for all conditions
2. **No cooling:** Droplet stays at injection temperature
3. **Different physics:** Pressure-driven vs thermal-limited
4. **Needs R0 and Ro:** Initial bubble/droplet radii must be stored

### 🔧 Required Parcel Variables (already exist!)
- ✅ `r_bubble[p_idx]` - Current bubble radius
- ✅ `v_bubble[p_idx]` - Wall velocity
- ✅ `r_bubble_0[p_idx]` - Initial bubble radius (R0)
- ✅ `r_drop_0[p_idx]` - Initial droplet radius (Ro)
- ✅ `temp[p_idx]` - Temperature (read-only for Song)
- ✅ `density[p_idx]`, `viscosity[p_idx]`, `surf_ten[p_idx]`

---

## Next Steps

1. **Create RPE_song.h** with structures and function declarations
2. **Create RPE_song.c** by adapting test_song_rpe.c
3. **Add USE_SONG_RPE flag** to user inputs
4. **Modify spray_drop_distort_NH3.c** to call Song solver
5. **Compile and test** with simple case
6. **Compare** thermal vs Song models
7. **Validate** against standalone test results

---

## Questions/Decisions Needed

1. ❓ **Sub-cycling:** Should Song solver sub-cycle like thermal model, or trust CONVERGE's dt_sub?
   - Recommendation: Trust dt_sub for now, CONVERGE already sub-cycles

2. ❓ **Termination behavior:** When void=0.99, should we:
   - Option A: Convert to child immediately (like thermal model)
   - Option B: Set pbt=0 and continue with evaporation only
   - Recommendation: Option A for consistency

3. ❓ **Initial conditions:** How is R0 set initially?
   - Current: From critical radius at injection
   - Song: Same approach should work

4. ❓ **Model selection:** Should it be:
   - Global flag (all parcels use same model)
   - Per-parcel flag (could switch mid-simulation)
   - Recommendation: Global flag for simplicity

---

**Ready to implement!** Start with creating RPE_song.h and RPE_song.c.
