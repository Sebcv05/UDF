# Adding USE_SONG_RPE Flag - Step-by-Step Guide

## Overview
This guide shows how to add a runtime flag to switch between thermal RPE and Song RPE models without recompiling.

---

## Pattern Used in Existing Code

The UDF uses a consistent pattern for user inputs:

1. **Declaration** in `include/globals.h` (with `extern`)
2. **Definition** in a `.c` file (usually `Breakup.c`) with default value
3. **Reading** in `src/read_input.c` from `user_inputs.in` file
4. **Usage** in other `.c` files (e.g., `spray_drop_distort_NH3.c`)

### Example: lk_correction_flag

**In `include/globals.h`:**
```c
extern CONVERGE_index_t lk_correction_flag;
```

**In `src/Breakup.c` (definition with default):**
```c
CONVERGE_index_t lk_correction_flag = 0;
```

**In `src/read_input.c` (reading from file):**
```c
// In UserInputs struct:
struct UserInputs {
    int lk_correction_flag;
    // ... other fields
};

// In read function, default value:
user_inputs->lk_correction_flag = 0;  // Off by default

// In read loop:
else if(strcmp(ktoken, "lk_correction_flag") == 0)
{
    user_inputs->lk_correction_flag = atoi(vtoken);
}

// After reading, assign to global:
lk_correction_flag = (CONVERGE_index_t)user_inputs->lk_correction_flag;

// Logging:
CONVERGE_logger_verbose("user_inputs->lk_correction_flag: %d", user_inputs->lk_correction_flag);

// Echo file:
CONVERGE_file_write(echo, "%-10d lk_correction_flag\n", user_inputs->lk_correction_flag);
```

**In `user_inputs.in` (user specifies value):**
```
0          lk_correction_flag     # 0=off, 1=on
```

---

## Step-by-Step: Adding USE_SONG_RPE Flag

### Step 1: Add to `include/globals.h`

**Location:** After the LK parameters (around line 24)

```c
// Langmuir-Knudsen evaporation model parameters
extern CONVERGE_index_t lk_correction_flag;
extern CONVERGE_index_t lk_diagnostic_flag;
extern CONVERGE_precision_t lk_chi_neq_min;
extern CONVERGE_precision_t lk_chi_neq_max;

// Song RPE model selection
extern CONVERGE_index_t use_song_rpe;    // ADD THIS LINE
```

**Full modified section:**
```c
//Global variables 

#ifndef GLOBALS_H
#define GLOBALS_H

#define MAX_NUM_CHILDREN 150  // Maximum allowed children per breakup event

// Global velocity variables
extern CONVERGE_precision_t user_child_velocity_x;
extern CONVERGE_precision_t user_child_velocity_y;
extern CONVERGE_precision_t user_child_velocity_z ;

// Breakup tuning parameter (default set in Breakup.c, overwritten by user input)
extern CONVERGE_precision_t breakup_velocity_scale;
extern CONVERGE_precision_t breakup_radius_scale;
extern CONVERGE_precision_t kb_threshold;
extern CONVERGE_index_t num_child_parcels;

// Langmuir-Knudsen evaporation model parameters
extern CONVERGE_index_t lk_correction_flag;
extern CONVERGE_index_t lk_diagnostic_flag;
extern CONVERGE_precision_t lk_chi_neq_min;
extern CONVERGE_precision_t lk_chi_neq_max;

// Song RPE model selection
extern CONVERGE_index_t use_song_rpe;

#endif
```

### Step 2: Define in `src/Breakup.c`

**Location:** After LK parameters (around line 33)

```c
// Langmuir-Knudsen evaporation model parameters (set by read_input.c)
CONVERGE_index_t lk_correction_flag = 0;
CONVERGE_index_t lk_diagnostic_flag = 0;
CONVERGE_precision_t lk_chi_neq_min = 0.0;
CONVERGE_precision_t lk_chi_neq_max = 0.9999;

// Song RPE model selection (default: thermal model)
CONVERGE_index_t use_song_rpe = 0;    // ADD THIS LINE
```

### Step 3: Add to `src/read_input.c`

**A. Add to UserInputs struct (around line 23):**
```c
struct UserInputs
{
   double breakup_velocity_scale;
   double breakup_radius_scale;
   double kb_threshold;
   double n_RR;  // Rosin-Rammler shape parameter
   int num_children;  // Number of children per breakup event
   
   // Langmuir-Knudsen evaporation model parameters
   int lk_correction_flag;
   int lk_diagnostic_flag;
   double lk_chi_neq_min;
   double lk_chi_neq_max;
   
   // Song RPE model selection
   int use_song_rpe;    // ADD THIS LINE
};
```

**B. Set default value in initialization (around line 74):**
```c
      // Default LK parameters
      user_inputs->lk_correction_flag = 0;     // Off by default
      user_inputs->lk_diagnostic_flag = 0;     // Off by default
      user_inputs->lk_chi_neq_min = 0.0;
      user_inputs->lk_chi_neq_max = 0.9999;
      
      // Default RPE model selection
      user_inputs->use_song_rpe = 0;           // ADD THIS LINE (0=thermal, 1=Song)
```

**C. Add reading logic in the parse loop (around line 154):**
```c
      else if(strcmp(ktoken, "lk_chi_neq_max") == 0)
      {
         user_inputs->lk_chi_neq_max = atof(vtoken);
      }
      // Song RPE model selection
      else if(strcmp(ktoken, "use_song_rpe") == 0)    // ADD THIS BLOCK
      {
         user_inputs->use_song_rpe = atoi(vtoken);
      }
```

**D. Assign to global variable (around line 199):**
```c
   // Set LK global variables
   lk_correction_flag = (CONVERGE_index_t)user_inputs->lk_correction_flag;
   lk_diagnostic_flag = (CONVERGE_index_t)user_inputs->lk_diagnostic_flag;
   lk_chi_neq_min = (CONVERGE_precision_t)user_inputs->lk_chi_neq_min;
   lk_chi_neq_max = (CONVERGE_precision_t)user_inputs->lk_chi_neq_max;
   
   // Set Song RPE model selection
   use_song_rpe = (CONVERGE_index_t)user_inputs->use_song_rpe;    // ADD THIS LINE
```

**E. Add logging (around line 212):**
```c
   CONVERGE_logger_verbose("user_inputs->lk_correction_flag: %d", user_inputs->lk_correction_flag);
   CONVERGE_logger_verbose("user_inputs->lk_diagnostic_flag: %d", user_inputs->lk_diagnostic_flag);
   CONVERGE_logger_verbose("user_inputs->lk_chi_neq_min: %f", user_inputs->lk_chi_neq_min);
   CONVERGE_logger_verbose("user_inputs->lk_chi_neq_max: %f", user_inputs->lk_chi_neq_max);
   CONVERGE_logger_verbose("user_inputs->use_song_rpe: %d", user_inputs->use_song_rpe);  // ADD THIS LINE
```

**F. Add to echo file (around line 225):**
```c
      CONVERGE_file_write(echo, "%-10d lk_correction_flag\n", user_inputs->lk_correction_flag);
      CONVERGE_file_write(echo, "%-10d lk_diagnostic_flag\n", user_inputs->lk_diagnostic_flag);
      CONVERGE_file_write(echo, "%-10.4f lk_chi_neq_min\n", user_inputs->lk_chi_neq_min);
      CONVERGE_file_write(echo, "%-10.4f lk_chi_neq_max\n", user_inputs->lk_chi_neq_max);
      CONVERGE_file_write(echo, "%-10d use_song_rpe\n", user_inputs->use_song_rpe);  // ADD THIS LINE
```

### Step 4: Add to `user_inputs.in` (in case directory)

**Example location:** `/home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev/user_inputs.in`

```
# Breakup parameters
0.92       breakup_velocity_scale   # aa parameter
10.0       breakup_radius_scale     # B parameter  
1.3        kb_threshold             # kb critical value
3.2        n_RR                     # Rosin-Rammler shape parameter
12         num_children             # Number of children per breakup

# Langmuir-Knudsen evaporation model
0          lk_correction_flag       # 0=off, 1=on
0          lk_diagnostic_flag       # 0=off, 1=on
0.0        lk_chi_neq_min          
0.9999     lk_chi_neq_max

# RPE model selection
0          use_song_rpe             # 0=thermal RPE (default), 1=Song isothermal RPE
```

### Step 5: Use in `src/spray_drop_distort_NH3.c`

**Location:** Around line 505 where RPE_euler_solver is called

**Add include at top:**
```c
#include <RPE_song.h>  // Add after #include <RPE_euler.h>
```

**Replace RPE call:**
```c
// OLD CODE:
RPE_euler_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                 hvap_table, cp_table, num_parcel_species);

// NEW CODE:
if (use_song_rpe) {
    // Song isothermal model
    RPE_song_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                    hvap_table, cp_table, num_parcel_species);
} else {
    // Current thermal model (default)
    RPE_euler_solver(&old_parcel_cloud, p_idx, P_amb, dt_sub,
                     hvap_table, cp_table, num_parcel_species);
}
```

**Optional: Add one-time logging at start of simulation:**
```c
// At the top of spray_distort_cell_NH3 function, add:
static int model_logged = 0;
if (!model_logged) {
    printf("[RPE_MODEL] Using %s RPE model\n", 
           use_song_rpe ? "Song isothermal" : "thermal");
    model_logged = 1;
}
```

---

## Complete File Modifications Summary

### Files to Modify:

1. **`include/globals.h`**
   - Add: `extern CONVERGE_index_t use_song_rpe;`

2. **`src/Breakup.c`**
   - Add: `CONVERGE_index_t use_song_rpe = 0;`

3. **`src/read_input.c`**
   - Add to struct: `int use_song_rpe;`
   - Set default: `user_inputs->use_song_rpe = 0;`
   - Read: `else if(strcmp(ktoken, "use_song_rpe") == 0) { user_inputs->use_song_rpe = atoi(vtoken); }`
   - Assign: `use_song_rpe = (CONVERGE_index_t)user_inputs->use_song_rpe;`
   - Log: `CONVERGE_logger_verbose("user_inputs->use_song_rpe: %d", ...);`
   - Echo: `CONVERGE_file_write(echo, "%-10d use_song_rpe\n", ...);`

4. **`src/spray_drop_distort_NH3.c`**
   - Add include: `#include <RPE_song.h>`
   - Replace RPE call with if/else switch

5. **`user_inputs.in` (case directory)**
   - Add line: `0          use_song_rpe`

---

## Testing the Flag

### Test 1: Verify Reading
```bash
cd /home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev
# Set use_song_rpe = 1 in user_inputs.in
./upc2.sh  # Compile
./run.sh   # Run (will fail if RPE_song.c not implemented yet)
# Check user_inputs.echo - should show: "1          use_song_rpe"
# Check outputs_original/converge.log - should show verbose log
```

### Test 2: Switch at Runtime
```bash
# Edit user_inputs.in: use_song_rpe = 0
./run.sh  # Uses thermal model

# Edit user_inputs.in: use_song_rpe = 1  
./run.sh  # Uses Song model (no recompilation needed!)
```

---

## Advantages of This Approach

✅ **No recompilation** - Just edit `user_inputs.in`
✅ **Clear documentation** - Value echoed to `user_inputs.echo`
✅ **Easy comparison** - Run same case with different models
✅ **Consistent pattern** - Follows existing code style
✅ **Logged output** - Shows which model is active
✅ **Default safe** - Defaults to 0 (thermal model)

---

## Error Handling

If user sets `use_song_rpe = 1` but `RPE_song.c` is not compiled:
- Linker error: "undefined reference to RPE_song_solver"
- Solution: Implement RPE_song.c first, or keep flag at 0

If `user_inputs.in` missing the flag:
- Default value (0) is used automatically
- No error, thermal model runs

---

## Next Steps

After implementing this flag system:
1. Create `include/RPE_song.h`
2. Create `src/RPE_song.c`
3. Modify files as listed above
4. Test compilation
5. Test switching between models

---

**This approach allows runtime model selection without any recompilation!**
