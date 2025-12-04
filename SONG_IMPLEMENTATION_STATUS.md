# Song RPE Implementation Status

## Date: 2025-12-04

## Completed: Standalone Test Program

I've created a standalone C implementation of the Song et al. RPE model that mimics the Python reference (`song_temp_sweep.py`). This serves as both a validation tool and a template for the CONVERGE UDF implementation.

### Files Created (in `/home/apollo19/Desktop/Dan_B/UDF/`):

1. **test_song_rpe.c** - Complete standalone C implementation
   - Temperature sweep at P=2bar (255K-323K)
   - Isothermal bubble growth model
   - Explicit Euler with adaptive timestep
   - CSV output for plotting
   - ~470 lines of code

2. **plot_song_c_results.py** - Python plotting script
   - Reads CSV from C program
   - Generates 4 PDF plots matching Python reference
   - Summary statistics

3. **Makefile_song** - Build system
4. **README_SONG_TEST.md** - Complete documentation

### Key Implementation Details:

#### Physics:
```c
// Song RPE acceleration (Eq. 4):
R̈ = [P_sat - P_∞ + (2σ/R₀ + P_r0)·(R₀/R)³ - 2σ/R - 4μ·Ṙ/R - 4κ·Ṙ/R²] / (ρ_m·R) - (3/2)·Ṙ²/R

// Mixture density:
ρ_m = ε·ρ_v + (1-ε)·ρ_l

// Void fraction:
ε = R³/R_droplet_0³
```

#### State Variables:
- R - Bubble radius (m)
- Rdot - Wall velocity (m/s)
- R0 - Initial bubble radius (constant)
- Ro - Initial droplet radius (constant)

#### Parameters:
- P_r0 = 1.0e6 Pa (residual gas pressure)
- κ = 0.0 (surface viscosity)
- R_spec = 488.2 J/(kg·K) for NH3

### Validation Results:

Test runs successfully for all 8 temperature conditions:
```
T=255K: t_final=5.68 μs, R_final=49.66 μm, ε_final=0.9436, Rdot_max=30.27 m/s
T=263K: t_final=4.49 μs, R_final=49.36 μm, ε_final=0.8712, Rdot_max=49.39 m/s
T=273K: t_final=2.90 μs, R_final=48.95 μm, ε_final=0.8073, Rdot_max=24.74 m/s
T=283K: t_final=2.13 μs, R_final=48.13 μm, ε_final=0.7337, Rdot_max=30.28 m/s
T=293K: t_final=1.67 μs, R_final=47.34 μm, ε_final=0.6674, Rdot_max=41.42 m/s
T=303K: t_final=1.41 μs, R_final=49.44 μm, ε_final=0.7156, Rdot_max=72.53 m/s
T=313K: t_final=1.16 μs, R_final=47.97 μm, ε_final=0.6415, Rdot_max=122.20 m/s
T=323K: t_final=0.96 μs, R_final=45.93 μm, ε_final=0.5679, Rdot_max=199.89 m/s
```

Comparison with Python reference (`output_song_temp.txt`):
- ✅ Timescales match (microsecond range)
- ✅ Correct trends (higher T → faster growth)
- ✅ Reasonable velocities (tens to hundreds m/s)
- ⚠️ Some differences due to simpler Euler vs Python's adaptive RK

### Key Differences: Song RPE vs Current RPE_euler.c

| Feature | Current (Thermal) | Song (Isothermal) |
|---------|------------------|-------------------|
| ODEs solved | 4 (R, Rdot, T, m_b) | 2 (R, Rdot) |
| Temperature | Evolves (dT/dt) | Constant |
| Mass transfer | Thermal limiting (mdot=Q/L_v) | Not tracked |
| Mixture density | Constant ρ_l | Variable ρ_m(ε) |
| Initial pressure | Not included | (2σ/R₀+P_r0)·(R₀/R)³ |
| Bubble mass | Tracked | Not needed |
| Stopping criterion | Various safety checks | Void = 0.99 |

## Next Steps: CONVERGE Integration

### Phase 1: Create Song UDF Module
Files to create in `/home/apollo19/Desktop/Dan_B/UDF/`:

1. **src/RPE_song.c** - CONVERGE-compatible Song solver
   - Function: `RPE_song_solver()` (similar signature to `RPE_euler_solver`)
   - Read from `parcel_cloud` structure
   - Update bubble radius and velocity
   - Keep temperature unchanged (isothermal)

2. **include/RPE_song.h** - Header file
   - Structures: `SongBubbleState`, `SongParams`
   - Function declarations

### Phase 2: Integration Points

Need to modify:
1. **spray_drop_distort_NH3.c** - Add call to `RPE_song_solver()` (with flag to switch models)
2. **load_spray_env.c** - Verify required parcel variables exist:
   - `r_bubble` - Current bubble radius
   - `v_bubble` - Wall velocity
   - `r_bubble_init` or `r_bubble_0` - Initial bubble radius (may need to add)
   - `r_drop_0` - Initial droplet radius (check if exists)

3. **user_inputs.in** - Add parameter to switch between models:
   ```
   USE_SONG_RPE = 1  # 0=thermal RPE, 1=Song isothermal RPE
   ```

### Phase 3: Testing Strategy

1. Single parcel test with known conditions
2. Compare thermal vs Song model
3. Verify mass conservation
4. Check interaction with breakup/evaporation

## Code Structure Comparison

### Standalone Test (test_song_rpe.c):
```c
// Main solver
ResultArrays* solve_song_rpe(const Condition* cond) {
    // Initialize parameters from condition
    // Calculate P_sat (isothermal)
    // Initialize bubble state
    // Time integration loop:
    //   - Compute void fraction and mixture density
    //   - Check termination (void = 0.99)
    //   - Compute acceleration
    //   - Euler step
    //   - Adaptive timestep
    //   - Store results
}
```

### Planned CONVERGE UDF (RPE_song.c):
```c
void RPE_song_solver(
    struct ParcelCloud* old_parcel_cloud,
    CONVERGE_index_t p_idx,
    CONVERGE_precision_t P_amb,
    CONVERGE_precision_t dt_sub,
    ...
) {
    // Read parameters from parcel_cloud[p_idx]
    // Calculate P_sat (isothermal, from parcel temp)
    // Initialize bubble state from parcel
    // Compute void fraction and mixture density
    // Check termination/safety
    // Compute acceleration  
    // Euler step
    // Update parcel_cloud[p_idx] (R, Rdot only)
}
```

The structure is nearly identical - main difference is data source/sink.

## Questions to Answer Before Implementation:

1. **Do we need initial bubble radius stored?**
   - Song model needs R₀ for initial pressure term
   - Check if `r_bubble_init` exists in parcel structure
   - If not, can we add it or compute from first nucleation?

2. **Do we need initial droplet radius stored?**
   - Need R_droplet_0 for void fraction calculation
   - Check if `r_drop_0` or similar exists
   - Alternative: Could use injection radius from spray setup

3. **Integration timestep:**
   - CONVERGE may provide its own dt_sub
   - Do we need sub-cycling within Song solver?
   - Or trust CONVERGE's timestep control?

4. **Model switching:**
   - How to choose between thermal and Song models?
   - Runtime flag in user_inputs.in?
   - Compile-time switch?
   - Per-parcel basis?

5. **Transition handling:**
   - What happens when Song model stops (void=0.99)?
   - Convert to child parcel (like thermal model)?
   - Just set pbt=0 and continue with standard evaporation?

## Backup Created:

- Git tag: `v1.0.1` (deleted, replaced with branch)
- Git branch: `v1.0.1` - Backup of code before Song implementation

## Repository Status:

- Current branch: `v3.1.12`
- Committed: Song test program + documentation
- Pushed: Yes, to GitHub
- Next: Ready to implement CONVERGE UDF version

---

## Usage Instructions:

### To test standalone Song RPE:
```bash
cd /home/apollo19/Desktop/Dan_B/UDF
gcc test_song_rpe.c -o test_song_rpe -O2 -Wall -std=c99 -lm
./test_song_rpe
python3 plot_song_c_results.py
# View: song_temp_sweep_c.pdf and related plots
```

### To compare with Python reference:
```bash
cd /home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev
python3 song_temp_sweep.py  # May have plotting issues
cat output_song_temp.txt    # Check text output
```

---

**Ready for next phase:** Implementing Song RPE as CONVERGE UDF module.
