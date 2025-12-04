# Song et al. RPE Test Program

## Overview

This directory contains a standalone C implementation of the Song et al. Rayleigh-Plesset Equation (RPE) model for isothermal bubble growth in superheated ammonia droplets. This serves as:

1. **Validation tool** - Compare C implementation against Python reference
2. **Template for CONVERGE** - The structure mirrors what will be implemented in CONVERGE UDF

## Files

- `test_song_rpe.c` - Standalone C test program (temperature sweep)
- `plot_song_c_results.py` - Python plotting script for C output
- `Makefile_song` - Build system
- `song_temp_sweep_c.csv` - Output data (generated)
- `song_temp_sweep_c*.pdf` - Result plots (generated)

## Building and Running

```bash
# Compile
gcc test_song_rpe.c -o test_song_rpe -O2 -Wall -std=c99 -lm

# Run
./test_song_rpe

# Plot results
python3 plot_song_c_results.py
```

Or use the Makefile:

```bash
make -f Makefile_song
make -f Makefile_song plot
```

## Test Conditions

Temperature sweep at P = 2 bar:
- T = 255, 263, 273, 283, 293, 303, 313, 323 K
- Initial droplet radius: 50 μm
- Initial bubble radius: 1.01 × R_c (critical radius)
- Target void fraction: 0.99

## Key Physics

### Song RPE (Equation 4):
```
(3/2)·ρ_m·Ṙ² + ρ_m·R·R̈ = P_sat - P_∞ + (2σ/R₀ + P_r0)·(R₀/R)³ - 2σ/R - 4μ·Ṙ/R - 4κ·Ṙ/R²
```

### Key Differences from Current RPE_euler.c:
1. **Isothermal** - Temperature held constant (no dT/dt)
2. **No thermal limiting** - Bubble growth driven purely by pressure
3. **Mixture density** - ρ_m = ε·ρ_v + (1-ε)·ρ_l varies with void fraction
4. **Initial pressure term** - Accounts for compressed gas in initial bubble
5. **Only 2 ODEs** - R and Rdot (vs 4 in thermal model)

### Implementation Details:

**State variables:**
- R - Bubble radius
- Rdot - Bubble wall velocity  
- R0 - Initial bubble radius (constant)
- Ro - Initial droplet radius (constant)

**Key parameters:**
- P_r0 = 1.0e6 Pa (residual gas pressure)
- κ = 0.0 (surface viscosity, negligible)
- R_spec = 488.2 J/(kg·K) for NH3

**Integration:**
- Explicit Euler with adaptive timestep
- CFL-like condition: dt ~ R/|Rdot|
- Termination at void fraction = 0.99

## Output

The program generates `song_temp_sweep_c.csv` with columns:
- condition - Temperature label (T=255K, etc.)
- time_s, time_us - Simulation time
- R_m, R_um - Bubble radius
- Rdot_ms - Wall velocity
- epsilon - Void fraction
- rho_m - Mixture density
- P_sat_Pa, P_sat_kPa - Saturation pressure
- film_um - Liquid film thickness

## Validation

Compare with reference Python implementation (`song_temp_sweep.py`):
- Similar timescales (μs range)
- Correct trends with temperature
- Reasonable velocities (tens of m/s)

Expected behavior:
- Higher temperature → Faster growth
- Higher superheat → Higher maximum velocity
- Smooth monotonic radius increase

## Next Steps

Once validated, this code structure will be adapted for CONVERGE UDF:
1. Replace standalone solver with UDF callback
2. Read state from `parcel_cloud` structure
3. Write updated R and Rdot back to parcel
4. Integrate with `spray_drop_distort_NH3.c`

## Notes

- Uses same Antoine coefficients as Python version
- Simple explicit Euler (CONVERGE may use its own integrator)
- Adaptive timestep helps stability
- Logs show: times, radii, velocities, void fractions

## Author

Created: 2025-12-04
Based on: Song et al. RPE formulation + euler_explicit.py reference
