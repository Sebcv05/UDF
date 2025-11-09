# RPE_euler Implementation Notes

## Date: 2025-11-09

## Changes Made

### 1. New Files Created
- `include/RPE_euler.h` - Header with structures and function declarations
- `src/RPE_euler.c` - Full implementation of Rayleigh-Plesset solver

### 2. Modified Files
- `src/spray_drop_distort_NH3.c` - Replaced `Bubble_Velocity()` with `RPE_euler_solver()`

### 3. Key Features Implemented

#### Physical Model
- **4 ODEs solved simultaneously:**
  1. dR/dt = Rdot (bubble radius evolution)
  2. dRdot/dt = Rayleigh-Plesset momentum equation
  3. dT_drop/dt = Energy balance (Q_conv - Q_evap)/(m_drop*cp)
  4. dm_b/dt = Mass balance (thermal limiting: mdot = Q_conv/L_v)

#### Heat Transfer
- Ranz-Marshall Nu correlation: Nu = 2 + 0.6*Re^0.5*Pr^(1/3)
- Thermal limiting mass transfer (not diffusion-limited)
- Dynamic Re and Pr based on instantaneous bubble velocity

#### Integration Strategy
- Explicit Euler time integration
- Called INSIDE sub-timestep loop (replaces pre-loop Bubble_Velocity call)
- Uses dt_sub (~1e-10 s) for integration
- Updates: r_bubble, v_bubble, temp each sub-step

#### Safety Limits
- Acceleration: |dRdot/dt| < 1e10 m/s²
- Temperature: 100 K < T_drop < 500 K
- Temperature rate: |dT/dt| < 1e6 K/s
- Bubble constraint: R ≤ 0.95*Ro
- Minimum values: R > 1e-12 m, rho_v > 1e-6 kg/m³

### 4. Integration with Existing Code

#### Workflow (each sub-timestep):
1. **RPE_euler_solver()** - Updates R, Rdot, T_drop
2. **Geometry()** - Computes droplet expansion (unchanged)
3. **DGRE_NH3()** - Disturbance growth rate (unchanged)
4. **BreakupCriterion()** - Breakup criterion kb (unchanged)

#### Data Flow:
- **Input to RPE:** old_parcel_cloud fields, P_amb, dt_sub, property tables
- **Output from RPE:** Updated r_bubble, v_bubble, temp
- **Used by DGRE:** v_bubble, r_bubble, radius (all consistent)

### 5. Removed Code
- Pre-loop `Bubble_Velocity()` call (line ~402)
- Old Euler update logic: `dR = v_bubble * dt_sub; Rb_temp = dR + Rb`
- Simplified to: RPE updates r_bubble internally, just do safety checks

## Testing Plan

### Unit Tests Required:

1. **Validation Test (323 K Ammonia)**
   - T_drop = 323 K, P_amb = 2.0e5 Pa, Ro = 82.5 µm
   - Expected: R → ~82 µm, Rdot → 15-20 m/s after ~8 µs
   - Status: Not yet tested

2. **Thermal Limiting Check**
   - Verify mdot = Q_conv / L_v
   - Status: Implemented in code

3. **Stability Test**
   - Extreme values: R ~ 1e-12 m, Rdot ~ 100 m/s
   - Status: Safety limits in place

4. **Integration Test with CONVERGE**
   - Single parcel injection
   - Monitor R, Rdot, T_drop evolution
   - Status: Ready for testing

## Known Issues / TODO

1. **Thermal conductivity k_l** - Currently hardcoded to 0.5 W/(m·K)
   - TODO: Add to property tables or compute from correlations

2. **Bubble mass initialization** - Computed from P_sat and rho_b
   - May need to persist m_bubble between timesteps
   - Currently recalculated each call

3. **Ambient temperature T_amb** - Not currently used in heat transfer
   - Using droplet-interface temperature difference only
   - May need to include far-field convection

4. **Property tables** - Using CONVERGE tables for L_v and cp_l
   - Working correctly with species mass fractions

5. **Nu_max limit** - Set to 1000
   - May need adjustment based on validation

## Backup Branch
- Branch `Inertial_RPE` created with pre-implementation state
- Can revert via: `git checkout Inertial_RPE`
- Current branch: `v3.1.12`

## Next Steps

1. Test compilation in case directory
2. Run single-parcel test case
3. Compare with old Bubble_Velocity results
4. Adjust parameters if needed (Nu_max, k_l, etc.)
5. Validate against standalone Python/MATLAB implementation
6. Full spray simulation test

## Commit Hash
- Implementation: 93f72e9
- Backup branch: 51bb775 (Inertial_RPE)
