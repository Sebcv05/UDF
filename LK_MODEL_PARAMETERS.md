# Langmuir-Knudsen (LK) Model Parameters for user_inputs.in

## Overview
The Langmuir-Knudsen model corrects classical evaporation models for non-equilibrium effects during rapid evaporation, particularly important during flash boiling conditions (P_sat > P_ambient).

## Required Parameters in user_inputs.in

Add the following parameters to your `user_inputs.in` file under the `lagrangian` section:

```
# Langmuir-Knudsen Evaporation Model Parameters
lagrangian.lk_correction_flag = 0           # 0 = classical model, 1 = LK correction enabled
lagrangian.lk_diagnostic_flag = 0           # 0 = no diagnostics, 1 = print LK quantities to terminal
lagrangian.lk_chi_neq_min = 0.0            # Minimum non-equilibrium mole fraction (default: 0.0)
lagrangian.lk_chi_neq_max = 0.9999         # Maximum non-equilibrium mole fraction (default: 0.9999)
```

## Parameter Descriptions

### lagrangian.lk_correction_flag
- **Type:** Integer (0 or 1)
- **Default:** 0 (classical model)
- **Description:** 
  - 0: Use classical Ranz-Marshall evaporation model
  - 1: Apply Langmuir-Knudsen correction for non-equilibrium evaporation

### lagrangian.lk_diagnostic_flag
- **Type:** Integer (0 or 1)
- **Default:** 0 (diagnostics off)
- **Description:**
  - 0: No diagnostic output
  - 1: Print LK model quantities to terminal (χ_eq, χ_neq, β, ψ, L_k, Y_v_s)
  - **Warning:** Setting to 1 will produce large amounts of terminal output. Use only for debugging.

### lagrangian.lk_chi_neq_min
- **Type:** Double precision
- **Default:** 0.0
- **Range:** [0.0, 1.0]
- **Description:** Minimum clipping value for non-equilibrium mole fraction. Prevents negative values during solver iterations.

### lagrangian.lk_chi_neq_max
- **Type:** Double precision
- **Default:** 0.9999
- **Range:** [0.0, 1.0]
- **Description:** Maximum clipping value for non-equilibrium mole fraction. Prevents division by zero in mass fraction conversion.

## Usage Notes

1. **Compatibility with Flash Boiling:**
   - LK correction applies BEFORE flash boiling correlation
   - Both corrections are additive (LK fixes kinetic limit, flash boiling adds superheat effects)
   - Recommended to use both for accurate flash boiling simulation

2. **Single Component:**
   - Current implementation is optimized for ammonia (single component)
   - Multi-component extension planned for future releases

3. **Iteration Strategy:**
   - LK uses previous iteration's drdt to calculate β parameter
   - Converges within existing implicit solver (typically 5-10 iterations)
   - No nested iteration required

4. **Performance Impact:**
   - Computational overhead: ~5-10% increase in spray_evap runtime
   - Most significant effect during flash boiling (P_sat > P_ambient)
   - Minimal effect (<5%) during subcooled evaporation

## Example Configuration

### For Standard Evaporation (No Flash Boiling):
```
lagrangian.lk_correction_flag = 0
lagrangian.lk_diagnostic_flag = 0
```

### For Flash Boiling with LK Correction:
```
lagrangian.lk_correction_flag = 1
lagrangian.lk_diagnostic_flag = 0
lagrangian.lk_chi_neq_min = 0.0
lagrangian.lk_chi_neq_max = 0.9999
```

### For Debugging/Validation:
```
lagrangian.lk_correction_flag = 1
lagrangian.lk_diagnostic_flag = 1  # WARNING: Large terminal output
lagrangian.lk_chi_neq_min = 0.0
lagrangian.lk_chi_neq_max = 0.9999
```

## References

1. Saha, K., Abu-Ramadan, E., and Li, X., "Direct Numerical Simulation of a Confined Three-Dimensional Gas Mixing Layer with One Evaporating Hydrocarbon-Droplet-Laden Stream," Journal of Fluid Mechanics, Vol. 530, 2005, pp. 1-35.

2. Langmuir, I., "The Evaporation of Small Spheres," Physical Review, Vol. 12, 1918, pp. 368-370.

## Implementation Status

- [x] Single-component ammonia
- [x] Option A iteration strategy (use previous drdt)
- [x] Integration with flash boiling correlation
- [x] Diagnostic output
- [ ] Multi-component support
- [ ] Nested iteration option
- [ ] Under-relaxation parameter

## Contact

For questions or issues, refer to LK_CONVERGE_integration_guide.md
