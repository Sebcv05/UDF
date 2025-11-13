# Hijacking Unused Parcel Variables for Diagnostics

## Strategy: Use Unused Model Variables

Instead of CSV output or cell-centered data, **write diagnostic data into unused parcel variables** that are already in the H5 output.

---

## Available Unused Variables

From your list, these are likely unused:

### Model Variables (Always Output):
- **`radiation_energy`** ← Perfect for is_child (0 or 1)
- **`distant`** ← Can use for r_bubble
- **`t_turb`** ← Can use for pbt flag
- **`time_turb_accum`** ← Can use for additional diagnostic

---

## Implementation

### Option 1: Set in Breakup.c (When Children Created)

When creating child parcels, set the diagnostic variables:

```c
// In Breakup.c, when creating children (around line 723):

// Set is_child flag in radiation_energy
new_parcel_cloud.radiation_energy[new_idx] = 1.0;  // 1 = child

// Set r_bubble in distant variable
new_parcel_cloud.distant[new_idx] = new_parcel_cloud.r_bubble[new_idx];

// Set pbt flag in t_turb variable  
new_parcel_cloud.t_turb[new_idx] = (CONVERGE_precision_t)new_parcel_cloud.pbt[new_idx];
```

For parent parcels (at injection):
```c
// In parcel_prop.c, parcel_inject function:

// Set is_child flag to 0 (parent)
parcel_cloud.radiation_energy[passed_parcel_idx] = 0.0;

// Set initial r_bubble
parcel_cloud.distant[passed_parcel_idx] = parcel_cloud.r_bubble[passed_parcel_idx];

// Set pbt flag
parcel_cloud.t_turb[passed_parcel_idx] = (CONVERGE_precision_t)parcel_cloud.pbt[passed_parcel_idx];
```

### Option 2: Update During RPE Solver (Continuous)

Update these diagnostic variables during bubble growth:

```c
// In spray_drop_distort_NH3.c, in the main parcel loop:

// Copy diagnostic data to output variables
parcel_cloud.radiation_energy[p_idx] = (CONVERGE_precision_t)parcel_cloud.is_child[p_idx];
parcel_cloud.distant[p_idx] = parcel_cloud.r_bubble[p_idx];
parcel_cloud.t_turb[p_idx] = (CONVERGE_precision_t)parcel_cloud.pbt[p_idx];
```

---

## Recommended Mapping

| Diagnostic Variable | Hijacked Output Variable | Units in Tecplot | Meaning |
|---------------------|-------------------------|------------------|---------|
| **is_child** | `radiation_energy` | J | 0=parent, 1=child |
| **r_bubble** | `distant` | m | Bubble radius (0-8e-5 m) |
| **pbt** | `t_turb` | s | 0=no breakup, 1=in breakup |

---

## Update post.in

Change from trying to use cell-centered user variables to using the model variables:

```yaml
parcels:
   liquid_parcels:
      physical:           [from_injector, from_nozzle, num_drop, radius, temp, velocity]
      model:              [tbreak_kh, tbreak_rt, radiation_energy, distant, t_turb]
```

**Remove the cell-centered user variables:**
```yaml
cells:
   geometry:              [grid_level, rank, idreg, xcen]
   general:               [density, mass, vol_frac, pressure, temperature, velocity, volume]
   turbulence:            [eps, tke, turb_viscosity, yplus]
   species_massfrac:      [NH3, O2, N2]
   species_molefrac:      [NH3, O2, N2]
   # REMOVE: user: [user_is_child, user_r_bubble, user_pbt]
```

---

## Implementation Code

### File 1: spray_drop_distort_NH3.c

Add at the end of the main parcel loop (after RPE solver, DGRE, etc.):

```c
// Copy diagnostic data to output variables (hijack unused fields)
// This makes is_child, r_bubble, pbt visible in Tecplot as parcel variables

for(CONVERGE_index_t p_idx = 0; p_idx < num_parcels; p_idx++)
{
   // is_child → radiation_energy
   parcel_cloud.radiation_energy[p_idx] = (CONVERGE_precision_t)parcel_cloud.is_child[p_idx];
   
   // r_bubble → distant  
   parcel_cloud.distant[p_idx] = parcel_cloud.r_bubble[p_idx];
   
   // pbt → t_turb
   parcel_cloud.t_turb[p_idx] = (CONVERGE_precision_t)parcel_cloud.pbt[p_idx];
}
```

### File 2: parcel_prop.c (parcel_inject)

Initialize for parent parcels:

```c
// Around line 175, after setting r_bubble:

// Hijack output variables for diagnostics
parcel_cloud.radiation_energy[passed_parcel_idx] = 0.0;  // 0 = parent
parcel_cloud.distant[passed_parcel_idx] = parcel_cloud.r_bubble[passed_parcel_idx];
parcel_cloud.t_turb[passed_parcel_idx] = 0.0;  // Not in breakup yet
```

### File 3: parcel_prop.c (parcel_child) 

Set for child parcels:

```c
// Around line 280, after setting other child properties:

// Hijack output variables for diagnostics  
parcel_cloud.radiation_energy[passed_child_parcel_idx] = 1.0;  // 1 = child
parcel_cloud.distant[passed_child_parcel_idx] = 0.0;  // Child has no bubble
parcel_cloud.t_turb[passed_child_parcel_idx] = 0.0;  // Not in breakup
```

---

## In Tecplot

After recompiling and running, you'll see in **Parcel Data**:

### radiation_energy
- 0 = parent parcel (injected, not yet broken up)
- 1 = child parcel (result of breakup)

**Filter for large parent parcels:**
```
radius > 7e-5 AND radiation_energy = 0
```

**Filter for large child parcels (RR tail):**
```
radius > 7e-5 AND radiation_energy = 1
```

### distant
- Bubble radius in meters
- 0 = no bubble or child parcel
- 1e-9 to 8e-5 = active bubble growing

**Filter for parcels about to break up:**
```
distant > 7e-5
```

### t_turb  
- 0 = not in thermal breakup
- 1 = actively breaking up

---

## Advantages of This Approach

✅ **Per-parcel data** - Every parcel has its own values  
✅ **No cell averaging** - Direct parcel-to-data mapping  
✅ **Native Tecplot** - Works with standard parcel visualization  
✅ **No CSV files** - All in H5 output  
✅ **Minimal code** - Just copy values to existing fields  
✅ **Already compiled** - parcel_cloud struct already has these fields

---

## Quick Implementation Summary

1. **Add to spray_drop_distort_NH3.c:** Copy diagnostic data to output variables at end of parcel loop
2. **Add to parcel_prop.c (inject):** Initialize for parent parcels  
3. **Add to parcel_prop.c (child):** Initialize for child parcels
4. **Update post.in:** Add `radiation_energy, distant, t_turb` to model output
5. **Compile:** Run `./upc2.sh`
6. **Run:** Normal CONVERGE run
7. **Tecplot:** Load H5, filter by `radiation_energy` to see parent vs child

---

**Estimated code changes:** ~15 lines total  
**Compilation time:** Same as normal  
**No new files needed:** Uses existing parcel_cloud struct

---

**Created:** 2025-11-13  
**Status:** Ready to implement
