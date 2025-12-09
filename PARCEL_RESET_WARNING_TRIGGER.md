# Where PARCEL_RESET_WARNING Occurs in Code - Execution Flow

## Location in Code Hierarchy

```
CONVERGE_udf_spray_drop_distort()  [spray_drop_distort_NH3.c]
  │
  ├─→ Loop over all parcels (line 187)
  │     │
  │     ├─→ Skip tiny parcels < 1 μm (line 190-198)
  │     │
  │     ├─→ Pre-checks and initialization (lines 200-410)
  │     │
  │     ├─→ Calculate P_sat at current temperature (line 410-411)
  │     │     Saturation_PressureNH3(T_drop, &P_sat_new)
  │     │
  │     ├─→ **SUPERHEAT CHECK** (line 415-428) ← YOUR WARNING ORIGINATES HERE
  │     │     │
  │     │     if (P_sat_new < P_amb) {
  │     │         // Parcel is SUBCOOLED (not superheated)
  │     │         printf("[THERMAL_ABORT] not superheated...");
  │     │         reset_parcel_to_child(..., "Not superheated"); ← CALLS THIS
  │     │         continue;
  │     │     }
  │     │
  │     └─→ Continue with thermal breakup if superheated...
  │
  └─→ reset_parcel_to_child()  [parcel_reset.c line 20]
        │
        ├─→ Get current state (line 25-27)
        │     R_old = parcel_cloud->radius[p_idx];
        │     R_new = parcel_cloud->r_drop_0[p_idx];  // "Injection" radius
        │     Nd_old = parcel_cloud->num_drop[p_idx];
        │
        ├─→ **SAFETY CHECK** (line 36-45) ← WARNING PRINTED HERE
        │     │
        │     if (R_new > R_old * 2.0) {  // 82.5 μm > 2 × 36 μm = FAIL
        │         printf("[PARCEL_RESET_WARNING] Invalid r_drop_0=%.3e...");
        │         R_new = R_old;  // Use current radius instead
        │         Nd_new = Nd_old;
        │     }
        │
        ├─→ Perform reset (line 48-56)
        │     parcel_cloud->radius[p_idx] = R_new;
        │     parcel_cloud->num_drop[p_idx] = Nd_new;
        │     parcel_cloud->r_bubble[p_idx] = 0.0;
        │     parcel_cloud->is_child[p_idx] = 1;
        │     parcel_cloud->thermal_breakup_flag[p_idx] = 999;
        │
        └─→ Log reset details (line 58-69)
              printf("[PARCEL_RESET] ...");
```

---

## Detailed Trigger Sequence

### Step 1: Parcel Loop Starts (spray_drop_distort_NH3.c, line 187)
```c
for (int p_idx = 0; p_idx < num_parcels_in_cloud; p_idx++) {
```

### Step 2: Calculate Saturation Pressure (line 410-411)
```c
CONVERGE_precision_t P_sat_new;
Saturation_PressureNH3(old_parcel_cloud.temp[p_idx], &P_sat_new);
```

**For your parcels:**
- `T_drop = 260-270 K` (parcel has cooled)
- `P_sat_new = Antoine(260-270K) ≈ 0.5-0.8 bar` (using NH3 Antoine equation)
- `P_amb = 1.0 bar` (ambient pressure in chamber)

### Step 3: Superheat Check (line 415-428) **← TRIGGER**
```c
// Superheat check: Parcel must be superheated to enter thermal breakup
// Condition: P_sat(T_drop) >= P_ambient
if (P_sat_new < P_amb)  // 0.7 bar < 1.0 bar = TRUE (subcooled!)
{
    // Not superheated - disable thermal breakup for this parcel
    static int not_superheated_count = 0;
    if (not_superheated_count < 10) {
        printf("[THERMAL_ABORT] p_idx=%li, not superheated: P_sat(T_drop)=%.3e < P_amb=%.3e Pa\n",
               p_idx, P_sat_new, P_amb);
        printf("                T_drop=%.2f K, lifetime=%.3e s, disabling thermal breakup\n",
               old_parcel_cloud.temp[p_idx], old_parcel_cloud.lifetime[p_idx]);
        not_superheated_count++;
    }
    reset_parcel_to_child(&old_parcel_cloud, p_idx, "Not superheated"); ← CALLS RESET
    continue;  // Skip rest of thermal breakup processing
}
```

**Why this triggers:**
- Parcel temperature has dropped due to evaporative cooling
- At lower temperature, saturation pressure is lower
- When `P_sat < P_amb`, the parcel is **subcooled** (not superheated)
- Subcooled parcels cannot undergo thermal breakup (no driving force for bubble growth)

### Step 4: Reset Function Called (parcel_reset.c, line 20)
```c
void reset_parcel_to_child(struct ParcelCloud* parcel_cloud, 
                            CONVERGE_index_t p_idx, 
                            const char* reason)
{
    // Get current state
    CONVERGE_precision_t R_old = parcel_cloud->radius[p_idx];      // e.g., 36 μm
    CONVERGE_precision_t R_new = parcel_cloud->r_drop_0[p_idx];    // e.g., 82.5 μm (parent!)
    CONVERGE_precision_t Nd_old = parcel_cloud->num_drop[p_idx];
```

### Step 5: Safety Check Detects Mismatch (line 36-45) **← WARNING PRINTED**
```c
// Safety check: if injection radius is invalid, use current radius
if (R_new <= 0.0 || R_new > R_old * 2.0) {
    // 82.5e-6 > 2 × 36e-6 = TRUE → INVALID
    static int warning_count = 0;
    if (warning_count < 5) {
        printf("[PARCEL_RESET_WARNING] Invalid r_drop_0=%.3e m, using current radius %.3e m (p_idx=%li, reason: %s)\n",
               R_new, R_old, p_idx, reason);
        //       82.5 μm       36 μm      1       "Not superheated"
        warning_count++;
    }
    R_new = R_old;   // Override: keep current radius (36 μm)
    Nd_new = Nd_old; // Keep current num_drop (mass conserved)
}
```

### Step 6: Reset Applied (line 48-56)
```c
// Perform reset
parcel_cloud->radius[p_idx] = R_new;          // = 36 μm (NOT 82.5 μm)
parcel_cloud->num_drop[p_idx] = Nd_new;       // unchanged
parcel_cloud->r_bubble[p_idx] = 0.0;          // zero bubble
parcel_cloud->r_bubble_0[p_idx] = 0.0;
parcel_cloud->v_bubble[p_idx] = 0.0;
parcel_cloud->pbt[p_idx] = 0;                 // disable pre-breakup tag
parcel_cloud->is_child[p_idx] = 1;            // mark as child
parcel_cloud->thermal_breakup_flag[p_idx] = 999;  // disable thermal breakup
```

---

## What Triggers This Specific Scenario?

### Parcel History Leading to Warning:

1. **t = 0 μs: Injection**
   - Parent injected at `R = 82.5 μm`, `T = 290 K` (hot)
   - `r_drop_0 = 82.5 μm` (stored)
   - `is_child = 0` (parent)

2. **t = 20 μs: KH-RT Breakup (Not Thermal)**
   - Parent breaks via standard KH-RT mechanism
   - Creates children: `R = 30-40 μm`
   - **Child inherits:** `r_drop_0 = 82.5 μm` (parent's radius!)
   - Child: `is_child = 1`, `radius = 36 μm`

3. **t = 50-100 μs: Evaporation**
   - Child evaporates in ambient gas
   - Loses mass → `radius = 36 μm → 33 μm → 30 μm`
   - Loses heat → `T = 290 K → 270 K → 265 K`

4. **t = 157 μs: Subcooled Check Triggers (YOUR LOG)**
   - Child has cooled: `T = 265 K`
   - Calculate: `P_sat(265K) ≈ 0.75 bar`
   - Compare: `0.75 bar < P_amb (1.0 bar)` → **SUBCOOLED**
   - Call `reset_parcel_to_child(..., "Not superheated")`

5. **Reset Function Detects Mismatch:**
   - Tries to use `r_drop_0 = 82.5 μm`
   - But current `radius = 26-36 μm`
   - **Safety check:** `82.5 > 2 × 36` → FAIL
   - **Warning printed:** `[PARCEL_RESET_WARNING]`
   - **Action:** Use current radius instead (36 μm)

---

## Why `r_drop_0` is Wrong for Children

From `parcel_prop.c` (child creation during KH-RT breakup):

```c
// Line 286: When creating child from parent
parcel_cloud.r_drop_0[passed_child_parcel_idx] = parcel_cloud.radius[passed_parent_parcel_idx];
//                                                 ↑
//                                                 This is PARENT radius at breakup time
//                                                 NOT the child's "injection" radius
```

**Result:** Children inherit parent's `r_drop_0`, which is meaningless for them.

---

## Code Locations Summary

| Event | File | Line | Code |
|-------|------|------|------|
| **Parcel loop** | `spray_drop_distort_NH3.c` | 187 | `for (int p_idx = 0; ...)` |
| **Superheat check** | `spray_drop_distort_NH3.c` | 415-428 | `if (P_sat_new < P_amb)` |
| **Reset called** | `spray_drop_distort_NH3.c` | 426 | `reset_parcel_to_child(..., "Not superheated")` |
| **Safety check** | `parcel_reset.c` | 36-45 | `if (R_new > R_old * 2.0)` |
| **Warning printed** | `parcel_reset.c` | 39 | `printf("[PARCEL_RESET_WARNING] ...")` |
| **Child inherits r_drop_0** | `parcel_prop.c` | 286 | `r_drop_0[child] = radius[parent]` |

---

## Physical Interpretation

**Your parcels are:**
1. ✅ Small children from KH-RT breakup (30-40 μm)
2. ✅ Evaporating and cooling normally
3. ✅ Correctly being disabled from thermal breakup (they're subcooled)
4. ✅ Safety check preventing artificial expansion to parent size

**This is expected behavior** for:
- Spray simulations with ambient T < boiling point
- Small droplets that evaporate and cool rapidly
- Post-breakup children that are no longer superheated

---

## How to See This in Action

Add this diagnostic right before the safety check in `parcel_reset.c` (line 35):

```c
// DIAGNOSTIC: Show why safety check triggers
printf("[RESET_DIAGNOSTIC] p_idx=%li, reason=%s\n", p_idx, reason);
printf("                   R_old=%.3e m, R_new(r_drop_0)=%.3e m, ratio=%.2f\n",
       R_old, R_new, R_new/R_old);
printf("                   is_child=%d (0=parent, 1=child)\n", 
       parcel_cloud->is_child[p_idx]);
```

This will show you that `is_child=1` (child) and the ratio `82.5/36 ≈ 2.3` exceeds the 2.0 limit.

---

## Summary

**Trigger:** Subcooled child parcels (from KH-RT breakup) being reset → safety check detects `r_drop_0` mismatch → warning printed

**Harmless:** Safety check prevents incorrect radius expansion, parcel behaves correctly

**Expected:** Normal behavior for cooling child parcels in spray simulation
