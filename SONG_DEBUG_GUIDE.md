# Song RPE Debugging and Testing Outputs

## Overview
The Song RPE implementation includes multiple levels of diagnostic output to validate behavior and compare with the thermal model.

---

## 1. Model Selection Confirmation (Startup)

**File:** `src/read_input.c`  
**Triggered:** When `user_inputs.in` is read at simulation start  
**Frequency:** Once per simulation

### Output with `use_song_rpe=0`:
```
========================================
RPE MODEL: Thermal (use_song_rpe=0)
========================================
```

### Output with `use_song_rpe=1`:
```
========================================
RPE MODEL: Song isothermal (use_song_rpe=1)
========================================
```

**Purpose:** Immediate confirmation of which model is active before any parcels run.

**How to check:**
```bash
grep "RPE MODEL" outputs_original/converge.log
```

---

## 2. Song Model Activation

**File:** `src/RPE_song.c` (lines 119-126)  
**Triggered:** First time Song solver is called (when first parcel enters thermal breakup)  
**Frequency:** Once per simulation

### Output:
```
========================================
[RPE_MODEL_SONG] Song isothermal RPE model active
========================================
```

**Purpose:** Confirms Song solver is actually being called, not just selected.

**How to check:**
```bash
grep "RPE_MODEL_SONG" outputs_original/converge.log
```

---

## 3. Detailed Step Logging

**File:** `src/RPE_song.c` (lines 193-198)  
**Triggered:** First 10 solver calls  
**Limited by:** `song_debug_count` static counter

### Output Format:
```
[SONG_STEP] p_idx=%li, R=%.3e m, Rdot=%.3e m/s, ε=%.4f, ρ_m=%.1f kg/m³, T=%.2f K
```

### Example:
```
[SONG_STEP] p_idx=42, R=1.234e-06 m, Rdot=5.678e-02 m/s, ε=0.0123, ρ_m=650.5 kg/m³, T=293.15 K
```

### What it shows:
- **p_idx**: Parcel index
- **R**: Current bubble radius (meters)
- **Rdot**: Bubble wall velocity (m/s)
- **ε (epsilon)**: Void fraction (0 to 1)
- **ρ_m**: Mixture density (kg/m³) - varies with void fraction
- **T**: Droplet temperature (K) - should remain constant for Song

**Purpose:** Track bubble growth evolution, verify physics calculations.

**How to check:**
```bash
grep "SONG_STEP" outputs_original/converge.log
```

---

## 4. Superheat Check Failure

**File:** `src/RPE_song.c` (lines 161-166)  
**Triggered:** When `P_sat ≤ P_amb` (liquid not superheated)  
**Limited:** First 10 occurrences per simulation

### Output:
```
[SONG_ABORT] Not superheated: P_sat=%.2e Pa, P_amb=%.2e Pa
```

### Example:
```
[SONG_ABORT] Not superheated: P_sat=1.95e+05 Pa, P_amb=2.00e+05 Pa
```

**Purpose:** Diagnose why parcels aren't entering or continuing thermal breakup.

**How to check:**
```bash
grep "SONG_ABORT" outputs_original/converge.log
```

---

## 5. Void Fraction Target Reached

**File:** `src/RPE_song.c` (lines 174-180)  
**Triggered:** When void fraction ≥ 0.99 (bubble growth complete)  
**Limited:** First 10 occurrences per simulation

### Output:
```
[SONG_COMPLETE] Void fraction target reached: ε=%.4f >= %.2f
```

### Example:
```
[SONG_COMPLETE] Void fraction target reached: ε=0.9901 >= 0.99
```

**Purpose:** Confirm successful completion of bubble growth.

**How to check:**
```bash
grep "SONG_COMPLETE" outputs_original/converge.log
```

---

## 6. Bubble Collapse Detection

**File:** `src/RPE_song.c` (lines 210-214)  
**Triggered:** When `Rdot < 0` (bubble shrinking)  
**Limited:** No limit (always prints)

### Output:
```
[SONG_COLLAPSE] Bubble collapsing: Rdot=%.3e m/s at p_idx=%li
```

### Example:
```
[SONG_COLLAPSE] Bubble collapsing: Rdot=-1.234e-02 m/s at p_idx=42
```

**Purpose:** Detect unphysical behavior. Bubbles shouldn't shrink in superheated liquid.  
**Action:** Investigate if this appears - may indicate numerical instability or incorrect parameters.

**How to check:**
```bash
grep "SONG_COLLAPSE" outputs_original/converge.log
```

---

## 7. Edge Capping

**File:** `src/RPE_song.c` (lines 217-221)  
**Triggered:** When bubble reaches droplet edge (`R > 0.999 * Ro`)  
**Limited:** No limit (always prints)

### Output:
```
[SONG_EDGE] Bubble capped at droplet edge: R=%.3e m, Ro=%.3e m
```

### Example:
```
[SONG_EDGE] Bubble capped at droplet edge: R=4.995e-05 m, Ro=5.000e-05 m
```

**Purpose:** Show when bubble growth is artificially limited by droplet size.

**How to check:**
```bash
grep "SONG_EDGE" outputs_original/converge.log
```

---

## 8. Error Messages (Invalid Initial Conditions)

**File:** `src/RPE_song.c` (lines 144-154)  
**Triggered:** Bad parcel initialization  
**Limited:** No limit (always prints)

### Output:
```
[SONG_ERROR] Droplet radius too small: Ro=%.3e m
[SONG_ERROR] Initial bubble radius too small: R0=%.3e m
```

**Purpose:** Catch bad parcel initialization that could cause crashes.

**How to check:**
```bash
grep "SONG_ERROR" outputs_original/converge.log
```

---

## Comparison: Song vs Thermal Debug Output

| Message Type | Song Model | Thermal Model |
|--------------|------------|---------------|
| Model selection | `RPE MODEL: Song...` | `RPE MODEL: Thermal...` |
| Activation | `[RPE_MODEL_SONG]` | `[RPE_MODEL_THERMAL]` |
| Step logging | `[SONG_STEP]` | Various thermal logs |
| Completion | `[SONG_COMPLETE]` | `[RPE_POST_RECOVERY]` |
| Superheat fail | `[SONG_ABORT]` | `[THERMAL_ABORT]` |
| Collapse | `[SONG_COLLAPSE]` | Recovery logic triggers |

---

## Testing Workflow

### 1. Quick Verification
```bash
cd /path/to/case
grep "RPE MODEL" outputs_original/converge.log
```
Expected: One line showing either "Thermal" or "Song isothermal"

### 2. Check Song Activation
```bash
grep "RPE_MODEL_SONG" outputs_original/converge.log
```
Expected: One line when first parcel enters thermal breakup

### 3. Monitor Bubble Growth
```bash
grep "SONG_STEP" outputs_original/converge.log
```
Expected: Up to 10 lines showing bubble evolution

### 4. Check Termination Reasons
```bash
grep -E "SONG_COMPLETE|SONG_ABORT|SONG_COLLAPSE|SONG_EDGE" outputs_original/converge.log
```
Expected: Messages explaining why bubbles stopped growing

### 5. Compare Models
```bash
# Run with use_song_rpe=0
grep -E "THERMAL|RPE" outputs_original/converge.log > thermal_run.txt

# Edit user_inputs.in: use_song_rpe=1
# Run again

grep -E "SONG|RPE" outputs_original/converge.log > song_run.txt

# Compare
diff thermal_run.txt song_run.txt
```

---

## Expected Behavior

### Healthy Song Run:
1. ✅ `RPE MODEL: Song isothermal` at startup
2. ✅ `[RPE_MODEL_SONG]` when first parcel arrives
3. ✅ `[SONG_STEP]` messages (up to 10) showing increasing R and ε
4. ✅ `[SONG_COMPLETE]` when ε reaches 0.99

### Problem Indicators:
- ⚠️ Many `[SONG_ABORT]` → Parcels not superheated (check T and P conditions)
- ⚠️ `[SONG_COLLAPSE]` → Numerical instability or bad parameters
- ⚠️ `[SONG_ERROR]` → Bad initialization (check parcel injection)
- ⚠️ No `[RPE_MODEL_SONG]` → Parcels not entering thermal breakup regime

---

## Additional Diagnostics (If Needed)

Can add more detailed logging by modifying `RPE_song.c`:

### Pressure Term Breakdown:
```c
printf("[SONG_PRESSURE] P_sat=%.2e, P_init=%.2e, P_surface=%.2e, P_viscous=%.2e\n",
       P_sat, P_init, -2.0*sigma/R, viscous_term);
```

### Timestep Monitoring:
```c
printf("[SONG_TIMESTEP] dt_sub=%.3e s, n_sub=%d\n", dt_sub, n_sub);
```

### Temperature Tracking:
```c
printf("[SONG_TEMP] T_drop=%.2f K (should be constant)\n", T_drop);
```

Add after line 198 in `RPE_song.c` if more detail needed.

---

## Output File Locations

- **Main log:** `outputs_original/converge.log`
- **Echo file:** `user_inputs.echo` (confirms flag value)
- **Post files:** `outputs_original/output/post*.h5` (parcel data)

---

## Summary

**Implemented Debug Levels:**
- ✅ Model selection (startup)
- ✅ Activation confirmation
- ✅ Step-by-step logging (first 10 calls)
- ✅ Termination reasons
- ✅ Error conditions
- ✅ Edge cases

**Format:** All messages use `[TAG]` prefix for easy grepping

**Limitation:** Some messages limited to first 10 occurrences to avoid log spam

**Comparison:** Easy to compare Song vs Thermal using grep patterns

Ready for Phase 3 testing!
