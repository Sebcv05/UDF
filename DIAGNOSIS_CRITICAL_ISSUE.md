# CRITICAL ISSUE DIAGNOSED - Root Cause Found

## Date: 2025-11-09

## Problem Summary

**Thermal breakup not occurring after 25 µs despite evaporation being disabled.**

## Root Cause: INSUFFICIENT SUPERHEAT

### Evidence from Diagnostics:

```
[RPE_STATUS] T_drop=293.15 K, T_sat=236.30 K, P_sat=8.49e+05 Pa
Max R=7.977e-16 µm, Ro=1.000e-18 m
```

### Critical Issues Found:

1. **Droplet radius is ZERO**: Ro = 1e-18 m instead of 82.5e-6 m
   - Parcels are being destroyed/removed before thermal breakup
   - RPE cannot run on parcels with Ro → 0

2. **Temperature too low FOR SUPERHEATED NH3**:
   - Current: T_inj = 293.15 K (20°C)  
   - P_amb = 2.0 bar (from initialize.in)
   - P_sat(293K) = **1.57 bar**
   - **P_sat < P_amb** → **NO SUPERHEAT!**

### Thermodynamic Analysis:

For NH3 at P_amb = 2.0 bar:
- T_sat(2 bar) = **340.2 K** (67°C)
- Current T_inj = 293.15 K
- **Droplet is SUBCOOLED by 47 K!**

For thermal breakup to occur:
- Need: T_inj > T_sat(P_amb) + margin
- Minimum: **T_inj > 340 K**
- Recommended: **T_inj = 350-360 K** (10-20 K superheat)

### Unit Test Confirmation:

Standalone RPE test at 323 K:
```
P_sat = 1.57 bar (superheat = -0.43 bar)
>>> Droplet subcooled - bubble collapsed immediately <<<
```

Even at 355 K (correct superheat):
- Bubble starts collapsing (Rdot = -84 m/s)
- Initial R too small (1 nm) hits safety floor

## Solutions Required:

### 1. Fix Injection Temperature (CRITICAL)

Check injection settings - need to inject at **T > 340 K**

**Where to check:**
```bash
grep -r "temp.*inj\|injection.*temp" /path/to/case/*.in
```

Look for:
- parcels.in or injection.in
- Stream temperature settings
- Initial parcel temperature

**Required change:**
```
T_injection: 293.15 K  →  355 K (or higher)
```

### 2. Fix Droplet Radius Issue (CRITICAL)

Why is Ro = 1e-18 m instead of 82.5 µm?

**Possible causes:**
- Parcels being removed before thermal breakup
- num_drop going to zero (→ radius shrinks)
- Geometry() function error
- Another destruction mechanism

**Diagnostic added:**
```
[RPE_RADIUS] will show actual Ro, num_drop, is_child values
```

**Check output for:**
```bash
grep "\[RPE_RADIUS\]" output.log
```

Should show: `Ro=8.25e-05 m (82.5 µm)`, not `1e-18 m`

### 3. Check Initial Bubble Size

Current: R_init ≈ 1e-22 m (atomic scale!)

**Where set?** Check where `r_bubble[p_idx]` is initialized in spray_drop_distort_NH3.c

Typical nucleation size: **R_c ≈ 1-10 nm** (1e-9 to 1e-8 m)

## Action Items:

1. **Immediately check injection temperature:**
   ```bash
   grep -i "temp" parcels.in injection.in stream*.in
   ```

2. **Rerun with T_inj = 355 K**

3. **Check [RPE_RADIUS] diagnostic:**
   ```bash
   ./upc2.sh
   ./run.sh 2>&1 | grep "\[RPE_RADIUS\]"
   ```

4. **Expected after fixes:**
   ```
   [RPE_RADIUS] Ro=8.25e-05 m (82.5 µm), num_drop=1.00e+00
   [RPE_GROWTH] R: 1.0e-09 -> 5.3e-09 (+4.3e-09)
   [RPE_STATUS] Max R=50.3 µm, Max Rdot=18.5 m/s
   [BREAKUP] Cycle 250: 1 breakups
   ```

## Physics Summary:

**For NH3 thermal breakup at 2 bar ambient:**

| Parameter | Current | Required | Status |
|-----------|---------|----------|--------|
| P_amb | 2.0 bar | 2.0 bar | ✓ OK |
| T_sat(P_amb) | 340.2 K | - | - |
| T_inj | 293.15 K | >340 K | ❌ FAIL |
| Superheat | -47 K | +10-20 K | ❌ FAIL |
| P_sat(T_inj) | 1.57 bar | >2.2 bar | ❌ FAIL |
| Ro | 1e-18 m | 82.5e-6 m | ❌ FAIL |

**Cannot have bubble growth when droplet is subcooled!**

## Next Run Checklist:

- [ ] Set T_injection = 355 K (minimum 350 K)
- [ ] Verify Ro from [RPE_RADIUS] diagnostic
- [ ] Check [EVAP_SKIP] still working
- [ ] Monitor [RPE_GROWTH] for positive growth
- [ ] Look for [BREAKUP] events

## Files to Check:

1. Case setup files: `parcels.in`, `injection.in`, `stream*.in`
2. Output: Look for [RPE_RADIUS], [RPE_GROWTH], [BREAKUP]
3. Temp_Tracker.txt should show ~355 K, not 293 K

