# FINAL: Safe Variables to Hijack

## The Safe Choice: KH/RT Breakup Timers

Since you're using **thermal breakup** (not KH or RT), these timer variables are SAFE:

### Final Mapping:

| Diagnostic | Hijacked Variable | Original Purpose | Why Safe |
|------------|------------------|------------------|----------|
| **is_child** | `tbreak_kh` | KH breakup timer | You're not using KH breakup |
| **r_bubble** | `tbreak_rt` | RT breakup timer | You're not using RT breakup |
| **pbt** | `t_turb` | Turbulent breakup time | Not used in thermal model |

## Why These Are Safe

1. **`tbreak_kh`** - Kelvin-Helmholtz breakup timer
   - Only active if KH breakup model is enabled
   - You're using thermal breakup instead
   - **100% safe to hijack**

2. **`tbreak_rt`** - Rayleigh-Taylor breakup timer
   - Only active if RT breakup model is enabled
   - You're using thermal breakup instead
   - **100% safe to hijack**

3. **`t_turb`** - Turbulent breakup accumulation time
   - Not used in your thermal breakup physics
   - **Safe to hijack**

## In Tecplot

After recompiling:

### `tbreak_kh` = is_child flag
- **0** = Parent parcel (injected)
- **1** = Child parcel (from breakup)

**Filter for large parent parcels:**
```
radius > 7e-5 AND tbreak_kh = 0
```

**Filter for large child parcels:**
```
radius > 7e-5 AND tbreak_kh = 1
```

### `tbreak_rt` = Bubble radius (meters)
- **0** = No bubble
- **1e-9 to 8e-5** = Active bubble

**Filter for parcels about to break:**
```
tbreak_rt > 7e-5
```

### `t_turb` = pbt flag
- **0** = Not breaking up
- **1** = Currently in breakup

---

## Files Modified

✅ `/home/apollo19/Desktop/Dan_B/UDF/src/parcel_prop.c` - Uses tbreak_kh, tbreak_rt, t_turb  
✅ `/home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev/post.in` - Outputs these 3 variables

## Ready to Compile

```bash
cd [CASE_DIR]
./upc2.sh
./run.sh
```

---

**Status:** ✅ These variables exist in ParcelCloud struct  
**Status:** ✅ Should compile successfully  
**Date:** 2025-11-13
