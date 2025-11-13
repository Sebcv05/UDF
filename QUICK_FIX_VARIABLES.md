# Fixed: Use Existing ParcelCloud Variables

## Problem
`radiation_energy` and `distant` don't exist in ParcelCloud struct.

## Solution
Use variables that **DO exist** and are already output:

### New Mapping:

| Diagnostic | Hijacked Variable | Original Purpose | New Use |
|------------|------------------|------------------|---------|
| **is_child** | `distort` | TAB distortion | 0=parent, 1=child |
| **r_bubble** | `distort_dot` | TAB distortion rate | Bubble radius (m) |
| **pbt** | `t_turb` | Turbulent breakup time | 0=no breakup, 1=in breakup |

## Why These Variables?

1. **`distort`** - Used for TAB (Taylor Analogy Breakup), which you're NOT using for thermal breakup
2. **`distort_dot`** - Rate of distortion change, also TAB-specific
3. **`t_turb`** - Turbulent breakup time, NOT used in your thermal breakup model

All three are safe to hijack!

## In Tecplot

After compilation, you'll see:

### `distort`
- 0 = Parent parcel
- 1 = Child parcel

**Filter for large parent parcels:**
```
radius > 7e-5 AND distort = 0
```

### `distort_dot`
- Bubble radius in meters (0 to 8e-5)

**Filter for parcels about to break:**
```
distort_dot > 7e-5
```

### `t_turb`
- 0 = Not in breakup
- 1 = Currently breaking up

---

## Files Modified

✅ `/home/apollo19/Desktop/Dan_B/UDF/src/parcel_prop.c`  
✅ `/home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev/post.in`  
✅ Documentation files

## Ready to Compile

```bash
cd [CASE_DIR]
./upc2.sh
./run.sh
```

**Status:** ✅ Should compile now!
