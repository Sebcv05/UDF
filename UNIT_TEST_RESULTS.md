# RPE Unit Test Results

## Date: 2025-11-09

## Test Conditions

Using **CORRECT Antoine coefficients** (Stull 1947):
- Temperature: 293.15 K (20°C)
- P_amb: 2.0 bar
- P_sat(293.15 K): **8.49 bar**
- **Superheat: 6.49 bar** ✓

## Results

**Bubble grew successfully from 1 µm to 74 µm in 2.78 µs**

| Metric | Value |
|--------|-------|
| Initial R | 1.0 µm |
| Final R | 74.25 µm |
| Growth | 73.25 µm |
| Max Rdot | 26.6 m/s |
| Time to fill | 2.78 µs |
| Final R/Ro | 90% |

## Conclusion

✅ **RPE physics is CORRECT**
✅ **Injection temperature (293 K) is CORRECT** - provides 6.49 bar superheat
✅ **Unit test shows rapid bubble growth** as expected

## Problem Identified in CONVERGE

From diagnostic output:
```
[RPE_STATUS] Max R=7.977e-16 µm, Ro=1.000e-18 m
```

❌ **Droplet radius is essentially ZERO** (Ro = 1e-18 m instead of 82.5e-6 m)

This prevents RPE from running:
- Bubble hits safety check immediately (R > Ro)
- No growth possible when droplet has no size

## Root Cause

Parcels are being **destroyed or having radius set to zero** before thermal breakup can occur.

Possible causes:
1. `num_drop[p_idx]` going to zero
2. Parcels being removed in spray_drop_distort_NH3.c
3. Geometry() function error
4. KH-RT or other breakup removing parcels
5. Mass conservation issue

## Next Steps

1. **Check [RPE_RADIUS] diagnostic:**
   ```bash
   grep "\[RPE_RADIUS\]" output.log
   ```
   Should show: `Ro=8.25e-05 m (82.5 µm)`, not `1e-18 m`

2. **Check if num_drop is zero:**
   - If num_drop → 0, then radius → 0 (mass conservation)

3. **Check for parcel removal:**
   - Search spray_drop_distort_NH3.c for where `radius[p_idx] = 0` or `num_drop[p_idx] = 0`
   - Check if temperature > 300 K removal is triggering (line 281)

4. **Verify parcels survive to thermal breakup:**
   - Check parcel lifetime values
   - Confirm is_child = 0 (parent parcels)

## Files Generated

- `bubble_growth.csv` - Raw data
- `rpe_unit_test_CORRECT.pdf` - Visualization
- `test_correct` - Standalone executable

## Test Can Be Run Independently

```bash
cd /home/apollo19/Desktop/Dan_B/UDF
./test_correct
# Output: bubble_growth.csv
```

No CONVERGE dependencies - pure C test of RPE physics.

