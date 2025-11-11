# CORRECTED ANALYSIS: Early Parcel Breakup Issue

## Date: 2025-11-11 (Revision)

## The Real Problem

After deeper analysis, the issue is **NOT** that the 1.5x check was blocking breakup.

### What's Actually Happening

**Tracked Parcel (p_idx=0):**
- Has 374,642 data points (extensive tracking)
- `thermal_breakup_flag` remains `-1` throughout
- kb values: 0.058 → 0.076 (growing slowly)
- **kb_threshold = 1.0** (from Breakup.c line 25)
- kb never approaches threshold!

### The 1.5x Expansion Check - Corrected Understanding

```c
Line 411: if (thermal_breakup_flag < 0 && pbt == 1)  // Entry condition
{
   Line 423: if(radius > 1.5*r_drop_0)  // INSIDE the block
   {
      tbt = 1;  // Triggers breakup
      thermal_breakup_flag = 6;
      continue;  // Skip to next parcel
   }
   
   // If radius < 1.5x, continues here:
   ... RPE solver ...
   ... DGRE calculation ...
   ... kb calculation ...
   
   Line 571: if (kb > kb_threshold)  // kb threshold check
   {
      tbt = 1;  // Also triggers breakup
      thermal_breakup_flag = 3;
   }
}
```

**The 1.5x check provides an ALTERNATIVE breakup trigger:**
- Path A: radius > 1.5x → immediate breakup (thermal_breakup_flag=6)
- Path B: kb > 1.0 → physics-based breakup (thermal_breakup_flag=3)

**Both paths lead to breakup!**

### So Why Are Early Parcels Persisting?

Given that the tracked parcel has:
- kb_max ≈ 0.076 << 1.0 threshold
- radius likely < 1.5x (subcooled, slow expansion)

**Neither condition is met**, so the parcel correctly persists without breakup!

## Two Possible Scenarios

### Scenario 1: This is Correct Physics
- Early parcels injected into 293K, ambient 298K
- **Subcooled** - actually losing superheat
- Bubble grows slowly, doesn't drive sufficient instability
- kb remains < 1.0 indefinitely
- Parcel should NOT break up (not physically ready)

### Scenario 2: kb Calculation is Wrong
- kb should be growing faster
- DGRE or BreakupCriterion has a bug
- Parcel physically should break up, but kb incorrectly calculated

## Questions to Answer

### 1. What is the expected kb behavior?
- Should kb reach 1.0 for a 293K parcel in 298K ambient?
- What superheat is needed for kb to exceed threshold?

### 2. Is the 1.5x check masking a kb calculation bug?
If we remove the 1.5x check:
- Parcels that would have triggered via expansion now must wait for kb
- If kb calculation is broken, they'll never break up
- **Removing 1.5x might expose the real problem**

### 3. What superheat do later parcels have?
- Later parcels presumably break up successfully
- Are they breaking via 1.5x expansion OR kb threshold?
- Check their `thermal_breakup_flag` values:
  - `=6` → broke via 1.5x expansion
  - `=3` → broke via kb threshold

## Action Items

### Check thermal_breakup_flag distribution in actual simulation:
```bash
grep "thermal_breakup_flag" log | awk '{print $NF}' | sort | uniq -c
```

Expected results:
- Mostly `=3` → kb threshold working, 1.5x rarely triggered
- Mostly `=6` → kb calculation may be broken, relying on expansion

### If mostly `=6`:
**Problem:** kb calculation not working properly  
**Solution:** Debug DGRE_NH3 and BreakupCriterion, not remove 1.5x check

### If mostly `=3`:
**Problem:** Early parcels physically shouldn't break up (too cold)  
**Solution:** Either adjust injection conditions OR lower kb_threshold

## Revised Recommendation

**DO NOT remove the 1.5x check yet.**

Instead:
1. Check distribution of thermal_breakup_flag values (3 vs 6 vs -1)
2. Verify kb calculation is correct
3. Check if early parcels are physically able to reach kb=1.0
4. If kb calculation is correct, early parcels may simply be too cold to break up

## Physical Justification

**Thermal breakup requires superheat:**
- Tinj = 293K, Tamb = 298K
- Parcel is initially **subcooled** relative to ambient
- Takes time to heat up from gas phase
- Early parcels may not receive enough heat before exiting domain
- kb < 1.0 is physically correct for insufficiently superheated droplets

**The 1.5x check is a safety mechanism:**
- If bubble grows extremely large (slow heating, large bubble)
- Droplet expands to 1.5x original size
- Even if kb < 1.0, geometric instability warrants breakup
- This is actually **physically reasonable**

## Conclusion

The 1.5x expansion check is **not the bug**. It's a valid alternative breakup criterion.

The real question is: **Why is kb not reaching 1.0 for early parcels?**

Possible answers:
1. kb calculation is broken → Fix DGRE/BreakupCriterion
2. Parcels are too cold → Adjust injection temperature or kb_threshold
3. Parcels exit domain before breaking up → Extend domain or reduce velocity

**Removing the 1.5x check without understanding why kb is low could make things worse.**

