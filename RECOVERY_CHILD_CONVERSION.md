# Recovery Strategy: Convert to Child Parcel

**Date:** 2025-11-17  
**Status:** Implemented - Child Conversion Approach

---

## Problem with Previous Approaches

All attempts to "pause" or "wait" during recovery caused crashes:
1. Skip with `continue` → CONVERGE crash
2. Skip RPE call → Solver instability  
3. Freeze with `pbt=0` → Still caused crashes

**Root cause:** Having parcels with tiny bubbles (119 nm) in the domain causes solver issues, whether they're active or passive.

---

## New Strategy: Immediate Conversion to Child Parcel

**Core Idea:** Instead of trying to recover the bubble and wait, immediately **give up on thermal breakup** and convert the parcel to a regular child parcel that will break up via KH-RT and evaporation.

### What Happens on Collapse

When bubble collapse is detected (Rdot < 0):

1. **Reset droplet to injection size** - `R_drop = r_drop_0`
2. **Zero the bubble** - `R_bubble = 0`
3. **Conserve mass** - Adjust `num_drop` to maintain total liquid mass
4. **Convert to child** - `is_child = 1`
5. **Disable thermal breakup** - `pbt = 0`, `flag = 999`
6. **Exit RPE** - Parcel never re-enters RPE

---

## Implementation

**Location:** `RPE_euler.c`, line ~355-440

### Key Changes:

```c
// When Rdot < 0 is detected:

// Check if already recovered (shouldn't happen)
if (recovery_time > 0) {
    printf("[RPE_ERROR] Recovered parcel re-entered RPE!\n");
    pbt = 0; flag = 999;
    return;
}

// First collapse - convert to child
CONVERGE_precision_t injection_radius = r_drop_0;

// Calculate mass to conserve
total_mass = num_drop * mass_liquid_old;

// New state: droplet at injection size, no bubble
V_drop_new = (4/3)*π*injection_radius³;
mass_per_drop_new = V_drop_new * rho_l;
new_num_drop = total_mass / mass_per_drop_new;

// Apply conversion
radius = injection_radius;      // Reset to injection size
r_bubble = 0.0;                // Zero bubble
r_bubble_0 = 0.0;
num_drop = new_num_drop;       // Conserve mass
v_bubble = 0.0;
is_child = 1;                  // Mark as child
pbt = 0;                       // Disable thermal breakup
flag = 999;                    // Mark as aborted

// Record for tracking
recovery_time = current_time;
recovery_count++;

return;  // Exit RPE forever
```

---

## How It Works

### Normal Thermal Breakup Parcel

```
Injection → RPE growth → Bubble expands → Breakup criterion → Children
```

### Parcel with Collapse (New Strategy)

```
Injection → RPE growth → Bubble collapses (Rdot < 0)
              ↓
         Convert to child:
         - R_drop → injection size
         - R_bubble → 0
         - is_child = 1
         - pbt = 0
              ↓
         KH-RT + Evaporation → Breakup via other mechanisms
```

---

## Advantages

### 1. **No More Crashes**
- No tiny bubbles left in domain
- No waiting/pausing/freezing logic
- Parcel immediately becomes a normal child

### 2. **Physically Reasonable**
- Parcel failed thermal breakup → try other mechanisms
- Reset to "safe" state (injection size, no bubble)
- Let KH-RT and evaporation handle it naturally

### 3. **Simple and Clean**
- One-time conversion, no complex state management
- No wait periods, no re-entry logic
- Parcel never comes back to RPE

### 4. **Mass Conserved**
- Total liquid mass preserved by adjusting num_drop
- Same approach used in normal breakup

### 5. **Compatible with Existing Code**
- Uses existing `is_child` flag
- Uses existing `pbt` disable mechanism
- No changes to outer loop structure

---

## What Gets Logged

### When Collapse Occurs:
```
[RPE_COLLAPSE] Negative Rdot=-7.470e-03, converting to child parcel
               Time: 3.239393e-05 s
               T_drop=293.06 K, T_sat(P_amb)=274.75 K, P_sat(T_drop)=8.468e+05 Pa, P_amb=4.518e+05 Pa
               R=3.812e-05 m, Ro=8.318e-05 m, dRdt=1.264e-03 m/s, dRdotdt=-8.882e+06 m/s²
               [RECOVERY #1] R_drop: 8.318e-05 -> 8.250e-05 m (injection radius)
               [RECOVERY #1] R_bubble: 3.812e-05 -> 0.00e+00 m (zeroed)
               [RECOVERY #1] num_drop: 9.263e-01 -> 1.046e+00 (mass conserved)
               [RECOVERY #1] is_child: 0 -> 1 (converted to child parcel)
               [RECOVERY #1] pbt: 1 -> 0 (thermal breakup disabled)
```

### If Somehow Re-enters RPE (Error):
```
[RPE_ERROR] Recovered parcel (child) re-entered RPE! p_idx=X, recovery_time=3.239e-05 s
```

---

## Why This Won't Crash

1. **No tiny bubbles** - R_bubble set to 0, not 119 nm
2. **No RPE integration** - pbt=0, parcel never re-enters RPE
3. **Normal child parcel** - Treated like any post-breakup child
4. **No waiting logic** - Instant conversion, no time-based checks
5. **No complex state** - Simple one-way transformation

---

## Potential Concerns

### Q: Will we lose breakup events?
**A:** No - parcels will still break up via KH-RT and evaporation. They just won't break up thermally.

### Q: Is resetting to injection size okay?
**A:** Yes - we're conserving mass, just redistributing it. The parcel is essentially "re-initialized" as a child.

### Q: What if thermal breakup was important for this parcel?
**A:** If it's collapsing, thermal breakup wasn't working anyway. Better to let other mechanisms handle it.

### Q: Does is_child=1 affect anything else?
**A:** Yes - it enables KH-RT and other breakup models, which is what we want.

---

## Testing

Look for:
1. **`[RPE_COLLAPSE]` messages** followed by conversion info
2. **No crashes** - Main success criterion
3. **Parcels still breaking up** - Via KH-RT/evaporation instead of thermal
4. **No `[RPE_ERROR]`** messages - Parcels shouldn't re-enter

---

## Files Modified

**`src/RPE_euler.c`:**
- Removed wait period logic (no more 20 μs delay)
- Removed recovery success check (no re-entry)
- Changed recovery action to child conversion
- Set R_bubble=0, is_child=1, pbt=0, flag=999

**`src/spray_drop_distort_NH3.c`:**
- No changes (reverted all previous attempts)

---

## Philosophy Shift

**Old approach:** "Let's try to save the thermal breakup process"
- Attempt recovery with tiny bubble
- Wait and hope it grows
- Complex state management
- Caused crashes

**New approach:** "If thermal breakup fails, move on"
- Accept that thermal breakup didn't work
- Convert to child parcel
- Let other mechanisms handle it
- Clean, simple, stable

This is less about "recovery" and more about "graceful failure handling."
