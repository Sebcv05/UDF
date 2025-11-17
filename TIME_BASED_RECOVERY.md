# Time-Based Bubble Collapse Recovery Implementation

**Date:** 2025-11-17  
**Status:** Ready for testing

---

## Overview

Replaced the cycle-based recovery counter (100 attempts) with a **time-based recovery system** using a 5e-5 second (50 microsecond) delay between recovery attempts.

---

## Motivation

The original cycle-based approach (counting to 100) wasn't working correctly - only one recovery attempt was observed instead of multiple. A time-based approach is more physically meaningful and avoids issues with sub-cycling loops.

---

## Implementation Details

### New Custom Parcel Fields

Added two new fields to the parcel cloud structure:

1. **`recovery_time`** (CONVERGE_DOUBLE)
   - Stores the simulation time (in seconds) when the last recovery occurred
   - Initialized to 0.0 for new parcels (both parent and child)
   - Reset to 0.0 when bubble recovers successfully

2. **`recovery_count`** (CONVERGE_INT)
   - Tracks the number of recovery attempts for debugging/monitoring
   - Initialized to 0 for new parcels
   - Incremented on each recovery attempt
   - NOT reset when recovery succeeds (kept for history tracking)

### Recovery Logic (RPE_euler.c)

When `Rdot < 0` (bubble collapsing):

1. **Check if in recovery period:**
   ```c
   time_since_recovery = current_time - recovery_time
   if (recovery_time > 0 && time_since_recovery < RECOVERY_PERIOD) {
       // Still waiting - keep bubble stable, return early
   }
   ```

2. **If not in recovery period, attempt recovery:**
   - Reset bubble to 1.1 × R_critical
   - Reset droplet to r_drop_0
   - Conserve mass by adjusting num_drop
   - **Record recovery time: `recovery_time = current_time`**
   - Increment `recovery_count`
   - Set flag=888 to break sub-loop

3. **On successful recovery (Rdot > 0):**
   - Reset `recovery_time = 0.0`
   - Keep `recovery_count` for history

### Recovery Period

```c
const CONVERGE_precision_t RECOVERY_PERIOD = 5.0e-5;  // 50 microseconds
```

This value can be easily adjusted if needed.

---

## Files Modified

### 1. `include/lagrangian/env.h`
- Added field ID declarations:
  ```c
  CONVERGE_id_t RECOVERY_TIME;
  CONVERGE_id_t RECOVERY_COUNT;
  ```
- Added pointers to ParcelCloud struct:
  ```c
  CONVERGE_precision_t* recovery_time;
  CONVERGE_int_t* recovery_count;
  ```

### 2. `src/load_spray_env.c`
- Registered new fields with CONVERGE:
  ```c
  CONVERGE_variable_register("recovery_time", CONVERGE_DOUBLE, ...);
  CONVERGE_variable_register("recovery_count", CONVERGE_INT, ...);
  ```
- Retrieved field IDs:
  ```c
  RECOVERY_TIME = CONVERGE_lagrangian_field_id("recovery_time");
  RECOVERY_COUNT = CONVERGE_lagrangian_field_id("recovery_count");
  ```
- **Added to `load_user_cloud()` function:**
  ```c
  parcel_cloud_loc->recovery_time = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, RECOVERY_TIME);
  parcel_cloud_loc->recovery_count = (CONVERGE_int_t *)CONVERGE_cloud_get_field_data(c, RECOVERY_COUNT);
  ```

### 3. `src/RPE_euler.c`
- Replaced `user_lag_var_v3[0]` with `recovery_time`
- Replaced `user_lag_var_i` increment with `recovery_count`
- Updated all logging messages
- Logic remains the same, just uses proper fields now

### 4. `src/parcel_prop.c`
- Initialize `recovery_time = 0.0` in `parcel_inject()`
- Initialize `recovery_count = 0` in `parcel_inject()`
- Initialize both fields in `parcel_child()`

---

## Advantages Over Cycle-Based Approach

1. **Physically meaningful:** Time delay is independent of sub-cycling/mesh resolution
2. **Simpler logic:** No need to track cycles or worry about sub-loop behavior
3. **More robust:** Works regardless of timestep structure
4. **Clearer code:** Custom fields with descriptive names instead of hijacking `user_lag_var_*`

---

## Expected Behavior

### First Collapse
```
[RPE_COLLAPSE] Negative Rdot=-7.470e-03, attempting recovery (bubble collapsing)
               Time: 1.234567e-05 s, Last recovery: 0.000000e+00 s (1.234567e-05 s ago)
               [RECOVERY #1] R_bubble: 3.812e-05 -> 1.191e-07 m (1.1*Rc=1.083e-07)
               [RECOVERY #1] Recovery time logged: 1.234567e-05 s, wait 5.000e-05 s before next
```

### During Recovery Wait Period (< 50 μs later)
```
[RPE_RECOVERY_WAIT] Still in recovery period: 1.234e-05 / 5.000e-05 s elapsed
                    Rdot=-2.345e-03 (negative but waiting), R=1.191e-07 m
```

### After Recovery Period (> 50 μs later, still collapsing)
```
[RPE_COLLAPSE] Negative Rdot=-1.234e-03, attempting recovery (bubble collapsing)
               Time: 6.789012e-05 s, Last recovery: 1.234567e-05 s (5.554445e-05 s ago)
               [RECOVERY #2] R_bubble: 1.191e-07 -> 1.191e-07 m (1.1*Rc=1.083e-07)
               [RECOVERY #2] Recovery time logged: 6.789012e-05 s, wait 5.000e-05 s before next
```

### Successful Recovery
```
[RPE_RECOVERY_SUCCESS] Bubble recovered! Rdot=1.234e-02 m/s, 5.678e-05 s after recovery
                       Resetting recovery timer (was 1.234567e-05 s)
```

---

## Testing Checklist

- [ ] Compile successfully with `upc2.sh`
- [ ] Check log for multiple recovery attempts on rogue parcels
- [ ] Verify time-based spacing (should see attempts ~50 μs apart)
- [ ] Confirm normal parcels (99%) unaffected
- [ ] Check that `recovery_count` increments properly
- [ ] Verify successful recovery resets `recovery_time` to 0

---

## Tuning Parameters

If 50 μs is too short/long, adjust:
```c
const CONVERGE_precision_t RECOVERY_PERIOD = 5.0e-5;  // seconds
```

Recommended range: 1e-5 to 1e-4 seconds (10-100 μs)

---

## Notes

- Old fields (`user_lag_var_i`, `user_lag_var_v3`) are still available but not used for recovery
- The counter-based code has been completely replaced
- No "give up after X attempts" logic - will keep trying indefinitely at 50 μs intervals
- If you want to add a maximum limit, you could check `recovery_count > MAX_ATTEMPTS`

---

## Compatibility

- Uses standard CONVERGE API for custom fields
- No hijacking of built-in fields
- Clean, readable code
- Easy to modify or extend in the future
