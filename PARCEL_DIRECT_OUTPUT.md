# Better Solution: Direct Parcel CSV Output

## Problem with Cell-Centered Approach

The cell-centered CONVERGE_POST approach has limitations:
- Only shows data for **first parcel** in each cell
- Multiple parcels in same cell → lose data for all but one
- No direct parcel-to-data mapping
- Can't distinguish individual parcels

## Solution: Write CSV Files Directly

Write parcel data directly to CSV files during simulation using `CONVERGE_OUTPUT`.

---

## Implementation Option 1: Timestep-Based CSV (Simple)

Add to your UDF a function that writes all parcel data to CSV at each output timestep.

**File:** `/home/apollo19/Desktop/Dan_B/UDF/src/parcel_output.c` (already created)

This creates files like:
```
parcel_diagnostics_0.000000e+00.csv
parcel_diagnostics_5.877070e-06.csv
parcel_diagnostics_1.147900e-05.csv
...
```

Each file contains ALL parcels at that timestep with columns:
- time, parcel_id, cell_id, x, y, z
- radius, temp, velocity_mag, num_drop
- **is_child**, **r_bubble**, **pbt**, from_injector

---

## Implementation Option 2: Track Specific Parcels (Targeted)

Only write data for parcels that meet certain criteria (e.g., large parcels).

**Add to your existing UDF:**

```c
// In spray_drop_distort_NH3.c or similar, during parcel loop:

// Track large parcels
if(parcel_cloud.radius[p_idx] > 7.0e-5) {  // > 70 µm
   
   // Write to diagnostic file (append mode)
   FILE* fp = fopen("large_parcels.csv", "a");
   
   fprintf(fp, "%.6e,%d,%.6e,%.6e,%d,%.6e,%d\n",
      CONVERGE_simulation_time(),
      p_idx,
      parcel_cloud.radius[p_idx],
      parcel_cloud.r_bubble[p_idx],
      parcel_cloud.is_child[p_idx],
      parcel_cloud.temp[p_idx],
      parcel_cloud.pbt[p_idx]
   );
   
   fclose(fp);
}
```

This creates a single file with only "interesting" parcels.

---

## Implementation Option 3: Parcel Tracking Array

Use CONVERGE's parcel tracking to follow specific parcels through time.

**Already exists in your code:**
- `parcel_cloud.from_injector[p_idx]` - injection event ID
- `parcel_cloud.is_child[p_idx]` - parent/child flag
- `parcel_cloud.p_idx` - parcel index

**Strategy:**
1. Identify mystery large parcel in Tecplot (note its position/time)
2. Add tracking when parcel created:
   ```c
   // In Breakup.c when creating children:
   static int breakup_event_counter = 0;
   
   for(int child_i = 0; child_i < num_children; child_i++) {
      // Tag this child with unique ID
      new_parcel_cloud.breakup_event_id[new_idx] = breakup_event_counter;
      new_parcel_cloud.child_number[new_idx] = child_i;
   }
   
   breakup_event_counter++;
   ```

3. Write CSV with tracking info
4. Post-process to follow individual parcels

---

## Quick Implementation (Recommended)

**Add CSV output to existing code with minimal changes:**

### Step 1: Add to spray_drop_distort_NH3.c

At the end of the main loop (after all parcel operations), add:

```c
// Once per output timestep
static CONVERGE_precision_t last_output_time = -1.0;
CONVERGE_precision_t current_time = CONVERGE_simulation_time();
CONVERGE_precision_t output_interval = 5.0e-6;  // Output every 5 µs

if(current_time - last_output_time >= output_interval) {
   
   int rank;
   CONVERGE_mpi_comm_rank(&rank);
   if(rank == 0) {  // Only rank 0 writes
      
      char filename[256];
      snprintf(filename, sizeof(filename), "parcels_t%.6e.csv", current_time);
      FILE* fp = fopen(filename, "w");
      
      fprintf(fp, "time,cell_id,x,y,z,radius,is_child,r_bubble,pbt,num_drop\n");
      
      // Loop through all parcels (use existing cloud loop)
      // Write each parcel's data
      
      fclose(fp);
      printf("[PARCEL_CSV] Wrote %s\n", filename);
   }
   
   last_output_time = current_time;
}
```

### Step 2: Post-Process CSV Files

```python
import pandas as pd
import glob

# Load all CSV files
csv_files = sorted(glob.glob('parcels_t*.csv'))
dfs = []

for file in csv_files:
    df = pd.read_csv(file)
    dfs.append(df)

all_parcels = pd.concat(dfs)

# Find large parcels
large_parcels = all_parcels[all_parcels['radius'] > 7e-5]

# Check if they're parents or children
print("Large parcels breakdown:")
print(large_parcels.groupby('is_child').size())

# Find specific large parcel at specific time/location
mystery_parcel = large_parcels[
    (large_parcels['time'] > 2.5e-5) & 
    (large_parcels['x'] > 0.01) &
    (large_parcels['radius'] > 7.5e-5)
]

print("\nMystery large parcel:")
print(mystery_parcel[['time', 'radius', 'is_child', 'r_bubble', 'pbt']])
```

---

## Comparison: H5 vs CSV Approaches

| Feature | H5 (cell-centered) | CSV (direct write) |
|---------|-------------------|-------------------|
| Per-parcel data | ❌ No (first only) | ✅ Yes (all parcels) |
| File size | Small | Larger |
| Tecplot integration | ✅ Native | ❌ Need post-process |
| Multiple parcels/cell | ❌ Lost | ✅ All captured |
| Real-time writing | ❌ No | ✅ Yes |
| Easy filtering | ✅ In Tecplot | ⚠️ In Python/script |
| Parcel tracking | ❌ Difficult | ✅ Easy |

---

## Recommendation

**Use both approaches:**

1. **H5 cell-centered** (current approach) - for Tecplot visualization
   - Good for general flow field + parcel distribution
   - Shows where parcels are

2. **CSV direct write** - for detailed parcel diagnostics
   - Write only at key timesteps (every 5-10 µs)
   - Captures ALL parcel data
   - Easy post-processing to answer specific questions

---

## Minimal CSV Implementation (Add to existing UDF)

Create `/home/apollo19/Desktop/Dan_B/UDF/src/parcel_output.c`:

Already created! Contains `CONVERGE_OUTPUT(write_parcel_data, ...)` function.

**To activate:**
1. It's automatically called at each CONVERGE output timestep
2. Creates `parcel_diagnostics_XXXXXX.csv` files
3. Contains ALL parcel data including is_child, r_bubble, pbt

**Usage:**
```bash
# After simulation:
ls parcel_diagnostics_*.csv

# Load in Python:
import pandas as pd
df = pd.read_csv('parcel_diagnostics_3.005860e-05.csv')

# Find large parcels:
large = df[df['radius'] > 7e-5]
print("Large parent parcels:", len(large[large['is_child'] == 0]))
print("Large child parcels:", len(large[large['is_child'] == 1]))
```

---

**Created:** 2025-11-13  
**Recommended:** Add parcel_output.c to your UDF for complete parcel tracking
