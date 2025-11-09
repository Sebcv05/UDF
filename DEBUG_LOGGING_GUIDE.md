# RPE Debug Logging Usage Guide

## Overview
The RPE_euler solver now includes automatic single-parcel tracking based on parcel lifetime. It tracks the **oldest parcel in thermal breakup** and logs detailed solver state every 50 sub-timesteps.

## How It Works

1. **Automatic Selection**: Tracks the parcel with longest lifetime that enters thermal breakup (pbt = 1)
2. **Lightweight Logging**: Only logs every 50th sub-timestep → ~200-500 data points per parcel
3. **Unique Identification**: Uses `lifetime` field, which is reliable pre-breakup
4. **Small Files**: Typically <100 KB per simulation

## Files Generated

- **`rpe_parcel_debug.csv`**: Debug output (created in run directory)
  - Contains: time, lifetime, R, Rdot, T_drop, m_b, Nu, Q_conv, mdot, pressures, derivatives, properties

## Usage

### Step 1: Run Simulation
```bash
cd /path/to/case/directory
./upc2.sh      # Compile UDF
./run.sh       # Run CONVERGE
```

During the run, you'll see:
```
[RPE_DEBUG] Opened log file: rpe_parcel_debug.csv
[RPE_DEBUG] Tracking parcel with lifetime = 1.234567e-04 s
```

### Step 2: Visualize Results
```bash
cd /path/to/case/directory
python /home/apollo19/Desktop/Dan_B/UDF/plot_rpe_debug.py
```

Or copy the script to your case directory:
```bash
cp /home/apollo19/Desktop/Dan_B/UDF/plot_rpe_debug.py .
python plot_rpe_debug.py
```

### Step 3: Analyze Output

The plot shows 12 subplots:

**Row 1: Basic Dynamics**
- Bubble radius R(t) vs droplet radius Ro(t)
- Bubble wall velocity Rdot(t)
- Bubble acceleration dRdot/dt(t)

**Row 2: Thermal**
- Pressures: Pb, P_sat, P_amb
- Droplet temperature T_drop(t)
- Cooling rate dT/dt(t)

**Row 3: Heat/Mass Transfer**
- Nusselt number Nu(t)
- Convective heat transfer Q_conv(t)
- Evaporation rate ṁ(t)

**Row 4: Diagnostics**
- Phase portrait: Rdot vs R
- Bubble vapor mass m_b(t)
- Driving pressure: Pb - P_amb

## Configuration

Edit `src/RPE_euler.c` (lines 16-18):

```c
#define RPE_DEBUG_LOGGING 1        // Set to 0 to disable
#define LOG_INTERVAL 50            // Log every Nth call (increase to reduce file size)
#define LOG_FILE_NAME "rpe_parcel_debug.csv"
```

## Validation Checks

### ✅ Expected Behavior:
- R grows monotonically (no oscillations)
- Rdot > 0 throughout growth
- T_drop decreases slightly (1-5 K cooling)
- Nu ∈ [2, 1000]
- Pb > P_amb (driving force)
- No NaN or negative values

### 🚨 Red Flags:
- R decreasing → check pressure difference
- Rdot < 0 → bubble collapsing (wrong physics)
- T_drop increasing → energy balance wrong
- Nu > 1000 → hitting limit (may need adjustment)
- Large temperature swings (>10 K) → timestep too large

## Example: Validation Test

**Expected for T_init = 323 K, P_amb = 2 bar, Ro = 82.5 µm:**

| Quantity | Expected Value |
|----------|---------------|
| Final R | ~80 µm (fills droplet) |
| Max Rdot | 15-20 m/s |
| ΔT_drop | 1-5 K (cooling) |
| Time to fill | 5-10 µs |
| Max Nu | 50-200 |

## Troubleshooting

**No log file created:**
- Check RPE_DEBUG_LOGGING = 1
- Verify parcels enter thermal breakup (pbt = 1)
- Check simulation runs long enough

**File is huge (>10 MB):**
- Increase LOG_INTERVAL (try 100 or 200)
- Check if multiple parcels are tracked (shouldn't happen)

**Plot script fails:**
- Install pandas and matplotlib: `pip install pandas matplotlib`
- Check CSV file exists and has data
- Verify Python 3 is available

**Values look wrong:**
- Check property inputs (k_l, Nu_max in RPE_euler.c)
- Verify Antoine coefficients for NH3
- Compare with standalone Python/MATLAB solver

## Disabling Debug Logging

For production runs (no logging overhead):

```c
// In src/RPE_euler.c, line 16:
#define RPE_DEBUG_LOGGING 0
```

Recompile with `./upc2.sh`

## Data Format

CSV columns:
```
time          - Simulation time (s)
lifetime      - Parcel lifetime (s)
R             - Bubble radius (m)
Rdot          - Bubble wall velocity (m/s)
T_drop        - Droplet temperature (K)
m_b           - Bubble vapor mass (kg)
Nu            - Nusselt number (-)
Q_conv        - Convective heat transfer (W)
mdot          - Mass transfer rate (kg/s)
Pb            - Bubble pressure (Pa)
P_sat         - Saturation pressure (Pa)
P_amb         - Ambient pressure (Pa)
Ro            - Droplet radius (m)
dRdt          - dR/dt (m/s)
dRdotdt       - d²R/dt² (m/s²)
dTdt          - dT/dt (K/s)
dmbdt         - dm_b/dt (kg/s)
rho_l         - Liquid density (kg/m³)
mu_l          - Liquid viscosity (Pa·s)
sigma         - Surface tension (N/m)
k_l           - Thermal conductivity (W/m/K)
cp_l          - Specific heat (J/kg/K)
L_v           - Latent heat (J/kg)
```

## Next Steps

1. Run simulation and generate debug file
2. Plot results with Python script
3. Validate against expected behavior
4. Adjust parameters if needed (k_l, Nu_max)
5. Compare with old Bubble_Velocity results
6. Disable logging for production runs
