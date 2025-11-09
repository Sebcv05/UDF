#!/usr/bin/env python3
"""
Plot RPE_euler debug output for single parcel validation
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def plot_rpe_debug(csv_file='rpe_parcel_debug.csv'):
    """Load and plot RPE debug data"""
    
    if not os.path.exists(csv_file):
        print(f"Error: {csv_file} not found!")
        print("Make sure you ran CONVERGE with RPE_DEBUG_LOGGING=1")
        return
    
    # Load data
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} data points")
    print(f"Time range: {df.time.min():.6e} to {df.time.max():.6e} s")
    print(f"Lifetime: {df.lifetime.iloc[-1]:.6e} s")
    
    # Convert to convenient units
    df['time_us'] = df['time'] * 1e6  # microseconds
    df['R_um'] = df['R'] * 1e6        # micrometers
    df['Ro_um'] = df['Ro'] * 1e6
    df['Pb_bar'] = df['Pb'] * 1e-5    # bar
    df['Psat_bar'] = df['P_sat'] * 1e-5
    df['Pamb_bar'] = df['P_amb'] * 1e-5
    df['Q_conv_kW'] = df['Q_conv'] * 1e-3  # kW
    df['mdot_mg_s'] = df['mdot'] * 1e6     # mg/s
    
    # Create comprehensive plot
    fig, axes = plt.subplots(4, 3, figsize=(16, 12))
    fig.suptitle(f'RPE Solver Validation: Single Parcel (lifetime={df.lifetime.iloc[-1]:.2e}s)', 
                 fontsize=14, fontweight='bold')
    
    # Row 1: Basic dynamics
    # R and Ro vs time
    ax = axes[0, 0]
    ax.plot(df.time_us, df.R_um, 'b-', linewidth=2, label='R (bubble)')
    ax.plot(df.time_us, df.Ro_um, 'r--', linewidth=2, label='Ro (droplet)')
    ax.axhline(df.Ro_um.iloc[0], color='k', linestyle=':', alpha=0.3, label='Initial Ro')
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Radius (µm)')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_title('Bubble & Droplet Radius')
    
    # Rdot vs time
    ax = axes[0, 1]
    ax.plot(df.time_us, df.Rdot, 'g-', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Bubble Velocity (m/s)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Bubble Wall Velocity')
    
    # Acceleration vs time
    ax = axes[0, 2]
    ax.plot(df.time_us, df.dRdotdt, 'm-', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Acceleration (m/s²)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Bubble Acceleration')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    
    # Row 2: Pressures and temperature
    # Pressures
    ax = axes[1, 0]
    ax.plot(df.time_us, df.Pb_bar, 'b-', linewidth=2, label='Pb (bubble)')
    ax.plot(df.time_us, df.Psat_bar, 'r--', linewidth=2, label='P_sat')
    ax.plot(df.time_us, df.Pamb_bar, 'k:', linewidth=2, label='P_amb')
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Pressure (bar)')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_title('Pressure Evolution')
    
    # Temperature
    ax = axes[1, 1]
    ax.plot(df.time_us, df.T_drop, 'r-', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Temperature (K)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Droplet Temperature')
    T_init = df.T_drop.iloc[0]
    T_final = df.T_drop.iloc[-1]
    ax.text(0.05, 0.95, f'ΔT = {T_init-T_final:.2f} K', 
            transform=ax.transAxes, va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Temperature rate
    ax = axes[1, 2]
    ax.plot(df.time_us, df.dTdt, 'orange', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('dT/dt (K/s)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Cooling Rate')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    
    # Row 3: Heat and mass transfer
    # Nusselt number
    ax = axes[2, 0]
    ax.plot(df.time_us, df.Nu, 'purple', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Nusselt Number')
    ax.grid(True, alpha=0.3)
    ax.set_title('Nusselt Number (Nu)')
    ax.text(0.05, 0.95, f'Range: {df.Nu.min():.1f} - {df.Nu.max():.1f}', 
            transform=ax.transAxes, va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Heat transfer rate
    ax = axes[2, 1]
    ax.plot(df.time_us, df.Q_conv_kW, 'red', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Q_conv (kW)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Convective Heat Transfer')
    
    # Mass transfer rate
    ax = axes[2, 2]
    ax.plot(df.time_us, df.mdot_mg_s, 'blue', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('ṁ (mg/s)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Evaporation Rate')
    
    # Row 4: Phase portraits and diagnostics
    # Phase portrait: Rdot vs R
    ax = axes[3, 0]
    ax.plot(df.R_um, df.Rdot, 'b-', linewidth=2)
    ax.plot(df.R_um.iloc[0], df.Rdot.iloc[0], 'go', markersize=10, label='Start')
    ax.plot(df.R_um.iloc[-1], df.Rdot.iloc[-1], 'ro', markersize=10, label='End')
    ax.set_xlabel('R (µm)')
    ax.set_ylabel('Rdot (m/s)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Phase Portrait')
    ax.legend()
    
    # Bubble mass
    ax = axes[3, 1]
    ax.plot(df.time_us, df.m_b * 1e9, 'cyan', linewidth=2)  # nanograms
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Bubble Mass (ng)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Bubble Vapor Mass')
    
    # Driving pressure
    ax = axes[3, 2]
    driving_pressure = df.Pb_bar - df.Pamb_bar
    ax.plot(df.time_us, driving_pressure, 'darkgreen', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Pb - P_amb (bar)')
    ax.grid(True, alpha=0.3)
    ax.set_title('Driving Pressure')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = csv_file.replace('.csv', '_validation.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    
    # Show plot
    # plt.show()
    
    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Initial R:       {df.R_um.iloc[0]:.6f} µm")
    print(f"Final R:         {df.R_um.iloc[-1]:.6f} µm")
    print(f"Initial Ro:      {df.Ro_um.iloc[0]:.6f} µm")
    print(f"Final Ro:        {df.Ro_um.iloc[-1]:.6f} µm")
    print(f"Max Rdot:        {df.Rdot.max():.3f} m/s")
    print(f"Initial T:       {df.T_drop.iloc[0]:.2f} K")
    print(f"Final T:         {df.T_drop.iloc[-1]:.2f} K")
    print(f"Temperature drop: {df.T_drop.iloc[0] - df.T_drop.iloc[-1]:.2f} K")
    print(f"Max Nu:          {df.Nu.max():.1f}")
    print(f"Min Nu:          {df.Nu.min():.1f}")
    print(f"Max Q_conv:      {df.Q_conv.max():.3e} W")
    print(f"Total time:      {df.time.iloc[-1]*1e6:.2f} µs")
    print("="*60)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        csv_file = 'rpe_parcel_debug.csv'
    
    plot_rpe_debug(csv_file)
