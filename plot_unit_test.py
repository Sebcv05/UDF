#!/usr/bin/env python3
"""
Visualize RPE unit test results
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_rpe_unit_test(csv_file='rpe_unit_test.csv'):
    """Plot RPE unit test results"""
    
    # Load data
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} data points")
    print(f"Time range: {df.time_us.min():.3f} to {df.time_us.max():.3f} µs")
    
    # Create plots
    fig, axes = plt.subplots(3, 3, figsize=(15, 10))
    fig.suptitle('RPE_euler Unit Test: NH3 Bubble Growth', fontsize=14, fontweight='bold')
    
    # Row 1: Basic dynamics
    ax = axes[0, 0]
    ax.plot(df.time_us, df.R_um, 'b-', linewidth=2)
    ax.axhline(82.5, color='r', linestyle='--', label='Ro (droplet)')
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('R (µm)')
    ax.set_title('Bubble Radius')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    ax = axes[0, 1]
    ax.plot(df.time_us, df.Rdot_m_s, 'g-', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Rdot (m/s)')
    ax.set_title('Bubble Velocity')
    ax.grid(True, alpha=0.3)
    
    ax = axes[0, 2]
    ax.plot(df.time_us, df.dRdotdt, 'm-', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Acceleration (m/s²)')
    ax.set_title('Bubble Acceleration')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    
    # Row 2: Thermal
    ax = axes[1, 0]
    ax.plot(df.time_us, df.Pb_bar, 'b-', linewidth=2, label='Pb')
    ax.plot(df.time_us, df.P_sat_bar, 'r--', linewidth=2, label='P_sat')
    ax.axhline(2.0, color='k', linestyle=':', linewidth=2, label='P_amb')
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Pressure (bar)')
    ax.set_title('Pressure Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 1]
    ax.plot(df.time_us, df.T_drop_K, 'r-', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('Droplet Temperature')
    ax.grid(True, alpha=0.3)
    T_init = df.T_drop_K.iloc[0]
    T_final = df.T_drop_K.iloc[-1]
    ax.text(0.05, 0.95, f'ΔT = {T_init-T_final:.2f} K', 
            transform=ax.transAxes, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax = axes[1, 2]
    ax.plot(df.time_us, df.Nu, 'purple', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Nusselt Number')
    ax.set_title('Nu (Heat Transfer)')
    ax.grid(True, alpha=0.3)
    
    # Row 3: Diagnostics
    ax = axes[2, 0]
    ax.plot(df.R_um, df.Rdot_m_s, 'b-', linewidth=2)
    ax.plot(df.R_um.iloc[0], df.Rdot_m_s.iloc[0], 'go', markersize=10, label='Start')
    ax.plot(df.R_um.iloc[-1], df.Rdot_m_s.iloc[-1], 'ro', markersize=10, label='End')
    ax.set_xlabel('R (µm)')
    ax.set_ylabel('Rdot (m/s)')
    ax.set_title('Phase Portrait')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    ax = axes[2, 1]
    ax.plot(df.time_us, df.m_b_ng, 'cyan', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Bubble Mass (ng)')
    ax.set_title('Vapor Mass')
    ax.grid(True, alpha=0.3)
    
    ax = axes[2, 2]
    driving_pressure = df.Pb_bar - 2.0  # P_amb = 2 bar
    ax.plot(df.time_us, driving_pressure, 'darkgreen', linewidth=2)
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Pb - P_amb (bar)')
    ax.set_title('Driving Pressure')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('rpe_unit_test.png', dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: rpe_unit_test.png")
    plt.show()
    
    # Summary statistics
    print("\n" + "="*70)
    print("UNIT TEST SUMMARY")
    print("="*70)
    print(f"Initial R:       {df.R_um.iloc[0]:.6f} µm")
    print(f"Final R:         {df.R_um.iloc[-1]:.3f} µm")
    print(f"Growth:          {df.R_um.iloc[-1] - df.R_um.iloc[0]:.3f} µm")
    print(f"Max Rdot:        {df.Rdot_m_s.max():.3f} m/s")
    print(f"Initial T:       {df.T_drop_K.iloc[0]:.2f} K")
    print(f"Final T:         {df.T_drop_K.iloc[-1]:.2f} K")
    print(f"Temperature drop: {df.T_drop_K.iloc[0] - df.T_drop_K.iloc[-1]:.2f} K")
    print(f"Max Nu:          {df.Nu.max():.1f}")
    print(f"Total time:      {df.time_us.iloc[-1]:.2f} µs")
    print(f"Final R/Ro:      {100*df.R_um.iloc[-1]/82.5:.1f}%")
    print("="*70)
    
    # Check success criteria
    print("\nSUCCESS CRITERIA:")
    R_growth = df.R_um.iloc[-1] / df.R_um.iloc[0]
    T_drop = df.T_drop_K.iloc[0] - df.T_drop_K.iloc[-1]
    
    print(f"  ✓ R grew by >100x:          {R_growth:.1f}x {'PASS' if R_growth > 100 else 'FAIL'}")
    print(f"  ✓ Max Rdot > 10 m/s:        {df.Rdot_m_s.max():.1f} m/s {'PASS' if df.Rdot_m_s.max() > 10 else 'FAIL'}")
    print(f"  ✓ Temperature drop 1-10 K:  {T_drop:.2f} K {'PASS' if 1 < T_drop < 10 else 'FAIL'}")
    print(f"  ✓ Nu in range [2, 1000]:    {df.Nu.min():.1f}-{df.Nu.max():.1f} {'PASS' if df.Nu.min() >= 2 and df.Nu.max() < 1000 else 'FAIL'}")
    print(f"  ✓ No NaN values:            {'PASS' if not df.isnull().any().any() else 'FAIL'}")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        csv_file = 'rpe_unit_test.csv'
    
    plot_rpe_unit_test(csv_file)
