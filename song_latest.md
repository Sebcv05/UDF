# Song RPE Implementation Progress - Latest Status

**Last Updated:** 2025-12-04 12:53 UTC

---

## 📊 Current Status: Phase 1 Complete ✅

### Overall Progress: 33% Complete

- ✅ **Phase 1: Runtime Flag** - COMPLETE AND VALIDATED
- ⏳ **Phase 2: Song Solver** - Ready to Begin
- ⏳ **Phase 3: Testing** - Pending Phase 2

---

## 🎯 Phase 1: Runtime Flag Implementation (COMPLETE)

### What Was Done:

**Code Implementation:**
1. ✅ Added `use_song_rpe` global flag in 3 core files
2. ✅ Reads from `user_inputs.in` without recompilation
3. ✅ Defaults to 0 (thermal model - backwards compatible)
4. ✅ Model selection switch added to `spray_drop_distort_NH3.c`

**Testing:**
- ✅ Test 1: `use_song_rpe = 0` → Thermal model (verified in user_inputs.echo)
- ✅ Test 2: `use_song_rpe = 1` → Song model flag (verified in user_inputs.echo)
- ✅ No recompilation needed to switch between models

**Files Modified:**
```
include/globals.h              - Declaration (extern use_song_rpe)
src/Breakup.c                  - Definition (= 0 default)
src/read_input.c               - Reading, logging, echo logic
src/RPE_euler.c                - Thermal model logging
src/spray_drop_distort_NH3.c   - Model selection if/else switch
```

**Git Status:**
- Branch: `v3.1.12`
- Latest Commit: `3cb8dfa`
- Pushed: Yes

### Reference Documents:

📄 **[PHASE1_COMPLETE.txt](PHASE1_COMPLETE.txt)** - Phase 1 completion summary

📄 **[PHASE1_TEST_RESULTS.txt](PHASE1_TEST_RESULTS.txt)** - Validation test results

📄 **[ADDING_USE_SONG_RPE_FLAG.md](ADDING_USE_SONG_RPE_FLAG.md)** - Step-by-step flag implementation guide

---

## 🔬 Standalone Validation (COMPLETE)

### Test Program Status: ✅ Validated

**Files:**
- `test_song_rpe.c` - Standalone C implementation (461 lines)
- `plot_song_c_results.py` - Plotting script
- `Makefile_song` - Build system

**Test Results:** (Temperature sweep at P=2 bar)
```
Condition    Time(μs)  R_final(μm)  ε_final  V_max(m/s)  Steps
────────────────────────────────────────────────────────────────
T=255K       5.87      49.89        0.989    20.35       589
T=263K       4.56      49.81        0.977    38.00       471
T=273K       2.79      49.75        0.967    94.53       299
T=283K       2.17      49.73        0.961    53.18       261
T=293K       1.72      49.62        0.949    47.82       238
T=303K       1.42      49.86        0.957    59.43       223
T=313K       1.20      49.80        0.948    68.40       201
T=323K       1.03      49.69        0.936    77.10       173
```

**Validation:** ✅ All conditions run successfully with smooth, physically realistic curves

**Key Parameters:**
- Initial timestep: 1e-10 s (10x reduction for smoothness)
- Safety factor: 0.05 (conservative)
- Void fraction target: 0.99
- Initial bubble radius: 1.01 × R_c (critical radius)

### Reference Documents:

📄 **[README_SONG_TEST.md](README_SONG_TEST.md)** - Test program usage guide

📄 **[SONG_IMPLEMENTATION_STATUS.md](SONG_IMPLEMENTATION_STATUS.md)** - Validation results and comparison

---

## 📋 Phase 2: Song Solver Implementation (NEXT)

### To Be Created:

**1. Header File: `include/RPE_song.h`**
- Structures: `SongParams`
- Function declarations
- Status: Not started

**2. Source File: `src/RPE_song.c`**
- Main solver: `RPE_song_solver()`
- Helper functions (adapted from `test_song_rpe.c`)
- Status: Not started

**3. Integration: `src/spray_drop_distort_NH3.c`**
- Replace stub with actual `RPE_song_solver()` call
- Status: Stub in place, ready for implementation

### Implementation Plan:

**Key Functions to Adapt:**
```c
// From test_song_rpe.c → RPE_song.c

compute_void_fraction()       - Keep as-is
compute_mixture_density()     - Keep as-is  
compute_song_acceleration()   - Adapt parameters
RPE_song_solver()            - Main changes here
```

**Data Flow:**
```
INPUT (from parcel_cloud[p_idx]):
  - r_bubble, v_bubble        (current state)
  - r_bubble_0, r_drop_0      (initial conditions)
  - temp, density, viscosity  (properties)

OUTPUT (to parcel_cloud[p_idx]):
  - r_bubble (updated)
  - v_bubble (updated)
  - temp (UNCHANGED - isothermal)
```

### Reference Documents:

📄 **[SONG_INTEGRATION_PLAN.md](SONG_INTEGRATION_PLAN.md)** - Complete integration guide (14KB, detailed)

📄 **[SONG_INTEGRATION_SUMMARY.txt](SONG_INTEGRATION_SUMMARY.txt)** - Quick reference checklist

📄 **[test_song_rpe.c](test_song_rpe.c)** - Template code (working standalone implementation)

---

## 🔑 Key Differences: Thermal vs Song Models

| Feature | Thermal RPE (Current) | Song RPE (New) |
|---------|----------------------|----------------|
| **ODEs** | 4 (R, Rdot, T, m_b) | 2 (R, Rdot) |
| **Temperature** | Evolves (dT/dt ≠ 0) | Constant (isothermal) |
| **Mass transfer** | Thermal limited (mdot=Q/L_v) | Not tracked |
| **Mixture density** | Constant ρ_l | Variable ρ_m(ε) |
| **Initial pressure** | Not included | (2σ/R₀+P_r0)·(R₀/R)³ |
| **Bubble mass** | Tracked | Not needed |
| **Stop criteria** | Multiple checks | void = 0.99 |
| **Physics** | Heat-limited growth | Pressure-driven growth |

---

## 🧮 Song Model Physics

### Core Equation:
```
R̈ = [P_sat - P_∞ + (2σ/R₀+P_r0)·(R₀/R)³ - 2σ/R - 4μ·Ṙ/R - 4κ·Ṙ/R²] / (ρ_m·R) - (3/2)·Ṙ²/R
```

### Key Parameters:
- **P_r0** = 1.0e6 Pa (residual gas pressure)
- **κ** = 0.0 (surface viscosity)
- **R_spec** = 488.2 J/(kg·K) for NH₃
- **Void fraction** = ε = R³/R_droplet_0³
- **Mixture density** = ρ_m = ε·ρ_v + (1-ε)·ρ_l

### Termination:
- Stop when void fraction ≥ 0.99
- Convert to child parcel (like thermal model)
- Temperature unchanged throughout

### Reference Documents:

📄 **[SONG_implement.md](../v3.1.12/Splitter/Dev/SONG_implement.md)** - Complete physics documentation (in Dev directory)

📄 **Physics Section in test_song_rpe.c** - Lines 140-172: Acceleration computation

---

## 📁 File Organization

### UDF Repository (`/home/apollo19/Desktop/Dan_B/UDF/`)

```
├── src/
│   ├── RPE_euler.c              ✅ Current thermal model
│   ├── RPE_song.c               ⏳ TO CREATE (Phase 2)
│   ├── spray_drop_distort_NH3.c ✅ Switch in place, stub ready
│   ├── read_input.c             ✅ Flag reading implemented
│   └── Breakup.c                ✅ Flag defined
│
├── include/
│   ├── RPE_euler.h              ✅ Current thermal header
│   ├── RPE_song.h               ⏳ TO CREATE (Phase 2)
│   └── globals.h                ✅ Flag declared
│
├── Test & Validation/
│   ├── test_song_rpe.c          ✅ Standalone template
│   ├── plot_song_c_results.py   ✅ Plotting tool
│   ├── Makefile_song            ✅ Build system
│   ├── song_temp_sweep_c.csv    ✅ Test results
│   └── song_temp_sweep_c*.pdf   ✅ Result plots
│
└── Documentation/
    ├── README_SONG_TEST.md              ✅ Test usage
    ├── SONG_IMPLEMENTATION_STATUS.md    ✅ Overall status
    ├── SONG_INTEGRATION_PLAN.md         ✅ Integration guide
    ├── SONG_INTEGRATION_SUMMARY.txt     ✅ Quick reference
    ├── ADDING_USE_SONG_RPE_FLAG.md      ✅ Flag guide
    ├── PHASE1_COMPLETE.txt              ✅ Phase 1 summary
    ├── PHASE1_TEST_RESULTS.txt          ✅ Test validation
    └── song_latest.md                   📍 THIS FILE
```

### Case Directory (`/home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev/`)

```
├── user_inputs.in       ✅ Flag added: use_song_rpe = 0
├── user_inputs.echo     ✅ Auto-generated verification
└── outputs_original/
    └── converge.log     ✅ Shows flag value on startup
```

---

## 🚀 Next Steps

### Immediate (Phase 2):

1. **Create `include/RPE_song.h`** (~10 minutes)
   - Copy structure from test_song_rpe.c
   - Adapt types to CONVERGE conventions
   - Add function declarations

2. **Create `src/RPE_song.c`** (~30 minutes)
   - Copy helper functions from test_song_rpe.c
   - Adapt main solver function
   - Read from parcel_cloud structure
   - Single timestep integration (not time loop)

3. **Update `src/spray_drop_distort_NH3.c`** (~5 minutes)
   - Replace stub printf with actual RPE_song_solver() call
   - Add include for RPE_song.h

4. **Compile and test** (~10 minutes)
   - Check for compilation errors
   - Run with use_song_rpe=0 (thermal) - should work as before
   - Run with use_song_rpe=1 (Song) - should call new solver

### Testing (Phase 3):

1. Single parcel test
2. Thermal vs Song comparison
3. Mass conservation check
4. Full spray validation

**Estimated Total Time:** ~1.5 hours for complete Phase 2 + 3

---

## 📚 Complete Documentation Index

### Implementation Guides:
1. **[SONG_INTEGRATION_PLAN.md](SONG_INTEGRATION_PLAN.md)** - START HERE for Phase 2
2. **[ADDING_USE_SONG_RPE_FLAG.md](ADDING_USE_SONG_RPE_FLAG.md)** - How flag system works
3. **[SONG_INTEGRATION_SUMMARY.txt](SONG_INTEGRATION_SUMMARY.txt)** - Quick checklist

### Physics & Theory:
4. **[SONG_implement.md](../v3.1.12/Splitter/Dev/SONG_implement.md)** - Complete physics equations
5. **[SONG_IMPLEMENTATION_STATUS.md](SONG_IMPLEMENTATION_STATUS.md)** - Thermal vs Song comparison

### Testing & Validation:
6. **[README_SONG_TEST.md](README_SONG_TEST.md)** - How to run standalone test
7. **[PHASE1_TEST_RESULTS.txt](PHASE1_TEST_RESULTS.txt)** - Phase 1 validation
8. **[test_song_rpe.c](test_song_rpe.c)** - Working template code

### Status Documents:
9. **[PHASE1_COMPLETE.txt](PHASE1_COMPLETE.txt)** - Phase 1 summary
10. **[song_latest.md](song_latest.md)** - THIS FILE - overall progress

---

## 💡 Key Success Factors

### ✅ Completed:
- Standalone Song solver validated against Python reference
- Runtime flag system working (no recompilation needed)
- All required parcel variables exist in framework
- Clean code structure following existing patterns
- Comprehensive documentation at multiple detail levels

### 🎯 For Phase 2:
- Use test_song_rpe.c as template (already working)
- Main change: time loop → single timestep
- Data source: test conditions → parcel_cloud
- Keep physics calculations identical
- Maintain same function structure

### ⚠️ Important Notes:
- Song model is **isothermal** - temperature does NOT change
- Void fraction calculated from **initial** droplet radius (R_drop_0)
- Initial bubble radius (R_bubble_0) needed for pressure term
- Stop at void = 0.99, then convert to child parcel
- Model selection at runtime via user_inputs.in

---

## 🔄 Git Repository Status

**Branch:** `v3.1.12`
**Latest Commit:** `3cb8dfa` - "Phase 1 Testing: Add model selection logging and validate flag"

**Recent Commits:**
```
3cb8dfa  Phase 1 Testing: Add model selection logging and validate flag
8bd8a28  Add Phase 1 completion summary  
c685ff0  Phase 1: Add use_song_rpe runtime flag
c370f71  Add detailed guide for implementing use_song_rpe runtime flag
81efbb0  Add detailed Song RPE integration plan and quick reference
aa1d213  Reduce timestep by 10x for smoother integration
949ea76  Add Song RPE standalone test program
```

**Backup:** Branch `v1.0.1` created before Song implementation started

---

## 📞 Quick Reference

### To Switch Models (No Recompilation):

```bash
# Edit user_inputs.in
cd /home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev

# For thermal model:
echo "0     use_song_rpe" >> user_inputs.in

# For Song model:
echo "1     use_song_rpe" >> user_inputs.in

# Run
./run.sh

# Verify which model was used
cat user_inputs.echo | grep use_song_rpe
```

### To Compile:

```bash
cd /home/apollo19/Desktop/Dan_B/v3.1.12/Splitter/Dev
./upc2.sh
```

### To Run Standalone Test:

```bash
cd /home/apollo19/Desktop/Dan_B/UDF
gcc test_song_rpe.c -o test_song_rpe -O2 -Wall -std=c99 -lm
./test_song_rpe
python3 plot_song_c_results.py
```

---

## 🎓 Learning Resources

**Understanding the Code:**
1. Read `test_song_rpe.c` - Complete working implementation
2. Compare with `RPE_euler.c` - See differences from thermal model
3. Review `SONG_INTEGRATION_PLAN.md` - Detailed mapping

**Physics Background:**
1. `SONG_implement.md` - Equation derivations
2. `test_song_rpe.c` lines 140-172 - Acceleration computation
3. Plots in `song_temp_sweep_c.pdf` - Expected behavior

---

## ✅ Validation Checklist

### Phase 1 (Complete):
- [x] Flag declared in globals.h
- [x] Flag defined in Breakup.c
- [x] Flag read in read_input.c
- [x] Flag echoed to user_inputs.echo
- [x] Model selection switch in spray_drop_distort_NH3.c
- [x] Test with flag=0 (thermal)
- [x] Test with flag=1 (Song stub)
- [x] Verify no recompilation needed

### Phase 2 (Pending):
- [ ] Create RPE_song.h
- [ ] Create RPE_song.c with helper functions
- [ ] Implement RPE_song_solver()
- [ ] Replace stub in spray_drop_distort_NH3.c
- [ ] Compile without errors
- [ ] Link without errors

### Phase 3 (Pending):
- [ ] Test single parcel
- [ ] Compare thermal vs Song
- [ ] Check mass conservation
- [ ] Full spray test
- [ ] Validate against standalone test
- [ ] Performance check

---

## 🏁 Summary

**Status:** Phase 1 Complete and Validated ✅

**Progress:** 33% (1 of 3 phases)

**Ready:** Phase 2 can begin immediately

**Confidence:** High - Standalone solver validated, flag system tested, all documentation in place

**Timeline:** ~1.5 hours to complete Phase 2 + 3

---

**For Questions:** Refer to detailed guides in documentation index above

**To Continue:** Start with SONG_INTEGRATION_PLAN.md, Section "Step 2: Create Source File"

**Last Updated:** 2025-12-04 12:53 UTC
