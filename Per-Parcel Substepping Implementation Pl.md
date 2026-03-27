Per-Parcel Substepping Implementation Plan for spray_evap.c (v3.1.12)

Scope

This document describes how to implement adaptive per-parcel substepping in the older spray_evap.c version (v3.1.12, the first file originally provided).

The goal is to advance stiff parcels with a smaller effective integration step inside one CFD timestep while leaving the global CFD timestep unchanged.

This plan is intentionally limited to the substepping procedure itself.

⸻

Objective

Add a parcel-local subcycling layer so that droplets which experience strong evaporation / strong thermal response can be integrated more accurately within the same CFD timestep.

The implementation should:
	•	keep the existing evaporation physics
	•	keep the existing implicit parcel solver logic
	•	subdivide the parcel update only when needed
	•	preserve the existing CFD timestep and gas-phase update structure

⸻

What Problem This Addresses

The current v3.1.12 routine advances each parcel once per CFD timestep.

For small or hot droplets, that can produce:
	•	large drdt
	•	large dm_dt
	•	large temperature drop within one CFD step
	•	poor convergence of the implicit parcel update
	•	unstable post-update parcel states

Per-parcel substepping reduces the size of the nonlinear update seen by the parcel solver.

⸻

High-Level Design

Old behavior

For each CFD timestep:
	1.	Read gas-cell fields
	2.	Loop over clouds
	3.	Loop over parcels
	4.	Run one implicit evaporation / temperature solve for each parcel
	5.	Accumulate source terms once

New behavior

For each CFD timestep:
	1.	Read gas-cell fields
	2.	Loop over clouds
	3.	Loop over parcels
	4.	Estimate parcel stiffness
	5.	Choose a parcel-specific number of substeps n_sub
	6.	Advance that parcel through n_sub substeps
	7.	Run the same implicit solver inside each substep
	8.	Accumulate the parcel contribution for the full CFD timestep

⸻

Core Rule

The implicit parcel solver stays in place.

Substepping does not replace the implicit solver.
It only changes the time interval over which that implicit solver is applied.

Each substep should still solve the same evaporation and temperature update logic, but over dt_sub = dt_cfd / n_sub.

⸻

Recommended Refactor Strategy

The current v3.1.12 code mixes several responsibilities inside one long parcel loop:
	•	property evaluation
	•	evaporation rate calculation
	•	implicit temperature iteration
	•	parcel state update
	•	source accumulation

For substepping, separate these responsibilities conceptually into:

1. Stiffness estimation

Determine whether a parcel needs substepping and how many substeps to use.

2. One parcel substep update

Advance a single parcel over dt_sub using the existing implicit evaporation logic.

3. Full-step accumulation

Sum parcel contributions over all substeps and apply them once for the CFD step.

This can be implemented either by extracting a helper function or by creating a repeated block inside the parcel loop.

⸻

Where to Modify the v3.1.12 File

In the first uploaded version, the main target is spray_evap_cell(CONVERGE_cloud_t cloud).

The insertion point is:
	•	after the parcel and gas fields are loaded
	•	after the transport-property setup
	•	before the current long implicit evaporation / temperature solve begins

In practice, the current single-pass parcel update should be wrapped in a substep loop.

⸻

Recommended Data Flow

For each parcel i_pc:

Step 1: Take a stiffness snapshot

Use the parcel state at the beginning of the CFD timestep to estimate the local thermal / evaporation response.

A practical estimate can be based on:
	•	current droplet radius
	•	current droplet temperature
	•	current estimated dm_dt or drdt
	•	latent heat
	•	liquid heat capacity

Step 2: Compute the number of substeps

A typical form is:

n_sub(i) = ceil(dt_cfd / (alpha * tau_local))

Where:
	•	tau_local is a local evaporation / thermal response timescale
	•	alpha is a safety factor less than 1
	•	n_sub(i) is clamped to a reasonable maximum

Step 3: March the parcel state in substeps

For j = 1..n_sub(i):
	•	set dt_sub = dt_cfd / n_sub(i)
	•	evaluate evaporation physics at the current substep state
	•	run the implicit solve
	•	update radius / temperature / mass fractions
	•	store local source increments

Step 4: Finish the parcel step

After the final substep:
	•	write the parcel’s final state back to the cloud arrays
	•	accumulate source terms for the full CFD timestep

⸻

Suggested Timescale for Substep Selection

The most useful trigger is a physics-based thermal/evaporation timescale.

One option is:

tau_evap ~ (m_drop * c_p_liquid) / (|dm_dt| * h_vap + eps)

Interpretation:
	•	if latent heat removal is fast relative to the droplet’s thermal mass, the parcel is stiff
	•	smaller tau_evap means more substeps are needed

If dm_dt is not yet reliable before the first estimate, a fallback geometric estimate can be used:

tau_R ~ R / (|drdt| + eps)

Then use the smaller of the two.

⸻

Substep Count Selection

Recommended structure:

n_sub = clamp(ceil(dt_cfd / (alpha * tau_local)), 1, n_sub_max)

Suggested starting values:
	•	alpha = 0.2 to 0.5
	•	n_sub_max = 5 initially
	•	n_sub_max = 10 only if needed later

Do not make substepping unconditional for every parcel unless absolutely necessary.

⸻

What Must Recompute Inside Each Substep

Inside every substep, recompute the same nonlinear quantities that the current solver already uses:
	•	droplet temperature-dependent properties
	•	vapor pressure
	•	latent heat
	•	mass transfer coefficient
	•	y1_star
	•	B_M
	•	drdt
	•	dm_dt
	•	parcel temperature update

Do not compute these once and reuse them for all substeps.

The point of substepping is to let the nonlinear state evolve between increments.

⸻

What Must Remain Frozen During Substepping

For the duration of one CFD timestep, the following should be treated as fixed inputs to the parcel substep loop:
	•	gas-cell density
	•	gas-cell pressure
	•	gas-cell temperature
	•	gas-cell species mass fractions
	•	gas-cell transport properties already retrieved for the CFD step

This is a parcel-local refinement only.
The gas phase is not advanced separately inside the parcel substeps.

⸻

State Variables to Preserve Carefully

The code currently uses both current-step and previous-step parcel values.

Keep as CFD-step history

Do not overwrite these during the substep loop:
	•	temp_tm1
	•	radius_tm1
	•	mfrac_tm1

Update during substeps

The following should evolve during substeps:
	•	temp
	•	radius
	•	mfrac
	•	drdt
	•	dm_dt
	•	temp_starm1 if the solver logic requires it as an iterate-local variable

The key rule is that _tm1 arrays remain tied to the CFD step, not the substep.

⸻

Implementation Shape

A good implementation shape is:

for (i_pc = first parcel; i_pc != -1; i_pc = next parcel)
{
    // 1. Estimate stiffness
    // 2. Compute n_sub
    // 3. dt_sub = dt / n_sub

    for (sub = 0; sub < n_sub; ++sub)
    {
        // 4. Recompute property / transfer terms
        // 5. Run the existing implicit parcel solve
        // 6. Update parcel state over dt_sub
        // 7. Accumulate substep-local source increments
    }

    // 8. Write final parcel state and parcel-level source contribution
}

The exact helper names are up to the implementation, but this structure should be preserved.

⸻

Helper Function Suggestion

A clean implementation is easier if one substep is isolated in a helper with input/output arguments.

Suggested helper purpose

Advance one parcel by one substep using the existing evaporation solver.

Suggested inputs
	•	parcel cloud handle
	•	parcel index
	•	frozen gas-cell state
	•	dt_sub
	•	species/table context
	•	flags needed by the solver

Suggested outputs
	•	updated parcel state
	•	local source increments for mass / energy / species
	•	convergence / failure status

The helper should not directly apply cell-wide sources.

⸻

Source Accumulation Strategy

The safest approach is to accumulate all substep increments locally and apply the full parcel contribution once.

That means tracking, per parcel:
	•	total evaporated mass over all substeps
	•	total latent / sensible energy contribution
	•	total species transfer over all substeps

At the end of the parcel’s substep loop, the final accumulated values are written into the existing source arrays.

This prevents double-counting and keeps the CFD-step bookkeeping consistent.

⸻

Interaction with the Existing Implicit Solver

The current implicit solver should remain the inner solver inside each substep.

That means each substep still performs the current iterative logic that solves for parcel temperature and evaporation state.

The only change is that the solver sees a smaller timestep and a smaller state jump.

This is the essential benefit of the refactor.

⸻

Debugging / Validation Requirements

Add diagnostics that confirm whether substepping is being activated correctly.

Recommended debug output for initial testing:
	•	parcel index
	•	cloud index
	•	chosen n_sub
	•	dt_sub
	•	initial and final parcel temperature over the CFD step
	•	initial and final radius
	•	number of implicit iterations per substep
	•	whether the solver recovered or failed

Suggested logic:
	•	print only for parcels with n_sub > 1
	•	print only for a small sample of parcels initially
	•	print the first few occurrences only

This keeps the log manageable.

⸻

Validation Tests

Test 1: Baseline equivalence

Set n_sub = 1 for all parcels.

Expected outcome:
	•	the new code should reproduce the current behavior as closely as possible
	•	any differences should be minimal and explainable

Test 2: Stiff small-droplet case

Run the smaller-diameter droplet case that previously produced many dT > 10 K events.

Expected outcome:
	•	n_sub should increase for the stiff parcels
	•	large single-step temperature jumps should reduce
	•	the implicit solver should converge more consistently

Test 3: Mixed cloud case

Run a case with both large and small droplets.

Expected outcome:
	•	only the stiff parcels should subcycle
	•	the rest should remain at n_sub = 1 or low values

Test 4: Conservation check

Compare total source terms against the non-subcycled baseline.

Expected outcome:
	•	total cell mass, energy, and species source accounting should remain consistent
	•	the substepping should change the path, not break bookkeeping

⸻

Performance Expectations

The additional cost should be localized.

Only parcels that exceed the stiffness threshold should subcycle.
The majority of parcels should remain at n_sub = 1 in normal cases.

This keeps the extra cost manageable.

⸻

Failure Modes to Watch For

1. Reusing one drdt across all substeps

This will defeat the purpose of substepping.

2. Updating _tm1 arrays inside the substep loop

This will corrupt the CFD-step history.

3. Double-applying cell source terms

Substep source increments must not be written repeatedly to the global arrays without accumulation logic.

4. Subcycling without recomputing nonlinear physics

This turns substepping into a cosmetic loop and does not improve stability.

5. Over-allocating substeps for every parcel

This can make the routine too expensive with little benefit.

⸻

Minimum Viable Implementation

The smallest useful version should do the following:
	1.	estimate a parcel stiffness metric
	2.	compute n_sub
	3.	set dt_sub = dt / n_sub
	4.	wrap the existing parcel implicit update in a substep loop
	5.	recompute the nonlinear evaporation quantities inside each substep
	6.	accumulate source increments once per parcel
	7.	keep the gas phase frozen over the CFD step

That is enough to test the concept.

⸻

Acceptance Criteria

The substepping implementation is acceptable if:
	•	parcels with large thermal response automatically subcycle
	•	the existing implicit solver still runs inside each substep
	•	n_sub = 1 reproduces the current behavior
	•	source bookkeeping remains consistent
	•	the stiff-droplet cases become less erratic

⸻

Notes for the CLI Tool

A code-generation or refactoring tool should treat this as a parcel-local integration layer.

The tool should preserve:
	•	existing field retrieval
	•	existing table loading / unloading
	•	existing parcel iteration structure
	•	existing gas coupling structure

The tool should add only:
	•	substep selection
	•	substep loop
	•	local substep state management
	•	local accumulation of substep increments

No global solver redesign is required for this task.