// ============================================================================
// FALLBACK: Uniform distribution (used if RR sampling fails)
// ============================================================================

void fallback_uniform_children(
    CONVERGE_precision_t parent_radius,
    CONVERGE_precision_t parent_num_drop,
    CONVERGE_precision_t R32_target,
    int N,
    CONVERGE_precision_t* child_radii,
    CONVERGE_precision_t* child_num_drop
) {
    // All children get same radius (R32_target) and same num_drop
    CONVERGE_precision_t parent_volume = CONVERGE_cube(parent_radius);
    CONVERGE_precision_t child_volume = CONVERGE_cube(R32_target);
    CONVERGE_precision_t uniform_num_drop = parent_num_drop * parent_volume / (N * child_volume);
    
    for (int i = 0; i < N; i++) {
        child_radii[i] = R32_target;
        child_num_drop[i] = uniform_num_drop;
    }
    
    static int fallback_warn_count = 0;
    if (fallback_warn_count < 3) {
        printf("[RR_FALLBACK] Using uniform distribution: all R=%.3e m, num_drop=%.3e\n",
               R32_target, uniform_num_drop);
        fallback_warn_count++;
    }
}

// ============================================================================
// ROSIN-RAMMLER CHILD SAMPLING FUNCTION
// 
// Purpose: Sample N child radii and num_drop values from RR distribution
//          while conserving D32 (Sauter Mean Diameter) and total volume
//
// Inputs:
//   - parent_radius: Parent droplet radius (m)
//   - parent_num_drop: Parent num_drop value
//   - R32_target: Target child radius (= calculated_radius from breakup model)
//   - N: Number of children (typically 12)
//   - n_RR: Rosin-Rammler shape parameter (2.5-4.0, default 3.2)
//   - gamma_ratio: Pre-computed tgamma(1+2/n) / tgamma(1+3/n)
//
// Outputs (arrays of size N):
//   - child_radii[N]: Sampled child radii (m)
//   - child_num_drop[N]: Corresponding num_drop values
//
// Conservation:
//   - D32 = sum(D^3) / sum(D^2) = 2 * R32_target (enforced exactly)
//   - Total volume = sum(num_drop[i] * R[i]^3) = parent_num_drop * parent_radius^3
//   - NO DENSITY used (pure volume conservation)
//
// Returns: 0 on success, -1 on failure (fallback to uniform distribution)
//
// ============================================================================

int sample_RR_children(
    CONVERGE_precision_t parent_radius,
    CONVERGE_precision_t parent_num_drop,
    CONVERGE_precision_t R32_target,
    int N,
    double n_RR,
    double gamma_ratio,
    CONVERGE_precision_t* child_radii,
    CONVERGE_precision_t* child_num_drop
) {
    const double PI = 3.14159265358979323846;
    const double U_MIN = 1.0e-12;  // Guard against log(0) or tiny u
    const double U_MAX = 1.0 - 1.0e-12;  // Guard against log(0)
    const double SCALE_MIN = 0.1;  // Minimum scale correction factor
    const double SCALE_MAX = 10.0;  // Maximum scale correction factor
    
    static int error_count = 0;
    const int MAX_ERRORS = 5;  // Limit error messages
    
    // Input validation
    if (parent_radius <= 0.0 || parent_num_drop <= 0.0 || R32_target <= 0.0 || N <= 0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Invalid inputs: R_p=%.3e, Nd_p=%.3e, R32=%.3e, N=%d\n",
                   parent_radius, parent_num_drop, R32_target, N);
            error_count++;
        }
        return -1;  // Signal failure
    }
    
    if (n_RR <= 0.0 || gamma_ratio <= 0.0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Invalid RR parameters: n_RR=%.3f, gamma_ratio=%.6f\n",
                   n_RR, gamma_ratio);
            error_count++;
        }
        return -1;
    }
    
    // Target D32 (diameter) from R32 (radius)
    CONVERGE_precision_t D32_target = 2.0 * R32_target;
    
    // Calculate RR scale parameter X using pre-computed gamma_ratio
    CONVERGE_precision_t X_RR = D32_target * gamma_ratio;
    
    // Precompute inverse of n_RR for efficiency
    double inv_n_RR = 1.0 / n_RR;
    
    // Step 1: Sample N diameters from RR distribution
    CONVERGE_precision_t child_diameters[N];
    CONVERGE_precision_t total_sampled_volume = 0.0;
    
    for (int i = 0; i < N; i++) {
        CONVERGE_precision_t u = CONVERGE_random_precision();  // Uniform [0,1]
        
        // Guard against extreme RNG values
        if (u < U_MIN) u = U_MIN;
        if (u > U_MAX) u = U_MAX;
        
        // Inverse CDF: D = X * (-ln(1-u))^(1/n)
        double log_arg = 1.0 - u;
        if (log_arg <= 0.0) {
            // Pathological case - should never happen with guards above
            if (error_count < MAX_ERRORS) {
                printf("[RR_ERROR] Invalid log argument: u=%.12e, log_arg=%.12e\n", u, log_arg);
                error_count++;
            }
            return -1;
        }
        
        child_diameters[i] = X_RR * pow(-log(log_arg), inv_n_RR);
        child_radii[i] = child_diameters[i] / 2.0;
        
        // Accumulate volume (normalized: no 4/3*pi factor)
        CONVERGE_precision_t V_i = CONVERGE_cube(child_radii[i]);
        total_sampled_volume += V_i;
    }
    
    // Guard against pathological sampling
    if (total_sampled_volume <= 0.0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Zero or negative sampled volume: %.3e. Aborting RR sampling.\n",
                   total_sampled_volume);
            error_count++;
        }
        return -1;
    }
    
    // Step 2: Calculate actual sampled D32 before correction
    CONVERGE_precision_t D2_sum = 0.0;
    CONVERGE_precision_t D3_sum = 0.0;
    for (int i = 0; i < N; i++) {
        CONVERGE_precision_t D = child_diameters[i];
        D2_sum += D * D;
        D3_sum += D * D * D;
    }
    
    // Guard against division by zero
    if (D2_sum <= 0.0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Zero D2_sum. Cannot compute D32. Aborting.\n");
            error_count++;
        }
        return -1;
    }
    
    CONVERGE_precision_t D32_sample = D3_sum / D2_sum;
    
    // Guard against pathological D32_sample
    if (D32_sample <= 0.0 || !isfinite(D32_sample)) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Invalid D32_sample: %.3e. Aborting.\n", D32_sample);
            error_count++;
        }
        return -1;
    }
    
    // Step 3: Apply per-parent correction to enforce exact D32
    CONVERGE_precision_t scale_correction = D32_target / D32_sample;
    
    // Guard against extreme scale corrections
    if (scale_correction < SCALE_MIN || scale_correction > SCALE_MAX) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_WARNING] Extreme scale correction: %.4f (clamping to [%.1f, %.1f]). D32_target=%.3e, D32_sample=%.3e\n",
                   scale_correction, SCALE_MIN, SCALE_MAX, D32_target, D32_sample);
            error_count++;
        }
        // Clamp to safe range
        if (scale_correction < SCALE_MIN) scale_correction = SCALE_MIN;
        if (scale_correction > SCALE_MAX) scale_correction = SCALE_MAX;
    }
    
    total_sampled_volume = 0.0;  // Recalculate after scaling
    for (int i = 0; i < N; i++) {
        child_diameters[i] *= scale_correction;
        child_radii[i] = child_diameters[i] / 2.0;
        total_sampled_volume += CONVERGE_cube(child_radii[i]);
    }
    
    // Guard against zero volume after scaling (shouldn't happen with proper clamping)
    if (total_sampled_volume <= 0.0) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Zero sampled volume after scaling. Aborting.\n");
            error_count++;
        }
        return -1;
    }
    
    // Step 4: Calculate uniform num_drop to conserve total volume
    // IMPORTANT: Volume convention is num_drop × R³ (NO 4/3*pi factor)
    // Parent volume (normalized): V_p = parent_num_drop * R_p^3
    // Child volume (normalized):  V_c = sum(Nd_i * R_i^3) = Nd * sum(R_i^3)
    // Conservation: V_c = V_p
    //   => Nd * sum(R_i^3) = parent_num_drop * R_p^3
    //   => Nd = parent_num_drop * R_p^3 / sum(R_i^3)
    
    CONVERGE_precision_t parent_volume = CONVERGE_cube(parent_radius);
    CONVERGE_precision_t base_num_drop = parent_num_drop * parent_volume / total_sampled_volume;
    
    // Guard against invalid num_drop
    if (base_num_drop <= 0.0 || !isfinite(base_num_drop)) {
        if (error_count < MAX_ERRORS) {
            printf("[RR_ERROR] Invalid base_num_drop: %.3e. Check inputs.\n", base_num_drop);
            error_count++;
        }
        return -1;
    }
    
    // Assign same num_drop to all children
    for (int i = 0; i < N; i++) {
        child_num_drop[i] = base_num_drop;
    }
    
    // Diagnostic output (rate-limited)
    static int rr_diag_count = 0;
    if (rr_diag_count < 3) {
        printf("[RR_SAMPLE] Parent: R=%.3e m, num_drop=%.3e\n", parent_radius, parent_num_drop);
        printf("[RR_SAMPLE] Target: R32=%.3e m (D32=%.3e m)\n", R32_target, D32_target);
        printf("[RR_SAMPLE] X_RR=%.3e m, n_RR=%.2f, gamma_ratio=%.6f\n", 
               X_RR, n_RR, gamma_ratio);
        printf("[RR_SAMPLE] Before correction: D32_sample=%.3e m (ratio=%.4f)\n", 
               D32_sample, D32_sample/D32_target);
        printf("[RR_SAMPLE] After correction: scale=%.4f (clamped to [%.1f, %.1f])\n", 
               scale_correction, SCALE_MIN, SCALE_MAX);
        printf("[RR_SAMPLE] Child radii (m): ");
        for (int i = 0; i < N; i++) {
            printf("%.3e ", child_radii[i]);
        }
        printf("\n");
        printf("[RR_SAMPLE] base_num_drop=%.3e (same for all children)\n", base_num_drop);
        
        // Verify conservation (IMPORTANT: No 4/3*pi factor)
        CONVERGE_precision_t total_child_volume = 0.0;
        for (int i = 0; i < N; i++) {
            total_child_volume += child_num_drop[i] * CONVERGE_cube(child_radii[i]);
        }
        CONVERGE_precision_t parent_total_volume = parent_num_drop * parent_volume;
        CONVERGE_precision_t volume_error = fabs(total_child_volume - parent_total_volume) / parent_total_volume;
        printf("[RR_SAMPLE] Volume check (num_drop*R^3): parent=%.3e, children=%.3e, error=%.2e%%\n",
               parent_total_volume, total_child_volume, volume_error * 100.0);
        printf("[RR_SAMPLE] NOTE: Volume = num_drop × R³ (NO 4/3*pi factor)\n");
        
        if (volume_error > 1.0e-6) {
            printf("[RR_WARNING] Volume error > 1 ppm. Check calculation.\n");
        }
        
        rr_diag_count++;
    }
    
    return 0;  // Success
}
