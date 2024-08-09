/*******************************************************************************
* CONVERGENT SCIENCE CONFIDENTIAL                                              *
* All rights reserved.                                                         *
* All information contained herein is the property of Convergent Science.      *
* The intellectual and technical concepts contained herein are                 *
* proprietary to Convergent Science.                                           *
* Dissemination of this information or reproduction of this material           *
* is strictly forbidden unless prior written permission is obtained from       *
* Convergent Science.                                                          *
*******************************************************************************/

#include "CONVERGE/udf_register.h"

#ifndef CONVERGE_SAMPLES_H_
#define CONVERGE_SAMPLES_H_

// General > Position
void before_main_loop();
void after_main_loop();
void before_time_step();
void after_time_step();
void before_chemistry();
void after_chemistry();
void before_dpm_source();
void after_dpm_source();
void before_transport();
void after_transport();
void before_turbulence();
void after_turbulence();
void before_property_update();
void after_property_update();
void before_grid_change();
void after_grid_change();
void before_dpm_move();
void after_dpm_move();
void before_write_file();
void after_write_file();
void before_write_post();
void after_write_post();
void before_write_restart();
void after_write_restart();
void before_write_transfer();
void after_write_transfer();
void before_source_main();
void after_source_main();
// void ga_merit();
void calc_statistical_variables_after();
void calc_statistical_variables_before();
// void post_write_points();

// Motion
// void piston_position(CONVERGE_PistonPositionInput_t* input, CONVERGE_PistonPositionOutput_t* output);
// void gradient_regression_setup(CONVERGE_MotionSetupInput_t* input);
// void gradient_regression_close(CONVERGE_MotionCloseInput_t* input);
// void gradient_regression_init(CONVERGE_MotionInitInput_t* input);
// void gradient_regression(CONVERGE_MotionInput_t* input, CONVERGE_MotionOutput_t* output);

// General > Other
void cyclic_event(CONVERGE_EventInput_t* input);
void cyclic_event_setup(CONVERGE_EventInput_t* input);
// void o2_init();
// void regions_flow();
// void mass_leakage(CONVERGE_SourceInput_t* input, CONVERGE_SourceOutput_t* output);
// void wiebe_law_heat_release_energy(CONVERGE_SourceInput_t* input, CONVERGE_SourceOutput_t* output);
// void wiebe_law_heat_release_energy_setup(CONVERGE_SourceInput_t* input);
// void per_unit_volume_time_species_source(CONVERGE_SourceInput_t* input, CONVERGE_SourceOutput_t* output);

// I/O
// void my_restart();
// void read_aniso_cond_input(CONVERGE_InputInput_t* input);
// void read_user(CONVERGE_InputInput_t* input);
// void read_source_input(CONVERGE_InputInput_t* input);
// void coupled_simulation_read_inputs(CONVERGE_InputInput_t* input);
// void soot_outfile_setup(CONVERGE_OutputSetupOutput_t* output);
// void soot_outfile(CONVERGE_OutputOutput_t* output);
// void soot_outfile_close(CONVERGE_OutputCloseOutput_t* output);
// void equiv_ratio_bin_setup(CONVERGE_OutputSetupOutput_t* output);
// void equiv_ratio_bin(CONVERGE_OutputOutput_t* output);
// void equiv_ratio_bin_close(CONVERGE_OutputCloseOutput_t* output);
// void mixing_output_setup(CONVERGE_OutputSetupOutput_t* output);
// void mixing_output(CONVERGE_OutputOutput_t* output);
// void mixing_output_close(CONVERGE_OutputCloseOutput_t* output);
// void react_ratio_bin_setup(CONVERGE_OutputSetupOutput_t* output);
// void react_ratio_bin(CONVERGE_OutputOutput_t* output);
// void react_ratio_bin_close(CONVERGE_OutputCloseOutput_t* output);
// void region_flow_output_setup(CONVERGE_OutputSetupOutput_t* output);
// void region_flow_output(CONVERGE_OutputOutput_t* output);
// void region_flow_output_close(CONVERGE_OutputCloseOutput_t* output);
// void outfile_method_setup(CONVERGE_OutputSetupOutput_t* output);
// void outfile_method(CONVERGE_OutputOutput_t* output);
// void outfile_method_close(CONVERGE_OutputCloseOutput_t* output);
// void wiebe_law_outfile_setup(CONVERGE_OutputSetupOutput_t* output);
// void wiebe_law_outfile(CONVERGE_OutputOutput_t* output);
// void wiebe_law_outfile_close(CONVERGE_OutputCloseOutput_t* output);
// void points_setup(CONVERGE_PointsSetupInput_t* input, CONVERGE_PointsOutput_t* output);
// void points(CONVERGE_PointsInput_t* input, CONVERGE_PointsOutput_t* output);
// void points_close(CONVERGE_PointsOutput_t* output);

// Properties
// void user_aniso_cond(CONVERGE_AnisoCondInput_t* input);
void gas_cp_cell(CONVERGE_CellPropertyInput_t* input, CONVERGE_CellPropertyOutput_t* output);
void gas_custom_prop(CONVERGE_CellPropertyInput_t* input, CONVERGE_CellPropertyOutput_t* output);
void gas_custom_nd_prop(CONVERGE_CellPropertyInput_t* input, CONVERGE_CellPropertyOutput_t* output);
void gas_mol_diffusivity(CONVERGE_CellPropertyInput_t* input, CONVERGE_CellPropertyOutput_t* output);
void absorption_coef_property(CONVERGE_AbsorptionCoefPropertyInput_t* input, CONVERGE_AbsorptionCoefPropertyOutput_t* output);
CONVERGE_shape* shape_nd;

// Simulation
// void constant_dt(CONVERGE_DtOutput_t* output);
void user_fsi(CONVERGE_FSIInput_t* input);
void user_fsi_dynamics(CONVERGE_FSIDynamicsInput_t* input);

// Lagrangian
void custom_drop_distort();
void custom_dynamic_spray_cone_angle(CONVERGE_DynamicSprayConeAngleInput_t* input, CONVERGE_DynamicSprayConeAngleOutput_t* output);
void custom_film_boiling(CONVERGE_FilmBoilingInput_t* input, CONVERGE_FilmBoilingOutput_t* output);
void custom_film_evap();
void custom_film_gradp(CONVERGE_FilmGradpInput_t* input);
void custom_film_jet(CONVERGE_FilmJetInput_t* input,CONVERGE_FilmJetOutput_t* output);
void custom_film_prop();
// void custom_film_splash_velocity(CONVERGE_FilmSplashVelocityInput_t* input, CONVERGE_FilmSplashVelocityOutput_t* output);
void custom_film_splash();
// void custom_film_strip();
void custom_film_velocity();
void custom_heat_trans_coeff_film(CONVERGE_HeatTransCoeffFilmInput_t* input,CONVERGE_HeatTransCoeffFilmOutput_t* output);
void custom_heat_trans_coeff_spray(CONVERGE_HeatTransCoeffSprayInput_t* input, CONVERGE_HeatTransCoeffSprayOutput_t* output);
void custom_ifpen_film_boiling(CONVERGE_IfpenFilmBoilingInput_t* input,CONVERGE_IfpenFilmBoilingOutput_t* output);
void custom_ifpen_film_thickness();
void custom_ifpen_spray_wall_heat_flux(CONVERGE_IfpenSprayWallHeatFluxInput_t* input, CONVERGE_IfpenSprayWallHeatFluxOutput_t* output);
// void custom_injector(CONVERGE_InjectorInput_t * input);
void custom_mass_trans_coeff_film(CONVERGE_MassTransCoeffFilmInput_t* input,CONVERGE_MassTransCoeffFilmOutput_t* output);
void custom_mass_trans_coeff_spray(CONVERGE_MassTransCoeffSprayInput_t* input, CONVERGE_MassTransCoeffSprayOutput_t* output);
void custom_nozzle(CONVERGE_NozzleInput_t* input, CONVERGE_NozzleOutput_t* output);
// void custom_nucleate_boiling(CONVERGE_NucleateBoilingInput_t* input, CONVERGE_NucleateBoilingOutput_t* output);
void custom_parcel_child(CONVERGE_ParcelChildInput_t* input);
void custom_parcel_inject(CONVERGE_ParcelInjectInput_t* input);
void custom_parcel_splash(CONVERGE_ParcelSplashInput_t* input);
void custom_parcel_strip(CONVERGE_ParcelStripInput_t* input);
void custom_rebound_velocity(CONVERGE_ReboundVelocityInput_t* input,CONVERGE_ReboundVelocityOutput_t* output);
void custom_regions_spray();
void custom_splash_crit(CONVERGE_SplashCritInput_t* input, CONVERGE_SplashCritOutput_t* output);
void custom_splash_mass(CONVERGE_SplashMassOutput_t* output);
void custom_splash_radius(CONVERGE_SplashRadiusInput_t* input,CONVERGE_SplashRadiusOutput_t* output);
void custom_spray_coalesce(CONVERGE_SprayCoalesceInput_t* input);
void custom_spray_collide(CONVERGE_SprayCollideInput_t* input);
void custom_spray_evap();
void custom_spray_gas_couple();
// void custom_spray_inject_profile(CONVERGE_SprayInjectProfileInput_t* input);
void custom_spray_kh(CONVERGE_SprayKHInput_t* input);
void custom_spray_rt(CONVERGE_SprayRTInput_t* input,CONVERGE_SprayRTOutput_t* output);
void custom_vof_boil_condensation();
// void custom_vof_cav_condensation();
// void custom_vof_contact_angle(CONVERGE_VofContactAngleInput_t* input, CONVERGE_VofContactAngleOutput_t* output);
void custom_vof_dissolved_gas();
// void spray_inject_custom(CONVERGE_SprayInjectInput_t* input, CONVERGE_SprayInjectOutput_t* output);
// void custom_film_sources();
// void custom_spray_boundary_inject(CONVERGE_SprayBoundaryInjectInput_t* input);
// void custom_spray_break();
// void custom_spray_main(CONVERGE_SprayMainInput_t* input);
// void custom_size_scale_flash_boiling(CONVERGE_SizeScaleFlashBoilingInput_t* input, CONVERGE_SizeScaleFlashBoilingOutput_t* output);
// void custom_liquid_spray_parcel_rebound_slide_velocity(CONVERGE_LiquidSprayParcelReboundSlideVelocityInput_t* input);
// void custom_solid_parcel_rebound_slide_velocity(CONVERGE_SolidParcelReboundSlideVelocityInput_t* input);
// void flux_limiter(CONVERGE_FluxLimiterInput_t* input, CONVERGE_FluxLimiterOutput_t* output);

// Turbulence
void custom_set_turb_prandtl_schmidt();
void custom_calc_wall_stress(CONVERGE_CalcWallStressInput_t* input);
void custom_calc_yplus(CONVERGE_CalcYPlusInput_t* input, CONVERGE_CalcYPlusOutput_t* output);
void turbulence_model_calc_eps_src();
void turbulence_model_calc_omega_src(CONVERGE_TurbModelCalcOmegaSrcInput_t* input);
void turbulence_model_calc_tke_src(CONVERGE_TurbModelCalcTkeSrcInput_t* input);
void turbulence_model_set_turb_visc();
void turbulence_model_setup();

// Combustion
void custom_calc_combust_ctc();
void custom_calc_combust_nox();
void custom_calc_combust_patch_rad();
void custom_calc_combust_patch_rad_2();
void custom_calc_combust_reinitialize_ctc();
void custom_calc_combust_shell();
void custom_calc_combust_soot();
void custom_calc_combust_tfm_efficiency();
void custom_calc_combust_tfm_thickening();
void custom_calc_equiv_ratio();
void custom_combust_model_setup();
void custom_combust_model();
void custom_combust_model_close();
void surface_rxn_2(CONVERGE_ReactionRateInput_t* input, CONVERGE_ReactionRateOutput_t* output);
void surface_rxn_5(CONVERGE_ReactionRateInput_t* input, CONVERGE_ReactionRateOutput_t* output);
void gas_rxn_617272(CONVERGE_ReactionRateInput_t* input, CONVERGE_ReactionRateOutput_t* output);
// void soot_reac_rate(CONVERGE_SootReactionRateInput_t* input, CONVERGE_SootReactionRateOutput_t* output);
void combust_adaptive_zone();
void combust_adaptive_zone_setup();
void combust_adaptive_zone_close();
void combust_rif();
void combust_rif_setup();
void user_mech_rate(CONVERGE_MechRateInput_t* input, CONVERGE_MechRateOutput_t* output);
void combust_sage_ode(CONVERGE_CombustSageOdeInput_t* input, CONVERGE_CombustSageOdeOutput_t* output);
void combust_sage_ode_init(CONVERGE_CombustSageOdeInitInput_t* input);
void combust_sage_ode_final();
void reaction_rate_multiplier();
// void custom_ignition_delay_init(CONVERGE_IgnitionDelayInitInput_t* input, CONVERGE_IgnitionDelayInitOutput_t* output);
// void custom_ignition_delay_final(CONVERGE_IgnitionDelayFinalOutput_t* output);
// void custom_ignition_delay(CONVERGE_IgnitionDelayInput_t* input, CONVERGE_IgnitionDelayOutput_t* output);
void custom_laminar_flame_speed_setup(CONVERGE_LaminarFlameSpeedOutput_t* output);
void custom_laminar_flame_speed(CONVERGE_LaminarFlameSpeedInput_t* input, CONVERGE_LaminarFlameSpeedOutput_t* output);
void custom_laminar_flame_speed_close(CONVERGE_LaminarFlameSpeedOutput_t* output);
// void combust_ecfm_flame_speed();
// void combust_ecfm_flame_speed_setup();
// void interpolated_chem_setup(CONVERGE_InterpolatedChemSetupInput_t* input, CONVERGE_InterpolatedChemSetupOutput_t* output);
// void interpolated_chem();
// void interpolated_chem_close();
// void interpolated_chem_get_cmean_eq(CONVERGE_InterpolatedChemGetCmeanEqInput_t* input, CONVERGE_InterpolatedChemGetCmeanEqOutput_t* output);
// void interpolated_chem_enthalpy_from_temp(CONVERGE_InterpolatedChemEnthalpyFromTempInput_t* input, CONVERGE_InterpolatedChemEnthalpyFromTempOutput_t* output);
// void interpolated_chem_interpolate_table(CONVERGE_InterpolatedChemInterpolateTableInput_t* input, CONVERGE_InterpolatedChemInterpolateTableOutput_t*);
// void interpolated_chem_var_index_map(CONVERGE_InterpolatedChemVarIndexMapInput_t* input, CONVERGE_InterpolatedChemVarIndexMapOutput_t* output);
void user_g_equation(CONVERGE_GEquationInput_t* input, CONVERGE_GEquationOutput_t* output);
void user_initialize_g_equation(CONVERGE_InitializeGEquationInput_t* input);
void user_reinitialize_g_equation(CONVERGE_ReInitializeGEquationInput_t* input, CONVERGE_ReInitializeGEquationOutput_t* output);
void user_set_grad_g_equal_one(CONVERGE_SetGradGEqualOneInput_t* input);
void user_g_eqn_spark_read(CONVERGE_GEquationSparkReadInput_t* input);
void user_g_eqn_spark_write(CONVERGE_GEquationSparkWriteInput_t* input);
void user_g_eqn_spark_read_restart(CONVERGE_GEquationSparkReadRestartInput_t* input);
void user_g_eqn_spark_write_restart(CONVERGE_GEquationSparkWriteRestartInput_t* input, CONVERGE_GEquationSparkWriteRestartOutput_t* output);
void user_test(CONVERGE_EcfmLftInput_t* input, CONVERGE_EcfmLftOutput_t* output);
void idaj(CONVERGE_EcfmSrcInput_t* input, CONVERGE_EcfmSrcOutput_t* output);

// Transport System
void generic_transport(CONVERGE_TransportSystemInput_t* input, CONVERGE_TransportSystemOutput_t* output);

// Heat Transfer
void amsden_1997(CONVERGE_HeatTransferInput_t* input, CONVERGE_HeatTransferOutput_t* output);
void custom_contact_resistance(CONVERGE_ContactResistanceInput_t* input, CONVERGE_ContactResistanceOutput_t* output);

// Profile
// void sine_profile_setup(CONVERGE_ProfileInput_t* input);
// void sine_profile(CONVERGE_ProfileInput_t* input);

// Electric potential
// void custom_epot_src();

// Radiation
void scatter_function(CONVERGE_ScatterFunctionInput_t* input, CONVERGE_ScatterFunctionOutput_t* output);

// Supporting onload/onclose functions
void custom_calc_combust_ctc_onload();
void custom_calc_combust_ctc_onclose();

void custom_calc_combust_shell_onload();
void custom_calc_combust_shell_onclose();

void custom_calc_combust_nox_onload();
void custom_calc_combust_nox_onclose();

void custom_calc_combust_soot_onload();
void custom_calc_combust_soot_onclose();

void custom_calc_equiv_ratio_onload();
void custom_calc_equiv_ratio_onclose();

void custom_calc_combust_tfm_thickening_onload();
void custom_calc_combust_tfm_thickening_onclose();

void custom_spray_evap_onload();
void custom_spray_evap_onclose();

void custom_set_turb_prandtl_schmidt_onload();
void custom_set_turb_prandtl_schmidt_onclose();

void custom_calc_wall_stress_onload();
void custom_calc_wall_stress_onclose();

void custom_calc_yplus_onload();
void custom_calc_yplus_onclose();

// void custom_vof_cav_condensation_onload();
// void custom_vof_cav_condensation_onclose();

void reg_sample_var_onload();

// void custom_epot_src_onload();
// void custom_epot_src_onclose();

void custom_film_sources_onload();
void custom_film_sources_onclose();

void release_shape_onclose();

void cvg_cloud_properties();
void spray_env();
// void custom_liquid_spray_parcel_rebound_slide_velocity_onload();
// void custom_liquid_spray_parcel_rebound_slide_velocity_onclose();

// void custom_solid_parcel_rebound_slide_velocity_onload();
// void custom_solid_parcel_rebound_slide_velocity_onclose();

void register_transport_var();

#endif
