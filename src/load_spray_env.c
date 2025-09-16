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
#include "lagrangian/env.h"
#include <counter.h>
#include <CONVERGE/udf.h>
#include <pthread.h>
//Global variables to index parcels and clouds - need to initialize here
 user_parcel_counter = 0; //First index
 user_cloud_counter = 0; //First index
 update_cloud_counter_flag = 0;   //Flag to tell program to tick cloud_counter
/** Load the spray environment
 */
CONVERGE_ONLOAD(spray_env, IN(CONVERGE_VOID))
{
   
   // Register a simple double data parcel field
   CONVERGE_variable_register("user_lag_var", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("r_bubble", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("v_bubble", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("v_bubble_tm1", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST); 
   CONVERGE_variable_register("v_drop", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST); 
   CONVERGE_variable_register("r_bubble_0", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("r_bubble_tm1",CONVERGE_DOUBLE,DEFAULT_PARCEL_VARIABLE_SETTINGS,END_ARG_LIST);
   CONVERGE_variable_register("r_drop_0", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("r_therm", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("omega", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("omega_tm1", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("int_omega", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("m0",        CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("child_index", CONVERGE_INT, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);

   CONVERGE_variable_register("eta_drop", CONVERGE_DOUBLE, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("thermal_breakup_flag", CONVERGE_INT, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);

   // Register a simple int data parcel field
   CONVERGE_variable_register("user_lag_var_i", CONVERGE_INT, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("dgre_cycle_count", CONVERGE_INT, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("parcel_index", CONVERGE_INT, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("cloud_index", CONVERGE_INT, DEFAULT_PARCEL_VARIABLE_SETTINGS, END_ARG_LIST);
   CONVERGE_variable_register("tbt",CONVERGE_INT,DEFAULT_PARCEL_VARIABLE_SETTINGS,END_ARG_LIST);
   CONVERGE_variable_register("is_child",CONVERGE_INT,DEFAULT_PARCEL_VARIABLE_SETTINGS,END_ARG_LIST);
   CONVERGE_variable_register("pbt",CONVERGE_INT,DEFAULT_PARCEL_VARIABLE_SETTINGS,END_ARG_LIST);

   // User defined component names, overrides automatic nameing for CONVERGE_VEC3
   const char *user_lag_var_v3_comp_names[] = {"user_lag_var0", "user_lag_var1", "user_lag_var2"};
   CONVERGE_variable_register(
      "user_lag_var_v3",
      CONVERGE_DOUBLE,
      DEFAULT_PARCEL_VARIABLE_SETTINGS,
      // Double wtih dimension 3 may also be simply CONVERGE_VEC3.
      // The purpose of specifying the dimension manually is to allow for component names to be specified separately
      "dimension",
      3,
      // Pass the component names after the "component_names" parameter
      "component_names",
      user_lag_var_v3_comp_names,
      END_ARG_LIST);
   CONVERGE_variable_register(
      "user_lag_var_v3b",
      // CONVERGE_VEC3 will automattically append _1/_2/_3 to the end of the variable name
      CONVERGE_VEC3,
      DEFAULT_PARCEL_VARIABLE_SETTINGS,
      // The purpose of specifying the dimension manually here is to demonstrate it is permitted to mix dimension with any CONVERGE_APIType
      "dimension",
      3,
      END_ARG_LIST);
      CONVERGE_variable_register(
         "child_uu",
         // CONVERGE_VEC3 will automattically append _1/_2/_3 to the end of the variable name
         CONVERGE_VEC3,
         DEFAULT_PARCEL_VARIABLE_SETTINGS,
         // The purpose of specifying the dimension manually here is to demonstrate it is permitted to mix dimension with any CONVERGE_APIType
         "dimension",
         3,
         END_ARG_LIST);
   // Get dynamic IDs to Lagrangian Cloud fields
   USER_LAG_VAR    = CONVERGE_lagrangian_field_id("user_lag_var");
   USER_LAG_VARi   = CONVERGE_lagrangian_field_id("user_lag_var_i");
   USER_LAG_VARv3  = CONVERGE_lagrangian_field_id("user_lag_var_v3");
   CHILD_UU = CONVERGE_lagrangian_field_id("child_uu");
   USER_LAG_VARv3b = CONVERGE_lagrangian_field_id("user_lag_var_v3b");
   R_BUBBLE = CONVERGE_lagrangian_field_id("r_bubble");
   V_BUBBLE = CONVERGE_lagrangian_field_id("v_bubble");
   V_B_TM1  = CONVERGE_lagrangian_field_id("v_bubble_tm1");
   V_DROP =    CONVERGE_lagrangian_field_id("v_drop");
   R_B_0    = CONVERGE_lagrangian_field_id("r_bubble_0");
   R_D_0    = CONVERGE_lagrangian_field_id("r_drop_0");
   R_B_TM1 = CONVERGE_lagrangian_field_id("r_bubble_tm1");
   R_THERM = CONVERGE_lagrangian_field_id("r_therm");
   OMEGA    = CONVERGE_lagrangian_field_id("omega");
   OMEGA_TM1    = CONVERGE_lagrangian_field_id("omega_tm1");
   INT_OMEGA    = CONVERGE_lagrangian_field_id("int_omega");
   ETA      = CONVERGE_lagrangian_field_id("eta_drop");   
   USER_LAG_VARi   = CONVERGE_lagrangian_field_id("user_lag_var_i");
   THERMAL_BREAKUP_FLAG = CONVERGE_lagrangian_field_id("thermal_breakup_flag");
   TBT      = CONVERGE_lagrangian_field_id("tbt");
   IS_CHILD = CONVERGE_lagrangian_field_id("is_child");
   CHILD_INDEX = CONVERGE_lagrangian_field_id("child_index");
   PBT      = CONVERGE_lagrangian_field_id("pbt");
   DGRE_COUNT = CONVERGE_lagrangian_field_id("dgre_cycle_count");
   CLOUD_INDEX = CONVERGE_lagrangian_field_id("cloud_index");
   PARCEL_INDEX = CONVERGE_lagrangian_field_id("parcel_index");
   USER_LAG_VARv3  = CONVERGE_lagrangian_field_id("user_lag_var_v3");
   USER_LAG_VARv3b = CONVERGE_lagrangian_field_id("user_lag_var_v3b");
   M0 = CONVERGE_lagrangian_field_id("m0");


   LAGRANGIAN_FROM_INJECTOR = CONVERGE_lagrangian_field_id("LAGRANGIAN_FROM_INJECTOR");
   LAGRANGIAN_FROM_INJECTOR_TYPE = CONVERGE_lagrangian_field_id("LAGRANGIAN_FROM_INJECTOR_TYPE");
   LAGRANGIAN_ON_TRIANGLE   = CONVERGE_lagrangian_field_id("LAGRANGIAN_ON_TRIANGLE");
   LAGRANGIAN_FROM_NOZZLE   = CONVERGE_lagrangian_field_id("LAGRANGIAN_FROM_NOZZLE");
   LAGRANGIAN_FILM_FLAG     = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_FLAG");

   LAGRANGIAN_JUST_HIT        = CONVERGE_lagrangian_field_id("LAGRANGIAN_JUST_HIT");
   LAGRANGIAN_JUST_HIT_LEIDEN = CONVERGE_lagrangian_field_id("LAGRANGIAN_JUST_HIT_LEIDEN");
   LAGRANGIAN_ISTHICK         = CONVERGE_lagrangian_field_id("LAGRANGIAN_ISTHICK");

   LAGRANGIAN_REL_VEL = CONVERGE_lagrangian_field_id("LAGRANGIAN_REL_VEL");
   LAGRANGIAN_UPRIME  = CONVERGE_lagrangian_field_id("LAGRANGIAN_UPRIME");
   LAGRANGIAN_UU      = CONVERGE_lagrangian_field_id("LAGRANGIAN_UU");
   LAGRANGIAN_UU_TM1  = CONVERGE_lagrangian_field_id("LAGRANGIAN_UU_TM1");
   LAGRANGIAN_XX      = CONVERGE_lagrangian_field_id("LAGRANGIAN_XX");
   LAGRANGIAN_XX_TM1  = CONVERGE_lagrangian_field_id("LAGRANGIAN_XX_TM1");

   LAGRANGIAN_V_NU               = CONVERGE_lagrangian_field_id("LAGRANGIAN_V_NU");
   LAGRANGIAN_V_SH               = CONVERGE_lagrangian_field_id("LAGRANGIAN_V_SH");
   LAGRANGIAN_TEMP               = CONVERGE_lagrangian_field_id("LAGRANGIAN_TEMP");
   LAGRANGIAN_TEMP_TM1           = CONVERGE_lagrangian_field_id("LAGRANGIAN_TEMP_TM1");
   LAGRANGIAN_TEMP_STARM1        = CONVERGE_lagrangian_field_id("LAGRANGIAN_TEMP_STARM1");
   LAGRANGIAN_REY_NUM            = CONVERGE_lagrangian_field_id("LAGRANGIAN_REY_NUM");
   LAGRANGIAN_REL_VEL_MAG        = CONVERGE_lagrangian_field_id("LAGRANGIAN_REL_VEL_MAG");
   LAGRANGIAN_RADIUS             = CONVERGE_lagrangian_field_id("LAGRANGIAN_RADIUS");
   LAGRANGIAN_RADIUS_TM1         = CONVERGE_lagrangian_field_id("LAGRANGIAN_RADIUS_TM1");
   LAGRANGIAN_PARENT             = CONVERGE_lagrangian_field_id("LAGRANGIAN_PARENT");
   LAGRANGIAN_DENSITY            = CONVERGE_lagrangian_field_id("LAGRANGIAN_DENSITY");
   LAGRANGIAN_DENSITY_TM1        = CONVERGE_lagrangian_field_id("LAGRANGIAN_DENSITY_TM1");
   LAGRANGIAN_GAS_DENSITY        = CONVERGE_lagrangian_field_id("LAGRANGIAN_GAS_DENSITY");
   LAGRANGIAN_MFRAC              = CONVERGE_lagrangian_field_id("LAGRANGIAN_MFRAC");
   LAGRANGIAN_MFRAC_TM1          = CONVERGE_lagrangian_field_id("LAGRANGIAN_MFRAC_TM1");
   LAGRANGIAN_NUM_DROP           = CONVERGE_lagrangian_field_id("LAGRANGIAN_NUM_DROP");
   LAGRANGIAN_SURF_TEMP          = CONVERGE_lagrangian_field_id("LAGRANGIAN_SURF_TEMP");
   LAGRANGIAN_TBREAK_KH          = CONVERGE_lagrangian_field_id("LAGRANGIAN_TBREAK_KH");
   LAGRANGIAN_SHED_NUM_DROP      = CONVERGE_lagrangian_field_id("LAGRANGIAN_SHED_NUM_DROP");
   LAGRANGIAN_SHED_MASS          = CONVERGE_lagrangian_field_id("LAGRANGIAN_SHED_MASS");
   LAGRANGIAN_SACTIVE            = CONVERGE_lagrangian_field_id("LAGRANGIAN_SACTIVE");
   LAGRANGIAN_SACTIVE_TM1        = CONVERGE_lagrangian_field_id("LAGRANGIAN_SACTIVE_TM1");
   LAGRANGIAN_SURF_TEN           = CONVERGE_lagrangian_field_id("LAGRANGIAN_SURF_TEN");
   LAGRANGIAN_VISCOSITY          = CONVERGE_lagrangian_field_id("LAGRANGIAN_VISCOSITY");
   LAGRANGIAN_DISTORT            = CONVERGE_lagrangian_field_id("LAGRANGIAN_DISTORT");
   LAGRANGIAN_DISTORT_DOT        = CONVERGE_lagrangian_field_id("LAGRANGIAN_DISTORT_DOT");
   LAGRANGIAN_DM_DT              = CONVERGE_lagrangian_field_id("LAGRANGIAN_DM_DT");
   LAGRANGIAN_DRDT               = CONVERGE_lagrangian_field_id("LAGRANGIAN_DRDT");
   LAGRANGIAN_WALL_HEAT_EXCHANGE = CONVERGE_lagrangian_field_id("LAGRANGIAN_WALL_HEAT_EXCHANGE");
   LAGRANGIAN_L_RR               = CONVERGE_lagrangian_field_id("LAGRANGIAN_L_RR");
   LAGRANGIAN_L_RC               = CONVERGE_lagrangian_field_id("LAGRANGIAN_L_RC");
   LAGRANGIAN_L_TEMP1            = CONVERGE_lagrangian_field_id("LAGRANGIAN_L_TEMP1");
   LAGRANGIAN_L_TEMP2            = CONVERGE_lagrangian_field_id("LAGRANGIAN_L_TEMP2");
   LAGRANGIAN_FILM_SHED          = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_SHED");
   LAGRANGIAN_TBREAK_RT          = CONVERGE_lagrangian_field_id("LAGRANGIAN_TBREAK_RT");
   LAGRANGIAN_SURF_TEMP_TM1      = CONVERGE_lagrangian_field_id("LAGRANGIAN_SURF_TEMP_TM1");
   LAGRANGIAN_NUM_DROP_TM1       = CONVERGE_lagrangian_field_id("LAGRANGIAN_NUM_DROP_TM1");
   LAGRANGIAN_FILM_ENERGY        = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_ENERGY");
   LAGRANGIAN_T_TURB             = CONVERGE_lagrangian_field_id("LAGRANGIAN_T_TURB");
   LAGRANGIAN_T_TURB_ACCUM       = CONVERGE_lagrangian_field_id("LAGRANGIAN_T_TURB_ACCUM");
   LAGRANGIAN_FILM_THICKNESS     = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_THICKNESS");
   LAGRANGIAN_AREA_REDUCTION     = CONVERGE_lagrangian_field_id("LAGRANGIAN_AREA_REDUCTION");
   LAGRANGIAN_TKE0               = CONVERGE_lagrangian_field_id("LAGRANGIAN_TKE0");
   LAGRANGIAN_EPS0               = CONVERGE_lagrangian_field_id("LAGRANGIAN_EPS0");
   LAGRANGIAN_LIFETIME           = CONVERGE_lagrangian_field_id("LAGRANGIAN_LIFETIME");

   LAGRANGIAN_FILM_THICKNESS_TM1 = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_THICKNESS_TM1");
   LAGRANGIAN_AREA_IN_FILM       = CONVERGE_lagrangian_field_id("LAGRANGIAN_AREA_IN_FILM");

   LAGRANGIAN_FILM_ACCUM_BIT_FLAG = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_ACCUM_BIT_FLAG");
   LAGRANGIAN_FILM_ACCUM_PLUS_BIT_FLAG = CONVERGE_lagrangian_field_id("LAGRANGIAN_FILM_ACCUM_PLUS_BIT_FLAG");

   //Get dynamic IDs for solid parcel cloud
   SOLID_PARCEL_FROM_INJECTOR = CONVERGE_solid_parcel_field_id("LAGRANGIAN_FROM_INJECTOR");
   SOLID_PARCEL_FROM_INJECTOR_TYPE = CONVERGE_solid_parcel_field_id("LAGRANGIAN_FROM_INJECTOR_TYPE");
   SOLID_PARCEL_ON_TRIANGLE   = CONVERGE_solid_parcel_field_id("LAGRANGIAN_ON_TRIANGLE");
   SOLID_PARCEL_FROM_NOZZLE   = CONVERGE_solid_parcel_field_id("LAGRANGIAN_FROM_NOZZLE");
   SOLID_PARCEL_FILM_FLAG     = CONVERGE_solid_parcel_field_id("LAGRANGIAN_FILM_FLAG");

   SOLID_PARCEL_JUST_HIT        = CONVERGE_solid_parcel_field_id("LAGRANGIAN_JUST_HIT");

   SOLID_PARCEL_REL_VEL = CONVERGE_solid_parcel_field_id("LAGRANGIAN_REL_VEL");
   SOLID_PARCEL_UPRIME  = CONVERGE_solid_parcel_field_id("LAGRANGIAN_UPRIME");
   SOLID_PARCEL_UU      = CONVERGE_solid_parcel_field_id("LAGRANGIAN_UU");
   SOLID_PARCEL_UU_TM1  = CONVERGE_solid_parcel_field_id("LAGRANGIAN_UU_TM1");
   SOLID_PARCEL_XX      = CONVERGE_solid_parcel_field_id("LAGRANGIAN_XX");
   SOLID_PARCEL_XX_TM1  = CONVERGE_solid_parcel_field_id("LAGRANGIAN_XX_TM1");

   SOLID_PARCEL_V_NU               = CONVERGE_solid_parcel_field_id("LAGRANGIAN_V_NU");
   SOLID_PARCEL_TEMP               = CONVERGE_solid_parcel_field_id("LAGRANGIAN_TEMP");
   SOLID_PARCEL_TEMP_TM1           = CONVERGE_solid_parcel_field_id("LAGRANGIAN_TEMP_TM1");
   SOLID_PARCEL_TEMP_STARM1        = CONVERGE_solid_parcel_field_id("LAGRANGIAN_TEMP_STARM1");
   SOLID_PARCEL_REY_NUM            = CONVERGE_solid_parcel_field_id("LAGRANGIAN_REY_NUM");
   SOLID_PARCEL_REL_VEL_MAG        = CONVERGE_solid_parcel_field_id("LAGRANGIAN_REL_VEL_MAG");
   SOLID_PARCEL_RADIUS             = CONVERGE_solid_parcel_field_id("LAGRANGIAN_RADIUS");
   SOLID_PARCEL_RADIUS_TM1         = CONVERGE_solid_parcel_field_id("LAGRANGIAN_RADIUS_TM1");
   SOLID_PARCEL_DENSITY            = CONVERGE_solid_parcel_field_id("LAGRANGIAN_DENSITY");
   SOLID_PARCEL_DENSITY_TM1        = CONVERGE_solid_parcel_field_id("LAGRANGIAN_DENSITY_TM1");
   SOLID_PARCEL_MFRAC              = CONVERGE_solid_parcel_field_id("LAGRANGIAN_MFRAC");
   SOLID_PARCEL_MFRAC_TM1          = CONVERGE_solid_parcel_field_id("LAGRANGIAN_MFRAC_TM1");
   SOLID_PARCEL_NUM_DROP           = CONVERGE_solid_parcel_field_id("LAGRANGIAN_NUM_DROP");
   SOLID_PARCEL_SURF_TEMP          = CONVERGE_solid_parcel_field_id("LAGRANGIAN_SURF_TEMP");
   SOLID_PARCEL_VISCOSITY          = CONVERGE_solid_parcel_field_id("LAGRANGIAN_VISCOSITY");
   SOLID_PARCEL_DISTORT            = CONVERGE_solid_parcel_field_id("LAGRANGIAN_DISTORT");
   SOLID_PARCEL_SURF_TEMP_TM1      = CONVERGE_solid_parcel_field_id("LAGRANGIAN_SURF_TEMP_TM1");
   SOLID_PARCEL_NUM_DROP_TM1       = CONVERGE_solid_parcel_field_id("LAGRANGIAN_NUM_DROP_TM1");
   SOLID_PARCEL_T_TURB             = CONVERGE_solid_parcel_field_id("LAGRANGIAN_T_TURB");
   SOLID_PARCEL_T_TURB_ACCUM       = CONVERGE_solid_parcel_field_id("LAGRANGIAN_T_TURB_ACCUM");

   // Get Dynamic IDs for Injectors parameters
   INJECTOR_INJECTED_PARCEL_INDEX        = CONVERGE_get_parameter_id("injector.injected_parcel_index");
   INJECTOR_INJECTED_PARCEL_TYPE         = CONVERGE_get_parameter_id("injector.injected_parcel_type");
   INJECTOR_INJECTED_NUM_PARCEL_SPECIES  = CONVERGE_get_parameter_id("injector.injected_num_parcel_species");
   INJECTOR_CONE_DISTRIBUTION_FLAG       = CONVERGE_get_parameter_id("injector.cone_distribution_flag");
   INJECTOR_DYNAMIC_SPRAY_CONE_ANGLE_FLAG= CONVERGE_get_parameter_id("injector.dynamic_spray_cone_angle_flag");
   INJECTOR_PENET_FRAC                   = CONVERGE_get_parameter_id("injector.penet_frac");
   INJECTOR_VAPOR_PENET_FRAC             = CONVERGE_get_parameter_id("injector.vapor_penet_frac");
   INJECTOR_PENET_BIN_SIZE               = CONVERGE_get_parameter_id("injector.penet_bin_size");
   INJECTOR_MFRAC                        = CONVERGE_get_parameter_id("injector.mfrac");
   INJECTOR_START_INJECT_IS_FILE         = CONVERGE_get_parameter_id("injector.start_inject_is_file");
   INJECTOR_DUR_INJECT_IS_FILE           = CONVERGE_get_parameter_id("injector.dur_inject_is_file");
   INJECTOR_MASS_INJECT_IS_FILE          = CONVERGE_get_parameter_id("injector.mass_inject_is_file");
   INJECTOR_NOZZLE_INIT_FLAG             = CONVERGE_get_parameter_id("injector.nozzle_init_flag");
   INJECTOR_SPRAY_INJECT_BC_FLAG         = CONVERGE_get_parameter_id("injector.spray_inject_bc_flag");
   INJECTOR_ELSA_FLAG                    = CONVERGE_get_parameter_id("injector.elsa_flag");
   INJECTOR_TEMP_FLAG                    = CONVERGE_get_parameter_id("injector.temp_flag");
   INJECTOR_TKE_FLAG                     = CONVERGE_get_parameter_id("injector.tke_flag");
   INJECTOR_EPS_FLAG                     = CONVERGE_get_parameter_id("injector.eps_flag");
   INJECTOR_INIT_CELL_TURB_FLAG          = CONVERGE_get_parameter_id("injector.init_cell_turb_flag");
   INJECTOR_KHACT_NOZZLE_FLOW_FLAG       = CONVERGE_get_parameter_id("injector.khact_nozzle_flow_flag");
   INJECTOR_KH_NEW_PARCEL_FLAG           = CONVERGE_get_parameter_id("injector.kh_new_parcel_flag");
   INJECTOR_KH_NO_ENLARGE_FLAG           = CONVERGE_get_parameter_id("injector.kh_no_enlarge_flag");
   INJECTOR_INJECT_DISTRIBUTION_FLAG     = CONVERGE_get_parameter_id("injector.inject_distribution_flag");
   INJECTOR_RT_FLAG                      = CONVERGE_get_parameter_id("injector.rt_flag");
   INJECTOR_RT_DISTRIBUTION_FLAG         = CONVERGE_get_parameter_id("injector.rt_distribution_flag");
   INJECTOR_TAB_FLAG                     = CONVERGE_get_parameter_id("injector.tab_flag");
   INJECTOR_TAB_DISTRIBUTION_FLAG        = CONVERGE_get_parameter_id("injector.tab_distribution_flag");
   INJECTOR_LISA_FLAG                    = CONVERGE_get_parameter_id("injector.lisa_flag");
   INJECTOR_DISCHARGE_COEFF_INPUT_FLAG   = CONVERGE_get_parameter_id("injector.discharge_coeff_input_flag");
   INJECTOR_LISA_DISTRIBUTION_FLAG       = CONVERGE_get_parameter_id("injector.lisa_distribution_flag");
   INJECTOR_RATE_SHAPE_INPUT_FLAG        = CONVERGE_get_parameter_id("injector.rate_shape_input_flag");
   INJECTOR_DISCHARGE_COEFF_FLAG         = CONVERGE_get_parameter_id("injector.discharge_coeff_flag");
   INJECTOR_INJECT_INDEX                 = CONVERGE_get_parameter_id("injector.inject_index");
   INJECTOR_NUM_PARCELS_PER_NOZZLE       = CONVERGE_get_parameter_id("injector.num_parcels_per_nozzle");
   INJECTOR_NUM_NOZZLES                  = CONVERGE_get_parameter_id("injector.num_nozzles");
   INJECTOR_STREAM_INDEX                 = CONVERGE_get_parameter_id("injector.stream_index");
   INJECTOR_NUM_INJECT_PARCELS           = CONVERGE_get_parameter_id("injector.num_inject_parcels");
   INJECTOR_CONE_FLAG                    = CONVERGE_get_parameter_id("injector.cone_flag");
   INJECTOR_TEMPORAL_TYPE                = CONVERGE_get_parameter_id("injector.temporal_type");
   INJECTOR_VELOCITY_COEFF_INPUT_FLAG    = CONVERGE_get_parameter_id("injector.velocity_coeff_input_flag");
   INJECTOR_POLAR_CPY_NUM                = CONVERGE_get_parameter_id("injector.polar_cpy_num");
   INJECTOR_Q_RR                         = CONVERGE_get_parameter_id("injector.q_rr");
   INJECTOR_GAMMA_RR_X_DIST              = CONVERGE_get_parameter_id("injector.gamma_rr_x_dist");
   INJECTOR_LISA_INJECTION_PRES          = CONVERGE_get_parameter_id("injector.lisa_injection_pressure");
   INJECTOR_LISA_KV                      = CONVERGE_get_parameter_id("injector.lisa_kv");
   INJECTOR_CYCLIC_PERIOD                = CONVERGE_get_parameter_id("injector.cyclic_period");
   INJECTOR_INJECT_START_TIME            = CONVERGE_get_parameter_id("injector.inject_start_time");
   INJECTOR_INJECT_DURATION              = CONVERGE_get_parameter_id("injector.inject_duration");
   INJECTOR_END_INJECT                   = CONVERGE_get_parameter_id("injector.end_inject");
   INJECTOR_INJECT_MASS                  = CONVERGE_get_parameter_id("injector.inject_mass");
   INJECTOR_INJECT_TEMP                  = CONVERGE_get_parameter_id("injector.inject_temperature");
   INJECTOR_AREA                         = CONVERGE_get_parameter_id("injector.tot_area");
   INJECTOR_KH_NEW_PARCEL_CUTOFF         = CONVERGE_get_parameter_id("injector.kh_new_parcel_cutoff");
   INJECTOR_KH_SHED_FACTOR               = CONVERGE_get_parameter_id("injector.kh_shed_factor");
   INJECTOR_KH_BALPHA                    = CONVERGE_get_parameter_id("injector.kh_balpha");
   INJECTOR_KH_CONST1                    = CONVERGE_get_parameter_id("injector.kh_const1");
   INJECTOR_KH_CONST2                    = CONVERGE_get_parameter_id("injector.kh_const2");
   INJECTOR_RT_LENGTH_CONST              = CONVERGE_get_parameter_id("injector.rt_length_const");
   INJECTOR_RT_CONST3                    = CONVERGE_get_parameter_id("injector.rt_const3");
   INJECTOR_RT_CONST2                    = CONVERGE_get_parameter_id("injector.rt_const2");
   INJECTOR_KHACT_TURB_KC                = CONVERGE_get_parameter_id("injector.khact_turb_kc");
   INJECTOR_KHACT_TURB_KE                = CONVERGE_get_parameter_id("injector.khact_turb_ke");
   INJECTOR_KHACT_TURB_S                 = CONVERGE_get_parameter_id("injector.khact_turb_s");
   INJECTOR_KHACT_C_TCAV                 = CONVERGE_get_parameter_id("injector.khact_c_tcav");
   INJECTOR_TIME_OFFSET                  = CONVERGE_get_parameter_id("injector.time_offset");
   INJECTOR_MASS_PER_PARCEL              = CONVERGE_get_parameter_id("injector.mass_per_parcel");
   INJECTOR_VELOCITY                     = CONVERGE_get_parameter_id("injector.velocity");
   INJECTOR_VELOCITY_TM1                 = CONVERGE_get_parameter_id("injector.velocity_tm1");
   INJECTOR_VELOCITY_OUT_OLD             = CONVERGE_get_parameter_id("injector.velocity_out_old");
   INJECTOR_VELOCITY_OUT_NEW             = CONVERGE_get_parameter_id("injector.velocity_out_new");
   INJECTOR_VELOCITY_COEFF               = CONVERGE_get_parameter_id("injector.velocity_coeff");
   INJECTOR_MASS_TM1                     = CONVERGE_get_parameter_id("injector.mass_tm1");
   INJECTOR_MASS_TM2                     = CONVERGE_get_parameter_id("injector.mass_tm2");
   INJECTOR_DISCHARGE_COEFF              = CONVERGE_get_parameter_id("injector.discharge_coeff");
   INJECTOR_RATE_SHAPE                   = CONVERGE_get_parameter_id("injector.rate_shape");
   INJECTOR_ANGLE_XY_INJ                 = CONVERGE_get_parameter_id("injector.angle_xy_inj");
   INJECTOR_ANGLE_XZ_INJ                 = CONVERGE_get_parameter_id("injector.angle_xz_inj");
   INJECTOR_ANGLE_CLOCK_INJ              = CONVERGE_get_parameter_id("injector.angle_clock_inj");
   INJECTOR_SWIRL_FRAC                   = CONVERGE_get_parameter_id("injector.swirl_frac");
   INJECTOR_TOT_INJECTED_MASS            = CONVERGE_get_parameter_id("injector.tot_inj_mass");
   INJECTOR_TOT_INJECTED_MASS_TM         = CONVERGE_get_parameter_id("injector.tot_inj_mass_tm");
   INJECTOR_AZIMUTH_ANGLE_START          = CONVERGE_get_parameter_id("injector.azimuth_angle_start");
   INJECTOR_AZIMUTH_ANGLE_END            = CONVERGE_get_parameter_id("injector.azimuth_angle_end");
   INJECTOR_SCALE_TEMP                   = CONVERGE_get_parameter_id("injector.scale_temperature");
   INJECTOR_OFFSET_TEMP                  = CONVERGE_get_parameter_id("injector.offset_temperature");
   INJECTOR_INJECT_TKE                   = CONVERGE_get_parameter_id("injector.inject_tke");
   INJECTOR_SCALE_TKE                    = CONVERGE_get_parameter_id("injector.scale_tke");
   INJECTOR_OFFSET_TKE                   = CONVERGE_get_parameter_id("injector.offset_tke");
   INJECTOR_INJECT_EPS                   = CONVERGE_get_parameter_id("injector.inject_epsilon");
   INJECTOR_SCALE_EPS                    = CONVERGE_get_parameter_id("injector.scale_epsilon");
   INJECTOR_OFFSET_EPS                   = CONVERGE_get_parameter_id("injector.offset_epsilon");
   INJECTOR_SWIRLER_MASS_FLOW_RATE       = CONVERGE_get_parameter_id("injector.swirler_mass_flow_rate");
   INJECTOR_SWIRLER_MEAN_ANGLE           = CONVERGE_get_parameter_id("injector.swirler_mean_angle");
   INJECTOR_SWIRLER_RADIUS               = CONVERGE_get_parameter_id("injector.swirler_radius");
   INJECTOR_VOF_SPRAY_MASS_PER_PARCEL    = CONVERGE_get_parameter_id("injector.vof_spray_mass_per_parcel");
   INJECTOR_VOF_SPRAY_LIQ_FRAC_THRESHOLD = CONVERGE_get_parameter_id("injector.vof_spray_liq_frac_threshold");
   INJECTOR_X_CEN                        = CONVERGE_get_parameter_id("injector.x_cen_inj");
   INJECTOR_AXI_VEC                      = CONVERGE_get_parameter_id("injector.axi_vec");
   INJECTOR_NORM_VEC                     = CONVERGE_get_parameter_id("injector.norm_vec");
   INJECTOR_OTHER_VEC                    = CONVERGE_get_parameter_id("injector.other_vec");

   // Get Dynamic IDs for Nozzle parameters
   NOZZLE_AXIAL_VEC                = CONVERGE_get_parameter_id("nozzle.axial_vec");
   NOZZLE_NORMAL_VEC               = CONVERGE_get_parameter_id("nozzle.normal_vec");
   NOZZLE_OTHER_VEC                = CONVERGE_get_parameter_id("nozzle.other_vec");
   NOZZLE_AREA                     = CONVERGE_get_parameter_id("nozzle.area_noz");
   NOZZLE_DIAMETER                 = CONVERGE_get_parameter_id("nozzle.diameter");
   NOZZLE_LENGTH                   = CONVERGE_get_parameter_id("nozzle.length");
   NOZZLE_SMD_DISTRIBUTION         = CONVERGE_get_parameter_id("nozzle.smd_distribution");
   NOZZLE_AMP_DISTORT              = CONVERGE_get_parameter_id("nozzle.amp_distort");
   NOZZLE_CONE_ANGLE               = CONVERGE_get_parameter_id("nozzle.cone_angle");
   NOZZLE_THICKNESS                = CONVERGE_get_parameter_id("nozzle.thickness");
   NOZZLE_RADIAL_DISTANCE          = CONVERGE_get_parameter_id("nozzle.radial_distance");
   NOZZLE_AXIAL_DISTANCE           = CONVERGE_get_parameter_id("nozzle.axial_distance");
   NOZZLE_POSITION_VEC             = CONVERGE_get_parameter_id("nozzle.position_vec");
   NOZZLE_XX_VEC                   = CONVERGE_get_parameter_id("nozzle.xx_vec");
   NOZZLE_YY_VEC                   = CONVERGE_get_parameter_id("nozzle.yy_vec");
   NOZZLE_ZZ_VEC                   = CONVERGE_get_parameter_id("nozzle.zz_vec");
   NOZZLE_POSITION                 = CONVERGE_get_parameter_id("nozzle.position");
   NOZZLE_X                        = CONVERGE_get_parameter_id("nozzle.x");
   NOZZLE_Y                        = CONVERGE_get_parameter_id("nozzle.y");
   NOZZLE_Z                        = CONVERGE_get_parameter_id("nozzle.z");
   NOZZLE_MAX_DISTANCE             = CONVERGE_get_parameter_id("nozzle.max_distance");
   NOZZLE_PENETRATION              = CONVERGE_get_parameter_id("nozzle.penetration");
   NOZZLE_VAPOR_PENETRATION        = CONVERGE_get_parameter_id("nozzle.vapor_penetration");
   NOZZLE_ECN_PENETRATION          = CONVERGE_get_parameter_id("nozzle.ecn_penetration");
   NOZZLE_TOTAL_MASS               = CONVERGE_get_parameter_id("nozzle.total_mass");
   NOZZLE_RADIUS_INJECT            = CONVERGE_get_parameter_id("nozzle.radius_inject");
   NOZZLE_VOL_COUNT                = CONVERGE_get_parameter_id("nozzle.vol_count");
   NOZZLE_CONE_INPUT_FLAG          = CONVERGE_get_parameter_id("nozzle.cone_noz_input_flag");
   NOZZLE_DIAMETER_INPUT_FLAG      = CONVERGE_get_parameter_id("nozzle.noz_diameter_input_flag");
   NOZZLE_RADIUS_INJECT_INPUT_FLAG = CONVERGE_get_parameter_id("nozzle.radius_inject_input_flag");
   NOZZLE_NOZ_THETA                    = CONVERGE_get_parameter_id("nozzle.noz_theta");
   NOZZLE_NOZ_ANGLE_XY             = CONVERGE_get_parameter_id("nozzle.noz_angle_xy");
   NOZZLE_NOZ_ANGLE_XZ             = CONVERGE_get_parameter_id("nozzle.noz_angle_xz");
   
   // Get Dynamic IDs for boundary injector parameters
   BOUNDARY_INJECTOR_CYCLIC_PERIOD = CONVERGE_get_parameter_id("boundary_injector.cyclic_period");
   BOUNDARY_INJECTOR_END_INJECT = CONVERGE_get_parameter_id("boundary_injector.end_inject");
   BOUNDARY_INJECTOR_INJECT_START_TIME = CONVERGE_get_parameter_id("boundary_injector.inject_start_time");
   BOUNDARY_INJECTOR_INJECT_DURATION = CONVERGE_get_parameter_id("boundary_injector.inject_duration");
   BOUNDARY_INJECTOR_VELOCITY = CONVERGE_get_parameter_id("boundary_injector.velocity");
   BOUNDARY_INJECTOR_VELOCITY_TM1 = CONVERGE_get_parameter_id("boundary_injector.velocity_tm1");
   BOUNDARY_INJECTOR_TOT_INJ_MASS = CONVERGE_get_parameter_id("boundary_injector.tot_inj_mass");
   BOUNDARY_INJECTOR_TOT_INJ_MASS_TM = CONVERGE_get_parameter_id("boundary_injector.tot_inj_mass_tm");
   BOUNDARY_INJECTOR_INJECT_TEMP = CONVERGE_get_parameter_id("boundary_injector.inject_temp");
   BOUNDARY_INJECTOR_Q_RR = CONVERGE_get_parameter_id("boundary_injector.q_rr");
   BOUNDARY_INJECTOR_RATE_SHAPE = CONVERGE_get_parameter_id("boundary_injector.rate_shape");
   BOUNDARY_INJECTOR_SMD = CONVERGE_get_parameter_id("boundary_injector.smd");
   BOUNDARY_INJECTOR_TOT_AREA = CONVERGE_get_parameter_id("boundary_injector.tot_area");
   BOUNDARY_INJECTOR_MASS_PER_PARCEL = CONVERGE_get_parameter_id("boundary_injector.mass_per_parcel");
   BOUNDARY_INJECTOR_MASS_TM1 = CONVERGE_get_parameter_id("boundary_injector.mass_tm1");
   BOUNDARY_INJECTOR_MASS_TM2 = CONVERGE_get_parameter_id("boundary_injector.mass_tm2");
   BOUNDARY_INJECTOR_INJECT_MASS = CONVERGE_get_parameter_id("boundary_injector.inject_mass");
   BOUNDARY_INJECTOR_VELOCITY_OLD_OUT = CONVERGE_get_parameter_id("boundary_injector.velocity_out_old");
   BOUNDARY_INJECTOR_VELOCITY_OLD_NEW = CONVERGE_get_parameter_id("boundary_injector.velocity_out_new");
   BOUNDARY_INJECTOR_PARCEL_RADIUS = CONVERGE_get_parameter_id("boundary_injector.parcel_radius");
   BOUNDARY_INJECTOR_AMP_DISTORT = CONVERGE_get_parameter_id("boundary_injector.amp_distort");
   BOUNDARY_INJECTOR_GAMMA_RR_X_DIST = CONVERGE_get_parameter_id("boundary_injector.gamma_rr_x_dist");
   BOUNDARY_INJECTOR_INJECT_RATIO_VALUE = CONVERGE_get_parameter_id("boundary_injector.inject_ratio_value");
   BOUNDARY_INJECTOR_MFRAC = CONVERGE_get_parameter_id("boundary_injector.mfrac");

   // Get Dynamic IDs for parcel setups
   PARCEL_IN_BREAKUP_FLAG          = CONVERGE_get_parameter_id("parcels_in.breakup_flag");
   PARCEL_IN_KH_FLAG               = CONVERGE_get_parameter_id("parcels_in.kh_flag");
   PARCEL_IN_KH_NEW_PARCEL_FLAG    = CONVERGE_get_parameter_id("parcels_in.kh_new_parcel_flag");
   PARCEL_IN_KH_NO_ENLARGE_FLAG    = CONVERGE_get_parameter_id("parcels_in.kh_no_enlarge_flag");
   PARCEL_IN_KH_NEW_PARCEL_CUTOFF  = CONVERGE_get_parameter_id("parcels_in.kh_new_parcel_cutoff");
   PARCEL_IN_KH_SHED_FACTOR        = CONVERGE_get_parameter_id("parcels_in.kh_shed_factor");
   PARCEL_IN_KH_VEL_CONST          = CONVERGE_get_parameter_id("parcels_in.kh_vel_const");
   PARCEL_IN_KH_SIZE_CONST         = CONVERGE_get_parameter_id("parcels_in.kh_size_const");
   PARCEL_IN_KH_TIME_CONST         = CONVERGE_get_parameter_id("parcels_in.kh_time_const");
   PARCEL_IN_KHACT_TURB_KC         = CONVERGE_get_parameter_id("parcels_in.khact_turb_kc");
   PARCEL_IN_KHACT_TURB_KE         = CONVERGE_get_parameter_id("parcels_in.khact_turb_ke");
   PARCEL_IN_KHACT_TURB_S          = CONVERGE_get_parameter_id("parcels_in.khact_turb_s");
   PARCEL_IN_KHACT_C_TCAV          = CONVERGE_get_parameter_id("parcels_in.khact_c_tcav");
   PARCEL_IN_RT_FLAG               = CONVERGE_get_parameter_id("parcels_in.rt_flag");
   PARCEL_IN_RT_DISTRIBUTION_FLAG  = CONVERGE_get_parameter_id("parcels_in.rt_distribution_flag");
   PARCEL_IN_RT_TIME_CONST         = CONVERGE_get_parameter_id("parcels_in.rt_time_const");
   PARCEL_IN_RT_SIZE_CONST         = CONVERGE_get_parameter_id("parcels_in.rt_size_const");
   PARCEL_IN_RT_LENGTH_CONST       = CONVERGE_get_parameter_id("parcels_in.rt_length_const");
   PARCEL_IN_TAB_FLAG              = CONVERGE_get_parameter_id("parcels_in.tab_flag");
   PARCEL_IN_TAB_DISTRIBUTION_FLAG = CONVERGE_get_parameter_id("parcels_in.tab_distribution_flag");
   PARCEL_IN_TAB_CSUBD             = CONVERGE_get_parameter_id("parcels_in.tab_csubd");
   PARCEL_IN_TAB_CSUBK             = CONVERGE_get_parameter_id("parcels_in.tab_csubk");
   PARCEL_IN_LISA_FLAG             = CONVERGE_get_parameter_id("parcels_in.lisa_flag");
   PARCEL_IN_LISA_DISTRIBUTION_FLAG= CONVERGE_get_parameter_id("parcels_in.lisa_distribution_flag");
   PARCEL_IN_LISA_LENGTH_CONST     = CONVERGE_get_parameter_id("parcels_in.lisa_length_const");
   PARCEL_IN_LISA_SIZE_CONST       = CONVERGE_get_parameter_id("parcels_in.lisa_size_const");

   // Get Dynamic IDs for solid parcel setups
   PARCELS_IN_SOLID_PARCELS_DRAG_MODEL = CONVERGE_get_parameter_id("parcels_in.solid_parcel_drag_flag");
   PARCELS_IN_SOLID_PARCELS_WALL_INTERACTION_FLAG = CONVERGE_get_parameter_id("parcels_in.solid_parcel_wall_interaction_flag");
   PARCELS_IN_SOLID_PARCELS_COEFF_OF_RESTITUTION_TYPE = CONVERGE_get_parameter_id("parcels_in.solid_parcel_coeff_of_restitution_type");
   PARCELS_IN_SOLID_PARCELS_COEFF_OF_RESTITUTION_NORMAL_DIR = CONVERGE_get_parameter_id("parcels_in.solid_parcel_coeff_of_restitution_normal_dir");
   PARCELS_IN_SOLID_PARCELS_COEFF_OF_RESTITUTION_TANG_DIR = CONVERGE_get_parameter_id("parcels_in.solid_parcel_coeff_of_restitution_tang_dir");
   PARCELS_IN_SOLID_PARCELS_STICK_RATE = CONVERGE_get_parameter_id("parcels_in.solid_parcel_stick_rate");
}

/** Load all of the cloud data from a CONVERGE_cloud_t to a wrapper structure
 */
void load_user_cloud(struct ParcelCloud *parcel_cloud_loc, CONVERGE_cloud_t c)
{
   // Load user defined field data
   parcel_cloud_loc->user_temp_starm1 = (double *)CONVERGE_cloud_get_field_data(c, USER_LAG_VAR);
   parcel_cloud_loc->r_bubble = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, R_BUBBLE);
   parcel_cloud_loc->v_bubble = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, V_BUBBLE);
   parcel_cloud_loc->v_drop = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c,V_DROP);
   parcel_cloud_loc->r_bubble_0 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, R_B_0);
   parcel_cloud_loc->r_bubble_tm1 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, R_B_TM1);
   parcel_cloud_loc->v_bubble_tm1 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, V_B_TM1);
   parcel_cloud_loc->r_drop_0 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, R_D_0);
   parcel_cloud_loc->r_therm = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, R_THERM);
   parcel_cloud_loc->omega = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, OMEGA);
   parcel_cloud_loc->omega_tm1 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, OMEGA_TM1);
   parcel_cloud_loc->int_omega = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, INT_OMEGA);
   parcel_cloud_loc->m0 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c,M0);

   parcel_cloud_loc->eta_drop  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, ETA);

   parcel_cloud_loc->user_lag_var_i   = (int *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARi);
   parcel_cloud_loc->tbt              = (int *)CONVERGE_cloud_get_field_data(c,TBT);
   parcel_cloud_loc->is_child = (int *)CONVERGE_cloud_get_field_data(c,IS_CHILD);
   parcel_cloud_loc->pbt              = (int *)CONVERGE_cloud_get_field_data(c,PBT);
   parcel_cloud_loc->child_index = (int *)CONVERGE_cloud_get_field_data(c,CHILD_INDEX);
   parcel_cloud_loc->thermal_breakup_flag = (int *)CONVERGE_cloud_get_field_data(c, THERMAL_BREAKUP_FLAG);
   parcel_cloud_loc->dgre_cycle_count = (int *)CONVERGE_cloud_get_field_data(c, DGRE_COUNT);
   parcel_cloud_loc->parcel_index = (int *)CONVERGE_cloud_get_field_data(c, PARCEL_INDEX);
   parcel_cloud_loc->cloud_index = (int *)CONVERGE_cloud_get_field_data(c, CLOUD_INDEX);
   parcel_cloud_loc->user_lag_var_v3  = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARv3);
   parcel_cloud_loc->user_lag_var_v3b = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARv3b);

   parcel_cloud_loc->user_lag_var_i   = (int *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARi);
   parcel_cloud_loc->child_uu = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, CHILD_UU);
   parcel_cloud_loc->user_lag_var_v3  = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARv3);
   parcel_cloud_loc->user_lag_var_v3b = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARv3b);

   parcel_cloud_loc->from_injector = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FROM_INJECTOR);
   parcel_cloud_loc->from_injector_type = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FROM_INJECTOR_TYPE);
   parcel_cloud_loc->on_triangle   = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_ON_TRIANGLE);
   parcel_cloud_loc->from_nozzle   = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FROM_NOZZLE);
   parcel_cloud_loc->film_flag     = (int *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_FLAG);

   parcel_cloud_loc->just_hit        = (short *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_JUST_HIT);
   parcel_cloud_loc->just_hit_leiden = (short *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_JUST_HIT_LEIDEN);
   parcel_cloud_loc->is_thick        = (char *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_ISTHICK);

   parcel_cloud_loc->rel_vel      = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_REL_VEL);
   parcel_cloud_loc->uprime       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_UPRIME);
   parcel_cloud_loc->uu           = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_UU);
   parcel_cloud_loc->uu_tm1       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_UU_TM1);
   parcel_cloud_loc->xx           = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_XX);
   parcel_cloud_loc->xx_tm1       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_XX_TM1);

   parcel_cloud_loc->v_nu          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_V_NU);
   parcel_cloud_loc->v_sh          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_V_SH);
   parcel_cloud_loc->temp          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TEMP);
   parcel_cloud_loc->temp_tm1      = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TEMP_TM1);
   parcel_cloud_loc->temp_starm1   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TEMP_STARM1);
   parcel_cloud_loc->rey_num       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_REY_NUM);
   parcel_cloud_loc->rel_vel_mag   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_REL_VEL_MAG);
   parcel_cloud_loc->radius        = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_RADIUS);
   parcel_cloud_loc->radius_tm1    = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_RADIUS_TM1);
   parcel_cloud_loc->parent        = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_PARENT);
   parcel_cloud_loc->density       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DENSITY);
   parcel_cloud_loc->density_tm1   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DENSITY_TM1);
   parcel_cloud_loc->gas_density   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_GAS_DENSITY);
   parcel_cloud_loc->mfrac         = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_MFRAC);
   parcel_cloud_loc->mfrac_tm1     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_MFRAC_TM1);
   parcel_cloud_loc->num_drop      = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_NUM_DROP);
   parcel_cloud_loc->surf_temp     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SURF_TEMP);
   parcel_cloud_loc->tbreak_kh     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TBREAK_KH);
   parcel_cloud_loc->shed_num_drop = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SHED_NUM_DROP);
   parcel_cloud_loc->shed_mass     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SHED_MASS);
   parcel_cloud_loc->sactive       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SACTIVE);
   parcel_cloud_loc->sactive_tm1   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SACTIVE_TM1);
   parcel_cloud_loc->surf_ten      = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SURF_TEN);
   parcel_cloud_loc->viscosity     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_VISCOSITY);
   parcel_cloud_loc->distort       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DISTORT);
   parcel_cloud_loc->distort_dot   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DISTORT_DOT);
   parcel_cloud_loc->dm_dt         = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DM_DT);
   parcel_cloud_loc->drdt          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_DRDT);
   parcel_cloud_loc->wall_heat_exchange =
      (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_WALL_HEAT_EXCHANGE);
   parcel_cloud_loc->l_rr          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_L_RR);
   parcel_cloud_loc->l_rc          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_L_RC);
   parcel_cloud_loc->l_temp1       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_L_TEMP1);
   parcel_cloud_loc->l_temp2       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_L_TEMP2);
   parcel_cloud_loc->film_shed     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_SHED);
   parcel_cloud_loc->tbreak_rt     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TBREAK_RT);
   parcel_cloud_loc->surf_temp_tm1 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_SURF_TEMP_TM1);
   parcel_cloud_loc->num_drop_tm1  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_NUM_DROP_TM1);
   parcel_cloud_loc->film_energy   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_ENERGY);
   parcel_cloud_loc->t_turb        = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_T_TURB);
   parcel_cloud_loc->t_turb_accum  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_T_TURB_ACCUM);
   parcel_cloud_loc->film_thickness =
      (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_THICKNESS);
   parcel_cloud_loc->area_reduction =
      (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_AREA_REDUCTION);
   parcel_cloud_loc->tke0     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_TKE0);
   parcel_cloud_loc->eps0     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_EPS0);
   parcel_cloud_loc->lifetime = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_LIFETIME);
   parcel_cloud_loc->film_accum_bit_flag =
      (unsigned long long*)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_ACCUM_BIT_FLAG);
   parcel_cloud_loc->film_accum_plus_bit_flag =
      (unsigned int*)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_ACCUM_PLUS_BIT_FLAG);

   if(CONVERGE_cloud_type(c) == LAGRANGIAN_FILM)
   {
      parcel_cloud_loc->film_thickness_tm1 =
         (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_FILM_THICKNESS_TM1);
      parcel_cloud_loc->area_in_film =
         (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, LAGRANGIAN_AREA_IN_FILM);
   }
   else
   {
      parcel_cloud_loc->film_thickness_tm1 = NULL;
      parcel_cloud_loc->area_in_film       = NULL;
   }
}

void load_user_solid_parcel_cloud(struct ParcelCloud *parcel_cloud_loc, CONVERGE_cloud_t c)
{
   // Load user defined field data
   parcel_cloud_loc->user_temp_starm1 = (double *)CONVERGE_cloud_get_field_data(c, USER_LAG_VAR);
   parcel_cloud_loc->user_lag_var_i   = (int *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARi);
   parcel_cloud_loc->user_lag_var_v3  = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARv3);
   parcel_cloud_loc->user_lag_var_v3b = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, USER_LAG_VARv3b);
   parcel_cloud_loc->child_uu =        (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, CHILD_UU);
   parcel_cloud_loc->child_index = (int *)CONVERGE_cloud_get_field_data(c,CHILD_INDEX);

   parcel_cloud_loc->from_injector = (int *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_FROM_INJECTOR);
   parcel_cloud_loc->from_injector_type = (int *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_FROM_INJECTOR_TYPE);
   parcel_cloud_loc->on_triangle   = (int *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_ON_TRIANGLE);
   parcel_cloud_loc->from_nozzle   = (int *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_FROM_NOZZLE);
   parcel_cloud_loc->film_flag     = (int *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_FILM_FLAG);

   parcel_cloud_loc->just_hit        = (short *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_JUST_HIT);

   parcel_cloud_loc->rel_vel      = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_REL_VEL);
   parcel_cloud_loc->uprime       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_UPRIME);
   parcel_cloud_loc->uu           = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_UU);
   parcel_cloud_loc->uu_tm1       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_UU_TM1);
   parcel_cloud_loc->xx           = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_XX);
   parcel_cloud_loc->xx_tm1       = (CONVERGE_vec3_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_XX_TM1);

   parcel_cloud_loc->v_nu          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_V_NU);
   parcel_cloud_loc->temp          = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_TEMP);
   parcel_cloud_loc->temp_tm1      = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_TEMP_TM1);
   parcel_cloud_loc->temp_starm1   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_TEMP_STARM1);
   parcel_cloud_loc->rey_num       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_REY_NUM);
   parcel_cloud_loc->rel_vel_mag   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_REL_VEL_MAG);
   parcel_cloud_loc->radius        = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_RADIUS);
   parcel_cloud_loc->radius_tm1    = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_RADIUS_TM1);
   parcel_cloud_loc->density       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_DENSITY);
   parcel_cloud_loc->density_tm1   = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_DENSITY_TM1);
   parcel_cloud_loc->mfrac         = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_MFRAC);
   parcel_cloud_loc->mfrac_tm1     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_MFRAC_TM1);
   parcel_cloud_loc->num_drop      = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_NUM_DROP);
   parcel_cloud_loc->surf_temp     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_SURF_TEMP);
   parcel_cloud_loc->viscosity     = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_VISCOSITY);
   parcel_cloud_loc->distort       = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_DISTORT);
   parcel_cloud_loc->surf_temp_tm1 = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_SURF_TEMP_TM1);
   parcel_cloud_loc->num_drop_tm1  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_NUM_DROP_TM1);
   parcel_cloud_loc->t_turb        = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_T_TURB);
   parcel_cloud_loc->t_turb_accum  = (CONVERGE_precision_t *)CONVERGE_cloud_get_field_data(c, SOLID_PARCEL_T_TURB_ACCUM);
}
