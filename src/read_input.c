#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <globals.h>
#include <Breakup.h>  // For init_RR_distribution()

// LK global variables are defined in Breakup.c, just set them here

struct UserInputs
{
   double breakup_velocity_scale;
   double breakup_radius_scale;
   double kb_threshold;
   double n_RR;  // Rosin-Rammler shape parameter
   int num_children;  // Number of children per breakup event
   
   // Langmuir-Knudsen evaporation model parameters
   int lk_correction_flag;
   int lk_diagnostic_flag;
   double lk_chi_neq_min;
   double lk_chi_neq_max;
   
   // Spray substepping parameters
   double substep_tau_coeff;
   int substep_max_n;
   
   // Song RPE model selection
   int use_song_rpe;
};

/**********************************************************************************************************/
/* macro: CONVERGE_INPUT                                                                                  */
/* name: user defined                                                                                     */
/* notes: This is called after all other inputs are read in CONVERGE and before the main CFD loop begins. */
/*        It is also called on re-read. It is up to the user to determine if re-read is required for      */
/*        their file.                                                                                     */
/*                                                                                                        */
/* inputs:                                                                                                */
/*          Mesh & Field variables                                                                        */
/*                                                                                                        */
/*          CONVERGE_path_t           passed_path                                                         */
/*           **The path that CONVERGE is currently looking for inputs at.                                 */
/*                                                                                                        */
/*          CONVERGE_flag_t           passed_reread_flag                                                  */
/*             A flag that denotes what type of re-read is occuring.                                      */
/*                0  -  Initial read.                                                                     */
/*                1  -  Re-reading during main loop (re-read)                                             */
/*              **2  -  Re-reading before main loop (override)                                            */
/*                                                                                                        */
/* outputs:                                                                                               */
/*          None                                                                                          */
/*                                                                                                        */
/*  ** This is related to a feature planned for 3.1 and is not required for UDFs implemented in 3.0       */
/**********************************************************************************************************/

CONVERGE_INPUT(read_user,
               IN(VALUE(CONVERGE_mesh_t, mesh),
                  VALUE(CONVERGE_path_t, passed_path),
                  VALUE(CONVERGE_flag_t, passed_reread_flag)))
{
   struct UserInputs* user_inputs;
   if(passed_reread_flag == 1)
   {
      user_inputs = (struct UserInputs*)CONVERGE_mesh_get_user_data(mesh, "user_inputs");
   }
   else
   {
      user_inputs = (struct UserInputs*)CONVERGE_mesh_set_user_data(mesh, "user_inputs", malloc(sizeof(struct UserInputs)));
      user_inputs->breakup_velocity_scale = NAN;
      user_inputs->breakup_radius_scale = NAN;
      user_inputs->kb_threshold = NAN;
      user_inputs->n_RR = 3.2;  // Default RR shape parameter (Kim & Park 2018)
      user_inputs->num_children = 12;  // Default: 12 children (Kim & Park 2018)
      
      // Default LK parameters
      user_inputs->lk_correction_flag = 0;     // Off by default
      user_inputs->lk_diagnostic_flag = 0;     // Off by default
      user_inputs->lk_chi_neq_min = 0.0;
      user_inputs->lk_chi_neq_max = 0.9999;
      
      // Default substepping parameters
      user_inputs->substep_tau_coeff = 0.005;
      user_inputs->substep_max_n = 100;
      
      // Default RPE model selection (set to sentinel value)
      user_inputs->use_song_rpe = 999;         // Sentinel: MUST be read from file
   }

   CONVERGE_logger_concise("reading %*s data from file %s", 16, "user_inputs.in", "user_inputs.in");

   // Open the input file to read
   CONVERGE_file_t file = CONVERGE_file_open_with_path("user_inputs.in", "r", passed_path);
   if(!CONVERGE_file_is_open(file))
   {
      if(passed_reread_flag == 0)
      {
         CONVERGE_logger_concise("Didn't find user_inputs.in, keeping default breakup_velocity_scale = %.3f",
                                 breakup_velocity_scale);
      }
      return;
   }

   int rank;
   CONVERGE_mpi_comm_rank(&rank);

   CONVERGE_file_t echo = CONVERGE_INVALID_HANDLE;
   if(rank == 0)
   {
      // Open the echo file to write
      echo = CONVERGE_file_open_with_path("user_inputs.echo", "w", passed_path);
      if(!CONVERGE_file_is_open(echo))
      {
         CONVERGE_logger_fatal("Couldn't create user_inputs.echo");
         return;
      }
   }

   char buffer[1024];
   CONVERGE_file_register_comment_delimiters(file, "#%");
   while(CONVERGE_file_read_line(file, buffer))
   {
      char* bookmark;
      // Use space, tab, carriage return, and newline as valid delimiters
      char* vtoken = strtok_r(buffer, " \t\r\n", &bookmark);
      char* ktoken = strtok_r(NULL, " \t\r\n", &bookmark);
      if(!vtoken || !ktoken)
      {
         continue;
      }
      // Allow reread/override of breakup velocity scale
      if(strstr(ktoken, "breakup_velocity_scale") || strstr(ktoken, "aa"))
      {
         user_inputs->breakup_velocity_scale = atof(vtoken);
      }
      else if(strstr(ktoken, "breakup_radius_scale") || strstr(ktoken, "B"))
      {
         user_inputs->breakup_radius_scale = atof(vtoken);
      }
      else if(strstr(ktoken, "kb_threshold") || strstr(ktoken, "kb"))
      {
         user_inputs->kb_threshold = atof(vtoken);
      }
      else if(strstr(ktoken, "n_RR") || strstr(ktoken, "n_rosin_rammler"))
      {
         user_inputs->n_RR = atof(vtoken);
      }
      else if(strstr(ktoken, "num_children") || strstr(ktoken, "n_children"))
      {
         user_inputs->num_children = atoi(vtoken);
      }
      // Langmuir-Knudsen parameters
      else if(strstr(ktoken, "lk_correction_flag"))
      {
         user_inputs->lk_correction_flag = atoi(vtoken);
      }
      else if(strstr(ktoken, "lk_diagnostic_flag"))
      {
         user_inputs->lk_diagnostic_flag = atoi(vtoken);
      }
      else if(strstr(ktoken, "lk_chi_neq_min"))
      {
         user_inputs->lk_chi_neq_min = atof(vtoken);
      }
      else if(strstr(ktoken, "lk_chi_neq_max"))
      {
         user_inputs->lk_chi_neq_max = atof(vtoken);
      }
      // Spray substepping parameters
      else if(strstr(ktoken, "substep_tau_coeff"))
      {
         user_inputs->substep_tau_coeff = atof(vtoken);
      }
      else if(strstr(ktoken, "substep_max_n"))
      {
         user_inputs->substep_max_n = atoi(vtoken);
      }
      // Song RPE model selection
      else if(strstr(ktoken, "use_song_rpe"))
      {
         user_inputs->use_song_rpe = atoi(vtoken);
      }
   }

   CONVERGE_bool_t error = CONVERGE_FALSE;
   if(isnan(user_inputs->breakup_velocity_scale))
   {
      CONVERGE_logger_err("user_inputs.in: breakup_velocity_scale missing or invalid");
      error = CONVERGE_TRUE;
   }
   if(isnan(user_inputs->breakup_radius_scale))
   {
      CONVERGE_logger_err("user_inputs.in: breakup_radius_scale missing or invalid");
      error = CONVERGE_TRUE;
   }
   if(isnan(user_inputs->kb_threshold))
   {
      CONVERGE_logger_err("user_inputs.in: kb_threshold missing or invalid");
      error = CONVERGE_TRUE;
   }
   if(user_inputs->num_children < 1)
   {
      CONVERGE_logger_err("user_inputs.in: num_children must be >= 1 (got %d)", user_inputs->num_children);
      error = CONVERGE_TRUE;
   }
   if(user_inputs->num_children > MAX_NUM_CHILDREN)
   {
      CONVERGE_logger_err("user_inputs.in: num_children=%d exceeds MAX_NUM_CHILDREN=%d", 
                          user_inputs->num_children, MAX_NUM_CHILDREN);
      error = CONVERGE_TRUE;
   }

   if(error)
   {
      CONVERGE_logger_fatal("Could not read user_inputs.in. Missing or invalid fields detected.");
   }

   breakup_velocity_scale = (CONVERGE_precision_t)user_inputs->breakup_velocity_scale;
   breakup_radius_scale = (CONVERGE_precision_t)user_inputs->breakup_radius_scale;
   kb_threshold = (CONVERGE_precision_t)user_inputs->kb_threshold;
   num_child_parcels = (CONVERGE_index_t)user_inputs->num_children;
   
   // Set LK global variables
   lk_correction_flag = (CONVERGE_index_t)user_inputs->lk_correction_flag;
   lk_diagnostic_flag = (CONVERGE_index_t)user_inputs->lk_diagnostic_flag;
   lk_chi_neq_min = (CONVERGE_precision_t)user_inputs->lk_chi_neq_min;
   lk_chi_neq_max = (CONVERGE_precision_t)user_inputs->lk_chi_neq_max;
   substep_tau_coeff = (CONVERGE_precision_t)user_inputs->substep_tau_coeff;
   substep_max_n = (int)user_inputs->substep_max_n;

   CONVERGE_logger_concise("========================================");
   CONVERGE_logger_concise("LK Model Parameters (Parsed from user_inputs.in):");
   CONVERGE_logger_concise("  lk_correction_flag = %d", user_inputs->lk_correction_flag);
   CONVERGE_logger_concise("  lk_diagnostic_flag = %d", user_inputs->lk_diagnostic_flag);
   CONVERGE_logger_concise("  lk_chi_neq_min = %.4f", user_inputs->lk_chi_neq_min);
   CONVERGE_logger_concise("  lk_chi_neq_max = %.4f", user_inputs->lk_chi_neq_max);
   CONVERGE_logger_concise("========================================");
   
   // VALIDATE: use_song_rpe MUST be read from file
   if(user_inputs->use_song_rpe == 999)
   {
      CONVERGE_logger_fatal("CRITICAL ERROR: use_song_rpe not found in user_inputs.in!");
      CONVERGE_logger_fatal("Add this line to user_inputs.in:");
      CONVERGE_logger_fatal("    0     use_song_rpe     # 0=thermal RPE, 1=Song isothermal");
   }
   
   // Validate range
   if(user_inputs->use_song_rpe != 0 && user_inputs->use_song_rpe != 1)
   {
      CONVERGE_logger_fatal("ERROR: use_song_rpe must be 0 or 1, got %d", user_inputs->use_song_rpe);
   }
   
   // Set Song RPE model selection
   use_song_rpe = (CONVERGE_index_t)user_inputs->use_song_rpe;
   
   // Log which model was selected
   if(use_song_rpe == 1)
   {
      CONVERGE_logger_concise("========================================");
      CONVERGE_logger_concise("RPE MODEL: Song isothermal (use_song_rpe=1)");
      CONVERGE_logger_concise("========================================");
   }
   else
   {
      CONVERGE_logger_concise("========================================");
      CONVERGE_logger_concise("RPE MODEL: Thermal (use_song_rpe=0)");
      CONVERGE_logger_concise("========================================");
   }
   
   // Initialize Rosin-Rammler distribution parameters
   init_RR_distribution(user_inputs->n_RR);
   
   CONVERGE_logger_verbose("user_inputs->breakup_velocity_scale: %f", user_inputs->breakup_velocity_scale);
   CONVERGE_logger_verbose("user_inputs->breakup_radius_scale: %f", user_inputs->breakup_radius_scale);
   CONVERGE_logger_verbose("user_inputs->kb_threshold: %f", user_inputs->kb_threshold);
   CONVERGE_logger_verbose("user_inputs->n_RR: %f", user_inputs->n_RR);
   CONVERGE_logger_verbose("user_inputs->num_children: %d", user_inputs->num_children);
   CONVERGE_logger_verbose("user_inputs->lk_correction_flag: %d", user_inputs->lk_correction_flag);
   CONVERGE_logger_verbose("user_inputs->lk_diagnostic_flag: %d", user_inputs->lk_diagnostic_flag);
   CONVERGE_logger_verbose("user_inputs->lk_chi_neq_min: %f", user_inputs->lk_chi_neq_min);
   CONVERGE_logger_verbose("user_inputs->lk_chi_neq_max: %f", user_inputs->lk_chi_neq_max);
   CONVERGE_logger_verbose("user_inputs->substep_tau_coeff: %f", user_inputs->substep_tau_coeff);
   CONVERGE_logger_verbose("user_inputs->substep_max_n: %d", user_inputs->substep_max_n);
   CONVERGE_logger_verbose("user_inputs->use_song_rpe: %d", user_inputs->use_song_rpe);

   // Write the echo file
   if(rank == 0)
   {
      CONVERGE_file_write(echo, "%-10.4f breakup_velocity_scale\n", user_inputs->breakup_velocity_scale);
      CONVERGE_file_write(echo, "%-10.4f breakup_radius_scale\n", user_inputs->breakup_radius_scale);
      CONVERGE_file_write(echo, "%-10.4f kb_threshold\n", user_inputs->kb_threshold);
      CONVERGE_file_write(echo, "%-10.4f n_RR\n", user_inputs->n_RR);
      CONVERGE_file_write(echo, "%-10d num_children\n", user_inputs->num_children);
      CONVERGE_file_write(echo, "%-10d lk_correction_flag\n", user_inputs->lk_correction_flag);
      CONVERGE_file_write(echo, "%-10d lk_diagnostic_flag\n", user_inputs->lk_diagnostic_flag);
      CONVERGE_file_write(echo, "%-10.4f lk_chi_neq_min\n", user_inputs->lk_chi_neq_min);
      CONVERGE_file_write(echo, "%-10.4f lk_chi_neq_max\n", user_inputs->lk_chi_neq_max);
      CONVERGE_file_write(echo, "%-10.4f substep_tau_coeff\n", user_inputs->substep_tau_coeff);
      CONVERGE_file_write(echo, "%-10d substep_max_n\n", user_inputs->substep_max_n);
      CONVERGE_file_write(echo, "%-10d use_song_rpe\n", user_inputs->use_song_rpe);
      CONVERGE_file_close(echo);
      CONVERGE_file_destroy(&echo);
   }

   CONVERGE_file_close(file);
   CONVERGE_file_destroy(&file);
}
