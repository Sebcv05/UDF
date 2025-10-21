#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <globals.h>

struct UserInputs
{
   double breakup_velocity_scale;
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

   char buffer[100];
   CONVERGE_file_register_comment_delimiters(file, "#%");
   while(CONVERGE_file_read_line(file, buffer))
   {
      char* bookmark;
      char* vtoken = strtok_r(buffer, " ", &bookmark);
      char* ktoken = strtok_r(bookmark, " ", &bookmark);
      if(!vtoken || !ktoken)
      {
         continue;
      }
      // Allow reread/override of breakup velocity scale
      if(strcmp(ktoken, "breakup_velocity_scale") == 0 || strcmp(ktoken, "aa") == 0)
      {
         user_inputs->breakup_velocity_scale = atof(vtoken);
      }
   }

   CONVERGE_bool_t error = CONVERGE_FALSE;
   if(isnan(user_inputs->breakup_velocity_scale))
   {
      CONVERGE_logger_err("user_inputs.in: breakup_velocity_scale missing or invalid");
      error = CONVERGE_TRUE;
   }

   if(error)
   {
      CONVERGE_logger_fatal("Could not read user_inputs.in. Missing or invalid fields detected.");
   }

   breakup_velocity_scale = (CONVERGE_precision_t)user_inputs->breakup_velocity_scale;
   CONVERGE_logger_verbose("user_inputs->breakup_velocity_scale: %f", user_inputs->breakup_velocity_scale);

   // Write the echo file
   if(rank == 0)
   {
      CONVERGE_file_write(echo, "%-10.4f breakup_velocity_scale\n", user_inputs->breakup_velocity_scale);
      CONVERGE_file_close(echo);
      CONVERGE_file_destroy(&echo);
   }

   CONVERGE_file_close(file);
   CONVERGE_file_destroy(&file);
}
