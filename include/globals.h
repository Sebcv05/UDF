//Global variables 

#ifndef GLOBALS_H
#define GLOBALS_H

// Global velocity variables
extern CONVERGE_precision_t user_child_velocity_x;
extern CONVERGE_precision_t user_child_velocity_y;
extern CONVERGE_precision_t user_child_velocity_z ;

// Breakup tuning parameter (default set in Breakup.c, overwritten by user input)
extern CONVERGE_precision_t breakup_velocity_scale;
extern CONVERGE_precision_t breakup_radius_scale;
extern CONVERGE_precision_t kb_threshold;

#endif
