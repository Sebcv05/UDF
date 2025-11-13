//Global variables 

#ifndef GLOBALS_H
#define GLOBALS_H

#define MAX_NUM_CHILDREN 50  // Maximum allowed children per breakup event

// Global velocity variables
extern CONVERGE_precision_t user_child_velocity_x;
extern CONVERGE_precision_t user_child_velocity_y;
extern CONVERGE_precision_t user_child_velocity_z ;

// Breakup tuning parameter (default set in Breakup.c, overwritten by user input)
extern CONVERGE_precision_t breakup_velocity_scale;
extern CONVERGE_precision_t breakup_radius_scale;
extern CONVERGE_precision_t kb_threshold;
extern CONVERGE_index_t num_child_parcels;

#endif
