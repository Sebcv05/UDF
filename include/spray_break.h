//Variables to reduce preamble in spray_break.c

// Local UDF Field Variable Arrays
// Declared here to simplify function signature for spray_break
static const CONVERGE_precision_t *global_density;
static const CONVERGE_precision_t *global_mol_cond;
static const CONVERGE_precision_t *global_csubp;
static const CONVERGE_precision_t *global_csubv;
static const CONVERGE_precision_t *global_mol_viscosity;
static const CONVERGE_precision_t *global_temperature;
static const CONVERGE_precision_t *global_volume;
static const CONVERGE_precision_t *global_pressure;
static const CONVERGE_precision_t *global_sensible_sie;
static const CONVERGE_precision_t *global_species_massfrac;
static const CONVERGE_precision_t *global_species_massfrac_tm1;
static const CONVERGE_precision_t *global_density_tm1;
static const CONVERGE_precision_t *global_mol_viscosity;
static CONVERGE_precision_t *global_src_ex_density;
static CONVERGE_precision_t *global_src_unburned_enth_evap_ex;
static CONVERGE_precision_t *global_src_sie_ex;
static CONVERGE_precision_t *global_src_species_ex;
static CONVERGE_precision_t *global_src_rif_fuel_evap_ex;
static CONVERGE_precision_t *global_ecfm_temp_evap_ex;
static const short *global_region_index;


static CONVERGE_species_t sp;
static CONVERGE_index_t num_gas_species;
static CONVERGE_index_t num_parcel_species;
static CONVERGE_index_t num_total_species;

   struct vect
   {
      CONVERGE_precision_t x;
      CONVERGE_precision_t y;
      CONVERGE_precision_t z;
   };


