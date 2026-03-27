import sys
import re
import math

with open('/home/apollo19/Desktop/Dan_B/UDF/src/spray_evap.c', 'r') as f:
    text = f.read()

idx = text.find('for(CONVERGE_index_t i_pc = CONVERGE_iterator_first(pc_it);')
if idx == -1:
    print("Cannot find parcel loop.")
    sys.exit(1)

pattern_setup = r'(for\(CONVERGE_index_t i_pc = CONVERGE_iterator_first\(pc_it\);\s*i_pc != -1;\s*i_pc = CONVERGE_iterator_next\(pc_it\)\)\n\s*\{)\n(\s*int inner_iter_flag)'
replacement_setup = r'''\1
         int n_sub = 1;
         CONVERGE_precision_t drdt_guess = 0.0;
         for (int isp = 0; isp < num_parcel_species; isp++) {
             drdt_guess += parcel_cloud.drdt[i_pc * num_parcel_species + isp];
         }
         if (fabs(drdt_guess) > 1.0e-12) {
             CONVERGE_precision_t tau_R = parcel_cloud.radius_tm1[i_pc] / fabs(drdt_guess);
             n_sub = (int)ceil(dt / (0.5 * tau_R));
         }
         if (n_sub < 1) n_sub = 1;
         if (n_sub > 10) n_sub = 10;
         
         CONVERGE_precision_t dt_sub = dt / n_sub;

         // Sub-step tm1 local variables
         CONVERGE_precision_t sub_temp_tm1 = parcel_cloud.temp_tm1[i_pc];
         CONVERGE_precision_t sub_radius_tm1 = parcel_cloud.radius_tm1[i_pc];
         CONVERGE_precision_t sub_density_tm1 = parcel_cloud.density_tm1[i_pc];
         CONVERGE_precision_t sub_mfrac_tm1[100]; // assuming max 100 species
         for(int isp=0; isp<num_parcel_species; ++isp) {
             sub_mfrac_tm1[isp] = parcel_cloud.mfrac_tm1[i_pc * num_parcel_species + isp];
         }
         CONVERGE_precision_t sub_mass_tm1 = sub_density_tm1 * (4.0 / 3.0) * M_PI * (sub_radius_tm1 * sub_radius_tm1 * sub_radius_tm1);
         
         // Accumulators for cell source sums over dt
         CONVERGE_precision_t acc_dm_dt[100];
         for(int isp=0; isp<num_parcel_species; ++isp) acc_dm_dt[isp] = 0.0;

         // Start sub-loop
         for (int sub_iter = 0; sub_iter < n_sub; ++sub_iter)
         {
\2'''

text1 = re.sub(pattern_setup, replacement_setup, text, count=1)

text2 = re.sub(
    r'(CONVERGE_precision_t tdrop\s*=\s*)parcel_cloud\.temp\[i_pc\];',
    r'\1 (sub_iter == 0) ? parcel_cloud.temp[i_pc] : sub_temp_tm1;',
    text1, count=1
)
text2 = re.sub(
    r'(CONVERGE_precision_t tdrop_starm1\s*=\s*)parcel_cloud\.temp\[i_pc\];',
    r'\1 (sub_iter == 0) ? parcel_cloud.temp[i_pc] : sub_temp_tm1;',
    text2, count=1
)
text2 = re.sub(
    r'(CONVERGE_precision_t temp_prev_timestep\s*=\s*)parcel_cloud\.temp_tm1\[i_pc\];',
    r'\1 sub_temp_tm1;',
    text2, count=1
)

text2 = re.sub(
    r'parcel_cloud\.temp_tm1\[i_pc\]\s*=\s*temp_prev_timestep;',
    r'sub_temp_tm1 = temp_prev_timestep;',
    text2, count=1
)

text3 = re.sub(
    r'mass_drop_tm1\s*=\s*parcel_cloud\.density_tm1\[i_pc\] \*\s*\(4\.0\s*/\s*3\.0\)\s*\*\s*PI\s*\*\s*\(CONVERGE_cube\(parcel_cloud\.radius_tm1\[i_pc\]\)\);',
    r'mass_drop_tm1 = sub_mass_tm1;',
    text2, count=1
)

text4 = text3.replace('parcel_cloud.mfrac_tm1[i_pc * num_parcel_species + isp]', 'sub_mfrac_tm1[isp]')

text6 = re.sub(r'csubp_liquid \*\s*mass_drop_new\s*\*\s*parcel_cloud\.temp\[i_pc\]', r'csubp_liquid * mass_drop_new * sub_temp_tm1', text4)

text7 = re.sub(r'parcel_cloud\.dm_dt\[(.*?)\]\s*\*\s*dt', r'parcel_cloud.dm_dt[\1] * dt_sub', text6)
text7 = re.sub(r'\(parcel_cloud\.radius\[i_pc\] - evap_min_radius\[isp\]\) / dt', r'(parcel_cloud.radius[i_pc] - evap_min_radius[isp]) / dt_sub', text7)
text7 = re.sub(r'evap_mass_drop_0\[isp\] / dt', r'evap_mass_drop_0[isp] / dt_sub', text7)
text7 = re.sub(r'mass_drop_tm1\) / dt', r'mass_drop_tm1) / dt_sub', text7)
text7 = re.sub(r'dt \* drop_area\b', r'dt_sub * drop_area', text7)
text7 = re.sub(r'dt \* drop_area_1\b', r'dt_sub * drop_area_1', text7)
text7 = re.sub(r'parcel_cloud\.drdt\[(.*?)\]\s*\*\s*dt', r'parcel_cloud.drdt[\1] * dt_sub', text7)


idx2 = text7.find('parcel_cloud.temp_starm1[i_pc] = tdrop;')
if idx2 != -1:
    close_subloop = r'''
            // Accumulate dm_dt
            for(int isp=0; isp<num_parcel_species; isp++) {
                acc_dm_dt[isp] += parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] * dt_sub;
            }

            // Update sub_tm1 variables for next sub-step
            sub_temp_tm1 = tdrop;
            sub_radius_tm1 = radius_new[i_pc];
            sub_density_tm1 = parcel_cloud.density[i_pc];
            sub_mass_tm1 = sub_density_tm1 * (4.0 / 3.0) * M_PI * (sub_radius_tm1 * sub_radius_tm1 * sub_radius_tm1);
            for(int isp=0; isp<num_parcel_species; ++isp) {
                sub_mfrac_tm1[isp] = parcel_cloud.mfrac[i_pc * num_parcel_species + isp];
            }
         } // end of sub_iter loop

         for(int isp=0; isp<num_parcel_species; isp++) {
             parcel_cloud.dm_dt[i_pc * num_parcel_species + isp] = acc_dm_dt[isp] / dt;
         }

         parcel_cloud.temp_starm1[i_pc] = tdrop;
'''
    text8 = text7[:idx2] + close_subloop + text7[idx2 + len('parcel_cloud.temp_starm1[i_pc] = tdrop;'):]

with open('/home/apollo19/Desktop/Dan_B/UDF/spray_evap_mod.c', 'w') as f:
    f.write(text8)

print("Modification complete.")
