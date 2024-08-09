#include <CONVERGE/udf.h>
#include <samples.h>

CONVERGE(configure)
{


    CONVERGE_register_onload("spray_env",&spray_env);
    // CONVERGE_register_onload("cvg_cloud_properties",&cvg_cloud_properties)
    CONVERGE_register_onload("custom_spray_evap_onload",&custom_spray_evap_onload);
    CONVERGE_register_spray_evap("custom_spray_evap",&custom_spray_evap);
    CONVERGE_register_onclose("custom_spray_evap_onclose:",&custom_spray_evap_onclose);


    CONVERGE_register_parcel_inject("custom_parcel_inject", &custom_parcel_inject);
    CONVERGE_register_parcel_child( "custom_parcel_child",  &custom_parcel_child);
    CONVERGE_register_parcel_splash( "custom_parcel_splash", &custom_parcel_splash);
    CONVERGE_register_parcel_strip( "custom_parcel_strip",  &custom_parcel_strip);
    CONVERGE_register_onclose("spray_env",&spray_env);
    //    CONVERGE_register_onclose("release_shape_onclose", &release_shape_onclose);

}



