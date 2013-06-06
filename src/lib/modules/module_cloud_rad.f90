
MODULE cloud_rad    

!     Cloud Radiation and Microphysics Parameters
      real, parameter :: bksct_eff_clwc = .07
      real, parameter :: bksct_eff_cice = .14
      real, parameter :: bksct_eff_rain = .07
      real, parameter :: bksct_eff_snow = .14

      real, parameter :: rholiq =    1e3 ! kilograms per cubic meter
      real, parameter :: rhosnow = .07e3 ! kilograms per cubic meter

      real, parameter :: reff_clwc = .000020 ! m
      real, parameter :: reff_cice = .000040 ! m
      real, parameter :: reff_rain = .001000 ! m
      real, parameter :: reff_snow = .004000 ! m

      PUBLIC albedo_to_clouds

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine albedo_to_clouds(albedo                                  & ! I
                                ,cloud_rad_trans,cloud_od,cloud_opacity)   ! O

        use mem_namelist, ONLY: r_missing_data

        clear_albedo = .2097063
        cloud_albedo = .4485300
!       cloud_albedo = .40

        if(albedo .eq. r_missing_data)then
            cloud_rad_trans = r_missing_data
            cloud_od = r_missing_data
            cloud_opacity = r_missing_data
            return
        endif

!       Uniform clouds are presently assumed in the grid box or cloud layer
        cloud_frac = 1.0

        arg = albedo

        call stretch2(clear_albedo,cloud_albedo,0.,1.,arg)

        cloud_albedo = min(arg,0.99999)

        cloud_rad_trans = 1.0 - cloud_albedo

        bksc_eff_od = -log(cloud_rad_trans)

        cloud_od = bksc_eff_od * 10.

        cloud_opacity = 1.0 - exp(cloud_od)

        return
   
   END subroutine albedo_to_clouds    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE cloud_rad   
