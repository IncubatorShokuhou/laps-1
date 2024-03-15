  subroutine saturate_lwc_points(sh, t, p, thresh, sh_m, rh_m)

     ! subroutine to saturate points containing a certain threshold
     ! of cloud liquid

     implicit none

     ! inputs:

     real, intent(in)    :: sh     ! specific humidity (kg/kg)
     real, intent(in)    :: t      ! temperature (k)
     real, intent(in)    :: p      ! pressure (pa)
     real, intent(in)    :: thresh ! saturation factor
     ! set thresh to 1.0 to convert cloud water up to
     ! vapor saturation.  1.1 will allow 110% rh, and so
     ! forth

     ! outputs:

     real, intent(out)   :: sh_m   ! adjusted specific humidity
     real, intent(out)   :: rh_m   ! adjusted rh (%)

     ! locals

     real :: shsat, mrmax, mr, mr_m, mrsat
     real, external :: ssh, make_rh

     ! set saturation specific humidity for this point

     shsat = ssh(p, t - 273.15)*0.001

     ! convert specific humidity to mixing ratio
     mrsat = shsat/(1.-shsat)
     mr = sh/(1.-sh)

     mrmax = mrsat*thresh

     rh_m = thresh*100.

     ! convert modified mixing ratio back to specific humidity
     sh_m = mrmax/(1.+mrmax)

     return
  end subroutine saturate_lwc_points
