  subroutine lwc2vapor(lwc, sh, t, p, thresh, lwc_m, sh_m, rh_m)

     ! subroutine to convert cloud water to vapor.

     implicit none

     ! inputs:

     real, intent(in)    :: lwc    ! cloud water mixing ratio (kg/kg)
     real, intent(in)    :: sh     ! specific humidity (kg/kg)
     real, intent(in)    :: t      ! temperature (k)
     real, intent(in)    :: p      ! pressure (pa)
     real, intent(in)    :: thresh ! saturation factor
     ! set thresh to 1.0 to convert cloud water up to
     ! vapor saturation.  1.1 will allow 110% rh, and so
     ! forth

     ! outputs:

     real, intent(out)   :: lwc_m  ! adjusted lwc
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

     ! create modified mixing ratio (mr_m) by adding cloud liquid

     mr_m = mr + lwc

     ! zero out the modified cloud water

     lwc_m = 0.

     ! if mr_m exceeds mrmax, convert the excess amount
     ! back to cloud water

     if (mr_m .gt. mrmax) then
        lwc_m = mr_m - mrmax
        mr_m = mrmax
     end if

     ! compute rh from modified mixing ratio
     rh_m = (mr_m/mrsat)*100.

     ! convert modified mixing ratio back to specific humidity
     sh_m = mr_m/(1.+mr_m)

     return
  end subroutine lwc2vapor
