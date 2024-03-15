  subroutine ice2vapor(ice, sh, t, p, thresh, ice_m, sh_m, rh_m)

     ! subroutine to convert cloud ice to vapor up to a saturation
     ! threshold (wrt ice)

     implicit none

     ! inputs:

     real, intent(in)    :: ice    ! cloud ice mixing ratio (kg/kg)
     real, intent(in)    :: sh     ! specific humidity (kg/kg)
     real, intent(in)    :: t      ! temperature (k)
     real, intent(in)    :: p      ! pressure (pa)
     real, intent(in)    :: thresh ! saturation factor
     ! set thresh to 1.0 to convert cloud ice up to
     ! ice saturation

     ! outputs:

     real, intent(out)   :: ice_m  ! adjusted lwc
     real, intent(out)   :: sh_m   ! adjusted specific humidity
     real, intent(out)   :: rh_m   ! adjusted rh (%)

     ! locals

     real :: shsat, mr, mrsat, mr_m, mrmax, tc
     real, external :: ssh2, make_rh

     ! set saturation specific humidity for ice for this point
     tc = t - 273.15
     shsat = ssh2(p, tc, tc, 0.)*0.001

     ! convert specific humidity to mixing ratio
     mrsat = shsat/(1.-shsat)
     mr = sh/(1.-sh)
     mrmax = mrsat*thresh

     ! create modified mr (mr_m) by adding cloud ice

     mr_m = mr + ice

     ! zero out the modified cloud ice

     ice_m = 0.

     ! if mr_m exceeds mrmax, convert the excess amount
     ! back to cloud ice

     if (mr_m .gt. mrmax) then
        ice_m = mr_m - mrmax
        mr_m = mrmax
     end if

     rh_m = (mr_m/mrsat)*100.

     ! convert mr_m to sh_m
     sh_m = mr_m/(1.+mr_m)
     return
  end subroutine ice2vapor
