  subroutine saturate_ice_points(t, p, thresh, sh_m, rh_m)

     ! subroutine to saturate grid boxes with respect to ice

     implicit none

     ! inputs:

     real, intent(in)    :: t      ! temperature (k)
     real, intent(in)    :: p      ! pressure (pa)
     real, intent(in)    :: thresh ! saturation factor
     ! set thresh to 1.0 to convert cloud ice up to
     ! ice saturation

     ! outputs:

     real, intent(out)   :: sh_m   ! adjusted specific humidity
     real, intent(out)   :: rh_m   ! adjusted rh (%)

     ! locals

     real :: shsat, mr, mrsat, mrmax, tc, esi, esw, e
     real, external :: ssh2, es, esice

     tc = t - 273.15

     ! determine ice satuaration vapor pressure
     esi = esice(tc)

     ! determine water saturation vapor pressure
     esw = es(tc)

     ! compute the e needed to reach thresh saturation wrt ice
     e = thresh*esi

     ! compute rh wrt liquid for this new e
     rh_m = e/esw*100.
     ! compute saturated sh by seting e = esi in the
     ! typical formula for q

     sh_m = (0.622*esi)/(p - 0.378*esi)

     return
  end subroutine saturate_ice_points
