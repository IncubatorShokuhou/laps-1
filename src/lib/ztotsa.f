
      real function ztotsa(z)

c*  this routine converts a height in meters into a temperature in kelvin

      if (z.lt.11000.) then
          ztotsa=288. - z/161.764              ! ramp to 220.
      else if (z.lt.20000.) then
          ztotsa=220.
      else if (z.lt.50000.) then
          ztotsa=220. + (z-20000.) * .00166666 ! ramp to 270.
      else if (z.lt.90000.) then
          ztotsa=270. - (z-50000.) * .002      ! ramp to 190.
      else
          ztotsa=190. + (z-90000.) * .00566666 ! ramp to 360. (120km)
      end if

      return
      end
