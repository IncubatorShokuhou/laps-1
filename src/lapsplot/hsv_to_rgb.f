
      subroutine hsv_to_rgb ( h, s, v, r, g, b )

!*****************************************************************************80
!
!! hsv_to_rgb converts hsv to rgb color coordinates.
!
!  discussion:
!
!    the hsv color system describes a color based on the three qualities
!    of hue, saturation, and value.  a given color will be represented
!    by three numbers, (h,s,v).  h, the value of hue, is an angle 
!    between 0 and 360 degrees, with 0 representing red.  s is the
!    saturation, and is between 0 and 1.  finally, v is the "value",
!    a measure of brightness, which goes from 0 for black, increasing 
!    to a maximum of 1 for the brightest colors.  the hsv color system 
!    is sometimes also called hsb, where the b stands for brightness.
!
!    the rgb color system describes a color based on the amounts of the 
!    base colors red, green, and blue.  thus, a particular color
!    has three coordinates, (r,g,b).  each coordinate must be between
!    0 and 1.  
!
!  licensing:
!
!    this code is distributed under the gnu lgpl license. 
!
!  modified:
!
!    29 august 1998
!
!  author:
!
!    john burkardt
!
!  reference:
!
!    james foley, andries van dam, steven feiner, john hughes,
!    computer graphics, principles and practice,
!    addison wesley, second edition, 1990.
!
!  parameters:
!
!    input, real ( kind = 8 ) h, s, v, the hsv color coordinates to 
!    be converted.
!
!    output, real ( kind = 8 ) r, g, b, the corresponding rgb color coordinates.
!
      implicit none

      real b
      real f
      real g
      real h
      real hue
      integer i
      real p
      real q
      real r
      real s
      real t
      real v

      if ( s == 0.0d+00 ) then

          r = v
          g = v
          b = v

      else
!
!     make sure hue lies between 0 and 360.0d+00
!
          hue = mod ( h, 360.0 )

          hue = hue / 60.0

          i = int ( hue )
          f = hue - real ( i, kind = 8 )
          p = v * ( 1.0 - s )
          q = v * ( 1.0 - s * f )
          t = v * ( 1.0 - s + s * f )

          if ( i == 0 ) then
              r = v
              g = t
              b = p
          else if ( i == 1 ) then
              r = q
              g = v
              b = p
          else if ( i == 2 ) then
              r = p
              g = v
              b = t
          else if ( i == 3 ) then
              r = p
              g = q
              b = v
          else if ( i == 4 ) then
              r = t
              g = p
              b = v
          else if ( i == 5 ) then
              r = v
              g = p
              b = q
          end if

      endif

      return
      end
