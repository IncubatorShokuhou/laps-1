
      subroutine xyz_to_polar_r(x,y,z,dec,ra,r)

      implicit real*8(a-z)

      atan3(x,y)=dmod((datan2(x,y)+6.2831853071796d0),6.2831853071796d0)

      r=dsqrt(x**2+y**2+z**2)
      dec=asin(z/r)
      ra=atan3(y,x)

      return
      end
