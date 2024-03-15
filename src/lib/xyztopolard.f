      subroutine xyz_to_polar_d(x,y,z,dec,ra,r)

      include 'trigd.inc'

      r=sqrt(x**2+y**2+z**2)
      dec=asind(z/r)
      ra=atan3d(y,x)

      return
      end
