      subroutine rotate_x(rx,ry,rz,o)

      include 'trigd.inc'

      ry1 = ry
      rz1 = rz
      ry=ry1*cosd(o)-rz1*sind(o)
      rz=ry1*sind(o)+rz1*cosd(o)

      return
      end

      subroutine rotate_y(rx,ry,rz,o)

      include 'trigd.inc'

      rz1 = rz
      rx1 = rx
      rz=rz1*cosd(o)-rx1*sind(o)
      rx=rz1*sind(o)+rz1*cosd(o)

      return
      end

      subroutine rotate_z(rx,ry,rz,o)

      include 'trigd.inc'

      rx1 = rx
      ry1 = ry
      rx=rx1*cosd(o)-ry1*sind(o)
      ry=rx1*sind(o)+ry1*cosd(o)

      return
      end
