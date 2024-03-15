        subroutine sidereal_time(ut1,lon,lst)

      implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
      real*8 lst,lon

      include '../../include/astparms.for'

!       ut1 is in days (jd)
!       lon is in degrees (west is negative)
!       lst is in radians

        parameter (rad_per_sec = 2d0 * pi / 86400d0)

      tu=(ut1-2451545.d0)/36525.d0

        gmst_0ut = 24110.54841d0 +
     1     tu * (8640184.812866d0  + tu * (.093104d0 - tu * 6.2d-6))

        gmst_0ut = gmst_0ut * rad_per_sec

        lst = gmst_0ut + (ut1+.5d0) * 2d0 * pi + lon * rpd

        lst = mod(lst,2d0*pi)

        return
        end
