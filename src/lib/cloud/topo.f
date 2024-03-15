
      subroutine topo(phi,lon,ut1,tx,ty,tz)
      implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
      real*8 lst,lon,jd

      include '../../include/astparms.for'

      parameter (r_e_au = r_e_km / km_per_au)

!       all arguments are real * 8
!       input phi = latitude in radians
!             lon = longitude in degrees (w is negative)
!             ut1 = julian date
!       output tx,ty,tz = equatorial coordinates of point on surface of earth

c     write(6,*)r_e_au

      ff=(1.d0-1.d0/298.257d0)**2
      cc=1./dsqrt(dcos(phi)**2+ff*dsin(phi)**2)
      xyg=r_e_au*cc*dcos(phi)

      call sidereal_time(ut1,lon,lst)

      tx=dcos(lst)*xyg
      ty=dsin(lst)*xyg
      tz=r_e_au*cc*dsin(phi)*ff

      return
      end

