
      subroutine phase(ox,oy,oz,px,py,pz,phase_angle,ill_frac) ! in radians
      implicit real*8 (a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
      real*8 ill_frac
      acos(x)=atan2(sqrt(1.-x*x),x)

      dx = ox - px
      dy = oy - py
      dz = oz - pz

      phase_angle = angle_vectors(-px,-py,-pz,dx,dy,dz)
      ill_frac = (1. + cos(phase_angle)) / 2.

      return

      end
