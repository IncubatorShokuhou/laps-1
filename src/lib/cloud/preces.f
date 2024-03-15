
      subroutine preces(ti,tf,x,y,z,mode)
      implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)

      atan3(x,y)=dmod((datan2(x,y)+6.2831853071796d0),6.2831853071796d0)

      t0=(ti-2415020.313417d0)/36524.21988d0
      tt=(tf-ti)/36524.21988d0
      zt=(.011171319d0+6.767999d-6*t0)*tt+1.4641d-6*tt*tt+ 8.68d-8*tt**3
      zz=                              zt+3.8349d-6*tt*tt
      th=(.009718973d0-4.135460d-6*t0)*tt-2.0653d-6*tt*tt-20.36d-8*tt**3

      if(mode.eq.1)then
          x0=x
          y0=y
          z0=z

      else ! mode .eq. 2
          x0=dcos(y)*dcos(x)
          y0=dsin(y)*dcos(x)
          z0=dsin(x)

      endif

      sinzt=sin(zt)
c     print*,'x0,y0,z0',x0,y0,z0
      coszt=cos(zt)
      sinth=sin(th)
      costh=cos(th)
      sinzz=sin(zz)
      coszz=cos(zz)
      xx=-sinzt*sinzz+coszt*costh*coszz
      xy=-coszt*sinzz-sinzt*costh*coszz
      xz=-sinth*coszz
      yx= sinzt*coszz+coszt*costh*sinzz
      yy= coszt*coszz-sinzt*costh*sinzz
      yz=-sinth*sinzz
      zx= sinth*coszt
      zy=-sinth*sinzt
      zz= costh
      x   =xx*x0+xy*y0+xz*z0
      y   =yx*x0+yy*y0+yz*z0
      zdum=zx*x0+zy*y0+zz*z0
c     print*,'x,ydum,zdum',x,ydum,zdum

      if(mode.eq.1)then
          z=zdum
c         print*,'preces - x,y,z',x,y,z

      else ! mode .eq. 2
          y=atan3(y,x)
          x=asin(zdum)
c         print*,'preces ra=',y,' dec = ',x

      endif

      return
      end
