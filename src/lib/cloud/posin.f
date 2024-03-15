
      subroutine posin(t,k,rx,ry,rz) ! 1950 coordinates
!     t is tdt (et)
      implicit real*8(a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
      real*8 i,i0,i1,i2,m,m0,m1,m2,m3
      dimension a0(1),a1(1),e0(1),e1(1),e2(1),m0(1),m1(1),m2(1),m3(1)
      dimension i0(1),i1(1),i2(1),ww0(1),ww1(1),ww2(1),w0(1),w1(1),w2(1)
      data a0/1.00000023d0/
      data a1/0.d0/
      data e0/0.0167301085d0/
      data e1/-.41926d-4/
      data e2/-.126d-6/
      data m0/358.000682d0/
      data m1/.9856002628d0/
      data m2/-.155d-3/
      data m3/-.333d-5/
      data i0/0.d0/
      data i1/.013076d0/
      data i2/-.9d-5/
      data ww0/174.40956d0/
      data ww1/-.24166d0/
      data ww2/.6d-4/
      data w0/287.67097d0/
      data w1/.56494d0/
      data w2/.9d-4/
      data pi/3.1415926535897932d0/,o/.4092062118d0/
      rpd=pi/180.
      td=t-2433282.5d0
      tc=td/36525.d0
      tcsq=tc*tc
      a=a0(k)+a1(k)*tc
      e=e0(k)+e1(k)*tc+e2(k)*tcsq
      m=(m0(k)+m1(k)*td+m2(k)*tcsq+m3(k)*tcsq*tc)*rpd
      i=(i0(k)+i1(k)*tc+i2(k)*tcsq)*rpd
      ww=(ww0(k)+ww1(k)*tc+ww2(k)*tcsq)*rpd
      w=(w0(k)+w1(k)*tc+w2(k)*tcsq)*rpd
      ee=m+e*sin(m)
8     arg=ee-(ee-e*dsin(ee)-m)/(1.d0-e*cos(ee))
c     print*,'m,e,ee',m,e,ee
      if(dabs(arg-ee)-1.0d-9)11,11,9
9     ee=arg
      goto8
11    ee=arg
      xw=a*(dcos(ee)-e)
      yw=a*dsqrt(1.d0-e*e)*sin(ee)
c     print*,'xw,yw',xw,yw
      sinw=dsin(w)
      cosw=dcos(w)
      sini=dsin(i)
      cosi=dcos(i)
      sinww=dsin(ww)
      cosww=dcos(ww)
      px= cosw*cosww-sinw*sinww*cosi
      py= cosw*sinww+sinw*cosww*cosi
      pz= sinw*sini
c     print*,'px,py,pz',px,py,pz
      qx=-sinw*cosww-cosw*sinww*cosi
      qy=-sinw*sinww+cosw*cosww*cosi
      qz= cosw*sini
c     print*,'qx,qy,qz',qx,qy,qz
      rx= xw*px+yw*qx
      ry1=xw*py+yw*qy
      rz1=xw*pz+yw*qz
      ry=ry1*cos(o)-rz1*sin(o)
      rz=ry1*sin(o)+rz1*cos(o)
      return
      end
