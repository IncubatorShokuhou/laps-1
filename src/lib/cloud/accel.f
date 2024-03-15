

      subroutine accel(iscrpt,nplan,nrel,nmass,npinms,t)
      implicit real*8 (a,b,c,d,e,f,g,h,o,p,q,r,s,t,u,v,w,x,y,z)
      double precision mass,ksq,ksdm15,lma
      dimension mass(13),x(13),y(13),z(13),xdd(13),ydd(13),zdd(13)
     .,ux(13),uy(13),uz(13)
      real*8 mx,my,mz
      common xp(13,59),yp(13,59),zp(13,59),xddp(13,59),yddp(13,59)
     .,zddp(13,59),xd(13),yd(13),zd(13),ksq
      data mass/3.040432924d-6,1.660136795d-7,2.447839598d-6
     .         ,3.227149362d-7,9.547861040d-4,2.858367872d-4
     .         ,4.372731645d-5,5.177591384d-5,3.333333333d-7
     .         ,5.9       d-10,1.3       d-10,1.2       d-10,0.d0/
      data fourm/3.948251502d-8/,cm2/3.335661215d-5/
      data emrat/81.3007d0/,f1/.9878494351d0/,f2/.0121505649/
      data clty,cltz/-.3978811766d0,.9174369566d0/
      data emratp1/82.3007d0/,t1950/2433282.423d0/

      do i=1,nplan
          x(i)=xp(i,iscrpt)
          y(i)=yp(i,iscrpt)
          z(i)=zp(i,iscrpt)
      enddo
c
c calculate mutual accelerations between sun and each planet
      do i=1,nplan
          rlcm1=1.d0/sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
          rlcm3=rlcm1*rlcm1*rlcm1

          if(i.eq.1)then

!             calculate newtonian solar acceleration of earth moon system

!             get earth moon vector for equinox of date
              call moon_brwn(t,mx,my,mz)

!             convert vector to 1950 coordinates
              call preces(t,t1950,mx,my,mz,1)

              px = -mx/emratp1
              py = -my/emratp1
              pz = -mz/emratp1

              xe=x(1)+px
              ye=y(1)+py
              ze=z(1)+pz

              xm=x(1)-px*emrat
              ym=y(1)-py*emrat
              zm=z(1)-pz*emrat

              arge=-ksq/sqrt(xe*xe+ye*ye+ze*ze)**3
              argm=-ksq/sqrt(xm*xm+ym*ym+zm*zm)**3

              xdd(i)=arge*xe*f1+argm*xm*f2
              ydd(i)=arge*ye*f1+argm*ym*f2
              zdd(i)=arge*ze*f1+argm*zm*f2

          else
              arg=-ksq*rlcm3
              xdd(i)=arg*x(i)
              ydd(i)=arg*y(i)
              zdd(i)=arg*z(i)

          endif ! i .eq. 1

30        if(nrel.ge.i)then

!             add perturbative acceleration due to general relativity
              rdrd=x(i)*xd(i)+y(i)*yd(i)+z(i)*zd(i)
              vsq=xd(i)*xd(i)+yd(i)*yd(i)+zd(i)*zd(i)
              a=fourm*rlcm1-vsq*cm2
              b=fourm*rdrd*rlcm3
              xdd(i)=xdd(i)*(1.d0-a)+b*xd(i)
              ydd(i)=ydd(i)*(1.d0-a)+b*yd(i)
              zdd(i)=zdd(i)*(1.d0-a)+b*zd(i)

          endif

          ux(i)=xdd(i)*mass(i)
          uy(i)=ydd(i)*mass(i)
          uz(i)=zdd(i)*mass(i)

      enddo ! n = 1,nplan

c
c add in mutual accelerations between pairs of planets
      do i=1,npinms
          i1=i+1
          do j=i1,nplan
              xdelt=x(j)-x(i)
              ydelt=y(j)-y(i)
              zdelt=z(j)-z(i)
              ksdm15=ksq/sqrt(xdelt*xdelt+ydelt*ydelt+zdelt*zdelt)**3
              arg=ksdm15*mass(j)
              xdd(i)=xdd(i)+arg*xdelt
              ydd(i)=ydd(i)+arg*ydelt
              zdd(i)=zdd(i)+arg*zdelt
              arg=-ksdm15*mass(i)
              xdd(j)=xdd(j)+arg*xdelt
              ydd(j)=ydd(j)+arg*ydelt
              zdd(j)=zdd(j)+arg*zdelt

          enddo

      enddo
c
c calculate total acceleration of sun
      sumux=ux(1)
      sumuy=uy(1)
      sumuz=uz(1)

      do i=2,nmass
          sumux=sumux+ux(i)
          sumuy=sumuy+uy(i)
          sumuz=sumuz+uz(i)

      enddo
c
c compute acc. relative to sun = total acc. - acc. of sun
      do i=1,nplan
          xdd(i)=xdd(i)+sumux
          ydd(i)=ydd(i)+sumuy
          zdd(i)=zdd(i)+sumuz
      enddo

c
c place accelerations in two dimensional array
      do i=1,nplan
          xddp(i,iscrpt)=xdd(i)
          yddp(i,iscrpt)=ydd(i)
          zddp(i,iscrpt)=zdd(i)
      enddo

      return
      end
