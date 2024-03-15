cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis

      subroutine plot_windob(dir,spd,ri,rj,lat,lon,imax,jmax,relsize
     1                      ,aspect,c4_rot)

      real lat(imax,jmax),lon(imax,jmax)
      character*4 c4_rot

      call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
!     write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
 1234 format(1x,4i5,4e12.4,i4)

!     this keeps about the size of barbs relative to the domain
!     du=(imax)/200. * relsize   

!     this tries to keep the same size of barbs relative to the grid points
      du = relsize

      call get_border(imax,jmax,x_1,x_2,y_1,y_2)
      call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),ltype)

      if(nint(ri).gt.imax.or.nint(rj).gt.jmax)then
          write(6,*)' plot_windob: skipping ob at ',ri,rj
          goto 1
      endif

      if(c4_rot .eq. 'true')then ! convert wind ob from true north to grid north
          rot = projrot_latlon(lat(nint(ri),nint(rj))
     1                        ,lon(nint(ri),nint(rj)),istatus) / 57.295
      else                       ! input is in grid north
          rot = 0.
      endif

!     convert ri and rj to x1 and y1 (u and v)
!     call supcon(alat,alon,x1,y1)
!     x1 = umin + (umax - umin) * (ri-1.) / float(imax-1)
!     y1 = vmin + (vmax - vmin) * (rj-1.) / float(jmax-1)

      if(dir .gt. -400.)then
          if(lat(nint(ri),nint(rj)) .ge. 0.)then
              sense = +1.0
          else
              sense = -1.0
          endif

          call barbs(spd,dir,ri,rj,du,rot,-1e10,+1e10,-1e10,+1e10
     1                                               ,sense,aspect)
      endif

    1 continue
      return
      end


      subroutine plot_vr(i,j,vel,imax,jmax,c1_plottype,n_plotted)

      character*1 c1_plottype

      save init
      data init/0/

!     call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
!     write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
!1234 format(1x,4i5,4e12.4,i4)
!     du=(umax-umin)/200.

!     call get_border(imax,jmax,x_1,x_2,y_1,y_2)
!     call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax))

!     du = (umax-umin)/imax
!     dv = (vmax-vmin)/jmax

      du = 1. / float(imax-1)
      dv = 1. / float(jmax-1)

      u = float(i-1) / float(imax-1)
      v = float(j-1) / float(jmax-1)

      if(.true.)then

      call get_border(imax,jmax,x_1,x_2,y_1,y_2)

      if(init .eq. 0)then ! avoid multiple calls for efficiency
          call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),1)
          call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
          write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
1234      format(1x,4i5,4e12.4,i4)
          init = 1
      endif

      u = i
      v = j
      du = 1.
      dv = 1.

      endif

      if(c1_plottype .eq. 'y')then

          icol = min(max(  110. + vel / 2.0   ,  100.  )  ,120.)
          write(6,*)i,j,vel,icol
          call setusv_dum(2hin,icol)
!         call setusv_dum(2hin,40)

          if(n_plotted .eq. 1)then
              do uu = u-du/2.,u+du/2.,du/25.
                  call line(uu,v-dv/2.,uu,v+dv/2.)
              enddo

          else ! 2 or more
              do uu = u-du/2.,u,du/25.
                  call line(uu,v-dv/2.,uu,v+dv/2.)
              enddo

          endif

      elseif(c1_plottype .eq. 'a')then ! plot arrows


      endif

      return
      end


      subroutine map(iproj,selat,selon,nwlat,nwlon)
      real nwlat,nwlon
      cenlat=90. ! (selat+nwlat)/2.
      cenlon=-105. ! (selon+nwlon)/2.
      rot=0.
      igrid=0
      idetail=4
      idot_pat=0
      call supmap(iproj,cenlat,cenlon,0.0,selat,selon,nwlat,nwlon,-2,
     1igrid,idetail,idot_pat,ier)
!     call savesup_gp('laps_sup.parms',istatus)
      return
      end


c
        subroutine barbs(spd,dir,u,v,du,projrot,umin,umax,vmin,vmax
     1                                              ,sense,aspect)
c
c---wind barb plotter from barbs_gp.
c
c       include 'sysdisk:[gudoc]edfvaxbox.for'
c       paul schultz    30-jun-1982     original version
c       hoagenson       02-sep-1982     changed name from barbs
c       schultz         25-oct-1984     changed zero-wind symbol
c
c---dir,spd     wind direction and speed
c---u,v         ncar coordinates of obs
c---du          arbitrarily chosen increment of screen space
c---projrot     rotation of barb due to map projection
c---aspect      ratio of du scaling divided by dv scaling (stretch in u dir)
c
        data degrad/.01745329/
c
c
c---return if missing speed or direction.
c
        if (dir .lt. 0. .or. spd .lt. 0.) return
        if (dir .gt. 360. .or. spd .gt. 300.) return

        thk_base = 1.0
c
c---directions:
c
        dr=dir*degrad+projrot
        dr1=(dir+(60.*sense))*degrad+projrot
        dr90=dr + 90.*degrad
        sind=sin(dr)
        sind1=sin(dr1)
        cosd=cos(dr)
        cosd1=cos(dr1)

        thk =sqrt(  sind**2          + (cosd *aspect)**2 ) ! staff
        call thk_transform(thk_base,aspect,dr,thk)

        thk1=sqrt(  sind1**2         + (cosd1*aspect)**2 ) ! barbs
        call thk_transform(thk_base,aspect,dr1,thk1)

        thkf=sqrt( (sind*aspect)**2  + cosd**2)            ! flag
        call thk_transform(thk_base,aspect,dr90,thkf)

!       write(6,11)aspect,dr*180./3.14159265,dr1*180./3.14159265
!    1            ,thk,thk1,thkf
!11     format(' aspect,dr,dr1,thk,thk1,thkf ',f6.3,2f8.1,3f9.3)
c
c---lengths:
c
        staff=du*4.                     !  multiplier arbitrarily chosen
        barb=staff*.5
        add=staff*.3
c
c---speed and counters:
c
        n50=0
        n10=0
c
!       if(u .le. umin .or. u .ge. umax .or. v .le. vmin .or. v .ge. vmax)
!       1                                                       return

        if (spd .lt. 1.0) then ! calm winds
!           call pwrit (u,v,'o',1,6,0,0)
!           call pcloqu(u,v,'o',du,angd,cntr)
            call plot_circle(u,v,du*0.8)
        return

        end if
        if (spd .lt. 2.5) then
        x1=u
        y1=v
        x2=x1+sind*staff*aspect
        y2=y1+cosd*staff

        if(x2 .lt. umin .or. x2 .gt. umax .or. y2 .lt. vmin .or. y2 .gt.
     1 vmax)then
            call gslwsc(thk_base)
            return
        endif

        call gslwsc(thk)
        call line(x1,y1,x2,y2)

        call gslwsc(thk_base)
        return

        endif
c
        sp=spd+2.5
c
        do while (sp .ge. 50.)
        n50=n50+1
        sp=sp-50.
        end do
c
        do while (sp .ge. 10.)
        n10=n10+1
        sp=sp-10.
        end do
c
c---draw staff
c
        x1=u
        y1=v
        x2=x1+sind*staff*aspect
        y2=y1+cosd*staff

        if(x2 .lt. umin .or. x2 .gt. umax .or. y2 .lt. vmin .or. y2 .gt.
     1 vmax)then
            return
        endif

        call gslwsc(thk)
        call line(x1,y1,x2,y2)
        call gslwsc(thk_base)
c
c---plot half-barb, if necessary
c
        if (sp .ge. 5.) then
        x1=x2+sind1*add*aspect
        y1=y2+cosd1*add

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        call gslwsc(thk1)
        call line(x1,y1,x2,y2)

        if (n50 .ne. 0 .or. n10 .ne. 0) go to 40
        x1=x2+sind*add*aspect
        y1=y2+cosd*add

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        call gslwsc(thk)
        call line(x1,y1,x2,y2)

        call gslwsc(thk_base)
        return

        end if
c
40      x1=x2
        y1=y2
c
c---plot barbs, if necessary
c
        do 50 i=1,n10
        x2=x1+sind*add*aspect
        y2=y1+cosd*add
        x3=x2+sind1*barb*aspect
        y3=y2+cosd1*barb

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

!       call frstpt(x1,y1)

!       if(x2 .le. umin .or. x2 .ge. umax .or. y2 .le. vmin .or. y2 .ge. vmax)
!       1                                                       return

!       call vector(x2,y2)

        call gslwsc(thk)
        call line(x1,y1,x2,y2)

!       if(x3 .le. umin .or. x3 .ge. umax .or. y3 .le. vmin .or. y3 .ge. vmax)
!       1                                                       return

        call gslwsc(thk1)
!       call vector(x3,y3)
        call line(x2,y2,x3,y3)
        x1=x2
        y1=y2
50      continue
c
c---plot flags, if necessary
c
        do 60 i=1,n50
        x2=x1+sind*add*aspect
        y2=y1+cosd*add
        x3=x2+sind1*barb*aspect
        y3=y2+cosd1*barb

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        call gslwsc(thk)
!       call frstpt(x1,y1)

!       if(x2 .le. umin .or. x2 .ge. umax .or. y2 .le. vmin .or. y2 .ge. vmax)
!       1                                                       return

!       call vector(x2,y2)
        call line(x1,y1,x2,y2)

!       if(x3 .le. umin .or. x3 .ge. umax .or. y3 .le. vmin .or. y3 .ge. vmax)
!       1                                                       return

        call gslwsc(thk1)
!       call vector(x3,y3)
        call line(x2,y2,x3,y3)

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        call gslwsc(thkf)
!       call vector(x1,y1)
        call line(x3,y3,x1,y1)
        x1=x2
        y1=y2
60      continue
c
        call gslwsc(thk_base)
        return
        end
c

        subroutine plot_circle(u,v,size)

        data degrad/.01745329/

        do iaz = 0,360,20
            az = float(iaz) * degrad
            x = u + size * sin(az)
            y = v + size * cos(az)
            if(iaz .eq. 0)then
                call frstpt(x,y)
            else
                call vector(x,y)
            endif
        enddo
 
        return
        end

        subroutine plot_circle_fill(u,v,radius,frac)

        data degrad/.01745329/

        if(frac .lt. 0.10)return
        
        if(frac .ge. 0.10 .and. frac .le. 0.50)then
!           add straight line
            call line(u,v-radius,u,v+radius)
        endif

        if(frac .gt. .50 .and. frac .lt. 0.9375)then ! bkn
!           overcast
            do iaz = -20,+20,40
                az = float(iaz) * degrad
                x1 = u + radius * sin(az)
                y1 = v - radius * cos(az)
                x2 = u + radius * sin(az)
                y2 = v + radius * cos(az)
                call line(x1,y1,x2,y2)
            enddo
        endif


        if(frac .ge. 0.9375)then ! overcast - fill full circle
            do iaz = 0,270,90
                az = float(iaz) * degrad
                x = u + radius * sin(az)
                y = v + radius * cos(az)
                call line(u,v,x,y)
            enddo
        endif

        return
        end

        subroutine thk_transform(t1,aspect,d1,t2)

        real d1     ! direction in radians
        real t1     ! reference thickness
        real aspect ! aspect ratio of projection stretching
        real t2     ! new thickness

        if(sin(d1) .eq. 0.)then
            t2 = t1 * aspect
            return
        endif

        if(cos(d1) .eq. 0)then
            t2 = t1
            return
        endif

        a1 = abs(t1/sin(d1))

        b1 = abs(t1/cos(d1))

        d2 = atan(aspect*tan(d1))

!       a2 = a1
!       b2 = b1 * aspect
!       t2 = a2 * sin(d2)

        t2 = abs(a1*sin(d2))
   
        return
        end
