cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
!       i4time = i4time_now_gg()
!       rlat = 40.
!       rlon = -105.
!       call solar_position(rlat,rlon,i4time,alt,dec,hrangle)
!       write(6,*)'hrangle,dec,alt',hrangle,dec,alt
!       end

        subroutine solar_position(rlat,rlon,i4time,alt,dec,hrangle)

c steve albers jan 1994
c argument      i/o     type                    description
c --------      ---     ----    -----------------------------------------------
c rlat           i      r*4     latitude (degrees)
c rlon           i      r*4     longitude (degrees)
c i4time         i      i*4     time (seconds since 1-1-1960)
c alt            o      r*4     solar elevation angle (degrees)
c dec            o      r*4     solar declination (degrees)
c hrangle        o      r*4     solar hour angle (degrees)
        include 'trigd.inc'
        character*9 asc9_time

        pi=3.14159265
        rpd=pi/180.
        call make_fnam_lp(i4time,asc9_time,istatus)

        read(asc9_time,1)jd
1       format(2x,i3)

        read(asc9_time,2)ih
2       format(5x,i2)

        read(asc9_time,3)im
3       format(7x,i2)

        is = i4time - ((i4time / 60) * 60)

        eqt=timeq2(i4time)/rpd          ! equation of time (degrees)
        dec=soldec2(i4time)/rpd         ! solar declination (degrees)
        hrangle = (ih-12)*15. + im/4. + is/240. + rlon + eqt

        coszen=sind(rlat)*sind(dec)+cosd(rlat)*cosd(dec)*cosd(hrangle)
        alt = 90. - acosd(coszen)       ! solar altitude (degrees)

        if(hrangle .lt. -180.)hrangle = hrangle + 360.
        if(hrangle .gt. +180.)hrangle = hrangle - 360.

!       write(6,*)'jd,ih,im',jd,ih,im
!       write(6,*)'hrangle,dec,alt',hrangle,dec,alt

        return
        end

        subroutine solalt(rlat,rlon,i4time,alt)

c steve albers jan 1994
c argument      i/o     type                    description
c --------      ---     ----    -----------------------------------------------
c rlat           i      r*4     latitude (degrees)
c rlng           i      r*4     longitude (degrees)
       include 'trigd.inc'
        character*9 asc9_time

        pi=3.14159265
        rpd=pi/180.
        call make_fnam_lp(i4time,asc9_time,istatus)

        read(asc9_time,1)jd
1       format(2x,i3)

        read(asc9_time,2)ih
2       format(5x,i2)

        read(asc9_time,3)im
3       format(7x,i2)

        eqt=timeq2(i4time)/rpd          ! equation of time (degrees)
        dec=soldec2(i4time)/rpd         ! solar declination (degrees)
        hrangle = (ih-12)*15. + im/4. + rlon + eqt

        coszen=sind(rlat)*sind(dec)+cosd(rlat)*cosd(dec)*cosd(hrangle)
        alt = 90. - acosd(coszen)       ! solar altitude (degrees)

!       write(6,*)'jd,ih,im',jd,ih,im
!       write(6,*)'hrangle,dec,alt',hrangle,dec,alt

        return
        end

        subroutine get_solalt_2d(rlat_a,rlon_a,i4time,ni,nj,alt)

c steve albers jan 1994
c argument      i/o     type                    description
c --------      ---     ----    -----------------------------------------------
c rlat           i      r*4     latitude (degrees)
c rlng           i      r*4     longitude (degrees)
        include 'trigd.inc'
        character*9 asc9_time

        real rlat_a(ni,nj)
        real rlon_a(ni,nj)
        real alt(ni,nj)

        pi=3.14159265
        rpd=pi/180.
        call make_fnam_lp(i4time,asc9_time,istatus)

        read(asc9_time,1)jd
1       format(2x,i3)

        read(asc9_time,2)ih
2       format(5x,i2)

        read(asc9_time,3)im
3       format(7x,i2)

        eqt=timeq2(i4time)/rpd          ! equation of time (degrees)
        dec=soldec2(i4time)/rpd         ! solar declination (degrees)

        do i = 1,ni
        do j = 1,nj

          rlat = rlat_a(i,j)
          rlon = rlon_a(i,j)

          hrangle = (ih-12)*15. + im/4. + rlon + eqt

          coszen=sind(rlat)*sind(dec)+cosd(rlat)*cosd(dec)*cosd(hrangle)
          alt(i,j) = 90. - acosd(coszen)       ! solar altitude (degrees)

        enddo ! j
        enddo ! i

!       write(6,*)'jd,ih,im',jd,ih,im
!       write(6,*)'hrangle,dec,alt',hrangle,dec,alt

        return
        end

        subroutine get_solaltaz_2d(rlat_a,rlon_a,i4time,ni,nj,alt,azi)

c steve albers jan 1994
c argument      i/o     type                    description
c --------      ---     ----    -----------------------------------------------
c rlat           i      r*4     latitude (degrees)
c rlng           i      r*4     longitude (degrees)
        include 'trigd.inc'
        character*9 asc9_time

        real rlat_a(ni,nj)
        real rlon_a(ni,nj)
        real alt(ni,nj)
        real azi(ni,nj)

        pi=3.14159265
        rpd=pi/180.
        call make_fnam_lp(i4time,asc9_time,istatus)

        read(asc9_time,1)jd
1       format(2x,i3)

        read(asc9_time,2)ih
2       format(5x,i2)

        read(asc9_time,3)im
3       format(7x,i2)

        eqt=timeq2(i4time)/rpd          ! equation of time (degrees)
        dec=soldec2(i4time)/rpd         ! solar declination (degrees)

        do i = 1,ni
        do j = 1,nj

          rlat = rlat_a(i,j)
          rlon = rlon_a(i,j)

          hrangle = (ih-12)*15. + im/4. + rlon + eqt

          coszen=sind(rlat)*sind(dec)+cosd(rlat)*cosd(dec)*cosd(hrangle)
          alt(i,j) = 90. - acosd(coszen)       ! solar altitude (degrees)

!         call equ_to_altaz_d(solar_dec,solar_ha,lat(i,j)
!    1                       ,altdum,azi(i,j))               

          phi = rlat
          ha = hrangle
          if(ha .lt. 0.)ha = ha + 360.

          sindec = sind(dec)
          cosdec = cosd(dec)
          sinphi = sind(phi)
          cosphi = cosd(phi)
          cosha  = cosd(ha)

!         alt=asind (sinphi*sindec+cosphi*cosdec*cosha)
          cosarg = (cosphi*sindec-sinphi*cosdec*cosha)/cosd(alt(i,j))
          cosarg = min(max(cosarg,-1.),+1.)
          az =acosd(cosarg)

!         if(i .eq. ni/2 .and. j .eq. nj/2)then
!           write(6,*)'ha,az',ha,az
!         endif

          if(ha .gt. 0. .and. ha .lt. 180.)az = 360.0 - az

!         if(i .eq. ni/2 .and. j .eq. nj/2)then
!           write(6,*)'ha,az',ha,az
!         endif

          azi(i,j) = az
          if(azi(i,j) .lt. 0.)azi(i,j) = azi(i,j) + 360.

          if(i .eq. ni/2 .and. j .eq. nj/2)then
            write(6,*)'rlat/rlon',rlat,rlon
            write(6,*)'jd,ih,im',jd,ih,im
            write(6,*)'ha,dec',ha,dec
            write(6,*)'alt,azi',alt(i,j),azi(i,j)
          endif

        enddo ! j
        enddo ! i

        return
        end

c       j. wakefield    28 jan 82       original version

c***these formulas are from paltridge and platt, 1976.  they reference spencer,
c***1971 for the solar declination and equation of time.
c------------------------------------------------------------------------------
        function radnorm(jd)
c***normalized earth-sun distance factor (r0/r)**2
c***jd is input julian day number
        dayang1=2.*3.14159265*(jd-1)/365.
        dayang2=2.*dayang1

        radnorm= 1.000110
     1          +0.034221*cos(dayang1)+0.001280*sin(dayang1)
     2          +0.000719*cos(dayang2)+0.000077*sin(dayang2)

        return
        end
c------------------------------------------------------------------------------
        function soldec(jd)
c***solar declination angle (radians)
c***jd is input julian day number
        dayang1=2.*3.14159265*(jd-1)/365.
        dayang2=2.*dayang1
        dayang3=3.*dayang1

        soldec=  0.006918
     1          -0.399912*cos(dayang1)+0.070257*sin(dayang1)
     2          -0.006758*cos(dayang2)+0.000907*sin(dayang2)
     3          -0.002697*cos(dayang3)+0.001480*sin(dayang3)

        return
        end
c------------------------------------------------------------------------------
        function soldec2(i4time)
c***solar declination angle (radians)
        double precision jd,by

!       compute astronomical julian date
        call i4time_to_jd(i4time,jd,istatus)

!       compute besselian year
        by = 1900.0 + (jd - 2415020.31352) / 365.242198781
        by_frac = by - int(by)

        dayang1=2.*3.14159265*by_frac     
        dayang2=2.*dayang1
        dayang3=3.*dayang1

        soldec2= 0.006918
     1          -0.399912*cos(dayang1)+0.070257*sin(dayang1)
     2          -0.006758*cos(dayang2)+0.000907*sin(dayang2)
     3          -0.002697*cos(dayang3)+0.001480*sin(dayang3)

        return
        end
c------------------------------------------------------------------------------
        function timeq(jd)
c***equation of time (radians)
c***jd is input julian day number
        dayang1=2.*3.14159265*(jd-1)/365.
        dayang2=2.*dayang1

        timeq=   0.000075
     1          +0.001868*cos(dayang1)-0.032077*sin(dayang1)
     2          -0.014615*cos(dayang2)-0.040849*sin(dayang2)

        return
        end
c------------------------------------------------------------------------------
        function timeq2(i4time)
c***equation of time (radians)
        double precision jd,by

!       compute astronomical julian date
        call i4time_to_jd(i4time,jd,istatus)

!       compute besselian year
        by = 1900.0 + (jd - 2415020.31352) / 365.242198781
        by_frac = by - int(by)

        dayang1=2.*3.14159265*by_frac    
        dayang2=2.*dayang1

        timeq2=  0.000075
     1          +0.001868*cos(dayang1)-0.032077*sin(dayang1)
     2          -0.014615*cos(dayang2)-0.040849*sin(dayang2)

        return
        end
