cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
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

        SUBROUTINE solar_position(rlat,rlon,i4time,alt,dec,hrangle)

C Steve Albers Jan 1994
C Argument      I/O     Type                    Description
C --------      ---     ----    -----------------------------------------------
C RLAT           I      R*4     Latitude (degrees)
C RLNG           I      R*4     Longitude (degrees)
        character*9 asc9_time

        PI=3.14159265
        RPD=PI/180.
        call make_fnam_lp(i4time,asc9_time,istatus)

        read(asc9_time,1)jd
1       format(2x,i3)

        read(asc9_time,2)ih
2       format(5x,i2)

        read(asc9_time,3)im
3       format(7x,i2)

        is = i4time - ((i4time / 60) * 60)

        EQT=TIMEQ(JD)/rpd               !Equation of Time (Degrees)
        DEC=SOLDEC(JD)/rpd              !Solar declination
        hrangle = (ih-12)*15. + im/4. + is/240. + rlon + EQT

        COSZEN=SIND(RLAT)*SIND(DEC)+COSD(RLAT)*COSD(DEC)*COSD(hrangle)
        alt = 90. - ACOSD(COSZEN)

        if(hrangle .lt. -180.)hrangle = hrangle + 360.
        if(hrangle .gt. +180.)hrangle = hrangle - 360.

!       write(6,*)'jd,ih,im',jd,ih,im
!       write(6,*)'hrangle,dec,alt',hrangle,dec,alt

        RETURN
        END

        SUBROUTINE solalt(rlat,rlon,i4time,alt)

C Steve Albers Jan 1994
C Argument      I/O     Type                    Description
C --------      ---     ----    -----------------------------------------------
C RLAT           I      R*4     Latitude (degrees)
C RLNG           I      R*4     Longitude (degrees)
        character*9 asc9_time

        PI=3.14159265
        RPD=PI/180.
        call make_fnam_lp(i4time,asc9_time,istatus)

        read(asc9_time,1)jd
1       format(2x,i3)

        read(asc9_time,2)ih
2       format(5x,i2)

        read(asc9_time,3)im
3       format(7x,i2)

        EQT=TIMEQ(JD)/rpd               !Equation of Time (Degrees)
        DEC=SOLDEC(JD)/rpd              !Solar declination
        HRANGLE = (ih-12)*15. + im/4. + rlon + EQT

        COSZEN=SIND(RLAT)*SIND(DEC)+COSD(RLAT)*COSD(DEC)*COSD(HRANGLE)
        alt = 90. - ACOSD(COSZEN)

!       write(6,*)'jd,ih,im',jd,ih,im
!       write(6,*)'hrangle,dec,alt',hrangle,dec,alt

        RETURN
        END

C       J. Wakefield    28 Jan 82       Original version

C***These formulas are from Paltridge and Platt, 1976.  They reference Spencer,
C***1971 for the solar declination and equation of time.
C------------------------------------------------------------------------------
        FUNCTION RADNORM(JD)
C***Normalized earth-sun distance factor (R0/R)**2
C***JD is input Julian day number
        DAYANG1=2.*3.14159265*(JD-1)/365.
        DAYANG2=2.*DAYANG1

        RADNORM= 1.000110
     1          +0.034221*COS(DAYANG1)+0.001280*SIN(DAYANG1)
     2          +0.000719*COS(DAYANG2)+0.000077*SIN(DAYANG2)

        RETURN
        END
C------------------------------------------------------------------------------
        FUNCTION SOLDEC(JD)
C***Solar declination angle (radians)
C***JD is input Julian day number
        DAYANG1=2.*3.14159265*(JD-1)/365.
        DAYANG2=2.*DAYANG1
        DAYANG3=3.*DAYANG1

        SOLDEC=  0.006918
     1          -0.399912*COS(DAYANG1)+0.070257*SIN(DAYANG1)
     2          -0.006758*COS(DAYANG2)+0.000907*SIN(DAYANG2)
     3          -0.002697*COS(DAYANG3)+0.001480*SIN(DAYANG3)

        RETURN
        END
C------------------------------------------------------------------------------
        FUNCTION TIMEQ(JD)
C***Equation of time (radians)
C***JD is input Julian day number
        DAYANG1=2.*3.14159265*(JD-1)/365.
        DAYANG2=2.*DAYANG1

        TIMEQ=   0.000075
     1          +0.001868*COS(DAYANG1)-0.032077*SIN(DAYANG1)
     2          -0.014615*COS(DAYANG2)-0.040849*SIN(DAYANG2)

        RETURN
        END
