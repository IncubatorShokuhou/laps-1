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
        subroutine   xy_to_met_xm( x,
     1                     y,
     1                     range,
     1                     dir,
     1                     istatus )

!  xy_to_met_xm  converts cartesian coordinates to meteorological (polar)
!                coordinates.

!       bob lipschutz   12-apr-1983     original version
!       windsor         15-aug-1985     fixed so differentiates btwn 0 and 180


        real   deg_per_rad
        parameter (deg_per_rad  = 180. / 3.14159)


        real
     1  x,              ! x-coordinate.
     1  y,              ! y-coordinate (in same units as x).
     1  range,          ! output:  range in same units as x and y).
     1  dir             ! output:  meteorological degrees.


        range  = sqrt( x*x + y*y )

        if ( x .ne. 0)  then

            dir  = atan2( y,x )
!       print *, 'after atan2 : ',dir

            dir  = dir * deg_per_rad
            if( dir .lt. 0 )  dir  = 360. + dir

!         ... convert to met. degrees.

            dir  = 450. - dir
            if( dir .gt. 360. )  dir  = dir - 360.

          else

            if (y .lt. 0) then

                dir = 180.0

             else

               dir  = 0.0

            endif

        end if

!       dir = float(nint(dir))

!       print *,' x,y : ',x,y,' deg : ',dir

        istatus   = 1
        return
        end
