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


        function sind(x)

        real x,sind,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        sind = sin(x*rpd)

        return
        end

        function cosd(x)

        real x,cosd,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        cosd = cos(x*rpd)

        return
        end

        function tand(x)

        real x,tand,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        tand = tan(x*rpd)

        return
        end

        function asind(x)

        real x,asind,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        asind = asin(x) / rpd

        return
        end

        function acosd(x)

        real x,acosd,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        acosd = acos(x) / rpd

        return
        end


        function atan2d(x,y)

        real x,y,atan2d,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        atan2_deg = atan2(x,y) / rpd

        return
        end








