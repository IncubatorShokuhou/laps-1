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
        subroutine getdudv_pol(lov,dx,dy,lat_orig,lon_orig,
     1	du,dv,u_orig,v_orig)


c	gets du dv of the polar grid
c	gets the origin (lowerl left) corner of the polar grid

        implicit none
        real dx,dy
        real du,dv,u1,u2,v1,v2,dlat,dlon
        real polat,polon,lov
        real pi,dummy,dummy2
        real u_orig,v_orig,lat_orig,lon_orig

        pi = acos(-1.0)
        polon = lov * pi/180.
        polat = 90. * pi/180.  ! by definition

c this routine computes the polar stereo grid given the input
c	determine increment in the across direction. (compute by 1/2 first)

        dlon =  dx/2. * 8.9992800576e-3  / cos(60.*pi/180.) ! km 2 nm 2 degrees
        dummy = (60.)*pi/180.
        dummy2 = polon +(dlon*pi/180.)
        call getuv_pol (polat,polon,dummy ,dummy2,u1,dummy)

        dummy = (60.)*pi/180.
        dummy2 = polon -(dlon*pi/180.)
        call getuv_pol (polat,polon,dummy ,dummy2,u2,dummy)

        du = u1 - u2

c	increment in up dir

        dlat =  dy/2. * 8.9992800576e-3
        dummy = (60.+dlat)*pi/180.
        call getuv_pol (polat,polon,dummy,polon,dummy,v1)

        dummy = (60.-dlat)*pi/180.
        call getuv_pol (polat,polon,dummy,polon,dummy,v2)

        dv = v1-v2

        call getuv_pol (polat,polon,
     1	lat_orig*pi/180.,lon_orig*pi/180.,u_orig,v_orig)


        return
        end
