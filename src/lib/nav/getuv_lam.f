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
	subroutine getuv_lam (polat,polon,lat,lon,u,v)

c 	derived for lambert conformal conic with one standard parallel
c 	routine assumes all angles in radians... 
c 	compare polat = theta_zero (ncar manual) 
c 	compare polon = phi_zero (ncar)

c	written by dan birkenheuer february 1994
c       j smart 6-96  modified the arguments for trig functions by computing
c 			them first and storing in variable term

	real*4 polat,polon,lat,lon,u,v
	real*8 n,pi
        real*8 term
        real*8 rtanterm,rexpterm

	pi = acos(-1.)
c
        term = pi/2.-polat
	n = cos(term)
        term = pi/4.-lat/2.
	rtanterm = tan(term)
        rexpterm = rtanterm**n
        term = n*(lon-polon)
	u = rexpterm* sin (term)
	v = (-1.)*rexpterm*cos(term)

	return
	end
