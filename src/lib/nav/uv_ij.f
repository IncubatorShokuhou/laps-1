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
	subroutine uv_ij (ny,u1,v1,du,dv,u,v,ri,rj)

c	this routine takes inital u1 v1 (lower left) map projection coords,
c	the increment in u and v (du and dv) (parent grid)
c	then given any set of u and v it generates real
c	representations of i and j in the parent grid corresponding to the 
c	u and v sets.

c	written by dan birkenheuer february 1994

	implicit none

	real u1,v1,du,dv,u,v,ri,rj
	integer ny

	ri = (u-u1)/du  +1.0
	rj = float(ny) - (v-v1)/dv

	return
	end
