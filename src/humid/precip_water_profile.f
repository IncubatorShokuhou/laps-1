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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
	subroutine precip_water_profile(q,p,za,n,pw)
c
c************************************************************************
c*									*
c*	module name:	precip_water_profile				*
c*									*
c*	language:	fortran-77	   library:			*
c*	version.rev:	1.0 30 april 96	programmer: kleespies		*
c*									*
c*	calling seq:	call precip_water_profile(q,p,za,n,pw)		*
c*									*
c*	description:	computes precipitable water profilefrom a mixing*
c*			ratio profile and a pressure profile, along	*
c*			a slant path.					*
c*									*
c*			this code largely lifted from radcom.		*
c*									*
c*	input args:	r*4 	q(n)	mixing ratio profile		*
c*			r*4	p(n)	pressure profile		*
c*			r*4	za	zenith angle			*
c*			r*4	n	number of points in profile	*
c*									*
c*	output args:	r*4	precipitable water profile		*
c*									*
c*	common blks:	none						*
c*									*
c*	include:	none						*
c*									*
c*	externals:	none						*
c*									*
c*	data files:	none						*
c*									*
c*	restrictions:	q,p,pw dimensioned in calling routine.		*
c*									*
c*	error codes:	none						*
c*									*
c*	error messages:	none						*
c*									*
c************************************************************************
c
	implicit none

*	input
	integer n
	real q(n),p(n),za  ! modified by dan birkenheuer 8/10/98

*	output
	real pw(n)
	real path, cgrav 
	data cgrav / 5.098581e-4 / ! .5/980.665
	real delp

	integer i	! local variables

	path = cgrav / cos(za*acos(-1.)/180.)

	pw(1) = 0.0
	
	do i = 2 , n
	 delp = abs(p(i) - p(i-1))
	 pw(i) = pw(i-1) + path*(q(i) + q(i-1))*delp
	enddo

	return
	end
