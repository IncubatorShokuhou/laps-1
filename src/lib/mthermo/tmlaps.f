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

	function tmlaps(thetae,p)
c
c   this function returns the temperature tmlaps (celsius) at pressure
c   p (millibars) along the moist adiabat corresponding to an equivalent
c   potential temperature thetae (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the algorithm was written by eric smith at colorado state
c   university.
	data crit/0.1/
c   cta = difference between kelvin and celsius temperatures.
c   crit = convergence criterion (degrees kelvin)
	eq0 = thetae
c   initial guess for solution
	tlev = 25.
c   compute the saturation equivalent potential temperature correspon-
c   ding to temperature tlev and pressure p.
	eq1 = ept(tlev,tlev,p)
	dif = abs(eq1-eq0)
	if (dif.lt.crit) go to 3
	if (eq1.gt.eq0) go to 1
c   dt is the initial stepping increment.
	dt = 10.
	i = -1
	go to 2
    1	dt = -10.
	i = 1
    2	tlev = tlev+dt
	eq1 = ept(tlev,tlev,p)
	dif = abs(eq1-eq0)
	if (dif.lt.crit) go to 3
	j = -1
	if (eq1.gt.eq0) j=1
	if (i.eq.j) go to 2
c   the solution has been passed. reverse the direction of search
c   and decrease the stepping increment.
	tlev = tlev-dt
	dt = dt/10.
	go to 2
    3	tmlaps = tlev
	return
	end
