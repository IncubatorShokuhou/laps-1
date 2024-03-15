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

	function tw(t,td,p)
c
c   this function returns the wet-bulb temperature tw (celsius)
c   given the temperature t (celsius), dew point td (celsius)
c   and pressure p (mb).  see p.13 in stipanuk (1973), referenced
c   above, for a description of the technique.
c
c	baker,schlatter	17-may-1982	original version
c
c   determine the mixing ratio line thru td and p.
	aw = w(td,p)
c
c   determine the dry adiabat thru t and p.
	ao = o(t,p)
	pi = p
c
c   iterate to locate pressure pi at the intersection of the two
c   curves .  pi has been set to p for the initial guess.
	do 4 i= 1,10
	   x= .02*(tmr(aw,pi)-tda(ao,pi))
	   if (abs(x).lt.0.01) go to 5
 4	   pi= pi*(2.**(x))
c   find the temperature on the dry adiabat ao at pressure pi.
 5	ti= tda(ao,pi)
c
c   the intersection has been located...now, find a saturation
c   adiabat thru this point. function os returns the equivalent 
c   potential temperature (c) of a parcel saturated at temperature
c   ti and pressure pi.
	aos= os(ti,pi)
c   function tsa returns the wet-bulb temperature (c) of a parcel at
c   pressure p whose equivalent potential temperature is aos.
	tw = tsa(aos,p)
	return
	end
