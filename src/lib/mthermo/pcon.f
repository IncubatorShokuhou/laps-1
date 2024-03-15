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

	function pcon(p,t,tc)

c   this function returns the pressure pcon (mb) at the lifted condensa-
c   tion level, given the initial pressure p (mb) and temperature t
c   (celsius) of the parcel and the temperature tc (celsius) at the 
c   lcl. the algorithm is exact.  it makes use of the formula for the
c   potential temperatures corresponding to t at p and tc at pcon.
c   these two potential temperatures are equal.

c	baker, schlatter  17-may-1982	  original version.

c
	data akapi/3.5037/

c   akapi = (specific heat at constant pressure for dry air) /
c	    (gas constant for dry air)

c   convert t and tc to kelvin temperatures.

	tk = t+273.15
	tck = tc+273.15
	pcon = p*(tck/tk)**akapi
	return
	end
