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

	function tlcl(t,td)
c   this function yields the temperature tlcl (celsius) of the lifting
c   condensation level, given the temperature t (celsius) and the
c   dew point td (celsius).  the formula used is in bolton, david,
c   1980: "the computation of equivalent potential temperature,"
c   monthly weather review, vol. 108, no. 7 (july), p. 1048, eq.(15).

c	baker, schlatter  17-may-1982	  original version.

c   convert from celsius to kelvin degrees.

	tk = t+273.15
	tdk = td+273.15
	a = 1./(tdk-56.)
	b = alog(tk/tdk)/800.
	tc = 1./(a+b)+56.
	tlcl = tc-273.15
	return
	end
