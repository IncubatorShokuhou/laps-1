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

	function dwpt(t,rh)
c
c   this function returns the dew point (celsius) given the temperature
c   (celsius) and relative humidity (%). 
c
c	baker,schlatter	17-may-1982	original version
c
c	the formula is used in the
c   processing of u.s. rawinsonde data and is referenced in parry, h.
c   dean, 1969: "the semiautomatic computation of rawinsondes,"
c   technical memorandum wbtm edl 10, u.s. department of commerce,
c   environmental science services administration, weather bureau,
c   office of systems development, equipment development laboratory,
c   silver spring, md (october), page 9 and page ii-4, line 460.
	x = 1.-0.01*rh
c   compute dew point depression.
	dpd =(14.55+0.114*t)*x+((2.5+0.007*t)*x)**3+(15.9+0.117*t)*x**14
	dwpt = t-dpd
	return
	end
