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

	function heatl(key,t)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the latent heat of
c		evaporation/condensation         for key=1
c		melting/freezing		 for key=2
c		sublimation/deposition		 for key=3
c   for water. the latent heat heatl (joules per kilogram) is a 
c   function of temperature t (celsius). the formulas are polynomial
c   approximations to the values in table 92, p. 343 of the smithsonian
c   meteorological tables, sixth revised edition, 1963 by roland list.
c   the approximations were developed by eric smith at colorado state
c   university.
c   polynomial coefficients

	data a0,a1,a2/ 3337118.5,-3642.8583, 2.1263947/
	data b0,b1,b2/-1161004.0, 9002.2648,-12.931292/
	data c0,c1,c2/ 2632536.8, 1726.9659,-3.6248111/
	hltnt = 0.
	tk = t+273.15
	if (key.eq.1) hltnt=a0+a1*tk+a2*tk*tk
	if (key.eq.2) hltnt=b0+b1*tk+b2*tk*tk
	if (key.eq.3) hltnt=c0+c1*tk+c2*tk*tk
	heatl = hltnt
	return
	end
