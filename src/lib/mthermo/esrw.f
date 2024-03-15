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

	function esrw(t)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the saturation vapor pressure over liquid
c   water esrw (millibars) given the temperature t (celsius). the
c   formula used is due to richards, j.m., 1971: simple expression
c   for the saturation vapour pressure of water in the range -50 to
c   140c, british journal of applied physics, vol. 4, pp.l15-l18.
c   the formula was quoted more recently by wigley, t.m.l.,1974:
c   comments on 'a simple but accurate formula for the saturation
c   vapor pressure over liquid water,' journal of applied meteorology,
c   vol. 13, no. 5 (august) p.606.

	data cta,ts,ews/273.15,373.15,1013.25/

c   cta = difference between kelvin and celsius temperature
c   ts = temperature of the boiling point of water (k)
c   ews = saturation vapor pressure over liquid water at 100c

	data c1,     c2,     c3,     c4
     1	/ 13.3185,-1.9760,-0.6445,-0.1299 /
	tk = t+cta
	x = 1.-ts/tk
	px = x*(c1+x*(c2+x*(c3+c4*x)))
	vp = ews*exp(px)
	if (vp.lt.0) vp = 0.
	esrw = vp
	return
	end
