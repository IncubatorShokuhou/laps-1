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

	function esice(t)
c
c   this function returns the saturation vapor pressure with respect to
c   ice esice (millibars) given the temperature t (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the formula used is based upon the integration of the clausius-
c   clapeyron equation by goff and gratch.  the formula appears on p.350
c   of the smithsonian meteorological tables, sixth revised edition,
c   1963.
c
	data cta,eis/273.16,6.1071/
c   cta = difference between kelvin and celsius temperature
c   eis = saturation vapor pressure (mb) over a water-ice mixture at 0c
	data c1,c2,c3/9.09718,3.56654,0.876793/
c   c1,c2,c3 = empirical coefficients in the goff-gratch formula
	if (t.le.0.) go to 5
	esice = 99999.
	write(6,3)esice
    3	format(' saturation vapor pressure for ice cannot be computed',
     1	       /' for temperature > 0c. esice =',f7.0)
	return
    5	continue
c   freezing point of water (k)
	tf = cta
	tk = t+cta
c   goff-gratch formula
	rhs = -c1*(tf/tk-1.)-c2*alog10(tf/tk)+c3*(1.-tk/tf)+alog10(eis)
	esi = 10.**rhs
	if (esi.lt.0.) esi = 0.
	esice = esi
	return
	end
