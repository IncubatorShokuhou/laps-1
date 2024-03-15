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

	function dpt(ew)
c
c   this function returns the dew point dpt (celsius), given the
c   water vapor pressure ew (millibars).
c
c	baker,schlatter	17-may-1982	original version
c
	data es0/6.1078/
c   es0 = saturation vapor pressure (mb) over water at 0c
c   return a flag value if the vapor pressure is out of range.
	if (ew.gt..06.and.ew.lt.1013.) go to 5
	dpt = 9999.
	return
    5	continue
c   approximate dew point by means of teten's formula.
c   the formula appears as eq.(8) in bolton, david, 1980:
c   "the computation of equivalent potential temperature,"
c   monthly weather review, vol 108, no. 7 (july), p.1047.
c   the formula is ew(t) = es0*10**(7.5*t/(t+237.3))
c            or    ew(t) = es0*exp(17.269388*t/(t+237.3))
c   the inverse formula is used below.
	x = alog(ew/es0)
	dnm = 17.269388-x
	t = 237.3*x/dnm
	fac = 1./(ew*dnm)
c   loop for iterative improvement of the estimate of dew point
   10	continue
c   get the precise vapor pressure corresponding to t.
	edp = esw(t)
c   estimate the change in temperature corresponding to (ew-edp)
c   assume that the derivative of temperature with respect to 
c   vapor pressure (dtdew) is given by the derivative of the
c   inverse teten formula.
	dtdew = (t+237.3)*fac
	dt = dtdew*(ew-edp)
	t = t+dt
	if (abs(dt).gt.1.e-04) go to 10
	dpt = t
	return
	end
