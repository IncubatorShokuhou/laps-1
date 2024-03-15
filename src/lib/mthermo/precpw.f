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

	function precpw(td,p,n)

c	baker, schlatter  17-may-1982	  original version.

c   this function computes total precipitable water precpw (cm) in a
c   vertical column of air based upon sounding data at n levels:
c	   td = dew point (celsius)
c	   p = pressure (millibars)
c   calculations are done in cgs units.

	dimension td(n),p(n)

c   g = acceleration due to the earth's gravity (cm/s**2)

	data g/980.616/

c   initialize value of precipitable water

	pw = 0.
	nl = n-1

c   calculate the mixing ratio at the lowest level.

	wbot = wmr(p(1),td(1))
	do 5 i=1,nl
	wtop = wmr(p(i+1),td(i+1))

c   calculate the layer-mean mixing ratio (g/kg).

	w = 0.5*(wtop+wbot)

c   make the mixing ratio dimensionless.

	wl = .001*w

c   calculate the specific humidity.

	ql = wl/(wl+1.)

c   the factor of 1000. below converts from millibars to dynes/cm**2.

	dp = 1000.*(p(i)-p(i+1))
	pw = pw+(ql/g)*dp
	wbot = wtop
    5	continue
	precpw = pw
	return
	end
