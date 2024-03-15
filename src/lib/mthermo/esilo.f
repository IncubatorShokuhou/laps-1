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

	function esilo(t)
c
c   this function returns the saturation vapor pressure over ice
c   esilo (millibars) given the temperature t (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the formula is due to lowe, paul r., 1977: an approximating polynomial for 
c   the computation of saturation vapor pressure, journal of applied
c   meteorology, vol. 16, no. 1 (january), pp. 100-103.
c   the polynomial coefficients are a0 through a6.
c
	data a0,a1,a2,a3,a4,a5,a6
     1	/6.109177956,     5.034698970e-01, 1.886013408e-02,
     2	 4.176223716e-04, 5.824720280e-06, 4.838803174e-08,
     3	 1.838826904e-10/
	if (t.le.0.) go to 5
	esilo = 9999.
	write(6,3)esilo
    3	format(' saturation vapor pressure over ice is undefined for',
     1	/' temperature > 0c. esilo =',f6.0)
	return
    5	continue
	esilo = a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t)))))
	return
	end
