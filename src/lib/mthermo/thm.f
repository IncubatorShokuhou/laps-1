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

	function thm(t,p)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the wet-bulb potential temperature thm
c   (celsius) corresponding to a parcel of air saturated at 
c   temperature t (celsius) and pressure p (millibars).

	f(x) =   1.8199427e+01+x*( 2.1640800e-01+x*( 3.0716310e-04+x*
     1	       (-3.8953660e-06+x*( 1.9618200e-08+x*( 5.2935570e-11+x*
     2	       ( 7.3995950e-14+x*(-4.1983500e-17)))))))
	thm = t
	if (p.eq.1000.) return

c   compute the potential temperature (celsius).

	thd = (t+273.15)*(1000./p)**.286-273.15
	thm = thd+6.071*(exp(t/f(t))-exp(thd/f(thd)))
	return
	end
