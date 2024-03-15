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

	function tv(t,td,p)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the virtual temperature tv (celsius) of
c   a parcel of air at temperature t (celsius), dew point td
c   (celsius), and pressure p (millibars). the equation appears
c   in most standard meteorological texts.

	data cta,eps/273.15,0.62197/

c   cta = difference between kelvin and celsius temperatures.
c   eps = ratio of the mean molecular weight of water (18.016 g/mole)
c	  to that of dry air (28.966 g/mole)

	tk = t+cta

c   calculate the dimensionless mixing ratio.

	w = .001*wmr(p,td)
	tv = tk*(1.+w/eps)/(1.+w)-cta
	return
	end
