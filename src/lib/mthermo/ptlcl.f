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

	subroutine ptlcl(p,t,td,pc,tc)

c   this subroutine estimates the pressure pc (mb) and the temperature
c   tc (celsius) at the lifted condensation level (lcl), given the
c   initial pressure p (mb), temperature t (celsius) and dew point
c   (celsius) of the parcel.  the approximation is that lines of 
c   constant potential temperature and constant mixing ratio are
c   straight on the skew t/log p chart.
c   teten's formula for saturation vapor pressure as a function of
c   pressure was used in the derivation of the formula below.  for
c   additional details, see math notes by t. schlatter dated 8 sep 81.
c   t. schlatter, noaa/erl/profs program office, boulder, colorado,
c   wrote this subroutine.

c	baker, schlatter  17-may-1982	  original version.

c   akap = (gas constant for dry air) / (specific heat at constant
c	   pressure for dry air)
c   cta = difference between kelvin and celsius temperatures

	data akap,cta/0.28541,273.15/
	c1 = 4098.026/(td+237.3)**2
	c2 = 1./(akap*(t+cta))
	pc = p*exp(c1*c2*(t-td)/(c2-c1))
	tc = t+c1*(t-td)/(c2-c1)
	return
	end
