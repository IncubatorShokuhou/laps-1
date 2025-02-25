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

	function alcl(t,td,p)

c	g.s. stipanuk     1973      	  original version.
c	reference stipanuk paper entitled:
c            "algorithms for generating a skew-t, log p
c	     diagram and computing selected meteorological
c	     quantities."
c	     atmospheric sciences laboratory
c	     u.s. army electronics command
c	     white sands missile range, new mexico 88002
c	     33 pages
c	baker, schlatter  17-may-1982	 

c   this function returns the pressure alcl (mb) of the lifting conden-
c   sation level (lcl) for a parcel initially at temperature t (celsius)
c   dew point td (celsius) and pressure p (millibars). alcl is computed
c   by an iterative procedure described by eqs. 8-12 in stipanuk (1973),
c   pp.13-14.
c   determine the mixing ratio line through td and p.

	aw = w(td,p)

c   determine the dry adiabat through t and p.

	ao = o(t,p)

c   iterate to locate pressure pi at the intersection of the two
c   curves. pi has been set to p for the initial guess.

 3	continue
	   pi = p
	   do 4 i= 1,10
	      x= .02*(tmr(aw,pi)-tda(ao,pi))
	      if (abs(x).lt.0.01) go to 5
 4	      pi= pi*(2.**(x))
 5	   alcl= pi
	return

	entry alclm(t,td,p)
c   for entry alclm only, t is the mean potential temperature (celsius)
c   and td is the mean mixing ratio (g/kg) of the layer containing the
c   parcel.

	aw = td
	ao = t
	go to 3
	end
