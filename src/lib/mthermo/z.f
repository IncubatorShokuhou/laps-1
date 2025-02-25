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

	function z(pt,p,t,td,n)

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

c   this function returns the thickness of a layer bounded by pressure
c   p(1) at the bottom and pressure pt at the top.
c   on input:
c	p = pressure (mb).  note that p(i).gt.p(i+1).
c	t = temperature (celsius)
c	td = dew point (celsius)
c	n = number of levels in the sounding and the dimension of
c	    p, t and td
c   on output:
c	z = geometric thickness of the layer (m)
c   the algorithm involves numerical integration of the hydrostatic
c   equation from p(1) to pt. it is described on p.15 of stipanuk
c   (1973).

	dimension t(n),p(n),td(n),tk(n)

c	c1 = .001*(1./eps-1.) where eps = .62197 is the ratio of the
c			      molecular weight of water to that of
c			      dry air. the factor 1000. converts the
c			      mixing ratio w from g/kg to a dimension-
c			      less ratio.
c	c2 = r/(2.*g) where r is the gas constant for dry air
c		      (287 kg/joule/deg k) and g is the acceleration
c		      due to the earth's gravity (9.8 m/s**2). the
c		      factor of 2 is used in averaging two virtual
c		      temperatures.

	data c1/.0006078/,c2/14.64285/
	do 5 i= 1,n
	   tk(i)= t(i)+273.15
    5	continue
	z= 0.0
	if (pt.lt.p(n)) go to 20
	i= 0
   10	i= i+1
	j= i+1
	if (pt.ge.p(j)) go to 15
	a1= tk(j)*(1.+c1*w(td(j),p(j)))
	a2= tk(i)*(1.+c1*w(td(i),p(i)))
	z= z+c2*(a1+a2)*(alog(p(i)/p(j)))
	go to 10
   15	continue
	a1= tk(j)*(1.+c1*w(td(j),p(j)))
	a2= tk(i)*(1.+c1*w(td(i),p(i)))
	z= z+c2*(a1+a2)*(alog(p(i)/pt))
	return
 20	z= -1.0
	return
	end
