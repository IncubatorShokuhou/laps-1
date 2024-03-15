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
cdis
cdis   
cdis
      program test
      implicit none
      real elev,slantkm,slant_range
      real height,range
      real height_grid,r_range
      real rheight_radar
c
c     doviak zrnic stuff
c
      real re,forthre,fthrsqd
      parameter (re=6371 000.,
     :           forthre=(4.*re/3.),
     :           fthrsqd=(forthre*forthre))
c
c     albers' stuff
c
      real rpd,mpd,radius_earth,radius_earth_8_thirds
      real hor_dist,curvature,height_factor
  10  print *, '  enter elev (deg), slant range (km): '
      read(5,*) elev,slantkm
      if(slantkm.eq.0.) stop
      slant_range=1000.*slantkm
      rheight_radar=0.
c
c     doviak zrnic calculation
c
      height=sqrt(slant_range*slant_range + fthrsqd +
     :           2.*forthre*slant_range*sind(elev)) -
     :           forthre
c
      range=forthre*asin((slant_range*cosd(elev)) /
     :                   (forthre + height))
      range=0.001*range
c
      height=height+rheight_radar
c
      print *, ' height, range (d & z ): ',height,range
c
c     albers' calculation
c
       rpd = 3.141592653589/180.
       mpd = 111194.
       radius_earth = 6371.e3
       radius_earth_8_thirds = 6371.e3 * 2.6666666

       hor_dist = slant_range * cosd(elev)

       curvature = hor_dist **2 / radius_earth_8_thirds
       height_grid =
     1  slant_range * sind(elev) + curvature + rheight_radar

       height_factor =
     1   (radius_earth + 0.5 * (rheight_radar + height_grid))
     1  /radius_earth

      r_range = hor_dist / height_factor * .001

      print *, ' height, range (albers): ',height_grid,r_range

      go to 10
      end
