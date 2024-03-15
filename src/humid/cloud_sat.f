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
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine cloud_sat (cg,delsat,data)

c     this routine is the basic cloud saturation model for lq3driver.x
c     this routine is called from the main routine when performing post
c     saturation computations.  it is now modified with test_cloud output.

c     dan birkenheuer  5/10/2002

      implicit none

      real :: cg                ! cloud fraction
      real :: delsat            ! output from test cloud, amount to sat
      real :: data              ! pre_analyzed ssh at location
      real :: data_old          ! data as input to routine (needs to be preserved)
      real :: computed_increment ! amount planned to increase data by, we can compare this to delsat to see if it is too high

      common /cloud_sat_insert/ max_cdelrh, cf_set
      real max_cdelrh
      real cf_set



c      if(cg > 0.6 .and. cg < 1.0) then !cloudy
c         
c         data = cg*sat + (1.- cg)*data
c         data = data + cg*delsat !ramp to saturation if above .6 cf

c      elseif (cg >= 1.0) then 
         
c         data = data + delsat
         
c      endif
c     new function goes here devised by steve albers 1/19/2012
      if (cg .gt. cf_set) then
      data_old = data
      computed_increment = (cg**0.2)*max_cdelrh  ! this is rh fraction
      computed_increment = computed_increment * (data+delsat)  ! this is sh fraction of saturation
      data = data + computed_increment   ! data now modify, must use data_old for comparision to saturation  
      data = min (data, (data_old+delsat) ) ! this caps data to 100% rh or saturated sh 
      endif

      return

      end subroutine cloud_sat
