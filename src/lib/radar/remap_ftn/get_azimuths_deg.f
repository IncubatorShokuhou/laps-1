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
      subroutine get_azimuths_deg ( i_tilt,    ! input
     :              max_rays,        ! input # rays dimensioned for remapper
     :              az_array,        ! output
     :              istatus )        ! output
c
c     purpose:
c       fill azimuth array with azimuths from the common buffer.
c
c     author:
c       steve albers, fsl
c     modifications:
c       clean-up and improve efficiency keith brwester, caps, june, 1994
c

      implicit none
c
c     input variables
c
      integer i_tilt
      integer max_rays
c
c     output variables
c
      real    az_array(max_rays)
      integer istatus
c
c     include files
c
      include 'remap_dims.inc'
      include 'remap_buffer.cmn'
c
c     misc internal variables
c
      integer i
c
      if( n_rays_cmn .gt. max_rays ) then
        write(6,805)
 805    format(' n_rays_cmn exceeds max_rays in get_azimuths_deg')
        write(6,810) n_rays_cmn,max_rays
 810    format(2i12)
        istatus = 0
        return
      end if
c
c     transfer common data to output array
c
      do 100 i = 1,n_rays_cmn
        az_array(i) = azim_cmn(i)
 100  continue
      istatus = 1
c
      return
      end
