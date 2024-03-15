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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine gen_bl_file (i4time,pb,igrid,jgrid,istatus)

c     this routine generates the boundary layer output file determined
c     for the sh mixing algorithm and subsequently for the retroactive
c     inclusion in the surface analysis.
      
c     inputs:
      
c     i4time
c     pb ! top of boundary layer pressure (mb)
      
c     outputs:   none!  this routine generates a file only where
c     filename (ascii time stamp for first record)
c     pb in a simple file (binary write)
      
      implicit none
      
      integer igrid,jgrid
      real pb(igrid,jgrid)
      integer i4time
      integer istatus, len
      character*9 filename
      character*200 pbfile

c     ----- begin exe  last modified 10/26/99 db
      
c     set istatus failure mode
      istatus = 0

      call get_directory('lpbl',pbfile,len)
      pbfile = pbfile(1:len)//'pbl_depth.lpbl'
      
      call make_fnam_lp (i4time,filename,istatus)
      
      if (istatus.eq.1) then
         istatus = 0
      else
         print*, 'error in routine make_fnam_lp'
         istatus = 0
         return
      endif
      
      open (unit=22, file=pbfile, status='replace')
      
      write (22,5) filename,i4time
      write (22,6) pb
 5    format (1x,a9,i12)
 6    format (1x, 7f10.3)
      
      close (22)
      
      istatus = 1
      
      return
      end
