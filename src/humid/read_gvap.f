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
      subroutine read_gvap (filename, nstations, path_to_gvap12,
     1     path_to_gvap10,time_diff,ihop_flag,
     1     lat,lon, wt,w1,w2,w3,gvap_p,
     1     nn, istatus)

c     routine to read goes watervapor for system testing
c     author dan birkenheuer
c     2/8/99
c     
c     paramter list
c     filename, name of the wisconsin file
c     nstations, dimension of variables
c     lat, latitude of each station
c     lon, longitude of each station
c     wt, total water
c     w1, layer 1 water
c     w2, layer 2 water
c     w3, layer 3 water
c     nn, number of stations with real data for consideration
c
      implicit none

c     input variables
      character*9 filename
      integer nstations,nn,istatus,idummy,idummy2
      character*256 path_to_gvap12,path_to_gvap10
      integer time_diff         !time allowed for latency (sec)
      integer ihop_flag
      real lat(nstations)
      real lon(nstations)
      real wt(nstations)
      real w1(nstations)
      real w2(nstations)
      real w3(nstations)
      real gvap_p(nstations)

c     internal variables

      integer i
      real dummy
      character*8 adummy
      integer ptg_index
      character*512 const_file
      integer cf
      character*9 filefound
      character*120 extension
      integer extension_index
      integer istatus_8,istatus_10
      integer file_name_length
      character*120 desired_ext
      integer de_index, n8


      file_name_length = 9
      desired_ext = 'tpw'
      istatus_8=0
      istatus_10 = 0

      call s_len(path_to_gvap12, ptg_index)


c     reading goes 12

c     get most recent file in directory

      write (6,*)
      write (6,*) 'acquiring goes west info'
      write (6,*)



      call get_newest_file (filename, time_diff,file_name_length,
     1     path_to_gvap12,ptg_index,filefound,desired_ext, de_index,
     1     extension, extension_index, istatus)

      if(istatus.ne.1) then     !failure in getting file
         return
      endif

      const_file = path_to_gvap12(1:ptg_index)//'20'//filefound
     1     //'.'//extension(1:extension_index)

      call s_len(const_file, cf)

      write(6,*) 'opening file ',const_file(1:cf)
      
c      if (filefound(1:2) .le. '09') then   ! big endian data
c         open(22, file=const_file(1:cf),form='formatted',
c     1        status='old',convert='swap',err = 668)
c      else ! little endian data
         open(22, file=const_file(1:cf),form='formatted',
     1        status='old',err = 668)
c      endif
      read(22,*,end=668,err=668) ! first header line is ignored
      do i = 1,nstations

         if(ihop_flag.eq.0) then
            read(22,*,end=665,err=665) idummy,idummy,lat(i),lon(i),
     1           idummy,idummy,wt(i), w1(i),w2(i),w3(i),idummy, 
     1           gvap_p(i)
         elseif (ihop_flag .eq. 1) then
            read(22,*,end=665,err=665) adummy,idummy,lat(i),lon(i),
     1           wt(i), w1(i),w2(i),w3(i),dummy,
     1           gvap_p(i)   
         endif  ! ihop special format


      enddo

 665  close (22)
      nn = i-1
      n8 = nn ! number of goes 12 stations read
      if (nn .eq. 0) go to 666
      write(6,*) nn, ' number of records read goes-east ',
     1     const_file(1:cf)
      istatus_8 = 1             !assign success to reading goes 8
c     write(6,*) (wt(i),i=1,nn)

      go to 669

 668  close (22)
      write(6,*) 'failed reading goes-east, ',const_file(1:cf)
      istatus_8 = 0


 669  continue

c     reading goes 10
c     skip reading goes 10 if ihop flag is set to on
      if(ihop_flag .eq. 1) then
         write (6,*) 'ihop activated, skipping goes-west'
         istatus_10 = 0
         go to 700
      endif

c     get most recent file in directory
      
      write (6,*)
      write (6,*) 'acquiring goes west gvap info'
      write (6,*)

      istatus_10 = 1

      call s_len(path_to_gvap10, ptg_index)

      desired_ext = 'tpw'

      call get_newest_file (filename, time_diff,file_name_length,
     1     path_to_gvap10,ptg_index,filefound,desired_ext, de_index,
     1     extension, extension_index, istatus)

      if(istatus.ne.1) then     !failure in getting file
         istatus_10 = 0
         go to 667
      endif

      const_file = path_to_gvap10(1:ptg_index)//'20'//filefound
     1     //'.'//extension(1:extension_index)

      call s_len(const_file, cf)

      write(6,*) 'opening file ',const_file(1:cf)
      
      open(23, file=const_file(1:cf),form='formatted',
     1     status='old',err = 666)
      read(23,*,end=666,err=666) ! first header line is ignored
      do i = nn,nstations
         read(23,*,end=667,err=667) idummy,idummy,lat(i),lon(i),
     1        idummy,idummy,wt(i), w1(i),w2(i),w3(i),idummy, 
     1        gvap_p(i)

      enddo

      write (6,*) 'loop reading goes west ended on dimension limit'
      write (6,*) 'consider expanding dimension size'
      write (6,*) 'current dimension limit is.... ', nstations


 667  close(23)

 700  if (istatus_8.eq.0) then
         write(6,*) 'goes east failed'
      endif
      if (istatus_10.eq.0) then
         write(6,*) 'goes west failed'
      endif

      if(istatus_8+istatus_10.eq.0) then
         return
      else
         write(6,*) 'got something,  processing....'
         istatus = 1
      endif
      
      nn = i-1
      if (nn .eq. 0) go to 666
      write(6,*) nn-n8, ' number of records read goes west'
      istatus = 1
c      write(6,*) (wt(i),i=1,nn)

      write (6,*) 'total gvap records read is ', nn

      return

 666  if(nn.eq.0) then
         istatus = 0
         write(6,*) 'process failure...  no available gvap data'
      else
         istatus = 1
      endif

      return

      end

