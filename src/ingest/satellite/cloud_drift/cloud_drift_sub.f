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

      subroutine get_cloud_drift_data(i4time_sys,i4_window
     1                                    ,nx_l,ny_l
     1                                    ,filename,istatus)

!     steve albers     jan-1998

!.............................................................................

      character*170 filename
      character*10 a10_time
      character*9 a9_timeobs,a10_to_a9 
      character*6 a6_time
      character*5 a5_time
      character*4 a4_time
      character*2 c2_sat_type

      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

      logical l_new_fmt, l_parse
      
      l_new_fmt = .true. ! .false.
      if(l_parse(filename,'goes11'))l_new_fmt = .true.
      if(l_parse(filename,'goes12'))l_new_fmt = .true.

      call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_laps_perimeter'
          return
      endif

      open(1,file=filename,status='old',err=990)
      i = 0

!     read the file...
      if(l_new_fmt)then
          write(6,*)' reading ascii file in new format ',filename
          nheader = 2
      else
          write(6,*)' reading ascii file in old format ',filename
          nheader = 1
      endif

      do iheader = 1,nheader
          read(1,*)               ! skip header lines
      enddo ! iheader

 40   if(.not. l_new_fmt)then
          read(1,11,err=50,end=990)c2_sat_type,a6_time,a4_time
     1                            ,rlat,rlon,ipres,spd,idir
 11       format(a2,15x,a6,2x,a4,f8.2,f9.2,i6,f7.1,i5)
      else ! new format
          read(1,12,err=50,end=990)c2_sat_type,a5_time,a4_time
     1                            ,rlat,rlon,ipres,spd,idir
 12       format(a2,17x,a5,4x,a4,2x,f10.0,f10.0,i10,f10.0,i10)
      endif

 50   continue

      i = i+1
      pres_pa = float(ipres) * 100. ! pascals
      dir = float(idir)
      rlon = -rlon

      if(rlat .le. rnorth .and. rlat .ge. south .and.
     1   rlon .ge. west   .and. rlon .le. east        )then       
          write(6,*)
          write(6,*)' cloud_drift #',i

      else ! outside lat/lon perimeter - reject
!         write(6,*)' lat/lon - reject'       
          goto 900
      endif

      if(c2_sat_type .ne. 'ir' .and. c2_sat_type .ne. 'vi'
     1                         .and. c2_sat_type .ne. 'wv')then       
          write(6,*)' bad sat type ',c2_sat_type
          goto 900
      endif

      if(l_new_fmt)then
          a9_timeobs = a5_time//a4_time
      else
          a10_time = a6_time//a4_time
          a9_timeobs = a10_to_a9(a10_time,istatus)
          if(istatus .ne. 1)then
              write(6,*)' bad observation time - reject ',a10_time       
              goto 900
          endif
      endif

      call cv_asc_i4time(a9_timeobs,i4time_ob)
      i4_resid = abs(i4time_ob - i4time_sys)
      if(i4_resid .gt. i4_window)then ! outside time window
          write(6,*)' time - reject ',a9_timeobs
     1                               ,i4_resid,i4_window
          goto 900        
      endif

      write(6 ,21)rlat,rlon,pres_pa,dir,spd,a9_timeobs
      write(11,21)rlat,rlon,pres_pa,dir,spd,a9_timeobs
 21   format(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

 900  continue
      
      go to 40 ! loop back to read next line

 990  close(1)

      return
      end
