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

      subroutine ingest_vad(istatus)

!     steve albers      dec-1997        

!     input file 
      character*170 filename_in
      character*9 a9_time
      character*150 dir_in
      character*255 c_filespec
      integer max_files
      parameter(max_files = 1000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)

!     output file
      character*31    ext

      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999

      call get_grid_dim_xy(nx_l,ny_l,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting horizontal domain dimensions'
          go to 999
      endif
 
      call get_laps_cycle_time(ilaps_cycle_time,istatus)
      if(istatus .eq. 1)then
          write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
      else
          write(6,*)' error getting laps_cycle_time'
          istatus = 0
          return
      endif

      iopen = 0

      lag_time_report = 4 * 3600

!     get list of input /public netcdf files
      call get_vad_parms(dir_in,istatus)
      if(istatus .eq. 1)then
          write(6,*)'dir_in = ',dir_in
      else
          write(6,*)' error getting dir_in'
          istatus = 0
          return
      endif

      call s_len(dir_in,len_dir_in)
      c_filespec = dir_in(1:len_dir_in)//'/'//'*0030r'
      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

!     loop through /public netcdf files and choose ones in time window
      write(6,*)' # of files on /public = ',i_nbr_files_ret
      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          filename_in = dir_in(1:len_dir_in)//'/'//a9_time//'0030r'

!         test whether we want the netcdf file for this time
          i4time_file_earliest = i4time_sys - (ilaps_cycle_time / 2)
          i4time_file_latest =   i4time_sys + (ilaps_cycle_time / 2) 
     1                                      + lag_time_report
          
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)' file is too early ',a9_time,' ',i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)' file is too late ',a9_time,' ',i
          else
              write(6,*)' file is in time window ',a9_time,i
              write(6,*)' input file ',filename_in

              if(iopen .eq. 0)then
!                 open output pro file
                  iopen = 1
                  ext = 'pro' ! 'vad'
                  call open_ext(1,i4time_sys,ext(1:3)
     1                                         ,istatus)
                  if(istatus .ne. 1)return
              endif

!             read from the netcdf pirep file and write to the opened pin file
              call get_vad_data(i4time_sys,ilaps_cycle_time
     1                                      ,nx_l,ny_l
     1                                      ,filename_in,istatus)
          endif
      enddo

 999  continue

      write(6,*)' end of ingest_vad subroutine'

      istatus = 1
      return
      end
 
 
       subroutine get_vad_parms(path_to_vad,istatus)

       character*150 path_to_vad
       namelist /vad_nl/ path_to_vad
 
       character*150 static_dir,filename
 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/vad.nl'
 
       open(21,file=filename,status='old',err=900)
       read(21,vad_nl,err=901)
       close(21)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading vad_nl in ',filename
       istatus = 0
       return

       end
