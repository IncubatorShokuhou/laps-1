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

      subroutine ingest_goessnd (path_to_data, i4time_sys, 
     1     lun_out, istatus)


      implicit none
      
      character*(*) path_to_data
      integer i4time_sys, lun_out, istatus


      integer max_files
      parameter (max_files = 3000)
      integer i4time_begin,i4time_end ! bracket of i4time of interest
      character*256 file_names(max_files) !file names found in path of interest
      character*256 files_to_process(max_files) ! files meeting time criteria
      character*9 a9_test !a9 filename to test for criteria
      integer     i4time_test   ! i4time of filename to test for time criteria
      integer num_files
      integer fname_len         ! length of file name found in path
      integer path_len          ! length of path_to_data
      integer i,j,ii,jj,lct
      real mdf               !missing data flag
      character*50 ext
      character *256 directory
      integer len_dir

      istatus = 0

      

c     get missing data flag
      call get_r_missing_data (mdf,istatus) 

c     define the spatial window for data
c     obtain dimensions of laps domain ii,jj
      call get_grid_dim_xy(ii,jj,istatus)


c     define time window to look for data
      i4time_begin = i4time_sys -4200
      i4time_end   = i4time_sys +3600

c     get the filenames to consider processing
      call s_len (path_to_data,path_len)
      write (6,*) path_to_data(1:path_len)
      call get_file_names (path_to_data(1:path_len), num_files,
     1     file_names,max_files,istatus)

      if(istatus.eq.1) then
         
c     trim file names removing path
         do i=1,num_files
            call s_len(file_names(i),fname_len)
            file_names(i) = file_names(i)(path_len+1:fname_len)
         enddo

      else
         write(6,*) 'error in get_file_names'
         return
      endif

      
c     go through list of filenames and test their valid times against window

      j = 0 ! counter to increment for good files to run through scanner
      do i = 1, num_files
         a9_test = file_names(i)(5:9)//file_names(i)(11:14)
         call i4time_fname_lp(a9_test,i4time_test,istatus)


c         i4time_begin = i4time_begin - 3800 ! inserted to expand time window

         if (i4time_test.ge.i4time_begin .and.
     1        i4time_test.le.i4time_end)then
            j = j+1
            files_to_process(j) = file_names(i)
         endif
      enddo                     !i                    

      if (j.eq.0) then
         write (6,*) 'no files to process'
         istatus = 0
         return
      endif

c     new call for steve albers check of lun
      call open_ext ( lun_out, i4time_sys, 'snd', istatus)

c     what is the new file name length (same for all files)

      call s_len(files_to_process(1),fname_len)

      do i = 1,j                !for all files to be processed

         write (6,*) 'processing data goes data files ',j
         write (6,*) 'filename is ',files_to_process(i)
         call process_goes_snd (path_to_data, path_len, 
     1        files_to_process(i),
     1        fname_len, ii,jj,
     1        i4time_begin,i4time_end,mdf,lun_out,istatus)

      enddo                     ! i

      istatus = 1
      return
      end




