cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

      program diff_test

        integer imax,jmax,kmax,nlvl

        include 'lapsparms.for'

!       parameter (imax=79)
!       parameter (jmax=73)
        parameter (imax=nx_l)
        parameter (jmax=ny_l)
        parameter (kmax=300)
        parameter (nlvl=300)
      
      integer*4	i4time,		!i4time of data
     1		lvl(nlvl),        !level of each field (4 digit max)
     1		lvl_avail(nlvl),        !level of each field (4 digit max)
     1		i,j,k,
     1		istatus
c
      real*4	data1(imax,jmax,kmax)	!raw data to be written
      real*4	data2(imax,jmax,kmax)	!raw data to be written
c
      character*70	dir_in		!directory to be written to
      character*70	dir_out		!directory to be written to
      character*31	ext		!file name ext (up to 31 chars)
      character*3	var(nlvl) 	        !3 letter id of each field
      character*3	laps_var_avail(nlvl) 	!3 letter id of each field
      character*19      var_avail(nlvl)
      character*4	lvl_coord(nlvl)	!vertical coordinate for each field
      character*10	units(nlvl)	!units of each field
!     character*250	comment1(nlvl)	!comments for each field
!     character*210	comment2(nlvl)	!comments for each field
      character*125	comment1(nlvl)	!comments for each field
      character*125	comment2(nlvl)	!comments for each field
      character*17	asctime
      
      character*4	version
      character*131	model 		!meteorological model in file
      character*131	origin		!location where file was created
      character*11    laps_dom_file   !name of domain file e.g. nest7grid
      logical         l_packed_data

      logical l_pass, l_is_vxx

      character*17 filename
      character*9 a9_time
      character*3 var_last

      real machine_factor

      var_last = '   '
      idiff_msg_flag = 0
      diff_max_all = 0.
      diff_max_all_rel = 0.
      n_files = 0
      ndiff_all = 0
      l_pass = .true.
      num_diff_field_thresh = 0

      write(6,*)' enter 1 if comparing different machines, 0 if same'
      read(5,*)machine_factor

      write(6,*)' filename?'
5     read(5,1)filename
1     format(a)

      if(filename(1:3) .eq. 'end')then
!         write(6,*)' end'
          goto 999
      endif

      read(5,1,err=999)dir_in
      read(5,1,err=999)dir_out

      write(6,*)' machine_factor = ',machine_factor

      if(filename(1:6) .eq. 'static')then
          if(.false.)then
              khmax = 24
              var(1) = 'lat'
              var(2) = 'lon'
              var(3) = 'avg'
              var(4) = 'std'
              var(5) = 'zin'
              var(6) = 'ldf'
              var(7) = 'lnd'
              var(8) = 'use'
              var(9) = 'sln'
              var(10) = 'slt'
              var(11) = 'stl'
              var(12) = 'sbl'
              var(13) = 'a01'
              var(14) = 'a02'
              var(15) = 'a03'
              var(16) = 'a04'
              var(17) = 'a05'
              var(18) = 'a06'
              var(19) = 'a07'
              var(20) = 'a08'
              var(21) = 'a09'
              var(22) = 'a10'
              var(23) = 'a11'
              var(24) = 'a12'
          else
              khmax = 38
              call get_gridgen_var(nlvl,khmax,var,comment1)
          endif
          ext = 'nest7grid'
          call rd_laps_static(dir_in,ext,imax,jmax,khmax,var,units
     1                     ,comment1,data1,grid_spacing_m,istatus)

          if(istatus .ne. 1)then
              write(6,*)' bad status reading 1st static file'
              stop
          endif

          call rd_laps_static(dir_out,ext,imax,jmax,khmax,var,units
     1                     ,comment2,data2,grid_spacing_m,istatus)

          if(istatus .ne. 1)then
              write(6,*)' bad status reading 2nd static file'
              stop
          endif

          thresh_write_pair = .01
          thresh_count_diff = .01

          ihmax = imax
          jhmax = jmax
 
      else

          call s_len(filename,lenf)
          if(lenf .eq. 13)then
              a9_time = filename(1:9)
              ext = filename(11:13)
          else
              a9_time = filename(1:9)
              ext = filename(15:17)
          endif
          call downcase(ext,ext)
          call cv_asc_i4time(a9_time,i4time)

          if(dir_in(1:3) .eq. 'end')then    
!             write(6,*)' end'
              goto 999
          endif

          do i = 1,nlvl
              var(i) = '   '
          enddo

          if(.true.)then
              ihmax = imax
              jhmax = jmax
              call get_laps_dimensions(nk,istatus)
              if(istatus .ne. 1)then
                  write(6,*)' bad status returned from get_laps_dims'
                  stop
              endif

              call rlh(ext,nk,var,lvl,khmax,istatus)
              if(istatus .ne. 1)then
                  write(6,*)' bad status returned from rlh'
                  stop
              endif

          else
!             call read_laps_header(
!    1                         i4time,dir_in,ext,ihmax,jhmax,khmax,
!    1                         laps_dom_file,asctime,version,
!    1                         model,origin,var,lvl,num_variables,
!    1                         var_avail,laps_var_avail,num_levels,
!    1                         lvl_avail,lvl_coord,units,
!    1                         comment1,l_packed_data,istatus)

          endif
 
!         for this extension, set default values of:
!         thresh_write_pair, thresh_count_diff, or num_diff_field_thresh

          thresh_write_pair = 1e-05
          thresh_count_diff = 0.
          num_diff_field_thresh = 10

!         for this particular extension, set values of: 
!         thresh_write_pair, thresh_count_diff, or num_diff_field_thresh

          if(ext(1:3) .eq. 'lc3')then
              do i = 1,42
                  var(i) = 'lc3'
                  lvl(i) = i
              enddo
              thresh_write_pair = .01
              thresh_count_diff = .01
              num_diff_field_thresh = 20

          elseif(ext(1:3) .eq. 'lcp')then
              thresh_write_pair = .01
              thresh_count_diff = .01
              num_diff_field_thresh = 20

          elseif(ext(1:3) .eq. 'lcv')then
              thresh_write_pair = .01
              thresh_count_diff = .01
              num_diff_field_thresh = 25

          elseif(ext(1:3) .eq. 'lco')then
              thresh_write_pair = .01
              thresh_count_diff = .01
              num_diff_field_thresh = 25

          elseif(ext(1:3) .eq. 'lcb')then
              thresh_write_pair = 1.0
              thresh_count_diff = 1.0

          elseif(ext(1:3) .eq. 'lwc')then
              thresh_write_pair = .000001
              thresh_count_diff = .000001

          elseif(ext(1:3) .eq. 'lw3')then
              thresh_write_pair = .01
              thresh_count_diff = .1

          elseif(ext(1:3) .eq. 'liw')then
              thresh_write_pair = .001
              thresh_count_diff = .01

          elseif(ext(1:3) .eq. 'lwm')then
              thresh_write_pair = .01
              thresh_count_diff = .1

          elseif(ext(1:3) .eq. 'lps')then
              thresh_write_pair = .01
              thresh_count_diff = .1

          elseif(ext(1:3) .eq. 'lhe')then
              thresh_write_pair = .001
              thresh_count_diff = .001

          elseif(ext(1:3) .eq. 'lil')then
              thresh_write_pair = .0001
              thresh_count_diff = .0001

          elseif(ext(1:3) .eq. 'lst')then
              thresh_write_pair = .05
              thresh_count_diff = .05

          elseif(ext(1:3) .eq. 'lvd')then
              thresh_write_pair = 1.
              thresh_count_diff = 1.

          elseif(ext(1:3) .eq. 'vrc')then
              thresh_write_pair = 0.1
              thresh_count_diff = 0.1

          elseif(ext(1:3) .eq. 'l1s')then
              do i = 1,4
                  lvl(i) = 0
              enddo
              var(1) = 'r01'
              var(2) = 'rto'
              var(3) = 's01'
              var(4) = 'sto'
              thresh_write_pair = .0001
              thresh_count_diff = .0001

          elseif(ext(1:3) .eq. 'lps')then
              do i = 1,21
                  lvl(i) = 1150 - 50 * i            
              enddo

          elseif(ext(1:3) .eq. 'lsx')then
              thresh_write_pair = .01
              thresh_count_diff = .01
              num_diff_field_thresh = 20

          endif

!         adjust values for this extension (dependent on machine) of: 
!         thresh_count_diff and num_diff_field_thresh

          if(machine_factor .eq. 0.)then
              thresh_count_diff = 0.
              num_diff_field_thresh = 0
          else
              thresh_count_diff = thresh_count_diff * machine_factor
          endif

          if(var(1) .eq. 'rh')var(1) = 'lhe'

!         read first file
          call s_len(dir_in,len_dir_in)

          if(lenf .eq. 13)then
              write(6,*)' reading: ',dir_in(1:len_dir_in),a9_time
     1                           ,'.',ext(1:3)

              call read_laps_data(i4time,dir_in,ext,ihmax,jhmax,khmax,
     1             khmax,var,lvl,lvl_coord,units,comment1,data1,
     1             istatus)
          else ! lga/lgb files
              write(6,*)' reading: ',dir_in(1:len_dir_in),a9_time
     1                           ,'.',ext(1:3)
              i4_valid = i4time + 7200
              call read_laps(i4time,i4_valid,dir_in,ext,
     1             ihmax,jhmax,khmax,khmax,var,lvl,lvl_coord,
     1             units,comment1,data1,
     1             istatus)
!             if(istatus .eq. 0)istatus = 1
          endif

          if(istatus .ne. 1)then
            if(.not. l_is_vxx(ext) )then
              if(l_pass)then
                  write(6,*)' read error: overall criteria failure'
                  l_pass = .false.     
              endif
            else
              write(6,*)' attempting compressed radar data read'
              do if = 1,3
                  write(6,*)' if = ',if
                  kp = 1 + ((if-1) * nk)
                  call read_laps_compressed(i4time,dir_in,ext
     1                               ,ihmax,jhmax,nk
     1                               ,var(kp),lvl(kp),lvl_coord(kp)
     1                               ,units(kp),comment1(kp)
     1                               ,data1(1,1,kp),istatus)
              enddo ! if
              if(istatus .ne. 1)then
                  write(6,*)' read error: overall criteria failure'
                  l_pass = .false.     
              endif
            endif
          endif

!         read second file
          call s_len(dir_out,len_dir_out)

          if(lenf .eq. 13)then  
              write(6,*)' reading: ',dir_out(1:len_dir_out),a9_time
     1                     ,'.',ext(1:3)
              call read_laps_data(i4time,dir_out,ext,ihmax,jhmax,khmax,
     1             khmax,var,lvl,lvl_coord,units,comment2,data2,
     1             istatus)
          else ! lga/lgb
              write(6,*)' reading: ',dir_out(1:len_dir_out),a9_time
     1                           ,'.',ext(1:3)
              i4_valid = i4time + 7200
              call read_laps(i4time,i4_valid,dir_out,ext,
     1             ihmax,jhmax,khmax,khmax,var,lvl,lvl_coord,
     1             units,comment2,data2,
     1             istatus)
!             if(istatus .eq. 0)istatus = 1
          endif

          if(istatus .ne. 1)then
            if(.not. l_is_vxx(ext) )then
              if(l_pass)then
                  write(6,*)' read error: overall criteria failure'
                  l_pass = .false.     
              endif
            else
              write(6,*)' attempting compressed radar data read'
              do if = 1,3
                  write(6,*)' if = ',if
                  kp = 1 + ((if-1) * nk)
                  call read_laps_compressed(i4time,dir_out,ext
     1                               ,ihmax,jhmax,nk
     1                               ,var(kp),lvl(kp),lvl_coord(kp)
     1                               ,units(kp),comment2(kp)
     1                               ,data2(1,1,kp),istatus)
              enddo ! if
              if(istatus .ne. 1)then
                  write(6,*)' read error: overall criteria failure'
                  l_pass = .false.     
              endif
            endif
          endif

      endif ! static file

!       write(6,*)' hit return to continue'
!	read(5,*)

        thresh_write_pair = thresh_write_pair * machine_factor

        diff_max_file_rel = 0.
        diff_max_file = 0.
        diff_max_var = 0.
        ndiff_file = 0
        nvar = 1
        n_levels = 0

        do k = 1,khmax

!       test whether we should switch variables within this file
        if(k .gt. 1)then
            if(var(k) .ne. var(k-1))then
                if(n_levels .gt. 1)then
                    write(6,*)' max diff for variable ',var(k-1)(1:3)
     1                      ,' =',diff_max_var
                    write(6,*)
                    diff_max_var = 0.
                    nvar = nvar + 1
                endif
                n_levels = 0
            endif
        endif

!       test if new variable
        if(var(k)(1:3) .ne. var_last)then

            write(6,*)
            write(6,*)' new var = ',var(k)(1:3)

!           for this variable, set values of:
!           thresh_write_pair, thresh_count_diff, or num_diff_field_thresh

            if    (ext(1:3) .eq. 'lvd')thresh_write_pair = 1.0

            if    (ext(1:3) .eq. 'lt1' .and. var(k)(1:2) .eq. 'ht')then       
                thresh_count_diff = 2.0 * machine_factor
            elseif(ext(1:3) .eq. 'lt1' .and. var(k)(1:2) .eq. 't3')then
                thresh_count_diff = .02 * machine_factor
            elseif(ext(1:3) .eq. 'lga' .and. var(k)(1:2) .eq. 'ht')then
                thresh_count_diff = .01 * machine_factor
            elseif(ext(1:3) .eq. 'lga' .and. var(k)(1:2) .eq. 't3')then
                thresh_count_diff = .001 * machine_factor
            elseif(ext(1:3) .eq. 'lgb' .and. var(k)(1:3) .eq. 'usf')then
                thresh_count_diff = .001 * machine_factor
            elseif(ext(1:3) .eq. 'lgb' .and. var(k)(1:3) .eq. 'vsf')then
                thresh_count_diff = .001 * machine_factor
            elseif(ext(1:3) .eq. 'lgb' .and. var(k)(1:3) .eq. 'tsf')then
                thresh_count_diff = .001 * machine_factor
            elseif(ext(1:3) .eq. 'lgb' .and. var(k)(1:3) .eq. 'dsf')then
                thresh_count_diff = .001 * machine_factor
            elseif(ext(1:3) .eq. 'lgb' .and. var(k)(1:3) .eq. 'tgd')then
                thresh_count_diff = .002 * machine_factor
            elseif(ext(1:3) .eq. 'lgb' .and. var(k)(1:3) .eq. 'slp')then
                thresh_count_diff = 15.0 * machine_factor
            elseif(ext(1:3) .eq. 'lgb' .and. var(k)(1:3) .eq. 'psf')then
                thresh_count_diff = 15.0 * machine_factor
            elseif(ext(1:3) .eq. 'lgb' .and. var(k)(1:3) .eq. 'p')then
                thresh_count_diff = 15.0 * machine_factor
            elseif(ext(1:3) .eq. 'lmt' .and. var(k)(1:3) .eq. 'lmt')then
                thresh_count_diff = 2.0 * machine_factor
            elseif(ext(1:3) .eq. 'lmt' .and. var(k)(1:3) .eq. 'llr')then
                thresh_count_diff = 2.0 * machine_factor
            elseif(ext(1:3) .eq. 'lsx' .and. var(k)(1:1) .eq. 'u')then
                thresh_count_diff = .02 * machine_factor
            elseif(ext(1:3) .eq. 'lsx' .and. var(k)(1:1) .eq. 'v')then
                thresh_count_diff = .02 * machine_factor
            elseif(ext(1:3) .eq. 'lsx' .and. var(k)(1:1) .eq. 'p')then
                thresh_count_diff = 15.0 * machine_factor
            elseif(ext(1:3) .eq. 'lsx' .and. var(k)(1:3) .eq. 'msl')then
                thresh_count_diff = 40.0 * machine_factor
            elseif(ext(1:3) .eq. 'lhe' .and. var(k)(1:3) .eq. 'lhe')then
                thresh_count_diff = 0.1 * machine_factor
            elseif(ext(1:3) .eq. 'lst' .and. var(k)(1:3) .eq. 'wb0')then
                thresh_count_diff = 5.0 * machine_factor
            elseif(ext(1:3) .eq. 'lvd' .and. var(k)(1:3) .eq. 'svn')then
                thresh_count_diff = 10.0 * machine_factor
            elseif(ext(1:3) .eq. 'lvd' .and. var(k)(1:3) .eq. 'alb')then
                thresh_count_diff = 0.01 * machine_factor
                thresh_write_pair = 0.01              !jsmart addition 11-1-96
!           elseif(ext(1:3) .eq. 'lw3' .and. var(k)(1:2) .ne. 'om')then
!               thresh_count_diff = .1 * machine_factor
            elseif(ext(1:9) .eq. 'nest7grid'
     1                                .and. var(k)(1:3) .eq. 'lat')then
                thresh_count_diff = .0005 * machine_factor
            elseif(ext(1:9) .eq. 'nest7grid'
     1                                .and. var(k)(1:3) .eq. 'lon')then
                thresh_count_diff = .0005 * machine_factor
            elseif(ext(1:9) .eq. 'nest7grid'
     1                                .and. var(k)(1:3) .eq. 'avg')then
                thresh_count_diff = 10.0 * machine_factor
            elseif(ext(1:9) .eq. 'nest7grid'
     1                                .and. var(k)(1:3) .eq. 'ldf')then
                thresh_count_diff = .005 * machine_factor
            endif

            write(6,*)
     1      ' threshold to write (first ten) grid point pairs = '
     1                                               ,thresh_write_pair
            write(6,*)
     1      ' threshold to count lvl grid point differences   = '
     1                                               ,thresh_count_diff
            write(6,*)
     1      ' max allowed count of lvl grid point differences = '
     1                                           ,num_diff_field_thresh
            write(6,*)


        endif

        n_levels = n_levels + 1

        diff_max_field = 0.
        diff_max_field_rel = 0.
        abs_value_max = 0.
        imaxd = 0
        jmaxd = 0
        iwrite = 0
        ndiff = 0
        inan = 0
        ndiff_msg = 0
        sumsq = 0.
        n_sq = 0

	do i = 1,ihmax
        do j = 1,jhmax

          call check_nan(data1(i,j,k),istatus_1)
          call check_nan(data2(i,j,k),istatus_2)

          if(istatus_1 .eq. 0 .or. istatus_2 .eq. 0)then
            iwrite = iwrite + 1
            if(iwrite .le. 10)then
                write(6,21)i,j,k,' nan'
            endif
            inan = inan + 1
          else
            diff     = abs(data1(i,j,k)-data2(i,j,k))

!           test if one of the points is missing and the other isn't
            if(   (data1(i,j,k) .eq. r_missing_data .or.
     1             data2(i,j,k) .eq. r_missing_data      )        
     1                          .and.
     1                     diff .gt. 0.                     )then

                ndiff_msg = ndiff_msg + 1
                idiff_msg_flag = 1

            else ! both data points are non-missing

                diff_max_file = max(diff_max_file,diff)
                diff_max_var = max(diff_max_var,diff)
                if(diff .gt. diff_max_field)then
                    diff_max_field = diff
                    imaxd = i
                    jmaxd = j
                endif

                sumsq = sumsq + diff**2
                n_sq = n_sq + 1

            endif

            if(data1(i,j,k) .ne. r_missing_data)then
                abs_value_max = max(abs_value_max,abs(data1(i,j,k)))
            endif

            if(data2(i,j,k) .ne. r_missing_data)then
                abs_value_max = max(abs_value_max,abs(data2(i,j,k)))
            endif


            if(diff .gt. thresh_count_diff)then
                ndiff = ndiff + 1
                ndiff_file = ndiff_file + 1
                ndiff_all = ndiff_all + 1
            endif

            if(diff .gt. thresh_write_pair)then
                iwrite = iwrite + 1 
                if(iwrite .le. 10)then
                    write(6,21,err=22)i,j,k,data1(i,j,k),data2(i,j,k)
     1                                                  ,diff
21                  format(1x,3i5,2f14.6,f12.6)
22              endif
            endif
          endif ! nan test
        enddo ! j
        enddo ! i

        if(n_sq .gt. 0)then
           rms = sqrt(sumsq / float(n_sq))
        else
           rms = 0.
        endif

        if(comment1(k)(1:80) .ne. comment2(k)(1:80))then
            write(6,*)' comments differ at level ',k
            write(6,*)'parallel    comment',trim(comment1(k)(1:80))
            write(6,*)'operational comment',trim(comment2(k)(1:80))
        else
            write(6,*)' comments similar at level ',k
            write(6,*)'parallel    comment',trim(comment1(k)(1:80))
            write(6,*)'operational comment',trim(comment2(k)(1:80))
        endif
        if(inan .gt. 0)write(6,*)' # of nans = ',inan

        if(abs_value_max .gt. 0.)then
!          if(ndiff_msg .eq. 0)then
               diff_max_field_rel = diff_max_field / abs_value_max
!          endif
        else
           diff_max_field_rel = 0.
        endif

        diff_max_file_rel  = max(diff_max_file_rel,diff_max_field_rel)

!	write(6,*)' df_mx - fld #',k,' ',var(k)
!    1  ,lvl(k),' abs/rel/rms/#',diff_max_field,diff_max_field_rel
!    1  ,rms,ndiff,'at',imaxd,jmaxd

 	write(6,101)k,var(k),lvl(k),diff_max_field,diff_max_field_rel
     1             ,rms,ndiff,imaxd,jmaxd

 101    format(' df_mx - fld #',i4,' ',a
     1  ,i6,' abs/rel/rms ',3e13.5,i8,' ndiff   max at',2i4)

        if(k .eq. khmax .and. nvar .gt. 1 
     1                  .and. n_levels .gt. 1)then
                write(6,*)
                write(6,*)' max diff for variable ',var(k)(1:3),' ='
     1                      ,diff_max_var
                nvar = nvar + 1
        endif

        if(ndiff_msg .gt. 0)then
            write(6,*)
            write(6,*)' warning: # of points differing '
     1                ,'wrt missing data = ',
     1                  ndiff_msg
        endif

        if(ndiff + ndiff_msg .gt. num_diff_field_thresh)then
          if(l_pass)then
            write(6,*)' overall criteria failure'
     1                           ,ndiff,ndiff_msg,num_diff_field_thresh
            l_pass = .false.     
          endif
        endif

        write(6,*)

        var_last = var(k)(1:3)

        enddo ! k


      if (istatus .ne. 1) write (6,*)'error in readlapsdata'

	   write(6,*)' overall file diff_max (',ext(1:3)
     1      ,') [abs/rel/#] = ',diff_max_file,diff_max_file_rel
     1               ,ndiff_file
           write(6,*)
           n_files = n_files + 1
           diff_max_all     = max(diff_max_all,diff_max_file)
           diff_max_all_rel = max(diff_max_all_rel,diff_max_file_rel)
           goto 5
      
999     continue

        if(n_files .gt. 1)then
          if(l_pass)then
             write(6,*)' max difference (all files)  [abs/rel/#] = '
     1              ,diff_max_all,diff_max_all_rel
     1              ,ndiff_all ! ,' passed'
             write(6,*)
             write(6,*)' passed'
          else
             write(6,*)' max difference (all files)  [abs/rel/#] = '
     1              ,diff_max_all,diff_max_all_rel
     1              ,ndiff_all ! ,' failed'
             write(6,*)
             write(6,*)' failed'
          endif
           write(6,*)
        else
          if(l_pass)then
            write(6,*)' passed'
            write(6,*)
          else
            write(6,*)' failed'
            write(6,*)
          endif
        endif
 
 
        if(idiff_msg_flag .eq. 1)then
            write(6,*)
            write(6,*)' warning: differences wrt missing data detected'
            write(6,*)
        endif

        stop

      end
      

      subroutine rlh(ext,nz_l,var,lvl,khmax,istatus)
!                     i   i    o   o    o      o
      
      character*(*) var(*)
      character*3 ext
      integer lvl(*)

      khmax = 0
      istatus = 1
      
      if    (ext .eq. 'lt1')then                                  ! temp
          call fill_3d_header('t3',var,lvl,nz_l,khmax)
          call fill_3d_header('ht',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'lsx')then
          khmax = 7
          var(1) = 't'
          var(2) = 'td'
          var(3) = 'u'
          var(4) = 'v'
          var(5) = 'p'
          var(6) = 'ps'
          var(7) = 'msl'
!         var(8) = 'vis'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0
          lvl(4) = 0
          lvl(5) = 0
          lvl(6) = 0
          lvl(7) = 0
!         lvl(8) = 0

      elseif(ext .eq. 'lst')then
          khmax = 8
          var(1) = 'pbe'
          var(2) = 'nbe'
          var(3) = 'li'
          var(4) = 'si'
          var(5) = 'tt'
          var(6) = 'lcl'
          var(7) = 'k'
          var(8) = 'wb0'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0
          lvl(4) = 0
          lvl(5) = 0
          lvl(6) = 0
          lvl(7) = 0
          lvl(8) = 0

      elseif(ext .eq. 'vrc')then
          khmax = 1
          var(1) = 'ref'
          lvl(1) = 0

      elseif(ext .eq. 'lw3')then                                  ! wind
          call fill_3d_header('u3',var,lvl,nz_l,khmax)
          call fill_3d_header('v3',var,lvl,nz_l,khmax)
          call fill_3d_header('om',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'lga')then                                  ! lga
          call fill_3d_header('ht',var,lvl,nz_l,khmax)
          call fill_3d_header('t3',var,lvl,nz_l,khmax)
          call fill_3d_header('sh',var,lvl,nz_l,khmax)
          call fill_3d_header('u3',var,lvl,nz_l,khmax)
          call fill_3d_header('v3',var,lvl,nz_l,khmax)
          call fill_3d_header('om',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'lgb')then
          khmax = 9  
          var(1) = 'usf'
          var(2) = 'vsf'
          var(3) = 'tsf'
          var(4) = 'tgd'
          var(5) = 'dsf'
          var(6) = 'slp'
          var(7) = 'psf'
          var(8) = 'rsf'
          var(9) = 'p'
!         var(10) = 'pcp'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0
          lvl(4) = 0
          lvl(5) = 0
          lvl(6) = 0
          lvl(7) = 0
          lvl(8) = 0
          lvl(9) = 0
!         lvl(10) = 0

      elseif(ext(1:2) .eq. 'v0')then                              ! v0x
          call fill_3d_header('ref',var,lvl,nz_l,khmax)
          call fill_3d_header('vel',var,lvl,nz_l,khmax)
          call fill_3d_header('nyq',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'vrz')then                                  
          call fill_3d_header('ref',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'lwm')then
          khmax = 2
          var(1) = 'su'
          var(2) = 'sv'
          lvl(1) = 0
          lvl(2) = 0

      elseif(ext .eq. 'liw')then
          khmax = 2
          var(1) = 'liw'
          var(2) = 'umf'
          lvl(1) = 0
          lvl(2) = 0

      elseif(ext .eq. 'lhe')then
          khmax = 3
          var(1) = 'lhe'
          var(2) = 'mu'
          var(3) = 'mv'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0

      elseif(ext .eq. 'lil')then
          khmax = 2
          var(1) = 'lil'
          var(2) = 'lic'
!         var(3) = 'cod'
!         var(4) = 'cla'
          lvl(1) = 0
          lvl(2) = 0
!         lvl(3) = 0
!         lvl(4) = 0

      elseif(ext .eq. 'lfr')then
          khmax = 4
          var(1) = 'vnt'
          var(2) = 'ham'
          var(3) = 'hah'
          var(4) = 'fwi'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0
          lvl(4) = 0

      elseif(ext .eq. 'lc3')then                                  ! cloud
          call fill_3d_header('lc3',var,lvl,42,khmax)

      elseif(ext .eq. 'lcp')then
          call fill_3d_header('lcp',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'lco')then
          call fill_3d_header('com',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'cty')then
          call fill_3d_header('cty',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'pty')then
          call fill_3d_header('pty',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'lwc')then
          call fill_3d_header('lwc',var,lvl,nz_l,khmax)
          call fill_3d_header('ice',var,lvl,nz_l,khmax)
          call fill_3d_header('pcn',var,lvl,nz_l,khmax)
          call fill_3d_header('rai',var,lvl,nz_l,khmax)
          call fill_3d_header('sno',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'lct')then                                
          khmax = 3
          var(1) = 'ptt'
          var(2) = 'pty'
          var(3) = 'sct'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0

      elseif(ext .eq. 'lcb')then                                
          khmax = 3
          var(1) = 'lcb'
          var(2) = 'lct'
          var(3) = 'cce'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0

      elseif(ext .eq. 'lcv')then                                
          khmax = 2
          var(1) = 'lcv'
          var(2) = 'csc'
          lvl(1) = 0
          lvl(2) = 0

      elseif(ext .eq. 'lmt')then                                
          khmax = 2
          var(1) = 'lmt'
          var(2) = 'llr'
          lvl(1) = 0
          lvl(2) = 0

      elseif(ext .eq. 'lmr')then                                
          khmax = 1
          var(1) = 'r'
          lvl(1) = 0

      elseif(ext .eq. 'lps')then
          call fill_3d_header('ref',var,lvl,nz_l,khmax)

      elseif(ext .eq. 'l1s')then                                  ! accum
          khmax = 4
          var(1) = 's01'
          var(2) = 'sto'
          var(3) = 'r01'
          var(4) = 'rto'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0
          lvl(4) = 0

      elseif(ext .eq. 'lvd')then                                  ! lvd
          khmax = 4
          var(1) = 'alb'
          var(2) = 'svs'
          var(3) = 'svn'
          var(4) = 's8a'
          lvl(1) = 0
          lvl(2) = 0
          lvl(3) = 0
          lvl(4) = 0

      else
          write(6,*)' ext is ',ext
          istatus = 0

      endif

      return
      end

      subroutine fill_3d_header(varin,var,lvl,nz_l,khmax)
!                                 i    o   o   i     o

      character*(*) varin,var(*)
      integer lvl(*)

      call get_r_missing_data(r_missing_data,istatus)
      call get_laps_dimensions(nk,istatus)

      do k = 1,nz_l
          khmax = khmax + 1
          var(khmax) = varin

          if(k .le. nk)then
              arg = pressure_of_level(k)
              lvl(khmax) = nint(pressure_of_level(k))/100
          endif
      enddo
 
      return
      end


