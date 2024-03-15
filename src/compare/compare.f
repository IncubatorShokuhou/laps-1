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

      include 'lapsparms.cmn'
c
       parameter (imax=125)
       parameter (jmax=105)
c      parameter (imax=nx_l_max)
c      parameter (jmax=ny_l_max)
      parameter (kmax=200)
      parameter (nlvl=200)
      
      integer	i4time,		!i4time of data
     1		kdim,		!k dimension of data array
     1		lvl(nlvl),        !level of each field (4 digit max)
     1		lvl_avail(nlvl),        !level of each field (4 digit max)
     1		flag,
     1		error(2),
     1		indx,
     1		i,j,k,start,
     1		istatus
c
      real	data1(imax,jmax,kmax)	!raw data to be written
      real	data2(imax,jmax,kmax)	!raw data to be written
c
      character*200	dir_in  !directory to be written to
      character*200	dir_out !directory to be written to
      character*31	ext		!file name ext (up to 31 chars)
      character*3	var(nlvl) 	        !3 letter id of each field
      character*3	laps_var_avail(nlvl) 	!3 letter id of each field
      character*19      var_avail(nlvl)
      character*4	lvl_coord(nlvl)	!vertical coordinate for each field
      character*10	units(nlvl)	!units of each field
      character*125	comment1(nlvl)	!comments for each field
      character*125	comment2(nlvl)	!comments for each field
      character*9	gtime
      character*91	file_name
      character*17	asctime
      
      character*4	version
      character*131	model 		!meteorological model in file
      character*131	origin		!location where file was created
      character*11    laps_dom_file   !name of domain file e.g. nest7grid
      logical         l_packed_data

      logical l_pass

      character*20 filename
      character*9 a9_time
      character*3 var_last


      var_last = '   '
      idiff_msg_flag = 0
      diff_max_all = 0.
      diff_max_all_rel = 0.
      n_files = 0
      ndiff_all = 0
      l_pass = .true.

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

      if(filename(1:6) .eq. 'static')then
          khmax = 4
          var(1) = 'lat'
          var(2) = 'lon'
          var(3) = 'avg'
          var(4) = 'ldf'
          ext = 'nest7grid'
          call rd_laps_static(dir_in,ext,imax,jmax,khmax,var,units
     1                     ,comment,data1,grid_spacing_m,istatus)

          if(istatus .ne. 1)then
              write(6,*)' bad status reading 1st static file'
              stop
          endif

          call rd_laps_static(dir_out,ext,imax,jmax,khmax,var,units
     1                     ,comment,data2,grid_spacing_m,istatus)

          if(istatus .ne. 1)then
              write(6,*)' bad status reading 2nd static file'
              stop
          endif

          thresh_write_pair = .01
          thresh_count_diff = .01

          ihmax = imax
          jhmax = jmax
 
      else
          n1 = index(filename,'.')
          n2 = index(filename,' ')

          a9_time = filename(1:9)
          ext = filename(n1+1:n2-1)

          call downcase(ext,ext)
          call cv_asc_i4time(a9_time,i4time)


          if(dir_in(1:3) .eq. 'end')then    
!             write(6,*)' end'
              goto 999
          endif

          do i = 1,nlvl
              var(i) = '   '
          enddo

          call read_laps_header(i4time,dir_in,ext,ihmax,jhmax,khmax,
     1                         laps_dom_file,asctime,version,
     1                         model,origin,var,lvl,num_variables,
     1                         var_avail,laps_var_avail,num_levels,
     1                         lvl_avail,lvl_coord,units,
     1                         comment1,l_packed_data,istatus)


 
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

          elseif(ext(1:3) .eq. 'lcp')then
              thresh_write_pair = .01
              thresh_count_diff = .01

          elseif(ext(1:3) .eq. 'lcv')then
              thresh_write_pair = .01
              thresh_count_diff = .01

          elseif(ext(1:3) .eq. 'lwc')then
              thresh_write_pair = .01
              thresh_count_diff = .01

          elseif(ext(1:3) .eq. 'lw3')then
              thresh_write_pair = .01
              thresh_count_diff = .1

          elseif(ext(1:3) .eq. 'lwm')then
              thresh_write_pair = .01
              thresh_count_diff = .1

          elseif(ext(1:3) .eq. 'lps')then
              thresh_write_pair = .01
              thresh_count_diff = .1

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
              num_diff_field_thresh = imax * jmax

          endif

!         adjust values for this extension (dependent on machine) of: 
!         thresh_count_diff and num_diff_field_thresh

          if(machine_factor .eq. 0)then
              thresh_count_diff = 0
              num_diff_field_thresh = 0
          endif

          if(var(1) .eq. 'rh')var(1) = 'lhe'

          write(6,*)' reading: ',dir_in,a9_time,'.',ext(1:3)

      
          call read_laps_data(i4time,dir_in,ext,ihmax,jhmax,khmax,
     1             khmax,var,lvl,lvl_coord,units,comment1,data1,
     1             istatus)

          if(istatus .ne. 1)then
            if(ext(1:3) .ne. 'lsx' .and. ext(1:3) .ne. 'liw')then
              if(l_pass)then
                  write(6,*)' read error: overall criteria failure'
                  l_pass = .false.     
              endif
            endif
          endif

          write(6,*)' reading: ',dir_out,a9_time,'.',ext(1:3)

      
          call read_laps_data(i4time,dir_out,ext,ihmax,jhmax,khmax,
     1             khmax,var,lvl,lvl_coord,units,comment2,data2,
     1             istatus)

          if(istatus .ne. 1)then
            if(ext(1:3) .ne. 'lsx' .and. ext(1:3) .ne. 'liw')then
              if(l_pass)then
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
                thresh_count_diff = .01 * machine_factor
            elseif(ext(1:3) .eq. 'lmt' .and. var(k)(1:3) .eq. 'lmt')then
                thresh_count_diff = 2.0 * machine_factor
            elseif(ext(1:3) .eq. 'lmt' .and. var(k)(1:3) .eq. 'llr')then
                thresh_count_diff = 2.0 * machine_factor
            elseif(ext(1:3) .eq. 'lvd' .and. var(k)(1:3) .eq. 'svn')then
                thresh_count_diff = 10.0 * machine_factor
            elseif(ext(1:3) .eq. 'lvd' .and. var(k)(1:3) .eq. 'alb')then
                thresh_count_diff = 0.1 * machine_factor
                thresh_write_pair = 0.1                                   !jsmart addition 11-1-96
!           elseif(ext(1:3) .eq. 'lw3' .and. var(k)(1:2) .ne. 'om')then
!               thresh_count_diff = .1 * machine_factor
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

	do i = 1,ihmax
           do j = 1,jhmax

c          if(
c!    1       data1(i,j,k) .le. r_min_normal()   .or.
c     1       data1(i,j,k) .ge. r_max_normal()   .or.         
c!    1       data2(i,j,k) .le. r_min_normal()   .or.         
c     1       data2(i,j,k) .ge. r_max_normal()            
c     1                                                    )then
              if(isnan(data1(i,j,k)).ne.0 .or. isnan(data2(i,j,k)).ne.0)
     +             then
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
21                  format(1x,3i5,2f12.6,f12.6)
22              endif
            endif
          endif ! nan test
        enddo ! j
        enddo ! i

        if(comment1(k) .ne. comment2(k))then
            write(6,*)comment1(k)(1:80)
            write(6,*)comment2(k)(1:80)
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

 	write(6,*)' df_mx - fld #',k,' ',var(k)
     1  ,lvl(k),' abs/rel/#',diff_max_field,diff_max_field_rel
     1  ,ndiff,imaxd,jmaxd

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
c           goto 5
      
999     continue

        if(n_files .gt. 1)then
          if(l_pass)then
             write(6,*)' max difference (all files)  [abs/rel/#] = '
     1              ,diff_max_all,diff_max_all_rel
     1              ,ndiff_all,' passed'
          else
             write(6,*)' max difference (all files)  [abs/rel/#] = '
     1              ,diff_max_all,diff_max_all_rel
     1              ,ndiff_all,' failed'
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
      
