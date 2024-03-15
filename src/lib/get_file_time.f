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

        subroutine get_file_time(c_filespec
     1          ,i4time_needed_in,i4time_nearest)

!       steve albers            1990
cdoc    given a directory path with a wildcard and a desired time, find the best
cdoc    matching file i4time in the given directory matching that wildcard.
cdoc    if there are very many files in the directory, you will need to
cdoc    use a directory name only, the ls command can't handle it otherwise.
!       john smart              1996
c       modified the routines within to handle filenames of the type
c       yyyymmdd_hhmm in addition to those of the form yyjjjmmhh.
c       the 13 character filenames are specific to the wfo platform.
c
!       john smart              1996
c       incorporated subroutine filter_non_numeric_fnames. this eliminates
c       c_fnames that are other that numeric so that i4time_fname_lp works.
c
!       steve albers            1998            slightly more generic
       
        integer max_files
        parameter (max_files = 20000)

        character*13 asc_time,asc_tim_nearest,asc13_tim_needed
        character*13 cvt_i4time_wfo_fname13
        character*9  asc9_tim_needed
        character*20 c20_type, a20_time
        character*(*) c_filespec                                 ! input
        character c_fnames(max_files)*150

        common /laps_diag/ no_laps_diag

        if(i4time_needed_in .eq. 0)then
            i4time_needed = i4time_now_gg()
        else
            i4time_needed = i4time_needed_in
        endif

!       call downcase(c_filespec,c_filespec)

        call    get_file_names(c_filespec,
     1                   i_nbr_files_ret,
     1                   c_fnames,
     1                   max_files,
     1                   i_status )

        call    filter_non_numeric_fnames(c_fnames,
     1                   i_nbr_files_ret,
     1                   i_nbr_files_out,
     1                   max_files,
     1                   istatus)

        min_diff = 1999999999

        if(i_nbr_files_out .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
            do i=1,i_nbr_files_out
                call get_filetime_type(c_fnames(i),c20_type,leni,lent)       
                call i4time_fname_lp(c_fnames(i),i4time_file,istatus)       

                i4_diff = abs(i4time_file - i4time_needed)

                if(i4_diff .lt. min_diff)then
                    min_diff = i4_diff
                    i_file_nearest = i
                    i4time_nearest = i4time_file
                    a20_time = c_fnames(i)(leni+1:leni+lent)
                endif

            enddo ! i
        endif ! i_nbr_files_out > 0

        if(i_nbr_files_out .gt. 0)then
            if(no_laps_diag .eq. 0)then
                if(c20_type .ne. 'yyyymmdd_hhmm')then
                    call make_fnam_lp(i4time_needed,asc9_tim_needed
     1                                            ,i_status)
                    asc13_tim_needed = asc9_tim_needed
                else
                    asc13_tim_needed = 
     1                         cvt_i4time_wfo_fname13(i4time_needed)       

                endif

                write(6,*)'    file time (needed/nearest) = '
     1                  ,asc13_tim_needed,' / ',a20_time(1:lent)

            endif

        else
            call s_len(c_filespec,len_spec)
            write(6,*)'  no files available - ',c_filespec(1:len_spec)
            i4time_nearest = 0

        endif ! i_nbr_files_out > 0

        return
        end


        subroutine get_latest_file_time(c_filespec,i4time_latest)

!       steve albers            1995
cdoc    given a directory path with a wildcard  find the latest
cdoc    matching file i4time in the given directory matching that wildcard.
cdoc    if there are very many files in the directory, you will need to
cdoc    use a directory name only, the ls command can't handle it otherwise.

!       steve albers            1998
c       more generic

        integer max_files
        parameter (max_files = 20000)

        character*9 asc_tim_latest
        character*(*) c_filespec                                   ! input
        character c_fnames(max_files)*150

        common /laps_diag/ no_laps_diag

!       call downcase(c_filespec,c_filespec)

        call    get_file_names(c_filespec,
     1                   i_nbr_files_ret,
     1                   c_fnames,
     1                   max_files,
     1                   i_status )

        call    filter_non_numeric_fnames(c_fnames,
     1                   i_nbr_files_ret,
     1                   i_nbr_files_out,
     1                   max_files,
     1                   istatus)

        if(i_nbr_files_out .gt. 0)then
            i = i_nbr_files_out
            call i4time_fname_lp(c_fnames(i),i4time_file,istatus)
            call make_fnam_lp(i4time_file,asc_tim_latest,istatus)

            i_file_latest = i
            i4time_latest = i4time_file
        endif


        if(i_nbr_files_ret .gt. 0)then
            if(no_laps_diag .eq. 0)then
                write(6,*)'    file time (latest) = ',asc_tim_latest
            endif
        else
            call s_len(c_filespec,len_spec)
            write(6,*)'  no files available - ',c_filespec(1:len_spec)       
            i4time_latest = 0
        endif

        return
        end

        subroutine get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_out,istatus)

!       steve albers            1996
cdoc    given a directory path with a wildcard find all files and
cdoc    return a list of the filenames and associated i4times
cdoc    if there are very many files in the directory, you will need to
cdoc    use a directory name only, the ls command can't handle it otherwise.

!       steve albers            1998
c       more generic

        integer max_files

        character*(*) c_filespec                           ! input
        character*(*) c_fnames(max_files)
        integer i4times(max_files)

        common /laps_diag/ no_laps_diag

        call    get_file_names(c_filespec,
     1                   i_nbr_files_ret,
     1                   c_fnames,
     1                   max_files,
     1                   istatus )

        call    filter_non_numeric_fnames(c_fnames,
     1                   i_nbr_files_ret,
     1                   i_nbr_files_out,
     1                   max_files,
     1                   istatus)

        if(i_nbr_files_out .gt. 0)then
            do i = 1,i_nbr_files_out
                call i4time_fname_lp(c_fnames(i),i4times(i),istatus)
            enddo ! i
        endif

        write(6,*)' get_file_times - # files = ',i_nbr_files_out

        return
        end
