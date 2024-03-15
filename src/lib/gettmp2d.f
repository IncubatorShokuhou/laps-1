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
        subroutine get_temp_2d(i4time_needed,i4time_tol,i4time_nearest
     1                          ,ilevel_mb,imax,jmax,temp_2d,istatus)

!       steve albers            1990

        logical ltest_vertical_grid

        character*150 directory
        character*31 ext

        character*125 comment_2d
        character*10 units_2d
        character*3 var_t
        integer lvl_2d
        character*4 lvl_coord_2d

        real temp_2d(imax,jmax)

        character*255 c_filespec

        write(6,*)
        write(6,*)' getting temperature analysis at ',ilevel_mb

        ext = 'lt1'
        call get_directory(ext,directory,len_dir)
        c_filespec = directory(1:len_dir)//'*.'//ext(1:3)

        call get_file_time(c_filespec
     1  ,i4time_needed,i4time_nearest)

        if(abs(i4time_needed - i4time_nearest) .gt. i4time_tol)then
            write(6,*)' no file available within requested time window'
            istatus = 0
            return
        endif

        units_2d  = 'k'
        if(ltest_vertical_grid('height'))then
            lvl_2d = zcoord_of_level(k)/10
            lvl_coord_2d = 'msl'
        elseif(ltest_vertical_grid('pressure'))then
            lvl_2d = ilevel_mb
            lvl_coord_2d = 'mb'
        else
            write(6,*)' error, vertical grid not supported,'
     1               ,' this routine supports pressure or height'
            istatus = 0
            return
        endif

        var_t = 't3'  ! newvar = 't3', oldvar = 't'

!       read in 2d t array
        call read_laps_data(i4time_nearest,directory,
     1          ext,imax,jmax,1,1,
     1          var_t,lvl_2d,lvl_coord_2d,units_2d,comment_2d,
     1          temp_2d,istatus)

!       test for bad temperatures
        if(temp_2d(29,29) .lt. 173.)istatus = 0

        return

        end

