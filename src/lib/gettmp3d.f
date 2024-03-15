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
        subroutine get_temp_3d(i4time_needed,i4time_nearest,iflag,
     1                          imax,jmax,kmax,field_3d,istatus)

cdoc    returns 3-d temperature (or ht) field within a time window. has options
cdoc    for theta and/or balanced fields depending on the 'iflag' input.

!       steve albers            2000

!       iflag = 0 (potential temperature k)
!       iflag = 1 (temperature k)
!       iflag = 2 (height m)
!       iflag = 3 (balanced pot'l temperature k)
!       iflag = 4 (balanced temperature k)
!       iflag = 5 (balanced height m)

        logical ltest_vertical_grid

        character*150 directory
        character*31 ext

        character*125 comment_3d(kmax)
        character*10 units_3d(kmax)
        character*3 var_t
        integer lvl_3d(kmax)
        character*4 lvl_coord_3d(kmax)

        real field_3d(imax,jmax,kmax)

        character*255 c_filespec

        o_k(t_k,p_pa)   =   o( t_k-273.15 , p_pa/100. )  + 273.15
!       tda_k(t_k,p_pa) = tda( t_k-273.15 , p_pa/100. )  + 273.15

        write(6,*)
        write(6,*)' subroutine get_temp_3d: iflag=',iflag

        ext = 'lt1'

        if(iflag .lt. 3)then ! 'temp.exe' analysis
            call get_directory(ext,directory,len_dir)
        else
            call get_directory('balance',directory,lend)
            directory=directory(1:lend)//'lt1/'
        endif

        if(iflag .eq. 2 .or. iflag .eq. 5)then
            var_t = 'ht'
        else
            var_t = 't3'
        endif

        call get_3dgrid_dname(directory
     1           ,i4time_needed,1000000,i4time_nearest
     1           ,ext,var_t,units_3d
     1           ,comment_3d,imax,jmax,kmax,field_3d,istatus)       


!       test for bad temperatures
        if(var_t .eq. 't3')then
            if(field_3d(29,29,17) .lt. 173.)istatus = 0
        endif

        if(istatus .ne. 1)then
            write(6,*)' error reading in 3d lt1 field ',var_t
        endif

!       convert from t to theta
        do i = 1,imax
        do j = 1,jmax

        if(iflag .eq. 0 .or. iflag .eq. 3)then
            do k = 1,kmax
                theta = o_k(field_3d(i,j,k),zcoord_of_level(k))
                field_3d(i,j,k) = theta
            enddo ! k
        endif

        enddo ! j
        enddo ! i

        return

        end

