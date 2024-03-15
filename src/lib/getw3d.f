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
        subroutine get_w_3d(i4time,imax,jmax,kmax,w_3d,ext,istatus)

!       returns vertical wind component from a 3d gridded laps dataset

!       steve albers            1990

        logical ltest_vertical_grid

        character*150 directory
        character*31 ext     ! input extension of file (normally 3 characters)

        character*125 comment_3d(kmax)
        character*10 units_3d(kmax)
        character*3 var_w(kmax)
        integer lvl_3d(kmax)
        character*4 lvl_coord_3d(kmax)

        logical l_convert

        common /lapsplot_omega/l_convert

        real w_3d(imax,jmax,kmax)

        call get_directory(ext,directory,len_dir)

        do k = 1,kmax
            units_3d(k)   = 'pa/s'
            comment_3d(k) = '3dwind'
            if(ltest_vertical_grid('height'))then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'msl'
            elseif(ltest_vertical_grid('pressure'))then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'hpa'
            else
                write(6,*)' error, vertical grid not supported,'
     1                   ,' this routine supports pressure or height'
                istatus = 0
                return
            endif

            if(ext(1:3) .eq. 'lco')then
                var_w(k) = 'com' ! newvar = 'com', oldvar = 'om'
            else
                var_w(k) = 'om'
            endif

        enddo ! k


!       read in 3d w array
        call read_laps_data(i4time,directory,
     1          ext,imax,jmax,kmax,kmax,
     1          var_w,lvl_3d,lvl_coord_3d,units_3d,comment_3d,
     1          w_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' sorry, file has not yet been generated this hour
     1'
        else
            write(6,*)' 3d - laps w analysis successfully read in'
        endif


        return
        end
