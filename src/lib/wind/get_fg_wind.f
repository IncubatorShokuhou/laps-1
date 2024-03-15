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

        subroutine get_fg_wind_new(
     1          i4time_lapswind,ilaps_cycle_time               ! input
     1          ,nx_l,ny_l,nz_l                                ! input
     1          ,ntmin,ntmax                                   ! input
     1          ,u_mdl_bkg_4d,v_mdl_bkg_4d                     ! output
     1          ,u_laps_fg,v_laps_fg                           ! output
     1          ,istatus                                       ! output
     1                                  )

!       1998    steve albers - fsl

        real u_mdl_bkg_4d(nx_l,ny_l,nz_l,ntmin:ntmax)
        real v_mdl_bkg_4d(nx_l,ny_l,nz_l,ntmin:ntmax)
        real u_laps_fg(nx_l,ny_l,nz_l),v_laps_fg(nx_l,ny_l,nz_l)

        character*3 var_2d

        write(6,*)
        write(6,*)' obtain wind tendency from model first guesses:'

        var_2d = 'u3'

        call get_fg_var(
     1          i4time_lapswind,ilaps_cycle_time               ! input
     1          ,nx_l,ny_l,nz_l                                ! input
     1          ,ntmin,ntmax                                   ! input
     1          ,var_2d                                        ! input
     1          ,u_mdl_bkg_4d                                  ! output
     1          ,istatus                                       ! output
     1                                  )

        if(istatus .ne. 1)then
            write(6,*)' error processing model first guess for ',var_2d       
            return
        endif

        var_2d = 'v3'

        call get_fg_var(
     1          i4time_lapswind,ilaps_cycle_time               ! input
     1          ,nx_l,ny_l,nz_l                                ! input
     1          ,ntmin,ntmax                                   ! input
     1          ,var_2d                                        ! input
     1          ,v_mdl_bkg_4d                                  ! output
     1          ,istatus                                       ! output
     1                                  )

        if(istatus .ne. 1)then
            write(6,*)' error processing model first guess for ',var_2d       
            return
        endif

        write(6,*)' using model winds for u/v_laps_fg'       

        call move_3d(u_mdl_bkg_4d(1,1,1,0),u_laps_fg,nx_l,ny_l,nz_l)       
        call move_3d(v_mdl_bkg_4d(1,1,1,0),v_laps_fg,nx_l,ny_l,nz_l)

        write(6,*)'u_laps_fg(nx_l/2+1,ny_l/2+1,1) = '
     1            ,u_laps_fg(nx_l/2+1,ny_l/2+1,1)
        write(6,*)'v_laps_fg(nx_l/2+1,ny_l/2+1,1) = '
     1            ,v_laps_fg(nx_l/2+1,ny_l/2+1,1)

        return
        end


        subroutine get_fg_var(
     1          i4time_lapswind,ilaps_cycle_time               ! input
     1          ,nx_l,ny_l,nz_l                                ! input
     1          ,ntmin,ntmax                                   ! input
     1          ,var_2d                                        ! input
     1          ,var_mdl_bkg_4d                                ! output
     1          ,istatus                                       ! output
     1                                  )

!       1999    steve albers - fsl

        dimension var_mdl_bkg_4d(nx_l,ny_l,nz_l,ntmin:ntmax)

        character*3 var_2d

        write(6,*)
        write(6,*)' obtain field tendency from model first guesses:'

!  ***  read in model data   ********************************************
        do nt = ntmin,ntmax

!         if(nt .ne. 1)then
          if(.true.)then

            write(6,*)' reading var_mdl_bkg_4d, var/nt= ',var_2d,nt
            call get_modelfg_3d(i4time_lapswind+ilaps_cycle_time*nt      
     1          ,var_2d,nx_l,ny_l,nz_l
     1          ,var_mdl_bkg_4d(1,1,1,nt),istatus)
            if(istatus .ne. 1)then
                write(6,*)' aborting from get_fg_var'
     1                   ,' - error reading model data ',var_2d,nt
                return
            endif

            call qc_field_3d('u3',var_mdl_bkg_4d(1,1,1,nt)
     1                     ,nx_l,ny_l,nz_l,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error: qc flag in background field'
                return
            endif

          else ! nt = 1

            write(6,*)' bogusing var_mdl_bkg_4d, var/nt= ',var_2d,nt

!           note that we can eventually read in the data for nt=1 directly
!           after all obs data is processed "the right way" wrt time
!           the current strategy allows for "null" testing of doing improved
!           time interpolation of obs.
            do k = 1,nz_l
            do j = 1,ny_l
            do i = 1,nx_l
                var_mdl_bkg_4d(i,j,k,1) = 2. * var_mdl_bkg_4d(i,j,k,0) 
     1                                -        var_mdl_bkg_4d(i,j,k,-1)        
            enddo ! i
            enddo ! j
            enddo ! k

          endif ! nt .ne. 1

          write(6,*)'var_mdl_bkg_4d(nx_l/2+1,ny_l/2+1,1,nt) = '
     1              ,var_mdl_bkg_4d(nx_l/2+1,ny_l/2+1,1,nt)

        enddo ! nt

        i4_elapsed = ishow_timer()

        return
        end

