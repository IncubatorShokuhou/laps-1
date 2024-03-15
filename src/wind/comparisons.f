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



        subroutine compare_wind (
     1                upass1,vpass1,cgrid,                                ! i
     1                istat_radar_vel,max_radars,grid_ra_vel,n_radars,    ! i
     1                nx_r,ny_r,ioffset,joffset,                          ! i
     1                rlat_radar,rlon_radar,rheight_radar,                ! i
     1                lat,lon,                                            ! i
     1                ni,nj,nk,r_missing_data,                            ! i
     1                obs_barnes_in,max_obs,ncnt_total_in,                ! i
     1                l_point_struct,l_withheld_only,                     ! i
     1                weight_pirep,weight_prof,weight_sfc,weight_cdw,     ! i
     1                grid_laps_u,grid_laps_v,grid_laps_wt,istatus)       ! i/o

c****************************************************************************
c
c  purpose: provide a single point out of lapswind_anal to call
c           diagnostic comparision routines.
c
c
c  inputs: upass1
c          vpass1
c          istat_radar_vel
c          grid_ra_vel
c          rlat_radar
c          rlon_radar
c          rheight_radar
c          n_radars
c
c  outputs: none
c
c*********************************************************************

c***************** declarations **************************************
        include 'barnesob.inc'
        type (barnesob) :: obs_barnes_in(max_obs)      
        type (barnesob) :: obs_barnes(max_obs)      

        integer istat_radar_vel
        integer l,n_radars,ni,nj,nk,max_radars
        integer ioffset(max_radars),joffset(max_radars)

        real rlat_radar(max_radars),rlon_radar(max_radars)
     1                     ,rheight_radar(max_radars)

        real lat(ni,nj),lon(ni,nj)

        real upass1(ni,nj,nk),vpass1(ni,nj,nk)
        real grid_ra_vel(nx_r,ny_r,nk,max_radars),r_missing_data
        real weight_pirep,weight_prof,weight_sfc,weight_cdw
        real grid_laps_u(ni,nj,nk),grid_laps_v(ni,nj,nk)
     1                                          ,grid_laps_wt(ni,nj,nk)

        integer max_obstypes
        parameter (max_obstypes=10)

        character*4 cgrid
        character*12 c_obstype_a(max_obstypes)
        logical l_parse, l_point_struct, l_withheld_only, l_compare_ob

c********************************************************************

        write(6,*)' subroutine compare_wind...',cgrid

!       copy obs structure into local structure depending on 'l_withheld_only'
        ncnt_total = 0
        do i = 1,ncnt_total_in
            if(l_withheld_only)then
                if(obs_barnes_in(i)%l_withhold)then
                    l_compare_ob = .true.
                    write(6,*)i,obs_barnes_in(i)%i,
     1                        obs_barnes_in(i)%j,obs_barnes_in(i)%type
                else
                    l_compare_ob = .false.
                endif
            else
                l_compare_ob = .true.
            endif

            if(l_compare_ob)then     
                ncnt_total = ncnt_total + 1
                obs_barnes(ncnt_total) = obs_barnes_in(i)
            endif
        enddo ! i      

        write(6,*)'l_withheld_only/ncnt_total_in/ncnt_total=',
     1             l_withheld_only,ncnt_total_in,ncnt_total

        if(ncnt_total .eq. 0)then
            write(6,*)' no obs detected, returning...'
            return
        endif

        call get_obstypes      (obs_barnes,max_obs,ncnt_total        ! i
     1                         ,c_obstype_a,max_obstypes,n_obstypes  ! i/o
     1                         ,istatus)                             ! o
        if(istatus .ne. 1)stop

!       n_obstypes = 4
!       c_obstype_a(1) = 'sfc '
!       c_obstype_a(2) = 'prof'
!       c_obstype_a(3) = 'pin '
!       c_obstype_a(4) = 'cdw '

        do i_obstype = 1,n_obstypes
            call s_len(c_obstype_a(i_obstype),len_obstype)
            write(6,*)

            if(c_obstype_a(i_obstype) .ne. 'radar')then

              if(l_withheld_only)then
                write(6,11)cgrid,c_obstype_a(i_obstype)(1:len_obstype)       
 11             format(1x,'  comparing ',a,' to ',a4
     1                   ,' obs (withheld - prior to qc)')    
              elseif(l_parse(cgrid,'fg'))then
                write(6,12)cgrid,c_obstype_a(i_obstype)(1:len_obstype)       
 12             format(1x,'  comparing ',a,' to ',a4
     1                   ,' obs (prior to qc)')    
              else
                write(6,13)cgrid,c_obstype_a(i_obstype)(1:len_obstype)
 13             format(1x,'  comparing ',a,' to ',a4
     1                   ,' obs (passing qc)')    
              endif

              call comp_grid_windobs(upass1,vpass1,ni,nj,nk
     1          ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_sfc
     1          ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1          ,cgrid,c_obstype_a(i_obstype),r_missing_data,rms)

            else
              write(6,*)' skip comparison for ',c_obstype_a(i_obstype)       

            endif

        enddo

        if(l_withheld_only)then
            write(6,*)' skip comparison for radial velocities'
            return
        endif

        do l = 1,n_radars
            write(6,*)
            write(6,*)'  comparing ',cgrid,' to radial velocities'       
     1                                                    ,' radar #',l
            if(istat_radar_vel .eq. 1)
     1        call comp_laps_vr(grid_ra_vel(1,1,1,l),upass1,vpass1
     1          ,ni,nj,nk,r_missing_data,cgrid,rms_fg_vr
     1          ,nx_r,ny_r,ioffset(l),joffset(l)                           ! i
     1          ,lat,lon,rlat_radar(l),rlon_radar(l),rheight_radar(l))

        enddo ! l

        return
        end



        subroutine comp_grid_windobs(u_3d,v_3d,ni,nj,nk
     1  ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_ob
     1  ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1  ,c_grid,c_obs,r_missing_data,rms)

        include 'barnesob.inc'
        type (barnesob) :: obs_barnes(max_obs)      

        real grid_laps_u(ni,nj,nk),grid_laps_v(ni,nj,nk)
        real grid_laps_wt(ni,nj,nk)
        real u_3d(ni,nj,nk),v_3d(ni,nj,nk) 

        character*4  c_grid
        character*12 c_ob_type,c_obs,c_obs_left,c_obs_right

        logical l_point_struct

        nobs = 0
        residualu = 0.
        residualv = 0.
        sumu = 0.
        sumv = 0.
        sumsp = 0.

        c_obs_left = c_obs
        call left_justify(c_obs_left)
        call s_len(c_obs_left,len_obstype)

        c_obs_right = c_obs
        call right_justify(c_obs_right)

        write(6,*)'comparing ',c_obs_left(1:len_obstype)
     1           ,' wind obs (passing qc) to ',c_grid,' grid'
        write(6,2)c_obs_right(4:12),c_grid
2       format(1x,'   i   j   k ',a,' ob (struct)    ob (array)'
     1                           ,a,' analysis        diff')       

        if(l_point_struct)then

!          checking the obstype is case insensitive
           call downcase(c_obs_left,c_obs_left)

           do iob = 1,ncnt_total
              c_ob_type = obs_barnes(iob)%type
              call downcase(c_ob_type,c_ob_type)
              call left_justify(c_ob_type)

              if(c_ob_type .eq. c_obs_left)then
                  nobs = nobs + 1
                  il = obs_barnes(iob)%i
                  jl = obs_barnes(iob)%j
                  k = obs_barnes(iob)%k
                  diffu = obs_barnes(iob)%valuef(1) - u_3d(il,jl,k) 
                  diffv = obs_barnes(iob)%valuef(2) - v_3d(il,jl,k) 

                  sumu = sumu + diffu
                  sumv = sumv + diffv

                  residualu = residualu + diffu ** 2
                  residualv = residualv + diffv ** 2

!                 calculate speed difference (ob - grid)
                  call uv_to_disp(obs_barnes(iob)%valuef(1),
     1                            obs_barnes(iob)%valuef(2),
     1                            di_dum,
     1                            grid_laps_sp)

                  call uv_to_disp(u_3d(il,jl,k),
     1                            v_3d(il,jl,k),
     1                            di_dum,
     1                            sp_3d)

                  sumsp = sumsp + (grid_laps_sp - sp_3d)
                  if(nobs .le. 200 .or. nobs .eq. (nobs/10)*10)then
                      write(6,101)il,jl,k
     1                ,obs_barnes(iob)%valuef(1)
     1                ,obs_barnes(iob)%valuef(2)       
     1                ,grid_laps_u(il,jl,k),grid_laps_v(il,jl,k)
     1                ,u_3d(il,jl,k),v_3d(il,jl,k)
     1                ,diffu,diffv
101                   format(1x,3i4,4(2x,2f7.1))
                  endif

              endif ! obstype match

           enddo ! iob

        else
           do jl = 1,nj

           do il = 1,ni

              do k = 1,nk

                if(grid_laps_wt(il,jl,k) .eq. weight_ob
     1             .and. u_3d(il,jl,k) .ne. r_missing_data )then
                  nobs = nobs + 1

!                 calculate u/v difference (ob - grid)
                  diffu = grid_laps_u(il,jl,k) - u_3d(il,jl,k) 
                  diffv = grid_laps_v(il,jl,k) - v_3d(il,jl,k) 

                  sumu = sumu + diffu
                  sumv = sumv + diffv

                  residualu = residualu + diffu ** 2
                  residualv = residualv + diffv ** 2

!                 calculate speed difference (ob - grid)
                  call uv_to_disp(grid_laps_u(il,jl,k),
     1                            grid_laps_v(il,jl,k),
     1                            di_dum,
     1                            grid_laps_sp)

                  call uv_to_disp(u_3d(il,jl,k),
     1                            v_3d(il,jl,k),
     1                            di_dum,
     1                            sp_3d)

                  sumsp = sumsp + (grid_laps_sp - sp_3d)

                  if(nobs .le. 200 .or. nobs .eq. (nobs/10)*10)then
                      write(6,101)il,jl,k
     1                  ,grid_laps_u(il,jl,k),grid_laps_v(il,jl,k)
     1                  ,u_3d(il,jl,k),v_3d(il,jl,k)
     1                  ,diffu,diffv
                  endif

                endif

              enddo ! k

           enddo ! il
           enddo ! jl

        endif ! l_point_struct

        if(nobs .gt. 0)then
            rmsu = sqrt(residualu/nobs)
            rmsv = sqrt(residualv/nobs)
            biasu = sumu / float(nobs)
            biasv = sumv / float(nobs)
            biassp = sumsp / float(nobs)

        else
            rmsu = 0.
            rmsv = 0.
            biasu = 0.
            biasv = 0.
            biassp = 0.

        endif

        rms  = sqrt(rmsu**2 + rmsv**2)

        if(nobs .gt. 0)then
            call upcase(c_obs_left,c_obs_left)
            write(6,102)c_obs_left(1:8),c_grid,nobs,biasu,biasv,biassp
     1                 ,rmsu,rmsv,rms      
102         format(' bias/rms between '
     1        ,a,' & ',a,' (n,biasu,biasv,biassp,rmsu,rmsv,rms) = '
     1        ,i5,6f6.1)
        endif

        return

        end


        subroutine get_obstypes(obs_barnes,max_obs,ncnt_total
     1                         ,c_obstype_a,max_obstypes,n_obstypes
     1                         ,istatus)

!       obtain list of obstypes contained within the obs data structure

        include 'barnesob.inc'
        type (barnesob) :: obs_barnes(max_obs)      

        character*12 c_obstype_a(max_obstypes)

        logical l_match_found

        n_obstypes = 0

        if(ncnt_total .eq. 0)then
            return
        endif        

        write(6,*)' subroutine get_obstypes, number of obs...'
     1            ,ncnt_total

        i = 1
        n_obstypes = 1
        write(6,*)n_obstypes,i,obs_barnes(i)%type
        c_obstype_a(i) = obs_barnes(i)%type

        if(ncnt_total .ge. 2)then
            do i = 1,ncnt_total
                l_match_found = .false.
                do j = 1,n_obstypes
                    if(obs_barnes(i)%type .eq. c_obstype_a(j))then
                        l_match_found = .true.
                    endif
                enddo ! j
                if(.not. l_match_found)then
                    n_obstypes = n_obstypes + 1
                    write(6,*)n_obstypes,i,obs_barnes(i)%type
                    if(n_obstypes .gt. max_obstypes)then
                        write(6,*)' error: too many obstypes'
                        istatus = 0
                        return
                    endif
                    c_obstype_a(n_obstypes) = obs_barnes(i)%type
                endif
            enddo ! i
        endif

        write(6,*)

!       determine how many obs have each obs type
        do itype = 1,n_obstypes
            n_match = 0
            do i = 1,ncnt_total
                if(obs_barnes(i)%type .eq. c_obstype_a(itype))then
                    n_match = n_match + 1
                endif
            enddo ! i
            write(6,*)n_match,' obs of type ',c_obstype_a(itype)
        enddo ! itype           

        istatus = 1

        return
        end 
