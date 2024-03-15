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



        subroutine compare_temp (
     1                tpass1,cgrid,
     1                ni,nj,nk,r_missing_data,                            ! i
     1                obs_barnes_in,max_obs,ncnt_total_in,                ! i
     1                l_withheld_only,                                    ! i
     1                weight_sfc,                                         ! i  
     1                istatus)                                            ! i/o

c****************************************************************************
c
c  purpose: provide a single point out of lapstemp_anal to call
c           diagnostic comparision routines.
c
c
c  inputs: tpass1
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

        real lat(ni,nj),lon(ni,nj)

        real tpass1(ni,nj,nk)
        real r_missing_data
        real weight_sfc

        integer max_obstypes
        parameter (max_obstypes=11)

        character*4 cgrid
        character*12 c_obstype_a(max_obstypes)
        logical l_parse, l_point_struct, l_withheld_only, l_compare_ob

c********************************************************************

        write(6,*)' subroutine compare_temp...',cgrid

        l_point_struct = .true.

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

        call get_temp_obstypes (obs_barnes,max_obs,ncnt_total        ! i
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

            call comp_grid_tempobs(tpass1,ni,nj,nk
     1          ,weight_sfc
     1          ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1          ,cgrid,c_obstype_a(i_obstype),r_missing_data,rms)

        enddo

        return
        end



        subroutine comp_grid_tempobs(t_3d,ni,nj,nk
     1  ,weight_ob
     1  ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1  ,c_grid,c_obs,r_missing_data,rms)

        include 'barnesob.inc'
        type (barnesob) obs_barnes(max_obs)      

        real t_3d(ni,nj,nk)

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
     1           ,' temp obs (passing qc) to ',c_grid,' grid'
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
                  difft = obs_barnes(iob)%value(1) - t_3d(il,jl,k) 

                  sumu = sumu + difft

                  residualu = residualu + difft ** 2

                  if(nobs .le. 200 .or. nobs .eq. (nobs/10)*10)then
                      write(6,101)il,jl,k
     1                ,obs_barnes(iob)%value(1)
     1                ,t_3d(il,jl,k)
     1                ,difft
101                   format(1x,3i4,3(2x,2f7.1))
                  endif

              endif ! obstype match

           enddo ! iob

        endif ! l_point_struct

        if(nobs .gt. 0)then
            rmsu = sqrt(residualu/nobs)
            bias_t = sumu / float(nobs)

        else
            rmsu = 0.
            bias_t = 0.

        endif

        rms  = sqrt(rmsu**2)

        if(nobs .gt. 0)then
            call upcase(c_obs_left,c_obs_left)
            write(6,102)c_obs_left(1:8),c_grid,nobs,bias_t,rmsu
102         format(' bias/rms between '
     1        ,a,' & ',a,' (n,bias_t,rms) = '
     1        ,i6,6f5.1)
        endif

        return

        end


        subroutine get_temp_obstypes(obs_barnes,max_obs,ncnt_total
     1                              ,c_obstype_a,max_obstypes,n_obstypes
     1                              ,istatus)

        include 'barnesob.inc'
        type (barnesob) obs_barnes(max_obs)      

        character*12 c_obstype_a(max_obstypes)

        logical l_match_found

        n_obstypes = 0

        if(ncnt_total .eq. 0)then
            return
        endif        

        write(6,*)' subroutine get_temp_obstypes, obstypes found...'

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

        istatus = 1

        return
        end 
