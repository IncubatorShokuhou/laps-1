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
        subroutine insert_pireps(i4time,cld_hts                          ! i
     1        ,default_clear_cover                                       ! i
     1        ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd          ! i/o
     1        ,lat,lon,ni,nj,nk,ix_low,ix_high,iy_low,iy_high,max_pireps ! i
     1        ,istatus)                                                  ! o

        real lat(ni,nj),lon(ni,nj)

        real cld_hts(nk)

!       arrays for cloud soundings
        real cld_snd(max_cld_snd,nk)
        real wt_snd(max_cld_snd,nk)
        integer i_snd(max_cld_snd)
        integer j_snd(max_cld_snd)

        character*9  asc9_tim_pirep,asc9_tim_rcvd
        character*80 string

        character*31    ext

        logical l_good_pirep

!       1 feb 1996     steve albers        remove spread( calls

        max_layers = 3  ! defined in pin file and netcdf raw data

        lun = 11
        ext = 'pin'
        call open_lapsprd_file_read(lun,i4time,ext,istatus)
        if(istatus .ne. 1)go to 998

        n_cloud = 0
        num_pirep = 0
        num_good_pirep = 0

5       read(11,101,end=900,err=50)string(1:6)
101     format(a6)

50      if(string(2:5) .eq. 'time')then
!           a9time = string(30:39)
            read(11,151)asc9_tim_pirep,asc9_tim_rcvd
151         format(1x,a9,2x,a9)
            write(6,*)
            write(6,151)asc9_tim_pirep,asc9_tim_rcvd
        endif

        if(string(2:4) .eq. 'lat')then
            read(11,201)rlat,rlon,ralt
201         format(2(f8.3,2x), f6.0,2i5)
            call latlon_to_rlapsgrid(rlat,rlon,lat,lon,ni,nj,ri,rj
     1                              ,istatus)
!           if(istatus.ne.1)return
            ilaps = nint(ri)
            jlaps = nint(rj)
            write(6,201)rlat,rlon,ralt,ilaps,jlaps
        endif

        if(string(2:5) .eq. 'clou')then

            if(                   ilaps .ge. ix_low   ! in bounds
     1                      .and. ilaps .le. ix_high  ! in bounds
     1                      .and. jlaps .ge. iy_low   ! in bounds
     1                      .and. jlaps .le. iy_high  ! in bounds
     1                                          )then ! in bounds

              n_cld_snd = n_cld_snd + 1  ! create a new souding
              num_pirep = num_pirep + 1

              if(num_pirep .gt. max_pireps)then
                  write(6,*)' insert_pireps: error, too many pireps,'
     1                     ,' check n_pirep parameter'
                  istatus = 0
                  return
              endif

              l_good_pirep = .false.

!             pin file input heights are feet msl
              do i = 1,max_layers
                read(11,203,err=500)cbase_ft,ctop_ft,icover
203             format (12x,2f8.0,i5)
                write(6,*,err=500)' got a cloud report',
     1                  nint(cbase_ft),nint(ctop_ft),icover

                if(icover .gt. 0)then ! good report (icover is in eighths)

                  if(cbase_ft .ge. 0. .and. ctop_ft .ge. 0.)then ! fill layer

                    cbase_m = cbase_ft / 3.281
                    ctop_m  = ctop_ft  / 3.281

                    cbuf_low = cbase_m - 500.
                    cbuf_high = ctop_m + 500.

                    cover_rpt = float(icover) / 8.0 ! 1.0

                    n_cloud = n_cloud + 1
                    write(6,*)' good layer report:      '
     1                       ,nint(cbase_m),nint(ctop_m),cover_rpt

                    do k=1,nk
                        cover = cover_rpt

!                       fill in cloud layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                     cld_hts(k) .lt. ctop_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = default_clear_cover

!                       fill in clear buffer under cloud layer
                        if(cld_hts(k) .ge. cbuf_low .and.
     1                     cld_hts(k) .lt. cbase_m                 )then       
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

!                       fill in clear buffer above cloud layer
                        if(cld_hts(k) .ge. ctop_m .and.
     1                     cld_hts(k) .lt. cbuf_high               )then       
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                    enddo ! k cld_hts

                  elseif(cbase_ft .gt. 0)then
                    cbase_m = cbase_ft / 3.281
                    ctop_m  = cbase_m + 1000.

                    cbuf_low = cbase_m - 500.

                    cover_rpt = float(icover) / 8.0 ! 1.0

                    write(6,*)' only a base reported:'
     1                       ,nint(cbase_m),nint(ctop_m),cover_rpt

                    n_cloud = n_cloud + 1

                    do k=1,nk
                        cover = cover_rpt

!                       fill in cloud layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                     cld_hts(k) .lt. ctop_m                  )then       
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = default_clear_cover

!                       fill in clear buffer under cloud layer
                        if(cld_hts(k) .ge. cbuf_low .and.
     1                     cld_hts(k) .lt. cbase_m                 )then       
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                    enddo ! k cld_hts

                  elseif(ctop_ft .gt. 0)then
                    ctop_m  = ctop_ft  / 3.281
                    cbase_m = ctop_m - 500.

                    cbuf_high = ctop_m + 500.

                    cover_rpt = float(icover) / 8.0 ! 1.0

                    write(6,*)' only a top reported:'
     1                       ,nint(cbase_m),nint(ctop_m),cover_rpt

                    n_cloud = n_cloud + 1

                    do k=1,nk
                        cover = cover_rpt

!                       fill in cloud layer
                        if(cld_hts(k) .ge. cbase_m .and.
     1                     cld_hts(k) .lt. ctop_m                  )then
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                        cover = default_clear_cover

!                       fill in clear buffer above cloud layer
                        if(cld_hts(k) .ge. ctop_m .and.
     1                     cld_hts(k) .lt. cbuf_high               )then      
                            call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                                  ,n_cld_snd,max_cld_snd
     1                                  ,nk,ilaps,jlaps,k,cover,1.)
                            write(6,*)' fill in k = ',k,cover,cld_hts(k)
                            l_good_pirep = .true.
                        endif

                    enddo ! k cld_hts

                  endif ! base and/or top reported

                endif ! cover > 0

              enddo ! i cloud layer

              if(l_good_pirep)num_good_pirep = num_good_pirep + 1

            else ! out of bounds
!             skip over cloud report
              do i = 1,max_layers
                read(11,203,err=500)cbase_ft,ctop_ft,icover
              enddo

              write(6,*)' out of domain perimeter'

            endif ! in bounds

        endif ! cloud report string

!       if(string(2:5) .eq. 'sky ')then
!           read(11,204,err=500)isky_cover
!204         format (40x,i4)
!            write(6,*)' sky cover = ',isky_cover
!        endif

500     goto5


900     write(6,*)
        write(6,*)' completed insertion of pireps'
        write(6,*)' num pireps/num good pireps/cloud layers = '
     1                  ,num_pirep,num_good_pirep,n_cloud
        close(11)

        istatus = 1
        return

998     write(6,*)' warning, could not find the pirep file'
        istatus = 1
        return

        end

