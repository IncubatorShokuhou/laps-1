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

        subroutine get_clouds_3dgrid(i4time_needed,i4time_nearest
     1          ,imax,jmax,kcloudin,ext ! ,b_cldcv
     1                          ,clouds_3d,cld_hts,cld_pres,istatus)

!       steve albers            1990
!       steve albers            1991     (b_cldcv now is an input with adj dim)
!       steve albers     30 oct 1991     cld_pres now has adjustable dims
!       steve albers      6 nov 1991     ext is passed in
!	linda wharton    26 oct 1998	 removed units_2d and comment_2d var not used

!       this is a newer version of 'get_clouds_3d' that reads the heights and
!       approximate pressures directly from the cloud grid

!       argument list
!       i4time_needed      input   desired i4time
!       i4time_nearest     output  actual i4time of nearest file
!       imax,jmax          input   laps horizontal grid dimensions
!       kcloudin           input   number of cloud vertical levels (dense grid)
!                                  this must be equal to value in laps_cloud.inc
!                                  if file is before i4time = 969476400
        character*31 ext ! input   file extension (normally 'lc3')
!       byte b_cldcv(imax,jmax,kcloudin)    ! output   clouds in bytes
        real clouds_3d(imax,jmax,kcloudin)! output 3d array of cloud cover (0-1.2)
        real cld_pres(kcloudin) ! output pressures of the vertical grid (approx)
                                  ! in pascals

        include 'laps_cloud.inc'

!       these two arrays deal with old cloud data files only
        real cld_hts_old(kcloud)

        data cld_hts_old/1200.,1300.,1400.,1500.,1600.,1700.,1800.,
     11900.,2000.,2200.,2400.,2500.,2600.,2800.,3000.,3200.,
     23400.,3600.,3800.,4000.,4500.,5000.,5500.,6000.,6500.,
     37000.,7500.,8000.,8500.,9000.,9500.,10000.,11000.,12000.,
     413000.,14000.,15000.,16000.,17000.,18000.,19000.,20000./


        real cld_hts_new(kcloud)

        data cld_hts_new/1200.,1300.,1400.,1500.,1600.,1700.,1800.,
     11900.,2000.,2100.,2200.,2400.,2600.,2800.,3000.,3200.,
     23400.,3600.,3800.,4000.,4500.,5000.,5500.,6000.,6500.,
     37000.,7500.,8000.,8500.,9000.,9500.,10000.,11000.,12000.,
     413000.,14000.,15000.,16000.,17000.,18000.,19000.,20000./

        character*255 c_filespec
        character*13 filename13
        character*150 directory

!       stuff for readlapsdata
        character*125 comment_3d(kcloud)
        character*10 units_3d(kcloud)

        character*3 var_3d(kcloud),var_2d
        data var_2d/'lc3'/

        integer lvl_3d(kcloud)
        character*4 lvl_coord_3d(kcloud)

        logical l_packed_data
        data l_packed_data /.false./

        common /laps_diag/ no_laps_diag

        if(no_laps_diag .eq. 0)write(6,*)' getting cloud analysis'

        call get_directory(ext,directory,len_dir)
        c_filespec = directory(1:len_dir)//'*.'//ext(1:3)

        call get_file_time(c_filespec,i4time_needed,i4time_nearest)

        if(i4time_nearest .eq. 0)then
            write(6,*)' error: no cloud (',ext(1:3),') files available'
            istatus = 0
            return
        endif

!       fill in cloud pressures array with standard atmosphere values
!       this is relavant for old files only

        i4time_switch = 969476400

        if(i4time_nearest .lt. i4time_switch)then
            if(kcloud .ne. kcloudin)then
                write(6,*)' error: incorrect vertical cloud dimension'
                istatus = 0
                return
            endif
            do i = 1,kcloud
                cld_hts(i) = cld_hts_old(i)
            enddo
        else
          if(kcloud .eq. kcloudin)then
            do i = 1,kcloud
                cld_hts(i) = cld_hts_new(i)
            enddo
          endif
        endif

        if(kcloud .eq. kcloudin)then
          do k=1,kcloud
            cld_pres(k) = ztopsa(cld_hts(k)) * 100.
          enddo ! k
        endif

!       read in cloud grid

        if(l_packed_data)then

            if(no_laps_diag .eq. 0)
     1       write(6,*)' reading in byte array from ',ext(1:3),' directo
     1ry'
            open(21,file=directory(1:len_dir)//filename13(i4time_nearest
     1,ext(1:3))
     1  ,access='sequential',form='unformatted',status='old',err=998)
!           read(21,err=998)(((b_cldcv(i,j,k),i=1,imax),j=1,jmax),k=1,kcloudin)
            read(21,err=100)(cld_hts(k),k=1,kcloudin),(cld_pres(k),k=1,k
     1cloudin)

            if(no_laps_diag .eq. 0)write(6,*)' updated cloud heights/pre
     1ssures'
            close(21)

            goto90

100         write(6,*)' using standard atmosphere for cloud pressure arr
     1ay'
            close(21)

!           convert cloud cover array from byte format
90          do k = 1,kcloudin
            do j = 1,jmax
            do i = 1,imax
!               if(b_cldcv(i,j,k) .ne. -1)then
!                   i_temp = b_cldcv(i,j,k)
!                   clouds_3d(i,j,k) = i_temp / 100.
!               else
!                   clouds_3d(i,j,k) = 0.0
!               endif
            enddo ! i
            enddo ! j
            enddo ! k

            istatus = 1
            return

998         write(6,*)' error reading cloud (lc3) file'
            istatus = 0
            return

        else ! use readlapsdata

            do k = 1,kcloudin
                units_3d(k) = '    '
                lvl_3d(k) = k
                lvl_coord_3d(k) = '   '
                var_3d(k) = var_2d
            enddo ! k

            call read_laps_data(i4time_nearest,directory,ext,imax,jmax,
     1  kcloudin,kcloud,var_3d,lvl_3d,lvl_coord_3d,units_3d,
     1                     comment_3d,clouds_3d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error reading cloud (lc3) file'
                return
            endif

!           decode cloud heights and pressures
            do k = 1,kcloudin
                read(comment_3d(k),1,err=999)cld_hts(k),cld_pres(k)
1               format(2e20.8)
            enddo ! k

        endif

        return ! normal return

999     write(6,*)' error reading comment field in lc3 file'
        istatus = 0
        return

        end

        subroutine get_clouds_3dpres(i4time_needed,imax,jmax,kmax
     1  ,kcloud,clouds_3d_pres,clouds_3d_height ! ,b_cldcv
     1  ,cld_hts,cld_pres_1d,cld_pres_3d,istatus)

        integer i4time_needed                      ! input
        integer imax,jmax,kmax                     ! input (laps dims)
        integer kcloud                             ! input (normally 42)
        real clouds_3d_pres(imax,jmax,kmax)        ! output
        real clouds_3d_height(imax,jmax,kcloud)    ! local dummy array
!       byte   b_cldcv(imax,jmax,kcloud)             ! local dummy array
        real cld_hts(kcloud)                       ! local dummy array
        real cld_pres_1d(kcloud)                   ! local dummy array
        real cld_pres_3d(imax,jmax,kcloud)         ! local dummy array
        integer istatus                            ! output

        character*31 ext

        ext = 'lc3'

        call get_clouds_3dgrid(i4time_needed,i4time_nearest
     1          ,imax,jmax,kcloud,ext ! ,b_cldcv
     1          ,clouds_3d_height,cld_hts,cld_pres_1d,istatus)

        if(istatus .ne. 1)return

        if(i4time_needed .ne. i4time_nearest)then
            istatus = 0
            return
        endif

!       create cld_pres_3d array
        do k = 1,kcloud
        do j = 1,jmax
        do i = 1,imax
            cld_pres_3d(i,j,k) = cld_pres_1d(k)
        enddo ! i
        enddo ! j
        enddo ! k

        call interp_height_pres(imax,jmax,kmax,kcloud
     1  ,clouds_3d_pres,clouds_3d_height,cld_pres_3d,istatus)

        return
        end


        subroutine interp_height_pres(imax,jmax,kmax,kcloud
     1  ,clouds_3d_pres,clouds_3d_height,cld_pres_3d,istatus)

        integer imax,jmax,kmax                     ! input (laps dims)
        integer kcloud                             ! input (normally 42)
        real clouds_3d_pres(imax,jmax,kmax)        ! output
        real clouds_3d_height(imax,jmax,kcloud)    ! input
        real cld_pres_3d(imax,jmax,kcloud)         ! input
        integer istatus                            ! output

!       interpolate from height grid to pressure grid (generating inputs can
!       lead to slowness or inaccuracies)
        do lvl = 1,kmax
            rlvl = pressure_of_level(lvl)

            do j = 1,jmax
            do i = 1,imax

!             efficiency speedup
              if(i .eq. 1 .or. i .eq. i/5*5)then

!               default if out of bounds
                kht_low = 1
                kht_high = 1
                frac_low = 1.0
                frac_high = 0.0

                do kht = kcloud-1,1,-1
                    if(cld_pres_3d(i,j,kht  ) .ge. rlvl .and.
     1               cld_pres_3d(i,j,kht+1) .le. rlvl           )then
                        kht_low = kht
                        kht_high = kht_low + 1
                        frac_high = (rlvl - cld_pres_3d(i,j,kht)) /
     1                        (cld_pres_3d(i,j,kht+1) - cld_pres_3d(i,j,
     1kht))
                        frac_low  = 1.0 - frac_high
                        goto 10 ! we found the cloud level, kick out of loop
                    endif
                enddo ! kht

10            endif

!d            write(6,1)lvl,kht_low,kht_high,frac_low,frac_high
!d      1            ,cld_pres_3d(i,j,kht_low),rlvl,cld_pres_3d(i,j,kht_high)
!d1           format(1x,' reading clouds lvl ',i5,3x,2i5,3x,2f8.3,2x,3f8.2)

              clouds_3d_pres(i,j,lvl)
     1                = clouds_3d_height(i,j,kht_low ) * frac_low
     1                      + clouds_3d_height(i,j,kht_high) * frac_high
            enddo ! i
            enddo ! j

        enddo ! lvl

        return
        end


        subroutine interp_height_pres_fast(imax,jmax,kmax,kcloud
     1  ,clouds_3d_pres,clouds_3d_height,heights_3d,cld_hts,istatus)

        integer imax,jmax,kmax                     ! input (laps dims)
        integer kcloud                             ! input (normally 42)
        real clouds_3d_pres(imax,jmax,kmax)        ! output
        real clouds_3d_height(imax,jmax,kcloud)    ! input
        real heights_3d(imax,jmax,kmax)            ! input
        real cld_hts (kcloud)
        integer istatus                            ! output

!       interpolate from height grid to pressure grid (this is fast and
!       accurate when heights are supplied)

        do j = 1,jmax
        do i = 1,imax

            kstart = 1

            do lvl = 1,kmax

                do kht = kstart,kcloud-1
                    if(cld_hts(kht+1) .ge. heights_3d(i,j,lvl)  )then
                        kht_low  = kht
                        kht_high = kht_low + 1

                        if(cld_hts(kht) .le. heights_3d(i,j,lvl))then
                            frac_high = (heights_3d(i,j,lvl) - cld_hts(k
     1ht)) /
     1                                (cld_hts(kht+1)      - cld_hts(kht
     1))
                            frac_low  = 1.0 - frac_high
                        else ! take care of lower boundary effects
                            kht_low = 1
                            kht_high = 1
                            frac_low = 1.0
                            frac_high = 0.0
                        endif

                        goto 10 ! we found the cloud level, kick out of loop
                    endif
                enddo ! kht

10            kstart = kht

              clouds_3d_pres(i,j,lvl)
     1                = clouds_3d_height(i,j,kht_low ) * frac_low
     1                + clouds_3d_height(i,j,kht_high) * frac_high

              if(i .eq. imax/2 .and. j .eq. jmax/2)then
                  write(6,1)lvl,kht_low,kht_high,frac_low,frac_high
     1                     ,clouds_3d_height(i,j,kht_low)
     1                     ,clouds_3d_height(i,j,kht_high)
     1                     ,heights_3d(i,j,lvl)
     1                     ,clouds_3d_pres(i,j,lvl)
 1                format(1x,' reading clouds lvl ',i5,3x,2i5,3x,2f8.3
     1                                            ,2x,2f8.2,f8.1,f8.2)
              endif

            enddo ! lvl
        enddo ! i
        enddo ! j

        return
        end
