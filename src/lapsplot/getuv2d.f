
        subroutine get_uv_2d(i4time,k_level,uv_2d,ext,imax,jmax
     1                      ,fcst_hhmm,c_model,istatus)

!       97-aug-17     ken dritz     used commenting to (temporarily) hardwire
!                                   vertical_grid to 'pressure' (without
!                                   accessing vertical_grid)
!       97-aug-17     ken dritz     removed include of lapsparms.for

        character*150 directory
        character*31 ext
        character*20 c_model

        character*125 comment_2d(2)
        character*10 units_2d(2)
        character*3 var(2)
        integer lvl_2d(2)
        character*4 lvl_coord_2d(2)
        character*5 fcst_hhmm
        character*9 a9time

        real uv_2d(imax,jmax,2)

        call get_directory(ext,directory,len_dir)
        call s_len(ext,lext)
        if(ext(1:lext).eq.'balance')then
           ext='lw3'
           directory=directory(1:len_dir)//'lw3'//'/'
        endif

        do k = 1,2
          units_2d(k)   = 'm/s'
          comment_2d(k) = '3dwind'

          if(k_level .gt. 0)then
!           if(vertical_grid .eq. 'height')then
!               lvl_2d(k) = zcoord_of_level(k_level)/10
!               lvl_coord_2d(k) = 'msl'
!           elseif(vertical_grid .eq. 'pressure')then
                lvl_2d(k) = zcoord_of_level(k_level)/100
                lvl_coord_2d(k) = 'hpa'
!           endif
            var(1) = 'u3' ! newvar 'u3', oldvar = 'u'
            var(2) = 'v3' ! newvar 'v3', oldvar = 'v'

          else
            var(1) = 'u'
            var(2) = 'v'
            lvl_2d(k) = 0
            lvl_coord_2d(k) = 'agl'

          endif

        enddo ! k

        if(ext(1:3) .eq. 'lga' .or. ext(1:3) .eq. 'fua')then
            call get_laps_cycle_time(laps_cycle_time,istatus)

            call input_background_info(
     1                              ext                     ! i
     1                             ,directory,c_model       ! o
     1                             ,i4time                  ! i
     1                             ,laps_cycle_time         ! i
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o
            if(istatus.ne.1)return

            call read_laps(i4_initial,i4_valid,directory,ext
     1                    ,imax,jmax,2,2,var,lvl_2d,lvl_coord_2d
     1                    ,units_2d,comment_2d,uv_2d,istatus)
            if(istatus .ne. 1)then
                write(6,*)
     1          ' sorry, unable to read in wind variables ',var          
                return
            else
                write(6,*)
     1          ' 2d - laps u and v analysis successfully read in'
     1          ,lvl_2d(1)        
            endif

            i4time = i4_valid

        else

!           read in 2d w array
            call read_laps_data(i4time,directory,ext,imax,jmax,2,2,
     1          var,lvl_2d,lvl_coord_2d,units_2d,comment_2d,
     1          uv_2d,istatus)
            if(istatus .ne. 1)then
                write(6,*)
     1      ' sorry, file has not yet been generated this hour'
                stop
            else
                write(6,*)
     1    ' 2d - laps u and v analysis successfully read in',lvl_2d(1)        
            endif

        endif

        return
        end

