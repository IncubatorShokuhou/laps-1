c
c**************new routine as adapted at fsl**************************

      subroutine fire_fields(ni,nj,nk,temp_3d,td_3d                 ! i
!    1                      ,heights_3d                             ! i
     1                      ,u_3d,v_3d                              ! i
     1                      ,t_sfc_k,p_sfc_pa                       ! i
     1                      ,rh_sfc,u_sfc,v_sfc                     ! i
     1                      ,r_missing_data,i4time                  ! i
     1                      ,istatus)                               ! o

      character*10  units_2d
      character*125 comment_2d
      character*3 var_2d

      integer max_fields
      parameter (max_fields = 7)
      real field_array(ni,nj,max_fields)
      character*3 var_a(max_fields)
      character*125 comment_a(max_fields)
      character*10  units_a(max_fields)
      character*31 ext

!     real heights_3d(ni,nj,nk)                                   ! i
      real temp_3d(ni,nj,nk)                                      ! i
      real td_3d(ni,nj,nk)                                        ! i
      real u_3d(ni,nj,nk)
      real v_3d(ni,nj,nk)

      real t_sfc_k(ni,nj)                                         ! i 
      real rh_sfc(ni,nj)                                          ! i 
      real p_sfc_pa(ni,nj)                                        ! i 
      real u_sfc(ni,nj)                                           ! i
      real v_sfc(ni,nj)                                           ! i

      real pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                   ! l
      real pres_3d_pa(ni,nj,nk)                                   ! l
      real haines_mid_2d(ni,nj)                                   ! l
      real haines_hi_2d(ni,nj)                                    ! l
      real vent_2d(ni,nj)                                         ! l
      real fosberg_2d(ni,nj)                                      ! l
      real umean_2d(ni,nj),vmean_2d(ni,nj)                        ! l
      real cfwi(ni,nj)                                            ! l

      write(6,*)' subroutine fire_fields (under construction)'

      ext = 'pbl'
      istat_pbl = 1

      write(6,*)' read in pressure of pbl top'
      var_2d = 'ptp'
      call get_laps_2dgrid(i4time,0,i4time_nearest
     1                    ,ext,var_2d,units_2d,comment_2d,ni,nj
     1                    ,pbl_top_pa,0,istatus)

      if(istatus .ne. 1)then
          write(6,*)' laps pbl top not available'
          istat_pbl = 0
      endif

      write(6,*)' read in pbl depth'
      var_2d = 'pdm'
      call get_laps_2dgrid(i4time,0,i4time_nearest
     1                    ,ext,var_2d,units_2d,comment_2d,ni,nj
     1                    ,pbl_depth_m,0,istatus)

      if(istatus .ne. 1)then
          write(6,*)' laps pbl depth not available'
          istat_pbl = 0
      endif

      call get_pres_3d(i4time,ni,nj,nk,pres_3d_pa,istatus)
      if(istatus .ne. 1)then
          write(6,*)' could not obtain pres_3d'
          return
      endif

      call cpt_fire_fields(ni,nj,nk,pres_3d_pa,temp_3d,td_3d        ! i
!    1                          ,heights_3d                         ! i
     1                          ,u_3d,v_3d                          ! i
     1                          ,pbl_top_pa,pbl_depth_m,istat_pbl   ! i
     1                          ,t_sfc_k                            ! i
     1                          ,rh_sfc                             ! i
     1                          ,u_sfc                              ! i
     1                          ,v_sfc                              ! i
     1                          ,p_sfc_pa                           ! i
     1                          ,r_missing_data                     ! i
     1                          ,i4time                             ! i
     1                          ,fosberg_2d                         ! o
     1                          ,haines_mid_2d                      ! o
     1                          ,haines_hi_2d                       ! o
     1                          ,vent_2d                            ! o
     1                          ,umean_2d,vmean_2d                  ! o
     1                          ,cfwi                               ! o
     1                          ,istatus)                           ! o

      if(istatus .eq. 1)then ! write out lfr fire fields file
          write(6,*)' write out lfr file'

          call move(vent_2d      ,field_array(1,1,1),ni,nj)
          call move(haines_mid_2d,field_array(1,1,2),ni,nj)
          call move(haines_hi_2d, field_array(1,1,3),ni,nj)
          call move(fosberg_2d   ,field_array(1,1,4),ni,nj)
          call move(umean_2d     ,field_array(1,1,5),ni,nj)
          call move(vmean_2d     ,field_array(1,1,6),ni,nj)
          call move(cfwi         ,field_array(1,1,7),ni,nj)

          ext = 'lfr'
          var_a(1) = 'vnt'
          var_a(2) = 'ham'
          var_a(3) = 'hah'
          var_a(4) = 'fwi'
          var_a(5) = 'upb'
          var_a(6) = 'vpb'
          var_a(7) = 'cwi'
          units_a(1) = 'm**2/s'
          units_a(2) = '      '
          units_a(3) = '      '
          units_a(4) = '      '
          units_a(5) = 'm/s'
          units_a(6) = 'm/s'
          units_a(7) = '      '
          comment_a(1) = 'ventilation index'
          comment_a(2) = 'haines index (850-700hpa)'
          comment_a(3) = 'haines index (700-500hpa)'
          comment_a(4) = 'fosberg fire wx index'
          comment_a(5) = 'boundary layer mean u component'
          comment_a(6) = 'boundary layer mean v component'
          comment_a(7) = 'critical fire weather index'
          call put_laps_multi_2d(i4time,ext,var_a,units_a
     1                          ,comment_a,field_array,ni,nj
     1                          ,7,istatus)

      else
          write(6,*)' skipping write of lfr file'

      endif

      i4_elapsed = ishow_timer()

      return
      end

       subroutine cpt_fire_fields(ni,nj,nk,pres_3d_pa,temp_3d,td_3d  ! i
!    1                           ,heights_3d                         ! i
     1                           ,u_3d,v_3d                          ! i
     1                           ,pbl_top_pa,pbl_depth_m,istat_pbl   ! i
     1                           ,t_sfc_k                            ! i
     1                           ,rh_sfc                             ! i
     1                           ,u_sfc                              ! i
     1                           ,v_sfc                              ! i
     1                           ,p_sfc_pa                           ! i
     1                           ,r_missing_data                     ! i
     1                           ,i4time                             ! i
     1                           ,fosberg_2d                         ! o
     1                           ,haines_mid_2d                      ! o
     1                           ,haines_hi_2d                       ! o
     1                           ,vent_2d                            ! o
     1                           ,umean_2d,vmean_2d                  ! o
     1                           ,cfwi                               ! o
     1                           ,istatus)                           ! o

       real pres_3d_pa(ni,nj,nk)                                   ! i
!      real heights_3d(ni,nj,nk)                                   ! i
       real temp_3d(ni,nj,nk)                                      ! i
       real td_3d(ni,nj,nk)                                        ! i
       real u_3d(ni,nj,nk)
       real v_3d(ni,nj,nk)

       real pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                   ! i
       real pbl_top_m(ni,nj)                                       ! l

       real t_sfc_k(ni,nj)
       real rh_sfc(ni,nj)
       real u_sfc(ni,nj)
       real v_sfc(ni,nj)
       real p_sfc_pa(ni,nj)                                        ! i
       real p_sfc_mb(ni,nj)

       real haines_mid_2d(ni,nj)                                   ! o
       real haines_hi_2d(ni,nj)                                    ! o
       real vent_2d(ni,nj)                                         ! o
       real umean_2d(ni,nj),vmean_2d(ni,nj)                        ! o
       real fosberg_2d(ni,nj)                                      ! o
       real cfwi(ni,nj)                                            ! o

       real pres_3d_mb(ni,nj,nk)                                   ! l

       write(6,*)' subroutine cpt_fire_fields...'

!      calculate haines index
       write(6,*)' calculate mid-level haines index'
       pres_3d_mb = pres_3d_pa / 100.
       call hainesindex(pres_3d_mb,temp_3d,td_3d,haines_mid_2d,ni,nj,nk       
     1                 ,850.,700.) 
       do i = 1,ni
       do j = 1,nj
           if(p_sfc_pa(i,j) .lt. 85000.)then
               haines_mid_2d(i,j) = r_missing_data
           endif
       enddo ! j
       enddo ! i

       call hainesindex(pres_3d_mb,temp_3d,td_3d,haines_hi_2d,ni,nj,nk
     1                 ,700.,500.) 
       do i = 1,ni
       do j = 1,nj
           if(p_sfc_pa(i,j) .lt. 70000.)then
               haines_hi_2d(i,j) = r_missing_data
           endif
       enddo ! j
       enddo ! i

!      calculate fosberg fireweather index
       write(6,*)' calculate fosberg fireweather index'
       p_sfc_mb = p_sfc_pa / 100.
!      argument list has been changed to use sfc inputs instead of 3d
       call fireweatherindex(t_sfc_k,rh_sfc,p_sfc_mb,u_sfc,v_sfc          ! i
     1                      ,ni,nj                                        ! i
     1                      ,fosberg_2d)                                  ! o

!      calculate ventilation index
       if(istat_pbl .eq. 1)then
           write(6,*)' calculate ventilation index'
           call ventilation_index(u_3d,v_3d,pbl_top_pa,pbl_depth_m        ! i
     1                           ,pres_3d_pa,p_sfc_pa                     ! i
     1                           ,ni,nj,nk                                ! i
     1                           ,r_missing_data                          ! i
!    1                           ,heights_3d                              ! i
     1                           ,umean_2d,vmean_2d                       ! o
     1                           ,vent_2d,istatus)                        ! o
       else
           write(6,*)' skip ventilation index due to no pbl'

       endif

!      calculate critical fire weather index.
!         (rh<15% and speed>20mph for any 3 consecutive hours during the
!          past 24 hours)
       write(6,*)' calculate critical fireweather index',i4time
       call critical_fwi(rh_sfc,u_sfc,v_sfc                               ! i
     1                  ,ni,nj,i4time                                     ! i
     1                  ,cfwi)                                            ! o

       return
       end

       subroutine ventilation_index(u_3d,v_3d,pbl_top_pa,pbl_depth_m      ! i
     1                             ,pres_3d_pa,p_sfc_pa                   ! i
     1                             ,ni,nj,nk                              ! i
     1                             ,r_missing_data                        ! i
!    1                             ,heights_3d                            ! i
     1                             ,umean_2d,vmean_2d                     ! o
     1                             ,vent_2d,istatus)                      ! o

!      real heights_3d(ni,nj,nk)                                        ! i
       real pres_3d_pa(ni,nj,nk)                                        ! i
       real u_3d(ni,nj,nk)                                              ! i
       real v_3d(ni,nj,nk)                                              ! i

       real pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                        ! i
       real umean_2d(ni,nj),vmean_2d(ni,nj)                             ! o
       real vent_2d(ni,nj)                                              ! o

       real topo(ni,nj)           ! switch to sfc_pres_pa?
       real p_sfc_pa(ni,nj)                                             ! i

       write(6,*)' subroutine ventilation_index'

       vent_2d = r_missing_data

!      calculate mean wind within the pbl
       call pbl_mean_wind(u_3d,v_3d,topo,pbl_top_pa,ni,nj,nk              ! i
     1                   ,pres_3d_pa,p_sfc_pa                             ! i
     1                   ,umean_2d,vmean_2d,istatus)                      ! o
       if(istatus .ne. 1)then
           write(6,*)' warning: bad status returned from pbl_mean_wind'       
           write(6,*)' returning vent_2d field as missing data'
           return
       endif

       write(6,*)' compute vi from mean wind speed and pbl depth'

!      multiply pbl depth by mean wind to obtain ventilation index
       do i = 1,ni
       do j = 1,nj
           if(abs(umean_2d(i,j)) .gt. 1000.)then
               write(6,*)' error, umean out of bounds',i,j,umean_2d(i,j)
               istatus = 0
               return
           endif

           if(abs(vmean_2d(i,j)) .gt. 1000.)then
               write(6,*)' error, vmean out of bounds',i,j,vmean_2d(i,j)       
               istatus = 0
               return
           endif

           spmean = sqrt(umean_2d(i,j)**2 + vmean_2d(i,j)**2)
           vent_2d(i,j) = pbl_depth_m(i,j) * spmean

       enddo ! j
       enddo ! i

       return
       end



        subroutine pbl_mean_wind(uanl,vanl,topo,pbl_top_pa        ! i
     1                          ,imax,jmax,kmax                   ! i
     1                          ,pres_3d_pa,p_sfc_pa              ! i
     1                          ,umean_2d,vmean_2d,istatus)       ! o

        logical ltest_vertical_grid

        real umean_2d(imax,jmax),vmean_2d(imax,jmax)            ! output
        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)          ! input
        real pres_3d_pa(imax,jmax,kmax)                         ! input

        real topo(imax,jmax)                                    ! input
        real p_sfc_pa(imax,jmax)                                ! input
        real pbl_top_pa(imax,jmax)                              ! input

        real sum(imax,jmax)                                     ! local
        real usum(imax,jmax)                                    ! local
        real vsum(imax,jmax)                                    ! local
        integer klow(imax,jmax)                                 ! local
        integer khigh(imax,jmax)                                ! local

!       topo = 0.             ! just for testing (also switch to pres_sfc_pa)?

!       calculate the mass weighted mean wind for the pbl. inputs are the
!       lower and upper bounds in terms of pressure. fractional levels are 
!       not accounted for - only whole levels between the pressure bounds
!       are integrated.

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        umean_2d = r_missing_data
        vmean_2d = r_missing_data

        do j = 1,jmax
          do i = 1,imax
             if(pbl_top_pa(i,j) .gt. p_sfc_pa(i,j))then
                 write(6,*)' error in pbl_mean_wind: pbl top > sfc p'
     1                    ,i,j,pbl_top_pa(i,j),p_sfc_pa(i,j)
                 istatus = 0
                 return
             endif

             khigh(i,j) = rlevel_of_field(pbl_top_pa(i,j),pres_3d_pa
     1                              ,imax,jmax,kmax,i,j,istatus)
             if(istatus .ne. 1)then
                 write(6,*)' mean_wind: error in rlevel_of_field'
                 return
             endif

             klow(i,j) = rlevel_of_field(p_sfc_pa(i,j),pres_3d_pa
     1                                  ,imax,jmax,kmax,i,j,istatus)
             if(istatus .ne. 1)then
                 write(6,*)' mean_wind: error in rlevel_of_field'
                 return
             endif

             sum(i,j) = 0.
             usum(i,j) = 0.
             vsum(i,j) = 0.
          enddo ! j
        enddo ! i

        do j = 1,jmax
          do i = 1,imax
            do k = 1,khigh(i,j)
              if(uanl(i,j,k) .ne. r_missing_data .and.
     1           vanl(i,j,k) .ne. r_missing_data
     1                        .and. k .ge. klow(i,j))then
                sum(i,j) = sum(i,j) + 1.
                usum(i,j) = usum(i,j) + uanl(i,j,k)
                vsum(i,j) = vsum(i,j) + vanl(i,j,k)
              endif
             enddo
          enddo
        enddo

        do j = 1,jmax
          do i = 1,imax
             if(sum(i,j) .gt. 0.)then ! mean wind through the layer
                 umean_2d(i,j) = usum(i,j) / sum(i,j)
                 vmean_2d(i,j) = vsum(i,j) / sum(i,j)

             else 
                 write(6,*)' warning: sum <= 0',i,j,klow(i,j),khigh(i,j)       

             endif
          enddo ! i
        enddo ! j

        return
        end


       subroutine critical_fwi(rh_sfc,u_sfc,v_sfc                         ! i
     1                        ,ni,nj,i4time                               ! i
     1                        ,cfwi)                                      ! o

       implicit none

       real spd_thresh
       parameter (spd_thresh = 79.90372)   ! 20mph squared in m/s
 
       integer ni,nj,i4time,pi4time,i,j,n,istatus

       real rh_sfc(ni,nj),u_sfc(ni,nj),v_sfc(ni,nj)
     .       ,rh(ni,nj),u(ni,nj),v(ni,nj)
     .       ,cfwi1(ni,nj,24),cfwi(ni,nj)
     .       ,speed

       character*31  ext
       character*10  units_2d
       character*125 comment_2d
       character*3   var_2d

       do n=1,24
          if (n .eq. 1) then

             do j=1,nj
             do i=1,ni
                u(i,j)=u_sfc(i,j)
                v(i,j)=v_sfc(i,j)
                rh(i,j)=rh_sfc(i,j)
             enddo
             enddo

          else

             pi4time = i4time - (n-1)*3600
             print *,'critical_fwi: i4time, pi4time',i4time,pi4time

             ext = 'lsx'

             var_2d = 'u'
             call get_laps_2d(pi4time,ext,var_2d,units_2d,comment_2d
     1                       ,ni,nj,u,istatus)
             var_2d = 'v'
             call get_laps_2d(pi4time,ext,var_2d,units_2d,comment_2d
     1                       ,ni,nj,v,istatus)
             var_2d = 'rh'
             call get_laps_2d(pi4time,ext,var_2d,units_2d,comment_2d
     1                       ,ni,nj,rh,istatus)

          endif

          do j=1,nj
          do i=1,ni
             if (rh(i,j) .gt. 0 .and. rh(i,j) .lt. 15. .and.
     1           u(i,j) .lt. 100. and. v(i,j) .lt. 100.) then
                speed = u(i,j)**2 + v(i,j)**2
                if (speed .gt. spd_thresh) then
                   cfwi1(i,j,n) = 1.
                else
                   cfwi1(i,j,n) = 0.
                endif
             else
                cfwi1(i,j,n) = 0.
             endif
          enddo
          enddo

       enddo

       do j=1,nj
       do i=1,ni
          cfwi(i,j) = 0.
       enddo
       enddo

       do n=3,24
       do j=1,nj
       do i=1,ni
          if (cfwi1(i,j,n)+cfwi1(i,j,n-1)+cfwi1(i,j,n-2) .gt. 0.) 
     1       cfwi(i,j) = 1.
       enddo
       enddo
       enddo

       return
       end
