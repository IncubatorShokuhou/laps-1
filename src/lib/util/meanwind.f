
        subroutine mean_wind_bunkers(uanl,vanl,topo,imax,jmax,kmax  ! i
     1                              ,heights_3d                     ! i
     1                              ,umean,vmean                    ! o
     1                              ,ushear,vshear                  ! o
     1                              ,ustorm,vstorm,istatus)         ! o

        logical ltest_vertical_grid

        real umean(imax,jmax),vmean(imax,jmax)                    ! o
        real ustorm(imax,jmax),vstorm(imax,jmax)                  ! o
        real ushear(imax,jmax),vshear(imax,jmax)                  ! o
        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)            ! i
        real heights_3d(imax,jmax,kmax)                           ! i

        real topo(imax,jmax)                                      ! i

        real sum(imax,jmax)                                       ! l
        real usum(imax,jmax)                                      ! l
        real vsum(imax,jmax)                                      ! l
        integer klow(imax,jmax)                                   ! l
        integer khigh(imax,jmax)                                  ! l

        write(6,*)
        write(6,*)' calculating mean wind (bsm)'

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        do j = 1,jmax
          do i = 1,imax
!            layer is 0-6 km agl, denoted from "sfc" to "top"

             klow(i,j) = nint(height_to_zcoord2(topo(i,j)      
     1                       ,heights_3d,imax,jmax,kmax,i,j,istatus))
             if(istatus .ne. 1)return

             khigh(i,j) = nint(height_to_zcoord2(topo(i,j)+6000.
     1                        ,heights_3d,imax,jmax,kmax,i,j,istatus))
             if(istatus .ne. 1)return

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
            enddo ! k
          enddo ! i
        enddo ! j

        do j = 1,jmax
          do i = 1,imax

c            compute storm motion vector (a la the ruc code)
c            it is defined as 7.5 m/s to the right of the 0-6 km mean
c            wind constrained along a line which is both perpendicular
c            to the 0-6 km mean vertical wind shear vector and passes
c            through the 0-6 km mean wind.  

!            mean wind through the layer
             umean(i,j) = usum(i,j) / sum(i,j)
             vmean(i,j) = vsum(i,j) / sum(i,j)

!            shear vector through the layer
             if(uanl(i,j,klow(i,j)) .ne. r_missing_data .and.
     1          vanl(i,j,klow(i,j)) .ne. r_missing_data       )then
                 ushear(i,j) = uanl(i,j,khigh(i,j))-uanl(i,j,klow(i,j))
                 vshear(i,j) = vanl(i,j,khigh(i,j))-vanl(i,j,klow(i,j))

                 shearspeed = sqrt(ushear(i,j)*ushear(i,j)
     1                            +vshear(i,j)*vshear(i,j))

                 ustorm(i,j) = umean(i,j) + (7.5*vshear(i,j)/shearspeed)
                 vstorm(i,j) = vmean(i,j) - (7.5*ushear(i,j)/shearspeed)

             else
                 write(6,*)' error in meanwind, missing low level wind'       

             endif


          enddo ! i
        enddo ! j

        return
        end

        subroutine mean_wind(uanl,vanl,topo,imax,jmax,kmax
     1                                  ,umean,vmean
     1                                  ,ustorm,vstorm,istatus)

        logical ltest_vertical_grid

        real umean(imax,jmax),vmean(imax,jmax)                  ! output
        real ustorm(imax,jmax),vstorm(imax,jmax)                ! output
        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)          ! input

        real topo(imax,jmax)                                    ! input

        real sum(imax,jmax)                                     ! local
        real usum(imax,jmax)                                    ! local
        real vsum(imax,jmax)                                    ! local
        integer klow(imax,jmax)                                 ! local

        write(6,*)
        write(6,*)' calculating mean wind (lsm)'

        if(ltest_vertical_grid('height'))then
            khigh = nint(height_to_zcoord(5000.,istatus))
        elseif(ltest_vertical_grid('pressure'))then
            pres_mb = 300.
            pres_pa = pres_mb * 100.
            khigh = nint(zcoord_of_pressure(pres_pa))
        else
            write(6,*)' mean_wind: unknown vertical grid'
            istatus = 0
            return
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        write(6,*)' top level of mean wind computation = ',khigh

!       mean wind (mass weighted) is calculated for whole levels within the 
!       range.
        do j = 1,jmax
          do i = 1,imax
             klow(i,j) =
     1            max(nint(height_to_zcoord(topo(i,j),istatus)),1)
             if(istatus .ne. 1)then
                 write(6,*)' mean_wind: error in height_to_zcoord'
                 return
             endif
             sum(i,j) = 0.
             usum(i,j) = 0.
             vsum(i,j) = 0.
          enddo ! j
        enddo ! i

        do k = 1,khigh
          do j = 1,jmax
            do i = 1,imax
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

!            mean wind through the layer
             umean(i,j) = usum(i,j) / sum(i,j)
             vmean(i,j) = vsum(i,j) / sum(i,j)

!            shear vector through the layer
             ushear = uanl(i,j,khigh) - uanl(i,j,klow(i,j))
             vshear = vanl(i,j,khigh) - vanl(i,j,klow(i,j))

!            estimate storm motion of a right moving storm
!            rotate shear vector by 90 deg, multiply by .15, add to mean wind
             ustorm(i,j) = umean(i,j) ! + 0.15 * vshear
             vstorm(i,j) = vmean(i,j) ! - 0.15 * ushear

          enddo ! i
        enddo ! j

        return
        end
