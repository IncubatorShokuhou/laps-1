
subroutine staggerxy_3d(vin, nx, ny, nz, onx, ony, onz, vout)

!==========================================================
!  this routine computes a stagger grid
!  using a given horizontal uniform grid
!
!  input:
!        vin:    data prepare to stagger
!        nx:        x grid point number before stagger
!        ny:        y grid point number before stagger
!        nz:        z grid point number before stagger
!       onx:    x grid point number after stagger
!       ony:    y grid point number after stagger
!       onz:    z grid point number after stagger
!
!  output:
!        vout:   data after stagger
!
!==========================================================

   implicit none

   integer, intent(in) :: nx, ny, nz, onx, ony, onz
   real, intent(in) ::    vin(nx, ny, nz)
   real, intent(out) ::   vout(onx, ony, onz)

   ! local variables:

   real :: vinx(nx - 1, ny, nz), viny(nx, ny - 1, nz), vinyy(nx - 1, ny - 1, nz)

   ! linear interpolation:

   if (nx /= onx .and. ny /= ony .and. onz /= 1) then

      vinx(1:nx - 1, 1:ny, 1:nz) = 0.5*( &
                                   vin(1:nx - 1, 1:ny, 1:nz) + &
                                   vin(2:nx, 1:ny, 1:nz))
      vinyy(1:nx - 1, 1:ny - 1, 1:nz) = 0.5*( &
                                        vinx(1:nx - 1, 1:ny - 1, 1:nz) + &
                                        vinx(1:nx - 1, 2:ny, 1:nz))

      vout(1:onx, 1:ony, 1:nz) = vinyy(1:nx - 1, 1:ny - 1, 1:nz)

   elseif (nx /= onx .and. ny /= ony .and. onz == 1) then

      vinx(1:nx - 1, 1:ny, 1) = 0.5*( &
                                vin(1:nx - 1, 1:ny, 1) + &
                                vin(2:nx, 1:ny, 1))
      vinyy(1:nx - 1, 1:ny - 1, 1) = 0.5*( &
                                     vinx(1:nx - 1, 1:ny - 1, 1) + &
                                     vinx(1:nx - 1, 2:ny, 1))

      vout(1:onx, 1:ony, 1) = vinyy(1:nx - 1, 1:ny - 1, 1)

   elseif (nx /= onx .and. ny == ony) then

      vinx(1:nx - 1, 1:ny, 1:nz) = 0.5*( &
                                   vin(1:nx - 1, 1:ny, 1:nz) + &
                                   vin(2:nx, 1:ny, 1:nz))

      vout(1:onx, 1:ony, 1:nz) = vinx(1:nx - 1, 1:ny, 1:nz)

   elseif (nx == onx .and. ny /= ony) then

      viny(1:nx, 1:ny - 1, 1:nz) = 0.5*( &
                                   vin(1:nx, 1:ny - 1, 1:nz) + &
                                   vin(1:nx, 2:ny, 1:nz))

      vout(1:onx, 1:ony, 1:nz) = viny(1:nx, 1:ny - 1, 1:nz)

   end if

end subroutine staggerxy_3d

subroutine staggerxy_2d(vin, nx, ny, onx, ony, vout)

!==========================================================
!  this routine computes a stagger grid
!  using a given horizontal uniform grid
!
!  input:
!        vin:    data prepare to stagger
!        nx:        x grid point number before stagger
!        ny:        y grid point number before stagger
!       onx:    x grid point number after stagger
!       ony:    y grid point number after stagger
!
!  output:
!        vout:   data after stagger
!
!==========================================================

   implicit none

   integer, intent(in) :: nx, ny, onx, ony
   real, intent(in) ::    vin(nx, ny)
   real, intent(out) ::   vout(onx, ony)

   ! local variables:

   real :: vinx(nx - 1, ny), viny(nx, ny - 1), vinyy(nx - 1, ny - 1)

   ! linear interpolation:

   if (nx /= onx .and. ny /= ony) then

      vinx(1:nx - 1, 1:ny) = 0.5*( &
                             vin(1:nx - 1, 1:ny) + &
                             vin(2:nx, 1:ny))
      vinyy(1:nx - 1, 1:ny - 1) = 0.5*( &
                                  vinx(1:nx - 1, 1:ny - 1) + &
                                  vinx(1:nx - 1, 2:ny))

      vout(1:onx, 1:ony) = vinyy(1:nx - 1, 1:ny - 1)

   elseif (nx /= onx .and. ny == ony) then

      vinx(1:nx - 1, 1:ny) = 0.5*( &
                             vin(1:nx - 1, 1:ny) + &
                             vin(2:nx, 1:ny))

      vout(1:onx, 1:ony) = vinx(1:nx - 1, 1:ny)

   elseif (nx == onx .and. ny /= ony) then

      viny(1:nx, 1:ny - 1) = 0.5*( &
                             vin(1:nx, 1:ny - 1) + &
                             vin(1:nx, 2:ny))

      vout(1:onx, 1:ony) = viny(1:nx, 1:ny - 1)

   end if

end subroutine staggerxy_2d

subroutine untaggerxy_3d(vin, nx, ny, nz, onx, ony, vout)

!==========================================================
!  this routine computes a unstagger grid
!  using a given horizontal stagger grid
!
!  input:
!        vin:    stagger data prepare to unstagger
!        nx:        number of stagger grid in x
!        ny:     number of stagger grid in y
!        nz:     z grid point number after unstagger
!       onx:    x grid point number after unstagger
!       ony:    y grid point number after unstagger
!
!  output:
!        vout:   data after unstagger
!
!==========================================================

   implicit none

   integer, intent(in) :: nx, ny, nz, onx, ony
   real, intent(in) ::    vin(nx, ny, nz)
   real, intent(out) ::   vout(onx, ony, nz)

   ! local variables:

   integer :: i, j, k
   real :: vinx(onx, ony, nz), viny(onx, ony, nz), vinyy(nx, ony, nz)

   ! linear interpolation:

   if (nx /= onx .and. ny /= ony) then

      do k = 1, nz
         do i = 1, nx
            do j = 2, ny
               vinyy(i, j, k) = &
                  0.5*(vin(i, j - 1, k) + vin(i, j, k))
            end do
            vinyy(i, 1, k) = &
               1.5*vin(i, 1, k) - &
               0.5*vin(i, 2, k)
            vinyy(i, ony, k) = &
               1.5*vin(i, ny, k) - &
               0.5*vin(i, ny - 1, k)
         end do
      end do

      do k = 1, nz
         do j = 1, ony
            do i = 2, nx
               vinx(i, j, k) = &
                  0.5*(vinyy(i - 1, j, k) + vinyy(i, j, k))
            end do
            vinx(1, j, k) = &
               1.5*vinyy(1, j, k) - &
               0.5*vinyy(2, j, k)
            vinx(onx, j, k) = &
               1.5*vinyy(nx, j, k) - &
               0.5*vinyy(nx - 1, j, k)
         end do
      end do

      vout(1:onx, 1:ony, 1:nz) = vinx(1:onx, 1:ony, 1:nz)

   elseif (nx /= onx .and. ny == ony) then

      do k = 1, nz
         do j = 1, ony
            do i = 2, nx
               vinx(i, j, k) = &
                  0.5*(vin(i - 1, j, k) + vin(i, j, k))
            end do
            vinx(1, j, k) = &
               1.5*vin(1, j, k) - &
               0.5*vin(2, j, k)
            vinx(onx, j, k) = &
               1.5*vin(nx, j, k) - &
               0.5*vin(nx - 1, j, k)
         end do
      end do

      vout(1:onx, 1:ony, 1:nz) = vinx(1:onx, 1:ony, 1:nz)

   elseif (nx == onx .and. ny /= ony) then

      do k = 1, nz
         do i = 1, onx
            do j = 2, ny
               viny(i, j, k) = &
                  0.5*(vin(i, j - 1, k) + vin(i, j, k))
            end do
            viny(i, 1, k) = &
               1.5*vin(i, 1, k) - &
               0.5*vin(i, 2, k)
            viny(i, ony, k) = &
               1.5*vin(i, ny, k) - &
               0.5*vin(i, ny - 1, k)
         end do
      end do

      vout(1:onx, 1:ony, 1:nz) = viny(1:onx, 1:ony, 1:nz)

   end if

end subroutine untaggerxy_3d

subroutine untaggerxy_2d(vin, nx, ny, onx, ony, vout)

!==========================================================
!  this routine computes a unstagger grid
!  using a given horizontal stagger grid
!
!  input:
!        vin:    stagger data prepare to unstagger
!        nx:        number of stagger grid in x
!        ny:     number of stagger grid in y
!       onx:    x grid point number after unstagger
!       ony:    y grid point number after unstagger
!
!  output:
!        vout:   data after unstagger
!
!==========================================================

   implicit none

   integer, intent(in) :: nx, ny, onx, ony
   real, intent(in) ::    vin(nx, ny)
   real, intent(out) ::   vout(onx, ony)

   ! local variables:

   integer :: i, j
   real :: vinx(onx, ony), viny(onx, ony), vinyy(nx, ony)

   ! linear interpolation:

   if (nx /= onx .and. ny /= ony) then

      do i = 1, nx
         do j = 2, ny
            vinyy(i, j) = &
               0.5*(vin(i, j - 1) + vin(i, j))
         end do
         vinyy(i, 1) = &
            1.5*vin(i, 1) - &
            0.5*vin(i, 2)
         vinyy(i, ony) = &
            1.5*vin(i, ny) - &
            0.5*vin(i, ny - 1)
      end do

      do j = 1, ony
         do i = 2, nx
            vinx(i, j) = &
               0.5*(vinyy(i - 1, j) + vinyy(i, j))
         end do
         vinx(1, j) = &
            1.5*vinyy(1, j) - &
            0.5*vinyy(2, j)
         vinx(onx, j) = &
            1.5*vinyy(nx, j) - &
            0.5*vinyy(nx - 1, j)
      end do

      vout(1:onx, 1:ony) = vinx(1:onx, 1:ony)

   elseif (nx /= onx .and. ny == ony) then

      do j = 1, ony
         do i = 2, nx
            vinx(i, j) = &
               0.5*(vin(i - 1, j) + vin(i, j))
         end do
         vinx(1, j) = &
            1.5*vin(1, j) - &
            0.5*vin(2, j)
         vinx(onx, j) = &
            1.5*vin(nx, j) - &
            0.5*vin(nx - 1, j)
      end do

      vout(1:onx, 1:ony) = vinx(1:onx, 1:ony)

   elseif (nx == onx .and. ny /= ony) then

      do i = 1, onx
         do j = 2, ny
            viny(i, j) = &
               0.5*(vin(i, j - 1) + vin(i, j))
         end do
         viny(i, 1) = &
            1.5*vin(i, 1) - &
            0.5*vin(i, 2)
         viny(i, ony) = &
            1.5*vin(i, ny) - &
            0.5*vin(i, ny - 1)
      end do

      vout(1:onx, 1:ony) = viny(1:onx, 1:ony)

   end if

end subroutine untaggerxy_2d

