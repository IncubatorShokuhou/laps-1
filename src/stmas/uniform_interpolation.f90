!>
!! this is a routine for an interpolation between two uniform grids over the same domain
!!
!! \author yuanfu xie
!! \b history: feb. 2014
!

subroutine uniform_interpolation2(nfrom, nto, vfrom, vto)

   implicit none

   integer, intent(in) :: nfrom(2), nto(2)
   real, intent(in) :: vfrom(nfrom(1), nfrom(2))
   real, intent(out) :: vto(nto(1), nto(2))

   ! local variables:
   integer :: indx(2, 2), i, j, ii, jj
   real    :: coef(2, 2)

   do j = 1, nto(2)
      coef(1, 2) = float(nfrom(2) - 1)*(j - 1)/float(nto(2) - 1) + 1.0 ! real position from (1,nfrom)
      indx(1, 2) = int(coef(1, 2))                              ! left grid
      coef(2, 2) = coef(1, 2) - indx(1, 2)                         ! distance to left: weight on right
      indx(2, 2) = min(indx(1, 2) + 1, nfrom(2))                   ! right grid
      coef(1, 2) = 1.0 - coef(2, 2)                               ! distance to right: weight to left
      do i = 1, nto(1)
         coef(1, 1) = float(nfrom(1) - 1)*(i - 1)/float(nto(1) - 1) + 1.0 ! real position from (1,nfrom)
         indx(1, 1) = int(coef(1, 1))                              ! left grid
         coef(2, 1) = coef(1, 1) - indx(1, 1)                         ! distance to left: weight on right
         indx(2, 1) = min(indx(1, 1) + 1, nfrom(1))                   ! right grid
         coef(1, 1) = 1.0 - coef(2, 1)                               ! distance to right: weight to left

         ! interpolation:
         vto(i, j) = 0.0
         do jj = 1, 2
         do ii = 1, 2
            vto(i, j) = vto(i, j) + coef(ii, 1)*coef(jj, 2)* &
                        vfrom(indx(ii, 1), indx(jj, 2))
         end do
         end do
      end do
   end do

end subroutine uniform_interpolation2

subroutine uniform_interpolation3(nfrom, nto, vfrom, vto)

   implicit none

   integer, intent(in) :: nfrom(3), nto(3)
   real, intent(in) :: vfrom(nfrom(1), nfrom(2), nfrom(3))
   real, intent(out) :: vto(nto(1), nto(2), nto(3))

   ! local variables:
   integer :: indx(2, 3), i, j, k, ii, jj, kk
   real    :: coef(2, 3)

   do k = 1, nto(3)
      coef(1, 3) = float(nfrom(3) - 1)*(k - 1)/float(nto(3) - 1) + 1.0 ! real position from (1,nfrom)
      indx(1, 3) = int(coef(1, 3))                              ! left grid
      coef(2, 3) = coef(1, 3) - indx(1, 3)                         ! distance to left: weight on right
      indx(2, 3) = min(indx(1, 3) + 1, nfrom(3))                   ! right grid
      coef(1, 3) = 1.0 - coef(2, 3)                               ! distance to right: weight to left
      do j = 1, nto(2)
         coef(1, 2) = float(nfrom(2) - 1)*(j - 1)/float(nto(2) - 1) + 1.0 ! real position from (1,nfrom)
         indx(1, 2) = int(coef(1, 2))                              ! left grid
         coef(2, 2) = coef(1, 2) - indx(1, 2)                         ! distance to left: weight on right
         indx(2, 2) = min(indx(1, 2) + 1, nfrom(2))                   ! right grid
         coef(1, 2) = 1.0 - coef(2, 2)                               ! distance to right: weight to left
         do i = 1, nto(1)
            coef(1, 1) = float(nfrom(1) - 1)*(i - 1)/float(nto(1) - 1) + 1.0 ! real position from (1,nfrom)
            indx(1, 1) = int(coef(1, 1))                              ! left grid
            coef(2, 1) = coef(1, 1) - indx(1, 1)                         ! distance to left: weight on right
            indx(2, 1) = min(indx(1, 1) + 1, nfrom(1))                   ! right grid
            coef(1, 1) = 1.0 - coef(2, 1)                               ! distance to right: weight to left

            ! interpolation:
            vto(i, j, k) = 0.0
            do kk = 1, 2
            do jj = 1, 2
            do ii = 1, 2
               vto(i, j, k) = vto(i, j, k) + coef(ii, 1)*coef(jj, 2)*coef(kk, 3)* &
                              vfrom(indx(ii, 1), indx(jj, 2), indx(kk, 3))
            end do
            end do
            end do
         end do
      end do
   end do

end subroutine uniform_interpolation3
