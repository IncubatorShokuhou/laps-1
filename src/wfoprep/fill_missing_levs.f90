subroutine fill_missing_levs(nx, ny, nz, plevs, data, missval, method)

! subroutine to fill in missing values in a 3d data array on
! pressure levels.  this routine assumes you have valid data
! on the top most and bottom most levels.

   implicit none
   integer, intent(in)                 :: nx
   integer, intent(in)                 :: ny
   integer, intent(in)                 :: nz
   real, intent(in)                    :: plevs(75)  ! press in pa
   real, intent(inout)                 :: data(nx, ny, nz)
   real, intent(in)                    :: missval
   integer, intent(in)                 :: method

   integer, parameter                  :: method_lin = 1
   integer, parameter                  :: method_log = 2

   integer                             :: k, kbot, ktop
   real                                :: weight_top, weight_bot
   logical                             :: goodlev(nz)

   goodlev(:) = .false.
   checklevs: do k = 1, nz
      if (maxval(data(:, :, k)) .ne. missval) then
         goodlev(k) = .true.
      end if
   end do checklevs

   fixlevs: do k = 1, nz
      if (goodlev(k)) cycle fixlevs

      ! else:

      ! find first level below that is "good"
      find_kbot: do kbot = k - 1, 1, -1
         if (goodlev(kbot)) exit find_kbot
      end do find_kbot

      ! find first level above that is good
      find_ktop: do ktop = k + 1, nz
         if (goodlev(ktop)) exit find_ktop
      end do find_ktop

      ! compute interp weights using plevs and method
      if (method .eq. method_lin) then
         ! linear interp in pressure
         weight_bot = (plevs(k) - plevs(ktop))/ &
                      (plevs(kbot) - plevs(ktop))
         weight_top = 1.-weight_bot
      else if (method .eq. method_log) then
         ! linear in ln(p)
         weight_bot = (alog(plevs(k)/plevs(ktop)))/ &
                      (alog(plevs(kbot)/plevs(ktop)))
         weight_top = 1.-weight_bot
      end if
      data(:, :, k) = weight_bot*data(:, :, kbot) + weight_top*data(:, :, ktop)
   end do fixlevs
   return
end subroutine fill_missing_levs
