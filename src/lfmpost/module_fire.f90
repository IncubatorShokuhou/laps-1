!dis
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis
!dis

module fire

   ! module to contain subroutines for computing various fire weather
   ! indices.

   implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ventilation(usig, vsig, zsig, pblhgt, topo, nx, ny, nz, &
                          upbl, vpbl, vent_ind)

      implicit none
      integer, intent(in)    :: nx
      integer, intent(in)    :: ny
      integer, intent(in)    :: nz
      real, intent(in)       :: usig(nx, ny, nz)  ! u on sigma
      real, intent(in)       :: vsig(nx, ny, nz)  ! v on sigma
      real, intent(in)       :: zsig(nx, ny, nz)  ! z on sigma
      real, intent(in)       :: pblhgt(nx, ny)
      real, intent(in)       :: topo(nx, ny)
      real, intent(out)      :: upbl(nx, ny)
      real, intent(out)      :: vpbl(nx, ny)
      real, intent(out)      :: vent_ind(nx, ny)

      integer                :: i, j, k, nbl
      real                   :: usum, vsum, umean, vmean, spmean

      print *, '---- subroutine ventilation ----'

      do j = 1, ny
         do i = 1, nx

            if (pblhgt(i, j) .gt. 0.) then

               ! compute mean wind within the pbl
               nbl = 0
               usum = 0.
               vsum = 0.
               umean = 0.
               vmean = 0.
               sum_pbl: do k = 1, nz
                  if (zsig(i, j, k) - topo(i, j) .le. pblhgt(i, j)) then
                     nbl = nbl + 1
                     usum = usum + usig(i, j, k)
                     vsum = vsum + vsig(i, j, k)
                  else
                     exit sum_pbl
                  end if
               end do sum_pbl
               if (nbl .gt. 0) then
                  umean = usum/float(nbl)
                  vmean = vsum/float(nbl)

                  ! compute mean wind speed for this layer
                  spmean = sqrt(umean**2 + vmean**2)

                  ! multiply mean speed by pbl depth to get index
                  vent_ind(i, j) = pblhgt(i, j)*spmean
                  upbl(i, j) = umean
                  vpbl(i, j) = vmean
               else
                  ! pbl height is lower than the lowest model level...use
                  ! lowest model wind
                  spmean = sqrt(usig(i, j, 1)**2 + vsig(i, j, 1)**2)
                  vent_ind(i, j) = pblhgt(i, j)*spmean
                  upbl(i, j) = usig(i, j, 1)
                  vpbl(i, j) = vsig(i, j, 1)
               end if
            else
               print *, 'warning:  pbl height <=0 in ventilation index'
               vent_ind(i, j) = 0.
               upbl(i, j) = 0.
               vpbl(i, j) = 0.
            end if
         end do
      end do
      return
   end subroutine ventilation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine haines_layer(p3d_mb, t3d_k, td3d_k, haines2d, nx, ny, nz, &
                           pmbbot, pmbtop)

      ! computes haines index for layer bounded by top and bottom pressure
      ! levels pmbtop and pmbbot (mb).

      implicit none

      integer, intent(in)      :: nx, ny, nz
      real, intent(in)         :: p3d_mb(nx, ny, nz)  ! 3d pressure in mb
      real, intent(in)         :: t3d_k(nx, ny, nz)  ! 3d temp in k
      real, intent(in)         :: td3d_k(nx, ny, nz)  ! 3d dewpoint in k
      real, intent(out)        :: haines2d(nx, ny)   ! 2d haines index
      real, intent(in)         :: pmbbot, pmbtop     ! bounding mb levels

      ! local variables

      integer :: i, j, k, kk, km1
      real :: tmkt, tmkb, tdkb, deltat, dpdep
      real :: factor1, factor2
      print *, '---- subroutine haines_layer ----'
      do j = 1, ny
         do i = 1, nx

            if (p3d_mb(i, j, 1) .lt. pmbbot) then
               haines2d(i, j) = 1e37  ! cannot be computed
            else

               do k = 2, nz

                  ! account for flipped vertical coordinate from original
                  ! algorithm from seattle

                  kk = nz + 1 - k
                  km1 = kk + 1

                  ! find temperature at the top
                  if (p3d_mb(i, j, kk) .gt. pmbtop .and. p3d_mb(i, j, km1) .le. pmbtop) then
                     tmkt = t3d_k(i, j, km1) + (t3d_k(i, j, kk) - t3d_k(i, j, km1))* &
                            (log(pmbtop) - log(p3d_mb(i, j, km1)))/ &
                            (log(p3d_mb(i, j, kk)) - log(p3d_mb(i, j, km1)))
                  end if

                  ! find temp/dewpoint at the bottom of the layer

                  if (p3d_mb(i, j, kk) .gt. pmbbot .and. p3d_mb(i, j, km1) .le. pmbbot) then
                     tmkb = t3d_k(i, j, km1) + (t3d_k(i, j, kk) - t3d_k(i, j, km1))* &
                            (log(pmbbot) - log(p3d_mb(i, j, km1)))/ &
                            (log(p3d_mb(i, j, kk)) - log(p3d_mb(i, j, km1)))
                     tdkb = td3d_k(i, j, km1) + (td3d_k(i, j, kk) - td3d_k(i, j, km1))* &
                            (log(pmbbot) - log(p3d_mb(i, j, km1)))/ &
                            (log(p3d_mb(i, j, kk)) - log(p3d_mb(i, j, km1)))
                  end if

               end do

               deltat = tmkb - tmkt
               dpdep = tmkb - tdkb

               if (nint(pmbbot) .eq. 700) then ! high haines
                  if (deltat .le. 17.5) then
                     factor1 = 1.
                  elseif (deltat .gt. 17.5 .and. deltat .le. 21.5) then
                     factor1 = 2.   ! deltat > 21.5
                  else
                     factor1 = 3.
                  end if

                  if (dpdep .le. 14.5) then
                     factor2 = 1.
                  elseif (dpdep .gt. 14.5 .and. dpdep .le. 20.5) then
                     factor2 = 2.
                  else    ! dpdep > 20.5
                     factor2 = 3.
                  end if

               elseif (nint(pmbbot) .eq. 850) then   ! mid-level haines
                  if (deltat .le. 5.5) then
                     factor1 = 1.
                  elseif (deltat .gt. 5.5 .and. deltat .le. 10.5) then
                     factor1 = 2.
                  else   ! deltat > 10.5
                     factor1 = 3.
                  end if

                  if (dpdep .le. 5.5) then
                     factor2 = 1.
                  elseif (dpdep .gt. 5.5 .and. dpdep .le. 12.5) then
                     factor2 = 2.
                  else   ! dpdep > 12.5
                     factor2 = 3.
                  end if

               else ! cannot determine mid or high

                  print *, 'haines_layer needs 850 or 700 as bottom layer'
                  print *, 'bottom level (mb) specified: ', pmbbot

                  stop 'bad_haines_layer'

               end if

               haines2d(i, j) = factor1 + factor2

            end if
         end do
      end do
      return
   end subroutine haines_layer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fosberg_fwi(t2k, rh2pct, p_sfc_mb, u10, v10, &
                          nx, ny, fwi)

      ! computes the fosberg fire weather index

      implicit none
      integer, intent(in)              :: nx, ny
      real, intent(in)                 :: t2k(nx, ny)   ! sfc temp (k)
      real, intent(in)                 :: rh2pct(nx, ny)   ! sfc rh (%)
      real, intent(in)                 :: p_sfc_mb(nx, ny) ! sfc press (mb)
      real, intent(in)                 :: u10(nx, ny)      ! 10 m u wind (m/s)
      real, intent(in)                 :: v10(nx, ny)      ! 10 m v wind (m/s)
      real, intent(out)                :: fwi(nx, ny)      ! fosberg index

      ! local variables

      integer :: i, j
      real :: uuu10, vvv10, m, n, t2f, rh2

      do j = 1, ny
         do i = 1, nx

            ! convert temperature from k to f
            t2f = 1.8*(t2k(i, j) - 273.15) + 32.0

            ! convert u/v from m/s to mph
            uuu10 = u10(i, j)*2.237
            vvv10 = v10(i, j)*2.237
            rh2 = rh2pct(i, j)

            if (rh2 .le. 10.5) then
               m = 0.03229 + (0.281073*rh2) - (0.000578*rh2*t2f)
            else if (rh2 .gt. 10.5 .and. rh2 .le. 50.5) then
               m = 2.22749 + (0.160107*rh2) - (0.014784*t2f)
            else if (rh2 .gt. 50.5 .and. rh2 .le. 100) then
               m = 21.0606 + (0.005565*rh2**2) - (0.00035*rh2*t2f) - &
                   (0.483199*rh2)
            else         !   rh2 > 100
               m = 21.0606 + (0.005565*100**2) - (0.00035*100*t2f) - &
                   (0.483199*100)
               print *, 'fwi calculation: rh2 > 100 ', j, i, rh2, t2f
            end if

            n = 1.0 - 2.0*(m/30.) + 1.5*(m/30.)**2 - 0.5*(m/30.)**3
            fwi(i, j) = (n*sqrt(1.0 + uuu10**2 + vvv10**2))/0.3002

         end do
      end do
      return
   end subroutine fosberg_fwi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module fire
