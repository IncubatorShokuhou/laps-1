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

  subroutine integrated_liquid(nx, ny, nz, cond_mr, vapor_mr, rho, height, topo, &
                               intliqwater, totpcpwater)

     ! computes integrated liquid water and total precip. water in a column.
     !
     !  adapted from usaf weather agency mm5 post processor
     !
     !  brent shaw, noaa forecast systems lab, dec 2000

     implicit none
     ! arguments

     integer, intent(in)                 :: nx
     integer, intent(in)                 :: ny
     integer, intent(in)                 :: nz
     real, intent(in)                 :: cond_mr(nx, ny, nz)
     real, intent(in)                 :: vapor_mr(nx, ny, nz)
     real, intent(in)                 :: rho(nx, ny, nz)
     real, intent(in)                 :: height(nx, ny, nz)
     real, intent(in)                 :: topo(nx, ny)
     real, intent(out)                :: intliqwater(nx, ny)
     real, intent(out)                :: totpcpwater(nx, ny)

     ! locals
     real                                :: height_top, height_bot
     real                                :: dz
     integer                             :: i, j, k

     do j = 1, ny
        do i = 1, nx
           intliqwater(i, j) = 0.0
           totpcpwater(i, j) = 0.0
           do k = 1, nz
              ! compute layer thickness
              if (k .eq. 1) then
                 height_bot = topo(i, j)
                 height_top = 0.5*(height(i, j, 1) + height(i, j, 2))
              else if (k .eq. nz) then
                 height_bot = 0.5*(height(i, j, nz - 1) + height(i, j, nz))
                 height_top = 2*height(i, j, nz) - height_bot
              else
                 height_bot = 0.5*(height(i, j, k - 1) + height(i, j, k))
                 height_top = 0.5*(height(i, j, k) + height(i, j, k + 1))
              end if
              dz = height_top - height_bot
              intliqwater(i, j) = intliqwater(i, j) + cond_mr(i, j, k)*rho(i, j, k)*dz
              totpcpwater(i, j) = totpcpwater(i, j) + vapor_mr(i, j, k)*rho(i, j, k)*dz
           end do
        end do
     end do
     return
  end subroutine integrated_liquid

