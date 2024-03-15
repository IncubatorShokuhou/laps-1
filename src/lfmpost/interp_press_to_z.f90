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

  subroutine interp_press_to_z(press_levels, heights, new_level, press_out, &
                               nx, ny, nz)

     ! given a 1d array of pressure levels and a 3d array of heights (m) at
     ! those levels, this routine interpolates the pressure to the desired
     ! new_level.

     ! pressures are in pa, heights are in m!

     implicit none
     integer, intent(in)               :: nx, ny, nz
     real, intent(in)                  :: press_levels(nz)
     real, intent(in)                  :: heights(nx, ny, nz)
     real, intent(in)                  :: new_level
     real, intent(out)                 :: press_out(nx, ny)
     integer                           :: i, j, k, ktop
     real                              :: ptop, pbot, ztop, zbot

     nsloop: do j = 1, ny
        ewloop: do i = 1, nx
           vertloop: do k = 2, nz
              if (heights(nx, ny, nz) .ge. new_level) then
                 ktop = k
                 exit vertloop
              end if
           end do vertloop
           ztop = heights(i, j, k)
           zbot = heights(i, j, k - 1)
           ptop = press_levels(k)*0.01    ! convert to mb!
           pbot = press_levels(k - 1)*0.01
           press_out(i, j) = exp( &
                             (new_level*alog(pbot/ptop) - &
                              ztop*alog(pbot) + &
                              zbot*alog(ptop))/ &
                             (zbot - ztop))*100.  ! convert back to pa
        end do ewloop
     end do nsloop
     return
  end subroutine interp_press_to_z

