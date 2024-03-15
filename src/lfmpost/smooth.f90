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

subroutine smooth(field, ix, iy, iz, smth)

!
! *** subprogram:  smooth - smooth a meteorological field.
!     author:  stan benjamin
!     date  :  90-06-15
!
! *** abstract:  shapiro smoother.
!
! *** program history log:
!        85-12-09  s. benjamin - original version
!        96-06-16  j. snook    - modified to do 3d rams fields
!                              - hold array is dynamically allocated
!        01-03-07  b. shaw     - adapted to free-form f90 as part o
!                                mm5post. includes use of array syntax.
!
! *** usage:  call smooth(field,ix,iy,iz,smth)
!
! *** input argument list:
!        field    - real array  field(ix,iy,iz)
!                               meteorological field
!        ix       - integer     x coordinates of field
!        iy       - integer     y coordinates of field
!        iz       - integer     z coordinates of field
!        smth     - real
!
! *** output argument list:
!        field    - real array  field(ix,iy,iz)
!                               smoothed meteorological field
!
! *** remarks:  reference:  shapiro, 1970: "smoothing, filtering, and
!        boundary effects", rev. geophys. sp. phys., 359-387.
!
!     this filter is of the type
!        z(i) = (1-s)z(i) + s(z(i+1)+z(i-1))/2
!     for a filter which is supposed to damp 2dx waves completely
!     but leave 4dx and longer with little damping,
!     it should be run with 2 passes using smth (or s) of 0.5
!     and -0.5
!----------------------------------------------------------------------------
   implicit none

   ! argument list

   integer, intent(in)         :: ix, iy, iz
   real, intent(inout)         :: field(ix, iy, iz)
   real, intent(in)            :: smth

   ! local variables
   integer                     :: i, j, k, i1, i2, it
   real                        :: hold(ix, 2)
   real                        :: smth1, smth2, smth3, smth4, &
                                  smth5, sum1, sum2
!----------------------------------------------------------------------------
   smth1 = 0.25*smth*smth
   smth2 = 0.50*smth*(1.-smth)
   smth3 = (1.-smth)*(1.-smth)
   smth4 = (1.-smth)
   smth5 = 0.5*smth

   do k = 1, iz

      hold(:, 1:2) = 0.
      i1 = 2
      i2 = 1

      do j = 2, iy - 1
         it = i1
         i1 = i2
         i2 = it
         do i = 2, ix - 1
            sum1 = field(i - 1, j + 1, k) + field(i - 1, j - 1, k) &
                   + field(i + 1, j + 1, k) + field(i + 1, j - 1, k)
            sum2 = field(i, j + 1, k) + field(i + 1, j, k) &
                   + field(i, j - 1, k) + field(i - 1, j, k)
            hold(i, i1) = smth1*sum1 + smth2*sum2 + smth3*field(i, j, k)
         end do
         if (j .ne. 2) field(2:ix - 1, j - 1, k) = hold(2:ix - 1, i2)
      end do

      field(2:ix - 1, iy - 1, k) = hold(2:ix - 1, i1)

      field(2:ix - 1, 1, k) = smth4*field(2:ix - 1, 1, k) + &
                              smth5*(field(1:ix - 2, 1, k) + field(3:ix, 1, k))
      field(2:ix - 1, iy, k) = smth4*field(2:ix - 1, iy, k) + &
                               smth5*(field(1:ix - 2, iy, k) + field(3:ix, iy, k))

      field(1, 2:iy - 1, k) = smth4*field(1, 2:iy - 1, k) + &
                              smth5*(field(1, 1:iy - 2, k) + field(1, 3:iy, k))

   end do
   return
end subroutine smooth
