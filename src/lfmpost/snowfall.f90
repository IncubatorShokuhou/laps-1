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

 subroutine snowfall(tsfc, prcpinc, preciptype, imax, jmax, &
                     snowinc, snowtot)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--
!--   name: snow accumulation algorithm
!--
!--   purpose
!--   =======
!--   this algorithm calculates incremental snow accumulation and
!--   total snow accumulation.
!--
!--   method
!--   ======
!--   - if precip type is snow
!--      - calculate surface temperature on dot
!--      - calculate incremental precip on dot
!--      - if surface temp is >= 10 f
!--        - use a 10:1 snow/liquid ratio
!--      - if surface temp is < 10 f
!--        - use a 15:1 snow/liquid ratio
!--   - if precip type is not snow
!--      - snow accumulation is zero
!--
!--   variables
!--   =========
!--   name             type      i/o     definition
!--   ----             ----      ---     ----------
!--   imax             integer  input   grid dimension, i direction
!--   jmax             integer  input   grid dimension, j direction
!--   prcpinc          real     input   3hr precip accum
!--                                     array(2-d)
!--   preciptype       integer  input   0 - is no precip
!--                                     1 - is rain
!--                                     2 - is freezing rain
!--                                     3 - is ice/mixed
!--                                     4 - is snow
!--   snowinc          real     output  3hr snow accum
!--                                         array(2-d)
!--   snowtot          real     output  total snow accum
!--                                        array(2-d)
!--   tsfc             real     input   surface temp array (2-d)
!--   tsfcf            real     local   surface temp in f
!--
!--   updates
!--   =======
!--   20 feb 98  initial version..................capt. john lewis/dnxt
!--   10 nov 98  changed preciptype flag for snow
!--              from 4 to 5; this corresponds
!--              with changes made to the precip
!--              type algorithm (wintprec.f)...capt david beberwyk/dnxt
!--    5 jan 99  removed interpolations to dot grid as mmpost now
!--              operates on the cross grid........................dnxm
!--    4 jan 01  adapted by fsl for use with laps.. b. shaw, noaa/fsl
!--
!----------------------------------------------------------------------
!----------------------------------------------------------------------

    use constants
    implicit none

    real, external        :: fahren
    integer                     :: i
    integer, intent(in)      :: imax
    integer                     :: j
    integer, intent(in)      :: jmax
    real, intent(in)      :: prcpinc(imax, jmax)
    integer, intent(in)      :: preciptype(imax, jmax)
    real, intent(out)     :: snowinc(imax, jmax)
    real, intent(inout)   :: snowtot(imax, jmax)
    real, intent(in)      :: tsfc(imax, jmax)
    real                        :: tsfcf

!--------------------------------------------------------------------------
!     -- begin the main double-do-loop
!-------------------------------------------------------------------------

    do j = 1, jmax
       do i = 1, imax

!-----------------------------------------------------------------------
!--       check if precipitation type is snow.
!-----------------------------------------------------------------------

          if (preciptype(i, j) == 5) then

             tsfcf = fahren(tsfc(i, j) - t0)

!-----------------------------------------------------------------------
!--         liquid equivalent of snow depends of surface temperature.
!-----------------------------------------------------------------------

             if (tsfcf >= 10.0) then

                snowinc(i, j) = 10.0*prcpinc(i, j)

             else

                snowinc(i, j) = 15.0*prcpinc(i, j)

             end if

!-----------------------------------------------------------------------
!--         if precip type is not snow then snow accum is zero.
!-----------------------------------------------------------------------

          else

             snowinc(i, j) = 0.0

          end if

       end do
    end do

!-------------------------------------------------------------------
!--   update snow total
!-------------------------------------------------------------------

    snowtot = snowtot + snowinc

 end subroutine snowfall
