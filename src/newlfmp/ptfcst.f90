subroutine lfm_initpts(lfmprd_dir, a24time, a9time)

!gfortran additions begin
   include 'trigd.inc'
!gfortran additions end
   use lfmgrid

   implicit none

   integer :: npt, nestpt, nc, chr, i, m
   integer :: ptpcptype
   real :: ptlat, ptlon, pti, ptj, ptzsfc, pttsfc, pttdsfc, ptrhsfc &
           , ptusfc, ptvsfc, ptspd, ptdir, ptupbl, ptvpbl, ptbspd, ptbdir &
           , ptceiling, ptvis, ptcld, ptpcp, ptsnow, ptvnt, ptpbl
   character(len=256) :: lfmprd_dir, ptfile, outfile
   character(len=80) :: customer
   character(len=32) :: cmtype
   character(len=24) :: a24time
   character(len=16) :: date_str
   character(len=10) :: ptid
   character(len=9) :: a9time
   character(len=8) :: ptwx
   character(len=3) :: domnum_str
   character(len=2) :: amonth
   logical :: back = .true., there

   character(len=3) :: amonths(12)
   data amonths/'jan', 'feb', 'mar', 'apr', 'may', 'jun' &
      , 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/

! prepare date/time string.

   do m = 1, 12
      if (a24time(4:6) == amonths(m)) exit
   end do
   write (amonth, '(i2.2)') m

   date_str = amonth//'/'//a24time(1:2)//'/'//a24time(8:11)//' '//a24time(13:17)

! read point forecast file.
! compute i,j locations of output point forecast locations in native grid.
! remove any stations outside of output domain.
! write point text file to model directory.

   nc = index(filename, '/', back)
   if (nc > 1) nc = index(filename(1:nc - 1), '/', back)
   ptfile = lfmprd_dir(1:nc)//'static/lfmpost_points.txt'
   inquire (file=trim(ptfile), exist=there)
   if (.not. there) then
      print *, ' '
      print *, 'could not find point forecast location file: ', trim(ptfile)
      npt = 0
      return
   end if

   if (verbose) then
      print *, ' '
      print *, 'reading points from: ', trim(ptfile)
   end if

   open (1, file=trim(ptfile), status='old', form='formatted', access='sequential')
   npt = 0
   do while (.true.)
      read (1, '(i2,1x,a10,1x,f8.4,1x,f9.4,1x,a)', end=1) nestpt, ptid, ptlat, ptlon, customer
      if (nestpt == domnum) then
         call latlon_to_ij(proj, ptlat, ptlon, pti, ptj)
         if (pti >= 1. .and. pti <= lx .and. &
             ptj >= 1. .and. pti <= ly) then
            npt = npt + 1
            call gdtost_lfm(zsfc, lx, ly, pti, ptj, ptzsfc)
            call gdtost_lfm(tsfc, lx, ly, pti, ptj, pttsfc)
            call gdtost_lfm(tdsfc, lx, ly, pti, ptj, pttdsfc)
            if (point_temp_units == 'f') then
               pttsfc = (pttsfc - 273.15)*1.8 + 32.
               pttdsfc = (pttdsfc - 273.15)*1.8 + 32.
            elseif (point_temp_units == 'c') then
               pttsfc = pttsfc - 273.15
               pttdsfc = pttdsfc - 273.15
            end if
            call gdtost_lfm(rhsfc, lx, ly, pti, ptj, ptrhsfc)
            call gdtost_lfm(usfc, lx, ly, pti, ptj, ptusfc)
            call gdtost_lfm(vsfc, lx, ly, pti, ptj, ptvsfc)
            ptspd = (ptusfc**2 + ptvsfc**2)**0.5
            ptdir = atan2d(-ptusfc, -ptvsfc)
            if (ptdir < 0.) ptdir = ptdir + 360.
            if (make_firewx) then
               call gdtost_lfm(upbl, lx, ly, pti, ptj, ptupbl)
               call gdtost_lfm(vpbl, lx, ly, pti, ptj, ptvpbl)
               ptbspd = (ptupbl**2 + ptvpbl**2)**0.5
               ptbdir = atan2d(-ptupbl, -ptvpbl)
               if (ptbdir < 0.) ptbdir = ptbdir + 360.
            else
               ptbspd = 0.
               ptbdir = 999.
            end if
            if (point_windspd_units == 'kts') then
               ptspd = ptspd*1.9425
               if (make_firewx) ptbspd = ptbspd*1.9425
            elseif (point_windspd_units == 'mph') then
               ptspd = ptspd*2.2369
               if (make_firewx) ptbspd = ptbspd*2.2369
            end if
            call gdtost_lfm(ceiling, lx, ly, pti, ptj, ptceiling)
            if (ptceiling < 50000.) then
               ptceiling = ptceiling*.032808  ! hundreds of feet
            else
               ptceiling = 999.
            end if
            call gdtost_lfm(visibility, lx, ly, pti, ptj, ptvis)
            ptvis = ptvis*0.00062317  ! miles
            call gdtost_lfm(cldamt, lx, ly, pti, ptj, ptcld)
            call gdtost_lfm(pcp_inc, lx, ly, pti, ptj, ptpcp)
            ptpcp = ptpcp*39.27
            call gdtost_lfm(snow_inc, lx, ly, pti, ptj, ptsnow)
            ptsnow = ptsnow*39.27
            call gdtost_lfm(vnt_index, lx, ly, pti, ptj, ptvnt)
            ptvnt = max(0., ptvnt)
            if (point_vent_units == 'kt-ft') ptvnt = ptvnt*6.3774
            call gdtost_lfm(pblhgt, lx, ly, pti, ptj, ptpbl)
            ptpbl = max(0., ptpbl)
            ptpbl = ptpbl*3.28

            ptpcptype = nint(pcptype_sfc(nint(pti), nint(ptj)))
            select case (ptpcptype)
            case (0)
               if (ptvis >= 7.) then
                  if (ptcld >= 0.75) then
                     ptwx = 'cloudy  '
                  elseif (ptcld < 0.75 .and. ptcld >= 0.25) then
                     ptwx = 'pt cldy '
                  else
                     ptwx = 'clear   '
                  end if
               elseif (ptvis < 7. .and. ptvis >= 3.) then
                  ptwx = 'haze    '
               else
                  ptwx = 'fog     '
               end if
            case (1)
               ptwx = 'rain    '
            case (2)
               ptwx = 'snow    '
            case (3)
               ptwx = 'frzrain '
            case (4)
               ptwx = 'sleet   '
            case (5)
               ptwx = 'hail    '
            case (6)
               ptwx = 'drizzle '
            case (7)
               ptwx = 'mixed   '
            case (8)
               ptwx = 'mixed   '
            case (9)
               ptwx = 'rain/ice'
            case default
               ptwx = 'unknown '
            end select

            write (domnum_str, '("d",i2.2)') domnum
            outfile = trim(lfmprd_dir)//'/'//domnum_str//'/points/'//trim(ptid)//'_'//a9time//'_fcst.txt'
            inquire (file=trim(outfile), exist=there)
            if (fcsttime == 0) there = .false.
            if (there) then
               open (2, file=trim(outfile), form='formatted', status='old', position='append')
            else
               cmtype = mtype
               do i = 1, len_trim(cmtype)
                  chr = ichar(cmtype(i:i))
                  if (chr > 96 .and. chr < 123) cmtype(i:i) = char(chr - 32)
               end do

               open (2, file=trim(outfile), form='formatted', status='new')
               write (2, '("**********************************************************************************************")')
               write (2, '("location: ",a,2x,"lat: ",f8.4,2x,"lon: ",f9.4,2x,"i: ",f7.2,2x,"j: ",f7.2)') &
                  ptid, ptlat, ptlon, pti, ptj
               write (2, '(a4,f4.1," km    forecast cycle: ",a,2x,"dom: ",i2,2x,"model elevation: ",i4)') &
                  cmtype(1:4), grid_spacing/1000., a9time, domnum, nint(ptzsfc)
               write (2, '("**********************************************************************************************")')
               write (2, '("date       time  tmp dpt rh  wind   cei vis  weather  precp snow vent   mixht pblwnd hm hh fbg")')
               write (2, '(a3,8x,a3,3x,a1,3x,a1,3x,"%",3x,"dg@",a3,1x,"hft mile",10x,"in",4x,"in",3x,a5,2x,"ftagl",1x,"dg@",a3)') &
                  point_tz_label, point_tz_label, point_temp_units &
                  , point_temp_units, point_windspd_units, point_vent_units, point_windspd_units
               write (2, '("---------- ----- --- --- --- ------ --- ---- -------- ----- ---- ------ ----- ------ -- -- ---")')
            end if

         write (2, '(a,3i4,1x,i3.3,"/",i2.2,1x,i3.3,1x,f4.1,1x,a,1x,f5.2,1x,f4.1,1x,i6,1x,i5,1x,i3.3,"/",i2.2,1x,i2,1x,i2,1x,i3)') &
               date_str, nint(pttsfc), nint(pttdsfc), nint(ptrhsfc), nint(ptdir/10.)*10, nint(ptspd) &
               , nint(ptceiling), ptvis, ptwx, ptpcp, ptsnow, nint(ptvnt), nint(ptpbl) &
               , nint(ptbdir/10.)*10, nint(ptbspd), nint(ham_index(nint(pti), nint(ptj))) &
               , nint(hah_index(nint(pti), nint(ptj))), nint(fwi_index(nint(pti), nint(ptj)))
            close (2)
         end if
      end if
   end do
1  continue

   return
end
