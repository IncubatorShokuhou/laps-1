      subroutine get_static_field_interp(ctype,i4time,nx,ny,data
     1,istatus)
   
! returns a time-interpolated (valid for time) 2d array of static climatological
! albedo using the monthly values in the static file.  the monthly values are 
! valid on the 15th day of each month.  this routine only interpolates to the
! nearest day and does not account for leap years, but this should not be any
! big deal.  
!
! j.smart 4-02: subroutine taken from wrfsi software (module_wrfsi_static.f)
!               and modified for laps use to get albedo.

      implicit none
      character(len=9)    :: ctime9
      character           :: ctype*(*)

      integer             :: nx,ny
      real                :: data(nx,ny)
    
      integer             :: midmonth_day(12)
      integer             :: valid_day
      integer             :: yyyyjjj
      integer             :: istatus
      real                :: sss
      integer             :: i4time
      integer             :: m, d1, d2, m1, m2
      character(len=3)    :: var_2d
      real, allocatable   :: data1(:,:),data2(:,:)
      real                :: w1, w2
      integer             :: int_file(9)
      integer             :: i
      integer             :: jday
!     integer             :: nyear,jday,nhour,min,month,nday

! midmonth_day is the julian day of the year corresponding to the 15th day
! of each month for a standard (non-leap) year

c     istatus = 0

      data midmonth_day / 15, 43, 74, 105, 135, 166, 196,
     &  227, 258, 288, 319, 349 /

! convert date string into integer yyyyjjj and sss
      call make_fnam_lp(i4time,ctime9,istatus)
      if(istatus.ne.1)then
         print*,'error returned: make_fnam_lp'
         return
      endif
      read(ctime9(3:5),'(i3)')valid_day

      print*,'time-interp monthly static field to day: ',valid_day

    ! find bounding months
      if ((valid_day .lt. midmonth_day(1)) .or.
     &    (valid_day .gt. midmonth_day(12))) then
    ! december and january are bounding months
           d1 = midmonth_day(12)
           d2 = midmonth_day(1)
           m1 = 12
           m2 = 1
      else
        find_bounds: do m = 1, 11
           d1 = midmonth_day(m)
           d2 = midmonth_day(m+1)
           if (valid_day .eq. d1) then
                d2 = d1
                m1 = m
                m2 = m1
                exit find_bounds
           else if (valid_day .eq. d2) then
                d1 = d2
                m1 = m + 1
                m2 = m1
                exit find_bounds
           else if ((valid_day .gt. d1).and.(valid_day .lt. d2)) then
                m1 = m
                m2 = m + 1
                exit find_bounds
           endif
        enddo find_bounds
      endif

! if d1 = d2, then we don't need any interpolation, just get that month's 
! data values
      if ( d1 .eq. d2) then
           if(ctype.eq.'albedo')then
              write(var_2d, '("a",i2.2)') m1
           elseif(ctype.eq.'green')then
              write(var_2d, '("g",i2.2)') m1
           endif
           call read_static_grid(nx,ny,var_2d,data,istatus)
           if(istatus.ne.1)then
c             print*,' error reading laps static: ',var_2d
              return
           endif
      else
! we need to get the two months of bounding data and time interpolate
         allocate(data1 (nx,ny))
         allocate(data2 (nx,ny))
         if(ctype.eq.'albedo')then
            write(var_2d, '("a",i2.2)') m1
         elseif(ctype.eq.'green')then
            write(var_2d, '("g",i2.2)') m1
         endif
         call read_static_grid(nx,ny,var_2d,data1,istatus)
         if(istatus .ne. 1)then
c           print*,' error reading laps static: ',var_2d
            return
         endif
         if(ctype.eq.'albedo')then
            write(var_2d, '("a",i2.2)') m2
         elseif(ctype.eq.'green')then
            write(var_2d, '("g",i2.2)') m2
         endif
         call read_static_grid(nx,ny,var_2d,data2,istatus)
         if(istatus .ne. 1)then
c           print*,' error reading laps static: ',var_2d
            return
         endif

!       compute weights
        if (d2 .gt. d1) then
            w1 = ( float(d2) - float(valid_day) ) / float(d2-d1)
        else ! we must be between dec 15 and jan 15
            if (valid_day .lt. midmonth_day(1)) then ! we are in january
                w1 = ( float(d2) - float(valid_day) ) / 31.
            else ! we are in december
                w1 = ( 366. - float(valid_day) +
     &                        float(midmonth_day(1))) / 31.
            endif
        endif
        w2 = 1. - w1
        data = w1*data1 + w2*data2
        deallocate(data1)
        deallocate(data2)
      endif

c     istatus=1
      return
      end subroutine get_static_field_interp
