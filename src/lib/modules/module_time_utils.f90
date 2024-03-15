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


module time_utils

  implicit none
 
  integer,private  :: days_in_month(12)
  data days_in_month / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine geth_idts (ndate, odate, idts)
   
      implicit none
      
      !  from 2 input mdates ('yyyy-mm-dd hh:mm:ss.ffff'), 
      !  compute the time difference.
      
      !  on entry     -  ndate  -  the new hdate.
      !                  odate  -  the old hdate.
      
      !  on exit      -  idts    -  the change in time in seconds.
      
      character (len=*) , intent(inout) :: ndate, odate
      integer           , intent(out)   :: idts
      
      !  local variables
      
      !  yrnew    -  indicates the year associated with "ndate"
      !  yrold    -  indicates the year associated with "odate"
      !  monew    -  indicates the month associated with "ndate"
      !  moold    -  indicates the month associated with "odate"
      !  dynew    -  indicates the day associated with "ndate"
      !  dyold    -  indicates the day associated with "odate"
      !  hrnew    -  indicates the hour associated with "ndate"
      !  hrold    -  indicates the hour associated with "odate"
      !  minew    -  indicates the minute associated with "ndate"
      !  miold    -  indicates the minute associated with "odate"
      !  scnew    -  indicates the second associated with "ndate"
      !  scold    -  indicates the second associated with "odate"
      !  i        -  loop counter
      !  mday     -  a list assigning the number of days in each month
      
      character (len=24) :: tdate
      integer :: olen, nlen
      integer :: yrnew, monew, dynew, hrnew, minew, scnew
      integer :: yrold, moold, dyold, hrold, miold, scold
      integer :: mday(12), i, newdys, olddys
      logical :: npass, opass
      integer :: isign
      
      if (odate.gt.ndate) then
         isign = -1
         tdate=ndate
         ndate=odate
         odate=tdate
      else
         isign = 1
      end if
      
      !  assign the number of days in a months
      
      mday( 1) = 31
      mday( 2) = 28
      mday( 3) = 31
      mday( 4) = 30
      mday( 5) = 31
      mday( 6) = 30
      mday( 7) = 31
      mday( 8) = 31
      mday( 9) = 30
      mday(10) = 31
      mday(11) = 30
      mday(12) = 31
      
      !  break down old hdate into parts
      
      hrold = 0
      miold = 0
      scold = 0
      olen = len(odate)
      
      read(odate(1:4),  '(i4)') yrold
      read(odate(6:7),  '(i2)') moold
      read(odate(9:10), '(i2)') dyold
      if (olen.ge.13) then
         read(odate(12:13),'(i2)') hrold
         if (olen.ge.16) then
            read(odate(15:16),'(i2)') miold
            if (olen.ge.19) then
               read(odate(18:19),'(i2)') scold
            end if
         end if
      end if
      
      !  break down new hdate into parts
      
      hrnew = 0
      minew = 0
      scnew = 0
      nlen = len(ndate)
      
      read(ndate(1:4),  '(i4)') yrnew
      read(ndate(6:7),  '(i2)') monew
      read(ndate(9:10), '(i2)') dynew
      if (nlen.ge.13) then
         read(ndate(12:13),'(i2)') hrnew
         if (nlen.ge.16) then
            read(ndate(15:16),'(i2)') minew
            if (nlen.ge.19) then
               read(ndate(18:19),'(i2)') scnew
            end if
         end if
      end if
      
      !  check that the dates make sense.
      
      npass = .true.
      opass = .true.
      
      !  check that the month of ndate makes sense.
      
      if ((monew.gt.12).or.(monew.lt.1)) then
         print*, 'geth_idts:  month of ndate = ', monew
         npass = .false.
      end if
      
      !  check that the month of odate makes sense.
      
      if ((moold.gt.12).or.(moold.lt.1)) then
         print*, 'geth_idts:  month of odate = ', moold
         opass = .false.
      end if
      
      !  check that the day of ndate makes sense.
      
      if (monew.ne.2) then
      ! ...... for all months but february
         if ((dynew.gt.mday(monew)).or.(dynew.lt.1)) then
            print*, 'geth_idts:  day of ndate = ', dynew
            npass = .false.
         end if
      else if (monew.eq.2) then
      ! ...... for february
         if ((dynew.gt.nfeb(yrnew)).or.(dynew.lt.1)) then
            print*, 'geth_idts:  day of ndate = ', dynew
            npass = .false.
         end if
      end if
      
      !  check that the day of odate makes sense.
      
      if (moold.ne.2) then
      ! ...... for all months but february
         if ((dyold.gt.mday(moold)).or.(dyold.lt.1)) then
            print*, 'geth_idts:  day of odate = ', dyold
            opass = .false.
         end if
      else if (moold.eq.2) then
      ! ....... for february
         if ((dyold.gt.nfeb(yrold)).or.(dyold.lt.1)) then
            print*, 'geth_idts:  day of odate = ', dyold
            opass = .false.
         end if
      end if
      
      !  check that the hour of ndate makes sense.
      
      if ((hrnew.gt.23).or.(hrnew.lt.0)) then
         print*, 'geth_idts:  hour of ndate = ', hrnew
         npass = .false.
      end if
      
      !  check that the hour of odate makes sense.
      
      if ((hrold.gt.23).or.(hrold.lt.0)) then
         print*, 'geth_idts:  hour of odate = ', hrold
         opass = .false.
      end if
      
      !  check that the minute of ndate makes sense.
      
      if ((minew.gt.59).or.(minew.lt.0)) then
         print*, 'geth_idts:  minute of ndate = ', minew
         npass = .false.
      end if
      
      !  check that the minute of odate makes sense.
      
      if ((miold.gt.59).or.(miold.lt.0)) then
         print*, 'geth_idts:  minute of odate = ', miold
         opass = .false.
      end if
      
      !  check that the second of ndate makes sense.
      
      if ((scnew.gt.59).or.(scnew.lt.0)) then
         print*, 'geth_idts:  second of ndate = ', scnew
         npass = .false.
      end if
      
      !  check that the second of odate makes sense.
      
      if ((scold.gt.59).or.(scold.lt.0)) then
         print*, 'geth_idts:  second of odate = ', scold
         opass = .false.
      end if
      
      if (.not. npass) then
         print*, 'screwy ndate: ', ndate(1:nlen)
         stop 'ndate_2'
      end if
      
      if (.not. opass) then
         print*, 'screwy odate: ', odate(1:olen)
         stop 'odate_1'
      end if
      
      !  date checks are completed.  continue.
      
      !  compute number of days from 1 january odate, 00:00:00 until ndate
      !  compute number of hours from 1 january odate, 00:00:00 until ndate
      !  compute number of minutes from 1 january odate, 00:00:00 until ndate
      
      newdys = 0
      do i = yrold, yrnew - 1
         newdys = newdys + (365 + (nfeb(i)-28))
      end do
      
      if (monew .gt. 1) then
         mday(2) = nfeb(yrnew)
         do i = 1, monew - 1
            newdys = newdys + mday(i)
         end do
         mday(2) = 28
      end if
      
      newdys = newdys + dynew-1
      
      !  compute number of hours from 1 january odate, 00:00:00 until odate
      !  compute number of minutes from 1 january odate, 00:00:00 until odate
      
      olddys = 0
      
      if (moold .gt. 1) then
         mday(2) = nfeb(yrold)
         do i = 1, moold - 1
            olddys = olddys + mday(i)
         end do
         mday(2) = 28
      end if
      
      olddys = olddys + dyold-1
      
      !  determine the time difference in seconds
      
      idts = (newdys - olddys) * 86400
      idts = idts + (hrnew - hrold) * 3600
      idts = idts + (minew - miold) * 60
      idts = idts + (scnew - scold)
      
      if (isign .eq. -1) then
         tdate=ndate
         ndate=odate
         odate=tdate
         idts = idts * isign
      end if
   
   end subroutine geth_idts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine geth_newdate (ndate, odate, idt)
   
      implicit none
      
      !  from old date ('yyyy-mm-dd hh:mm:ss.ffff') and 
      !  delta-time, compute the new date.
   
      !  on entry     -  odate  -  the old hdate.
      !                  idt    -  the change in time
   
      !  on exit      -  ndate  -  the new hdate.
      
      integer , intent(in)           :: idt
      character (len=*) , intent(out) :: ndate
      character (len=*) , intent(in)  :: odate
      
       
      !  local variables
       
      !  yrold    -  indicates the year associated with "odate"
      !  moold    -  indicates the month associated with "odate"
      !  dyold    -  indicates the day associated with "odate"
      !  hrold    -  indicates the hour associated with "odate"
      !  miold    -  indicates the minute associated with "odate"
      !  scold    -  indicates the second associated with "odate"
       
      !  yrnew    -  indicates the year associated with "ndate"
      !  monew    -  indicates the month associated with "ndate"
      !  dynew    -  indicates the day associated with "ndate"
      !  hrnew    -  indicates the hour associated with "ndate"
      !  minew    -  indicates the minute associated with "ndate"
      !  scnew    -  indicates the second associated with "ndate"
       
      !  mday     -  a list assigning the number of days in each month
      
      !  i        -  loop counter
      !  nday     -  the integer number of days represented by "idt"
      !  nhour    -  the integer number of hours in "idt" after taking out
      !              all the whole days
      !  nmin     -  the integer number of minutes in "idt" after taking out
      !              all the whole days and whole hours.
      !  nsec     -  the integer number of minutes in "idt" after taking out
      !              all the whole days, whole hours, and whole minutes.
       
      integer :: nlen, olen
      integer :: yrnew, monew, dynew, hrnew, minew, scnew, frnew
      integer :: yrold, moold, dyold, hrold, miold, scold, frold
      integer :: mday(12), nday, nhour, nmin, nsec, nfrac, i, ifrc
      logical :: opass
      character (len=10) :: hfrc
      character (len=1) :: sp
      ! integer, external :: nfeb  ! in the same module now
      
      !  assign the number of days in a months
      
      mday( 1) = 31
      mday( 2) = 28
      mday( 3) = 31
      mday( 4) = 30
      mday( 5) = 31
      mday( 6) = 30
      mday( 7) = 31
      mday( 8) = 31
      mday( 9) = 30
      mday(10) = 31
      mday(11) = 30
      mday(12) = 31
      
      !  break down old hdate into parts
      
      hrold = 0
      miold = 0
      scold = 0
      frold = 0
      olen = len(odate)
      if (olen.ge.11) then
         sp = odate(11:11)
      else
         sp = ' '
      end if
      
      !  use internal read statements to convert the character string
      !  date into integer components.
   
      read(odate(1:4),  '(i4)') yrold
      read(odate(6:7),  '(i2)') moold
      read(odate(9:10), '(i2)') dyold
      if (olen.ge.13) then
         read(odate(12:13),'(i2)') hrold
         if (olen.ge.16) then
            read(odate(15:16),'(i2)') miold
            if (olen.ge.19) then
               read(odate(18:19),'(i2)') scold
               if (olen.gt.20) then
                  read(odate(21:olen),'(i2)') frold
               end if
            end if
         end if
      end if
      
      !  set the number of days in february for that year.
      
      mday(2) = nfeb(yrold)
      
      !  check that odate makes sense.
      
      opass = .true.
      
      !  check that the month of odate makes sense.
      
      if ((moold.gt.12).or.(moold.lt.1)) then
         write(*,*) 'geth_newdate:  month of odate = ', moold
         opass = .false.
      end if
      
      !  check that the day of odate makes sense.
      
      if ((dyold.gt.mday(moold)).or.(dyold.lt.1)) then
         write(*,*) 'geth_newdate:  day of odate = ', dyold
         opass = .false.
      end if
      
      !  check that the hour of odate makes sense.
      
      if ((hrold.gt.23).or.(hrold.lt.0)) then
         write(*,*) 'geth_newdate:  hour of odate = ', hrold
         opass = .false.
      end if
      
      !  check that the minute of odate makes sense.
      
      if ((miold.gt.59).or.(miold.lt.0)) then
         write(*,*) 'geth_newdate:  minute of odate = ', miold
         opass = .false.
      end if
      
      !  check that the second of odate makes sense.
      
      if ((scold.gt.59).or.(scold.lt.0)) then
         write(*,*) 'geth_newdate:  second of odate = ', scold
         opass = .false.
      end if
      
      !  check that the fractional part  of odate makes sense.
      
      
      if (.not.opass) then
         write(*,*) 'geth_newdate: crazy odate: ', odate(1:olen), olen
         stop 'odate_3'
      end if
      
      !  date checks are completed.  continue.
      
      
      !  compute the number of days, hours, minutes, and seconds in idt
      
      if (olen.gt.20) then !idt should be in fractions of seconds
         ifrc = olen-20
         ifrc = 10**ifrc
         nday   = abs(idt)/(86400*ifrc)
         nhour  = mod(abs(idt),86400*ifrc)/(3600*ifrc)
         nmin   = mod(abs(idt),3600*ifrc)/(60*ifrc)
         nsec   = mod(abs(idt),60*ifrc)/(ifrc)
         nfrac = mod(abs(idt), ifrc)
      else if (olen.eq.19) then  !idt should be in seconds
         ifrc = 1
         nday   = abs(idt)/86400 ! integer number of days in delta-time
         nhour  = mod(abs(idt),86400)/3600
         nmin   = mod(abs(idt),3600)/60
         nsec   = mod(abs(idt),60)
         nfrac  = 0
      else if (olen.eq.16) then !idt should be in minutes
         ifrc = 1
         nday   = abs(idt)/1440 ! integer number of days in delta-time
         nhour  = mod(abs(idt),1440)/60
         nmin   = mod(abs(idt),60)
         nsec   = 0
         nfrac  = 0
      else if (olen.eq.13) then !idt should be in hours
         ifrc = 1
         nday   = abs(idt)/24 ! integer number of days in delta-time
         nhour  = mod(abs(idt),24)
         nmin   = 0
         nsec   = 0
         nfrac  = 0
      else if (olen.eq.10) then !idt should be in days
         ifrc = 1
         nday   = abs(idt)/24 ! integer number of days in delta-time
         nhour  = 0
         nmin   = 0
         nsec   = 0
         nfrac  = 0
      else
         write(*,'(''geth_newdate: strange length for odate: '', i3)') &
              olen
         write(*,*) odate(1:olen)
         stop 'odate_4'
      end if
      
      if (idt.ge.0) then
      
         frnew = frold + nfrac
         if (frnew.ge.ifrc) then
            frnew = frnew - ifrc
            nsec = nsec + 1
         end if
      
         scnew = scold + nsec
         if (scnew .ge. 60) then
            scnew = scnew - 60
            nmin  = nmin + 1
         end if
      
         minew = miold + nmin
         if (minew .ge. 60) then
            minew = minew - 60
            nhour  = nhour + 1
         end if
      
         hrnew = hrold + nhour
         if (hrnew .ge. 24) then
            hrnew = hrnew - 24
            nday  = nday + 1
         end if
      
         dynew = dyold
         monew = moold
         yrnew = yrold
         do i = 1, nday
            dynew = dynew + 1
            if (dynew.gt.mday(monew)) then
               dynew = dynew - mday(monew)
               monew = monew + 1
               if (monew .gt. 12) then
                  monew = 1
                  yrnew = yrnew + 1
                  ! if the year changes, recompute the number of days in february
                  mday(2) = nfeb(yrnew)
               end if
            end if
         end do
      
      else if (idt.lt.0) then
      
         frnew = frold - nfrac
         if (frnew .lt. 0) then
            frnew = frnew + ifrc
            nsec = nsec - 1
         end if
      
         scnew = scold - nsec
         if (scnew .lt. 00) then
            scnew = scnew + 60
            nmin  = nmin + 1
         end if
      
         minew = miold - nmin
         if (minew .lt. 00) then
            minew = minew + 60
            nhour  = nhour + 1
         end if
      
         hrnew = hrold - nhour
         if (hrnew .lt. 00) then
            hrnew = hrnew + 24
            nday  = nday + 1
         end if
      
         dynew = dyold
         monew = moold
         yrnew = yrold
         do i = 1, nday
            dynew = dynew - 1
            if (dynew.eq.0) then
               monew = monew - 1
               if (monew.eq.0) then
                  monew = 12
                  yrnew = yrnew - 1
                  ! if the year changes, recompute the number of days in february
                  mday(2) = nfeb(yrnew)
               end if
               dynew = mday(monew)
            end if
         end do
      end if
      
      !  now construct the new mdate
      
      nlen = len(ndate)
      
      if (nlen.gt.20) then
         write(ndate(1:19),19) yrnew, monew, dynew, hrnew, minew, scnew
         write(hfrc,'(i10)') frnew+1000000000
         ndate = ndate(1:19)//'.'//hfrc(31-nlen:10)
      
      else if (nlen.eq.19.or.nlen.eq.20) then
         write(ndate(1:19),19) yrnew, monew, dynew, hrnew, minew, scnew
      19   format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':',i2.2,':',i2.2)
         if (nlen.eq.20) ndate = ndate(1:19)//'.'
      
      else if (nlen.eq.16) then
         write(ndate,16) yrnew, monew, dynew, hrnew, minew
      16   format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':',i2.2)
      
      else if (nlen.eq.13) then
         write(ndate,13) yrnew, monew, dynew, hrnew
      13   format(i4,'-',i2.2,'-',i2.2,'_',i2.2)
      
      else if (nlen.eq.10) then
         write(ndate,10) yrnew, monew, dynew
      10   format(i4,'-',i2.2,'-',i2.2)
      
      end if
      
      if (olen.ge.11) ndate(11:11) = sp
   
   end subroutine geth_newdate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function nfeb ( year ) result (num_days)
   
      ! compute the number of days in february for the given year
   
      implicit none
   
      integer :: year
      integer :: num_days
   
      num_days = 28 ! by default, february has 28 days ...
      if (mod(year,4).eq.0) then  
         num_days = 29  ! but every four years, it has 29 days ...
         if (mod(year,100).eq.0) then
            num_days = 28  ! except every 100 years, when it has 28 days ...
            if (mod(year,400).eq.0) then
               num_days = 29  ! except every 400 years, when it has 29 days.
            end if
         end if
      end if
   
   end function nfeb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine split_date_char ( date , century_year , month , day , hour , minute , second )
     
      implicit none
   
      !  input data.
   
      character(len=19) , intent(in) :: date 
   
      !  output data.
   
      integer , intent(out) :: century_year , month , day , hour , minute , second
      
      read(date,fmt='(    i4.4)') century_year
      read(date,fmt='( 5x,i2.2)') month
      read(date,fmt='( 8x,i2.2)') day
      read(date,fmt='(11x,i2.2)') hour
      read(date,fmt='(14x,i2.2)') minute
      read(date,fmt='(17x,i2.2)') second
   
   end subroutine split_date_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function compute_day_of_year(century_year, month, day) result (day_of_year)

    implicit none
    integer                     :: century_year, month, day
    integer                     :: day_of_year
    integer                     :: m

    days_in_month(2) = nfeb(century_year)           
    if ( month .eq. 1) then
       day_of_year = day
    else  
       day_of_year = 0
       do m = 1, month-1
         day_of_year = day_of_year + days_in_month(m)
       enddo
       day_of_year = day_of_year + day
    endif
 
  end function compute_day_of_year
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function hms_to_wrf_time (hour,minute,second) result (seconds_utc)
    implicit none
    integer                     :: hour, minute, second
    real                        :: seconds_utc
    seconds_utc = float(hour*3600 + minute*60 + second)
  end function hms_to_wrf_time     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wrf_date_to_ymd (wrf_date, century_year,month,day)
  
    implicit none
    integer , intent(in)                :: wrf_date
    integer , intent(out)               :: century_year
    integer,  intent(out)               :: month
    integer,  intent(out)               :: day
    integer :: day_of_year, m, total_days

    century_year = wrf_date/1000
    day_of_year = mod(wrf_date,1000)
    days_in_month(2) = nfeb(century_year)
    total_days = 0
    month_loop: do m = 1,12
      total_days = total_days + days_in_month(m)
      if (total_days .lt. day_of_year) then
         cycle month_loop
      else
        month = m
        day = days_in_month(m) - (total_days-day_of_year)
        exit month_loop
      endif
    end do month_loop
  end subroutine wrf_date_to_ymd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wrf_time_to_hms(wrf_time, hour,minute,second)
    implicit none
    real, intent(in)                    :: wrf_time
    integer, intent(out)                :: hour, minute, second
    integer                             :: leftover_seconds
    hour = floor(wrf_time)/3600
    leftover_seconds = mod(floor(wrf_time),3600)
    minute = leftover_seconds/60           
    second = mod(leftover_seconds,60)
  end subroutine wrf_time_to_hms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wrf_to_mm5_date (wrf_date, wrf_time, mm5_date)

    implicit none
    integer , intent(in)                :: wrf_date
    real ,    intent(in)                :: wrf_time
    character (len=19), intent(out)     :: mm5_date
    integer :: century_year, month, day, hour, minute, second

    call wrf_date_to_ymd(wrf_date,century_year,month,day)
    call wrf_time_to_hms(wrf_time, hour, minute, second)
    write (mm5_date, 91) century_year,month,day,hour,minute,second
    91 format(i4.4,'-',i2.2,'-',i2.2,'_',i2.2,':',i2.2,':',i2.2) 

  end subroutine wrf_to_mm5_date

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mm5_to_wrf_date (mm5_date, wrf_date, wrf_time)

    implicit none
    character (len=19), intent(in)      :: mm5_date
    integer, intent(out)                :: wrf_date
    real, intent(out)                   :: wrf_time
    integer :: century_year, month, day, hour, minute, second
    integer :: day_of_year

    call split_date_char(mm5_date, century_year, month, day, hour, &
            minute, second)
    wrf_time = hms_to_wrf_time(hour,minute,second)
    day_of_year = compute_day_of_year(century_year,month,day)
    wrf_date = century_year*1000+day_of_year

  end subroutine mm5_to_wrf_date 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mm5_to_lapstime (mm5_date, lapstime)
  
    implicit none
    character(len=24), intent(in)   :: mm5_date
    integer, intent(out)            :: lapstime
    integer                         :: year,month,day,hour,minute,second
    integer                         :: jday
    integer                         :: total_years,leap_days, total_days
    integer                         :: total_hours
    
    call split_date_char(mm5_date(1:19),year,month,day,hour,minute,second)
    jday = compute_day_of_year(year,month,day)

    total_years = year - 1960
    leap_days = nint((float(total_years)+1.)/4.)
    total_days = total_years*365 + leap_days + jday - 1
    total_hours = total_days * 24 + hour
    lapstime = total_hours * 3600 + minute * 60
    return
  end subroutine mm5_to_lapstime
end module time_utils


