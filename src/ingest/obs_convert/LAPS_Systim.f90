subroutine read_systim(ddmmyyhm,timeofch,yymmddhm,exstatus)

!==============================================================================
!  this routine decodes an ascii time string into integers of day, month, year,
!  hour, minute.
!
!  history:
!	creation:	yuanfu xie	jun 2007
!==============================================================================

  implicit none

  character, intent(in)  :: ddmmyyhm*16
  character, intent(out) :: timeofch*14
  integer,   intent(out) :: yymmddhm(5)
  integer,   intent(out) :: exstatus

  ! local variable:
  integer :: zero

  zero = ichar('0')
  exstatus = 0

  ! second:
  timeofch(13:14) = '00'

  ! day:
  yymmddhm(3) = (ichar(ddmmyyhm(1:1))-zero)*10+ichar(ddmmyyhm(2:2))-zero
  timeofch(7:8) = ddmmyyhm(1:2)

  ! month:
  select case (ddmmyyhm(4:6))
  case ('jan') 
    yymmddhm(2) =  1
    timeofch(5:6) = '01'
  case ('feb') 
    yymmddhm(2) =  2
    timeofch(5:6) = '02'
  case ('mar')
    yymmddhm(2) =  3
    timeofch(5:6) = '03'
  case ('apr')
    yymmddhm(2) =  4
    timeofch(5:6) = '04'
  case ('may')
    yymmddhm(2) =  5
    timeofch(5:6) = '05'
  case ('jun')
    yymmddhm(2) =  6
    timeofch(5:6) = '06'
  case ('jul')
    yymmddhm(2) =  7
    timeofch(5:6) = '07'
  case ('aug')
    yymmddhm(2) =  8
    timeofch(5:6) = '08'
  case ('sep')
    yymmddhm(2) =  9
    timeofch(5:6) = '09'
  case ('oct')
    yymmddhm(2) = 10
    timeofch(5:6) = '10'
  case ('nov')
    yymmddhm(2) = 11
    timeofch(5:6) = '11'
  case ('dec')
    yymmddhm(2) = 12
    timeofch(5:6) = '12'
  end select

  ! year:
  yymmddhm(1) = (ichar(ddmmyyhm(8 : 8))-zero)*1000+ &
                (ichar(ddmmyyhm(9 : 9))-zero)*100 + &
                (ichar(ddmmyyhm(10:10))-zero)*10+ &
                 ichar(ddmmyyhm(11:11))-zero
  timeofch(1:4) = ddmmyyhm(8:11)

  yymmddhm(4) = (ichar(ddmmyyhm(13:13))-zero)*10+&
                 ichar(ddmmyyhm(14:14))-zero
  timeofch(9:10) = ddmmyyhm(13:14)

  yymmddhm(5) = (ichar(ddmmyyhm(15:15))-zero)*10+&
                 ichar(ddmmyyhm(16:16))-zero
  timeofch(11:12) = ddmmyyhm(15:16)

end subroutine read_systim
