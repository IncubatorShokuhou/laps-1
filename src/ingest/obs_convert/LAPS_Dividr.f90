subroutine laps_divider

!==============================================================================
!  this routine prints a line divider to format an output.
!
!  history:
!	creation:	yuanfu xie	jun 2007
!==============================================================================

  character :: symbol
  integer   :: loopvr,ncount

  symbol = '='
  ncount = 70

  write(6,*) (symbol,loopvr=1,ncount)

end subroutine laps_divider
