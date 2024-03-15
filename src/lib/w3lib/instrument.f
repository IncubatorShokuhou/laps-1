!-----------------------------------------------------------------------
      subroutine instrument(k,kall,ttot,tmin,tmax)
!$$$  subprogram documentation  block
!                .      .    .                                       .
! subprogram:  instrument    monitor wall-clock times, etc.
!   prgmmr: iredell          org: np23        date:1998-07-16
!
! abstract: this subprogram is useful in instrumenting a code
!   by monitoring the number of times each given section
!   of a program is invoked as well as the minimum, maximum
!   and total wall-clock time spent in the given section.
!
! program history log:
!   1998-07-16  iredell
!
! usage:    call instrument(k,kall,ttot,tmin,tmax)
!   input argument list:
!     k        - integer positive section number
!                or maximum section number in the first invocation
!                or zero to reset all wall-clock statistics
!                or negative section number to skip monitoring
!                and just return statistics.
!
!   output argument list:
!     kall     - integer number of times section is called
!     ttot     - real total seconds spent in section
!     tmin     - real minimum seconds spent in section
!     tmax     - real maximum seconds spent in section
!
! subprograms called:
!   w3utcdat     return the utc date and time
!   w3difdat     return a time interval between two dates
!
! remarks:
!   this subprogram should not be invoked from a multitasking region.
!   normally, time spent inside this subprogram is not counted.
!   wall-clock times are kept to the nearest millisecond.
!
!   example.
!     call instrument(2,kall,ttot,tmin,tmax)    ! keep stats for 2 subs
!     do k=1,n
!       call sub1
!       call instrument(1,kall,ttot,tmin,tmax)  ! accum stats for sub1
!       call sub2
!       call instrument(2,kall,ttot,tmin,tmax)  ! accum stats for sub2
!     enddo
!     print *,'sub2 stats: ',kall,ttot,tmin,tmax
!     call instrument(-1,kall,ttot,tmin,tmax)   ! return stats for sub1
!     print *,'sub1 stats: ',kall,ttot,tmin,tmax
!
! attributes:
!   language: fortran 90
!
!$$$
        implicit none
        integer,intent(in):: k
        integer,intent(out):: kall
        real,intent(out):: ttot,tmin,tmax
        integer,save:: kmax=0
        integer,dimension(:),allocatable,save:: kalls
        real,dimension(:),allocatable,save:: ttots,tmins,tmaxs
        integer,dimension(8),save:: idat
        integer,dimension(8):: jdat
        real,dimension(5):: rinc
        integer:: ka
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ka=abs(k)
!  allocate monitoring arrays if initial invocation
        if(kmax.eq.0) then
          kmax=k
          allocate(kalls(kmax))
          allocate(ttots(kmax))
          allocate(tmins(kmax))
          allocate(tmaxs(kmax))
          kalls=0
          ka=0
!  or reset all statistics back to zero
        elseif(k.eq.0) then
          kalls=0
!  or count time since last invocation against this section
        elseif(k.gt.0) then
          call w3utcdat(jdat)
          call w3difdat(jdat,idat,4,rinc)
          kalls(k)=kalls(k)+1
          if(kalls(k).eq.1) then
            ttots(k)=rinc(4)
            tmins(k)=rinc(4)
            tmaxs(k)=rinc(4)
          else
            ttots(k)=ttots(k)+rinc(4)
            tmins(k)=min(tmins(k),rinc(4))
            tmaxs(k)=max(tmaxs(k),rinc(4))
          endif
        endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  return statistics
        if(ka.ge.1.and.ka.le.kmax.and.kalls(ka).gt.0) then
          kall=kalls(ka)
          ttot=ttots(ka)
          tmin=tmins(ka)
          tmax=tmaxs(ka)
        else
          kall=0
          ttot=0
          tmin=0
          tmax=0
        endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  keep current time for next invocation
        if(k.ge.0) call w3utcdat(idat)
      end subroutine instrument
