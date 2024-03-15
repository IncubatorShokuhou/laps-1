      subroutine mkieee(a,rieee,num)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    mkieee 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-09
!
! abstract: this subroutine stores a list of real values in 
!   32-bit ieee floating point format.
!
! program history log:
! 2000-05-09  gilbert
!
! usage:    call mkieee(a,rieee,num)
!   input argument list:
!     a        - input array of floating point values.
!     num      - number of floating point values to convert.
!
!   output argument list:      
!     rieee    - output array of floating point values in 32-bit ieee format.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      real,intent(in) :: a(num)
      real(4),intent(out) :: rieee(num)
      integer,intent(in) :: num

      integer(4) :: ieee 

      real,save :: two23
      real,save :: two126
      integer,save :: once=0

      if ( once .eq. 0 ) then
         once=1
         two23=scale(1.0,23)
         two126=scale(1.0,126)
      endif

      alog2=alog(2.0)

      do j=1,num
        ieee=0

        if (a(j).eq.0.) then
          ieee=0
          rieee(j)=transfer(ieee,rieee(j))
!       write(6,fmt='(f20.10,5x,b32)') a,a
!       write(6,fmt='(f20.10,5x,b32)') rieee,rieee
          cycle
        endif
        
!
!  set sign bit (bit 31 - leftmost bit)
!
        if (a(j).lt.0.0) then
          ieee=ibset(ieee,31)
          atemp=abs(a(j))
        else
          ieee=ibclr(ieee,31)
          atemp=a(j)
        endif
!
!  determine exponent n with base 2
!
        n=floor(alog(atemp)/alog2)
        iexp=n+127
        if (n.gt.127) iexp=255     ! overflow
        if (n.lt.-127) iexp=0
        !      set exponent bits ( bits 30-23 )
        call mvbits(iexp,0,8,ieee,23)
!
!  determine mantissa
! 
        if (iexp.ne.255) then
          if (iexp.ne.0) then
            atemp=(atemp/(2.0**n))-1.0
          else
            atemp=atemp*two126
          endif
          imant=nint(atemp*two23)
        else
          imant=0
        endif
        !      set mantissa bits ( bits 22-0 )
        call mvbits(imant,0,23,ieee,0)
!
!  transfer ieee bit string to real variable
!
        rieee(j)=transfer(ieee,rieee(j))
!       write(6,fmt='(f20.10,5x,b32)') a,a
!       write(6,fmt='(f20.10,5x,b32)') rieee,rieee

      enddo

      return
      end

