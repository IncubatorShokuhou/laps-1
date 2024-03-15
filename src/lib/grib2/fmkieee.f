      real function fmkieee(a)
 !$$$ subprogram documentation block
 ! . . . .
 ! subprogram: mkieee 
 ! prgmmr: gilbert org: w/np11 date: 2000-05-09
 !
 ! abstract: this subroutine stores a list of real values in 
 ! 32-bit ieee floating point format.
 !
 ! program history log:
 ! 2000-05-09 gilbert
 !
 ! usage: call mkieee(a)
 ! input argument list:
 ! a - input floating point value.
 !
 ! output argument list: none.
 !
 ! remarks: none
 !
 ! attributes:
 ! language: fortran 90
 ! machine: ibm sp
 !
 !$$$

      parameter (two23=2.**23)
      parameter (two126=2.**126)
c
      equivalence(rtemp,ieee)
c
      alog2=alog(2.0)

      ieee=0

      if (a.ne.0.) then
 !
 ! set sign bit (bit 31 - leftmost bit)
 !
         if (a.lt.0.0) then
           ieee=ibset(ieee,31)
           atemp=abs(a)
         else
           ieee=ibclr(ieee,31)
           atemp=a
         endif
 !
 ! determine exponent n with base 2
 !
         n=int(flr(alog(atemp)/alog2))
         iexp=n+127
         if (n.gt.127) iexp=255 ! overflow
         if (n.lt.-127) iexp=0
         ! set exponent bits ( bits 30-23 )
         call mvbits(iexp,0,8,ieee,23)
 !
 ! determine mantissa
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
         ! set mantissa bits ( bits 22-0 )
         call mvbits(imant,0,23,ieee,0)
      endif
c
      fmkieee=rtemp
      return
      end
