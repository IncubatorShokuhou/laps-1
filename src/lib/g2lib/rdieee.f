      subroutine rdieee(rieee,a,num)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    rdieee 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-09
!
! abstract: this subroutine reads a list of real values in 
!   32-bit ieee floating point format.
!
! program history log:
! 2000-05-09  gilbert
!
! usage:    call rdieee(rieee,a,num)
!   input argument list:
!     rieee    - input array of floating point values in 32-bit ieee format.
!     num      - number of floating point values to convert.
!
!   output argument list:      
!     a        - output array of real values.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      real(4),intent(in) :: rieee(num)
      real,intent(out) :: a(num)
      integer,intent(in) :: num

      integer(4) :: ieee

      real,save :: two23
      real,save :: two126
      integer,save :: once=0

      if ( once .eq. 0 ) then
         once=1
         two23=scale(1.0,-23)
         two126=scale(1.0,-126)
      endif

      do j=1,num
!
!  transfer ieee bit string to integer variable
!
        ieee=transfer(rieee(j),ieee)
!
!  extract sign bit, exponent, and mantissa
!
        isign=ibits(ieee,31,1)
        iexp=ibits(ieee,23,8)
        imant=ibits(ieee,0,23)
        sign=1.0
        if (isign.eq.1) sign=-1.0
        
        if ( (iexp.gt.0).and.(iexp.lt.255) ) then
          temp=2.0**(iexp-127)
          a(j)=sign*temp*(1.0+(two23*real(imant)))

        elseif ( iexp.eq.0 ) then
          if ( imant.ne.0 ) then
            a(j)=sign*two126*two23*real(imant)
          else
            a(j)=sign*0.0
          endif

        elseif ( iexp.eq.255 ) then
          a(j)=sign*huge(a(j))

        endif

      enddo

      return
      end

