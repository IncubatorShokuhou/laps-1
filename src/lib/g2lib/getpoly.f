      subroutine getpoly(csec3,lcsec3,jj,kk,mm)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getpoly 
!   prgmmr: gilbert         org: w/np11    date: 2002-12-11
!
! abstract: this subroutine returns the j, k, and m pentagonal resolution
!   parameters specified in a grib grid definition section used
!   spherical harmonic coefficients using gdt 5.50 through 5.53
!
! program history log:
! 2002-12-11  gilbert
!
! usage:    call getpoly(csec3,lcsec3,jj,kk,mm)
!   input argument list:
!     csec3    - character array that contains the packed grib2 gds
!    lcsec3    - length (in octets) of section 3
!
!   output argument list:      
!         jj   = j - pentagonal resolution parameter
!         kk   = k - pentagonal resolution parameter
!         mm   = m - pentagonal resolution parameter
!
! remarks:  returns jj, kk, and mm set to zero, if grid template
!           not recognized.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
!      use grib_mod
    
      character(len=1),intent(in) :: csec3(*)
      integer,intent(in) :: lcsec3
      integer,intent(out) :: jj,kk,mm
      
      integer,pointer,dimension(:) :: igdstmpl,list_opt
      integer :: igds(5)
      integer iofst,igdtlen,num_opt,jerr

      interface
         subroutine gf_unpack3(cgrib,lcgrib,iofst,igds,igdstmpl,
     &                         mapgridlen,ideflist,idefnum,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib
            integer,intent(inout) :: iofst
            integer,pointer,dimension(:) :: igdstmpl,ideflist
            integer,intent(out) :: igds(5)
            integer,intent(out) :: ierr,idefnum
         end subroutine gf_unpack3
      end interface

      nullify(igdstmpl,list_opt)
        !
      iofst=0       ! set offset to beginning of section
      call gf_unpack3(csec3,lcsec3,iofst,igds,igdstmpl,
     &                 igdtlen,list_opt,num_opt,jerr)
      if (jerr.eq.0) then
         selectcase( igds(5) )     !  template number
           case (50:53)   ! spherical harmonic coefficients
              jj=igdstmpl(1)
              kk=igdstmpl(2)
              mm=igdstmpl(3)
           case default
              jj=0
              kk=0
              mm=0
         end select
      else
         jj=0
         kk=0
         mm=0
      endif
        !
      if (associated(igdstmpl)) deallocate(igdstmpl)
      if (associated(list_opt)) deallocate(list_opt)

      return
      end
