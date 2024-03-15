      subroutine getdim(csec3,lcsec3,width,height,iscan)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    getdim 
!   prgmmr: gilbert         org: w/np11    date: 2002-12-11
!
! abstract: this subroutine returns the dimensions and scanning mode of 
!   a grid definition packed in grib2 grid definition section 3 format.
!
! program history log:
! 2002-12-11  gilbert
!
! usage:    call getdim(csec3,lcsec3,width,height,iscan)
!   input argument list:
!     csec3    - character array that contains the packed grib2 gds
!    lcsec3    - length (in octets) of section 3
!
!   output argument list:      
!     width    - x (or i) dimension of the grid.
!     height   - y (or j) dimension of the grid.
!     iscan    - scanning mode ( see code table 3.4 )
!
! remarks:  returns width and height set to zero, if grid template
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
      integer,intent(out) :: width,height,iscan
      
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
           case (0:3)   ! lat/lon
              width=igdstmpl(8)
              height=igdstmpl(9)
              iscan=igdstmpl(19)
           case (10)   ! mercator
              width=igdstmpl(8)
              height=igdstmpl(9)
              iscan=igdstmpl(16)
           case (20)   ! polar stereographic
              width=igdstmpl(8)
              height=igdstmpl(9)
              iscan=igdstmpl(18)
           case (30)   ! lambert conformal
              width=igdstmpl(8)
              height=igdstmpl(9)
              iscan=igdstmpl(18)
           case (40:43)   ! gaussian
              width=igdstmpl(8)
              height=igdstmpl(9)
              iscan=igdstmpl(19)
           case (90)   ! space view/orthographic
              width=igdstmpl(8)
              height=igdstmpl(9)
              iscan=igdstmpl(17)
           case (110)   ! equatorial azimuthal
              width=igdstmpl(8)
              height=igdstmpl(9)
              iscan=igdstmpl(16)
           case default
              width=0
              height=0
              iscan=0
         end select
      else
         width=0
         height=0
      endif
        !
      if (associated(igdstmpl)) deallocate(igdstmpl)
      if (associated(list_opt)) deallocate(list_opt)

      return
      end
