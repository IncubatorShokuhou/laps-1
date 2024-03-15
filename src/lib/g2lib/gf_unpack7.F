      subroutine gf_unpack7(cgrib,lcgrib,iofst,igdsnum,igdstmpl,
     &                      idrsnum,idrstmpl,ndpts,fld,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_unpack7 
!   prgmmr: gilbert         org: w/np11    date: 2002-01-24
!
! abstract: this subroutine unpacks grib2 section 7 (data section).
!
! program history log:
! 2002-01-24  gilbert
! 2002-12-17  gilbert  - added support for new templates using 
!                        png and jpeg2000 algorithms/templates.
! 2004-12-29  gilbert  - added check on comunpack return code.
!
! usage:    call gf_unpack7(cgrib,lcgrib,iofst,igdsnum,igdstmpl,
!    &                      idrsnum,idrstmpl,ndpts,fld,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 7.
!     igdsnum  - grid definition template number ( see code table 3.0)
!                (only required to unpack drt 5.51)
!     igdstmpl - pointer to an integer array containing the data values for
!                the specified grid definition
!                template ( n=igdsnum ).  each element of this integer
!                array contains an entry (in the order specified) of grid
!                definition template 3.n
!                (only required to unpack drt 5.51)
!     idrsnum  - data representation template number ( see code table 5.0)
!     idrstmpl - pointer to an integer array containing the data values for
!                the specified data representation
!                template ( n=idrsnum ).  each element of this integer
!                array contains an entry (in the order specified) of data
!                representation template 5.n
!     ndpts    - number of data points unpacked and returned.
!
!   output argument list:      
!     iofst    - bit offset at the end of section 7, returned.
!     fld()    - pointer to a real array containing the unpacked data field.
!     ierr     - error return code.
!                0 = no error
!                4 = unrecognized data representation template
!                5 = one of gdt 3.50 through 3.53 required to unpack drt 5.51
!                6 = memory allocation error
!                7 = corrupt section 7.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib,ndpts,igdsnum,idrsnum
      integer,intent(inout) :: iofst
      integer,pointer,dimension(:) :: igdstmpl,idrstmpl
      integer,intent(out) :: ierr
      real,pointer,dimension(:) :: fld


      ierr=0
      nullify(fld)

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32    
      iofst=iofst+8     ! skip section number

      ipos=(iofst/8)+1
      istat=0
      allocate(fld(ndpts),stat=istat)
      if (istat.ne.0) then
         ierr=6
         return
      endif

      if (idrsnum.eq.0) then
        call simunpack(cgrib(ipos),lensec-5,idrstmpl,ndpts,fld)
      elseif (idrsnum.eq.2.or.idrsnum.eq.3) then
        call comunpack(cgrib(ipos),lensec-5,lensec,idrsnum,idrstmpl,
     &                 ndpts,fld,ier)
        if ( ier .ne. 0 ) then
           ierr=7
           return
        endif
      elseif (idrsnum.eq.50) then      !  spectral simple
        call simunpack(cgrib(ipos),lensec-5,idrstmpl,ndpts-1,
     &                 fld(2))
        ieee=idrstmpl(5)
        call rdieee(ieee,fld(1),1)
      elseif (idrsnum.eq.51) then      !  spectral complex
        if (igdsnum.ge.50.and.igdsnum.le.53) then
          call specunpack(cgrib(ipos),lensec-5,idrstmpl,ndpts,
     &                    igdstmpl(1),igdstmpl(2),igdstmpl(3),fld)
        else
          print *,'gf_unpack7: cannot use gdt 3.',igdsnum,
     &            ' to unpack data section 5.51.'
          ierr=5
          nullify(fld)
          return
        endif
#ifdef use_jpeg2000
      elseif (idrsnum.eq.40 .or. idrsnum.eq.40000) then
        call jpcunpack(cgrib(ipos),lensec-5,idrstmpl,ndpts,fld)
#endif  /* use_jpeg2000 */
#ifdef use_png
      elseif (idrsnum.eq.41 .or. idrsnum.eq.40010) then
        call pngunpack(cgrib(ipos),lensec-5,idrstmpl,ndpts,fld)
#endif  /* use_png */
      else
        print *,'gf_unpack7: data representation template ',idrsnum,
     &          ' not yet implemented.'
        ierr=4
        nullify(fld)
        return
      endif

      iofst=iofst+(8*lensec)
      
      return    ! end of section 7 processing
      end

