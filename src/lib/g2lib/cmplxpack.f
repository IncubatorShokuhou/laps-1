      subroutine cmplxpack(fld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    cmplxpack
!   prgmmr: gilbert          org: w/np11    date: 2004-08-27
!
! abstract: this subroutine packs up a data field using a complex
!   packing algorithm as defined in the grib2 documention.  it
!   supports grib2 complex packing templates with or without
!   spatial differences (i.e. drts 5.2 and 5.3).
!   it also fills in grib2 data representation template 5.2 or 5.3 
!   with the appropriate values.
!
! program history log:
! 2004-08-27  gilbert
!
! usage:    call cmplxpack(fld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
!   input argument list:
!     fld()    - contains the data values to pack
!     ndpts    - the number of data values in array fld()
!     idrsnum  - data representation template number 5.n
!                must equal 2 or 3.
!     idrstmpl - contains the array of values for data representation
!                template 5.2 or 5.3
!                (1) = reference value - ignored on input
!                (2) = binary scale factor
!                (3) = decimal scale factor
!                    .
!                    .
!                (7) = missing value management
!                (8) = primary missing value
!                (9) = secondary missing value
!                    .
!                    .
!               (17) = order of spatial differencing  ( 1 or 2 )
!                    .
!                    .
!
!   output argument list: 
!     idrstmpl - contains the array of values for data representation
!                template 5.3
!                (1) = reference value - set by compack routine.
!                (2) = binary scale factor - unchanged from input
!                (3) = decimal scale factor - unchanged from input
!                    .
!                    .
!     cpack    - the packed data field (character*1 array)
!     lcpack   - length of packed field cpack().
!
! remarks: none
!
! attributes:
!   language: xl fortran 90
!   machine:  ibm sp
!
!$$$

      integer,intent(in) :: ndpts,idrsnum
      real,intent(in) :: fld(ndpts)
      character(len=1),intent(out) :: cpack(*)
      integer,intent(inout) :: idrstmpl(*)
      integer,intent(out) :: lcpack

      

      if ( idrstmpl(7) .eq. 0 ) then       ! no internal missing values
         call compack(fld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
      elseif ( idrstmpl(7).eq.1 .or. idrstmpl(7).eq.2) then
         call misspack(fld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
      else
         print *,'cmplxpack: don:t recognize missing value option.'
         lcpack=-1
      endif

      return
      end
