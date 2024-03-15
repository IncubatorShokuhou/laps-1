      subroutine gribcreate(cgrib,lcgrib,listsec0,listsec1,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gribcreate 
!   prgmmr: gilbert         org: w/np11    date: 2000-04-28
!
! abstract: this subroutine initializes a new grib2 message and packs
!   grib2 sections 0 (indicator section) and 1 (identification section).
!   this routine is used with routines "addlocal", "addgrid", "addfield",
!   and "gribend" to create a complete grib2 message.  subroutine
!   gribcreate must be called first to initialize a new grib2 message.
!   also, a call to gribend is required to complete grib2 message
!   after all fields have been added.
!
! program history log:
! 2000-04-28  gilbert
!
! usage:    call gribcreate(cgrib,lcgrib,listsec0,listsec1,ierr)
!   input argument list:
!     cgrib    - character array to contain the grib2 message
!     lcgrib   - maximum length (bytes) of array cgrib.
!     listsec0 - contains information needed for grib indicator section 0.
!                must be dimensioned >= 2.
!                listsec0(1)=discipline-grib master table number
!                            (see code table 0.0)
!                listsec0(2)=grib edition number (currently 2)
!     listsec1 - contains information needed for grib identification section 1.
!                must be dimensioned >= 13.
!                listsec1(1)=id of orginating centre (common code table c-1)
!                listsec1(2)=id of orginating sub-centre (local table)
!                listsec1(3)=grib master tables version number (code table 1.0)
!                listsec1(4)=grib local tables version number (code table 1.1)
!                listsec1(5)=significance of reference time (code table 1.2)
!                listsec1(6)=reference time - year (4 digits)
!                listsec1(7)=reference time - month
!                listsec1(8)=reference time - day
!                listsec1(9)=reference time - hour
!                listsec1(10)=reference time - minute
!                listsec1(11)=reference time - second
!                listsec1(12)=production status of data (code table 1.3)
!                listsec1(13)=type of processed data (code table 1.4)
!
!   output argument list:      
!     cgrib    - character array to contain the grib2 message
!     ierr     - error return code.
!                0 = no error
!                1 = tried to use for version other than grib edition 2
!
! remarks: this routine is intended for use with routines "addlocal", 
!          "addgrid", "addfield", and "gribend" to create a complete 
!          grib2 message.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=1),intent(inout) :: cgrib(lcgrib)
      integer,intent(in) :: listsec0(*),listsec1(*)
      integer,intent(in) :: lcgrib
      integer,intent(out) :: ierr
      
      character(len=4),parameter :: grib='grib'
      integer,parameter :: zero=0,one=1
      integer,parameter :: mapsec1len=13
      integer,parameter :: 
     &        mapsec1(mapsec1len)=(/ 2,2,1,1,1,2,1,1,1,1,1,1,1 /)
      integer lensec0,iofst,ibeg

      ierr=0
!
!  currently handles only grib edition 2.
!  
      if (listsec0(2).ne.2) then
        print *,'gribcreate: can only code grib edition 2.'
        ierr=1
        return
      endif
!
!  pack section 0 - indicator section 
!  ( except for total length of grib message )
!
!      cgrib=' '
      cgrib(1)=grib(1:1)                     ! beginning of grib message
      cgrib(2)=grib(2:2)   
      cgrib(3)=grib(3:3)   
      cgrib(4)=grib(4:4)   
      call sbyte(cgrib,zero,32,16)           ! reserved for future use
      call sbyte(cgrib,listsec0(1),48,8)     ! discipline
      call sbyte(cgrib,listsec0(2),56,8)     ! grib edition number
      lensec0=16      ! bytes (octets)
!
!  pack section 1 - identification section
!
      ibeg=lensec0*8        !   calculate offset for beginning of section 1
      iofst=ibeg+32         !   leave space for length of section
      call sbyte(cgrib,one,iofst,8)     ! store section number ( 1 )
      iofst=iofst+8
      !
      !   pack up each input value in array listsec1 into the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mapsec1.
      !
      do i=1,mapsec1len
        nbits=mapsec1(i)*8
        call sbyte(cgrib,listsec1(i),iofst,nbits)
        iofst=iofst+nbits
      enddo
      !
      !   calculate length of section 1 and store it in octets
      !   1-4 of section 1.
      !
      lensec1=(iofst-ibeg)/8
      call sbyte(cgrib,lensec1,ibeg,32)
!
!  put current byte total of message into section 0
!
      call sbyte(cgrib,zero,64,32)
      call sbyte(cgrib,lensec0+lensec1,96,32)

      return
      end
