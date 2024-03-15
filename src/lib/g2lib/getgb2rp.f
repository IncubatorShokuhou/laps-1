c-----------------------------------------------------------------------
      subroutine getgb2rp(lugb,cindex,extract,gribm,leng,iret)
c$$$  subprogram documentation block
c
c subprogram: getgb2rp       extracts a grib message from a file
c   prgmmr: gilbert          org: w/nmc23     date: 2003-12-31
c
c abstract: find and extracts a grib message from a file given the 
c   index for the requested field.
c   the grib message returned can contain only the requested field
c   (extract=.true.). or the complete grib message originally containing
c   the desired field can be returned (extract=.false.) even if other
c   fields were included in the grib message.
c   if the grib field is not found, then the return code will be nonzero.
c
c program history log:
c 2003-12-31  gilbert
c
c usage:    call getgb2rp(lugb,cindex,extract,gribm,leng,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib data file.
c                  file must be opened with baopen or baopenr before calling 
c                  this routine.
c     cindex       index record of the grib file  ( see docblock of
c                  subroutine ixgb2 for description of an index record.)
c     extract       logical value indicating whether to return a grib2 
c                   message with just the requested field, or the entire
c                   grib2 message containing the requested field.
c                  .true. = return grib2 message containing only the requested
c                           field.
c                  .false. = return entire grib2 message containing the
c                            requested field.
c
c   output arguments:
c     gribm         returned grib message.
c     leng         length of returned grib message in bytes.
c     iret         integer return code
c                    0      all ok
c                    97     error reading grib file
c
c subprograms called:
c   baread          byte-addressable read
c
c remarks: none 
c
c attributes:
c   language: fortran 90
c
c$$$

      integer,intent(in) :: lugb
      character(len=1),intent(in) :: cindex(*)
      logical,intent(in) :: extract
      integer,intent(out) :: leng,iret
      character(len=1),pointer,dimension(:) :: gribm
 
      integer,parameter :: zero=0
      character(len=1),allocatable,dimension(:) :: csec2,csec6,csec7
      character(len=4) :: ctemp

      iret=0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  extract grib message from file
      if ( extract ) then
         len0=16
         len8=4
         call gbyte(cindex,iskip,4*8,4*8)    ! bytes to skip in file
         call gbyte(cindex,iskp2,8*8,4*8)    ! bytes to skip for section 2
         if ( iskp2 .gt. 0 ) then
            call baread(lugb,iskip+iskp2,4,lread,ctemp)
            call gbyte(ctemp,len2,0,4*8)      ! length of section 2
            allocate(csec2(len2))
            call baread(lugb,iskip+iskp2,len2,lread,csec2)
         else
            len2=0
         endif
         call gbyte(cindex,len1,44*8,4*8)      ! length of section 1
         ipos=44+len1
         call gbyte(cindex,len3,ipos*8,4*8)      ! length of section 3
         ipos=ipos+len3
         call gbyte(cindex,len4,ipos*8,4*8)      ! length of section 4
         ipos=ipos+len4
         call gbyte(cindex,len5,ipos*8,4*8)      ! length of section 5
         ipos=ipos+len5
         call gbyte(cindex,len6,ipos*8,4*8)      ! length of section 6
         ipos=ipos+5
         call gbyte(cindex,ibmap,ipos*8,1*8)      ! bitmap indicator
         if ( ibmap .eq. 254 ) then
            call gbyte(cindex,iskp6,24*8,4*8)    ! bytes to skip for section 6
            call baread(lugb,iskip+iskp6,4,lread,ctemp)
            call gbyte(ctemp,len6,0,4*8)      ! length of section 6
         endif
         !
         !  read in section 7 from file
         !
         call gbyte(cindex,iskp7,28*8,4*8)    ! bytes to skip for section 7
         call baread(lugb,iskip+iskp7,4,lread,ctemp)
         call gbyte(ctemp,len7,0,4*8)      ! length of section 7
         allocate(csec7(len7))
         call baread(lugb,iskip+iskp7,len7,lread,csec7)

         leng=len0+len1+len2+len3+len4+len5+len6+len7+len8
         if (.not. associated(gribm)) allocate(gribm(leng))

         ! create section 0
         !
         gribm(1)='g'
         gribm(2)='r'
         gribm(3)='i'
         gribm(4)='b'
         gribm(5)=char(0)
         gribm(6)=char(0)
         gribm(7)=cindex(42)
         gribm(8)=cindex(41)
         gribm(9)=char(0)
         gribm(10)=char(0)
         gribm(11)=char(0)
         gribm(12)=char(0)
         call sbyte(gribm,leng,12*8,4*8)
         !
         ! copy section 1
         !
         gribm(17:16+len1)=cindex(45:44+len1)
         lencur=16+len1
         ipos=44+len1
         !
         ! copy section 2, if necessary
         !
         if ( iskp2 .gt. 0 ) then
           gribm(lencur+1:lencur+len2)=csec2(1:len2)
           lencur=lencur+len2
         endif
         !
         ! copy sections 3 through 5
         !
         gribm(lencur+1:lencur+len3+len4+len5)=
     &                      cindex(ipos+1:ipos+len3+len4+len5)
         lencur=lencur+len3+len4+len5
         ipos=ipos+len3+len4+len5
         !
         ! copy section 6
         !
         if ( len6 .eq. 6 .and. ibmap .ne. 254 ) then
            gribm(lencur+1:lencur+len6)=cindex(ipos+1:ipos+len6)
            lencur=lencur+len6
         else
            call gbyte(cindex,iskp6,24*8,4*8)    ! bytes to skip for section 6
            call baread(lugb,iskip+iskp6,4,lread,ctemp)
            call gbyte(ctemp,len6,0,4*8)      ! length of section 6
            allocate(csec6(len6))
            call baread(lugb,iskip+iskp6,len6,lread,csec6)
            gribm(lencur+1:lencur+len6)=csec6(1:len6)
            lencur=lencur+len6
            if ( allocated(csec6)) deallocate(csec6)
         endif
         !
         ! copy section 7
         !
         gribm(lencur+1:lencur+len7)=csec7(1:len7)
         lencur=lencur+len7
         !
         ! section 8
         !
         gribm(lencur+1)='7'
         gribm(lencur+2)='7'
         gribm(lencur+3)='7'
         gribm(lencur+4)='7'

         !  clean up
         !
         if ( allocated(csec2)) deallocate(csec2)
         if ( allocated(csec7)) deallocate(csec7)

      else    ! do not extract field from message :  get entire message

         call gbyte(cindex,iskip,4*8,4*8)    ! bytes to skip in file
         call gbyte(cindex,leng,36*8,4*8)      ! length of grib message
         if (.not. associated(gribm)) allocate(gribm(leng))
         call baread(lugb,iskip,leng,lread,gribm)
         if ( leng .ne. lread ) then
            deallocate(gribm)
            nullify(gribm)
            iret=97
            return
         endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
