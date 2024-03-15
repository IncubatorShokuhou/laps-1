      subroutine jpcpack(fld,width,height,idrstmpl,cpack,lcpack)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    jpcpack
!   prgmmr: gilbert          org: w/np11    date: 2002-12-17
!
! abstract: this subroutine packs up a data field into a jpeg2000 code stream.
!   after the data field is scaled, and the reference value is subtracted out,
!   it is treated as a grayscale image and passed to a jpeg2000 encoder.
!   it also fills in grib2 data representation template 5.40 or 5.40000 with the
!   appropriate values.
!
! program history log:
! 2002-12-17  gilbert
! 2004-07-19  gilbert - added check on whether the jpeg2000 encoding was
!                       successful.  if not, try again with different encoder
!                       options.
!
! usage:    call jpcpack(fld,width,height,idrstmpl,cpack,lcpack)
!   input argument list:
!     fld()    - contains the data values to pack
!     width    - number of points in the x direction
!     height   - number of points in the y direction
!     idrstmpl - contains the array of values for data representation
!                template 5.40 or 5.40000
!                (1) = reference value - ignored on input
!                (2) = binary scale factor
!                (3) = decimal scale factor
!                (4) = number of bits for each data value - ignored on input
!                (5) = original field type - currently ignored on input
!                      data values assumed to be reals.
!                (6) = 0 - use lossless compression
!                    = 1 - use lossy compression
!                (7) = desired compression ratio, if idrstmpl(6)=1.
!                      set to 255, if idrstmpl(6)=0.
!     lcpack   - size of array cpack().
!
!   output argument list: 
!     idrstmpl - contains the array of values for data representation
!                template 5.0
!                (1) = reference value - set by jpcpack routine.
!                (2) = binary scale factor - unchanged from input
!                (3) = decimal scale factor - unchanged from input
!                (4) = number of bits containing each grayscale pixel value
!                (5) = original field type - currently set = 0 on output.
!                      data values assumed to be reals.
!                (6) = 0 - use lossless compression
!                    = 1 - use lossy compression
!                (7) = desired compression ratio, if idrstmpl(6)=1
!     cpack    - the packed data field (character*1 array)
!     lcpack   - length of packed field in cpack().
!
! remarks: none
!
! attributes:
!   language: xl fortran 90
!   machine:  ibm sp
!
!$$$

      integer,intent(in) :: width,height
      real,intent(in) :: fld(width*height)
      character(len=1),intent(out) :: cpack(*)
      integer,intent(inout) :: idrstmpl(*)
      integer,intent(inout) :: lcpack

      real(4) :: ref
      integer(4) :: iref
      integer :: ifld(width*height),retry
      integer,parameter :: zero=0
      integer :: enc_jpeg2000
      character(len=1),allocatable :: ctemp(:)
#ifdef use_jpeg2000
      
      ndpts=width*height
      bscale=2.0**real(-idrstmpl(2))
      dscale=10.0**real(idrstmpl(3))
!
!  find max and min values in the data
!
      rmax=fld(1)
      rmin=fld(1)
      do j=2,ndpts
        if (fld(j).gt.rmax) rmax=fld(j)
        if (fld(j).lt.rmin) rmin=fld(j)
      enddo
      if (idrstmpl(2).eq.0) then
         maxdif=nint(rmax*dscale)-nint(rmin*dscale)
      else
         maxdif=nint((rmax-rmin)*dscale*bscale)
      endif
!
!  if max and min values are not equal, pack up field.
!  if they are equal, we have a constant field, and the reference
!  value (rmin) is the value for each point in the field and
!  set nbits to 0.
!
      if (rmin.ne.rmax .and. maxdif.ne.0) then
        !
        !  determine which algorithm to use based on user-supplied 
        !  binary scale factor and number of bits.
        !
        if (idrstmpl(2).eq.0) then
           !
           !  no binary scaling and calculate minimum number of 
           !  bits in which the data will fit.
           !
           imin=nint(rmin*dscale)
           imax=nint(rmax*dscale)
           maxdif=imax-imin
           temp=alog(real(maxdif+1))/alog(2.0)
           nbits=ceiling(temp)
           rmin=real(imin)
           !   scale data
           do j=1,ndpts
             ifld(j)=nint(fld(j)*dscale)-imin
           enddo
        else
           !
           !  use binary scaling factor and calculate minimum number of 
           !  bits in which the data will fit.
           !
           rmin=rmin*dscale
           rmax=rmax*dscale
           maxdif=nint((rmax-rmin)*bscale)
           temp=alog(real(maxdif+1))/alog(2.0)
           nbits=ceiling(temp)
           !   scale data
           do j=1,ndpts
             ifld(j)=nint(((fld(j)*dscale)-rmin)*bscale)
           enddo
        endif
        !
        !  pack data into full octets, then do jpeg2000 encode.
        !  and calculate the length of the packed data in bytes
        !
        retry=0
        nbytes=(nbits+7)/8
        nsize=lcpack      ! needed for input to enc_jpeg2000
        allocate(ctemp(nbytes*ndpts))
        call sbytes(ctemp,ifld,0,nbytes*8,0,ndpts)
        lcpack=enc_jpeg2000(ctemp,width,height,nbits,idrstmpl(6),
     &                      idrstmpl(7),retry,cpack,nsize)
        if (lcpack.le.0) then
           print *,'jpcpack: error packing jpc=',lcpack
           if (lcpack.eq.-3) then
              retry=1
              print *,'jpcpack: retrying....'
              lcpack=enc_jpeg2000(ctemp,width,height,nbits,idrstmpl(6),
     &                         idrstmpl(7),retry,cpack,nsize)
              if (lcpack.le.0) then
                 print *,'jpcpack: retry failed.'
              else
                 print *,'jpcpack: retry successful.'
              endif
           endif
        endif
        deallocate(ctemp)

      else
        nbits=0
        lcpack=0
      endif

!
!  fill in ref value and number of bits in template 5.0
!
      call mkieee(rmin,ref,1)   ! ensure reference value is ieee format
!      call gbyte(ref,idrstmpl(1),0,32)
      iref=transfer(ref,iref)
      idrstmpl(1)=iref
      idrstmpl(4)=nbits
      idrstmpl(5)=0         ! original data were reals
      if (idrstmpl(6).eq.0) idrstmpl(7)=255       ! lossy not used

#endif /* use_jpeg2000 */
      return
      end
