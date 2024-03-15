      subroutine gf_unpack5(cgrib,lcgrib,iofst,ndpts,idrsnum,idrstmpl,
     &                   mapdrslen,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_unpack5 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 5 (data representation section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
! 2002-01-24  gilbert  - changed to dynamically allocate arrays
!                        and to pass pointers to those arrays through
!                        the argument list.
!
! usage:    call gf_unpack5(cgrib,lcgrib,iofst,ndpts,idrsnum,idrstmpl,
!                        mapdrslen,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 5.
!
!   output argument list:      
!     iofst    - bit offset at the end of section 5, returned.
!     ndpts    - number of data points unpacked and returned.
!     idrsnum  - data representation template number ( see code table 5.0)
!     idrstmpl - pointer to an integer array containing the data values for 
!                the specified data representation
!                template ( n=idrsnum ).  each element of this integer
!                array contains an entry (in the order specified) of data
!                representation template 5.n
!     mapdrslen- number of elements in idrstmpl().  i.e. number of entries
!                in data representation template 5.n  ( n=idrsnum ).
!     ierr     - error return code.
!                0 = no error
!                6 = memory allocation error
!                7 = "grib" message contains an undefined data
!                    representation template.
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      use drstemplates
      use re_alloc        !  needed for subroutine realloc

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(inout) :: iofst
      integer,intent(out) :: ndpts,idrsnum
      integer,pointer,dimension(:) :: idrstmpl
      integer,intent(out) :: ierr

      integer,allocatable :: mapdrs(:)
      integer :: mapdrslen
      logical needext

      ierr=0
      nullify(idrstmpl)

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32
      iofst=iofst+8     ! skip section number
      allocate(mapdrs(lensec))

      call gbyte(cgrib,ndpts,iofst,32)    ! get num of data points
      iofst=iofst+32
      call gbyte(cgrib,idrsnum,iofst,16)     ! get data rep template num.
      iofst=iofst+16
      !   gen data representation template
      call getdrstemplate(idrsnum,mapdrslen,mapdrs,needext,iret)
      if (iret.ne.0) then
        ierr=7
        if( allocated(mapdrs) ) deallocate(mapdrs)
        return
      endif
      !
      !   unpack each value into array ipdstmpl from the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mappds.
      !
      istat=0
      if (mapdrslen.gt.0) allocate(idrstmpl(mapdrslen),stat=istat)
      if (istat.ne.0) then
         ierr=6
         nullify(idrstmpl)
         if( allocated(mapdrs) ) deallocate(mapdrs)
         return
      endif
      do i=1,mapdrslen
        nbits=iabs(mapdrs(i))*8
        if ( mapdrs(i).ge.0 ) then
          call gbyte(cgrib,idrstmpl(i),iofst,nbits)
        else
          call gbyte(cgrib,isign,iofst,1)
          call gbyte(cgrib,idrstmpl(i),iofst+1,nbits-1)
          if (isign.eq.1) idrstmpl(i)=-idrstmpl(i)
        endif
        iofst=iofst+nbits
      enddo
      !
      !   check to see if the data representation template needs to be
      !   extended.
      !   the number of values in a specific template may vary
      !   depending on data specified in the "static" part of the
      !   template.
      !
      if ( needext ) then
        call extdrstemplate(idrsnum,idrstmpl,newmapdrslen,mapdrs)
        call realloc(idrstmpl,mapdrslen,newmapdrslen,istat)
        !   unpack the rest of the data representation template
        do i=mapdrslen+1,newmapdrslen
          nbits=iabs(mapdrs(i))*8
          if ( mapdrs(i).ge.0 ) then
            call gbyte(cgrib,idrstmpl(i),iofst,nbits)
          else
            call gbyte(cgrib,isign,iofst,1)
            call gbyte(cgrib,idrstmpl(i),iofst+1,nbits-1)
            if (isign.eq.1) idrstmpl(i)=-idrstmpl(i)
          endif
          iofst=iofst+nbits
        enddo
        mapdrslen=newmapdrslen
      endif
      if( allocated(mapdrs) ) deallocate(mapdrs)

      return    ! end of section 5 processing
      end

