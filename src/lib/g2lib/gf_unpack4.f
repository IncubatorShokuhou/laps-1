      subroutine gf_unpack4(cgrib,lcgrib,iofst,ipdsnum,ipdstmpl,
     &                      mappdslen,coordlist,numcoord,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_unpack4 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 4 (product definition section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
! 2002-01-24  gilbert  - changed to dynamically allocate arrays
!                        and to pass pointers to those arrays through
!                        the argument list.
!
! usage:    call gf_unpack4(cgrib,lcgrib,iofst,ipdsnum,ipdstmpl,mappdslen,
!    &                   coordlist,numcoord,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 4.
!
!   output argument list:      
!     iofst    - bit offset of the end of section 4, returned.
!     ipdsnum  - product definition template number ( see code table 4.0)
!     ipdstmpl - pointer to integer array containing the data values for 
!                the specified product definition
!                template ( n=ipdsnum ).  each element of this integer
!                array contains an entry (in the order specified) of product
!                defintion template 4.n
!     mappdslen- number of elements in ipdstmpl().  i.e. number of entries
!                in product defintion template 4.n  ( n=ipdsnum ).
!     coordlist- pointer to real array containing floating point values 
!                intended to document
!                the vertical discretisation associated to model data
!                on hybrid coordinate vertical levels.  (part of section 4)
!     numcoord - number of values in array coordlist.
!     ierr     - error return code.
!                0 = no error
!                5 = "grib" message contains an undefined product definition
!                    template.
!                6 = memory allocation error
!
! remarks: uses fortran 90 module pdstemplates and module re_alloc.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      use pdstemplates
      use re_alloc        !  needed for subroutine realloc

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(inout) :: iofst
      real,pointer,dimension(:) :: coordlist
      integer,pointer,dimension(:) :: ipdstmpl
      integer,intent(out) :: ipdsnum
      integer,intent(out) :: ierr,numcoord

      real(4),allocatable :: coordieee(:)
      integer,allocatable :: mappds(:)
      integer :: mappdslen
      logical needext

      ierr=0
      nullify(ipdstmpl,coordlist)

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32
      iofst=iofst+8     ! skip section number
      allocate(mappds(lensec))

      call gbyte(cgrib,numcoord,iofst,16)    ! get num of coordinate values
      iofst=iofst+16
      call gbyte(cgrib,ipdsnum,iofst,16)    ! get prod. def template num.
      iofst=iofst+16
      !   get product definition template
      call getpdstemplate(ipdsnum,mappdslen,mappds,needext,iret)
      if (iret.ne.0) then
        ierr=5
        if( allocated(mappds) ) deallocate(mappds)
        return
      endif
      !
      !   unpack each value into array ipdstmpl from the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mappds.
      !
      istat=0
      if (mappdslen.gt.0) allocate(ipdstmpl(mappdslen),stat=istat)
      if (istat.ne.0) then
         ierr=6
         nullify(ipdstmpl)
         if( allocated(mappds) ) deallocate(mappds)
         return
      endif
      do i=1,mappdslen
        nbits=iabs(mappds(i))*8
        if ( mappds(i).ge.0 ) then
          call gbyte(cgrib,ipdstmpl(i),iofst,nbits)
        else
          call gbyte(cgrib,isign,iofst,1)
          call gbyte(cgrib,ipdstmpl(i),iofst+1,nbits-1)
          if (isign.eq.1) ipdstmpl(i)=-ipdstmpl(i)
        endif
        iofst=iofst+nbits
      enddo
      !
      !   check to see if the product definition template needs to be
      !   extended.
      !   the number of values in a specific template may vary
      !   depending on data specified in the "static" part of the
      !   template.
      !
      if ( needext ) then
        call extpdstemplate(ipdsnum,ipdstmpl,newmappdslen,mappds)
        call realloc(ipdstmpl,mappdslen,newmappdslen,istat)
        !   unpack the rest of the product definition template
        do i=mappdslen+1,newmappdslen
          nbits=iabs(mappds(i))*8
          if ( mappds(i).ge.0 ) then
            call gbyte(cgrib,ipdstmpl(i),iofst,nbits)
          else
            call gbyte(cgrib,isign,iofst,1)
            call gbyte(cgrib,ipdstmpl(i),iofst+1,nbits-1)
            if (isign.eq.1) ipdstmpl(i)=-ipdstmpl(i)
          endif
          iofst=iofst+nbits
        enddo
        mappdslen=newmappdslen
      endif
      if( allocated(mappds) ) deallocate(mappds)
      !
      !   get optional list of vertical coordinate values
      !   after the product definition template, if necessary.
      !
      nullify(coordlist)
      if ( numcoord .ne. 0 ) then
         allocate (coordieee(numcoord),stat=istat1)
         allocate(coordlist(numcoord),stat=istat)
         if ((istat1+istat).ne.0) then
            ierr=6
            nullify(coordlist)
            if( allocated(coordieee) ) deallocate(coordieee)
            return
         endif
        call gbytes(cgrib,coordieee,iofst,32,0,numcoord)
        call rdieee(coordieee,coordlist,numcoord)
        deallocate (coordieee)
        iofst=iofst+(32*numcoord)
      endif
      
      return    ! end of section 4 processing
      end

