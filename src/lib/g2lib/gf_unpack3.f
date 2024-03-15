      subroutine gf_unpack3(cgrib,lcgrib,iofst,igds,igdstmpl,
     &                   mapgridlen,ideflist,idefnum,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_unpack3 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine unpacks section 3 (grid definition section)
!   starting at octet 6 of that section.  
!
! program history log:
! 2000-05-26  gilbert
! 2002-01-24  gilbert  - changed to dynamically allocate arrays
!                        and to pass pointers to those arrays through
!                        the argument list.
!
! usage:    call gf_unpack3(cgrib,lcgrib,lensec,iofst,igds,igdstmpl,
!    &                   mapgridlen,ideflist,idefnum,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     iofst    - bit offset of the beginning of section 3.
!
!   output argument list:      
!     iofst    - bit offset at the end of section 3, returned.
!     igds     - contains information read from the appropriate grib grid 
!                definition section 3 for the field being returned.
!                must be dimensioned >= 5.
!                igds(1)=source of grid definition (see code table 3.0)
!                igds(2)=number of grid points in the defined grid.
!                igds(3)=number of octets needed for each 
!                            additional grid points definition.  
!                            used to define number of
!                            points in each row ( or column ) for
!                            non-regular grids.  
!                            = 0, if using regular grid.
!                igds(4)=interpretation of list for optional points 
!                            definition.  (code table 3.11)
!                igds(5)=grid definition template number (code table 3.1)
!     igdstmpl - pointer to integer array containing the data values for 
!                the specified grid definition
!                template ( nn=igds(5) ).  each element of this integer 
!                array contains an entry (in the order specified) of grid
!                defintion template 3.nn
!     mapgridlen- number of elements in igdstmpl().  i.e. number of entries
!                in grid defintion template 3.nn  ( nn=igds(5) ).
!     ideflist - (used if igds(3) .ne. 0)  pointer to integer array containing
!                the number of grid points contained in each row ( or column ).
!                (part of section 3)
!     idefnum  - (used if igds(3) .ne. 0)  the number of entries
!                in array ideflist.  i.e. number of rows ( or columns )
!                for which optional grid points are defined.
!     ierr     - error return code.
!                0 = no error
!                5 = "grib" message contains an undefined grid definition
!                    template.
!                6 = memory allocation error
!
! remarks: uses fortran 90 module gridtemplates and module re_alloc.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      use gridtemplates
      use re_alloc        !  needed for subroutine realloc

      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib
      integer,intent(inout) :: iofst
      integer,pointer,dimension(:) :: igdstmpl,ideflist
      integer,intent(out) :: igds(5)
      integer,intent(out) :: ierr,idefnum

      integer,allocatable :: mapgrid(:)
      integer :: mapgridlen,ibyttem
      logical needext

      ierr=0
      nullify(igdstmpl,ideflist)

      call gbyte(cgrib,lensec,iofst,32)        ! get length of section
      iofst=iofst+32
      iofst=iofst+8     ! skip section number

      call gbyte(cgrib,igds(1),iofst,8)     ! get source of grid def.
      iofst=iofst+8
      call gbyte(cgrib,igds(2),iofst,32)    ! get number of grid pts.
      iofst=iofst+32
      call gbyte(cgrib,igds(3),iofst,8)     ! get num octets for opt. list
      iofst=iofst+8
      call gbyte(cgrib,igds(4),iofst,8)     ! get interpret. for opt. list
      iofst=iofst+8
      call gbyte(cgrib,igds(5),iofst,16)    ! get grid def template num.
      iofst=iofst+16
!      if (igds(1).eq.0) then
      if (igds(1).eq.0.or.igds(1).eq.255) then  ! for ecmwf test only
        allocate(mapgrid(lensec))
        !   get grid definition template
        call getgridtemplate(igds(5),mapgridlen,mapgrid,needext,
     &                       iret)
        if (iret.ne.0) then
          ierr=5
          if( allocated(mapgrid) ) deallocate(mapgrid)
          return
        endif
      else
!        igdstmpl=-1
        mapgridlen=0
        needext=.false.
      endif
      !
      !   unpack each value into array igdstmpl from the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mapgrid.
      !
      istat=0
      if (mapgridlen.gt.0) allocate(igdstmpl(mapgridlen),stat=istat)
      if (istat.ne.0) then
         ierr=6
         nullify(igdstmpl)
         if( allocated(mapgrid) ) deallocate(mapgrid)
         return
      endif
      ibyttem=0
      do i=1,mapgridlen
        nbits=iabs(mapgrid(i))*8
        if ( mapgrid(i).ge.0 ) then
          call gbyte(cgrib,igdstmpl(i),iofst,nbits)
        else
          call gbyte(cgrib,isign,iofst,1)
          call gbyte(cgrib,igdstmpl(i),iofst+1,nbits-1)
          if (isign.eq.1) igdstmpl(i)=-igdstmpl(i)
        endif
        iofst=iofst+nbits
        ibyttem=ibyttem+iabs(mapgrid(i))
      enddo
      !
      !   check to see if the grid definition template needs to be
      !   extended.
      !   the number of values in a specific template may vary
      !   depending on data specified in the "static" part of the
      !   template.
      !
      if ( needext ) then
        call extgridtemplate(igds(5),igdstmpl,newmapgridlen,mapgrid)
        !   unpack the rest of the grid definition template
        call realloc(igdstmpl,mapgridlen,newmapgridlen,istat)
        do i=mapgridlen+1,newmapgridlen
          nbits=iabs(mapgrid(i))*8
          if ( mapgrid(i).ge.0 ) then
            call gbyte(cgrib,igdstmpl(i),iofst,nbits)
          else
            call gbyte(cgrib,isign,iofst,1)
            call gbyte(cgrib,igdstmpl(i),iofst+1,nbits-1)
            if (isign.eq.1) igdstmpl(i)=-igdstmpl(i)
          endif
          iofst=iofst+nbits
          ibyttem=ibyttem+iabs(mapgrid(i))
        enddo
        mapgridlen=newmapgridlen
      endif
      if( allocated(mapgrid) ) deallocate(mapgrid)
      !
      !   unpack optional list of numbers defining number of points
      !   in each row or column, if included.  this is used for non regular
      !   grids.
      !
      if ( igds(3).ne.0 ) then
         nbits=igds(3)*8
         idefnum=(lensec-14-ibyttem)/igds(3)
         istat=0
         if (idefnum.gt.0) allocate(ideflist(idefnum),stat=istat)
         if (istat.ne.0) then
            ierr=6
            nullify(ideflist)
            return
         endif
         call gbytes(cgrib,ideflist,iofst,nbits,0,idefnum)
         iofst=iofst+(nbits*idefnum)
      else
         idefnum=0
         nullify(ideflist)
      endif
      
      return    ! end of section 3 processing
      end
