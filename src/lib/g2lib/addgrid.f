      subroutine addgrid(cgrib,lcgrib,igds,igdstmpl,igdstmplen,
     &                   ideflist,idefnum,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    addgrid 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-01
!
! abstract: this subroutine packs up a grid definition section (section 3) 
!   and adds it to a grib2 message.
!   this routine is used with routines "gribcreate", "addlocal", "addfield",
!   and "gribend" to create a complete grib2 message.  subroutine
!   gribcreate must be called first to initialize a new grib2 message.
!
! program history log:
! 2000-05-01  gilbert
!
! usage:    call addgrid(cgrib,lcgrib,igds,igdstmpl,igdstmplen,
!                        ideflist,idefnum,ierr)
!   input argument list:
!     cgrib    - character array to contain the grib2 message
!     lcgrib   - maximum length (bytes) of array cgrib.
!     igds     - contains information needed for grib grid definition section 3.
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
!     igdstmpl - contains the data values for the specified grid definition
!                template ( nn=igds(5) ).  each element of this integer 
!                array contains an entry (in the order specified) of grid
!                defintion template 3.nn
!   igdstmplen - max dimension of igdstmpl()
!     ideflist - (used if igds(3) .ne. 0)  this array contains the
!                number of grid points contained in each row ( or column )
!      idefnum - (used if igds(3) .ne. 0)  the number of entries
!                in array ideflist.  i.e. number of rows ( or columns )
!                for which optional grid points are defined.
!
!   output argument list:      
!     cgrib    - character array to contain the grib2 message
!     ierr     - error return code.
!                0 = no error
!                1 = grib message was not initialized.  need to call
!                    routine gribcreate first.
!                2 = grib message already complete.  cannot add new section.
!                3 = sum of section byte counts doesn't add to total byte count.
!                4 = previous section was not 1, 2 or 7.
!                5 = could not find requested grid definition template.
!
! remarks: note that the local use section ( section 2 ) can only follow
!          section 1 or section 7 in a grib2 message.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      use gridtemplates

      character(len=1),intent(inout) :: cgrib(lcgrib)
      integer,intent(in) :: igds(*),igdstmpl(*),ideflist(idefnum)
      integer,intent(in) :: lcgrib,idefnum,igdstmplen
      integer,intent(out) :: ierr
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4):: ctemp
      integer:: mapgrid(igdstmplen)
      integer,parameter :: one=1,three=3
      integer lensec3,iofst,ibeg,lencurr,len,mapgridlen
      logical needext
 
      ierr=0
!
!  check to see if beginning of grib message exists
!
      ctemp=cgrib(1)//cgrib(2)//cgrib(3)//cgrib(4)
      if ( ctemp.ne.grib ) then
        print *,'addgrid: grib not found in given message.'
        print *,'addgrid: call to routine gribcreate required',
     &          ' to initialize grib messge.'
        ierr=1
        return
      endif
!
!  get current length of grib message
!  
      call gbyte(cgrib,lencurr,96,32)
!
!  check to see if grib message is already complete
!  
      ctemp=cgrib(lencurr-3)//cgrib(lencurr-2)//cgrib(lencurr-1)
     &      //cgrib(lencurr)
      if ( ctemp.eq.c7777 ) then
        print *,'addgrid: grib message already complete.  cannot',
     &          ' add new section.'
        ierr=2
        return
      endif
!
!  loop through all current sections of the grib message to
!  find the last section number.
!
      len=16    ! length of section 0
      do 
      !    get section number and length of next section
        iofst=len*8
        call gbyte(cgrib,ilen,iofst,32)
        iofst=iofst+32
        call gbyte(cgrib,isecnum,iofst,8)
        len=len+ilen
      !    exit loop if last section reached
        if ( len.eq.lencurr ) exit
      !    if byte count for each section doesn't match current
      !    total length, then there is a problem.
        if ( len.gt.lencurr ) then
          print *,'addgrid: section byte counts don''t add to total.'
          print *,'addgrid: sum of section byte counts = ',len
          print *,'addgrid: total byte count in section 0 = ',lencurr
          ierr=3
          return
        endif
      enddo
!
!  section 3 can only be added after sections 1, 2 and 7.
!
      if ( (isecnum.ne.1) .and. (isecnum.ne.2) .and. 
     &     (isecnum.ne.7) ) then
        print *,'addgrid: section 3 can only be added after section',
     &          ' 1, 2 or 7.'
        print *,'addgrid: section ',isecnum,' was the last found in',
     &          ' given grib message.'
        ierr=4
        return
      endif
!
!  add section 3  - grid definition section
!
      ibeg=lencurr*8        !   calculate offset for beginning of section 3
      iofst=ibeg+32         !   leave space for length of section
      call sbyte(cgrib,three,iofst,8)     ! store section number ( 3 )
      iofst=iofst+8
      call sbyte(cgrib,igds(1),iofst,8)     ! store source of grid def.
      iofst=iofst+8
      call sbyte(cgrib,igds(2),iofst,32)    ! store number of data pts.
      iofst=iofst+32
      call sbyte(cgrib,igds(3),iofst,8)     ! store number of extra octets.
      iofst=iofst+8
      call sbyte(cgrib,igds(4),iofst,8)     ! store interp. of extra octets.
      iofst=iofst+8
      !   if octet 6 is not equal to zero, grid definition template may
      !   not be supplied.
      if ( igds(1).eq.0 ) then
        call sbyte(cgrib,igds(5),iofst,16)  ! store grid def template num.
      else
        call sbyte(cgrib,65535,iofst,16)   ! store missing value as grid def template num.
      endif
      iofst=iofst+16
      !
      !   get grid definition template
      !
      if (igds(1).eq.0) then
        call getgridtemplate(igds(5),mapgridlen,mapgrid,needext,
     &                       iret)
        if (iret.ne.0) then
          ierr=5
          return
        endif
        !
        !   extend the grid definition template, if necessary.
        !   the number of values in a specific template may vary
        !   depending on data specified in the "static" part of the
        !   template.
        !
        if ( needext ) then
          call extgridtemplate(igds(5),igdstmpl,mapgridlen,mapgrid)
        endif
      else
        mapgridlen=0
      endif
      !
      !   pack up each input value in array igdstmpl into the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mapgrid.
      !
      do i=1,mapgridlen
        nbits=iabs(mapgrid(i))*8
        if ( (mapgrid(i).ge.0).or.(igdstmpl(i).ge.0) ) then
          call sbyte(cgrib,igdstmpl(i),iofst,nbits)
        else
          call sbyte(cgrib,one,iofst,1)
          call sbyte(cgrib,iabs(igdstmpl(i)),iofst+1,nbits-1)
        endif
        iofst=iofst+nbits
      enddo
      !
      !   if requested,
      !   insert optional list of numbers defining number of points
      !   in each row or column.  this is used for non regular
      !   grids.
      !
      if ( igds(3).ne.0 ) then
         nbits=igds(3)*8
         call sbytes(cgrib,ideflist,iofst,nbits,0,idefnum)
         iofst=iofst+(nbits*idefnum)
      endif
      !
      !   calculate length of section 3 and store it in octets
      !   1-4 of section 3.
      !
      lensec3=(iofst-ibeg)/8
      call sbyte(cgrib,lensec3,ibeg,32)

!
!  update current byte total of message in section 0
!
      call sbyte(cgrib,lencurr+lensec3,96,32)

      return
      end

