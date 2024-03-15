      subroutine addfield(cgrib,lcgrib,ipdsnum,ipdstmpl,ipdstmplen,
     &                    coordlist,numcoord,idrsnum,idrstmpl,
     &                    idrstmplen,fld,ngrdpts,ibmap,bmap,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    addfield 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-02
!
! abstract: this subroutine packs up sections 4 through 7 for a given field
!   and adds them to a grib2 message.  they are product definition section,
!   data representation section, bit-map section and data section, 
!   respectively.
!   this routine is used with routines "gribcreate", "addlocal", "addgrid",
!   and "gribend" to create a complete grib2 message.  subroutine
!   gribcreate must be called first to initialize a new grib2 message.
!   also, subroutine addgrid must be called after gribcreate and
!   before this routine to add the appropriate grid description to
!   the grib2 message.   also, a call to gribend is required to complete 
!   grib2 message after all fields have been added.
!
! program history log:
! 2000-05-02  gilbert
! 2002-12-17  gilbert  - added support for new templates using
!                        png and jpeg2000 algorithms/templates.
! 2004-06-22  gilbert  - added check to determine if packing algorithm failed.
!
! usage:    call addfield(cgrib,lcgrib,ipdsnum,ipdstmpl,ipdstmplen,
!                         coordlist,numcoord,idrsnum,idrstmpl,
!                         idrstmplen,fld,ngrdpts,ibmap,bmap,ierr)
!   input argument list:
!     cgrib    - character array to contain the grib2 message
!     lcgrib   - maximum length (bytes) of array cgrib.
!     ipdsnum  - product definition template number ( see code table 4.0)
!     ipdstmpl - contains the data values for the specified product definition
!                template ( n=ipdsnum ).  each element of this integer 
!                array contains an entry (in the order specified) of product
!                defintion template 4.n
!   ipdstmplen - max dimension of ipdstmpl()
!     coordlist- array containg floating point values intended to document
!                the vertical discretisation associated to model data
!                on hybrid coordinate vertical levels.
!     numcoord - number of values in array coordlist.
!     idrsnum  - data representation template number ( see code table 5.0 )
!     idrstmpl - contains the data values for the specified data representation
!                template ( n=idrsnum ).  each element of this integer 
!                array contains an entry (in the order specified) of data
!                representation template 5.n
!                note that some values in this template (eg. reference
!                values, number of bits, etc...) may be changed by the
!                data packing algorithms.
!                use this to specify scaling factors and order of
!                spatial differencing, if desired.
!   idrstmplen - max dimension of idrstmpl()
!     fld()    - array of data points to pack.
!     ngrdpts  - number of data points in grid.
!                i.e.  size of fld and bmap.
!     ibmap    - bitmap indicator ( see code table 6.0 )
!                0 = bitmap applies and is included in section 6.
!                1-253 = predefined bitmap applies
!                254 = previously defined bitmap applies to this field
!                255 = bit map does not apply to this product.
!     bmap()   - logical*1 array containing bitmap to be added. 
!                ( if ibmap=0 or ibmap=254)
!
!   output argument list:      
!     cgrib    - character array to contain the grib2 message
!     ierr     - error return code.
!                0 = no error
!                1 = grib message was not initialized.  need to call
!                    routine gribcreate first.
!                2 = grib message already complete.  cannot add new section.
!                3 = sum of section byte counts does not add to total 
!                    byte count.
!                4 = previous section was not 3 or 7.
!                5 = could not find requested product definition template.
!                6 = section 3 (gds) not previously defined in message
!                7 = tried to use unsupported data representationi template
!                8 = specified use of a previously defined bitmap, but one
!                    does not exist in the grib message.
!                9 = gdt of one of 5.50 through 5.53 required to pack
!                    using drt 5.51
!               10 = error packing data field.
!
! remarks: note that the local use section ( section 2 ) can only follow
!          section 1 or section 7 in a grib2 message.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      use pdstemplates
      use drstemplates

      character(len=1),intent(inout) :: cgrib(lcgrib)
      integer,intent(in) :: ipdsnum,ipdstmpl(*)
      integer,intent(in) :: idrsnum,numcoord,ipdstmplen,idrstmplen
      integer,intent(in) :: lcgrib,ngrdpts,ibmap
      real,intent(in) :: coordlist(numcoord)
      real,target,intent(in) :: fld(ngrdpts)
      integer,intent(out) :: ierr
      integer,intent(inout) :: idrstmpl(*)
      logical*1,intent(in) :: bmap(ngrdpts)
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4):: ctemp
      character(len=1),allocatable :: cpack(:)
      real,pointer,dimension(:) :: pfld
      real(4) :: coordieee(numcoord),re00
      integer(4) :: ire00,allones
      integer :: mappds(ipdstmplen),intbmap(ngrdpts),mapdrs(idrstmplen)
      integer,parameter :: zero=0,one=1,four=4,five=5,six=6,seven=7
      integer iofst,ibeg,lencurr,len,mappdslen,mapdrslen,lpos3
      integer width,height,ndpts
      integer lensec3,lensec4,lensec5,lensec6,lensec7
      logical issec3,needext,isprevbmap
 
      ierr=0
      do jj=0,31
         allones=ibset(allones,jj)
      enddo
!
!  check to see if beginning of grib message exists
!
      ctemp=cgrib(1)//cgrib(2)//cgrib(3)//cgrib(4)
      if ( ctemp.ne.grib ) then
        print *,'addfield: grib not found in given message.'
        print *,'addfield: call to routine gribcreate required',
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
        print *,'addfield: grib message already complete.  cannot',
     &          ' add new section.'
        ierr=2
        return
      endif
!
!  loop through all current sections of the grib message to
!  find the last section number.
!
      issec3=.false.
      isprevbmap=.false.
      len=16    ! length of section 0
      do 
      !    get number and length of next section
        iofst=len*8
        call gbyte(cgrib,ilen,iofst,32)
        iofst=iofst+32
        call gbyte(cgrib,isecnum,iofst,8)
        iofst=iofst+8
      !  check if previous section 3 exists and save location of
      !  the section 3 in case needed later.
        if (isecnum.eq.3) then
           issec3=.true.
           lpos3=len+1
           lensec3=ilen
        endif
      !  check if a previous defined bitmap exists
        if (isecnum.eq.6) then
          call gbyte(cgrib,ibmprev,iofst,8)
          iofst=iofst+8
          if ((ibmprev.ge.0).and.(ibmprev.le.253)) isprevbmap=.true.
        endif
        len=len+ilen
      !    exit loop if last section reached
        if ( len.eq.lencurr ) exit
      !    if byte count for each section does not match current
      !    total length, then there is a problem.
        if ( len.gt.lencurr ) then
          print *,'addfield: section byte counts don''t add to total.'
          print *,'addfield: sum of section byte counts = ',len
          print *,'addfield: total byte count in section 0 = ',lencurr
          ierr=3
          return
        endif
      enddo
!
!  sections 4 through 7 can only be added after section 3 or 7.
!
      if ( (isecnum.ne.3) .and. (isecnum.ne.7) ) then
        print *,'addfield: sections 4-7 can only be added after',
     &          ' section 3 or 7.'
        print *,'addfield: section ',isecnum,' was the last found in',
     &          ' given grib message.'
        ierr=4
        return
!
!  sections 4 through 7 can only be added if section 3 was previously defined.
!
      elseif (.not.issec3) then
        print *,'addfield: sections 4-7 can only be added if section',
     &          ' 3 was previously included.'
        print *,'addfield: section 3 was not found in',
     &          ' given grib message.'
        print *,'addfield: call to routine addgrid required',
     &          ' to specify grid definition.'
        ierr=6
        return
      endif
!
!  add section 4  - product definition section
!
      ibeg=lencurr*8        !   calculate offset for beginning of section 4
      iofst=ibeg+32         !   leave space for length of section
      call sbyte(cgrib,four,iofst,8)     ! store section number ( 4 )
      iofst=iofst+8
      call sbyte(cgrib,numcoord,iofst,16)   ! store num of coordinate values
      iofst=iofst+16
      call sbyte(cgrib,ipdsnum,iofst,16)    ! store prod def template num.
      iofst=iofst+16
      !
      !   get product definition template
      !
      call getpdstemplate(ipdsnum,mappdslen,mappds,needext,iret)
      if (iret.ne.0) then
        ierr=5
        return
      endif
      !
      !   extend the product definition template, if necessary.
      !   the number of values in a specific template may vary
      !   depending on data specified in the "static" part of the
      !   template.
      !
      if ( needext ) then
        call extpdstemplate(ipdsnum,ipdstmpl,mappdslen,mappds)
      endif
      !
      !   pack up each input value in array ipdstmpl into the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mappds.
      !
      do i=1,mappdslen
        nbits=iabs(mappds(i))*8
        if ( (mappds(i).ge.0).or.(ipdstmpl(i).ge.0) ) then
          call sbyte(cgrib,ipdstmpl(i),iofst,nbits)
        else
          call sbyte(cgrib,one,iofst,1)
          call sbyte(cgrib,iabs(ipdstmpl(i)),iofst+1,nbits-1)
        endif
        iofst=iofst+nbits
      enddo
      !
      !   add optional list of vertical coordinate values
      !   after the product definition template, if necessary.
      !
      if ( numcoord .ne. 0 ) then
        call mkieee(coordlist,coordieee,numcoord)
        call sbytes(cgrib,coordieee,iofst,32,0,numcoord)
        iofst=iofst+(32*numcoord)
      endif
      !
      !   calculate length of section 4 and store it in octets
      !   1-4 of section 4.
      !
      lensec4=(iofst-ibeg)/8
      call sbyte(cgrib,lensec4,ibeg,32)
!
!  pack data using appropriate algorithm
!
      !
      !   get data representation template
      !
      call getdrstemplate(idrsnum,mapdrslen,mapdrs,needext,iret)
      if (iret.ne.0) then
        ierr=5
        return
      endif
      !
      !  contract data field, removing data at invalid grid points,
      !  if bit-map is provided with field.
      !
      if ( ibmap.eq.0 .or. ibmap.eq.254 ) then
         allocate(pfld(ngrdpts))
         ndpts=0;
         do jj=1,ngrdpts
             intbmap(jj)=0
             if ( bmap(jj) ) then
                intbmap(jj)=1
                ndpts=ndpts+1
                pfld(ndpts)=fld(jj);
             endif
         enddo
      else 
         ndpts=ngrdpts;
         pfld=>fld;
      endif
      lcpack=0
      allocate(cpack(ndpts*4),stat=istat)
      if (idrsnum.eq.0) then      !  simple packing
        call simpack(pfld,ndpts,idrstmpl,cpack,lcpack)
      elseif (idrsnum.eq.2.or.idrsnum.eq.3) then      !  complex packing
        call cmplxpack(pfld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
      elseif (idrsnum.eq.50) then      !  sperical harmonic simple packing
        call simpack(pfld(2),ndpts-1,idrstmpl,cpack,lcpack)
        call mkieee(real(pfld(1)),re00,1)  ! ensure re(0,0) value is ieee format
        !call gbyte(re00,idrstmpl(5),0,32)
        ire00=transfer(re00,ire00)
        idrstmpl(5)=ire00
      elseif (idrsnum.eq.51) then      !  sperical harmonic complex packing
           call getpoly(cgrib(lpos3),lensec3,jj,kk,mm)
           if (jj.ne.0 .and. kk.ne.0 .and. mm.ne.0) then
             call specpack(pfld,ndpts,jj,kk,mm,idrstmpl,cpack,lcpack)
           else
             print *,'addfield: cannot pack drt 5.51.'
             ierr=9
             return
           endif
#ifdef use_jpeg2000
      elseif (idrsnum.eq.40 .or. idrsnum.eq.40000) then     !  jpeg2000 encoding
        if (ibmap.eq.255) then
           call getdim(cgrib(lpos3),lensec3,width,height,iscan)
           if ( width.eq.0 .or. height.eq.0 ) then
              width=ndpts
              height=1
           elseif ( width.eq.allones .or. height.eq.allones ) then
              width=ndpts
              height=1
           elseif ( ibits(iscan,5,1) .eq. 1) then   ! scanning mode: bit 3
              itemp=width
              width=height
              height=itemp
           endif
        else
           width=ndpts
           height=1
        endif
        call jpcpack(pfld,width,height,idrstmpl,cpack,lcpack)
#endif  /* use_jpeg2000 */
#ifdef use_png
      elseif (idrsnum.eq.41 .or. idrsnum.eq.40010) then      !  png encoding
        if (ibmap.eq.255) then
           call getdim(cgrib(lpos3),lensec3,width,height,iscan)
           if ( width.eq.0 .or. height.eq.0 ) then
              width=ndpts
              height=1
           elseif ( width.eq.allones .or. height.eq.allones ) then
              width=ndpts
              height=1
           elseif ( ibits(iscan,5,1) .eq. 1) then   ! scanning mode: bit 3
              itemp=width
              width=height
              height=itemp
           endif
        else
           width=ndpts
           height=1
        endif
        call pngpack(pfld,width,height,idrstmpl,cpack,lcpack)
#endif   /* use_png */
      else
        print *,'addfield: data representation template 5.',idrsnum,
     *          ' not yet implemented.'
        ierr=7
        return
      endif
      if ( ibmap.eq.0 .or. ibmap.eq.254 ) then
         deallocate(pfld)
      endif
      if ( lcpack .lt. 0 ) then
        if( allocated(cpack) )deallocate(cpack)
        ierr=10
        return
      endif

!
!  add section 5  - data representation section
!
      ibeg=iofst            !   calculate offset for beginning of section 5
      iofst=ibeg+32         !   leave space for length of section
      call sbyte(cgrib,five,iofst,8)     ! store section number ( 5 )
      iofst=iofst+8
      call sbyte(cgrib,ndpts,iofst,32)    ! store num of actual data points
      iofst=iofst+32
      call sbyte(cgrib,idrsnum,iofst,16)    ! store data repr. template num.
      iofst=iofst+16
      !
      !   pack up each input value in array idrstmpl into the
      !   the appropriate number of octets, which are specified in
      !   corresponding entries in array mapdrs.
      !
      do i=1,mapdrslen
        nbits=iabs(mapdrs(i))*8
        if ( (mapdrs(i).ge.0).or.(idrstmpl(i).ge.0) ) then
          call sbyte(cgrib,idrstmpl(i),iofst,nbits)
        else
          call sbyte(cgrib,one,iofst,1)
          call sbyte(cgrib,iabs(idrstmpl(i)),iofst+1,nbits-1)
        endif
        iofst=iofst+nbits
      enddo
      !
      !   calculate length of section 5 and store it in octets
      !   1-4 of section 5.
      !
      lensec5=(iofst-ibeg)/8
      call sbyte(cgrib,lensec5,ibeg,32)

!
!  add section 6  - bit-map section
!
      ibeg=iofst            !   calculate offset for beginning of section 6
      iofst=ibeg+32         !   leave space for length of section
      call sbyte(cgrib,six,iofst,8)     ! store section number ( 6 )
      iofst=iofst+8
      call sbyte(cgrib,ibmap,iofst,8)    ! store bit map indicator
      iofst=iofst+8
      !
      !  store bitmap, if supplied
      !
      if (ibmap.eq.0) then
        call sbytes(cgrib,intbmap,iofst,1,0,ngrdpts)    ! store bitmap
        iofst=iofst+ngrdpts
      endif
      !
      !  if specifying a previously defined bit-map, make sure
      !  one already exists in the current grib message.
      !
      if ((ibmap.eq.254).and.(.not.isprevbmap)) then
        print *,'addfield: requested previously defined bitmap, ',
     &        ' but one does not exist in the current grib message.'
        ierr=8
        return
      endif
      !
      !   calculate length of section 6 and store it in octets
      !   1-4 of section 6.  pad to end of octect, if necessary.
      !
      left=8-mod(iofst,8)
      if (left.ne.8) then
        call sbyte(cgrib,zero,iofst,left)     ! pad with zeros to fill octet
        iofst=iofst+left
      endif
      lensec6=(iofst-ibeg)/8
      call sbyte(cgrib,lensec6,ibeg,32)

!
!  add section 7  - data section
!
      ibeg=iofst            !   calculate offset for beginning of section 7
      iofst=ibeg+32         !   leave space for length of section
      call sbyte(cgrib,seven,iofst,8)     ! store section number ( 7 )
      iofst=iofst+8
      !      store packed binary data values, if non-constant field
      if (lcpack.ne.0) then
        ioctet=iofst/8           
        cgrib(ioctet+1:ioctet+lcpack)=cpack(1:lcpack)
        iofst=iofst+(8*lcpack)
      endif
      !
      !   calculate length of section 7 and store it in octets
      !   1-4 of section 7.  
      !
      lensec7=(iofst-ibeg)/8
      call sbyte(cgrib,lensec7,ibeg,32)

      if( allocated(cpack) )deallocate(cpack)
!
!  update current byte total of message in section 0
!
      newlen=lencurr+lensec4+lensec5+lensec6+lensec7
      call sbyte(cgrib,newlen,96,32)

      return
      end

