      subroutine gf_getfld(cgrib,lcgrib,ifldnum,unpack,expand,gfld,ierr)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_getfld 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine returns the grid definition, product definition,
!   bit-map ( if applicable ), and the unpacked data for a given data
!   field.  all of the information returned is stored in a derived
!   type variable, gfld.  gfld is of type gribfield, which is defined
!   in module grib_mod, so users of this routine will need to include
!   the line "use grib_mod" in their calling routine.  each component of the 
!   gribfield type is described in the output argument list section below.
!
!   since there can be multiple data fields packed into a grib2
!   message, the calling routine indicates which field is being requested
!   with the ifldnum argument.
!
! program history log:
! 2000-05-26  gilbert
! 2002-01-24  gilbert  - changed to pass back derived type gribfield
!                        variable through argument list, instead of
!                        having many different arguments.
! 2004-05-20  gilbert  - added check to see if previous a bit-map is specified,
!                        but none was found.
!
! usage:    call gf_getfld(cgrib,lcgrib,ifldnum,unpack,expand,gfld,ierr)
!   input argument list:
!     cgrib    - character array that contains the grib2 message
!     lcgrib   - length (in bytes) of grib message array cgrib.
!     ifldnum  - specifies which field in the grib2 message to return.
!     unpack   - logical value indicating whether to unpack bitmap/data
!                .true. = unpack bitmap and data values
!                .false. = do not unpack bitmap and data values
!     expand   - boolean value indicating whether the data points should be
!                expanded to the correspond grid, if a bit-map is present.
!                1 = if possible, expand data field to grid, inserting zero
!                    values at gridpoints that are bitmapped out.
!                    (see remarks2)
!                0 = do not expand data field, leaving it an array of
!                    consecutive data points for each "1" in the bitmap.
!                this argument is ignored if unpack == 0 or if the
!                returned field does not contain a bit-map.
!
!   output argument list:      
!     gfld - derived type gribfield ( defined in module grib_mod )
!            ( note: see remarks section )
!        gfld%version = grib edition number ( currently 2 )
!        gfld%discipline = message discipline ( see code table 0.0 )
!        gfld%idsect() = contains the entries in the identification
!                        section ( section 1 )
!                        this element is actually a pointer to an array
!                        that holds the data.
!            gfld%idsect(1)  = identification of originating centre 
!                                    ( see common code table c-1 )
!                             7 - us national weather service
!            gfld%idsect(2)  = identification of originating sub-centre
!            gfld%idsect(3)  = grib master tables version number
!                                    ( see code table 1.0 )
!                             0 - experimental
!                             1 - initial operational version number
!            gfld%idsect(4)  = grib local tables version number
!                                    ( see code table 1.1 )
!                             0     - local tables not used
!                             1-254 - number of local tables version used
!            gfld%idsect(5)  = significance of reference time (code table 1.2)
!                             0 - analysis
!                             1 - start of forecast
!                             2 - verifying time of forecast
!                             3 - observation time
!            gfld%idsect(6)  = year ( 4 digits )
!            gfld%idsect(7)  = month
!            gfld%idsect(8)  = day
!            gfld%idsect(9)  = hour
!            gfld%idsect(10)  = minute
!            gfld%idsect(11)  = second
!            gfld%idsect(12)  = production status of processed data
!                                    ( see code table 1.3 )
!                              0 - operational products
!                              1 - operational test products
!                              2 - research products
!                              3 - re-analysis products
!            gfld%idsect(13)  = type of processed data ( see code table 1.4 )
!                              0  - analysis products
!                              1  - forecast products
!                              2  - analysis and forecast products
!                              3  - control forecast products
!                              4  - perturbed forecast products
!                              5  - control and perturbed forecast products
!                              6  - processed satellite observations
!                              7  - processed radar observations
!        gfld%idsectlen = number of elements in gfld%idsect().
!        gfld%local() = pointer to character array containing contents
!                       of local section 2, if included
!        gfld%locallen = length of array gfld%local()
!        gfld%ifldnum = field number within grib message
!        gfld%griddef = source of grid definition (see code table 3.0)
!                      0 - specified in code table 3.1
!                      1 - predetermined grid defined by originating centre
!        gfld%ngrdpts = number of grid points in the defined grid.
!        gfld%numoct_opt = number of octets needed for each 
!                          additional grid points definition.  
!                          used to define number of
!                          points in each row ( or column ) for
!                          non-regular grids.  
!                          = 0, if using regular grid.
!        gfld%interp_opt = interpretation of list for optional points 
!                          definition.  (code table 3.11)
!        gfld%igdtnum = grid definition template number (code table 3.1)
!        gfld%igdtmpl() = contains the data values for the specified grid 
!                         definition template ( nn=gfld%igdtnum ).  each 
!                         element of this integer array contains an entry (in 
!                         the order specified) of grid defintion template 3.nn
!                         this element is actually a pointer to an array
!                         that holds the data.
!        gfld%igdtlen = number of elements in gfld%igdtmpl().  i.e. number of
!                       entries in grid defintion template 3.nn  
!                       ( nn=gfld%igdtnum ).
!        gfld%list_opt() = (used if gfld%numoct_opt .ne. 0)  this array 
!                          contains the number of grid points contained in 
!                          each row ( or column ).  (part of section 3)
!                          this element is actually a pointer to an array
!                          that holds the data.  this pointer is nullified
!                          if gfld%numoct_opt=0.
!        gfld%num_opt = (used if gfld%numoct_opt .ne. 0)  the number of entries
!                       in array ideflist.  i.e. number of rows ( or columns )
!                       for which optional grid points are defined.  this value
!                       is set to zero, if gfld%numoct_opt=0.
!        gfdl%ipdtnum = product definition template number (see code table 4.0)
!        gfld%ipdtmpl() = contains the data values for the specified product 
!                         definition template ( n=gfdl%ipdtnum ).  each element
!                         of this integer array contains an entry (in the 
!                         order specified) of product defintion template 4.n.
!                         this element is actually a pointer to an array
!                         that holds the data.
!        gfld%ipdtlen = number of elements in gfld%ipdtmpl().  i.e. number of
!                       entries in product defintion template 4.n  
!                       ( n=gfdl%ipdtnum ).
!        gfld%coord_list() = real array containing floating point values 
!                            intended to document the vertical discretisation
!                            associated to model data on hybrid coordinate
!                            vertical levels.  (part of section 4)
!                            this element is actually a pointer to an array
!                            that holds the data.
!        gfld%num_coord = number of values in array gfld%coord_list().
!        gfld%ndpts = number of data points unpacked and returned.
!        gfld%idrtnum = data representation template number 
!                       ( see code table 5.0)
!        gfld%idrtmpl() = contains the data values for the specified data 
!                         representation template ( n=gfld%idrtnum ).  each 
!                         element of this integer array contains an entry 
!                         (in the order specified) of product defintion 
!                         template 5.n.
!                         this element is actually a pointer to an array
!                         that holds the data.
!        gfld%idrtlen = number of elements in gfld%idrtmpl().  i.e. number 
!                       of entries in data representation template 5.n 
!                       ( n=gfld%idrtnum ).
!        gfld%unpacked = logical value indicating whether the bitmap and
!                        data values were unpacked.  if false, 
!                        gfld%bmap and gfld%fld pointers are nullified.
!        gfld%expanded = logical value indicating whether the data field
!                         was expanded to the grid in the case where a
!                         bit-map is present.  if true, the data points in
!                         gfld%fld match the grid points and zeros were
!                         inserted at grid points where data was bit-mapped
!                         out.  if false, the data values in gfld%fld were
!                         not expanded to the grid and are just a consecutive
!                         array of data points corresponding to each value of
!                         "1" in gfld%bmap.
!        gfld%ibmap = bitmap indicator ( see code table 6.0 )
!                     0 = bitmap applies and is included in section 6.
!                     1-253 = predefined bitmap applies
!                     254 = previously defined bitmap applies to this field
!                     255 = bit map does not apply to this product.
!        gfld%bmap() = logical*1 array containing decoded bitmap, 
!                      if ibmap=0 or ibap=254.  otherwise nullified.
!                      this element is actually a pointer to an array
!                      that holds the data.
!        gfld%fld() = array of gfld%ndpts unpacked data points.
!                     this element is actually a pointer to an array
!                     that holds the data.
!     ierr     - error return code.
!                0 = no error
!                1 = beginning characters "grib" not found.
!                2 = grib message is not edition 2.
!                3 = the data field request number was not positive.
!                4 = end string "7777" found, but not where expected.
!                6 = grib message did not contain the requested number of
!                    data fields.
!                7 = end string "7777" not found at end of message.
!                8 = unrecognized section encountered.
!                9 = data representation template 5.nn not yet implemented.
!               15 = error unpacking section 1.
!               16 = error unpacking section 2.
!               10 = error unpacking section 3.
!               11 = error unpacking section 4.
!               12 = error unpacking section 5.
!               13 = error unpacking section 6.
!               14 = error unpacking section 7.
!               17 = previous bitmap specified, but none exists.
!
! remarks: note that derived type gribfield contains pointers to many
!          arrays of data.  the memory for these arrays is allocated
!          when the values in the arrays are set, to help minimize
!          problems with array overloading.  because of this users
!          are encouraged to free up this memory, when it is no longer
!          needed, by an explicit call to subroutine gf_free.
!          ( i.e.   call gf_free(gfld) )
!
!          subroutine gb_info can be used to first determine
!          how many data fields exist in a given grib message.
!
! remarks2: it may not always be possible to expand a bit-mapped data field.
!           if a pre-defined bit-map is used and not included in the grib2
!           message itself, this routine would not have the necessary
!           information to expand the data.  in this case, gfld%expanded would
!           would be set to 0 (false), regardless of the value of input
!           argument expand.
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
      use grib_mod
    
      character(len=1),intent(in) :: cgrib(lcgrib)
      integer,intent(in) :: lcgrib,ifldnum
      logical,intent(in) :: unpack,expand
      type(gribfield),intent(out) :: gfld
      integer,intent(out) :: ierr
!      integer,intent(out) :: igds(*),igdstmpl(*),ideflist(*)
!      integer,intent(out) :: ipdsnum,ipdstmpl(*)
!      integer,intent(out) :: idrsnum,idrstmpl(*)
!      integer,intent(out) :: ndpts,ibmap,idefnum,numcoord
!      logical*1,intent(out) :: bmap(*)
!      real,intent(out) :: fld(*),coordlist(*)
      
      character(len=4),parameter :: grib='grib',c7777='7777'
      character(len=4) :: ctemp
      real,pointer,dimension(:) :: newfld
      integer:: listsec0(2),igds(5)
      integer iofst,ibeg,istart
      integer(4) :: ieee
      logical*1,pointer,dimension(:) :: bmpsave
      logical have3,have4,have5,have6,have7

      interface
         subroutine gf_unpack1(cgrib,lcgrib,iofst,ids,idslen,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib
            integer,intent(inout) :: iofst
            integer,pointer,dimension(:) :: ids
            integer,intent(out) :: ierr,idslen
         end subroutine gf_unpack1
         subroutine gf_unpack2(cgrib,lcgrib,iofst,lencsec2,csec2,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib
            integer,intent(inout) :: iofst
            integer,intent(out) :: lencsec2
            integer,intent(out) :: ierr
            character(len=1),pointer,dimension(:) :: csec2
         end subroutine gf_unpack2
         subroutine gf_unpack3(cgrib,lcgrib,iofst,igds,igdstmpl,
     &                         mapgridlen,ideflist,idefnum,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib
            integer,intent(inout) :: iofst
            integer,pointer,dimension(:) :: igdstmpl,ideflist
            integer,intent(out) :: igds(5)
            integer,intent(out) :: ierr,idefnum
         end subroutine gf_unpack3
         subroutine gf_unpack4(cgrib,lcgrib,iofst,ipdsnum,ipdstmpl,
     &                      mappdslen,coordlist,numcoord,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib
            integer,intent(inout) :: iofst
            real,pointer,dimension(:) :: coordlist
            integer,pointer,dimension(:) :: ipdstmpl
            integer,intent(out) :: ipdsnum
            integer,intent(out) :: ierr,numcoord
         end subroutine gf_unpack4
         subroutine gf_unpack5(cgrib,lcgrib,iofst,ndpts,idrsnum,
     &                         idrstmpl,mapdrslen,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib
            integer,intent(inout) :: iofst
            integer,intent(out) :: ndpts,idrsnum
            integer,pointer,dimension(:) :: idrstmpl
            integer,intent(out) :: ierr
         end subroutine gf_unpack5
         subroutine gf_unpack6(cgrib,lcgrib,iofst,ngpts,ibmap,bmap,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib,ngpts
            integer,intent(inout) :: iofst
            integer,intent(out) :: ibmap
            integer,intent(out) :: ierr
            logical*1,pointer,dimension(:) :: bmap
         end subroutine gf_unpack6
         subroutine gf_unpack7(cgrib,lcgrib,iofst,igdsnum,igdstmpl,
     &                         idrsnum,idrstmpl,ndpts,fld,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib,ndpts,idrsnum,igdsnum
            integer,intent(inout) :: iofst
            integer,pointer,dimension(:) :: idrstmpl,igdstmpl
            integer,intent(out) :: ierr
            real,pointer,dimension(:) :: fld
         end subroutine gf_unpack7
      end interface

      have3=.false.
      have4=.false.
      have5=.false.
      have6=.false.
      have7=.false.
      ierr=0
      numfld=0
      gfld%locallen=0
      nullify(gfld%list_opt,gfld%igdtmpl,gfld%ipdtmpl)
      nullify(gfld%coord_list,gfld%idrtmpl,gfld%bmap,gfld%fld)
!
!  check for valid request number
!  
      if (ifldnum.le.0) then
        print *,'gf_getfld: request for field number must be positive.'
        ierr=3
        return
      endif
!
!  check for beginning of grib message in the first 100 bytes
!
      istart=0
      do j=1,100
        ctemp=cgrib(j)//cgrib(j+1)//cgrib(j+2)//cgrib(j+3)
        if (ctemp.eq.grib ) then
          istart=j
          exit
        endif
      enddo
      if (istart.eq.0) then
        print *,'gf_getfld:  beginning characters grib not found.'
        ierr=1
        return
      endif
!
!  unpack section 0 - indicator section 
!
      iofst=8*(istart+5)
      call gbyte(cgrib,listsec0(1),iofst,8)     ! discipline
      iofst=iofst+8
      call gbyte(cgrib,listsec0(2),iofst,8)     ! grib edition number
      iofst=iofst+8
      iofst=iofst+32
      call gbyte(cgrib,lengrib,iofst,32)        ! length of grib message
      iofst=iofst+32
      lensec0=16
      ipos=istart+lensec0
!
!  currently handles only grib edition 2.
!  
      if (listsec0(2).ne.2) then
        print *,'gf_getfld: can only decode grib edition 2.'
        ierr=2
        return
      endif
!
!  loop through the remaining sections keeping track of the 
!  length of each.  also keep the latest grid definition section info.
!  unpack the requested field number.
!
      do
        !    check to see if we are at end of grib message
        ctemp=cgrib(ipos)//cgrib(ipos+1)//cgrib(ipos+2)//cgrib(ipos+3)
        if (ctemp.eq.c7777 ) then
          ipos=ipos+4
          !    if end of grib message not where expected, issue error
          if (ipos.ne.(istart+lengrib)) then
            print *,'gf_getfld: "7777" found, but not where expected.'
            ierr=4
            return
          endif
          exit
        endif
        !     get length of section and section number
        iofst=(ipos-1)*8
        call gbyte(cgrib,lensec,iofst,32)        ! get length of section
        iofst=iofst+32
        call gbyte(cgrib,isecnum,iofst,8)         ! get section number
        iofst=iofst+8
        !print *,' lensec= ',lensec,'    secnum= ',isecnum
        !
        !  check to see if section number is valid
        !
        if ( (isecnum.lt.1).or.(isecnum.gt.7) ) then
          print *,'gf_getfld: unrecognized section encountered=',isecnum     
          ierr=8
          return
        endif
        !
        !   if found section 1, decode elements in identification section
        !
        if (isecnum.eq.1) then
          iofst=iofst-40       ! reset offset to beginning of section
          call gf_unpack1(cgrib,lcgrib,iofst,gfld%idsect,
     &                    gfld%idsectlen,jerr)
          if (jerr.ne.0) then
            ierr=15
            return
          endif
        endif
        !
        !   if found section 2, grab local section
        !   save in case this is the latest one before the requested field.
        !
        if (isecnum.eq.2) then
          iofst=iofst-40       ! reset offset to beginning of section
          if (associated(gfld%local)) deallocate(gfld%local)
          call gf_unpack2(cgrib,lcgrib,iofst,gfld%locallen,
     &                    gfld%local,jerr)
          if (jerr.ne.0) then
            ierr=16
            return
          endif
        endif
        !
        !   if found section 3, unpack the gds info using the 
        !   appropriate template.  save in case this is the latest
        !   grid before the requested field.
        !
        if (isecnum.eq.3) then
          iofst=iofst-40       ! reset offset to beginning of section
          if (associated(gfld%igdtmpl)) deallocate(gfld%igdtmpl)
          if (associated(gfld%list_opt)) deallocate(gfld%list_opt)
          call gf_unpack3(cgrib,lcgrib,iofst,igds,gfld%igdtmpl,
     &                 gfld%igdtlen,gfld%list_opt,gfld%num_opt,jerr)
          if (jerr.eq.0) then
            have3=.true.
            gfld%griddef=igds(1)
            gfld%ngrdpts=igds(2)
            gfld%numoct_opt=igds(3)
            gfld%interp_opt=igds(4)
            gfld%igdtnum=igds(5)
          else
            ierr=10
            return
          endif
        endif
        !
        !   if found section 4, check to see if this field is the
        !   one requested.
        !
        if (isecnum.eq.4) then
          numfld=numfld+1
          if (numfld.eq.ifldnum) then
            gfld%discipline=listsec0(1)
            gfld%version=listsec0(2)
            gfld%ifldnum=ifldnum
            gfld%unpacked=unpack
            gfld%expanded=.false.
            iofst=iofst-40       ! reset offset to beginning of section
            call gf_unpack4(cgrib,lcgrib,iofst,gfld%ipdtnum,
     &                      gfld%ipdtmpl,gfld%ipdtlen,gfld%coord_list,
     &                      gfld%num_coord,jerr)
            if (jerr.eq.0) then
              have4=.true.
            else
              ierr=11
              return
            endif
          endif
        endif
        !
        !   if found section 5, check to see if this field is the
        !   one requested.
        !
        if ((isecnum.eq.5).and.(numfld.eq.ifldnum)) then
          iofst=iofst-40       ! reset offset to beginning of section
          call gf_unpack5(cgrib,lcgrib,iofst,gfld%ndpts,gfld%idrtnum,
     &                    gfld%idrtmpl,gfld%idrtlen,jerr)
          if (jerr.eq.0) then
            have5=.true.
          else
            ierr=12
            return
          endif
        endif
        !
        !   if found section 6, unpack bitmap.
        !   save in case this is the latest
        !   bitmap before the requested field.
        !
        if (isecnum.eq.6) then
          if (unpack) then   ! unpack bitmap
            iofst=iofst-40       ! reset offset to beginning of section
            bmpsave=>gfld%bmap      ! save pointer to previous bitmap
            call gf_unpack6(cgrib,lcgrib,iofst,gfld%ngrdpts,gfld%ibmap,
     &                   gfld%bmap,jerr)
            if (jerr.eq.0) then
              have6=.true.
              if (gfld%ibmap .eq. 254) then    ! use previously specified bitmap
                 if ( associated(bmpsave) ) then
                    gfld%bmap=>bmpsave
                 else
                    print *,'gf_getfld:  previous bit-map specified,',
     &                       ' but none exists,'
                    ierr=17
                    return
                 endif
              else                             ! get rid of it
                 if ( associated(bmpsave) ) deallocate(bmpsave)
              endif
            else
              ierr=13
              return
            endif
          else    ! do not unpack bitmap
            call gbyte(cgrib,gfld%ibmap,iofst,8)      ! get bitmap indicator
            have6=.true.
          endif
        endif
        !
        !   if found section 7, check to see if this field is the
        !   one requested.
        !
        if ((isecnum.eq.7).and.(numfld.eq.ifldnum).and.unpack) then
          iofst=iofst-40       ! reset offset to beginning of section
          call gf_unpack7(cgrib,lcgrib,iofst,gfld%igdtnum,
     &                    gfld%igdtmpl,gfld%idrtnum,
     &                    gfld%idrtmpl,gfld%ndpts,
     &                    gfld%fld,jerr)
          if (jerr.eq.0) then
            have7=.true.
            !  if bitmap is used with this field, expand data field
            !  to grid, if possible.
            if ( gfld%ibmap .ne. 255 .and. associated(gfld%bmap) ) then
               if ( expand ) then
                  allocate(newfld(gfld%ngrdpts))
                  !newfld(1:gfld%ngrdpts)=0.0
                  !newfld=unpack(gfld%fld,gfld%bmap,newfld)
                  n=1
                  do j=1,gfld%ngrdpts
                      if ( gfld%bmap(j) ) then
                        newfld(j)=gfld%fld(n)
                        n=n+1
                      else
                        newfld(j)=0.0
                      endif
                  enddo
                  deallocate(gfld%fld);
                  gfld%fld=>newfld;
                  gfld%expanded=.true.
               else 
                  gfld%expanded=.false.
               endif
            else 
               gfld%expanded=.true.
            endif
          else
            print *,'gf_getfld: return from gf_unpack7 = ',jerr
            ierr=14
            return
          endif
        endif
        !
        !   check to see if we read pass the end of the grib
        !   message and missed the terminator string '7777'.
        !
        ipos=ipos+lensec                 ! update beginning of section pointer
        if (ipos.gt.(istart+lengrib)) then
          print *,'gf_getfld: "7777"  not found at end of grib message.'
          ierr=7
          return
        endif
        !
        !  if unpacking requested, return when all sections have been
        !  processed
        !
        if (unpack.and.have3.and.have4.and.have5.and.have6.and.have7)
     &      return
        !
        !  if unpacking is not requested, return when sections 
        !  3 through 6 have been processed
        !
        if ((.not.unpack).and.have3.and.have4.and.have5.and.have6)
     &      return
        
      enddo

!
!  if exited from above loop, the end of the grib message was reached
!  before the requested field was found.
!
      print *,'gf_getfld: grib message contained ',numlocal,
     &        ' different fields.'
      print *,'gf_getfld: the request was for the ',ifldnum,
     &        ' field.'
      ierr=6

      return
      end

