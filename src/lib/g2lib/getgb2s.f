c-----------------------------------------------------------------------
      subroutine getgb2s(cbuf,nlen,nnum,j,jdisc,jids,jpdtn,jpdt,jgdtn,
     &                   jgdt,k,gfld,lpos,iret)
c$$$  subprogram documentation block
c
c subprogram: getgb2s        finds a grib message
c   prgmmr: gilbert          org: w/np11     date: 02-01-15
c
c abstract: find a grib message.
c   find in the index file a reference to the grib field requested.
c   the grib field request specifies the number of messages to skip
c   and the unpacked identification section, grid definition template and 
c   product defintion section parameters.  (a requested parameter
c   of -9999 means to allow any value of this parameter to be found.)
c
c           each index record has the following form:
c       byte 001 - 004: length of index record
c       byte 005 - 008: bytes to skip in data file before grib message
c       byte 009 - 012: bytes to skip in message before lus (local use)
c                       set = 0, if no local use section in grib2 message.
c       byte 013 - 016: bytes to skip in message before gds
c       byte 017 - 020: bytes to skip in message before pds
c       byte 021 - 024: bytes to skip in message before drs
c       byte 025 - 028: bytes to skip in message before bms
c       byte 029 - 032: bytes to skip in message before data section
c       byte 033 - 040: bytes total in the message
c       byte 041 - 041: grib version number ( currently 2 )
c       byte 042 - 042: message discipline
c       byte 043 - 044: field number within grib2 message
c       byte 045 -  ii: identification section (ids)
c       byte ii+1-  jj: grid definition section (gds)
c       byte jj+1-  kk: product definition section (pds)
c       byte kk+1-  ll: the data representation section (drs)
c       byte ll+1-ll+6: first 6 bytes of the bit map section (bms)
c
c   most of the decoded information for the selected grib field
c   is returned in a derived type variable, gfld.  
c   gfld is of type gribfield, which is defined
c   in module grib_mod, so users of this routine will need to include
c   the line "use grib_mod" in their calling routine.  each component of the
c   gribfield type is described in the output argument list section below.
c   only the unpacked bitmap and data field components are not set by this 
c   routine.
c
c program history log:
c   95-10-31  iredell
c 2002-01-02  gilbert   modified from getg1s to work with grib2
c
c usage:    call getgb2s(cbuf,nlen,nnum,j,jdisc,jids,jpdtn,jpdt,jgdtn,
c    &                   jgdt,k,gfld,lpos,iret)
c   input arguments:
c     cbuf         character*1 (nlen) buffer containing index data
c     nlen         integer total length of all index records
c     nnum         integer number of index records
c     j            integer number of messages to skip
c                  (=0 to search from beginning)
c     jdisc        grib2 discipline number of requested field
c                  ( if = -1, accept any discipline)
c                  ( see code table 0.0 )
c                  0 - meteorological products
c                  1 - hydrological products
c                  2 - land surface products
c                  3 - space products
c                  10 - oceanographic products
c     jids()       integer array of values in the identification section
c                  (=-9999 for wildcard)
c            jids(1)   = identification of originating centre
c                         ( see common code table c-1 )
c            jids(2)   = identification of originating sub-centre
c            jids(3)   = grib master tables version number
c                         ( see code table 1.0 )
c                       0 - experimental
c                       1 - initial operational version number
c            jids(4)   = grib local tables version number
c                         ( see code table 1.1 )
c                       0     - local tables not used
c                       1-254 - number of local tables version used
c            jids(5)   = significance of reference time (code table 1.2)
c                       0 - analysis
c                       1 - start of forecast
c                       2 - verifying time of forecast
c                       3 - observation time
c            jids(6)   = year ( 4 digits )
c            jids(7)   = month
c            jids(8)   = day
c            jids(9)   = hour
c            jids(10)  = minute
c            jids(11)  = second
c            jids(12)  = production status of processed data
c                         ( see code table 1.3 )
c                       0 - operational products
c                       1 - operational test products
c                       2 - research products
c                       3 - re-analysis products
c            jids(13)  = type of processed data ( see code table 1.4 )
c                       0  - analysis products
c                       1  - forecast products
c                       2  - analysis and forecast products
c                       3  - control forecast products
c                       4  - perturbed forecast products
c                       5  - control and perturbed forecast products
c                       6  - processed satellite observations
c                       7  - processed radar observations
c     jpdtn        integer product definition template number (n)
c                  ( if = -1, don't bother matching pdt )
c     jpdt()       integer array of values defining the product definition 
c                  template 4.n of the field for which to search
c                  (=-9999 for wildcard)
c     jgdtn        integer grid definition template number (m)
c                  ( if = -1, don't bother matching gdt )
c     jgdt()       integer array of values defining the grid definition
c                  template 3.m of the field for which to search
c                  (=-9999 for wildcard)
c   output arguments:
c     k            integer message number found
c                  (can be same as j in calling program
c                  in order to facilitate multiple searches)
c     gfld - derived type gribfield ( defined in module grib_mod )
c            ( note: see remarks section )
c        gfld%version = grib edition number ( currently 2 )
c        gfld%discipline = message discipline ( see code table 0.0 )
c        gfld%idsect() = contains the entries in the identification
c                        section ( section 1 )
c                        this element is actually a pointer to an array
c                        that holds the data.
c            gfld%idsect(1)  = identification of originating centre
c                                    ( see common code table c-1 )
c                             7 - us national weather service
c            gfld%idsect(2)  = identification of originating sub-centre
c            gfld%idsect(3)  = grib master tables version number
c                                    ( see code table 1.0 )
c                             0 - experimental
c                             1 - initial operational version number
c            gfld%idsect(4)  = grib local tables version number
c                                    ( see code table 1.1 )
c                             0     - local tables not used
c                             1-254 - number of local tables version used
c            gfld%idsect(5)  = significance of reference time (code table 1.2)
c                             0 - analysis
c                             1 - start of forecast
c                             2 - verifying time of forecast
c                             3 - observation time
c            gfld%idsect(6)  = year ( 4 digits )
c            gfld%idsect(7)  = month
c            gfld%idsect(8)  = day
c            gfld%idsect(9)  = hour
c            gfld%idsect(10)  = minute
c            gfld%idsect(11)  = second
c            gfld%idsect(12)  = production status of processed data
c                                    ( see code table 1.3 )
c                              0 - operational products
c                              1 - operational test products
c                              2 - research products
c                              3 - re-analysis products
c            gfld%idsect(13)  = type of processed data ( see code table 1.4 )
c                              0  - analysis products
c                              1  - forecast products
c                              2  - analysis and forecast products
c                              3  - control forecast products
c                              4  - perturbed forecast products
c                              5  - control and perturbed forecast products
c                              6  - processed satellite observations
c                              7  - processed radar observations
c        gfld%idsectlen = number of elements in gfld%idsect().
c        gfld%local() = pointer to character array containing contents
c                       of local section 2, if included
c        gfld%locallen = length of array gfld%local()
c        gfld%ifldnum = field number within grib message
c        gfld%griddef = source of grid definition (see code table 3.0)
c                      0 - specified in code table 3.1
c                      1 - predetermined grid defined by originating centre
c        gfld%ngrdpts = number of grid points in the defined grid.
c        gfld%numoct_opt = number of octets needed for each
c                          additional grid points definition.
c                          used to define number of
c                          points in each row ( or column ) for
c                          non-regular grids.
c                          = 0, if using regular grid.
c        gfld%interp_opt = interpretation of list for optional points
c                          definition.  (code table 3.11)
c        gfld%igdtnum = grid definition template number (code table 3.1)
c        gfld%igdtmpl() = contains the data values for the specified grid
c                         definition template ( nn=gfld%igdtnum ).  each
c                         element of this integer array contains an entry (in
c                         the order specified) of grid defintion template 3.nn
c                         this element is actually a pointer to an array
c                         that holds the data.
c        gfld%igdtlen = number of elements in gfld%igdtmpl().  i.e. number of
c                       entries in grid defintion template 3.nn
c                       ( nn=gfld%igdtnum ).
c        gfld%list_opt() = (used if gfld%numoct_opt .ne. 0)  this array
c                          contains the number of grid points contained in
c                          each row ( or column ).  (part of section 3)
c                          this element is actually a pointer to an array
c                          that holds the data.  this pointer is nullified
c                          if gfld%numoct_opt=0.
c        gfld%num_opt = (used if gfld%numoct_opt .ne. 0)  the number of entries
c                       in array ideflist.  i.e. number of rows ( or columns )
c                       for which optional grid points are defined.  this value
c                       is set to zero, if gfld%numoct_opt=0.
c        gfdl%ipdtnum = product definition template number (see code table 4.0)
c        gfld%ipdtmpl() = contains the data values for the specified product
c                         definition template ( n=gfdl%ipdtnum ).  each element
c                         of this integer array contains an entry (in the
c                         order specified) of product defintion template 4.n.
c                         this element is actually a pointer to an array
c                         that holds the data.
c        gfld%ipdtlen = number of elements in gfld%ipdtmpl().  i.e. number of
c                       entries in product defintion template 4.n
c                       ( n=gfdl%ipdtnum ).
c        gfld%coord_list() = real array containing floating point values
c                            intended to document the vertical discretisation
c                            associated to model data on hybrid coordinate
c                            vertical levels.  (part of section 4)
c                            this element is actually a pointer to an array
c                            that holds the data.
c        gfld%num_coord = number of values in array gfld%coord_list().
c        gfld%ndpts = number of data points unpacked and returned.
c        gfld%idrtnum = data representation template number
c                       ( see code table 5.0)
c        gfld%idrtmpl() = contains the data values for the specified data
c                         representation template ( n=gfld%idrtnum ).  each
c                         element of this integer array contains an entry
c                         (in the order specified) of product defintion
c                         template 5.n.
c                         this element is actually a pointer to an array
c                         that holds the data.
c        gfld%idrtlen = number of elements in gfld%idrtmpl().  i.e. number
c                       of entries in data representation template 5.n
c                       ( n=gfld%idrtnum ).
c        gfld%unpacked = logical value indicating whether the bitmap and
c                        data values were unpacked.  if false,
c                        gfld%bmap and gfld%fld pointers are nullified.
c                        note: this routine sets this component to .false.
c        gfld%ibmap = bitmap indicator ( see code table 6.0 )
c                     0 = bitmap applies and is included in section 6.
c                     1-253 = predefined bitmap applies
c                     254 = previously defined bitmap applies to this field
c                     255 = bit map does not apply to this product.
c        gfld%bmap() = logical*1 array containing decoded bitmap,
c                      if ibmap=0 or ibap=254.  otherwise nullified.
c                      this element is actually a pointer to an array
c                      that holds the data.
c                      note: this component is not set by this routine.
c        gfld%fld() = array of gfld%ndpts unpacked data points.
c                     this element is actually a pointer to an array
c                     that holds the data.
c                      note: this component is not set by this routine.
c     lpos         starting position of the found index record within
c                  the complete index buffer, cbuf.
c                  = 0, if request not found
c     iret         integer return code
c                    0      all ok
c                    1      request not found
c
c remarks: 
c   this subprogram is intended for private use by getgb2 routines only.
c
c   note that derived type gribfield contains pointers to many
c   arrays of data.  the memory for these arrays is allocated
c   when the values in the arrays are set, to help minimize
c   problems with array overloading.  because of this users
c   are encouraged to free up this memory, when it is no longer
c   needed, by an explicit call to subroutine gf_free.
c   ( i.e.   call gf_free(gfld) )
c
c subprograms called:
c   gbyte            unpack bytes
c   gf_unpack1          unpack ids
c   gf_unpack4          unpack pds
c   gf_unpack3          unpack gds
c
c attributes:
c   language: fortran 90
c
c$$$
      use grib_mod

!      character(len=1),pointer,dimension(:) :: cbuf
      character(len=1),intent(in) :: cbuf(nlen)
      integer,intent(in) :: nlen,nnum,j,jdisc,jpdtn,jgdtn
      integer,dimension(:) :: jids(*),jpdt(*),jgdt(*)
      integer,intent(out) :: k,lpos,iret
      type(gribfield),intent(out) :: gfld

      integer :: kgds(5)
      logical :: match1,match3,match4
!      integer,pointer,dimension(:) :: kids,kpdt,kgdt
!      integer,pointer,dimension(:) :: idef
!      real,pointer,dimension(:) :: coord

      interface
         subroutine gf_unpack1(cgrib,lcgrib,iofst,ids,idslen,ierr)
            character(len=1),intent(in) :: cgrib(lcgrib)
            integer,intent(in) :: lcgrib
            integer,intent(inout) :: iofst
            integer,pointer,dimension(:) :: ids
            integer,intent(out) :: ierr,idslen
         end subroutine gf_unpack1
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
      end interface
      
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  initialize
      k=0
      lpos=0
      iret=1
      ipos=0
      nullify(gfld%list_opt,gfld%igdtmpl,gfld%ipdtmpl)
      nullify(gfld%coord_list,gfld%idrtmpl,gfld%bmap,gfld%fld)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  search for request
      dowhile(iret.ne.0.and.k.lt.nnum)
        k=k+1
        call gbyte(cbuf,inlen,ipos*8,4*8)    ! get length of current
                                              ! index record
        if ( k.le.j ) then           ! skip this index
           ipos=ipos+inlen
           cycle
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  check if grib2 discipline is a match
        call gbyte(cbuf,gfld%discipline,(ipos+41)*8,1*8)
        if ( (jdisc.ne.-1).and.(jdisc.ne.gfld%discipline) ) then
           ipos=ipos+inlen
           cycle
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  check if identification section is a match
        match1=.false.
        call gbyte(cbuf,lsec1,(ipos+44)*8,4*8)  ! get length of ids 
        iof=0
        call gf_unpack1(cbuf(ipos+45),lsec1,iof,gfld%idsect,
     &                  gfld%idsectlen,icnd)
        if ( icnd.eq.0 ) then
           match1=.true.
           do i=1,gfld%idsectlen
              if ( (jids(i).ne.-9999).and.
     &             (jids(i).ne.gfld%idsect(i)) ) then
                 match1=.false.
                 exit
              endif
           enddo
        endif
        if ( .not. match1 ) then
           deallocate(gfld%idsect)
           ipos=ipos+inlen
           cycle
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  check if grid definition template is a match
        jpos=ipos+44+lsec1
        match3=.false.
        call gbyte(cbuf,lsec3,jpos*8,4*8)  ! get length of gds 
        if ( jgdtn.eq.-1 ) then
           match3=.true.
        else
           call gbyte(cbuf,numgdt,(jpos+12)*8,2*8)  ! get gdt template no.
           if ( jgdtn.eq.numgdt ) then
              iof=0
              call gf_unpack3(cbuf(jpos+1),lsec3,iof,kgds,gfld%igdtmpl,
     &                     gfld%igdtlen,gfld%list_opt,gfld%num_opt,icnd)
              if ( icnd.eq.0 ) then
                 match3=.true.
                 do i=1,gfld%igdtlen
                    if ( (jgdt(i).ne.-9999).and.
     &                   (jgdt(i).ne.gfld%igdtmpl(i)) ) then
                       match3=.false.
                       exit
                    endif
                 enddo
c                 where ( jgdt(1:gfld%igdtlen).ne.-9999 ) 
c     &              match3=all(jgdt(1:gfld%igdtlen).eq.gfld%igdtmpl(1:gfld%igdtlen))
              endif
           endif
        endif
        if ( .not. match3 ) then
           if (associated(gfld%igdtmpl)) deallocate(gfld%igdtmpl)
           if (associated(gfld%list_opt)) deallocate(gfld%list_opt)
           ipos=ipos+inlen
           cycle
        else
           gfld%griddef=kgds(1)
           gfld%ngrdpts=kgds(2)
           gfld%numoct_opt=kgds(3)
           gfld%interp_opt=kgds(4)
           gfld%igdtnum=kgds(5)
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  check if product definition template is a match
        jpos=jpos+lsec3
        match4=.false.
        call gbyte(cbuf,lsec4,jpos*8,4*8)  ! get length of pds 
        if ( jpdtn.eq.-1 ) then
           match4=.true.
        else
           call gbyte(cbuf,numpdt,(jpos+7)*8,2*8)  ! get pdt template no.
           if ( jpdtn.eq.numpdt ) then
              iof=0
              call gf_unpack4(cbuf(jpos+1),lsec4,iof,gfld%ipdtnum,
     &                        gfld%ipdtmpl,gfld%ipdtlen,
     &                        gfld%coord_list,gfld%num_coord,icnd)
              if ( icnd.eq.0 ) then
                 match4=.true.
                 do i=1,gfld%ipdtlen
                    if ( (jpdt(i).ne.-9999).and.
     &                   (jpdt(i).ne.gfld%ipdtmpl(i)) ) then
                       match4=.false.
                       exit
                    endif
                 enddo
c                 where ( jpdt.ne.-9999) 
c     &              match4=all( jpdt(1:gfld%ipdtlen) .eq. gfld%ipdtmpl(1:gfld%ipdtlen) )
              endif
           endif
        endif
        if ( .not. match4 ) then
           if (associated(gfld%ipdtmpl)) deallocate(gfld%ipdtmpl)
           if (associated(gfld%coord_list)) deallocate(gfld%coord_list)
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  if request is found
c  set values for derived type gfld and return
        if(match1.and.match3.and.match4) then
           lpos=ipos+1
           call gbyte(cbuf,gfld%version,(ipos+40)*8,1*8)
           call gbyte(cbuf,gfld%ifldnum,(ipos+42)*8,2*8)
           gfld%unpacked=.false.
           jpos=ipos+44+lsec1
           if ( jgdtn.eq.-1 ) then     ! unpack gds, if not done before
              iof=0
              call gf_unpack3(cbuf(jpos+1),lsec3,iof,kgds,gfld%igdtmpl,
     &                     gfld%igdtlen,gfld%list_opt,gfld%num_opt,icnd)
              gfld%griddef=kgds(1)
              gfld%ngrdpts=kgds(2)
              gfld%numoct_opt=kgds(3)
              gfld%interp_opt=kgds(4)
              gfld%igdtnum=kgds(5)
           endif
           jpos=jpos+lsec3
           if ( jpdtn.eq.-1 ) then     ! unpack pds, if not done before
              iof=0
              call gf_unpack4(cbuf(jpos+1),lsec4,iof,gfld%ipdtnum,
     &                        gfld%ipdtmpl,gfld%ipdtlen,
     &                        gfld%coord_list,gfld%num_coord,icnd)
           endif
           jpos=jpos+lsec4
           call gbyte(cbuf,lsec5,jpos*8,4*8)  ! get length of drs 
           iof=0
           call gf_unpack5(cbuf(jpos+1),lsec5,iof,gfld%ndpts,
     &                     gfld%idrtnum,gfld%idrtmpl,
     &                     gfld%idrtlen,icnd)
           jpos=jpos+lsec5
           call gbyte(cbuf,gfld%ibmap,(jpos+5)*8,1*8)  ! get ibmap
           iret=0
        else      ! pdt did not match
           ipos=ipos+inlen
        endif
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
