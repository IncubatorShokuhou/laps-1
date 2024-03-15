c-----------------------------------------------------------------------
      subroutine getgb2r(lugb,cindex,gfld,iret)
c$$$  subprogram documentation block
c
c subprogram: getgb2r        reads and unpacks a grib field
c   prgmmr: gilbert          org: w/np11     date: 02-01-15
c
c abstract: read and unpack sections 6 and 7 from a grib2 message.
c
c   this routine assumes that the "metadata" for this field
c   already exists in derived type gribfield.  specifically,
c   it requires gfld%ibmap,gfld%ngrdpts,gfld%idrtnum,gfld%idrtmpl,
c   and gfld%ndpts.
c   
c   the decoded information for the selected grib field
c   is returned in a derived type variable, gfld.
c   gfld is of type gribfield, which is defined
c   in module grib_mod, so users of this routine will need to include
c   the line "use grib_mod" in their calling routine.  each component of the
c   gribfield type is described in the output argument list section below.
c
c program history log:
c   95-10-31  iredell
c 2002-01-11  gilbert     modified from getgb1r to work with grib2
c
c usage:    call getgb2r(lugb,cindex,gfld,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib data file
c     cindex       index record of the grib field  ( see docblock of
c                  subroutine ixgb2 for description of an index record.)
c   output arguments:
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
c        gfld%expanded = logical value indicating whether the data field
c                         was expanded to the grid in the case where a
c                         bit-map is present.  if true, the data points in
c                         gfld%fld match the grid points and zeros were
c                         inserted at grid points where data was bit-mapped
c                         out.  if false, the data values in gfld%fld were
c                         not expanded to the grid and are just a consecutive
c                         array of data points corresponding to each value of
c                         "1" in gfld%bmap.
c        gfld%ibmap = bitmap indicator ( see code table 6.0 )
c                     0 = bitmap applies and is included in section 6.
c                     1-253 = predefined bitmap applies
c                     254 = previously defined bitmap applies to this field
c                     255 = bit map does not apply to this product.
c        gfld%bmap() = logical*1 array containing decoded bitmap,
c                      if ibmap=0 or ibap=254.  otherwise nullified.
c                      this element is actually a pointer to an array
c                      that holds the data.
c        gfld%fld() = array of gfld%ndpts unpacked data points.
c                     this element is actually a pointer to an array
c                     that holds the data.
c     iret         integer return code
c                    0      all ok
c                    97     error reading grib file
c                    other  gf_getfld grib unpacker return code
c
c subprograms called:
c   baread         byte-addressable read
c   gf_unpack6     unapcks bit_map section
c   gf_unpack7     unapcks data section 
c
c remarks: 
c   do not engage the same logical unit from more than one processor.
c   this subprogram is intended for private use by getgb2 routines only.
c
c   note that derived type gribfield contains pointers to many
c   arrays of data.  the memory for these arrays is allocated
c   when the values in the arrays are set, to help minimize
c   problems with array overloading.  because of this, users
c   are encouraged to free up this memory, when it is no longer
c   needed, by an explicit call to subroutine gf_free.
c   ( i.e.   call gf_free(gfld) )
c
c attributes:
c   language: fortran 90
c
c$$$
      use grib_mod

      integer,intent(in) :: lugb
      character(len=1),intent(in) :: cindex(*)
      integer,intent(out) :: iret
      type(gribfield) :: gfld

      integer :: lskip,skip6,skip7
      character(len=1):: csize(4)
      character(len=1),allocatable :: ctemp(:)
      real,pointer,dimension(:) :: newfld

      interface
         subroutine gf_unpack6(cgrib,lcgrib,iofst,ngpts,ibmap,
     &                         bmap,ierr)
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
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  get info
      nullify(gfld%bmap,gfld%fld)
      iret=0
      call gbyte(cindex,lskip,4*8,4*8)
      call gbyte(cindex,skip6,24*8,4*8)
      call gbyte(cindex,skip7,28*8,4*8)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  read and unpack bit_map, if present
      if ( gfld%ibmap.eq.0.or.gfld%ibmap.eq.254 ) then
         iskip=lskip+skip6
         call baread(lugb,iskip,4,lread,csize)    ! get length of section
         call gbyte(csize,ilen,0,32)
         allocate(ctemp(ilen))
         call baread(lugb,iskip,ilen,lread,ctemp)  ! read in section
         if (ilen.ne.lread) then
            iret=97
            deallocate(ctemp)
            return
         endif
         iofst=0
         call gf_unpack6(ctemp,ilen,iofst,gfld%ngrdpts,idum,
     &                   gfld%bmap,ierr)
         if (ierr.ne.0) then
            iret=98
            deallocate(ctemp)
            return
         endif
         deallocate(ctemp)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  read and unpack data field 
      iskip=lskip+skip7
      call baread(lugb,iskip,4,lread,csize)    ! get length of section
      call gbyte(csize,ilen,0,32)
      allocate(ctemp(ilen))
      call baread(lugb,iskip,ilen,lread,ctemp)  ! read in section
      if (ilen.ne.lread) then
         iret=97
         deallocate(ctemp)
         return
      endif
      iofst=0
      call gf_unpack7(ctemp,ilen,iofst,gfld%igdtnum,gfld%igdtmpl,
     &                   gfld%idrtnum,gfld%idrtmpl,gfld%ndpts,
     &                   gfld%fld,ierr)
      if (ierr.ne.0) then
         iret=98
         deallocate(ctemp)
         return
      endif
      deallocate(ctemp)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  if bitmap is used with this field, expand data field
      !  to grid, if possible.
      if ( gfld%ibmap .ne. 255 .and. associated(gfld%bmap) ) then
            allocate(newfld(gfld%ngrdpts))
            !newfld=0.0
            !newfld=unpack(lgfld%fld,lgfld%bmap,newfld)
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
         gfld%expanded=.true.
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
