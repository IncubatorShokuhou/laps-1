c-----------------------------------------------------------------------
      subroutine getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     &                  unpack,k,gfld,iret)
c$$$  subprogram documentation block
c
c subprogram: getgb2         finds and unpacks a grib message
c   prgmmr: iredell          org: w/nmc23     date: 94-04-01
c
c abstract: find and unpack a grib message.
c   read a grib index file (or optionally the grib file itself)
c   to get the index buffer (i.e. table of contents) for the grib file.
c   find in the index buffer a reference to the grib field requested.
c   the grib field request specifies the number of fields to skip
c   and the unpacked identification section, grid definition template and
c   product defintion section parameters.  (a requested parameter
c   of -9999 means to allow any value of this parameter to be found.)
c   if the requested grib field is found, then it is read from the
c   grib file and unpacked.  its number is returned along with
c   the associated unpacked parameters.  the bitmap (if any),
c   and the data values are unpacked only if argument "unpack" is set to
c   true.  if the grib field is not found, then the
c   return code will be nonzero.
c
c   the decoded information for the selected grib field
c   is returned in a derived type variable, gfld.
c   gfld is of type gribfield, which is defined
c   in module grib_mod, so users of this routine will need to include
c   the line "use grib_mod" in their calling routine.  each component of the
c   gribfield type is described in the output argument list section below.
c
c program history log:
c   94-04-01  iredell
c   95-10-31  iredell     modularized portions of code into subprograms
c                         and allowed for unspecified index file
c 2002-01-11  gilbert     modified from getgb and getgbm to work with grib2
c
c usage:    call getgb2(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
c    &                  unpack,k,gfld,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib data file.
c                  file must be opened with baopen or baopenr before calling 
c                  this routine.
c     lugi         integer unit of the unblocked grib index file.
c                  if nonzero, file must be opened with baopen baopenr before 
c                  calling this routine.
c                  >0 - read index from index file lugi, if index doesn"t
c                       already exist.
c                  =0 - to get index buffer from the grib file, if index
c                       doesn"t already exist.
c                  <0 - force reread of index from index file abs(lugi).
c                  =lugb - force regeneration of index from grib2 file lugb.
c     j            integer number of fields to skip
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
c                  ( if = -1, don't bother matching pdt - accept any )
c     jpdt()       integer array of values defining the product definition
c                  template 4.n of the field for which to search
c                  (=-9999 for wildcard)
c     jgdtn        integer grid definition template number (m)
c                  ( if = -1, don't bother matching gdt - accept any )
c     jgdt()       integer array of values defining the grid definition
c                  template 3.m of the field for which to search
c                  (=-9999 for wildcard)
c     unpack       logical value indicating whether to unpack bitmap/data
c                  .true. = unpack bitmap and data values
c                  .false. = do not unpack bitmap and data values
c
c   output arguments:
c     k            integer field number unpacked
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
c                    96     error reading index
c                    97     error reading grib file
c                    99     request not found
c                    other  gf_getfld grib2 unpacker return code
c
c subprograms called:
c   getidx         get index
c   getgb2s        search index records
c   getgb2r        read and unpack grib record
c   gf_free        frees memory used by gfld  ( see remarks )
c
c remarks: specify an index file if feasible to increase speed.
c   do not engage the same logical unit from more than one processor.
c
c   note that derived type gribfield contains pointers to many
c   arrays of data.  the memory for these arrays is allocated
c   when the values in the arrays are set, to help minimize
c   problems with array overloading.  because of this users
c   are encouraged to free up this memory, when it is no longer
c   needed, by an explicit call to subroutine gf_free.
c   ( i.e.   call gf_free(gfld) )
c
c attributes:
c   language: fortran 90
c
c$$$
      use grib_mod

      integer,intent(in) :: lugb,lugi,j,jdisc,jpdtn,jgdtn
      integer,dimension(:) :: jids(*),jpdt(*),jgdt(*)
      logical,intent(in) :: unpack
      integer,intent(out) :: k,iret
      type(gribfield),intent(out) :: gfld

      character(len=1),pointer,dimension(:) :: cbuf

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  declare interfaces (required for cbuf pointer)
      interface
         subroutine getidx(lugb,lugi,cbuf,nlen,nnum,irgi)
            character(len=1),pointer,dimension(:) :: cbuf
            integer,intent(in) :: lugb,lugi
            integer,intent(out) :: nlen,nnum,irgi
         end subroutine getidx
      end interface
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  determine whether index buffer needs to be initialized
      irgi=0
      call getidx(lugb,lugi,cbuf,nlen,nnum,irgi)
      if(irgi.gt.1) then
        iret=96
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  search index buffer
      call getgb2s(cbuf,nlen,nnum,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     &             jk,gfld,lpos,irgs)
      if(irgs.ne.0) then
        iret=99
        call gf_free(gfld)
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  read local use section, if available
      call getgb2l(lugb,cbuf(lpos),gfld,iret)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  read and unpack grib record
      if (unpack) then
    !    numfld=gfld%ifldnum
    !    call gf_free(gfld)
        call getgb2r(lugb,cbuf(lpos),gfld,iret)
      endif
      k=jk
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
