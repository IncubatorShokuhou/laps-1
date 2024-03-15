c-----------------------------------------------------------------------
      subroutine putgb2(lugb,gfld,iret)
c$$$  subprogram documentation block
c
c subprogram: putgb2         packs and writes a grib2 message
c   prgmmr: gilbert          org: w/np11     date: 2002-04-22
c
c abstract: packs a single field into a grib2 message
c   and writes out that message to the file associated with unit lugb.
c   note that file/unit lugb should be opened woth a call to 
c   subroutine baopenw before this routine is called.
c
c   the information to be packed into the grib field
c   is stored in a derived type variable, gfld.
c   gfld is of type gribfield, which is defined
c   in module grib_mod, so users of this routine will need to include
c   the line "use grib_mod" in their calling routine.  each component of the
c   gribfield type is described in the input argument list section below.
c
c program history log:
c 2002-04-22  gilbert  
c 2005-02-28  gilbert   - changed dimension of array cgrib to be a multiple
c                         of gfld%ngrdpts instead of gfld%ndpts.
c
c usage:    call putgb2(lugb,gfld,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib data file.
c                  file must be opened with baopen or baopenw before calling 
c                  this routine.
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
c
c   output arguments:
c     iret         integer return code
c                    0      all ok
c                    2      memory allocation error
c                    10     no section 1 info available
c                    11     no grid definition template info available
c                    12     missing some required data field info
c
c subprograms called:
c   gribcreate     start a new grib2 message
c   addlocal       add local section to a grib2 message
c   addgrid        add grid info to a grib2 message
c   addfield       add data field to a grib2 message
c   gribend        end grib2 message
c
c remarks: 
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

      integer,intent(in) :: lugb
      type(gribfield),intent(in) :: gfld
      integer,intent(out) :: iret

      character(len=1),allocatable,dimension(:) :: cgrib
      integer :: listsec0(2)=(/0,2/)
      integer :: igds(5)=(/0,0,0,0,0/)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  allocate array for grib2 field
      lcgrib=gfld%ngrdpts*4
      allocate(cgrib(lcgrib),stat=is)
      if ( is.ne.0 ) then
         print *,'putgb2: cannot allocate memory. ',is
         iret=2
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  create new message
      listsec0(1)=gfld%discipline
      listsec0(2)=gfld%version
      if ( associated(gfld%idsect) ) then
         call gribcreate(cgrib,lcgrib,listsec0,gfld%idsect,ierr)
         if (ierr.ne.0) then
            write(6,*) 'putgb2: error creating new grib2 field = ',ierr
         endif
      else
         print *,'putgb2: no section 1 info available. '
         iret=10
         deallocate(cgrib)
         return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  add local use section to grib2 message
      if ( associated(gfld%local).and.gfld%locallen.gt.0 ) then
         call addlocal(cgrib,lcgrib,gfld%local,gfld%locallen,ierr)
         if (ierr.ne.0) then
            write(6,*) 'putgb2: error adding local info = ',ierr
         endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  add grid to grib2 message
      igds(1)=gfld%griddef
      igds(2)=gfld%ngrdpts
      igds(3)=gfld%numoct_opt
      igds(4)=gfld%interp_opt
      igds(5)=gfld%igdtnum
      if ( associated(gfld%igdtmpl) ) then
         call addgrid(cgrib,lcgrib,igds,gfld%igdtmpl,gfld%igdtlen,
     &                gfld%list_opt,gfld%num_opt,ierr)
         if (ierr.ne.0) then
            write(6,*) 'putgb2: error adding grid info = ',ierr
         endif
      else
         print *,'putgb2: no gdt info available. '
         iret=11
         deallocate(cgrib)
         return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  add data field to grib2 message
      if ( associated(gfld%ipdtmpl).and.
     &     associated(gfld%idrtmpl).and.
     &     associated(gfld%fld) ) then
         call addfield(cgrib,lcgrib,gfld%ipdtnum,gfld%ipdtmpl,
     &                 gfld%ipdtlen,gfld%coord_list,gfld%num_coord,
     &                 gfld%idrtnum,gfld%idrtmpl,gfld%idrtlen,
     &                 gfld%fld,gfld%ngrdpts,gfld%ibmap,gfld%bmap,
     &                 ierr)
         if (ierr.ne.0) then
            write(6,*) 'putgb2: error adding data field = ',ierr
         endif
      else
         print *,'putgb2: missing some field info. '
         iret=12
         deallocate(cgrib)
         return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  close grib2 message and write to file
      call gribend(cgrib,lcgrib,lengrib,ierr)
      call wryte(lugb,lengrib,cgrib)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      deallocate(cgrib)
      return
      end
