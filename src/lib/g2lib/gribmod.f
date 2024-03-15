      module grib_mod
!$$$  subprogram documentation block
!                .      .    .                                       .
! module:    grib_mod 
!   prgmmr: gilbert         org: w/np11    date: 2002-01-23
!
! abstract: this fortran module contains the declaration
!   of derived type gribfield.
!   if variable gfld is declared of type gribfield 
!   ( i.e. type(gribfield) :: gfld ), it would have the following componenets:
!
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
!
!
! program history log:
! 2002-01-23  gilbert
!
! usage:    use grib_mod
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      character(len=12) :: g2_version="g2lib-1.0.7"

      type gribfield
          integer :: version,discipline
          integer,pointer,dimension(:) :: idsect
          integer :: idsectlen
          character(len=1),pointer,dimension(:) :: local
          integer :: locallen
          integer :: ifldnum
          integer :: griddef,ngrdpts
          integer :: numoct_opt,interp_opt,num_opt
          integer,pointer,dimension(:) :: list_opt
          integer :: igdtnum,igdtlen
          integer,pointer,dimension(:) :: igdtmpl
          integer :: ipdtnum,ipdtlen
          integer,pointer,dimension(:) :: ipdtmpl
          integer :: num_coord
          real,pointer,dimension(:) :: coord_list
          integer :: ndpts,idrtnum,idrtlen
          integer,pointer,dimension(:) :: idrtmpl
          logical :: unpacked
          logical :: expanded
          integer :: ibmap
          logical*1,pointer,dimension(:) :: bmap
          real,pointer,dimension(:) :: fld
      end type gribfield

      end module
