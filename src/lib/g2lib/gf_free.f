      subroutine gf_free(gfld)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    gf_free 
!   prgmmr: gilbert         org: w/np11    date: 2000-05-26
!
! abstract: this subroutine frees up memory that was used to store
!   array values in derived type gribfield.
!
! program history log:
! 2000-05-26  gilbert
!
! usage:    call gf_free(gfld)
!   input argument list:
!     gfld - derived type gribfield ( defined in module grib_mod )
!
!   output argument list:      
!     gfld - derived type gribfield ( defined in module grib_mod )
!        gfld%version = grib edition number
!        gfld%discipline = message discipline ( see code table 0.0 )
!        gfld%idsect() = contains the entries in the identification
!                        section ( section 1 )
!                        this element is actually a pointer to an array
!                        that holds the data.
!            gfld%idsect(1)  = identification of originating centre
!                                    ( see common code table c-1 )
!            gfld%idsect(2)  = identification of originating sub-centre
!            gfld%idsect(3)  = grib master tables version number
!                                    ( see code table 1.0 )
!            gfld%idsect(4)  = grib local tables version number
!                                    ( see code table 1.1 )
!            gfld%idsect(5)  = significance of reference time (code table 1.2)
!            gfld%idsect(6)  = year ( 4 digits )
!            gfld%idsect(7)  = month
!            gfld%idsect(8)  = day
!            gfld%idsect(9)  = hour
!            gfld%idsect(10)  = minute
!            gfld%idsect(11)  = second
!            gfld%idsect(12)  = production status of processed data
!                                    ( see code table 1.3 )
!            gfld%idsect(13)  = type of processed data ( see code table 1.4 )
!        gfld%idsectlen = number of elements in gfld%idsect().
!        gfld%local() = pointer to character array containing contents
!                       of local section 2, if included
!        gfld%locallen = length of array gfld%local()
!        gfld%ifldnum = field number within grib message
!        gfld%griddef = source of grid definition (see code table 3.0)
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
!                        data values were unpacked.  if false, gfld%ndpts
!                        is set to zero, and gfld%bmap and gfld%fld 
!                        pointers are nullified.
!        gfld%ibmap = bitmap indicator ( see code table 6.0 )
!                     0 = bitmap applies and is included in section 6.
!                     1-253 = predefined bitmap applies
!                     254 = previously defined bitmap applies to this field
!                     255 = bit map does not apply to this product.
!        gfld%bmap() - logical*1 array containing decoded bitmap, 
!                      if ibmap=0 or ibap=254.  otherwise nullified.
!                      this element is actually a pointer to an array
!                      that holds the data.
!        gfld%fld() = array of gfld%ndpts unpacked data points.
!                     this element is actually a pointer to an array
!                     that holds the data.
!
! remarks: 
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
      use grib_mod
    
      type(gribfield) :: gfld

      if (associated(gfld%idsect)) then
         deallocate(gfld%idsect)
         !deallocate(gfld%idsect,stat=is)
         !print *,'gfld%idsect: ',is
      endif
      nullify(gfld%idsect)

      if (associated(gfld%local)) then
         deallocate(gfld%local)
         !deallocate(gfld%local,stat=is)
         !print *,'gfld%local: ',is
      endif
      nullify(gfld%local)

      if (associated(gfld%list_opt)) then
         deallocate(gfld%list_opt)
         !deallocate(gfld%list_opt,stat=is)
         !print *,'gfld%list_opt: ',is
      endif
      nullify(gfld%list_opt)

      if (associated(gfld%igdtmpl)) then
         deallocate(gfld%igdtmpl)
         !deallocate(gfld%igdtmpl,stat=is)
         !print *,'gfld%igdtmpl: ',is
      endif
      nullify(gfld%igdtmpl)

      if (associated(gfld%ipdtmpl)) then
         deallocate(gfld%ipdtmpl)
         !deallocate(gfld%ipdtmpl,stat=is)
         !print *,'gfld%ipdtmpl: ',is
      endif
      nullify(gfld%ipdtmpl)

      if (associated(gfld%coord_list)) then
         deallocate(gfld%coord_list)
         !deallocate(gfld%coord_list,stat=is)
         !print *,'gfld%coord_list: ',is
      endif
      nullify(gfld%coord_list)

      if (associated(gfld%idrtmpl)) then
         deallocate(gfld%idrtmpl)
         !deallocate(gfld%idrtmpl,stat=is)
         !print *,'gfld%idrtmpl: ',is
      endif
      nullify(gfld%idrtmpl)

      if (associated(gfld%bmap)) then
         deallocate(gfld%bmap)
         !deallocate(gfld%bmap,stat=is)
         !print *,'gfld%bmap: ',is
      endif
      nullify(gfld%bmap)

      if (associated(gfld%fld)) then
         deallocate(gfld%fld)
         !deallocate(gfld%fld,stat=is)
         !print *,'gfld%fld: ',is
      endif
      nullify(gfld%fld)

      return
      end
