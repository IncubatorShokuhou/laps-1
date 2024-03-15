c-----------------------------------------------------------------------
      subroutine getgb2p(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     &                   extract,k,gribm,leng,iret)
c$$$  subprogram documentation block
c
c subprogram: getgb2p        finds and extracts a grib message
c   prgmmr: iredell          org: w/nmc23     date: 94-04-01
c
c abstract: find and extracts a grib message from a file.
c   read a grib index file (or optionally the grib file itself)
c   to get the index buffer (i.e. table of contents) for the grib file.
c   find in the index buffer a reference to the grib field requested.
c   the grib field request specifies the number of fields to skip
c   and the unpacked identification section, grid definition template and
c   product defintion section parameters.  (a requested parameter
c   of -9999 means to allow any value of this parameter to be found.)
c   if the requested grib field is found, then it is read from the
c   grib file and returned. 
c   if the grib field is not found, then the return code will be nonzero.
c
c program history log:
c   94-04-01  iredell
c   95-10-31  iredell     modularized portions of code into subprograms
c                         and allowed for unspecified index file
c 2002-01-11  gilbert     modified from getgb and getgbm to work with grib2
c 2003-12-17  gilbert     modified from getgb2 to return packed grib2 message.
c
c usage:    call getgb2p(lugb,lugi,j,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
c    &                  extract,k,gribm,leng,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib data file.
c                  file must be opened with baopen or baopenr before calling 
c                  this routine.
c     lugi         integer unit of the unblocked grib index file.
c                  if nonzero, file must be opened with baopen baopenr before 
c                  calling this routine.
c                  (=0 to get index buffer from the grib file)
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
c     extract       logical value indicating whether to return a grib2 
c                   message with just the requested field, or the entire
c                   grib2 message containing the requested field.
c                  .true. = return grib2 message containing only the requested
c                           field.
c                  .false. = return entire grib2 message containing the
c                            requested field.
c
c   output arguments:
c     k            integer field number returned.
c     gribm         returned grib message.
c     leng         length of returned grib message in bytes.
c     iret         integer return code
c                    0      all ok
c                    96     error reading index file
c                    97     error reading grib file
c                    99     request not found
c
c subprograms called:
c   getg2i          read index file
c   getg2ir         read index buffer from grib file
c   getgb2s        search index records
c   getgb2rp        read a packed grib record
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
      logical,intent(in) :: extract
      integer,intent(out) :: k,iret,leng
      character(len=1),pointer,dimension(:) :: gribm

      type(gribfield) :: gfld

      character(len=1),pointer,dimension(:) :: cbuf
      parameter(msk1=32000,msk2=4000)

      save cbuf,nlen,nnum
      data lux/0/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  declare interfaces (required for cbuf pointer)
      interface
         subroutine getg2i(lugi,cbuf,nlen,nnum,iret)
            character(len=1),pointer,dimension(:) :: cbuf
            integer,intent(in) :: lugi
            integer,intent(out) :: nlen,nnum,iret
         end subroutine getg2i
         subroutine getg2ir(lugb,msk1,msk2,mnum,cbuf,nlen,nnum,
     &                      nmess,iret)
            character(len=1),pointer,dimension(:) :: cbuf
            integer,intent(in) :: lugb,msk1,msk2,mnum
            integer,intent(out) :: nlen,nnum,nmess,iret
         end subroutine getg2ir
         subroutine getgb2rp(lugb,cindex,extract,gribm,leng,iret)
            integer,intent(in) :: lugb
            character(len=1),intent(in) :: cindex(*)
            logical,intent(in) :: extract
            integer,intent(out) :: leng,iret
            character(len=1),pointer,dimension(:) :: gribm
         end subroutine getgb2rp
      end interface

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  determine whether index buffer needs to be initialized
      irgi=0
      if(lugi.gt.0.and.lugi.ne.lux) then
        call getg2i(lugi,cbuf,nlen,nnum,irgi)
        lux=lugi
      elseif(lugi.le.0.and.lugb.ne.lux) then
        mskp=0
        call getg2ir(lugb,msk1,msk2,mskp,cbuf,nlen,nnum,nmess,irgi)
        lux=lugb
      endif
      if(irgi.gt.1) then
        iret=96
        lux=0
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
c  extract grib message from file
      call getgb2rp(lugb,cbuf(lpos:),extract,gribm,leng,iret)
!      if ( extract ) then
!         print *,'not supposed to be here.'
!      else
!         ipos=(lpos+3)*8
!         call gbyte(cbuf,iskip,ipos,32)     ! bytes to skip in file
!         ipos=ipos+(32*8)
!         call gbyte(cbuf,leng,ipos,32)      ! length of grib message
!         if (.not. associated(gribm)) allocate(gribm(leng))
!         call baread(lugb,iskip,leng,lread,gribm)
!         if ( leng .ne. lread ) then
!            iret=97
!            call gf_free(gfld)
!            return
!         endif
!      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      k=jk
      call gf_free(gfld)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
