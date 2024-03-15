c-----------------------------------------------------------------------
      subroutine getidx(lugb,lugi,cindex,nlen,nnum,iret)
c$$$  subprogram documentation block
c
c subprogram: getidx         finds, reads or generates a grib2 index 
c   prgmmr: gilbert          org: w/np11     date: 2005-03-15
c
c abstract: finds, reads or generates a grib2 index for the grib2 file
c  associated with unit lugb.  if the index already exists, it is returned.
c  otherwise, the index is (1) read from an existing indexfile associated with
c  unit lugi. or (2) generated from the grib2file lugb ( if lugi=0 ). 
c  users can force a regeneration of an index.  if lugi equals lugb, the index
c  will be regenerated from the data in file lugb.  if lugi is less than
c  zero, then the index is re read from index file abs(lugi).  
c
c program history log:
c 2005-03-15  gilbert
c
c usage:    call getidx(lugb,lugi,cindex,nlen,nnum,iret)
c
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
c
c   output arguments:
c     cindex       character*1 pointer to a buffer that contains index records.
c     nlen         integer total length of all index records
c     nnum         integer number of index records
c     iret         integer return code
c                    0      all ok
c                    90     unit number out of range
c                    96     error reading/creating index file
c
c subprograms called:
c   getg2i          read index file
c   getg2ir         read index buffer from grib file
c
c remarks: 
c
c
c attributes:
c   language: fortran 90
c
c$$$

      integer,intent(in) :: lugb,lugi
      integer,intent(out) :: nlen,nnum,iret
      character(len=1),pointer,dimension(:) :: cindex

      integer,parameter :: maxidx=100
      integer,parameter :: msk1=32000,msk2=4000
 
      type gindex
         integer :: nlen
         integer :: nnum
         character(len=1),pointer,dimension(:) :: cbuf
      end type gindex
     
      type(gindex),save :: idxlist(100)

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
      end interface

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  determine whether index buffer needs to be initialized
      lux=0
      iret=0
      if ( lugb.le.0 .and. lugb.gt.100 ) then
         iret=90
         return
      endif
      if (lugi.eq.lugb) then      ! force regeneration of index from grib2 file
         if ( associated( idxlist(lugb)%cbuf ) ) 
     &                  deallocate(idxlist(lugb)%cbuf)
         nullify(idxlist(lugb)%cbuf)
         idxlist(lugb)%nlen=0
         idxlist(lugb)%nnum=0
         lux=0
      endif
      if (lugi.lt.0) then      ! force re-read of index from indexfile
                               ! associated with unit abs(lugi)
         if ( associated( idxlist(lugb)%cbuf ) ) 
     &                  deallocate(idxlist(lugb)%cbuf)
         nullify(idxlist(lugb)%cbuf)
         idxlist(lugb)%nlen=0
         idxlist(lugb)%nnum=0
         lux=abs(lugi)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  check if index already exists in memory
      if ( associated( idxlist(lugb)%cbuf ) ) then
         cindex => idxlist(lugb)%cbuf
         nlen = idxlist(lugb)%nlen
         nnum = idxlist(lugb)%nnum
         return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      irgi=0
      if(lux.gt.0) then
        call getg2i(lux,idxlist(lugb)%cbuf,nlen,nnum,irgi)
      elseif(lux.le.0) then
        mskp=0
        call getg2ir(lugb,msk1,msk2,mskp,idxlist(lugb)%cbuf,
     &               nlen,nnum,nmess,irgi)
      endif
      if(irgi.eq.0) then
         cindex => idxlist(lugb)%cbuf
         idxlist(lugb)%nlen = nlen
         idxlist(lugb)%nnum = nnum
      else
         nlen = 0
         nnum = 0
         iret=96
         return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
