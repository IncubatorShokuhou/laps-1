c-----------------------------------------------------------------------
      subroutine getg2i(lugi,cbuf,nlen,nnum,iret)
c$$$  subprogram documentation block
c
c subprogram: getg2i          reads a grib2 index file
c   prgmmr: iredell          org: w/nmc23     date: 95-10-31
c
c abstract: read a grib2 index file and return its contents.
c   version 1 of the index file has the following format:
c     81-byte s.lord header with 'gb2ix1' in columns 42-47 followed by
c     81-byte header with number of bytes to skip before index records,
c     total length in bytes of the index records, number of index records,
c     and grib file basename written in format ('ix1form:',3i10,2x,a40).
c     each following index record corresponds to a grib message
c     and has the internal format:
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
c program history log:
c   95-10-31  iredell
c   96-10-31  iredell   augmented optional definitions to byte 320
c 2002-01-03  gilbert   modified from getgi to work with grib2 
c
c usage:    call getg2i(lugi,cbuf,nlen,nnum,iret)
c   input arguments:
c     lugi         integer unit of the unblocked grib index file
c   output arguments:
c     cbuf         character*1 pointer to a buffer that contains index records.
c                  users should free memory that cbuf points to
c                  using deallocate(cbuf) when cbuf is no longer needed.
c     nlen         integer total length of all index records
c     nnum         integer number of index records
c     iret         integer return code
c                    0      all ok
c                    2      not enough memory to hold index buffer
c                    3      error reading index file buffer
c                    4      error reading index file header
c
c subprograms called:
c   baread         byte-addressable read
c
c remarks: subprogram can be called from a multiprocessing environment.
c   do not engage the same logical unit from more than one processor.
c
c attributes:
c   language: fortran 90
c
c$$$
      character(len=1),pointer,dimension(:) :: cbuf
      integer,intent(in) :: lugi
      integer,intent(out) :: nlen,nnum,iret
      character chead*162
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (associated(cbuf)) nullify(cbuf)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      nlen=0
      nnum=0
      iret=4
      call baread(lugi,0,162,lhead,chead)
      if(lhead.eq.162.and.chead(42:47).eq.'gb2ix1') then
        read(chead(82:162),'(8x,3i10,2x,a40)',iostat=ios) nskp,nlen,nnum
        if(ios.eq.0) then
          
          allocate(cbuf(nlen),stat=istat)    ! allocate space for cbuf
          if (istat.ne.0) then
             iret=2
             return
          endif
          iret=0
          call baread(lugi,nskp,nlen,lbuf,cbuf)
          if(lbuf.ne.nlen) iret=3

        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
