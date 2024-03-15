c-----------------------------------------------------------------------
      subroutine getgi(lugi,mnum,mbuf,cbuf,nlen,nnum,iret)
c$$$  subprogram documentation block
c
c subprogram: getgi          reads a grib index file
c   prgmmr: iredell          org: w/nmc23     date: 95-10-31
c
c abstract: read a grib index file and return its contents.
c   version 1 of the index file has the following format:
c     81-byte s.lord header with 'gb1ix1' in columns 42-47 followed by
c     81-byte header with number of bytes to skip before index records,
c     number of bytes in each index record, number of index records,
c     and grib file basename written in format ('ix1form:',3i10,2x,a40).
c     each following index record corresponds to a grib message
c     and has the internal format:
c       byte 001-004: bytes to skip in data file before grib message
c       byte 005-008: bytes to skip in message before pds
c       byte 009-012: bytes to skip in message before gds (0 if no gds)
c       byte 013-016: bytes to skip in message before bms (0 if no bms)
c       byte 017-020: bytes to skip in message before bds
c       byte 021-024: bytes total in the message
c       byte 025-025: grib version number
c       byte 026-053: product definition section (pds)
c       byte 054-095: grid definition section (gds) (or nulls)
c       byte 096-101: first part of the bit map section (bms) (or nulls)
c       byte 102-112: first part of the binary data section (bds)
c       byte 113-172: (optional) bytes 41-100 of the pds
c       byte 173-184: (optional) bytes 29-40 of the pds
c       byte 185-320: (optional) bytes 43-178 of the gds
c
c program history log:
c   95-10-31  iredell
c   96-10-31  iredell   augmented optional definitions to byte 320
c
c usage:    call getgi(lugi,mnum,mbuf,cbuf,nlen,nnum,iret)
c   input arguments:
c     lugi         integer unit of the unblocked grib index file
c     mnum         integer number of index records to skip (usually 0)
c     mbuf         integer length of cbuf in bytes
c   output arguments:
c     cbuf         character*1 (mbuf) buffer to receive index data
c     nlen         integer length of each index record in bytes
c     nnum         integer number of index records
c     iret         integer return code
c                    0      all ok
c                    1      cbuf too small to hold index buffer
c                    2      error reading index file buffer
c                    3      error reading index file header
c
c subprograms called:
c   baread         byte-addressable read
c
c remarks: subprogram can be called from a multiprocessing environment.
c   do not engage the same logical unit from more than one processor.
c
c attributes:
c   language: fortran 77
c   machine:  cray, workstations
c
c$$$
      character cbuf(mbuf)
      character chead*162
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      nlen=0
      nnum=0
      iret=3
      call baread(lugi,0,162,lhead,chead)
      if(lhead.eq.162.and.chead(42:47).eq.'gb1ix1') then
        read(chead(82:162),'(8x,3i10,2x,a40)',iostat=ios) nskp,nlen,nnum
        if(ios.eq.0) then
          nskp=nskp+mnum*nlen
          nnum=nnum-mnum
          nbuf=nnum*nlen
          iret=0
          if(nbuf.gt.mbuf) then
            nnum=mbuf/nlen
            nbuf=nnum*nlen
            iret=1
          endif
          if(nbuf.gt.0) then
            call baread(lugi,nskp,nbuf,lbuf,cbuf)
            if(lbuf.ne.nbuf) iret=2
          endif
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
