c-----------------------------------------------------------------------
      subroutine getgir(lugb,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
c$$$  subprogram documentation block
c
c subprogram: getgir         reads a grib index file
c   prgmmr: iredell          org: w/nmc23     date: 95-10-31
c
c abstract: read a grib file and return its index contents.
c   the index buffer returned contains index records with the internal format:
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
c usage:    call getgir(lugb,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib file
c     msk1         integer number of bytes to search for first message
c     msk2         integer number of bytes to search for other messages
c     mnum         integer number of index records to skip (usually 0)
c     mbuf         integer length of cbuf in bytes
c   output arguments:
c     cbuf         character*1 (mbuf) buffer to receive index data
c     nlen         integer length of each index record in bytes
c     nnum         integer number of index records
c                  (=0 if no grib messages are found)
c     iret         integer return code
c                    0      all ok
c                    1      cbuf too small to hold index data
c
c subprograms called:
c   skgb           seek next grib message
c   ixgb           make index record
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
      parameter(mindex=320)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  search for first grib message
      iseek=0
      call skgb(lugb,iseek,msk1,lskip,lgrib)
      if(lgrib.gt.0.and.mindex.le.mbuf) then
        call ixgb(lugb,lskip,lgrib,mindex,1,nlen,cbuf)
      else
        nlen=mindex
      endif
      do m=1,mnum
        if(lgrib.gt.0) then
          iseek=lskip+lgrib
          call skgb(lugb,iseek,msk2,lskip,lgrib)
        endif
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  make an index record for every grib record found
      nnum=0
      iret=0
      dowhile(iret.eq.0.and.lgrib.gt.0)
        if(nlen*(nnum+1).le.mbuf) then
          nnum=nnum+1
          call ixgb(lugb,lskip,lgrib,nlen,nnum,mlen,cbuf)
          iseek=lskip+lgrib
          call skgb(lugb,iseek,msk2,lskip,lgrib)
        else
          iret=1
        endif
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
