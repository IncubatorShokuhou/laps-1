c-----------------------------------------------------------------------
      subroutine ixgb(lugb,lskip,lgrib,nlen,nnum,mlen,cbuf)
c$$$  subprogram documentation block
c
c subprogram: ixgb           make index record
c   prgmmr: iredell          org: w/nmc23     date: 95-10-31
c
c abstract: this subprogram makes one index record.
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
c usage:    call wrgi1r(lugb,lskip,lgrib,lugi)
c   input arguments:
c     lugb         integer logical unit of input grib file
c     lskip        integer number of bytes to skip before grib message
c     lgrib        integer number of bytes in grib message
c     nlen         integer length of each index record in bytes
c     nnum         integer index record number to make
c   output arguments:
c     mlen         integer actual valid length of index record
c     cbuf         character*1 (mbuf) buffer to receive index data
c
c subprograms called:
c   gbyte        get integer data from bytes
c   sbyte        store integer data in bytes
c   baread       byte-addressable read
c
c attributes:
c   language: cray fortran
c
c$$$
      character cbuf(*)
      parameter(lindex=112,mindex=320)
      parameter(ixskp=0,ixspd=4,ixsgd=8,ixsbm=12,ixsbd=16,ixlen=20,
     &          ixver=24,ixpds=25,ixgds=53,ixbms=95,ixbds=101,
     &          ixpdx=112,ixpdw=172,ixgdx=184)
      parameter(mxskp=4,mxspd=4,mxsgd=4,mxsbm=4,mxsbd=4,mxlen=4,
     &          mxver=1,mxpds=28,mxgds=42,mxbms=6,mxbds=11,
     &          mxpdx=60,mxpdw=12,mxgdx=136)
      character cbread(mindex),cindex(mindex)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  initialize index record and read grib message
      mlen=lindex
      cindex=char(0)
      call sbyte(cindex,lskip,8*ixskp,8*mxskp)
      call sbyte(cindex,lgrib,8*ixlen,8*mxlen)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  put pds in index record
      iskpds=8
      ibskip=lskip
      ibread=iskpds+mxpds
      call baread(lugb,ibskip,ibread,lbread,cbread)
      if(lbread.ne.ibread) return
      cindex(ixver+1)=cbread(8)
      call sbyte(cindex,iskpds,8*ixspd,8*mxspd)
      call gbyte(cbread,lenpds,8*iskpds,8*3)
      call gbyte(cbread,incgds,8*iskpds+8*7+0,1)
      call gbyte(cbread,incbms,8*iskpds+8*7+1,1)
      ilnpds=min(lenpds,mxpds)
      cindex(ixpds+1:ixpds+ilnpds)=cbread(iskpds+1:iskpds+ilnpds)
      isktot=iskpds+lenpds
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  put pds extension in index record
      if(lenpds.gt.mxpds) then
        iskpdw=iskpds+mxpds
        ilnpdw=min(lenpds-mxpds,mxpdw)
        ibskip=lskip+iskpdw
        ibread=ilnpdw
        call baread(lugb,ibskip,ibread,lbread,cbread)
        if(lbread.ne.ibread) return
        cindex(ixpdw+1:ixpdw+ilnpdw)=cbread(1:ilnpdw)
        iskpdx=iskpds+(mxpds+mxpdw)
        ilnpdx=min(lenpds-(mxpds+mxpdw),mxpdx)
        ibskip=lskip+iskpdx
        ibread=ilnpdx
        call baread(lugb,ibskip,ibread,lbread,cbread)
        if(lbread.ne.ibread) return
        cindex(ixpdx+1:ixpdx+ilnpdx)=cbread(1:ilnpdx)
        mlen=max(mlen,ixpdw+ilnpdw,ixpdx+ilnpdx)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  put gds in index record
      if(incgds.ne.0) then
        iskgds=isktot
        ibskip=lskip+iskgds
        ibread=mxgds
        call baread(lugb,ibskip,ibread,lbread,cbread)
        if(lbread.ne.ibread) return
        call sbyte(cindex,iskgds,8*ixsgd,8*mxsgd)
        call gbyte(cbread,lengds,0,8*3)
        ilngds=min(lengds,mxgds)
        cindex(ixgds+1:ixgds+ilngds)=cbread(1:ilngds)
        isktot=iskgds+lengds
        if(lengds.gt.mxgds) then
          iskgdx=iskgds+mxgds
          ilngdx=min(lengds-mxgds,mxgdx)
          ibskip=lskip+iskgdx
          ibread=ilngdx
          call baread(lugb,ibskip,ibread,lbread,cbread)
          if(lbread.ne.ibread) return
          cindex(ixgdx+1:ixgdx+ilngdx)=cbread(1:ilngdx)
          mlen=max(mlen,ixgdx+ilngdx)
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  put bms in index record
      if(incbms.ne.0) then
        iskbms=isktot
        ibskip=lskip+iskbms
        ibread=mxbms
        call baread(lugb,ibskip,ibread,lbread,cbread)
        if(lbread.ne.ibread) return
        call sbyte(cindex,iskbms,8*ixsbm,8*mxsbm)
        call gbyte(cbread,lenbms,0,8*3)
        ilnbms=min(lenbms,mxbms)
        cindex(ixbms+1:ixbms+ilnbms)=cbread(1:ilnbms)
        isktot=iskbms+lenbms
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  put bds in index record
      iskbds=isktot
      ibskip=lskip+iskbds
      ibread=mxbds
      call baread(lugb,ibskip,ibread,lbread,cbread)
      if(lbread.ne.ibread) return
      call sbyte(cindex,iskbds,8*ixsbd,8*mxsbd)
      call gbyte(cbread,lenbds,0,8*3)
      ilnbds=min(lenbds,mxbds)
      cindex(ixbds+1:ixbds+ilnbds)=cbread(1:ilnbds)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  store index record
      mlen=min(mlen,nlen)
      nskip=nlen*(nnum-1)
      cbuf(nskip+1:nskip+mlen)=cindex(1:mlen)
      cbuf(nskip+mlen+1:nskip+nlen)=char(0)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
