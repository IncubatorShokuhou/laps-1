c-----------------------------------------------------------------------
      subroutine ixgb2(lugb,lskip,lgrib,cbuf,numfld,mlen,iret)
c$$$  subprogram documentation block
c
c subprogram: ixgb2          make index records for fields in a grib2 message
c   prgmmr: gilbert          org: w/np11      date: 2001-12-10
c
c abstract: this subprogram generates an index record for each field in a
c           grib2 message.  the index records are written to index buffer
c           pointed to by cbuf.
c
c           each index record has the following form:
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
c 2001-12-10  gilbert   modified from ixgb to create grib2 indexes
c 2002-01-31  gilbert   added identification section to index record
c
c usage:    call ixgb2(lugb,lskip,lgrib,cbuf,numfld,mlen,iret)
c   input arguments:
c     lugb         integer logical unit of input grib file
c     lskip        integer number of bytes to skip before grib message
c     lgrib        integer number of bytes in grib message
c   output arguments:
c     cbuf         character*1 pointer to a buffer that contains index records.
c                  users should free memory that cbuf points to
c                  using deallocate(cbuf) when cbuf is no longer needed.
c     numfld       integer number of index records created.
c                  = 0, if problems
c     mlen         integer total length of all index records
c     iret         integer return code
c                  =0, all ok
c                  =1, not enough memory to allocate initial index buffer
c                  =2, i/o error in read
c                  =3, grib message is not edition 2
c                  =4, not enough memory to allocate extent to index buffer
c                  =5, unidentified grib section encountered...problem 
c                      somewhere.
c
c subprograms called:
c   gbyte        get integer data from bytes
c   sbyte        store integer data in bytes
c   baread       byte-addressable read
c   realloc      re-allocates more memory
c
c attributes:
c   language: fortran 90
c
c$$$
      use re_alloc          ! needed for subroutine realloc
      character(len=1),pointer,dimension(:) :: cbuf
      parameter(linmax=5000,init=50000,next=10000)
      parameter(ixskp=4,ixlus=8,ixsgd=12,ixspd=16,ixsdr=20,ixsbm=24,
     &          ixds=28,ixlen=36,ixfld=42,ixids=44)
      parameter(mxskp=4,mxlus=4,mxsgd=4,mxspd=4,mxsdr=4,mxsbm=4,
     &          mxds=4,mxlen=4,mxfld=2,mxbms=6)
      character cbread(linmax),cindex(linmax)
      character cver,cdisc
      character cids(linmax),cgds(linmax),cbms(6)
      character(len=4) :: ctemp
      integer loclus,locgds,lengds,locbms
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      loclus=0
      iret=0
      mlen=0
      numfld=0
      if (associated(cbuf)) nullify(cbuf)
      mbuf=init
      allocate(cbuf(mbuf),stat=istat)    ! allocate initial space for cbuf
      if (istat.ne.0) then
         iret=1
         return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  read sections 0 and 1 for versin number and discipline
      ibread=min(lgrib,linmax)
      call baread(lugb,lskip,ibread,lbread,cbread)
      if(lbread.ne.ibread) then
         iret=2
         return
      endif
      if(cbread(8).ne.char(2)) then          !  not grib edition 2
         iret=3
         return
      endif
      cver=cbread(8)
      cdisc=cbread(7)
      call gbyte(cbread,lensec1,16*8,4*8)
      lensec1=min(lensec1,ibread)
      cids(1:lensec1)=cbread(17:16+lensec1)
      ibskip=lskip+16+lensec1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  loop through remaining sections creating an index for each field
      ibread=max(5,mxbms)
      do
         call baread(lugb,ibskip,ibread,lbread,cbread)     
         ctemp=cbread(1)//cbread(2)//cbread(3)//cbread(4)
         if (ctemp.eq.'7777') return        ! end of message found
         if(lbread.ne.ibread) then
            iret=2
            return
         endif
         call gbyte(cbread,lensec,0*8,4*8)
         call gbyte(cbread,numsec,4*8,1*8)

         if (numsec.eq.2) then                 ! save local use location
            loclus=ibskip-lskip
         elseif (numsec.eq.3) then                 ! save gds info
            lengds=lensec
            cgds=char(0)
            call baread(lugb,ibskip,lengds,lbread,cgds)     
            if(lbread.ne.lengds) then
               iret=2
               return
            endif
            locgds=ibskip-lskip
         elseif (numsec.eq.4) then                 ! found pds
            cindex=char(0)
            call sbyte(cindex,lskip,8*ixskp,8*mxskp)    ! bytes to skip
            call sbyte(cindex,loclus,8*ixlus,8*mxlus)   ! location of local use
            call sbyte(cindex,locgds,8*ixsgd,8*mxsgd)   ! location of gds
            call sbyte(cindex,ibskip-lskip,8*ixspd,8*mxspd)  ! location of pds
            call sbyte(cindex,lgrib,8*ixlen,8*mxlen)    ! len of grib2
            cindex(41)=cver
            cindex(42)=cdisc
            call sbyte(cindex,numfld+1,8*ixfld,8*mxfld)   ! field num
            cindex(ixids+1:ixids+lensec1)=cids(1:lensec1)
            lindex=ixids+lensec1
            cindex(lindex+1:lindex+lengds)=cgds(1:lengds)
            lindex=lindex+lengds
            ilnpds=lensec
            call baread(lugb,ibskip,ilnpds,lbread,cindex(lindex+1))     
            if(lbread.ne.ilnpds) then
               iret=2
               return
            endif
            !   cindex(lindex+1:lindex+ilnpds)=cbread(1:ilnpds)
            lindex=lindex+ilnpds
         elseif (numsec.eq.5) then                 ! found drs
            call sbyte(cindex,ibskip-lskip,8*ixsdr,8*mxsdr)  ! location of drs
            ilndrs=lensec
            call baread(lugb,ibskip,ilndrs,lbread,cindex(lindex+1))     
            if(lbread.ne.ilndrs) then
               iret=2
               return
            endif
            !   cindex(lindex+1:lindex+ilndrs)=cbread(1:ilndrs)
            lindex=lindex+ilndrs
         elseif (numsec.eq.6) then                 ! found bms
            indbmp=mova2i(cbread(6))
            if ( indbmp.lt.254 ) then
               locbms=ibskip-lskip
               call sbyte(cindex,locbms,8*ixsbm,8*mxsbm)  ! loc. of bms
            elseif ( indbmp.eq.254 ) then
               call sbyte(cindex,locbms,8*ixsbm,8*mxsbm)  ! loc. of bms
            elseif ( indbmp.eq.255 ) then
               call sbyte(cindex,ibskip-lskip,8*ixsbm,8*mxsbm)  ! loc. of bms
            endif
            cindex(lindex+1:lindex+mxbms)=cbread(1:mxbms)
            lindex=lindex+mxbms
            call sbyte(cindex,lindex,0,8*4)    ! num bytes in index record
         elseif (numsec.eq.7) then                 ! found data section
            call sbyte(cindex,ibskip-lskip,8*ixds,8*mxds)   ! loc. of data sec.
            numfld=numfld+1
            if ((lindex+mlen).gt.mbuf) then        ! allocate more space if
                                                   ! necessary
               newsize=max(mbuf+next,mbuf+lindex)
               call realloc(cbuf,mlen,newsize,istat)
               if ( istat .ne. 0 ) then
                  numfld=numfld-1
                  iret=4
                  return
               endif
               mbuf=newsize
            endif
            cbuf(mlen+1:mlen+lindex)=cindex(1:lindex)
            mlen=mlen+lindex
         else                           ! unrecognized section
            iret=5
            return
         endif
         ibskip=ibskip+lensec
      enddo

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
