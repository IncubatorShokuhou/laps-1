c-----------------------------------------------------------------------
      subroutine getg2ir(lugb,msk1,msk2,mnum,cbuf,nlen,nnum,nmess,iret)
c$$$  subprogram documentation block
c
c subprogram: getg2ir        creates an index of a grib2 file
c   prgmmr: gilbert          org: w/np11      date: 2002-01-02
c
c abstract: read a grib file and return its index contents.
c   the index buffer returned contains index records with the internal format:
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
c 2002-01-02  gilbert   modified from getgir to create grib2 indexes
c
c usage:    call getg2ir(lugb,msk1,msk2,mnum,cbuf,nlen,nnum,nmess,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib file
c     msk1         integer number of bytes to search for first message
c     msk2         integer number of bytes to search for other messages
c     mnum         integer number of grib messages to skip (usually 0)
c   output arguments:
c     cbuf         character*1 pointer to a buffer that contains index records.
c                  users should free memory that cbuf points to
c                  using deallocate(cbuf) when cbuf is no longer needed.
c     nlen         integer total length of index record buffer in bytes
c     nnum         integer number of index records
c                  (=0 if no grib messages are found)
c     nmess        last grib message in file successfully processed
c     iret         integer return code
c                    0      all ok
c                    1      not enough memory available to hold full index 
c                           buffer
c                    2      not enough memory to allocate initial index buffer
c
c subprograms called:
c   skgb           seek next grib message
c   ixgb2          make index record
c
c remarks: subprogram can be called from a multiprocessing environment.
c   do not engage the same logical unit from more than one processor.
c
c attributes:
c   language: fortran 90
c
c$$$
      use re_alloc          ! needed for subroutine realloc
      parameter(init=50000,next=10000)
      character(len=1),pointer,dimension(:) :: cbuf
      integer,intent(in) :: lugb,msk1,msk2,mnum
      integer,intent(out) :: nlen,nnum,nmess,iret
      character(len=1),pointer,dimension(:) :: cbuftmp
      interface      ! required for cbuf pointer
         subroutine ixgb2(lugb,lskip,lgrib,cbuf,numfld,mlen,iret)
           integer,intent(in) :: lugb,lskip,lgrib
           character(len=1),pointer,dimension(:) :: cbuf
           integer,intent(out) :: numfld,mlen,iret
         end subroutine ixgb2
      end interface
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  initialize
      iret=0
      if (associated(cbuf)) nullify(cbuf)
      mbuf=init
      allocate(cbuf(mbuf),stat=istat)    ! allocate initial space for cbuf
      if (istat.ne.0) then
         iret=2
         return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  search for first grib message
      iseek=0
      call skgb(lugb,iseek,msk1,lskip,lgrib)
      do m=1,mnum
        if(lgrib.gt.0) then
          iseek=lskip+lgrib
          call skgb(lugb,iseek,msk2,lskip,lgrib)
        endif
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  get index records for every grib message found
      nlen=0
      nnum=0
      nmess=mnum
      dowhile(iret.eq.0.and.lgrib.gt.0)
        call ixgb2(lugb,lskip,lgrib,cbuftmp,numfld,nbytes,iret1)
        if (iret1.ne.0) print *,' sagt ',numfld,nbytes,iret1
        if((nbytes+nlen).gt.mbuf) then             ! allocate more space, if
                                                   ! necessary
           newsize=max(mbuf+next,mbuf+nbytes)
           call realloc(cbuf,nlen,newsize,istat)
           if ( istat .ne. 0 ) then
              iret=1
              return
           endif
           mbuf=newsize
        endif
        !
        !  if index records were returned in cbuftmp from ixgb2,
        !  copy cbuftmp into cbuf, then deallocate cbuftmp when done
        !
        if ( associated(cbuftmp) ) then
           cbuf(nlen+1:nlen+nbytes)=cbuftmp(1:nbytes)
           deallocate(cbuftmp,stat=istat)
           if (istat.ne.0) then
             print *,' deallocating cbuftmp ... ',istat
             stop 99
           endif
           nullify(cbuftmp)
           nnum=nnum+numfld
           nlen=nlen+nbytes
           nmess=nmess+1
        endif
        !      look for next grib message
        iseek=lskip+lgrib
        call skgb(lugb,iseek,msk2,lskip,lgrib)
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
