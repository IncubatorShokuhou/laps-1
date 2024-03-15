c-----------------------------------------------------------------------
      subroutine getgb1s(cbuf,nlen,nnum,j,jpds,jgds,jens,
     &                   k,kpds,kgds,kens,lskip,lgrib,iret)
c$$$  subprogram documentation block
c
c subprogram: getgb1s        finds a grib message
c   prgmmr: iredell          org: w/nmc23     date: 95-10-31
c
c abstract: find a grib message.
c   find in the index file a reference to the grib message requested.
c   the grib message request specifies the number of messages to skip
c   and the unpacked pds and gds parameters.  (a requested parameter
c   of -1 means to allow any value of this parameter to be found.)
c
c program history log:
c   95-10-31  iredell
c
c usage:    call getgb1s(cbuf,nlen,nnum,j,jpds,jgds,jens,
c    &                   k,kpds,kgds,kens,lskip,lgrib,iret)
c   input arguments:
c     cbuf         character*1 (nlen*nnum) buffer containing index data
c     nlen         integer length of each index record in bytes
c     nnum         integer number of index records
c     j            integer number of messages to skip
c                  (=0 to search from beginning)
c     jpds         integer (200) pds parameters for which to search
c                  (=-1 for wildcard)
c     jgds         integer (200) gds parameters for which to search
c                  (only searched if jpds(3)=255)
c                  (=-1 for wildcard)
c     jens         integer (200) ensemble pds parms for which to search
c                  (only searched if jpds(23)=2)
c                  (=-1 for wildcard)
c   output arguments:
c     k            integer message number found
c                  (can be same as j in calling program
c                  in order to facilitate multiple searches)
c     kpds         integer (200) unpacked pds parameters
c     kgds         integer (200) unpacked gds parameters
c     kens         integer (200) unpacked ensemble pds parms
c     lskip        integer number of bytes to skip
c     lgrib        integer number of bytes to read
c     iret         integer return code
c                    0      all ok
c                    1      request not found
c
c remarks: subprogram can be called from a multiprocessing environment.
c   this subprogram is intended for private use by getgb routines only.
c
c subprograms called:
c   gbyte          unpack bytes
c   fi632          unpack pds
c   fi633          unpack gds
c   pdseup         unpack pds extension
c
c attributes:
c   language: fortran 77
c   machine:  cray, workstations
c
c$$$
      character cbuf(nlen*nnum)
      integer jpds(200),jgds(200),jens(200)
      integer kpds(200),kgds(200),kens(200)
      parameter(lpds=23,lgds=22,lens=5)     ! actual search ranges
      character cpds(400)*1,cgds(400)*1
      integer kptr(200)
      integer ipdsp(lpds),jpdsp(lpds)
      integer igdsp(lgds),jgdsp(lgds)
      integer iensp(lens),jensp(lens)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  compress request lists
      k=j
      lskip=0
      lgrib=0
      iret=1
c  compress pds request
      lpdsp=0
      do i=1,lpds
        if(jpds(i).ne.-1) then
          lpdsp=lpdsp+1
          ipdsp(lpdsp)=i
          jpdsp(lpdsp)=jpds(i)
        endif
      enddo
c  compress gds request
      lgdsp=0
      if(jpds(3).eq.255) then
        do i=1,lgds
          if(jgds(i).ne.-1) then
            lgdsp=lgdsp+1
            igdsp(lgdsp)=i
            jgdsp(lgdsp)=jgds(i)
          endif
        enddo
      endif
c  compress ens request
      lensp=0
      if(jpds(23).eq.2) then
        do i=1,lens
          if(jens(i).ne.-1) then
            lensp=lensp+1
            iensp(lensp)=i
            jensp(lensp)=jens(i)
          endif
        enddo
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  search for request
      dowhile(iret.ne.0.and.k.lt.nnum)
        k=k+1
        lt=0
c  search for pds request
        if(lpdsp.gt.0) then
          cpds=char(0)
          cpds(1:28)=cbuf((k-1)*nlen+26:(k-1)*nlen+53)
          nless=max(184-nlen,0)
          cpds(29:40-nless)=cbuf((k-1)*nlen+173:(k-1)*nlen+184-nless)
          kptr=0
          call gbyte(cbuf,kptr(3),(k-1)*nlen*8+25*8,3*8)
          kpds(18)=1
          call gbyte(cpds,kpds(4),7*8,8)
          call fi632(cpds,kptr,kpds,iret)
          do i=1,lpdsp
            ip=ipdsp(i)
            lt=lt+abs(jpds(ip)-kpds(ip))
          enddo
        endif
c  search for gds request
        if(lt.eq.0.and.lgdsp.gt.0) then
          cgds=char(0)
          cgds(1:42)=cbuf((k-1)*nlen+54:(k-1)*nlen+95)
          nless=max(320-nlen,0)
          cgds(43:178-nless)=cbuf((k-1)*nlen+185:(k-1)*nlen+320-nless)
          kptr=0
          call fi633(cgds,kptr,kgds,iret)
          do i=1,lgdsp
            ip=igdsp(i)
            lt=lt+abs(jgds(ip)-kgds(ip))
          enddo
        endif
c  search for ens request
        if(lt.eq.0.and.lensp.gt.0) then
          nless=max(172-nlen,0)
          cpds(41:100-nless)=cbuf((k-1)*nlen+113:(k-1)*nlen+172-nless)
          call pdseup(kens,kprob,xprob,kclust,kmembr,45,cpds)
          do i=1,lensp
            ip=iensp(i)
            lt=lt+abs(jens(ip)-kens(ip))
          enddo
        endif
c  return if request is found
        if(lt.eq.0) then
          call gbyte(cbuf,lskip,(k-1)*nlen*8,4*8)
          call gbyte(cbuf,lgrib,(k-1)*nlen*8+20*8,4*8)
          if(lpdsp.eq.0) then
            cpds=char(0)
            cpds(1:28)=cbuf((k-1)*nlen+26:(k-1)*nlen+53)
            nless=max(184-nlen,0)
            cpds(29:40-nless)=cbuf((k-1)*nlen+173:(k-1)*nlen+184-nless)
            kptr=0
            call gbyte(cbuf,kptr(3),(k-1)*nlen*8+25*8,3*8)
            kpds(18)=1
            call gbyte(cpds,kpds(4),7*8,8)
            call fi632(cpds,kptr,kpds,iret)
          endif
          if(lgdsp.eq.0) then
            cgds=char(0)
            cgds(1:42)=cbuf((k-1)*nlen+54:(k-1)*nlen+95)
            nless=max(320-nlen,0)
            cgds(43:178-nless)=cbuf((k-1)*nlen+185:(k-1)*nlen+320-nless)
            kptr=0
            call fi633(cgds,kptr,kgds,iret)
          endif
          if(kpds(23).eq.2.and.lensp.eq.0) then
            nless=max(172-nlen,0)
            cpds(41:100-nless)=cbuf((k-1)*nlen+113:(k-1)*nlen+172-nless)
            call pdseup(kens,kprob,xprob,kclust,kmembr,45,cpds)
          endif
          iret=0
        endif
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
