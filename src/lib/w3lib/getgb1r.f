c-----------------------------------------------------------------------
      subroutine getgb1r(lugb,lskip,lgrib,kf,kpds,kgds,kens,lb,f,nbitss
     +      ,iret)
c$$$  subprogram documentation block
c
c subprogram: getgb1r        reads and unpacks a grib message
c   prgmmr: iredell          org: w/nmc23     date: 95-10-31
c
c abstract: read and unpack a grib message.
c
c program history log:
c   95-10-31  iredell
c   04-07-22  chuang  add packing bit number nbitss in the argument 
c            list because eta grib files need it to repack grib file  
c usage:    call getgb1r(lugb,lskip,lgrib,kf,kpds,kgds,kens,lb,f,iret)
c   input arguments:
c     lugb         integer unit of the unblocked grib data file
c     lskip        integer number of bytes to skip
c     lgrib        integer number of bytes to read
c   output arguments:
c     kf           integer number of data points unpacked
c     kpds         integer (200) unpacked pds parameters
c     kgds         integer (200) unpacked gds parameters
c     kens         integer (200) unpacked ensemble pds parms
c     lb           logical*1 (kf) unpacked bitmap if present
c     f            real (kf) unpacked data
c     iret         integer return code
c                    0      all ok
c                    97     error reading grib file
c                    other  w3fi63 grib unpacker return code
c
c subprograms called:
c   baread         byte-addressable read
c   w3fi63         unpack grib
c   pdseup         unpack pds extension
c
c remarks: there is no protection against unpacking too much data.
c   subprogram can be called from a multiprocessing environment.
c   do not engage the same logical unit from more than one processor.
c   this subprogram is intended for private use by getgb routines only.
c
c attributes:
c   language: fortran 77
c   machine:  cray, workstations
c
c$$$
      integer kpds(200),kgds(200),kens(200)
      logical*1 lb(*)
      real f(*)
      integer kptr(200)
      character grib(lgrib)*1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  read grib record
      call baread(lugb,lskip,lgrib,lread,grib)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  unpack grib record
      if(lread.eq.lgrib) then
        call w3fi63(grib,kpds,kgds,lb,f,kptr,iret)
        if(iret.eq.0.and.kpds(23).eq.2) then
          call pdseup(kens,kprob,xprob,kclust,kmembr,45,grib(9))
        endif
      else
        iret=97
      endif
      nbitss=kptr(20) 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  return number of points
      if(iret.eq.0) then
        kf=kptr(10)
      else
        kf=0
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
