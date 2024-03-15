c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    pdseup.f    unpacks grib pds extension 41- for ensemble
c   prgmmr: richard wobus    org: w/np20     date: 98-09-28
c
c abstract: unpacks grib pds extension starting on byte 41 for ensemble
c	forecast products. for format of pds extension, see nmc office note 38
c
c program history log:
c   95-03-14  zoltan toth and mark iredell
c   95-10-31  iredell     removed saves and prints
c   98-09-28  wobus       corrected member extraction
c
c usage:    call pdsens.f(kens,kprob,xprob,kclust,kmembr,ilast,msga)
c   input argument list:
c     ilast    - last byte to be unpacked (if greater/equal to first byt
c                in any of four sections below, whole section is packed.
c     msga     - full pds section, including new ensemble extension
c
c   output argument list:      (including work arrays)
c     kens(5)  - bytes 41-45 (general section, always present.)
c     kprob(2) - bytes 46-47 (probability section, present only if neede
c     xprob(2) - bytes 48-51&52-55 (probability section, if needed.)
c     kclust(16)-bytes 61-76 (clustering section, if needed.)
c     kmembr(80)-bytes 77-86 (cluster membership section, if needed.)
c
c remarks: use pdsens.f for packing pds ensemble extension.
c   subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: cf77 fortran
c   machine:  cray, workstations
c
c$$$
c
	  subroutine pdseup(kens,kprob,xprob,kclust,kmembr,ilast,msga)
	  integer kens(5),kprob(2),kclust(16),kmembr(80)
	  dimension xprob(2)
	  integer kref
	  character*1 msga(100)
	  real refnce
	  character*1 ckref(8)
	  equivalence   (ckref(1),kref,refnce)
c	checking total number of bytes in pds (ibytes)
	  call gbyte(msga, ibytes, 0,24)
	  if(ilast.gt.ibytes) then
c	  ilast=ibytes
	  go to 333
	  endif
	  if(ilast.lt.41) then
	  go to 333
	  endif
c	unpacking first section (general information)
      call gbytes(msga,kens,40*8,8,0,5)
c	unpacking 2nd section (probability section)
      if(ilast.ge.46) then
      call gbytes(msga,kprob,45*8,8,0,2)
c
c
      call gbyte (msga,kref,47*8,32)
      call w3fi01(lw)
      if (lw.eq.4) then
        call gbyte (ckref,jsgn,0,1)
        call gbyte (ckref,jexp,1,7)
        call gbyte (ckref,ifr,8,24)
      else
        call gbyte (ckref,jsgn,32,1)
        call gbyte (ckref,jexp,33,7)
        call gbyte (ckref,ifr,40,24)
      endif
      if (ifr.eq.0) then
          refnce  = 0.0
      else if (jexp.eq.0.and.ifr.eq.0) then
          refnce  = 0.0
      else
          refnce  = float(ifr) * 16.0 ** (jexp - 64 - 6)
          if (jsgn.ne.0) refnce = - refnce
      end if
	  xprob(1)=refnce
c
      call gbyte (msga,kref,51*8,32)
      call w3fi01(lw)
      if (lw.eq.4) then
        call gbyte (ckref,jsgn,0,1)
        call gbyte (ckref,jexp,1,7)
        call gbyte (ckref,ifr,8,24)
      else
        call gbyte (ckref,jsgn,32,1)
        call gbyte (ckref,jexp,33,7)
        call gbyte (ckref,ifr,40,24)
      endif
      if (ifr.eq.0) then
          refnce  = 0.0
      else if (jexp.eq.0.and.ifr.eq.0) then
          refnce  = 0.0
      else
          refnce  = float(ifr) * 16.0 ** (jexp - 64 - 6)
          if (jsgn.ne.0) refnce = - refnce
      end if
	  xprob(2)=refnce
	  endif
c
c	unpacking 3rd section (clustering information)
      if(ilast.ge.61) call gbytes(msga,kclust,60*8,8,0,16)
c	unpacking 4th section (clustermembership information)
      if(ilast.ge.77) call gbytes(msga,kmembr,76*8,1,0,80)
c
 333  continue
	  return
	  end
