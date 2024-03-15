c-----------------------------------------------------------------------
      module bacio_module
c$$$  f90-module documentation block
c
c f90-module: bacio_module   byte-addressable i/o module
c   prgmmr: iredell          org: np23        date: 98-06-04
c
c abstract: module to share file descriptors
c   in the byte-addessable i/o package.
c
c program history log:
c   98-06-04  iredell
c
c attributes:
c   language: fortran 90
c
c$$$
      integer,external:: bacio
      integer,dimension(999),save:: fd=999*0
      integer,dimension(20),save:: baopts=0
      include 'baciof.h'
      end
c-----------------------------------------------------------------------
      subroutine baseto(nopt,vopt)
c$$$  subprogram documentation block
c
c subprogram: baseto         byte-addressable set options
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: set options for byte-addressable i/o.
c   all options default to 0.
c   option 1: blocked reading option
c             if the option value is 1, then the reading is blocked
c             into four 4096-byte buffers.  this may be efficient if
c             the reads will be requested in much smaller chunks.
c             otherwise, each call to baread initiates a physical read.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call baseto(nopt,vopt)
c   input arguments:
c     nopt         integer option number
c     vopt         integer option value
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      integer nopt,vopt
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(nopt.ge.1.and.nopt.le.20) baopts(nopt)=vopt
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine baopen(lu,cfn,iret)
c$$$  subprogram documentation block
c
c subprogram: baopen         byte-addressable open
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: open a byte-addressable file.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call baopen(lu,cfn,iret)
c   input arguments:
c     lu           integer unit to open
c     cfn          character filename to open
c                  (consisting of nonblank printable characters)
c   output arguments:
c     iret         integer return code
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      character cfn*(*)
      character(80) cmsg
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(lu.lt.001.or.lu.gt.999) then
        iret=6
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iret=bacio(bacio_openrw,ib,jb,1,nb,ka,fd(lu),cfn,a)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine baopenr(lu,cfn,iret)
c$$$  subprogram documentation block
c
c subprogram: baopenr        byte-addressable open
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: open a byte-addressable file for read only.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call baopenr(lu,cfn,iret)
c   input arguments:
c     lu           integer unit to open
c     cfn          character filename to open
c                  (consisting of nonblank printable characters)
c   output arguments:
c     iret         integer return code
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      character cfn*(*)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(lu.lt.001.or.lu.gt.999) then
        iret=6
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iret=bacio(bacio_openr,ib,jb,1,nb,ka,fd(lu),cfn,a)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine baopenw(lu,cfn,iret)
c$$$  subprogram documentation block
c
c subprogram: baopenw        byte-addressable open
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: open a byte-addressable file for write only.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call baopenw(lu,cfn,iret)
c   input arguments:
c     lu           integer unit to open
c     cfn          character filename to open
c                  (consisting of nonblank printable characters)
c   output arguments:
c     iret         integer return code
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      character cfn*(*)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(lu.lt.001.or.lu.gt.999) then
        iret=6
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iret=bacio(bacio_openw,ib,jb,1,nb,ka,fd(lu),cfn,a)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine baopenwt(lu,cfn,iret)
c$$$  subprogram documentation block
c
c subprogram: baopenwt       byte-addressable open
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: open a byte-addressable file for write only with truncation.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call baopenwt(lu,cfn,iret)
c   input arguments:
c     lu           integer unit to open
c     cfn          character filename to open
c                  (consisting of nonblank printable characters)
c   output arguments:
c     iret         integer return code
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      character cfn*(*)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(lu.lt.001.or.lu.gt.999) then
        iret=6
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iret=bacio(bacio_openwt,ib,jb,1,nb,ka,fd(lu),cfn,a)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine baopenwa(lu,cfn,iret)
c$$$  subprogram documentation block
c
c subprogram: baopenwa       byte-addressable open
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: open a byte-addressable file for write only with append.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call baopenwa(lu,cfn,iret)
c   input arguments:
c     lu           integer unit to open
c     cfn          character filename to open
c                  (consisting of nonblank printable characters)
c   output arguments:
c     iret         integer return code
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      character cfn*(*)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(lu.lt.001.or.lu.gt.999) then
        iret=6
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iret=bacio(bacio_openwa,ib,jb,1,nb,ka,fd(lu),cfn,a)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine baclose(lu,iret)
c$$$  subprogram documentation block
c
c subprogram: baclose        byte-addressable close
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: close a byte-addressable file.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call baclose(lu,iret)
c   input arguments:
c     lu           integer unit to close
c   output arguments:
c     iret         integer return code
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c remarks:  a baopen must have already been called.
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(lu.lt.001.or.lu.gt.999) then
        iret=6
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iret=bacio(bacio_close,ib,jb,1,nb,ka,fd(lu),cfn,a)
      if(iret.eq.0) fd(lu)=0
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine baread(lu,ib,nb,ka,a)
c$$$  subprogram documentation block
c
c subprogram: baread         byte-addressable read
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: read a given number of bytes from an unblocked file,
c   skipping a given number of bytes.
c   the physical i/o is blocked into four 4096-byte buffers
c   if the byte-addressable option 1 has been set to 1 by baseto.
c   this buffered reading is incompatible with no-seek reading.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call baread(lu,ib,nb,ka,a)
c   input arguments:
c     lu           integer unit to read
c     ib           integer number of bytes to skip
c                  (if ib<0, then the file is accessed with no seeking)
c     nb           integer number of bytes to read
c   output arguments:
c     ka           integer number of bytes actually read
c     a            character*1 (nb) data read
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c remarks:  a baopen must have already been called.
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      character a(nb)
      character cfn
      parameter(ny=4096,my=4)
      integer ns(my),nn(my)
      character y(ny,my)
      data lux/0/
      save jy,ns,nn,y,lux
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(fd(lu).le.0) then
        ka=0
        return
      endif
      if(ib.lt.0.and.baopts(1).eq.1) then
        ka=0
        return
      endif
      if(nb.le.0) then
        ka=0
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  unbuffered i/o
      if(baopts(1).ne.1) then
        if(ib.ge.0) then
          iret=bacio(bacio_read,ib,jb,1,nb,ka,fd(lu),cfn,a)
        else
          iret=bacio(bacio_read+bacio_noseek,0,jb,1,nb,ka,fd(lu),cfn,a)
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  buffered i/o
c  get data from previous call if possible
      else
        ka=0
        if(lux.ne.lu) then
          jy=0
          ns=0
          nn=0
        else
          do i=1,my
            iy=mod(jy+i-1,my)+1
            ky=ib+ka-ns(iy)
            if(ka.lt.nb.and.ky.ge.0.and.ky.lt.nn(iy)) then
              k=min(nb-ka,nn(iy)-ky)
              a(ka+1:ka+k)=y(ky+1:ky+k,iy)
              ka=ka+k
            endif
          enddo
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  set position and read buffer and get data
        if(ka.lt.nb) then
          lux=abs(lu)
          jy=mod(jy,my)+1
          ns(jy)=ib+ka
          iret=bacio(bacio_read,ns(jy),jb,1,ny,nn(jy),
     &               fd(lux),cfn,y(1,jy))
          if(nn(jy).gt.0) then
            k=min(nb-ka,nn(jy))
            a(ka+1:ka+k)=y(1:k,jy)
            ka=ka+k
          endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  continue to read buffer and get data
          dowhile(nn(jy).eq.ny.and.ka.lt.nb)
            jy=mod(jy,my)+1
            ns(jy)=ns(jy)+nn(jy)
            iret=bacio(bacio_read+bacio_noseek,ns(jy),jb,1,ny,nn(jy),
     &                 fd(lux),cfn,y(1,jy))
            if(nn(jy).gt.0) then
              k=min(nb-ka,nn(jy))
              a(ka+1:ka+k)=y(1:k,jy)
              ka=ka+k
            endif
          enddo
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine bawrite(lu,ib,nb,ka,a)
c$$$  subprogram documentation block
c
c subprogram: bawrite        byte-addressable write
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: write a given number of bytes to an unblocked file,
c   skipping a given number of bytes.
c
c program history log:
c   1998-06-04  iredell
c
c usage:    call bawrite(lu,ib,nb,ka,a)
c   input arguments:
c     lu           integer unit to write
c     ib           integer number of bytes to skip
c                  (if ib<0, then the file is accessed with no seeking)
c     nb           integer number of bytes to write
c     a            character*1 (nb) data to write
c   output arguments:
c     ka           integer number of bytes actually written
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c remarks:  a baopen must have already been called.
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      character a(nb)
      character cfn
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(fd(lu).le.0) then
        ka=0
        return
      endif
      if(nb.le.0) then
        ka=0
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(ib.ge.0) then
        iret=bacio(bacio_write,ib,jb,1,nb,ka,fd(lu),cfn,a)
      else
        iret=bacio(bacio_write+bacio_noseek,0,jb,1,nb,ka,fd(lu),cfn,a)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine wryte(lu,nb,a)
c$$$  subprogram documentation block
c
c subprogram: wryte          write data out by bytes
c   prgmmr: iredell          org: w/nmc23     date: 1998-06-04
c
c abstract: write a given number of bytes to an unblocked file.
c
c program history log:
c   92-10-31  iredell
c   95-10-31  iredell     workstation version
c   1998-06-04  iredell   bacio version
c
c usage:    call wryte(lu,nb,a)
c   input arguments:
c     lu           integer unit to which to write
c     nb           integer number of bytes to write
c     a            character*1 (nb) data to write
c
c modules used:
c   bacio_module   byte-addressable i/o fortran interface
c
c subprograms called:
c   bacio          byte-addressable i/o c package
c
c remarks:  a baopen must have already been called.
c
c attributes:
c   language: fortran 90
c
c$$$
      use bacio_module
      character a(nb)
      character cfn
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(fd(lu).le.0) then
        return
      endif
      if(nb.le.0) then
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iret=bacio(bacio_write+bacio_noseek,0,jb,1,nb,ka,fd(lu),cfn,a)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
