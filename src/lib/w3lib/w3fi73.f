      subroutine w3fi73 (ibflag,ibmap,iblen,bms,lenbms,ier)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:  w3fi73        construct grib bit map section (bms)
c   prgmmr: farley           org: nmc421      date:92-11-16
c
c abstract: this subroutine constructs a grib bit map section.
c
c program history log:
c   92-07-01  m. farley      original author
c   94-02-14  cavanaugh      recoded
c   95-10-31  iredell        removed saves and prints
c
c usage:    call w3fi73 (ibflag, ibmap, iblen, bms, lenbms, ier)
c   input argument list:
c     ibflag      - 0, if bit map supplied by user
c                 - #, number of predefined center bit map
c     ibmap       - integer array containing user bit map
c     iblen       - length of bit map
c
c   output argument list:
c     bms       - completed grib bit map section
c     lenbms    - length of bit map section
c     ier       - 0 normal exit, 8 = ibmap values are all zero
c
c remarks: subprogram can be called from a multiprocessing environment.
c
c attributes:
c   language: ibm370 vs fortran 77, cray cft77 fortran
c   machine:  hds, cray c916-128, cray y-mp8/864, cray y-mp el2/256
c
c$$$
c
      integer       ibmap(*)
      integer       lenbms
      integer       iblen
      integer       ibflag
c
      character*1   bms   (*)
c
      ier   = 0
c
c
      iz  = 0
      do 20 i = 1, iblen
          if (ibmap(i).eq.0) iz  = iz + 1
   20 continue
      if (iz.eq.iblen) then
c
c                         at this point all bit map positions are zero
c
              ier = 8
              return
      end if
c
c                          bit map is a combination of ones and zeros
c                          or      bit map all ones
c
c                     construct bit map field of bit map section
c
      call sbytes (bms,ibmap,48,1,0,iblen)
c
      if (mod(iblen,16).ne.0) then
          nleft  = 16 - mod(iblen,16)
      else
          nleft  = 0
      end if
c
      num  = 6 + (iblen+nleft) / 8
c
c
c                          construct bms from collected data
c
c                          size into first three bytes
      call sbyte (bms,num,0,24)
c                          number of fill bits into byte 4
      call sbyte (bms,nleft,24,8)
c                          octet 5-6 to contain info from ibflag
      call sbyte (bms,ibflag,32,16)
c
c                          bit map may be all ones or a combination
c                          of ones and zeros
c
c                          actual bits of bit map placed all ready
c
c                          install fill positions if needed
      if (nleft.ne.0) then
          nleft  = 16 - nleft
c                          zero fill positions
          call sbyte (bms,0,iblen+48,nleft)
      end if
c
c     store num in lenbms  (length of bms section)
c
      lenbms = num
c     print *,'w3fi73 - bms len =',num,lenbms
c
      return
      end
