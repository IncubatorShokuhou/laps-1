      subroutine r63w72(kpds,kgds,ipds,igds)
c$$$  subprogram documentation block
c
c subprogram:    r63w72      convert w3fi63 parms to w3fi72 parms
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31
c
c abstract: determines the integer pds and gds parameters
c           for the grib1 packing routine w3fi72 given the parameters
c           returned from the grib1 unpacking routine w3fi63.
c
c program history log:
c   91-10-31  mark iredell
c   96-05-03  mark iredell  corrected some level types and
c                           some data representation types
c   97-02-14  mark iredell  only altered ipds(26:27) for extended pds
c   98-06-01  chris caruso  y2k fix for year of century
c 2005-05-06  diane stokes  recognize level 236
c
c usage:    call r63w72(kpds,kgds,ipds,igds)
c
c   input argument list:
c     kpds     - integer (200) pds parameters from w3fi63
c     kgds     - integer (200) gds parameters from w3fi63
c
c   output argument list:
c     ipds     - integer (200) pds parameters for w3fi72
c     igds     - integer (200) gds parameters for w3fi72
c
c remarks: kgds and igds extend beyond their dimensions here
c          if pl parameters are present.
c
c attributes:
c   language: cray fortran
c
c$$$
      dimension kpds(200),kgds(200),ipds(200),igds(200)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  determine product definition section (pds) parameters
      if(kpds(23).ne.2) then
        ipds(1)=28                      ! length of pds
      else
        ipds(1)=45                      ! length of pds
      endif
      ipds(2)=kpds(19)                  ! parameter table version
      ipds(3)=kpds(1)                   ! originating center
      ipds(4)=kpds(2)                   ! generating model
      ipds(5)=kpds(3)                   ! grid definition
      ipds(6)=mod(kpds(4)/128,2)        ! gds flag
      ipds(7)=mod(kpds(4)/64,2)         ! bms flag
      ipds(8)=kpds(5)                   ! parameter indicator
      ipds(9)=kpds(6)                   ! level type
      if(kpds(6).eq.101.or.kpds(6).eq.104.or.kpds(6).eq.106.or.
     &   kpds(6).eq.108.or.kpds(6).eq.110.or.kpds(6).eq.112.or.
     &   kpds(6).eq.114.or.kpds(6).eq.116.or.kpds(6).eq.121.or.
     &   kpds(6).eq.128.or.kpds(6).eq.141.or.kpds(6).eq.236)  then
        ipds(10)=mod(kpds(7)/256,256)   ! level value 1
        ipds(11)=mod(kpds(7),256)       ! level value 2
      else
        ipds(10)=0                      ! level value 1
        ipds(11)=kpds(7)                ! level value 2
      endif
      ipds(12)=kpds(8)                  ! year of century
      ipds(13)=kpds(9)                  ! month
      ipds(14)=kpds(10)                 ! day
      ipds(15)=kpds(11)                 ! hour
      ipds(16)=kpds(12)                 ! minute
      ipds(17)=kpds(13)                 ! forecast time unit
      ipds(18)=kpds(14)                 ! time range 1
      ipds(19)=kpds(15)                 ! time range 2
      ipds(20)=kpds(16)                 ! time range indicator
      ipds(21)=kpds(17)                 ! number in average
      ipds(22)=kpds(20)                 ! number missing in average
      ipds(23)=kpds(21)                 ! century
      ipds(24)=kpds(23)                 ! subcenter
      ipds(25)=kpds(22)                 ! decimal scaling
      if(ipds(1).gt.28) then
        ipds(26)=0                      ! pds byte 29
        ipds(27)=0                      ! pds byte 30
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  determine grid definition section (gds) parameters
      igds(1)=kgds(19)                  ! number of vertical coordinates
      igds(2)=kgds(20)                  ! vertical coordinates
      igds(3)=kgds(1)                   ! data representation
      igds(4)=kgds(2)                   ! (unique to representation)
      igds(5)=kgds(3)                   ! (unique to representation)
      igds(6)=kgds(4)                   ! (unique to representation)
      igds(7)=kgds(5)                   ! (unique to representation)
      igds(8)=kgds(6)                   ! (unique to representation)
      igds(9)=kgds(7)                   ! (unique to representation)
      igds(10)=kgds(8)                  ! (unique to representation)
      igds(11)=kgds(9)                  ! (unique to representation)
      igds(12)=kgds(10)                 ! (unique to representation)
      igds(13)=kgds(11)                 ! (unique to representation)
      igds(14)=kgds(12)                 ! (unique to representation)
      igds(15)=kgds(13)                 ! (unique to representation)
      igds(16)=kgds(14)                 ! (unique to representation)
      igds(17)=kgds(15)                 ! (unique to representation)
      igds(18)=kgds(16)                 ! (unique to representation)
c  exceptions for latlon or gaussian
      if(kgds(1).eq.0.or.kgds(1).eq.4) then
        igds(11)=kgds(10)
        igds(12)=kgds(9)
c  exceptions for mercator
      elseif(kgds(1).eq.1) then
        igds(11)=kgds(13)
        igds(12)=kgds(12)
        igds(13)=kgds(9)
        igds(14)=kgds(11)
c  exceptions for lambert conformal
      elseif(kgds(1).eq.3) then
        igds(15)=kgds(12)
        igds(16)=kgds(13)
        igds(17)=kgds(14)
        igds(18)=kgds(15)
      endif
c  extension for pl parameters
      if(kgds(1).eq.0.and.kgds(19).eq.0.and.kgds(20).ne.255) then
        do j=1,kgds(3)
          igds(18+j)=kgds(21+j)
        enddo
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
