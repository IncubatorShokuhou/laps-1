c-----------------------------------------------------------------------
      function lengds(kgds)
c$$$  subprogram documentation block
c
c subprogram:    lengds      return the length of a grid
c   prgmmr: iredell          org: w/nmc23     date: 96-07-19
c
c abstract: given a grid description section (in w3fi63 format),
c   return its size in terms of number of data points.
c
c program history log:
c   96-07-19  iredell
c
c usage:    call lengds(kgds)
c   input arguments:
c     kgds         integer (200) gds parameters in w3fi63 format
c   output arguments:
c     lengds       integer size of grid
c
c attributes:
c   language: fortran
c
c$$$
      integer kgds(200)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  special case of staggered eta
      if(kgds(1).eq.201) then
        lengds=kgds(7)*kgds(8)-kgds(8)/2
c  special case of filled eta
      elseif(kgds(1).eq.202) then
        lengds=kgds(7)*kgds(8)
c  special case of thinned wafs
      elseif(kgds(19).eq.0.and.kgds(20).ne.255) then
        lengds=kgds(21)
c  general case
      else
        lengds=kgds(2)*kgds(3)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
