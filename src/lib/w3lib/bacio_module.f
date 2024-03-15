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

