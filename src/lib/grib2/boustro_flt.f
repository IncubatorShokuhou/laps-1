      subroutine boustro_flt(a,nx,ny)
c
c        april   2000   lawrence   called by prepr
c        january 2001   glahn      comments; changed algorithm to
c                                  to use 2-d array rather than 1-d
c
c        purpose
c            scans and reorders an array of data boustrophedonically.
c            boustrophedonic ordering is a process in which
c            an array of data is scanned starting from the jy = 1
c            row and working up to the jy = ny row. the rows are
c            alternatively scanned from left to right and from
c            right to left. this scanning method is advantageous
c            when used to process a field of data that will be
c            packed with the grib2 complex packing method (with
c            or without spatial differences).  it generally
c            produces a data array which can be better broken
c            down into groups, a process which is core to the
c            complex packing method.
c
c        data set use
c           none
c
c        variables
c            a(ix,jy) = contains the data field to be
c                       boustrophedonically reordered
c                       (ix=1,nx) (jy=1,ny).  (input/output)
c                  nx = the number of columns in the array.  (input)
c                  ny = the number of rows in the array.  (input)
c
c        local variables
c                nxp1 = nx + 1.  this is done up front in case some
c                       compilers won't lift it out of the loop.
c                temp = saves the value of a( , ) currently being swapped
c                       so that it is not overwritten.
c
c         non system subroutines called
c            none
c
      dimension a(nx,ny)
c
      nxp1=nx+1
c
      do 125 jy=1,ny
c
         if(mod(jy,2).eq.0)then
c
            do 123 ix=1,nx/2
               temp=a(ix,jy)
               a(ix,jy)=a(nxp1-ix,jy)
               a(nxp1-ix,jy)=temp
 123        continue
c
         endif
c  
 125  continue
c
      return
      end
