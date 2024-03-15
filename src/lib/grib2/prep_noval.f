      subroutine prep_noval(ib,nxy,ibitmap,ipkopt,ier)
c
c        may      2000   lawrence  original coding
c        january  2001   glahn     comments
c        february 2002   glahn     ipkopt set = 2 when a bit must
c                                  be pakced 
c        march    2002   glahn     comments
c
c        purpose
c            to handle the situation where there are no values
c            in the grid. three possibilities are checked by this
c            routine:
c
c            1) there are no values in the grid, but there is a
c               bit-map and it indicates that there should be. this
c               is treated as an error.
c            2) there are no values in the grid, and the bit-map
c               indicates that there are no values in the grid.
c               the user is notified and the field is
c               packed using the simple method.
c            3) there are no values in the grid, and there is no
c               accompanying bit-map. this is not considered to
c               be a reasonable situation. it is treated as an
c               error.
c
c        data set use
c           kfildo - unit number for output (print) file. (output)
c
c        variables
c               ib(k) = the bit map when one is used.  (input/output)
c                       it can be input or in can be calculated if
c                       the simple method is used (k=1,nxy).
c                       complex and spatial differencing do not
c                       use a bit map, but will accept one and insert
c                       the missing values.
c                 nxy = dimension of ib( ).  (input)
c             ibitmap = 1 when there is a bitmap in ib( ).
c                       0 otherwise.
c                       (input)
c              ipkopt = packing indicator:
c                       0 = error, don't pack
c                       1 = pack ia( ), simple
c                       2 = pack ia( ) and ib( ), simple
c                       3 = pack complex or spatial differencing
c                       4 = pack complex.
c                       (output)
c                 ier = return status codes. values can come from
c                       subroutines; otherwise:
c                       0 = good return.
c                     902 = there are no "good" values in the grid
c                           and the bit-map indicates that.  not treated
c                           as fatal in prepr.
c                     903 = there are no values in the grid and the
c                           bit-map indicates that there should be.
c                           treated as fatal in prepr.
c                     904 = there are no values in the grid, and there
c                           is no bit-map.  treated as fatal in prepr.
c                    
c        local variables
c                     k = a loop index variable.
c
c        non system subroutines called
c           none
c
      dimension ib(nxy)
c
c        this section when there are no values in
c        the data field. that is, nval = 0.
      if(ibitmap.eq.1)then
c
         do 10 k=1,nxy
c
            if(ib(k).eq.1)then
c                 there is a bit map.
c                 there are no values in the grid
c                 and the bit map indicates there
c                 should be. don't pack.
               ier=903
               ipkopt=0
               go to 900
            endif
c
 10      continue
c
c           there is a bit map.
c           there are no values in the grid
c           and the bit map indicates no values.
c           notify user and pack.  even though all values
c           in the bit map are zero, the wmo standards
c           evidently do not allow for not packing the 
c           bit map.
         ier=902
         ipkopt=2
c
      else
c           there are no values in grid and no
c           bit map.  not a reasonable situation.
c           don't pack.
         ier=904
         ipkopt=0
         go to 900
      endif
c
 900  return
      end
