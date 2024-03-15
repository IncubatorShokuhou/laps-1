      subroutine dlaswp( n, a, lda, k1, k2, ipiv, incx )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      integer            incx, k1, k2, lda, n
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      double precision   a( lda, * )
*     ..
*
*  purpose
*  =======
*
*  dlaswp performs a series of row interchanges on the matrix a.
*  one row interchange is initiated for each of rows k1 through k2 of a.
*
*  arguments
*  =========
*
*  n       (input) integer
*          the number of columns of the matrix a.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the matrix of column dimension n to which the row
*          interchanges will be applied.
*          on exit, the permuted matrix.
*
*  lda     (input) integer
*          the leading dimension of the array a.
*
*  k1      (input) integer
*          the first element of ipiv for which a row interchange will
*          be done.
*
*  k2      (input) integer
*          the last element of ipiv for which a row interchange will
*          be done.
*
*  ipiv    (input) integer array, dimension (m*abs(incx))
*          the vector of pivot indices.  only the elements in positions
*          k1 through k2 of ipiv are accessed.
*          ipiv(k) = l implies rows k and l are to be interchanged.
*
*  incx    (input) integer
*          the increment between successive values of ipiv.  if ipiv
*          is negative, the pivots are applied in reverse order.
*
* =====================================================================
*
*     .. local scalars ..
      integer            i, ip, ix
*     ..
*     .. external subroutines ..
      external           dswap
*     ..
*     .. executable statements ..
*
*     interchange row i with row ipiv(i) for each of rows k1 through k2.
*
      if( incx.eq.0 )
     $   return
      if( incx.gt.0 ) then
         ix = k1
      else
         ix = 1 + ( 1-k2 )*incx
      end if
      if( incx.eq.1 ) then
         do 10 i = k1, k2
            ip = ipiv( i )
            if( ip.ne.i )
     $         call dswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
   10    continue
      else if( incx.gt.1 ) then
         do 20 i = k1, k2
            ip = ipiv( ix )
            if( ip.ne.i )
     $         call dswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
            ix = ix + incx
   20    continue
      else if( incx.lt.0 ) then
         do 30 i = k2, k1, -1
            ip = ipiv( ix )
            if( ip.ne.i )
     $         call dswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
            ix = ix + incx
   30    continue
      end if
*
      return
*
*     end of dlaswp
*
      end
