      subroutine dgetf2( m, n, a, lda, ipiv, info )
*
*  -- lapack routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     june 30, 1992
*
*     .. scalar arguments ..
      integer            info, lda, m, n
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      double precision   a( lda, * )
*     ..
*
*  purpose
*  =======
*
*  dgetf2 computes an lu factorization of a general m-by-n matrix a
*  using partial pivoting with row interchanges.
*
*  the factorization has the form
*     a = p * l * u
*  where p is a permutation matrix, l is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and u is upper
*  triangular (upper trapezoidal if m < n).
*
*  this is the right-looking level 2 blas version of the algorithm.
*
*  arguments
*  =========
*
*  m       (input) integer
*          the number of rows of the matrix a.  m >= 0.
*
*  n       (input) integer
*          the number of columns of the matrix a.  n >= 0.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the m by n matrix to be factored.
*          on exit, the factors l and u from the factorization
*          a = p*l*u; the unit diagonal elements of l are not stored.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,m).
*
*  ipiv    (output) integer array, dimension (min(m,n))
*          the pivot indices; for 1 <= i <= min(m,n), row i of the
*          matrix was interchanged with row ipiv(i).
*
*  info    (output) integer
*          = 0: successful exit
*          < 0: if info = -k, the k-th argument had an illegal value
*          > 0: if info = k, u(k,k) is exactly zero. the factorization
*               has been completed, but the factor u is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      integer            j, jp
*     ..
*     .. external functions ..
      integer            idamax
      external           idamax
*     ..
*     .. external subroutines ..
      external           dger, dscal, dswap, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          max, min
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      if( m.lt.0 ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, m ) ) then
         info = -4
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgetf2', -info )
         return
      end if
*
*     quick return if possible
*
      if( m.eq.0 .or. n.eq.0 )
     $   return
*
      do 10 j = 1, min( m, n )
*
*        find pivot and test for singularity.
*
         jp = j - 1 + idamax( m-j+1, a( j, j ), 1 )
         ipiv( j ) = jp
         if( a( jp, j ).ne.zero ) then
*
*           apply the interchange to columns 1:n.
*
            if( jp.ne.j )
     $         call dswap( n, a( j, 1 ), lda, a( jp, 1 ), lda )
*
*           compute elements j+1:m of j-th column.
*
            if( j.lt.m )
     $         call dscal( m-j, one / a( j, j ), a( j+1, j ), 1 )
*
         else if( info.eq.0 ) then
*
            info = j
         end if
*
         if( j.lt.min( m, n ) ) then
*
*           update trailing submatrix.
*
            call dger( m-j, n-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda,
     $                 a( j+1, j+1 ), lda )
         end if
   10 continue
      return
*
*     end of dgetf2
*
      end
