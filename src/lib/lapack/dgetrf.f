      subroutine dgetrf( m, n, a, lda, ipiv, info )
*
*  -- lapack routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     march 31, 1993
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
*  dgetrf computes an lu factorization of a general m-by-n matrix a
*  using partial pivoting with row interchanges.
*
*  the factorization has the form
*     a = p * l * u
*  where p is a permutation matrix, l is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and u is upper
*  triangular (upper trapezoidal if m < n).
*
*  this is the right-looking level 3 blas version of the algorithm.
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
*          on entry, the m-by-n matrix to be factored.
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
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*          > 0:  if info = i, u(i,i) is exactly zero. the factorization
*                has been completed, but the factor u is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. parameters ..
      double precision   one
      parameter          ( one = 1.0d+0 )
*     ..
*     .. local scalars ..
      integer            i, iinfo, j, jb, nb
*     ..
*     .. external subroutines ..
      external           dgemm, dgetf2, dlaswp, dtrsm, xerbla
*     ..
*     .. external functions ..
      integer            ilaenv
      external           ilaenv
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
         call xerbla( 'dgetrf', -info )
         return
      end if
*
*     quick return if possible
*
      if( m.eq.0 .or. n.eq.0 )
     $   return
*
*     determine the block size for this environment.
*
      nb = ilaenv( 1, 'dgetrf', ' ', m, n, -1, -1 )
      if( nb.le.1 .or. nb.ge.min( m, n ) ) then
*
*        use unblocked code.
*
         call dgetf2( m, n, a, lda, ipiv, info )
      else
*
*        use blocked code.
*
         do 20 j = 1, min( m, n ), nb
            jb = min( min( m, n )-j+1, nb )
*
*           factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            call dgetf2( m-j+1, jb, a( j, j ), lda, ipiv( j ), iinfo )
*
*           adjust info and the pivot indices.
*
            if( info.eq.0 .and. iinfo.gt.0 )
     $         info = iinfo + j - 1
            do 10 i = j, min( m, j+jb-1 )
               ipiv( i ) = j - 1 + ipiv( i )
   10       continue
*
*           apply interchanges to columns 1:j-1.
*
            call dlaswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )
*
            if( j+jb.le.n ) then
*
*              apply interchanges to columns j+jb:n.
*
               call dlaswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1,
     $                      ipiv, 1 )
*
*              compute block row of u.
*
               call dtrsm( 'left', 'lower', 'no transpose', 'unit', jb,
     $                     n-j-jb+1, one, a( j, j ), lda, a( j, j+jb ),
     $                     lda )
               if( j+jb.le.m ) then
*
*                 update trailing submatrix.
*
                  call dgemm( 'no transpose', 'no transpose', m-j-jb+1,
     $                        n-j-jb+1, jb, -one, a( j+jb, j ), lda,
     $                        a( j, j+jb ), lda, one, a( j+jb, j+jb ),
     $                        lda )
               end if
            end if
   20    continue
      end if
      return
*
*     end of dgetrf
*
      end
