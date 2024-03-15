      subroutine dgetri( n, a, lda, ipiv, work, lwork, info )
*
*  -- lapack routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     september 30, 1994
*
*     .. scalar arguments ..
      integer            info, lda, lwork, n
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      double precision   a( lda, * ), work( lwork )
*     ..
*
*  purpose
*  =======
*
*  dgetri computes the inverse of a matrix using the lu factorization
*  computed by dgetrf.
*
*  this method inverts u and then computes inv(a) by solving the system
*  inv(a)*l = inv(u) for inv(a).
*
*  arguments
*  =========
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the factors l and u from the factorization
*          a = p*l*u as computed by dgetrf.
*          on exit, if info = 0, the inverse of the original matrix a.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  ipiv    (input) integer array, dimension (n)
*          the pivot indices from dgetrf; for 1<=i<=n, row i of the
*          matrix was interchanged with row ipiv(i).
*
*  work    (workspace/output) double precision array, dimension (lwork)
*          on exit, if info=0, then work(1) returns the optimal lwork.
*
*  lwork   (input) integer
*          the dimension of the array work.  lwork >= max(1,n).
*          for optimal performance lwork >= n*nb, where nb is
*          the optimal blocksize returned by ilaenv.
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*          > 0:  if info = i, u(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. parameters ..
      double precision   zero, one
      parameter          ( zero = 0.0d+0, one = 1.0d+0 )
*     ..
*     .. local scalars ..
      integer            i, iws, j, jb, jj, jp, ldwork, nb, nbmin, nn
*     ..
*     .. external functions ..
      integer            ilaenv
      external           ilaenv
*     ..
*     .. external subroutines ..
      external           dgemm, dgemv, dswap, dtrsm, dtrtri, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          max, min
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      work( 1 ) = max( n, 1 )
      if( n.lt.0 ) then
         info = -1
      else if( lda.lt.max( 1, n ) ) then
         info = -3
      else if( lwork.lt.max( 1, n ) ) then
         info = -6
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgetri', -info )
         return
      end if
*
*     quick return if possible
*
      if( n.eq.0 )
     $   return
*
*     form inv(u).  if info > 0 from dtrtri, then u is singular,
*     and the inverse is not computed.
*
      call dtrtri( 'upper', 'non-unit', n, a, lda, info )
      if( info.gt.0 )
     $   return
*
*     determine the block size for this environment.
*
      nb = ilaenv( 1, 'dgetri', ' ', n, -1, -1, -1 )
      nbmin = 2
      ldwork = n
      if( nb.gt.1 .and. nb.lt.n ) then
         iws = max( ldwork*nb, 1 )
         if( lwork.lt.iws ) then
            nb = lwork / ldwork
            nbmin = max( 2, ilaenv( 2, 'dgetri', ' ', n, -1, -1, -1 ) )
         end if
      else
         iws = n
      end if
*
*     solve the equation inv(a)*l = inv(u) for inv(a).
*
      if( nb.lt.nbmin .or. nb.ge.n ) then
*
*        use unblocked code.
*
         do 20 j = n, 1, -1
*
*           copy current column of l to work and replace with zeros.
*
            do 10 i = j + 1, n
               work( i ) = a( i, j )
               a( i, j ) = zero
   10       continue
*
*           compute current column of inv(a).
*
            if( j.lt.n )
     $         call dgemv( 'no transpose', n, n-j, -one, a( 1, j+1 ),
     $                     lda, work( j+1 ), 1, one, a( 1, j ), 1 )
   20    continue
      else
*
*        use blocked code.
*
         nn = ( ( n-1 ) / nb )*nb + 1
         do 50 j = nn, 1, -nb
            jb = min( nb, n-j+1 )
*
*           copy current block column of l to work and replace with
*           zeros.
*
            do 40 jj = j, j + jb - 1
               do 30 i = jj + 1, n
                  work( i+( jj-j )*ldwork ) = a( i, jj )
                  a( i, jj ) = zero
   30          continue
   40       continue
*
*           compute current block column of inv(a).
*
            if( j+jb.le.n )
     $         call dgemm( 'no transpose', 'no transpose', n, jb,
     $                     n-j-jb+1, -one, a( 1, j+jb ), lda,
     $                     work( j+jb ), ldwork, one, a( 1, j ), lda )
            call dtrsm( 'right', 'lower', 'no transpose', 'unit', n, jb,
     $                  one, work( j ), ldwork, a( 1, j ), lda )
   50    continue
      end if
*
*     apply column interchanges.
*
      do 60 j = n - 1, 1, -1
         jp = ipiv( j )
         if( jp.ne.j )
     $      call dswap( n, a( 1, j ), 1, a( 1, jp ), 1 )
   60 continue
*
      work( 1 ) = iws
      return
*
*     end of dgetri
*
      end
