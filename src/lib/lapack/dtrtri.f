      subroutine dtrtri( uplo, diag, n, a, lda, info )
*
*  -- lapack routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     march 31, 1993
*
*     .. scalar arguments ..
      character          diag, uplo
      integer            info, lda, n
*     ..
*     .. array arguments ..
      double precision   a( lda, * )
*     ..
*
*  purpose
*  =======
*
*  dtrtri computes the inverse of a real upper or lower triangular
*  matrix a.
*
*  this is the level 3 blas version of the algorithm.
*
*  arguments
*  =========
*
*  uplo    (input) character*1
*          = 'u':  a is upper triangular;
*          = 'l':  a is lower triangular.
*
*  diag    (input) character*1
*          = 'n':  a is non-unit triangular;
*          = 'u':  a is unit triangular.
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the triangular matrix a.  if uplo = 'u', the
*          leading n-by-n upper triangular part of the array a contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of a is not referenced.  if uplo = 'l', the
*          leading n-by-n lower triangular part of the array a contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of a is not referenced.  if diag = 'u', the
*          diagonal elements of a are also not referenced and are
*          assumed to be 1.
*          on exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  info    (output) integer
*          = 0: successful exit
*          < 0: if info = -i, the i-th argument had an illegal value
*          > 0: if info = i, a(i,i) is exactly zero.  the triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      logical            nounit, upper
      integer            j, jb, nb, nn
*     ..
*     .. external functions ..
      logical            lsame
      integer            ilaenv
      external           lsame, ilaenv
*     ..
*     .. external subroutines ..
      external           dtrmm, dtrsm, dtrti2, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          max, min
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      upper = lsame( uplo, 'u' )
      nounit = lsame( diag, 'n' )
      if( .not.upper .and. .not.lsame( uplo, 'l' ) ) then
         info = -1
      else if( .not.nounit .and. .not.lsame( diag, 'u' ) ) then
         info = -2
      else if( n.lt.0 ) then
         info = -3
      else if( lda.lt.max( 1, n ) ) then
         info = -5
      end if
      if( info.ne.0 ) then
         call xerbla( 'dtrtri', -info )
         return
      end if
*
*     quick return if possible
*
      if( n.eq.0 )
     $   return
*
*     check for singularity if non-unit.
*
      if( nounit ) then
         do 10 info = 1, n
            if( a( info, info ).eq.zero )
     $         return
   10    continue
         info = 0
      end if
*
*     determine the block size for this environment.
*
      nb = ilaenv( 1, 'dtrtri', uplo // diag, n, -1, -1, -1 )
      if( nb.le.1 .or. nb.ge.n ) then
*
*        use unblocked code
*
         call dtrti2( uplo, diag, n, a, lda, info )
      else
*
*        use blocked code
*
         if( upper ) then
*
*           compute inverse of upper triangular matrix
*
            do 20 j = 1, n, nb
               jb = min( nb, n-j+1 )
*
*              compute rows 1:j-1 of current block column
*
               call dtrmm( 'left', 'upper', 'no transpose', diag, j-1,
     $                     jb, one, a, lda, a( 1, j ), lda )
               call dtrsm( 'right', 'upper', 'no transpose', diag, j-1,
     $                     jb, -one, a( j, j ), lda, a( 1, j ), lda )
*
*              compute inverse of current diagonal block
*
               call dtrti2( 'upper', diag, jb, a( j, j ), lda, info )
   20       continue
         else
*
*           compute inverse of lower triangular matrix
*
            nn = ( ( n-1 ) / nb )*nb + 1
            do 30 j = nn, 1, -nb
               jb = min( nb, n-j+1 )
               if( j+jb.le.n ) then
*
*                 compute rows j+jb:n of current block column
*
                  call dtrmm( 'left', 'lower', 'no transpose', diag,
     $                        n-j-jb+1, jb, one, a( j+jb, j+jb ), lda,
     $                        a( j+jb, j ), lda )
                  call dtrsm( 'right', 'lower', 'no transpose', diag,
     $                        n-j-jb+1, jb, -one, a( j, j ), lda,
     $                        a( j+jb, j ), lda )
               end if
*
*              compute inverse of current diagonal block
*
               call dtrti2( 'lower', diag, jb, a( j, j ), lda, info )
   30       continue
         end if
      end if
*
      return
*
*     end of dtrtri
*
      end
