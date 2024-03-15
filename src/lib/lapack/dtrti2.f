      subroutine dtrti2( uplo, diag, n, a, lda, info )
*
*  -- lapack routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     february 29, 1992
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
*  dtrti2 computes the inverse of a real upper or lower triangular
*  matrix.
*
*  this is the level 2 blas version of the algorithm.
*
*  arguments
*  =========
*
*  uplo    (input) character*1
*          specifies whether the matrix a is upper or lower triangular.
*          = 'u':  upper triangular
*          = 'l':  lower triangular
*
*  diag    (input) character*1
*          specifies whether or not the matrix a is unit triangular.
*          = 'n':  non-unit triangular
*          = 'u':  unit triangular
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the triangular matrix a.  if uplo = 'u', the
*          leading n by n upper triangular part of the array a contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of a is not referenced.  if uplo = 'l', the
*          leading n by n lower triangular part of the array a contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of a is not referenced.  if diag = 'u', the
*          diagonal elements of a are also not referenced and are
*          assumed to be 1.
*
*          on exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  info    (output) integer
*          = 0: successful exit
*          < 0: if info = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. parameters ..
      double precision   one
      parameter          ( one = 1.0d+0 )
*     ..
*     .. local scalars ..
      logical            nounit, upper
      integer            j
      double precision   ajj
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. external subroutines ..
      external           dscal, dtrmv, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          max
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
         call xerbla( 'dtrti2', -info )
         return
      end if
*
      if( upper ) then
*
*        compute inverse of upper triangular matrix.
*
         do 10 j = 1, n
            if( nounit ) then
               a( j, j ) = one / a( j, j )
               ajj = -a( j, j )
            else
               ajj = -one
            end if
*
*           compute elements 1:j-1 of j-th column.
*
            call dtrmv( 'upper', 'no transpose', diag, j-1, a, lda,
     $                  a( 1, j ), 1 )
            call dscal( j-1, ajj, a( 1, j ), 1 )
   10    continue
      else
*
*        compute inverse of lower triangular matrix.
*
         do 20 j = n, 1, -1
            if( nounit ) then
               a( j, j ) = one / a( j, j )
               ajj = -a( j, j )
            else
               ajj = -one
            end if
            if( j.lt.n ) then
*
*              compute elements j+1:n of j-th column.
*
               call dtrmv( 'lower', 'no transpose', diag, n-j,
     $                     a( j+1, j+1 ), lda, a( j+1, j ), 1 )
               call dscal( n-j, ajj, a( j+1, j ), 1 )
            end if
   20    continue
      end if
*
      return
*
*     end of dtrti2
*
      end
