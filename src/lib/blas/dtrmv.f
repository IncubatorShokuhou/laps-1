      subroutine dtrmv ( uplo, trans, diag, n, a, lda, x, incx )
*     .. scalar arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
*     .. array arguments ..
      double precision   a( lda, * ), x( * )
*     ..
*
*  purpose
*  =======
*
*  dtrmv  performs one of the matrix-vector operations
*
*     x := a*x,   or   x := a'*x,
*
*  where x is an n element vector and  a is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  parameters
*  ==========
*
*  uplo   - character*1.
*           on entry, uplo specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              uplo = 'u' or 'u'   a is an upper triangular matrix.
*
*              uplo = 'l' or 'l'   a is a lower triangular matrix.
*
*           unchanged on exit.
*
*  trans  - character*1.
*           on entry, trans specifies the operation to be performed as
*           follows:
*
*              trans = 'n' or 'n'   x := a*x.
*
*              trans = 't' or 't'   x := a'*x.
*
*              trans = 'c' or 'c'   x := a'*x.
*
*           unchanged on exit.
*
*  diag   - character*1.
*           on entry, diag specifies whether or not a is unit
*           triangular as follows:
*
*              diag = 'u' or 'u'   a is assumed to be unit triangular.
*
*              diag = 'n' or 'n'   a is not assumed to be unit
*                                  triangular.
*
*           unchanged on exit.
*
*  n      - integer.
*           on entry, n specifies the order of the matrix a.
*           n must be at least zero.
*           unchanged on exit.
*
*  a      - double precision array of dimension ( lda, n ).
*           before entry with  uplo = 'u' or 'u', the leading n by n
*           upper triangular part of the array a must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           a is not referenced.
*           before entry with uplo = 'l' or 'l', the leading n by n
*           lower triangular part of the array a must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           a is not referenced.
*           note that when  diag = 'u' or 'u', the diagonal elements of
*           a are not referenced either, but are assumed to be unity.
*           unchanged on exit.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program. lda must be at least
*           max( 1, n ).
*           unchanged on exit.
*
*  x      - double precision array of dimension at least
*           ( 1 + ( n - 1 )*abs( incx ) ).
*           before entry, the incremented array x must contain the n
*           element vector x. on exit, x is overwritten with the
*           tranformed vector x.
*
*  incx   - integer.
*           on entry, incx specifies the increment for the elements of
*           x. incx must not be zero.
*           unchanged on exit.
*
*
*  level 2 blas routine.
*
*  -- written on 22-october-1986.
*     jack dongarra, argonne national lab.
*     jeremy du croz, nag central office.
*     sven hammarling, nag central office.
*     richard hanson, sandia national labs.
*
*
*     .. parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
*     .. local scalars ..
      double precision   temp
      integer            i, info, ix, j, jx, kx
      logical            nounit
*     .. external functions ..
      logical            lsame
      external           lsame
*     .. external subroutines ..
      external           xerbla
*     .. intrinsic functions ..
      intrinsic          max
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call xerbla( 'dtrmv ', info )
         return
      end if
*
*     quick return if possible.
*
      if( n.eq.0 )
     $   return
*
      nounit = lsame( diag, 'n' )
*
*     set up the start point in x if the increment is not unity. this
*     will be  ( n - 1 )*incx  too small for descending loops.
*
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
*
*     start the operations. in this version the elements of a are
*     accessed sequentially with one pass through a.
*
      if( lsame( trans, 'n' ) )then
*
*        form  x := a*x.
*
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*a( i, j )
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, i = 1, j - 1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx + incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*a( i, j )
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, i = n, j + 1, -1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx - incx
   80          continue
            end if
         end if
      else
*
*        form  x := a'*x.
*
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 100, j = n, 1, -1
                  temp = x( j )
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 90, i = j - 1, 1, -1
                     temp = temp + a( i, j )*x( i )
   90             continue
                  x( j ) = temp
  100          continue
            else
               jx = kx + ( n - 1 )*incx
               do 120, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 110, i = j - 1, 1, -1
                     ix   = ix   - incx
                     temp = temp + a( i, j )*x( ix )
  110             continue
                  x( jx ) = temp
                  jx      = jx   - incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = 1, n
                  temp = x( j )
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 130, i = j + 1, n
                     temp = temp + a( i, j )*x( i )
  130             continue
                  x( j ) = temp
  140          continue
            else
               jx = kx
               do 160, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 150, i = j + 1, n
                     ix   = ix   + incx
                     temp = temp + a( i, j )*x( ix )
  150             continue
                  x( jx ) = temp
                  jx      = jx   + incx
  160          continue
            end if
         end if
      end if
*
      return
*
*     end of dtrmv .
*
      end
