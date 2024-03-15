      subroutine dtrmm ( side, uplo, transa, diag, m, n, alpha, a, lda,
     $                   b, ldb )
*     .. scalar arguments ..
      character*1        side, uplo, transa, diag
      integer            m, n, lda, ldb
      double precision   alpha
*     .. array arguments ..
      double precision   a( lda, * ), b( ldb, * )
*     ..
*
*  purpose
*  =======
*
*  dtrmm  performs one of the matrix-matrix operations
*
*     b := alpha*op( a )*b,   or   b := alpha*b*op( a ),
*
*  where  alpha  is a scalar,  b  is an m by n matrix,  a  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of
*
*     op( a ) = a   or   op( a ) = a'.
*
*  parameters
*  ==========
*
*  side   - character*1.
*           on entry,  side specifies whether  op( a ) multiplies b from
*           the left or right as follows:
*
*              side = 'l' or 'l'   b := alpha*op( a )*b.
*
*              side = 'r' or 'r'   b := alpha*b*op( a ).
*
*           unchanged on exit.
*
*  uplo   - character*1.
*           on entry, uplo specifies whether the matrix a is an upper or
*           lower triangular matrix as follows:
*
*              uplo = 'u' or 'u'   a is an upper triangular matrix.
*
*              uplo = 'l' or 'l'   a is a lower triangular matrix.
*
*           unchanged on exit.
*
*  transa - character*1.
*           on entry, transa specifies the form of op( a ) to be used in
*           the matrix multiplication as follows:
*
*              transa = 'n' or 'n'   op( a ) = a.
*
*              transa = 't' or 't'   op( a ) = a'.
*
*              transa = 'c' or 'c'   op( a ) = a'.
*
*           unchanged on exit.
*
*  diag   - character*1.
*           on entry, diag specifies whether or not a is unit triangular
*           as follows:
*
*              diag = 'u' or 'u'   a is assumed to be unit triangular.
*
*              diag = 'n' or 'n'   a is not assumed to be unit
*                                  triangular.
*
*           unchanged on exit.
*
*  m      - integer.
*           on entry, m specifies the number of rows of b. m must be at
*           least zero.
*           unchanged on exit.
*
*  n      - integer.
*           on entry, n specifies the number of columns of b.  n must be
*           at least zero.
*           unchanged on exit.
*
*  alpha  - double precision.
*           on entry,  alpha specifies the scalar  alpha. when  alpha is
*           zero then  a is not referenced and  b need not be set before
*           entry.
*           unchanged on exit.
*
*  a      - double precision array of dimension ( lda, k ), where k is m
*           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
*           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
*           upper triangular part of the array  a must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           a is not referenced.
*           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
*           lower triangular part of the array  a must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           a is not referenced.
*           note that when  diag = 'u' or 'u',  the diagonal elements of
*           a  are not referenced either,  but are assumed to be  unity.
*           unchanged on exit.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program.  when  side = 'l' or 'l'  then
*           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
*           then lda must be at least max( 1, n ).
*           unchanged on exit.
*
*  b      - double precision array of dimension ( ldb, n ).
*           before entry,  the leading  m by n part of the array  b must
*           contain the matrix  b,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  ldb    - integer.
*           on entry, ldb specifies the first dimension of b as declared
*           in  the  calling  (sub)  program.   ldb  must  be  at  least
*           max( 1, m ).
*           unchanged on exit.
*
*
*  level 3 blas routine.
*
*  -- written on 8-february-1989.
*     jack dongarra, argonne national laboratory.
*     iain duff, aere harwell.
*     jeremy du croz, numerical algorithms group ltd.
*     sven hammarling, numerical algorithms group ltd.
*
*
*     .. external functions ..
      logical            lsame
      external           lsame
*     .. external subroutines ..
      external           xerbla
*     .. intrinsic functions ..
      intrinsic          max
*     .. local scalars ..
      logical            lside, nounit, upper
      integer            i, info, j, k, nrowa
      double precision   temp
*     .. parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      lside  = lsame( side  , 'l' )
      if( lside )then
         nrowa = m
      else
         nrowa = n
      end if
      nounit = lsame( diag  , 'n' )
      upper  = lsame( uplo  , 'u' )
*
      info   = 0
      if(      ( .not.lside                ).and.
     $         ( .not.lsame( side  , 'r' ) )      )then
         info = 1
      else if( ( .not.upper                ).and.
     $         ( .not.lsame( uplo  , 'l' ) )      )then
         info = 2
      else if( ( .not.lsame( transa, 'n' ) ).and.
     $         ( .not.lsame( transa, 't' ) ).and.
     $         ( .not.lsame( transa, 'c' ) )      )then
         info = 3
      else if( ( .not.lsame( diag  , 'u' ) ).and.
     $         ( .not.lsame( diag  , 'n' ) )      )then
         info = 4
      else if( m  .lt.0               )then
         info = 5
      else if( n  .lt.0               )then
         info = 6
      else if( lda.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldb.lt.max( 1, m     ) )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'dtrmm ', info )
         return
      end if
*
*     quick return if possible.
*
      if( n.eq.0 )
     $   return
*
*     and when  alpha.eq.zero.
*
      if( alpha.eq.zero )then
         do 20, j = 1, n
            do 10, i = 1, m
               b( i, j ) = zero
   10       continue
   20    continue
         return
      end if
*
*     start the operations.
*
      if( lside )then
         if( lsame( transa, 'n' ) )then
*
*           form  b := alpha*a*b.
*
            if( upper )then
               do 50, j = 1, n
                  do 40, k = 1, m
                     if( b( k, j ).ne.zero )then
                        temp = alpha*b( k, j )
                        do 30, i = 1, k - 1
                           b( i, j ) = b( i, j ) + temp*a( i, k )
   30                   continue
                        if( nounit )
     $                     temp = temp*a( k, k )
                        b( k, j ) = temp
                     end if
   40             continue
   50          continue
            else
               do 80, j = 1, n
                  do 70 k = m, 1, -1
                     if( b( k, j ).ne.zero )then
                        temp      = alpha*b( k, j )
                        b( k, j ) = temp
                        if( nounit )
     $                     b( k, j ) = b( k, j )*a( k, k )
                        do 60, i = k + 1, m
                           b( i, j ) = b( i, j ) + temp*a( i, k )
   60                   continue
                     end if
   70             continue
   80          continue
            end if
         else
*
*           form  b := alpha*b*a'.
*
            if( upper )then
               do 110, j = 1, n
                  do 100, i = m, 1, -1
                     temp = b( i, j )
                     if( nounit )
     $                  temp = temp*a( i, i )
                     do 90, k = 1, i - 1
                        temp = temp + a( k, i )*b( k, j )
   90                continue
                     b( i, j ) = alpha*temp
  100             continue
  110          continue
            else
               do 140, j = 1, n
                  do 130, i = 1, m
                     temp = b( i, j )
                     if( nounit )
     $                  temp = temp*a( i, i )
                     do 120, k = i + 1, m
                        temp = temp + a( k, i )*b( k, j )
  120                continue
                     b( i, j ) = alpha*temp
  130             continue
  140          continue
            end if
         end if
      else
         if( lsame( transa, 'n' ) )then
*
*           form  b := alpha*b*a.
*
            if( upper )then
               do 180, j = n, 1, -1
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 150, i = 1, m
                     b( i, j ) = temp*b( i, j )
  150             continue
                  do 170, k = 1, j - 1
                     if( a( k, j ).ne.zero )then
                        temp = alpha*a( k, j )
                        do 160, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  160                   continue
                     end if
  170             continue
  180          continue
            else
               do 220, j = 1, n
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 190, i = 1, m
                     b( i, j ) = temp*b( i, j )
  190             continue
                  do 210, k = j + 1, n
                     if( a( k, j ).ne.zero )then
                        temp = alpha*a( k, j )
                        do 200, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  200                   continue
                     end if
  210             continue
  220          continue
            end if
         else
*
*           form  b := alpha*b*a'.
*
            if( upper )then
               do 260, k = 1, n
                  do 240, j = 1, k - 1
                     if( a( j, k ).ne.zero )then
                        temp = alpha*a( j, k )
                        do 230, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  230                   continue
                     end if
  240             continue
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( k, k )
                  if( temp.ne.one )then
                     do 250, i = 1, m
                        b( i, k ) = temp*b( i, k )
  250                continue
                  end if
  260          continue
            else
               do 300, k = n, 1, -1
                  do 280, j = k + 1, n
                     if( a( j, k ).ne.zero )then
                        temp = alpha*a( j, k )
                        do 270, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  270                   continue
                     end if
  280             continue
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( k, k )
                  if( temp.ne.one )then
                     do 290, i = 1, m
                        b( i, k ) = temp*b( i, k )
  290                continue
                  end if
  300          continue
            end if
         end if
      end if
*
      return
*
*     end of dtrmm .
*
      end
