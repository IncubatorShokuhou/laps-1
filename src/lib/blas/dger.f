      subroutine dger  ( m, n, alpha, x, incx, y, incy, a, lda )
*     .. scalar arguments ..
      double precision   alpha
      integer            incx, incy, lda, m, n
*     .. array arguments ..
      double precision   a( lda, * ), x( * ), y( * )
*     ..
*
*  purpose
*  =======
*
*  dger   performs the rank 1 operation
*
*     a := alpha*x*y' + a,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and a is an m by n matrix.
*
*  parameters
*  ==========
*
*  m      - integer.
*           on entry, m specifies the number of rows of the matrix a.
*           m must be at least zero.
*           unchanged on exit.
*
*  n      - integer.
*           on entry, n specifies the number of columns of the matrix a.
*           n must be at least zero.
*           unchanged on exit.
*
*  alpha  - double precision.
*           on entry, alpha specifies the scalar alpha.
*           unchanged on exit.
*
*  x      - double precision array of dimension at least
*           ( 1 + ( m - 1 )*abs( incx ) ).
*           before entry, the incremented array x must contain the m
*           element vector x.
*           unchanged on exit.
*
*  incx   - integer.
*           on entry, incx specifies the increment for the elements of
*           x. incx must not be zero.
*           unchanged on exit.
*
*  y      - double precision array of dimension at least
*           ( 1 + ( n - 1 )*abs( incy ) ).
*           before entry, the incremented array y must contain the n
*           element vector y.
*           unchanged on exit.
*
*  incy   - integer.
*           on entry, incy specifies the increment for the elements of
*           y. incy must not be zero.
*           unchanged on exit.
*
*  a      - double precision array of dimension ( lda, n ).
*           before entry, the leading m by n part of the array a must
*           contain the matrix of coefficients. on exit, a is
*           overwritten by the updated matrix.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program. lda must be at least
*           max( 1, m ).
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
      integer            i, info, ix, j, jy, kx
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
      if     ( m.lt.0 )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, m ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'dger  ', info )
         return
      end if
*
*     quick return if possible.
*
      if( ( m.eq.0 ).or.( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
*
*     start the operations. in this version the elements of a are
*     accessed sequentially with one pass through a.
*
      if( incy.gt.0 )then
         jy = 1
      else
         jy = 1 - ( n - 1 )*incy
      end if
      if( incx.eq.1 )then
         do 20, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               do 10, i = 1, m
                  a( i, j ) = a( i, j ) + x( i )*temp
   10          continue
            end if
            jy = jy + incy
   20    continue
      else
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( m - 1 )*incx
         end if
         do 40, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               ix   = kx
               do 30, i = 1, m
                  a( i, j ) = a( i, j ) + x( ix )*temp
                  ix        = ix        + incx
   30          continue
            end if
            jy = jy + incy
   40    continue
      end if
*
      return
*
*     end of dger  .
*
      end
