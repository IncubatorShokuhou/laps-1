      subroutine dlacon( n, v, x, isgn, est, kase )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     february 29, 1992
*
*     .. scalar arguments ..
      integer            kase, n
      double precision   est
*     ..
*     .. array arguments ..
      integer            isgn( * )
      double precision   v( * ), x( * )
*     ..
*
*  purpose
*  =======
*
*  dlacon estimates the 1-norm of a square, real matrix a.
*  reverse communication is used for evaluating matrix-vector products.
*
*  arguments
*  =========
*
*  n      (input) integer
*         the order of the matrix.  n >= 1.
*
*  v      (workspace) double precision array, dimension (n)
*         on the final return, v = a*w,  where  est = norm(v)/norm(w)
*         (w is not returned).
*
*  x      (input/output) double precision array, dimension (n)
*         on an intermediate return, x should be overwritten by
*               a * x,   if kase=1,
*               a' * x,  if kase=2,
*         and dlacon must be re-called with all the other parameters
*         unchanged.
*
*  isgn   (workspace) integer array, dimension (n)
*
*  est    (output) double precision
*         an estimate (a lower bound) for norm(a).
*
*  kase   (input/output) integer
*         on the initial call to dlacon, kase should be 0.
*         on an intermediate return, kase will be 1 or 2, indicating
*         whether x should be overwritten by a * x  or a' * x.
*         on the final return from dlacon, kase will again be 0.
*
*  further details
*  ======= =======
*
*  contributed by nick higham, university of manchester.
*  originally named sonest, dated march 16, 1988.
*
*  reference: n.j. higham, "fortran codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  acm trans. math. soft., vol. 14, no. 4, pp. 381-396, december 1988.
*
*  =====================================================================
*
*     .. parameters ..
      integer            itmax
      parameter          ( itmax = 5 )
      double precision   zero, one, two
      parameter          ( zero = 0.0d+0, one = 1.0d+0, two = 2.0d+0 )
*     ..
*     .. local scalars ..
      integer            i, iter, j, jlast, jump
      double precision   altsgn, estold, temp
*     ..
*     .. external functions ..
      integer            idamax
      double precision   dasum
      external           idamax, dasum
*     ..
*     .. external subroutines ..
      external           dcopy
*     ..
*     .. intrinsic functions ..
      intrinsic          abs, dble, nint, sign
*     ..
*     .. save statement ..
      save
*     ..
*     .. executable statements ..
*
      if( kase.eq.0 ) then
         do 10 i = 1, n
            x( i ) = one / dble( n )
   10    continue
         kase = 1
         jump = 1
         return
      end if
*
      go to ( 20, 40, 70, 110, 140 )jump
*
*     ................ entry   (jump = 1)
*     first iteration.  x has been overwritten by a*x.
*
   20 continue
      if( n.eq.1 ) then
         v( 1 ) = x( 1 )
         est = abs( v( 1 ) )
*        ... quit
         go to 150
      end if
      est = dasum( n, x, 1 )
*
      do 30 i = 1, n
         x( i ) = sign( one, x( i ) )
         isgn( i ) = nint( x( i ) )
   30 continue
      kase = 2
      jump = 2
      return
*
*     ................ entry   (jump = 2)
*     first iteration.  x has been overwritten by trandpose(a)*x.
*
   40 continue
      j = idamax( n, x, 1 )
      iter = 2
*
*     main loop - iterations 2,3,...,itmax.
*
   50 continue
      do 60 i = 1, n
         x( i ) = zero
   60 continue
      x( j ) = one
      kase = 1
      jump = 3
      return
*
*     ................ entry   (jump = 3)
*     x has been overwritten by a*x.
*
   70 continue
      call dcopy( n, x, 1, v, 1 )
      estold = est
      est = dasum( n, v, 1 )
      do 80 i = 1, n
         if( nint( sign( one, x( i ) ) ).ne.isgn( i ) )
     $      go to 90
   80 continue
*     repeated sign vector detected, hence algorithm has converged.
      go to 120
*
   90 continue
*     test for cycling.
      if( est.le.estold )
     $   go to 120
*
      do 100 i = 1, n
         x( i ) = sign( one, x( i ) )
         isgn( i ) = nint( x( i ) )
  100 continue
      kase = 2
      jump = 4
      return
*
*     ................ entry   (jump = 4)
*     x has been overwritten by trandpose(a)*x.
*
  110 continue
      jlast = j
      j = idamax( n, x, 1 )
      if( ( x( jlast ).ne.abs( x( j ) ) ) .and. ( iter.lt.itmax ) ) then
         iter = iter + 1
         go to 50
      end if
*
*     iteration complete.  final stage.
*
  120 continue
      altsgn = one
      do 130 i = 1, n
         x( i ) = altsgn*( one+dble( i-1 ) / dble( n-1 ) )
         altsgn = -altsgn
  130 continue
      kase = 1
      jump = 5
      return
*
*     ................ entry   (jump = 5)
*     x has been overwritten by a*x.
*
  140 continue
      temp = two*( dasum( n, x, 1 ) / dble( 3*n ) )
      if( temp.gt.est ) then
         call dcopy( n, x, 1, v, 1 )
         est = temp
      end if
*
  150 continue
      kase = 0
      return
*
*     end of dlacon
*
      end
