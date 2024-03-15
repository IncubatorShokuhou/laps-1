      subroutine drscl( n, sa, sx, incx )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     september 30, 1994
*
*     .. scalar arguments ..
      integer            incx, n
      double precision   sa
*     ..
*     .. array arguments ..
      double precision   sx( * )
*     ..
*
*  purpose
*  =======
*
*  drscl multiplies an n-element real vector x by the real scalar 1/a.
*  this is done without overflow or underflow as long as
*  the final result x/a does not overflow or underflow.
*
*  arguments
*  =========
*
*  n       (input) integer
*          the number of components of the vector x.
*
*  sa      (input) double precision
*          the scalar a which is used to divide each component of x.
*          sa must be >= 0, or the subroutine will divide by zero.
*
*  sx      (input/output) double precision array, dimension
*                         (1+(n-1)*abs(incx))
*          the n-element vector x.
*
*  incx    (input) integer
*          the increment between successive values of the vector sx.
*          > 0:  sx(1) = x(1) and sx(1+(i-1)*incx) = x(i),     1< i<= n
*
* =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      logical            done
      double precision   bignum, cden, cden1, cnum, cnum1, mul, smlnum
*     ..
*     .. external functions ..
      double precision   dlamch
      external           dlamch
*     ..
*     .. external subroutines ..
      external           dlabad, dscal
*     ..
*     .. intrinsic functions ..
      intrinsic          abs
*     ..
*     .. executable statements ..
*
*     quick return if possible
*
      if( n.le.0 )
     $   return
*
*     get machine parameters
*
      smlnum = dlamch( 's' )
      bignum = one / smlnum
      call dlabad( smlnum, bignum )
*
*     initialize the denominator to sa and the numerator to 1.
*
      cden = sa
      cnum = one
*
   10 continue
      cden1 = cden*smlnum
      cnum1 = cnum / bignum
      if( abs( cden1 ).gt.abs( cnum ) .and. cnum.ne.zero ) then
*
*        pre-multiply x by smlnum if cden is large compared to cnum.
*
         mul = smlnum
         done = .false.
         cden = cden1
      else if( abs( cnum1 ).gt.abs( cden ) ) then
*
*        pre-multiply x by bignum if cden is small compared to cnum.
*
         mul = bignum
         done = .false.
         cnum = cnum1
      else
*
*        multiply x by cnum / cden and return.
*
         mul = cnum / cden
         done = .true.
      end if
*
*     scale the vector x by mul
*
      call dscal( n, mul, sx, incx )
*
      if( .not.done )
     $   go to 10
*
      return
*
*     end of drscl
*
      end
