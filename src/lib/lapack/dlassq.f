      subroutine dlassq( n, x, incx, scale, sumsq )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      integer            incx, n
      double precision   scale, sumsq
*     ..
*     .. array arguments ..
      double precision   x( * )
*     ..
*
*  purpose
*  =======
*
*  dlassq  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = x( 1 + ( i - 1 )*incx ). the value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in scale and sumsq and
*  scl and smsq are overwritten on scale and sumsq respectively.
*
*  the routine makes only one pass through the vector x.
*
*  arguments
*  =========
*
*  n       (input) integer
*          the number of elements to be used from the vector x.
*
*  x       (input) double precision
*          the vector for which a scaled sum of squares is computed.
*             x( i )  = x( 1 + ( i - 1 )*incx ), 1 <= i <= n.
*
*  incx    (input) integer
*          the increment between successive values of the vector x.
*          incx > 0.
*
*  scale   (input/output) double precision
*          on entry, the value  scale  in the equation above.
*          on exit, scale is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  sumsq   (input/output) double precision
*          on entry, the value  sumsq  in the equation above.
*          on exit, sumsq is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. parameters ..
      double precision   zero
      parameter          ( zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      integer            ix
      double precision   absxi
*     ..
*     .. intrinsic functions ..
      intrinsic          abs
*     ..
*     .. executable statements ..
*
      if( n.gt.0 ) then
         do 10 ix = 1, 1 + ( n-1 )*incx, incx
            if( x( ix ).ne.zero ) then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi ) then
                  sumsq = 1 + sumsq*( scale / absxi )**2
                  scale = absxi
               else
                  sumsq = sumsq + ( absxi / scale )**2
               end if
            end if
   10    continue
      end if
      return
*
*     end of dlassq
*
      end
