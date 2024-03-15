      double precision function dlange( norm, m, n, a, lda, work )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      character          norm
      integer            lda, m, n
*     ..
*     .. array arguments ..
      double precision   a( lda, * ), work( * )
*     ..
*
*  purpose
*  =======
*
*  dlange  returns the value of the one norm,  or the frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real matrix a.
*
*  description
*  ===========
*
*  dlange returns the value
*
*     dlange = ( max(abs(a(i,j))), norm = 'm' or 'm'
*              (
*              ( norm1(a),         norm = '1', 'o' or 'o'
*              (
*              ( normi(a),         norm = 'i' or 'i'
*              (
*              ( normf(a),         norm = 'f', 'f', 'e' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normi  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normf  denotes the  frobenius norm of a matrix (square root of sum of
*  squares).  note that  max(abs(a(i,j)))  is not a  matrix norm.
*
*  arguments
*  =========
*
*  norm    (input) character*1
*          specifies the value to be returned in dlange as described
*          above.
*
*  m       (input) integer
*          the number of rows of the matrix a.  m >= 0.  when m = 0,
*          dlange is set to zero.
*
*  n       (input) integer
*          the number of columns of the matrix a.  n >= 0.  when n = 0,
*          dlange is set to zero.
*
*  a       (input) double precision array, dimension (lda,n)
*          the m by n matrix a.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(m,1).
*
*  work    (workspace) double precision array, dimension (lwork),
*          where lwork >= m when norm = 'i'; otherwise, work is not
*          referenced.
*
* =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      integer            i, j
      double precision   scale, sum, value
*     ..
*     .. external subroutines ..
      external           dlassq
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. intrinsic functions ..
      intrinsic          abs, max, min, sqrt
*     ..
*     .. executable statements ..
*
      if( min( m, n ).eq.0 ) then
         value = zero
      else if( lsame( norm, 'm' ) ) then
*
*        find max(abs(a(i,j))).
*
         value = zero
         do 20 j = 1, n
            do 10 i = 1, m
               value = max( value, abs( a( i, j ) ) )
   10       continue
   20    continue
      else if( ( lsame( norm, 'o' ) ) .or. ( norm.eq.'1' ) ) then
*
*        find norm1(a).
*
         value = zero
         do 40 j = 1, n
            sum = zero
            do 30 i = 1, m
               sum = sum + abs( a( i, j ) )
   30       continue
            value = max( value, sum )
   40    continue
      else if( lsame( norm, 'i' ) ) then
*
*        find normi(a).
*
         do 50 i = 1, m
            work( i ) = zero
   50    continue
         do 70 j = 1, n
            do 60 i = 1, m
               work( i ) = work( i ) + abs( a( i, j ) )
   60       continue
   70    continue
         value = zero
         do 80 i = 1, m
            value = max( value, work( i ) )
   80    continue
      else if( ( lsame( norm, 'f' ) ) .or. ( lsame( norm, 'e' ) ) ) then
*
*        find normf(a).
*
         scale = zero
         sum = one
         do 90 j = 1, n
            call dlassq( m, a( 1, j ), 1, scale, sum )
   90    continue
         value = scale*sqrt( sum )
      end if
*
      dlange = value
      return
*
*     end of dlange
*
      end
