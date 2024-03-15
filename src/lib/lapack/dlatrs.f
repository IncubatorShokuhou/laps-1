      subroutine dlatrs( uplo, trans, diag, normin, n, a, lda, x, scale,
     $                   cnorm, info )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     june 30, 1992
*
*     .. scalar arguments ..
      character          diag, normin, trans, uplo
      integer            info, lda, n
      double precision   scale
*     ..
*     .. array arguments ..
      double precision   a( lda, * ), cnorm( * ), x( * )
*     ..
*
*  purpose
*  =======
*
*  dlatrs solves one of the triangular systems
*
*     a *x = s*b  or  a'*x = s*b
*
*  with scaling to prevent overflow.  here a is an upper or lower
*  triangular matrix, a' denotes the transpose of a, x and b are
*  n-element vectors, and s is a scaling factor, usually less than
*  or equal to 1, chosen so that the components of x will be less than
*  the overflow threshold.  if the unscaled problem will not cause
*  overflow, the level 2 blas routine dtrsv is called.  if the matrix a
*  is singular (a(j,j) = 0 for some j), then s is set to 0 and a
*  non-trivial solution to a*x = 0 is returned.
*
*  arguments
*  =========
*
*  uplo    (input) character*1
*          specifies whether the matrix a is upper or lower triangular.
*          = 'u':  upper triangular
*          = 'l':  lower triangular
*
*  trans   (input) character*1
*          specifies the operation applied to a.
*          = 'n':  solve a * x = s*b  (no transpose)
*          = 't':  solve a'* x = s*b  (transpose)
*          = 'c':  solve a'* x = s*b  (conjugate transpose = transpose)
*
*  diag    (input) character*1
*          specifies whether or not the matrix a is unit triangular.
*          = 'n':  non-unit triangular
*          = 'u':  unit triangular
*
*  normin  (input) character*1
*          specifies whether cnorm has been set or not.
*          = 'y':  cnorm contains the column norms on entry
*          = 'n':  cnorm is not set on entry.  on exit, the norms will
*                  be computed and stored in cnorm.
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  a       (input) double precision array, dimension (lda,n)
*          the triangular matrix a.  if uplo = 'u', the leading n by n
*          upper triangular part of the array a contains the upper
*          triangular matrix, and the strictly lower triangular part of
*          a is not referenced.  if uplo = 'l', the leading n by n lower
*          triangular part of the array a contains the lower triangular
*          matrix, and the strictly upper triangular part of a is not
*          referenced.  if diag = 'u', the diagonal elements of a are
*          also not referenced and are assumed to be 1.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max (1,n).
*
*  x       (input/output) double precision array, dimension (n)
*          on entry, the right hand side b of the triangular system.
*          on exit, x is overwritten by the solution vector x.
*
*  scale   (output) double precision
*          the scaling factor s for the triangular system
*             a * x = s*b  or  a'* x = s*b.
*          if scale = 0, the matrix a is singular or badly scaled, and
*          the vector x is an exact or approximate solution to a*x = 0.
*
*  cnorm   (input or output) double precision array, dimension (n)
*
*          if normin = 'y', cnorm is an input argument and cnorm(j)
*          contains the norm of the off-diagonal part of the j-th column
*          of a.  if trans = 'n', cnorm(j) must be greater than or equal
*          to the infinity-norm, and if trans = 't' or 'c', cnorm(j)
*          must be greater than or equal to the 1-norm.
*
*          if normin = 'n', cnorm is an output argument and cnorm(j)
*          returns the 1-norm of the offdiagonal part of the j-th column
*          of a.
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -k, the k-th argument had an illegal value
*
*  further details
*  ======= =======
*
*  a rough bound on x is computed; if that is less than overflow, dtrsv
*  is called, otherwise, specific code is used which checks for possible
*  overflow or divide-by-zero at every operation.
*
*  a columnwise scheme is used for solving a*x = b.  the basic algorithm
*  if a is lower triangular is
*
*       x[1:n] := b[1:n]
*       for j = 1, ..., n
*            x(j) := x(j) / a(j,j)
*            x[j+1:n] := x[j+1:n] - x(j) * a[j+1:n,j]
*       end
*
*  define bounds on the components of x after j iterations of the loop:
*     m(j) = bound on x[1:j]
*     g(j) = bound on x[j+1:n]
*  initially, let m(0) = 0 and g(0) = max{x(i), i=1,...,n}.
*
*  then for iteration j+1 we have
*     m(j+1) <= g(j) / | a(j+1,j+1) |
*     g(j+1) <= g(j) + m(j+1) * | a[j+2:n,j+1] |
*            <= g(j) ( 1 + cnorm(j+1) / | a(j+1,j+1) | )
*
*  where cnorm(j+1) is greater than or equal to the infinity-norm of
*  column j+1 of a, not counting the diagonal.  hence
*
*     g(j) <= g(0) product ( 1 + cnorm(i) / | a(i,i) | )
*                  1<=i<=j
*  and
*
*     |x(j)| <= ( g(0) / |a(j,j)| ) product ( 1 + cnorm(i) / |a(i,i)| )
*                                   1<=i< j
*
*  since |x(j)| <= m(j), we use the level 2 blas routine dtrsv if the
*  reciprocal of the largest m(j), j=1,..,n, is larger than
*  max(underflow, 1/overflow).
*
*  the bound on x(j) is also used to determine when a step in the
*  columnwise method can be performed without fear of overflow.  if
*  the computed bound is greater than a large constant, x is scaled to
*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*  1, and scale to 0, and a non-trivial solution to a*x = 0 is found.
*
*  similarly, a row-wise scheme is used to solve a'*x = b.  the basic
*  algorithm for a upper triangular is
*
*       for j = 1, ..., n
*            x(j) := ( b(j) - a[1:j-1,j]' * x[1:j-1] ) / a(j,j)
*       end
*
*  we simultaneously compute two bounds
*       g(j) = bound on ( b(i) - a[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
*       m(j) = bound on x(i), 1<=i<=j
*
*  the initial values are g(0) = 0, m(0) = max{b(i), i=1,..,n}, and we
*  add the constraint g(j) >= g(j-1) and m(j) >= m(j-1) for j >= 1.
*  then the bound on x(j) is
*
*       m(j) <= m(j-1) * ( 1 + cnorm(j) ) / | a(j,j) |
*
*            <= m(0) * product ( ( 1 + cnorm(i) ) / |a(i,i)| )
*                      1<=i<=j
*
*  and we can safely call dtrsv if 1/m(n) and 1/g(n) are both greater
*  than max(underflow, 1/overflow).
*
*  =====================================================================
*
*     .. parameters ..
      double precision   zero, half, one
      parameter          ( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0 )
*     ..
*     .. local scalars ..
      logical            notran, nounit, upper
      integer            i, imax, j, jfirst, jinc, jlast
      double precision   bignum, grow, rec, smlnum, sumj, tjj, tjjs,
     $                   tmax, tscal, uscal, xbnd, xj, xmax
*     ..
*     .. external functions ..
      logical            lsame
      integer            idamax
      double precision   dasum, ddot, dlamch
      external           lsame, idamax, dasum, ddot, dlamch
*     ..
*     .. external subroutines ..
      external           daxpy, dscal, dtrsv, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          abs, max, min
*     ..
*     .. executable statements ..
*
      info = 0
      upper = lsame( uplo, 'u' )
      notran = lsame( trans, 'n' )
      nounit = lsame( diag, 'n' )
*
*     test the input parameters.
*
      if( .not.upper .and. .not.lsame( uplo, 'l' ) ) then
         info = -1
      else if( .not.notran .and. .not.lsame( trans, 't' ) .and. .not.
     $         lsame( trans, 'c' ) ) then
         info = -2
      else if( .not.nounit .and. .not.lsame( diag, 'u' ) ) then
         info = -3
      else if( .not.lsame( normin, 'y' ) .and. .not.
     $         lsame( normin, 'n' ) ) then
         info = -4
      else if( n.lt.0 ) then
         info = -5
      else if( lda.lt.max( 1, n ) ) then
         info = -7
      end if
      if( info.ne.0 ) then
         call xerbla( 'dlatrs', -info )
         return
      end if
*
*     quick return if possible
*
      if( n.eq.0 )
     $   return
*
*     determine machine dependent parameters to control overflow.
*
      smlnum = dlamch( 'safe minimum' ) / dlamch( 'precision' )
      bignum = one / smlnum
      scale = one
*
      if( lsame( normin, 'n' ) ) then
*
*        compute the 1-norm of each column, not including the diagonal.
*
         if( upper ) then
*
*           a is upper triangular.
*
            do 10 j = 1, n
               cnorm( j ) = dasum( j-1, a( 1, j ), 1 )
   10       continue
         else
*
*           a is lower triangular.
*
            do 20 j = 1, n - 1
               cnorm( j ) = dasum( n-j, a( j+1, j ), 1 )
   20       continue
            cnorm( n ) = zero
         end if
      end if
*
*     scale the column norms by tscal if the maximum element in cnorm is
*     greater than bignum.
*
      imax = idamax( n, cnorm, 1 )
      tmax = cnorm( imax )
      if( tmax.le.bignum ) then
         tscal = one
      else
         tscal = one / ( smlnum*tmax )
         call dscal( n, tscal, cnorm, 1 )
      end if
*
*     compute a bound on the computed solution vector to see if the
*     level 2 blas routine dtrsv can be used.
*
      j = idamax( n, x, 1 )
      xmax = abs( x( j ) )
      xbnd = xmax
      if( notran ) then
*
*        compute the growth in a * x = b.
*
         if( upper ) then
            jfirst = n
            jlast = 1
            jinc = -1
         else
            jfirst = 1
            jlast = n
            jinc = 1
         end if
*
         if( tscal.ne.one ) then
            grow = zero
            go to 50
         end if
*
         if( nounit ) then
*
*           a is non-unit triangular.
*
*           compute grow = 1/g(j) and xbnd = 1/m(j).
*           initially, g(0) = max{x(i), i=1,...,n}.
*
            grow = one / max( xbnd, smlnum )
            xbnd = grow
            do 30 j = jfirst, jlast, jinc
*
*              exit the loop if the growth factor is too small.
*
               if( grow.le.smlnum )
     $            go to 50
*
*              m(j) = g(j-1) / abs(a(j,j))
*
               tjj = abs( a( j, j ) )
               xbnd = min( xbnd, min( one, tjj )*grow )
               if( tjj+cnorm( j ).ge.smlnum ) then
*
*                 g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
*
                  grow = grow*( tjj / ( tjj+cnorm( j ) ) )
               else
*
*                 g(j) could overflow, set grow to 0.
*
                  grow = zero
               end if
   30       continue
            grow = xbnd
         else
*
*           a is unit triangular.
*
*           compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
*
            grow = min( one, one / max( xbnd, smlnum ) )
            do 40 j = jfirst, jlast, jinc
*
*              exit the loop if the growth factor is too small.
*
               if( grow.le.smlnum )
     $            go to 50
*
*              g(j) = g(j-1)*( 1 + cnorm(j) )
*
               grow = grow*( one / ( one+cnorm( j ) ) )
   40       continue
         end if
   50    continue
*
      else
*
*        compute the growth in a' * x = b.
*
         if( upper ) then
            jfirst = 1
            jlast = n
            jinc = 1
         else
            jfirst = n
            jlast = 1
            jinc = -1
         end if
*
         if( tscal.ne.one ) then
            grow = zero
            go to 80
         end if
*
         if( nounit ) then
*
*           a is non-unit triangular.
*
*           compute grow = 1/g(j) and xbnd = 1/m(j).
*           initially, m(0) = max{x(i), i=1,...,n}.
*
            grow = one / max( xbnd, smlnum )
            xbnd = grow
            do 60 j = jfirst, jlast, jinc
*
*              exit the loop if the growth factor is too small.
*
               if( grow.le.smlnum )
     $            go to 80
*
*              g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
*
               xj = one + cnorm( j )
               grow = min( grow, xbnd / xj )
*
*              m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
*
               tjj = abs( a( j, j ) )
               if( xj.gt.tjj )
     $            xbnd = xbnd*( tjj / xj )
   60       continue
            grow = min( grow, xbnd )
         else
*
*           a is unit triangular.
*
*           compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
*
            grow = min( one, one / max( xbnd, smlnum ) )
            do 70 j = jfirst, jlast, jinc
*
*              exit the loop if the growth factor is too small.
*
               if( grow.le.smlnum )
     $            go to 80
*
*              g(j) = ( 1 + cnorm(j) )*g(j-1)
*
               xj = one + cnorm( j )
               grow = grow / xj
   70       continue
         end if
   80    continue
      end if
*
      if( ( grow*tscal ).gt.smlnum ) then
*
*        use the level 2 blas solve if the reciprocal of the bound on
*        elements of x is not too small.
*
         call dtrsv( uplo, trans, diag, n, a, lda, x, 1 )
      else
*
*        use a level 1 blas solve, scaling intermediate results.
*
         if( xmax.gt.bignum ) then
*
*           scale x so that its components are less than or equal to
*           bignum in absolute value.
*
            scale = bignum / xmax
            call dscal( n, scale, x, 1 )
            xmax = bignum
         end if
*
         if( notran ) then
*
*           solve a * x = b
*
            do 110 j = jfirst, jlast, jinc
*
*              compute x(j) = b(j) / a(j,j), scaling x if necessary.
*
               xj = abs( x( j ) )
               if( nounit ) then
                  tjjs = a( j, j )*tscal
               else
                  tjjs = tscal
                  if( tscal.eq.one )
     $               go to 100
               end if
               tjj = abs( tjjs )
               if( tjj.gt.smlnum ) then
*
*                    abs(a(j,j)) > smlnum:
*
                  if( tjj.lt.one ) then
                     if( xj.gt.tjj*bignum ) then
*
*                          scale x by 1/b(j).
*
                        rec = one / xj
                        call dscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     end if
                  end if
                  x( j ) = x( j ) / tjjs
                  xj = abs( x( j ) )
               else if( tjj.gt.zero ) then
*
*                    0 < abs(a(j,j)) <= smlnum:
*
                  if( xj.gt.tjj*bignum ) then
*
*                       scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
*                       to avoid overflow when dividing by a(j,j).
*
                     rec = ( tjj*bignum ) / xj
                     if( cnorm( j ).gt.one ) then
*
*                          scale by 1/cnorm(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        rec = rec / cnorm( j )
                     end if
                     call dscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  end if
                  x( j ) = x( j ) / tjjs
                  xj = abs( x( j ) )
               else
*
*                    a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to a*x = 0.
*
                  do 90 i = 1, n
                     x( i ) = zero
   90             continue
                  x( j ) = one
                  xj = one
                  scale = zero
                  xmax = zero
               end if
  100          continue
*
*              scale x if necessary to avoid overflow when adding a
*              multiple of column j of a.
*
               if( xj.gt.one ) then
                  rec = one / xj
                  if( cnorm( j ).gt.( bignum-xmax )*rec ) then
*
*                    scale x by 1/(2*abs(x(j))).
*
                     rec = rec*half
                     call dscal( n, rec, x, 1 )
                     scale = scale*rec
                  end if
               else if( xj*cnorm( j ).gt.( bignum-xmax ) ) then
*
*                 scale x by 1/2.
*
                  call dscal( n, half, x, 1 )
                  scale = scale*half
               end if
*
               if( upper ) then
                  if( j.gt.1 ) then
*
*                    compute the update
*                       x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
*
                     call daxpy( j-1, -x( j )*tscal, a( 1, j ), 1, x,
     $                           1 )
                     i = idamax( j-1, x, 1 )
                     xmax = abs( x( i ) )
                  end if
               else
                  if( j.lt.n ) then
*
*                    compute the update
*                       x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
*
                     call daxpy( n-j, -x( j )*tscal, a( j+1, j ), 1,
     $                           x( j+1 ), 1 )
                     i = j + idamax( n-j, x( j+1 ), 1 )
                     xmax = abs( x( i ) )
                  end if
               end if
  110       continue
*
         else
*
*           solve a' * x = b
*
            do 160 j = jfirst, jlast, jinc
*
*              compute x(j) = b(j) - sum a(k,j)*x(k).
*                                    k<>j
*
               xj = abs( x( j ) )
               uscal = tscal
               rec = one / max( xmax, one )
               if( cnorm( j ).gt.( bignum-xj )*rec ) then
*
*                 if x(j) could overflow, scale x by 1/(2*xmax).
*
                  rec = rec*half
                  if( nounit ) then
                     tjjs = a( j, j )*tscal
                  else
                     tjjs = tscal
                  end if
                  tjj = abs( tjjs )
                  if( tjj.gt.one ) then
*
*                       divide by a(j,j) when scaling x if a(j,j) > 1.
*
                     rec = min( one, rec*tjj )
                     uscal = uscal / tjjs
                  end if
                  if( rec.lt.one ) then
                     call dscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  end if
               end if
*
               sumj = zero
               if( uscal.eq.one ) then
*
*                 if the scaling needed for a in the dot product is 1,
*                 call ddot to perform the dot product.
*
                  if( upper ) then
                     sumj = ddot( j-1, a( 1, j ), 1, x, 1 )
                  else if( j.lt.n ) then
                     sumj = ddot( n-j, a( j+1, j ), 1, x( j+1 ), 1 )
                  end if
               else
*
*                 otherwise, use in-line code for the dot product.
*
                  if( upper ) then
                     do 120 i = 1, j - 1
                        sumj = sumj + ( a( i, j )*uscal )*x( i )
  120                continue
                  else if( j.lt.n ) then
                     do 130 i = j + 1, n
                        sumj = sumj + ( a( i, j )*uscal )*x( i )
  130                continue
                  end if
               end if
*
               if( uscal.eq.tscal ) then
*
*                 compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
*                 was not used to scale the dotproduct.
*
                  x( j ) = x( j ) - sumj
                  xj = abs( x( j ) )
                  if( nounit ) then
                     tjjs = a( j, j )*tscal
                  else
                     tjjs = tscal
                     if( tscal.eq.one )
     $                  go to 150
                  end if
*
*                    compute x(j) = x(j) / a(j,j), scaling if necessary.
*
                  tjj = abs( tjjs )
                  if( tjj.gt.smlnum ) then
*
*                       abs(a(j,j)) > smlnum:
*
                     if( tjj.lt.one ) then
                        if( xj.gt.tjj*bignum ) then
*
*                             scale x by 1/abs(x(j)).
*
                           rec = one / xj
                           call dscal( n, rec, x, 1 )
                           scale = scale*rec
                           xmax = xmax*rec
                        end if
                     end if
                     x( j ) = x( j ) / tjjs
                  else if( tjj.gt.zero ) then
*
*                       0 < abs(a(j,j)) <= smlnum:
*
                     if( xj.gt.tjj*bignum ) then
*
*                          scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
*
                        rec = ( tjj*bignum ) / xj
                        call dscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     end if
                     x( j ) = x( j ) / tjjs
                  else
*
*                       a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
*                       scale = 0, and compute a solution to a'*x = 0.
*
                     do 140 i = 1, n
                        x( i ) = zero
  140                continue
                     x( j ) = one
                     scale = zero
                     xmax = zero
                  end if
  150             continue
               else
*
*                 compute x(j) := x(j) / a(j,j)  - sumj if the dot
*                 product has already been divided by 1/a(j,j).
*
                  x( j ) = x( j ) / tjjs - sumj
               end if
               xmax = max( xmax, abs( x( j ) ) )
  160       continue
         end if
         scale = scale / tscal
      end if
*
*     scale the column norms by 1/tscal for return.
*
      if( tscal.ne.one ) then
         call dscal( n, one / tscal, cnorm, 1 )
      end if
*
      return
*
*     end of dlatrs
*
      end
