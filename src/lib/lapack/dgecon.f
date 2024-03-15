      subroutine dgecon( norm, n, a, lda, anorm, rcond, work, iwork,
     $                   info )
*
*  -- lapack routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     february 29, 1992
*
*     .. scalar arguments ..
      character          norm
      integer            info, lda, n
      double precision   anorm, rcond
*     ..
*     .. array arguments ..
      integer            iwork( * )
      double precision   a( lda, * ), work( * )
*     ..
*
*  purpose
*  =======
*
*  dgecon estimates the reciprocal of the condition number of a general
*  real matrix a, in either the 1-norm or the infinity-norm, using
*  the lu factorization computed by dgetrf.
*
*  an estimate is obtained for norm(inv(a)), and the reciprocal of the
*  condition number is computed as
*     rcond = 1 / ( norm(a) * norm(inv(a)) ).
*
*  arguments
*  =========
*
*  norm    (input) character*1
*          specifies whether the 1-norm condition number or the
*          infinity-norm condition number is required:
*          = '1' or 'o':  1-norm;
*          = 'i':         infinity-norm.
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  a       (input) double precision array, dimension (lda,n)
*          the factors l and u from the factorization a = p*l*u
*          as computed by dgetrf.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  anorm   (input) double precision
*          if norm = '1' or 'o', the 1-norm of the original matrix a.
*          if norm = 'i', the infinity-norm of the original matrix a.
*
*  rcond   (output) double precision
*          the reciprocal of the condition number of the matrix a,
*          computed as rcond = 1/(norm(a) * norm(inv(a))).
*
*  work    (workspace) double precision array, dimension (4*n)
*
*  iwork   (workspace) integer array, dimension (n)
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      logical            onenrm
      character          normin
      integer            ix, kase, kase1
      double precision   ainvnm, scale, sl, smlnum, su
*     ..
*     .. external functions ..
      logical            lsame
      integer            idamax
      double precision   dlamch
      external           lsame, idamax, dlamch
*     ..
*     .. external subroutines ..
      external           dlacon, dlatrs, drscl, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          abs, max
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      onenrm = norm.eq.'1' .or. lsame( norm, 'o' )
      if( .not.onenrm .and. .not.lsame( norm, 'i' ) ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, n ) ) then
         info = -4
      else if( anorm.lt.zero ) then
         info = -5
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgecon', -info )
         return
      end if
*
*     quick return if possible
*
      rcond = zero
      if( n.eq.0 ) then
         rcond = one
         return
      else if( anorm.eq.zero ) then
         return
      end if
*
      smlnum = dlamch( 'safe minimum' )
*
*     estimate the norm of inv(a).
*
      ainvnm = zero
      normin = 'n'
      if( onenrm ) then
         kase1 = 1
      else
         kase1 = 2
      end if
      kase = 0
   10 continue
      call dlacon( n, work( n+1 ), work, iwork, ainvnm, kase )
      if( kase.ne.0 ) then
         if( kase.eq.kase1 ) then
*
*           multiply by inv(l).
*
            call dlatrs( 'lower', 'no transpose', 'unit', normin, n, a,
     $                   lda, work, sl, work( 2*n+1 ), info )
*
*           multiply by inv(u).
*
            call dlatrs( 'upper', 'no transpose', 'non-unit', normin, n,
     $                   a, lda, work, su, work( 3*n+1 ), info )
         else
*
*           multiply by inv(u').
*
            call dlatrs( 'upper', 'transpose', 'non-unit', normin, n, a,
     $                   lda, work, su, work( 3*n+1 ), info )
*
*           multiply by inv(l').
*
            call dlatrs( 'lower', 'transpose', 'unit', normin, n, a,
     $                   lda, work, sl, work( 2*n+1 ), info )
         end if
*
*        divide x by 1/(sl*su) if doing so will not cause overflow.
*
         scale = sl*su
         normin = 'y'
         if( scale.ne.one ) then
            ix = idamax( n, work, 1 )
            if( scale.lt.abs( work( ix ) )*smlnum .or. scale.eq.zero )
     $         go to 20
            call drscl( n, scale, work, 1 )
         end if
         go to 10
      end if
*
*     compute the estimate of the reciprocal condition number.
*
      if( ainvnm.ne.zero )
     $   rcond = ( one / ainvnm ) / anorm
*
   20 continue
      return
*
*     end of dgecon
*
      end
