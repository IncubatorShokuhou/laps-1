      subroutine dlabad( small, large )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      double precision   large, small
*     ..
*
*  purpose
*  =======
*
*  dlabad takes as input the values computed by slamch for underflow and
*  overflow, and returns the square root of each of these values if the
*  log of large is sufficiently large.  this subroutine is intended to
*  identify machines with a large exponent range, such as the crays, and
*  redefine the underflow and overflow limits to be the square roots of
*  the values computed by dlamch.  this subroutine is needed because
*  dlamch does not compensate for poor arithmetic in the upper half of
*  the exponent range, as is found on a cray.
*
*  arguments
*  =========
*
*  small   (input/output) double precision
*          on entry, the underflow threshold as computed by dlamch.
*          on exit, if log10(large) is sufficiently large, the square
*          root of small, otherwise unchanged.
*
*  large   (input/output) double precision
*          on entry, the overflow threshold as computed by dlamch.
*          on exit, if log10(large) is sufficiently large, the square
*          root of large, otherwise unchanged.
*
*  =====================================================================
*
*     .. intrinsic functions ..
      intrinsic          log10, sqrt
*     ..
*     .. executable statements ..
*
*     if it looks like we're on a cray, take the square root of
*     small and large to avoid overflow and underflow problems.
*
      if( log10( large ).gt.2000.d0 ) then
         small = sqrt( small )
         large = sqrt( large )
      end if
*
      return
*
*     end of dlabad
*
      end
