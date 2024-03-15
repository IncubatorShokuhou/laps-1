      subroutine xerbla( srname, info )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     september 30, 1994
*
*     .. scalar arguments ..
      character*6        srname
      integer            info
*     ..
*
*  purpose
*  =======
*
*  xerbla  is an error handler for the lapack routines.
*  it is called by an lapack routine if an input parameter has an
*  invalid value.  a message is printed and execution stops.
*
*  installers may consider modifying the stop statement in order to
*  call system-specific exception-handling facilities.
*
*  arguments
*  =========
*
*  srname  (input) character*6
*          the name of the routine which called xerbla.
*
*  info    (input) integer
*          the position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. executable statements ..
*
      write( *, fmt = 9999 )srname, info
*
      stop
*
 9999 format( ' ** on entry to ', a6, ' parameter number ', i2, ' had ',
     $      'an illegal value' )
*
*     end of xerbla
*
      end
