      integer          function ilaenv( ispec, name, opts, n1, n2, n3,
     $                 n4 )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     september 30, 1994
*
*     .. scalar arguments ..
      character*( * )    name, opts
      integer            ispec, n1, n2, n3, n4
*     ..
*
*  purpose
*  =======
*
*  ilaenv is called from the lapack routines to choose problem-dependent
*  parameters for the local environment.  see ispec for a description of
*  the parameters.
*
*  this version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  this routine will not function correctly if it is converted to all
*  lower case.  converting it to all upper case is allowed.
*
*  arguments
*  =========
*
*  ispec   (input) integer
*          specifies the parameter to be returned as the value of
*          ilaenv.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for n less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ilaenv(2,...) and m by ilaenv(5,...)
*          = 6: the crossover point for the svd (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a qr factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift qr and qz methods
*               for nonsymmetric eigenvalue problems.
*
*  name    (input) character*(*)
*          the name of the calling subroutine, in either upper case or
*          lower case.
*
*  opts    (input) character*(*)
*          the character options to the subroutine name, concatenated
*          into a single character string.  for example, uplo = 'u',
*          trans = 't', and diag = 'n' for a triangular routine would
*          be specified as opts = 'utn'.
*
*  n1      (input) integer
*  n2      (input) integer
*  n3      (input) integer
*  n4      (input) integer
*          problem dimensions for the subroutine name; these may not all
*          be required.
*
* (ilaenv) (output) integer
*          >= 0: the value of the parameter specified by ispec
*          < 0:  if ilaenv = -k, the k-th argument had an illegal value.
*
*  further details
*  ===============
*
*  the following conventions have been used when calling ilaenv from the
*  lapack routines:
*  1)  opts is a concatenation of all of the character options to
*      subroutine name, in the same order that they appear in the
*      argument list for name, even if they are not used in determining
*      the value of the parameter specified by ispec.
*  2)  the problem dimensions n1, n2, n3, n4 are specified in the order
*      that they appear in the argument list for name.  n1 is used
*      first, n2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  the parameter value returned by ilaenv is checked for validity in
*      the calling subroutine.  for example, ilaenv is used to retrieve
*      the optimal blocksize for strtri as follows:
*
*      nb = ilaenv( 1, 'strtri', uplo // diag, n, -1, -1, -1 )
*      if( nb.le.1 ) nb = max( 1, n )
*
*  =====================================================================
*
*     .. local scalars ..
      logical            cname, sname
      character*1        c1
      character*2        c2, c4
      character*3        c3
      character*6        subnam
      integer            i, ic, iz, nb, nbmin, nx
*     ..
*     .. intrinsic functions ..
      intrinsic          char, ichar, int, min, real
*     ..
*     .. executable statements ..
*
      go to ( 100, 100, 100, 400, 500, 600, 700, 800 ) ispec
*
*     invalid value for ispec
*
      ilaenv = -1
      return
*
  100 continue
*
*     convert name to upper case if the first character is lower case.
*
      ilaenv = 1
      subnam = name
      ic = ichar( subnam( 1:1 ) )
      iz = ichar( 'z' )
      if( iz.eq.90 .or. iz.eq.122 ) then
*
*        ascii character set
*
         if( ic.ge.97 .and. ic.le.122 ) then
            subnam( 1:1 ) = char( ic-32 )
            do 10 i = 2, 6
               ic = ichar( subnam( i:i ) )
               if( ic.ge.97 .and. ic.le.122 )
     $            subnam( i:i ) = char( ic-32 )
   10       continue
         end if
*
      else if( iz.eq.233 .or. iz.eq.169 ) then
*
*        ebcdic character set
*
         if( ( ic.ge.129 .and. ic.le.137 ) .or.
     $       ( ic.ge.145 .and. ic.le.153 ) .or.
     $       ( ic.ge.162 .and. ic.le.169 ) ) then
            subnam( 1:1 ) = char( ic+64 )
            do 20 i = 2, 6
               ic = ichar( subnam( i:i ) )
               if( ( ic.ge.129 .and. ic.le.137 ) .or.
     $             ( ic.ge.145 .and. ic.le.153 ) .or.
     $             ( ic.ge.162 .and. ic.le.169 ) )
     $            subnam( i:i ) = char( ic+64 )
   20       continue
         end if
*
      else if( iz.eq.218 .or. iz.eq.250 ) then
*
*        prime machines:  ascii+128
*
         if( ic.ge.225 .and. ic.le.250 ) then
            subnam( 1:1 ) = char( ic-32 )
            do 30 i = 2, 6
               ic = ichar( subnam( i:i ) )
               if( ic.ge.225 .and. ic.le.250 )
     $            subnam( i:i ) = char( ic-32 )
   30       continue
         end if
      end if
*
      c1 = subnam( 1:1 )
      sname = c1.eq.'s' .or. c1.eq.'d'
      cname = c1.eq.'c' .or. c1.eq.'z'
      if( .not.( cname .or. sname ) )
     $   return
      c2 = subnam( 2:3 )
      c3 = subnam( 4:6 )
      c4 = c3( 2:3 )
*
      go to ( 110, 200, 300 ) ispec
*
  110 continue
*
*     ispec = 1:  block size
*
*     in these examples, separate code is provided for setting nb for
*     real and complex.  we assume that nb will take the same value in
*     single or double precision.
*
      nb = 1
*
      if( c2.eq.'ge' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         else if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or.
     $            c3.eq.'qlf' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'tri' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'po' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         else if( sname .and. c3.eq.'trd' ) then
            nb = 1
         else if( sname .and. c3.eq.'gst' ) then
            nb = 64
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trf' ) then
            nb = 64
         else if( c3.eq.'trd' ) then
            nb = 1
         else if( c3.eq.'gst' ) then
            nb = 64
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1:1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nb = 32
            end if
         else if( c3( 1:1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nb = 32
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1:1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nb = 32
            end if
         else if( c3( 1:1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nb = 32
            end if
         end if
      else if( c2.eq.'gb' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               if( n4.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            else
               if( n4.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            end if
         end if
      else if( c2.eq.'pb' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               if( n2.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            else
               if( n2.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            end if
         end if
      else if( c2.eq.'tr' ) then
         if( c3.eq.'tri' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'la' ) then
         if( c3.eq.'uum' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( sname .and. c2.eq.'st' ) then
         if( c3.eq.'ebz' ) then
            nb = 1
         end if
      end if
      ilaenv = nb
      return
*
  200 continue
*
*     ispec = 2:  minimum block size
*
      nbmin = 2
      if( c2.eq.'ge' ) then
         if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or.
     $       c3.eq.'qlf' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'tri' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nbmin = 8
            else
               nbmin = 8
            end if
         else if( sname .and. c3.eq.'trd' ) then
            nbmin = 2
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trd' ) then
            nbmin = 2
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1:1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nbmin = 2
            end if
         else if( c3( 1:1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nbmin = 2
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1:1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nbmin = 2
            end if
         else if( c3( 1:1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nbmin = 2
            end if
         end if
      end if
      ilaenv = nbmin
      return
*
  300 continue
*
*     ispec = 3:  crossover point
*
      nx = 0
      if( c2.eq.'ge' ) then
         if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or.
     $       c3.eq.'qlf' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( sname .and. c3.eq.'trd' ) then
            nx = 1
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trd' ) then
            nx = 1
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1:1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nx = 128
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1:1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or.
     $          c4.eq.'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or.
     $          c4.eq.'br' ) then
               nx = 128
            end if
         end if
      end if
      ilaenv = nx
      return
*
  400 continue
*
*     ispec = 4:  number of shifts (used by xhseqr)
*
      ilaenv = 6
      return
*
  500 continue
*
*     ispec = 5:  minimum column dimension (not used)
*
      ilaenv = 2
      return
*
  600 continue 
*
*     ispec = 6:  crossover point for svd (used by xgelss and xgesvd)
*
      ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
      return
*
  700 continue
*
*     ispec = 7:  number of processors (not used)
*
      ilaenv = 1
      return
*
  800 continue
*
*     ispec = 8:  crossover point for multishift (used by xhseqr)
*
      ilaenv = 50
      return
*
*     end of ilaenv
*
      end
