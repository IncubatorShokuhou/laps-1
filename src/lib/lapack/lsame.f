      logical          function lsame( ca, cb )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     september 30, 1994
*
*     .. scalar arguments ..
      character          ca, cb
*     ..
*
*  purpose
*  =======
*
*  lsame returns .true. if ca is the same letter as cb regardless of
*  case.
*
*  arguments
*  =========
*
*  ca      (input) character*1
*  cb      (input) character*1
*          ca and cb specify the single characters to be compared.
*
* =====================================================================
*
*     .. intrinsic functions ..
      intrinsic          ichar
*     ..
*     .. local scalars ..
      integer            inta, intb, zcode
*     ..
*     .. executable statements ..
*
*     test if the characters are equal
*
      lsame = ca.eq.cb
      if( lsame )
     $   return
*
*     now test for equivalence if both characters are alphabetic.
*
      zcode = ichar( 'z' )
*
*     use 'z' rather than 'a' so that ascii can be detected on prime
*     machines, on which ichar returns a value with bit 8 set.
*     ichar('a') on prime machines returns 193 which is the same as
*     ichar('a') on an ebcdic machine.
*
      inta = ichar( ca )
      intb = ichar( cb )
*
      if( zcode.eq.90 .or. zcode.eq.122 ) then
*
*        ascii is assumed - zcode is the ascii code of either lower or
*        upper case 'z'.
*
         if( inta.ge.97 .and. inta.le.122 ) inta = inta - 32
         if( intb.ge.97 .and. intb.le.122 ) intb = intb - 32
*
      else if( zcode.eq.233 .or. zcode.eq.169 ) then
*
*        ebcdic is assumed - zcode is the ebcdic code of either lower or
*        upper case 'z'.
*
         if( inta.ge.129 .and. inta.le.137 .or.
     $       inta.ge.145 .and. inta.le.153 .or.
     $       inta.ge.162 .and. inta.le.169 ) inta = inta + 64
         if( intb.ge.129 .and. intb.le.137 .or.
     $       intb.ge.145 .and. intb.le.153 .or.
     $       intb.ge.162 .and. intb.le.169 ) intb = intb + 64
*
      else if( zcode.eq.218 .or. zcode.eq.250 ) then
*
*        ascii is assumed, on prime machines - zcode is the ascii code
*        plus 128 of either lower or upper case 'z'.
*
         if( inta.ge.225 .and. inta.le.250 ) inta = inta - 32
         if( intb.ge.225 .and. intb.le.250 ) intb = intb - 32
      end if
      lsame = inta.eq.intb
*
*     return
*
*     end of lsame
*
      end
