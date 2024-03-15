! this analysis program is based on the idea of dr. xie in xie et al. (2005) paper         !
! 'a sequential variational analysis approach for mesoscale data assimilation', in         !
! which an innovation of data assimilation has been proposed that data assimilation        !
! should be decomposed to two step, with the first step to catch information provided      !
! by observations and the second step to combine background and observations to make       !
! the final optimal analysis. this code is to implement the first step which can provide   !
! reasonable analysis field. for high efficiency, we choose the multi-grid method to       !
! implement this idea. the previous version of this code is the stmas surface analysis     !
! system coded by dr. xie himself. some applicaiton of this surface stmas in ocean can     !
! be refered to li et al. (2007) and he et al. (2007). in this four-dimensional version,   !
! we use some penalty term following the suggestion in xie et al. (2002) paper 'impact     !
! of formulation ofcost function and constraionts on three-dimensional variational data    !
! assimialtion'. so far, this four-dimensional version can deal with analysis from 1d to   !
! 4d, multi-variable anaysis using balance and general vertical coordinate.                !

! name of local variable               has  1 or 2 characters
! physical scale variable              has     3   characters
! name of i/o device unit              has     5   characters
! name of integer global varialbe      has     7   characters
! name of parameter                    has     7   characters
! name of global array                 has     7   characters
! name of observation or grid variable has     8   characters
! name of subroutine                   has     9   characters
! name of module                       has    12   characters

program main

   use input_bg_obs, only: bkgrndobs
   use prmtrs_stmas, only: if_test
   use stmas4d_core, only: mganalyss
   use output_anals, only: tmpoutput, outptlaps, outputana

   call bkgrndobs
   call mganalyss
   if (if_test .eq. 0) then
      call outptlaps
   else
      call outputana
   end if
   print *, 'end of analysis'

end program main
