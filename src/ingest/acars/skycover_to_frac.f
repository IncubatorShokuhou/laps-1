

      subroutine skycover_to_frac(c_skycover,fraction,istatus)

      character*(*) c_skycover              ! input 
      real        fraction                ! output

      call s_len(c_skycover,ilen)
      ilen = max(ilen,1)

      istatus = 0
      fraction = 0.

      if(ilen .gt. 4)goto 999

!     calculate position of beginning of 3 and 4 character words
      ibeg3 = max(ilen-2,1)
      ibeg4 = max(ilen-3,1)

!     if(c_skycover(ibeg3:ilen) .eq. 'clr') fraction =  .01      
      if(c_skycover(ibeg3:ilen) .eq. 'sct') fraction =  .25      
      if(c_skycover(ibeg3:ilen) .eq. 'bkn') fraction =  .70     
      if(c_skycover(ibeg3:ilen) .eq. 'ovc') fraction = 1.00
      if(c_skycover(ibeg4:ilen) .eq. '-sct')fraction =  .15      
      if(c_skycover(ibeg4:ilen) .eq. '-bkn')fraction =  .40      
      if(c_skycover(ibeg4:ilen) .eq. '-ovc')fraction =  .60      
       
      if(fraction .ne. 0.)istatus = 1

 999  return
      end
