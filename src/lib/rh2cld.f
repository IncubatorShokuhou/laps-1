

        function rh_to_cldcv(rh)

!       1999 steve albers fsl

!       convert fractional rh into fractional cloud cover. this is defined
!       as a total cloud cover through a particular cloud layer.
!       it is recommended that the input rh be calculated with subroutine
!       'make_rh' using a 't_ref' between 0c and -10c.

        ramp_thresh = 0.80

        if(rh .lt. ramp_thresh)then
            rh_to_cldcv = 0.0
        elseif(rh .le. 1.00)then
            rh_to_cldcv = (rh - ramp_thresh) / (1.00 - ramp_thresh)  
        else
            rh_to_cldcv = 1.00
        endif

        return
        end
