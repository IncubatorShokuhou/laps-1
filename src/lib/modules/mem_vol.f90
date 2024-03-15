module mem_vol  

integer, allocatable, dimension(:,:,:) :: reflectivity, reflectivity_hi, radialvelocity, radialvelocity_hi

real, allocatable, dimension(:,:) :: elevationr, elevationr_hi, elevationv, elevationv_hi          

real, allocatable, dimension(:,:) :: azimuthr, azimuthr_hi, azimuthv, azimuthv_hi          

real, allocatable, dimension(:) :: distancer, distancer_hi, distancev, distancev_hi, nyquistvelocityv, nyquistvelocityv_hi          

integer :: gater, gater_hi, gatev, gatev_hi, radialr, radialr_hi, &
           radialv, radialv_hi, scanr, scanr_hi, scanv, &
           scanv_hi,nf_fid, nf_vid, nf_status  

character*8 :: c8_fname_format

end module
