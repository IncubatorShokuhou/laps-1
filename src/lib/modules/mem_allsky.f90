
module mem_allsky

!     input arrays on model grid
      real, allocatable, dimension(:,:,:) :: pres_3d
      real, allocatable, dimension(:,:,:) :: heights_3d
      real, allocatable, dimension(:,:,:) :: clwc_3d    ! kg/m**3
      real, allocatable, dimension(:,:,:) :: cice_3d    ! kg/m**3
      real, allocatable, dimension(:,:,:) :: rain_3d    ! kg/m**3
      real, allocatable, dimension(:,:,:) :: snow_3d    ! kg/m**3
      real, allocatable, dimension(:,:,:) :: aod_3d

!     local arrays on model grid
      real, allocatable, dimension(:,:,:) :: transm_3d
      real, allocatable, dimension(:,:,:,:) :: transm_4d
      real, allocatable, dimension(:,:,:,:) :: uprad_4d 
      real, allocatable, dimension(:,:,:) :: upxrad_3d 
      real, allocatable, dimension(:,:,:) :: upyrad_3d 

!     2d arrays on sky grid
      real, allocatable, dimension(:,:) :: aod_ill_opac
      real, allocatable, dimension(:,:) :: aod_ill_opac_potl

!     various non-gridded variables
      parameter (nc = 3)
      parameter (day_int0 = 3e9) ! solar relative brightness scaling constant

!     https://www.goes-r.gov/products/atbds/baseline/imagery_v2.0_no_color.pdf
      real nl_2_reflectance
      parameter (nl_2_reflectance = 1. / (4. * day_int0))
      real ghi_sim
      real ext_g(nc)
      integer mode_aero_cld /1/ ! treat aerosols more as clouds [1,2,3]
      integer mil,mih,mjl,mjh

      public alloc_allsky, dealloc_allsky

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine alloc_allsky(ni,nj,nk,nc,istatus)   ! i/o

!       allocate some though not all arrays mentioned above

        allocate(pres_3d(ni,nj,nk))
        allocate(heights_3d(ni,nj,nk))
        allocate(clwc_3d(ni,nj,nk))
        allocate(cice_3d(ni,nj,nk))
        allocate(rain_3d(ni,nj,nk))
        allocate(snow_3d(ni,nj,nk))
        allocate(aod_3d(ni,nj,nk))
        allocate(transm_3d(ni,nj,nk))
        allocate(transm_4d(ni,nj,nk,nc))
        allocate(uprad_4d(ni,nj,nk,nc))
        allocate(upxrad_3d(ni,nj,nk))
        allocate(upyrad_3d(ni,nj,nk))

        write(6,*)' allsky successfully allocated'

        return
   
     end subroutine alloc_allsky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine dealloc_allsky()             

!       deallocate some though not all arrays mentioned above

        deallocate(pres_3d)
        deallocate(heights_3d)
        deallocate(clwc_3d)
        deallocate(cice_3d)
        deallocate(rain_3d)
        deallocate(snow_3d)
        deallocate(aod_3d)
        deallocate(transm_3d)
        deallocate(transm_4d)
        deallocate(uprad_4d)
        deallocate(upxrad_3d)
        deallocate(upyrad_3d)

        write(6,*)' allsky successfully deallocated'

        return
   
     end subroutine dealloc_allsky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mem_allsky
