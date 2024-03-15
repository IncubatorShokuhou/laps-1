      module module_sfc_structure

      type :: lbsi
      
      real :: lat ! deg
      real :: lon ! deg
      real :: sfc_temp ! k
      real :: sfc_pres ! mb
      real :: secsola ! no units
      real :: secza(3)  !no units, 1=goese, 2=goesw, 3=goes 9
      real :: sfc_emiss (25) ! unknown units
      real :: sfc_refl (25) ! unknown units

      end type lbsi

      end module module_sfc_structure
