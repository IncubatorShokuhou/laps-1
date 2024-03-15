!
! copyright (c) 2007  ; all rights reserved ; atmet, llc
! 
! this file is free software; you can redistribute it and/or modify it under the
! terms of the gnu general public license as published by the free software 
! foundation; either version 2 of the license, or (at your option) any later version.
! 
! this software is distributed in the hope that it will be useful, but without any 
! warranty; without even the implied warranty of merchantability or fitness for a 
! particular purpose.  see the gnu general public license for more details.
!
! you should have received a copy of the gnu general public license along with this 
! program; if not, write to the free software foundation, inc., 
! 59 temple place - suite 330, boston, ma 02111-1307, usa.
!======================================================================================

module mem_temp



type temp_fields
real, pointer, dimension(:,:,:)  :: t3, ht, p3
end type

type(temp_fields) :: temp

! pointers for renaming arrays
real, pointer, dimension(:,:,:) :: &
         temp_3d,heights_3d,pres_3d_pa

integer num_temp_obs

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine alloc_temp_arrays (nxl, nyl, nzl)

implicit none
integer :: nxl, nyl, nzl

integer :: nt

   allocate(temp%t3 (nxl,nyl,nzl))
   allocate(temp%ht (nxl,nyl,nzl))
   allocate(temp%p3 (nxl,nyl,nzl))

  
return
end subroutine alloc_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine point_temp_arrays ()

implicit none

temp_3d =>  temp%t3   
heights_3d =>  temp%ht  
pres_3d_pa  =>  temp%p3  
  
return
end subroutine point_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deallocate_temp_arrays()

implicit none


integer :: nt

   if (associated(temp%t3)) deallocate(temp%t3 )
   if (associated(temp%ht)) deallocate(temp%ht )
   if (associated(temp%p3)) deallocate(temp%p3 )

return
end subroutine deallocate_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nullify_temp_arrays ()

implicit none

integer :: nt

   if (associated(temp%t3)) nullify(temp%t3 )
   if (associated(temp%ht)) nullify(temp%ht )
   if (associated(temp%p3)) nullify(temp%p3 )

return
end subroutine nullify_temp_arrays
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module
