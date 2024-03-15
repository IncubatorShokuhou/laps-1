subroutine get_file_unit(unit)

   ! simple subroutine to get the next available file unit number

   implicit none
   integer, parameter   :: min_unit = 7
   integer, parameter   :: max_unit = 1023
   integer, intent(out) :: unit
   logical              :: unit_used
   logical              :: unit_found

   unit_found = .false.
   unit_loop: do unit = min_unit, max_unit
      inquire (unit=unit, opened=unit_used)
      if (.not. unit_used) then
         unit_found = .true.
         exit unit_loop
      end if
   end do unit_loop
   if (.not. unit_found) then
      print '(a)', 'no available unit numbers!'
      stop
   end if
   return
end subroutine get_file_unit
