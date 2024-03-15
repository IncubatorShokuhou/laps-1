! searches a directory for wrf v2.1 arw files in netcdf format
!  (1 output time per file) that either match or surround (temporally)
!  a desired laps i4time
!
!  if it finds a perfect match, it returns nfiles=1 with the filename.
!  else, if there are two bounding files, nfiles=2
!  otherwise, nfiles=0

subroutine get_acceptable_wrf(bgpath,i4time_needed,nfiles, &
           i4times_found,filenames)

  use time_utils

  implicit none
  character(len=*), intent(in)  :: bgpath
  integer,intent(in)            :: i4time_needed
  integer,intent(out)           :: nfiles
  integer,intent(out)           :: i4times_found(2)
  character(len=256)            :: filenames(3) ! plus 1 dim. by wei-ting (130312) to put previous time
  character(len=256)            :: thisfile
  character(len=256)            :: all_files(20000)

  ! local vars
  integer  :: delta_t, delta_t_1, delta_t_2,total_files,i, istatus
  integer, parameter  :: max_files = 500
  integer, parameter  :: dtmax = 86400
  integer, parameter  :: dtmin = -86400
  character(len=11),parameter   :: filter = "wrfout_d01_"
  character(len=24)             :: wrftime
  integer                       :: i4time
  integer                       :: fname_len

  delta_t_1 = dtmax
  delta_t_2 = dtmin  
  call get_file_names(bgpath,total_files,all_files,max_files,istatus)

  do i = 1, total_files 

    thisfile = all_files(i)
    fname_len = len_trim(all_files(i))
    if (thisfile(fname_len-29:fname_len-19) .eq. filter) then
      wrftime = thisfile(fname_len-18:fname_len) // ".0000" 
      call mm5_to_lapstime(wrftime,i4time)
      if (i4time .eq. i4time_needed) then
         print *, "file matches time needed: ",trim(thisfile)
         nfiles = 1
         i4times_found(1) = i4time
         filenames(1)     = thisfile
         if (i .gt. 1) then
            filenames(2)  = all_files(i-1) ! modified by wei-ting (130326) to get previous time
         endif
         return
      else
       ! delta_t = i4time_needed - i4time
        delta_t = i4time - i4time_needed
        if (delta_t .gt. 0) then
          if (delta_t .lt. delta_t_1) then
            delta_t_1 = delta_t
            i4times_found(1) = i4time
            filenames(1) = thisfile
          endif
        else
          if (delta_t .gt. delta_t_2) then
            delta_t_2 = delta_t
            i4times_found(2) = i4time
            filenames(2) = thisfile
            if (i .gt. 1) then
               filenames(3) = all_files(i-1) ! added by wei-ting (130326) to get previous time
            endif
          endif
        endif
      endif
    endif
  enddo
  
  ! if we are here, then we did not find an exact
  ! match.  we need to check if we got appropriate
  ! bounding files
 
  if ( (i4times_found(1)-i4time_needed .lt. dtmax).and. & ! add -i4time_needed by wei-ting (130326)
       (i4times_found(2)-i4time_needed .gt. dtmin) ) then ! add -i4time_needed by wei-ting (130326)
     nfiles = 2
  else
     nfiles = 0
  endif
  return
end subroutine get_acceptable_wrf
