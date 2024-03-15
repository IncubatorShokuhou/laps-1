!  output in netcdf format (assumes 1 forecast time per
!  file)
!
!  brent l. shaw, weathernews americas, inc., 2005

subroutine lga_driver_wrfarw(nx_laps,ny_laps,nz_laps,bgpath, cmodel, & ! add nx_laps,ny_laps,nz_laps by wei-ting(130326)
                             use_analysis,forecast_length, &
                             luse_sfc_bkgd, &
                             i4time_now, smooth_fields,lga_status)

  use wrf_netcdf
  use time_utils
  use wrf_lga
  implicit none
  include 'netcdf.inc'
  ! === add variables by wei-ting(130326) ===
  integer, intent(in)          :: nx_laps,ny_laps,nz_laps ! laps grid dimensions
  real                         :: ht_1d(nz_laps), &       !laps vert grid (sigma_ht grid only)
                                  pr1d_pa(nz_laps), &     !laps pressures (pa)
                                  pr1d_mb(nz_laps)        !laps pressures (mb)
  character(len=9)             :: wfo_fname13_to_fname9
  character(len=9)             :: fname9_reftime          ! reftime julian date
  integer                      :: i4time_reftime,i4time_fcst(2)
  character(len=256)           :: outdir
  character(len=31)            :: ext
  integer                      :: iloop,len_dir,hr_now,min_now
  ! === end of add by wei-ting(130326) ===
  character(len=*),intent(in)  :: bgpath
  character(len=12)            :: cmodel
  logical, intent(in)          :: use_analysis
  integer, intent(in)          :: forecast_length ! hours
  logical, intent(in)          :: luse_sfc_bkgd
  integer, intent(in)          :: i4time_now 
  logical, intent(in)          :: smooth_fields
  integer, intent(out)         :: lga_status ! 1 = success

  ! local variables
  integer              :: nfiles,cdf,istatus
  integer              :: files_i4time(2)  
  character(len=256)   :: filenames(3) ! plus 1 dim. by wei-ting (130312) to put previous time
  character(len=19)    :: reftime
  real                 :: dt
  integer              :: itimestep, tau_hr, tau_min, tau_sec 
  logical              :: already_done
  lga_status = 1  ! if we encounter a failure, change to 0


  ! get list of acceptable files
  call get_acceptable_wrf(bgpath,i4time_now,nfiles, &
     files_i4time, filenames)


  ! do we have an exact time match?  
  if (nfiles .eq. 1) then
    ! yes, get its attributes
    print *, "found exact match: ", trim(filenames(1))
    ! make sure file is ready
    call wrfio_wait(filenames(1),300)
    call open_wrfnc(filenames(1),cdf,istatus)
    if (istatus .ne. 0) then
      print *, "problem opening ",trim(filenames(1))
      lga_status = 0
      return 
    endif
    ! get time information
    call get_wrf2_timeinfo(cdf,reftime,dt,itimestep,tau_hr, &
           tau_min, tau_sec, istatus)
    call close_wrfnc(cdf)
    if (istatus .ne. 0) then
      print*, "problem getting time info from ", trim(filenames(1))
      lga_status = 0
      return
    endif
    ! is it a new enough forecast?
    if (tau_hr .le. forecast_length) then
      print *, "forecast hour = ",tau_hr

      ! only process if we don't already have an lga output file from wrfarw
      ! that matches this reftime+tau
      print *, "calling check_wrf_lga for:",cmodel,reftime
      call check_wrf_lga(cmodel,reftime,tau_hr,tau_min,already_done) 
      print *, "already done: ",already_done
      if (.not. already_done) then
        print *, "calling wrf2lga"
        call wrf2lga(filenames(1:2),i4time_now,cmodel,istatus) ! modified by wei-ting (130326) to insert previous time ( filename(1) -> filname(1:2) )
        if (istatus .ne. 1) then
          print *, "failure in wrf2lga"
          lga_status = 0
          return
        endif
      else
        print *, "file already processed previously"
        lga_status = 1
        return
      endif
    else
      print *, "forecast is too old :", tau_hr, forecast_length
      lga_status = 0
      return
    endif
  elseif(nfiles .eq. 2) then
    ! no, but we do have 2 files that bound this time
    print *, "found bounding files: "
    print *, trim(filenames(2)) ! eariler time : exchange two lines by wei-ting(130326)
    print *, trim(filenames(1)) ! later time   : exchange two lines by wei-ting(130326)
    print *, " "
    ! does latest bound exceed forecast limit?
             ! no...read both in, then time interpolate
                ! get domain info
                ! allocate three sets of arrays (t-1,t,t+1)
                ! get destaggered variables for t-1 and t-1
                ! time interpolate to t
                ! deallocate t-1 and t+1

             ! yes...exit with error
        ! no
          ! exit with error
    
    ! ===== below codes is added to complete        =====
    ! ===== time interpolation by wei-ting (130326) =====
    do iloop = 2,1,-1
       print *, "processing :", trim(filenames(iloop))
       ! make sure file is ready
       call wrfio_wait(filenames(iloop),300)
       call open_wrfnc(filenames(iloop),cdf,istatus)
       if (istatus .ne. 0) then
         print *, "problem opening ",trim(filenames(iloop))
         lga_status = 0
         return 
       endif
       ! get time information
       call get_wrf2_timeinfo(cdf,reftime,dt,itimestep,tau_hr, &
              tau_min, tau_sec, istatus)
       call close_wrfnc(cdf)
       if (istatus .ne. 0) then
         print*, "problem getting time info from ", trim(filenames(iloop))
         lga_status = 0
         return
       endif
       ! is it a new enough forecast?
       if (tau_hr .le. forecast_length) then
         print *, "forecast hour = ",tau_hr

         ! only process if we don't already have an lga output file from wrfarw
         ! that matches this reftime+tau
         print *, "calling check_wrf_lga for:",cmodel,reftime
         call check_wrf_lga(cmodel,reftime,tau_hr,tau_min,already_done) 
         print *, "already done: ",already_done
         if (.not. already_done) then
           print *, "calling wrf2lga"
           call wrf2lga(filenames(iloop:iloop+1),i4time_now,cmodel,istatus)
           if (istatus .ne. 1) then
             print *, "failure in wrf2lga"
             lga_status = 0
             return
           endif
         else
           print *, "file already processed previously"
           lga_status = 1
           ! return
           continue
         endif
       else
         print *, "forecast is too old :", tau_hr, forecast_length
         lga_status = 0
         return
       endif
    enddo
    ! ===== time interpolation =====
    fname9_reftime = wfo_fname13_to_fname9&
                  (reftime(1:4)//reftime(6:7)//reftime(9:13)//reftime(15:16))
    call i4time_fname_lp(fname9_reftime,i4time_reftime,istatus)
    i4time_fcst(1) = files_i4time(1)-i4time_reftime
    i4time_fcst(2) = files_i4time(2)-i4time_reftime
    hr_now = (i4time_now-i4time_reftime)/3600
    min_now = mod(i4time_now-i4time_reftime,3600)/60
    print *, "calling check_wrf_lga for:",cmodel,reftime
    call check_wrf_lga(cmodel,reftime,hr_now,min_now,already_done)
    print *, "already done: ",already_done
    if (.not. already_done) then
       print*,'get 1d pressures'
       call get_pres_1d(i4time_now,nz_laps,pr1d_pa,istatus)
       if(istatus.ne.1)then
          print*,'error returned from get_pres_1d'
          print*,'check pressures.nl or nk_laps in nest7grid.parms'
          stop
       endif
       pr1d_mb(:)=pr1d_pa(:)/100.  ! pa to mb
       print*,i4time_now,i4time_reftime,i4time_reftime, &
              i4time_fcst(2),i4time_fcst(1), &
              files_i4time(2),files_i4time(1)
       ext = 'lga'
       call get_directory(ext,outdir,len_dir)
       print*,outdir,ext,nz_laps
       ! interp 3d fields
       call time_interp(outdir,ext, &
                   nx_laps,ny_laps,nz_laps,6,pr1d_mb,ht_1d, &
                   files_i4time(2),files_i4time(1), &
                   i4time_now,i4time_reftime,i4time_fcst(1), &
                   i4time_reftime,i4time_fcst(2))
       ext = 'lgb'
       call get_directory(ext,outdir,len_dir)
       print*,outdir,ext
       ! interp 2d fields
       call time_interp(outdir,ext, &
                   nx_laps,ny_laps,1,10,pr1d_mb(1),ht_1d(1), &
                   files_i4time(2),files_i4time(1), &
                   i4time_now,i4time_reftime,i4time_fcst(1), &
                   i4time_reftime,i4time_fcst(2))
       
    else
       print *, "file already processed previously"
       lga_status = 1
       return
    endif
    ! ===== end of modification by wei-ting (130326) =====
  else
    print *, "no acceptable wrf files found!"
    lga_status = 0
    return
  endif

  return
end subroutine lga_driver_wrfarw

subroutine  check_wrf_lga(cmodel,reftime,tau_hr,tau_min,already_done)
  implicit none
  include 'netcdf.inc'
  character(len=12),intent(in)         :: cmodel
  character(len=19), intent(in)        :: reftime
  integer,intent(in)                   :: tau_hr,tau_min
  logical,intent(out)                  :: already_done

  character(len=9),external            :: wfo_fname13_to_fname9
  character(len=9)                     :: lga_reftime
  character(len=13)                    :: reftime13
  character(len=200)                   :: lgadir
  character(len=13)                    :: lgatime 
  character(len=200)                   :: lgafile
  character(len=132)                   :: ht_comment
  integer                              :: lendir,cdf,nstatus,varid,cmodel_len
  integer                              :: start(3), count(3)
  already_done = .false.
  print *, "in check_wrf_lga:"
  print *, "cmodel = ", cmodel
  print *, "reftime = ", reftime
  print *, "tau_hr / tau_min = " , tau_hr, tau_min   
  reftime13 = reftime(1:4)//reftime(6:7)//reftime(9:13)//reftime(15:16)
  print *, reftime13
  lga_reftime = wfo_fname13_to_fname9(reftime13)
  print *, reftime13,lga_reftime
  write(lgatime, '(a9,i2.2,i2.2)') lga_reftime,tau_hr,tau_min
  call get_directory('lga',lgadir,lendir)
  lgafile = lgadir(1:lendir)//lgatime//'.lga'
  print *, "lgafile ", lgafile
  inquire(file=trim(lgafile),exist=already_done)
  if (already_done) then

    ! check the comment field to see if this is wrfarw
    nstatus = nf_open(lgafile,nf_nowrite,cdf)
    if (nstatus .ne. nf_noerr) then
      already_done = .false.
      return
    endif
    nstatus = nf_inq_varid(cdf, 'ht_comment', varid)
    if (nstatus .ne. nf_noerr) then
      print *, "problem getting ht_comment"
      nstatus = nf_close(cdf)
      already_done = .false.
      return
    endif
    start = (/1,1,1/)
    count = (/132,1,1/)
    nstatus = nf_get_vara_text(cdf,varid,start,count,ht_comment)
    nstatus = nf_close(cdf)
    if (cmodel(1:6) .eq. ht_comment(1:6)) then
      print *, "file already processed for ",cmodel
      already_done = .true.
    else
      print *, cmodel(1:6),"|",ht_comment(1:6)
      print *, "different cmodel : ", ht_comment(1:6) 
      already_done = .false.
    endif
  endif
  return
end subroutine check_wrf_lga
