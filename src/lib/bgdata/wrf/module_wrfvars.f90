module basicvars

!  contains the arrays and routines needed for allocating
!  and populating the basic variables from wrf.

  use constants
  use wrf_netcdf
  use wrfpost_config
  implicit none

  public
  ! 3d
  real, allocatable        :: p3 (:,:,:)
  real, allocatable        :: t3 (:,:,:)
  real, allocatable        :: z3 (:,:,:)
  real, allocatable        :: z3f(:,:,:)
  real, allocatable        :: qv3(:,:,:)
  real, allocatable        :: rh3(:,:,:)
  real, allocatable        :: u3 (:,:,:)
  real, allocatable        :: v3 (:,:,:)
  real, allocatable        :: w3 (:,:,:)
 
  ! 2d
  real, allocatable        :: p_sfc(:,:)
  real, allocatable        :: mslp(:,:)
  real, allocatable        :: z_1000(:,:)
  real, allocatable        :: tv_sl(:,:)
  real, allocatable        :: t_2m (:,:)
  real, allocatable        :: q_2m (:,:)
  real, allocatable        :: td_2m(:,:)
  real, allocatable        :: rh_2m(:,:)
  real, allocatable        :: u_10m (:,:)
  real, allocatable        :: v_10m (:,:)
  real, allocatable        :: wdir_10m (:,:)
  real, allocatable        :: wspd_10m (:,:)
  real, allocatable        :: qpf_sum(:,:)
  real, allocatable        :: qpf_inc(:,:)
  real, allocatable        :: pcprate(:,:)
  integer, allocatable     :: conv_flag(:,:) 
  real, allocatable        :: xlat(:,:),xlon(:,:),topo(:,:)
  real, allocatable        :: landmask(:,:)
  real, allocatable        :: t_skin(:,:)
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrfpost_get_basic()

  implicit none
  include 'netcdf.inc'
  real, allocatable               :: d1(:,:,:)
  real, allocatable               :: d2(:,:,:)
  real, allocatable               :: d1f(:,:,:)
  real, allocatable               :: d2f(:,:,:)
  real, allocatable               :: dum2d (:,:)
  real, external                  :: relhum, dewpt2, wspd, wdir,mixsat
  integer :: status, i, j, k, rcode
 
  allocate (d1 (nx,ny,nz) )
  allocate (d2 (nx,ny,nz) )
  allocate (d1f (nx,ny,nz+1) )
  allocate (d2f (nx,ny,nz+1) )

  if (.not. allocated(p3)) allocate ( p3 (nx,ny,nz) )
  if (.not. allocated(z3)) allocate ( z3 (nx,ny,nz) )
  if (.not. allocated(z3f)) allocate ( z3f (nx,ny,nz+1) )
  if (.not. allocated(t3)) allocate ( t3 (nx,ny,nz) )
  if (.not. allocated(qv3)) allocate ( qv3(nx,ny,nz) )
  if (.not. allocated(u3))allocate ( u3 (nx,ny,nz) )
  if (.not. allocated(v3)) allocate ( v3 (nx,ny,nz) )
  if (.not. allocated(w3)) allocate ( w3 (nx,ny,nz) )
  if (.not. allocated(p_sfc)) allocate ( p_sfc (nx,ny) )
  if (.not. allocated(t_2m)) allocate ( t_2m(nx,ny) )
  if (.not. allocated(q_2m))allocate ( q_2m (nx,ny) )
  if (.not. allocated(u_10m)) allocate ( u_10m (nx,ny) )
  if (.not. allocated(v_10m)) allocate ( v_10m (nx,ny) )
  if (.not. allocated(qpf_sum)) allocate ( qpf_sum (nx,ny) )
  if (.not. allocated(pcprate)) allocate ( pcprate (nx,ny) )
  if (.not. allocated(xlat)) allocate ( xlat (nx,ny) )
  if (.not. allocated(xlon)) allocate ( xlon (nx,ny) )
  if (.not. allocated(topo)) allocate ( topo (nx,ny) )
  if (.not. allocated(t_skin)) allocate ( t_skin (nx,ny) )
  if (.not. allocated(landmask)) allocate ( landmask (nx,ny) )
  if (.not. allocated(conv_flag)) allocate(conv_flag (nx,ny) )
  conv_flag(:,:) = 0
  ! get pressures
  call get_wrfnc_3d(ncfile,"pb","t",nx,ny,nz,1,d1,status)
  if (status .gt. 0) then
    print *, "+++critical:  could not get base pressure!"
    call abort
  endif
  p3 = d1

  call get_wrfnc_3d(ncfile,"p","t",nx,ny,nz,1,d1,status)
  if (status .gt. 0) then
    print *, "+++critical:  could not get pert pressure!"
    call abort
  endif
  p3 = p3 + d1

  
  ! get heights
  call get_wrfnc_3d(ncfile,"phb","t",nx,ny,nz+1,1,d1f,status)
  if (status.ne.0) then
    print *, 'could not properly obtain wrf base-state geopotential.'
    call abort
  endif
  d2f = d1f
  call get_wrfnc_3d(ncfile,"ph","t",nx,ny,nz+1,1,d1f,status)
  if (status.ne.0) then
    print *, 'could not properly obtain wrf geopotential.'
    call abort
  endif
  d2f = (d2f + d1f) / grav
  z3f = d2f
  do k = 1,nz
    z3(:,:,k) = 0.5 * (z3f(:,:,k) + z3f(:,:,k+1))
  enddo

  ! get theta and convert to temperature
  call get_wrfnc_3d(ncfile, "t","t",nx,ny,nz,1,d1,status)
  if (status.ne.0) then
    print *, 'could not properly obtain wrf perturbation theta.'
    call abort
  endif
  d1 = d1 + 300.
  do k = 1, nz
    do j = 1, ny
      do i = 1,nx
        t3(i,j,k) = d1(i,j,k)/ ((100000./p3(i,j,k))**kappa)
      enddo
    enddo
  enddo
 
   ! get q on sigma
   call get_wrfnc_3d(ncfile, "qvapor","t",nx,ny,nz,1,qv3,status)
   if (status.ne.0) then
     print *, 'could not properly obtain wrf mixing ratio.'
     call abort
   endif

   ! get u on sigma
   call get_wrfnc_3d(ncfile, "u","t",nx,ny,nz,1,u3,status)
   if (status.ne.0) then
     print *, 'could not properly obtain wrf u-comp.'
     call abort
   endif

   ! get u on sigma
   call get_wrfnc_3d(ncfile, "v","t",nx,ny,nz,1,v3,status)
   if (status.ne.0) then
     print *, 'could not properly obtain wrf v-comp.'
     call abort
   endif

   ! get w on sigma
   call get_wrfnc_3d(ncfile, "w","t",nx,ny,nz+1,1,d1f,status)
   if (status.ne.0) then
     print *, 'could not properly obtain wrf w-comp.'
     call abort
   endif
   do k = 1,nz
     w3(:,:,k) = 0.5*(d1f(:,:,k)+d1f(:,:,k+1))
   enddo

  
  deallocate(d1)
  deallocate(d2)
  deallocate(d1f)
  deallocate(d2f)

  ! get the 2d stuff
  
  ! pressure at surface
  p_sfc(:,:) = 0.
  call get_wrfnc_2d(ncfile, "psfc","t",nx,ny,1,p_sfc,status)
  if ((status .ne. 0).or.(maxval(p_sfc) .lt. 10000.))then
    print *, "could not get psfc, using lowest sigma level"
    p_sfc = p3(:,:,1)
  endif
  
  ! 2m temp
  t_2m(:,:) = 0.
  call get_wrfnc_2d(ncfile, "t2","t",nx,ny,1,t_2m,status)
  if ((status .ne. 0).or.(maxval(t_2m) .lt. 100.))then
    print *, "could not get t2, using lowest sigma level"
    t_2m = t3(:,:,1)
  endif
  
  ! qvapor at 2m
  q_2m(:,:) = 0.
  call get_wrfnc_2d(ncfile, "q2","t",nx,ny,1,q_2m,status)
  if ((status .ne. 0).or.(maxval(q_2m) .lt. 0.0001))then
    print *, "could not get q2, using lowest sigma level"
    q_2m = qv3(:,:,1)
  endif
  ! because 2m qv and t are derived from the pbl scheme and
  ! the wrf is apparently not checking for saturation, clean this
  ! up now
  do j = 1, ny
    do i= 1, nx
       q_2m(i,j) = min(q_2m(i,j),mixsat(t_2m(i,j),p_sfc(i,j)))
    enddo
  enddo
  
  ! 10m wind
  !u_10m(:,:) = 0. 
  !call get_wrfnc_2d(ncfile, "u10","t",nx,ny,1,u_10m,status) 
  !if ((status .ne. 0).or.(maxval(u_10m) .eq. 0.))then
   ! print *, "could not get u10, using lowest sigma level"
    u_10m = u3(:,:,1)
  !endif
  
  !v_10m(:,:) = 0.
  !call get_wrfnc_2d(ncfile, "v10","t",nx,ny,1,v_10m,status)
  !if ((status .ne. 0).or.(maxval(v_10m) .eq. 0.))then
  !  print *, "could not get v10, using lowest sigma level"
    v_10m = v3(:,:,1)
  !endif  
  
  ! precip
  allocate(dum2d (nx,ny))
  qpf_sum(:,:) = 0.
  call get_wrfnc_2d(ncfile, "rainnc","t",nx,ny,1,dum2d,status)
  if (status .ne. 0)then
    print *, "could not get rainnc, assuming 0."
  else
    qpf_sum = dum2d
    call get_wrfnc_2d(ncfile, "rainc","t",nx,ny,1,dum2d,status)
    if (status .ne. 0)then
       print *, "could not get rainc, assuming 0."
    else
      ! set convective flag for any point that received parameterized precip
       do j = 1, ny
         do i = 1, nx
           if(dum2d(i,j) .gt. 0) conv_flag(i,j) = 1
            ! zero out small values
           !if (dum2d(i,j) .lt. 1e-4) dum2d(i,j) = 0
          enddo
        enddo
       qpf_sum = qpf_sum + dum2d
    endif
  endif
  where(qpf_sum .lt. 1.e-4) qpf_sum = 0.
  ! lets get incremental precipitation
  dum2d(:,:) = 0.
  if (.not. allocated (qpf_inc)) allocate(qpf_inc(nx,ny))
  if (have_prev) then
    if (pretime_str .ne. reftime_str) then
      ! there is a previous output time available
      ! and it is not the 0-h tau, so let's compute
      ! get the qpf from the previous time and subtract it out
      call open_wrfnc(filename_prev, ncfile_prev, status)
      if (status .eq. 0) then
        call get_wrfnc_2d(ncfile_prev, "rainnc","t",nx,ny,1,dum2d,status)
	if (status .ne. 0) dum2d(:,:) = 0.
	qpf_inc = dum2d
	call get_wrfnc_2d(ncfile_prev, "rainc","t",nx,ny,1,dum2d,status)
	if (status .eq. 0) then 
          ! if we had convective precip at any grid point in the previous time
          ! go ahead and set the flag for convection
          do j = 1, ny
            do i = 1, nx
              if(dum2d(i,j) .gt. 0) conv_flag(i,j) = 1
              ! zero out small values
              !if (dum2d(i,j) .lt. 1e-4) dum2d(i,j)= 0.
            enddo
          enddo
          qpf_inc = qpf_inc + dum2d
        endif
        where(qpf_inc .lt. 1.e-4) qpf_inc = 0.
        qpf_inc = qpf_sum - qpf_inc
        call close_wrfnc(ncfile_prev)
      else 
        print *, "error opening previous file"
        qpf_inc = 0.
      endif
    else
      qpf_inc = qpf_sum 
    endif
  else
    qpf_inc(:,:) = 0.
  endif
  
  print *, "   min/max qpf total: ", minval(qpf_sum),maxval(qpf_sum)
  print *, "   min/max qpf incre: ", minval(qpf_inc),maxval(qpf_inc)

  ! preciptation rate (get the timestep resolved and convective
  ! precip for this period, divide by the timestep (seconds) and
  ! multiply by 3600 to get mm per hour.  if we cannot obtain
  ! it, then deallocate it.
  print *, '  -- precip rate'
  allocate(dum2d(nx,ny))
  call get_wrfnc_2d(ncfile, "rainncv", "t",nx,ny,1,dum2d,status)
  if (status.gt.0) then
    print *, "did not get rainncv"
    pcprate = 0.
  else
    where(dum2d .lt. 1.e-07) dum2d = 0.
    pcprate = dum2d
    call get_wrfnc_2d(ncfile, "raincv", "t",nx,ny,1,dum2d,status)
    if (status .eq. 0) then
      where(dum2d .lt. 1.e-06) dum2d =0.
      pcprate = pcprate + dum2d
      do j = 1, ny
        do i = 1, nx
          if(dum2d(i,j) .gt. 0) conv_flag(i,j) = 1
        enddo
      enddo
    endif
  endif
  rcode = nf_get_att_real(ncfile,0,"dt",dt)
  if (rcode .eq. nf_noerr) then
    pcprate = pcprate / dt   ! converts to mm per sec
    where (pcprate .lt. 1.e-8) pcprate = 0.
  else
    print *, "could not get timestep to compute pcprate!"
    deallocate(pcprate)
  endif
  deallocate(dum2d)

  call get_wrfnc_2d(ncfile, "xlat","t",nx,ny,1,xlat,status)
  call get_wrfnc_2d(ncfile, "xlong","t",nx,ny,1,xlon,status)
  call get_wrfnc_2d(ncfile, "hgt","t",nx,ny,1,topo,status)
  call get_wrfnc_2d(ncfile, "landmask","t",nx,ny,1,landmask,status)
  call get_wrfnc_2d(ncfile, "tsk", "t", nx,ny,1,t_skin,status)


  ! compute a few more things
  ! rh and dewpoint at surface are nice to have
  print *, "  computing td_2m and rh_2m"
  allocate(rh_2m(nx,ny))
  allocate(td_2m(nx,ny))
  allocate(rh3(nx,ny,nz))
 
  do k= 1, nz
    do j = 1, ny
      do i = 1, nx 
        rh3(i,j,k) = relhum(t3(i,j,k),qv3(i,j,k),p3(i,j,k)) * 100.
      enddo
    enddo
  enddo
  do j = 1, ny
    do i = 1, nx
      rh_2m(i,j) = relhum(t_2m(i,j),q_2m(i,j),p_sfc(i,j))
      if (rh_2m(i,j) .gt. 1.0) then
        print *,"i/j/t/p/q/rh = ",i,j,t_2m(i,j),p_sfc(i,j),q_2m(i,j),rh_2m(i,j)
      endif
      td_2m(i,j) = dewpt2(t_2m(i,j),rh_2m(i,j))
      rh_2m(i,j) = rh_2m(i,j) * 100.
    enddo
  enddo
  print *, '    min/max td_2m = ', minval(td_2m),maxval(td_2m)
  print *, '    min/max rh_2m = ', minval(rh_2m),maxval(rh_2m)

  ! lets also get a surface true wind direction and speed
  allocate(wdir_10m(nx,ny))
  allocate(wspd_10m(nx,ny))
  do j = 1, ny
    do i = 1, nx
      wdir_10m(i,j) = wdir(u_10m(i,j),v_10m(i,j),xlon(i,j),grid%stdlon, &
                     grid%cone)
      wspd_10m(i,j) = wspd(u_10m(i,j),v_10m(i,j))
    enddo
  enddo

  ! mean sea-level pressure
  allocate(mslp(nx,ny))
  allocate(tv_sl(nx,ny))
  allocate(z_1000(nx,ny))
  print *, "wrfpost_basicvars:  calling derive_mslp"
  call derive_mslp(nx,ny,nz,t3,qv3,p3,p_sfc,topo,mslp,tv_sl,z_1000)

  return

end subroutine wrfpost_get_basic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module basicvars
