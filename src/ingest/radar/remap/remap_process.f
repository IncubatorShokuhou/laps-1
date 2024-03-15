cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis
      subroutine remap_process(
     :         i_tilt,                                     ! integer (input)
     :         i_last_scan,                                ! integer (input)
     :         i_first_scan,                               ! integer (input)
     :         grid_rvel,grid_rvel_sq,grid_nyq,ngrids_vel,n_pot_vel, ! (output)
     :         grid_ref,ngrids_ref,n_pot_ref,                        ! (output)
     :         nx_l,ny_l,nz_l,nx_r,ny_r,                   ! integer   (input)
     :         l_offset_radar,ioffset,joffset,             !           (input)
     1         lat,lon,topo,                               !           (input)
     1         i_scan_mode,                                !           (input)
     :         slant_ranges_m,                             !           (input)
     :         n_rays,                                     !           (input)
     :         n_gates,                                    !           (input)
     1         velocity,reflect,                           !           (input)
     1         az_array,max_ray_tilt,elevation_deg,        !           (input)
     1         vel_nyquist,                                !           (input)
     :         ref_min,min_ref_samples,min_vel_samples,dgr,! integer (input)
     :         laps_radar_ext,c3_radar_subdir,             ! char      (input)
     :         path_to_vrc,                                ! char      (input)
     :         namelist_parms,                             ! struct    (input)
     :         i_product_i4time,                           ! integer (input)
     :         i_num_finished_products,                    ! integer (output)
     :         i_status_tilt,i_status)                     ! integer (output)
c
c     subroutine remap_process
c
c     purpose:
c       main process routine for the remapping algorithm
c
c **************************** history section ****************************
c
c       windsor, c. r.  10-jul-1985     original version
c       albers, steve    7-apr-1986     update for version 03
c       albers, steve      jun-1987     map velocities to laps grid
c       albers, steve      mar-1988     streamlined and converted to 87 data
c       albers, steve      dec-1988     further conversions for rt87 cartesian
c       albers, steve      may-1992     turn off range unfolding for velocities
c       albers, steve      feb-1993     min#, 40% frac qc for reflectivity added
c       albers, steve      may-1994     88d version for sun risc box
c       brewster, keith    aug-1994     clean-out of artifacts
c       brewster, keith    apr-1995     added initial_gate parameter
c                                       modified volume nyquist determination.
c       brewster, keith    sep-1995     added point-by-point nyquist calc.
c       albers, steve   19-dec-1995     changed gate_spacing_m variable to 
c                                       gate_spacing_m_ret to prevent 
c                                       reassigning a parameter value. the 
c                                       location is the call to read_data_88d 
c                                       and a new declaration.
c                                       environment variable evaluations added
c                                       for ftping and purging the output.
c                                       new streamlined purging function.
c       albers, steve      feb-1996     linear reflectivity averaging (via lut)
c       albers, steve      may-1996     new igate_lut to reduce processing 
c       albers, steve          1998     more flexibility added

*********************** declaration section **************************
c
      include 'trigd.inc'
      implicit none
c
c     input variables
c
      integer i_tilt, i_status_tilt
      integer i_last_scan
      integer i_first_scan
      integer i_product_i4time
      integer max_ray_tilt
      integer nx_l,ny_l,nz_l,nx_r,ny_r
      integer ioffset,joffset,io,jo,iomin,iomax,jomin,jomax

      integer min_ref_samples,min_vel_samples

      real ref_min,dgr
c
c     laps grid dimensions
c
      include 'remap_constants.dat'
      include 'remap.cmn'
c
c     velocity obs
c
      real grid_rvel(nx_r,ny_r,nz_l)  !  radial radar velocities
      real grid_rvel_sq(nx_r,ny_r,nz_l)
      real grid_nyq(nx_r,ny_r,nz_l)
      integer ngrids_vel(nx_r,ny_r,nz_l)
      integer n_pot_vel(nx_r,ny_r,nz_l)
c
c     reflectivity obs
c
      real grid_ref (nx_r,ny_r,nz_l)  !  radar reflectivities
      integer ngrids_ref (nx_r,ny_r,nz_l)
      integer n_pot_ref (nx_r,ny_r,nz_l)
c
c     output variables
c
      integer i_num_finished_products
      integer i_status
c
c     processing parameters
c
      real re43
      parameter (re43 = 8503700.) ! 4/3 radius of the earth in meters
c
c
      integer max_fields
      parameter (max_fields = 10)
c
c     variables for netcdf i/o
c
      character*150 dir

      character*9 a9time
      character*31 ext,ext_in
      character*3 var_a(max_fields)
      character*125 comment_a(max_fields)
      character*10  units_a(max_fields)
      character*4 laps_radar_ext
      character*3 c3_radar_subdir
      character*(*) path_to_vrc

c
c     functions
c
      real height_to_zcoord
      real height_of_level
c
c     misc local variables
c
      integer igate,i_scan_mode,jray,end_ext,ilut_ref
c
      real  slant_ranges_m (max_gates),
     :        elevation_deg,
     :        az_array(max_ray_tilt),
     :        velocity(max_gates,max_ray_tilt),
     :        reflect(max_gates,max_ray_tilt)

      real, allocatable, dimension(:,:,:,:) :: out_array_4d

      real r_missing_data

      real lat(nx_l,ny_l)      
      real lon(nx_l,ny_l)      
      real topo(nx_l,ny_l)     
      real dum_2d(nx_l,ny_l)   ! local
      integer k_eff(nx_l,ny_l) ! local
c
      logical l_unfold, l_compress_output, l_domain_read, l_offset_radar       
      save l_domain_read
      data l_domain_read /.false./
c 
      real avgvel,vel_nyquist,vel_value,ref_value,lat_dum,lon_dum
      real v_nyquist_tilt(max_tilts)
      real v_nyquist_vol
      real gate_spacing_m_ret,grid_spacing_cen_m
      real height_grid,range_dum,range_new,azimuth,elevation_dum
      real height_guess
c
      integer i,j,k,k_low,ielev,igate_lut,iter,klut,idebug
      integer nazi,iran
      integer num_sweeps,n_rays,n_gates,n_obs_vel,n_output_data,nf
      integer igate_max
      integer igate_interval
      integer n_vel_grids_final,n_vel_grids_prelim
      integer n_ref_grids,n_ref_grids_qc_fail,nycor
      integer istatus,istatus_qc,istat_alloc,len_ext
      integer ishow_timer,i4_elapsed
      integer i_purge
      integer init_ref_gate_hyb,init_ref_gate_actual
      integer mingate_valid_ref,maxgate_valid_ref

      real rvel,azimuth_interval
      real rmax,height_max,rlat_radar,rlon_radar,rheight_radar
      real vknt,rknt,variance,hybrid_range

      character*4 c4_radarname
      character*7 c7_laps_xmit
      character*7 c7_laps_purge
      character*7 c7_laps_sleep

      save height_max, k_low, n_obs_vel, n_output_data
c
c     beginning of executable code
c
      idebug = 0 ! more verbose output (0-2 range)

      write(6,*)
      write(6,805) i_first_scan,i_last_scan,i_tilt
  805 format(' remap_process > ifirst,ilast,tilt',4i5)

      call get_l_compress_radar(l_compress_output,i_status)
      if(i_status .ne. 1)then
          write(6,*)' error in obtaining l_compress_radar'
          return
      endif

      rlat_radar = rlat_radar_cmn
      rlon_radar = rlon_radar_cmn
      rheight_radar = rheight_radar_cmn
      c4_radarname = c4_radarname_cmn

      if(rlat_radar .eq. 0.)then
          write(6,*)' error: no radar coords in remap_process'
          i_status = 0
          return
      endif

      i_num_finished_products = 0

      call get_r_missing_data(r_missing_data, i_status)
      if(i_status .ne. 1)then
          write(6,*)' error in obtaining r_missing_data'
          return
      endif
c
c     for first scan, initialize sums and counters to zero.
c
      if (i_first_scan .eq. 1 .or. i_first_scan .eq. 999) then

        i4_elapsed = ishow_timer()

        write(6,806)
  806   format
     1  (' remap_process > 1st sweep - initializing vel/ref arrays')

        n_obs_vel = 0
        n_output_data = 0

        grid_rvel(:,:,:) = 0.
        grid_rvel_sq(:,:,:) = 0.
        grid_nyq(:,:,:) = 0.
        grid_ref(:,:,:) = 0.
        ngrids_vel(:,:,:) = 0
        ngrids_ref(:,:,:) = 0
        n_pot_vel(:,:,:) = 0
        n_pot_ref(:,:,:) = 0
c
c       compute maximum height of data needed.
c
        height_max = height_of_level(nz_l)
c
c       define lower limit of radar coverage in laps grid
c
        k_low = int(height_to_zcoord(rheight_radar,i_status))
        k_low = max(k_low,1)

!       true will map one tilt in each laps level (for testing only)
        if(namelist_parms%l_ppi_mode)then 
            k_low = 1 
        endif

      end if ! initialize for 1st scan

c     former location of 'read_data_88d' call
      if (i_status_tilt .ne. 1) go to 998 ! abnormal return
c
      v_nyquist_tilt(i_tilt) = vel_nyquist
c
      write(6,*)' remap_process > vel_nyquist for this tilt = '
     :        ,i_tilt,vel_nyquist
c
c     compute max range from elevation angle
c
      rmax = -re43 * sind(elevation_deg)
     :  + sqrt(re43*re43*sind(elevation_deg)**2 +
     :          height_max * (2.*re43 + height_max))


      print *, ' rmax,height_max= ',rmax,height_max

      print *, ' gate_spacing_m,gate_interval= ',gate_spacing_m,
     :           gate_interval

      igate_max = min(int(rmax/gate_spacing_m) , n_gates)

      write(6,809) i_scan_mode,n_gates,igate_max,elevation_deg
 809  format
     :(' remap_process > i_scan_mode,n_gates,igate_max,elevation = '
     :                                                 ,i3,2i5,f5.1)

!     calculate effective range/height at grid point centers (done iteratively)
      
!     first read domain grid info if needed
!     if(.not. l_domain_read)then
      if(.false.)then
          write(6,*)' remap_process > call get_laps_domain_95'
          call get_laps_domain_95(nx_l,ny_l,lat,lon,topo
     1                           ,dum_2d,grid_spacing_cen_m
     1                           ,istatus)
          if(istatus .ne. 1)then
              write(6,*)' error return from get_laps_domain_95'
              return
          endif
          l_domain_read = .true.
      endif

      i4_elapsed = ishow_timer()

      write(6,*)' remap_process > calculate k_eff array'

      height_guess = height_max

      do j = 1,ny_l
      do i = 1,nx_l
!         use guessed height to obtain initial range value
          call latlon_to_radar(lat(i,j),lon(i,j),height_guess           ! i
     1                        ,azimuth,range_new,elevation_dum          ! o
     1                        ,rlat_radar,rlon_radar,rheight_radar)     ! i  

          iter = 0

          if(range_new .le. rmax)then

            range_dum = r_missing_data

!           converge on height where radar beam hits the grid point
            do while (abs(range_dum-range_new).gt.10. .and. iter.lt.10)          
              range_dum = range_new

!             obtain height where radar beam hits grid point
              call radar_to_latlon(lat_dum,lon_dum,height_grid          ! o
     1                            ,azimuth,range_dum,elevation_deg      ! i
     1                            ,rlat_radar,rlon_radar,rheight_radar) ! i

!             use corrected height to obtain more accurate range
              call latlon_to_radar(lat(i,j),lon(i,j),height_grid        ! i
     1                            ,azimuth,range_new,elevation_dum      ! o
     1                            ,rlat_radar,rlon_radar,rheight_radar) ! i  

              iter = iter + 1

              height_guess = height_grid

            enddo 

            k_eff(i,j) = nint(height_to_zcoord(height_grid,i_status))

          else
            k_eff(i,j) = 0

          endif

      enddo
      enddo

      i4_elapsed = ishow_timer()
c
c     find elevation index in look up table
c
      ielev = nint(((elevation_deg-min_elev) * lut_elevs)
     1                  /(max_elev-min_elev)                 )
      ielev = max(ielev,0)
      ielev = min(ielev,lut_elevs)
      write(6,*)' remap_process > elev index = ',ielev

      i4_elapsed = ishow_timer()

      write(6,*)' remap_process > looping through rays and gates'

      azimuth_interval = 360. / float(lut_azimuths)
      write(6,*)' azimuth_interval = ',azimuth_interval

      write(6,*)' first azimuth = ',az_array(1)

      mingate_valid_ref = 99999
      maxgate_valid_ref = 0

      do 200 jray=1, n_rays

        if(az_array(jray) .ne. r_missing_data)then
            nazi = nint(az_array(jray) / azimuth_interval)
            nazi = mod(nazi,lut_azimuths)
        else
            goto200
        endif

        if(namelist_parms%l_hybrid_first_gate)then
          if(elevation_deg .lt. 1.0)then     ! < 1.0
              hybrid_range = 100000. ! 80000.
          elseif(elevation_deg .lt. 2.0)then ! between 1.0 and 2.0
              hybrid_range = 60000. ! 60000.
          elseif(elevation_deg .lt. 3.0)then ! between 2.0 and 3.0
              hybrid_range = 28000. ! 30000.
          else
              hybrid_range = 0.
          endif

          init_ref_gate_hyb = hybrid_range / gate_spacing_m

          if(abs(az_array(jray) - 250.) .lt. 80.)then ! higher terrain
              init_ref_gate_actual = max(initial_ref_gate
     1                                  ,init_ref_gate_hyb)       
          else
              init_ref_gate_actual = initial_ref_gate
          endif

!         write(6,*)
!    1        ' l_hybrid_first_gate flag is set, first range/gate = '      
!    1        ,hybrid_range,init_ref_gate_actual

        else
          init_ref_gate_actual = initial_ref_gate

        endif

        igate_interval=1

        do 180 igate=init_ref_gate_actual,igate_max,igate_interval       

          if(lgate_lut(igate))then ! we'll process this gate, it may have data

            igate_lut = igate/gate_interval

            iran = gate_elev_to_projran_lut(igate_lut,ielev)

            i = azran_to_igrid_lut(nazi,iran)
            j = azran_to_jgrid_lut(nazi,iran)

            if (i .eq. 0 .or. j .eq. 0) go to 180

            io = i - ioffset
            jo = j - joffset 

            if(k_eff(i,j) .ne. 0)then
                k = k_eff(i,j)
                klut = gate_elev_to_z_lut(igate_lut,ielev)
                if(jray .eq. 1 .and. idebug .ge. 1)then
                    write(6,5)igate,igate_lut,iran,i,j,k,klut
5                   format(' igate,igate_lut,iran,i,j,k,klut= ',7i6)
                endif
            else
                k = gate_elev_to_z_lut(igate_lut,ielev)
            endif

            if (k .eq. 0) go to 180

            if (k .lt. k_low .and. elevation_deg .ge. 0.)then
                write(6,*)' error: inconsistent k values - ',k,k_low
                stop
            endif

!           true will map one tilt in each laps level (for testing only)
            if(namelist_parms%l_ppi_mode)then 
                k = i_tilt 
            endif

            if( lgate_vel_lut(igate) ) then

!           if( igate .ge. initial_vel_gate) then
c
c      velocity data
c
              n_pot_vel(io,jo,k) = n_pot_vel(io,jo,k) + 1
c
c      map velocity if data present and abs value of velocity is
c      more than 2 ms-1.
c
              vel_value = velocity(igate,jray)

              if (abs(vel_value) .lt. vel_mis_check .and.
     :            abs(vel_value) .gt. namelist_parms%abs_vel_min ) then

                if(ngrids_vel(io,jo,k).eq.0 .or. 
     1             vel_nyquist .eq. r_missing_data) then
                  rvel =  vel_value

                elseif(vel_nyquist .ne. r_missing_data) then ! .and. ngrids > 0
                  avgvel=grid_rvel(io,jo,k)/float(ngrids_vel(io,jo,k))
                  nycor=nint(0.5*(avgvel-vel_value)/
     :                     vel_nyquist)
                  rvel=vel_value+((2*nycor)*vel_nyquist)

                end if

                n_obs_vel = n_obs_vel + 1
!               if(n_obs_vel .le. 100)then
!                   write(6,*)'n_obs_vel = ',n_obs_vel,i,j,k
!               endif
                grid_rvel(io,jo,k) = grid_rvel(io,jo,k) + rvel
                grid_rvel_sq(io,jo,k) =
     :          grid_rvel_sq(io,jo,k) + rvel*rvel

                ngrids_vel(io,jo,k) = ngrids_vel(io,jo,k) + 1

                if(vel_nyquist .ne. r_missing_data)then
                    grid_nyq(io,jo,k)=grid_nyq(io,jo,k)+vel_nyquist
                endif

              end if

            end if
c
c     map reflectivity
c
            if( lgate_ref_lut(igate) ) then

              n_pot_ref(io,jo,k) = n_pot_ref(io,jo,k) + 1

              ref_value = reflect(igate,jray)

              if (abs(ref_value) .lt. ref_mis_check) then ! datum is present

c               grid_ref(io,jo,k) =
c    :          grid_ref(io,jo,k) + ref_value

                if(jray .eq. 1 .and. idebug .ge. 2)then
                    write(6,170)ref_value
 170                format(45x,'ref = ',f8.2)
                endif

                ilut_ref = nint(ref_value * 10.) ! tenths of a dbz
                grid_ref(io,jo,k) =
     :          grid_ref(io,jo,k) + dbz_to_z_lut(ilut_ref)
                ngrids_ref(io,jo,k) = ngrids_ref(io,jo,k) + 1

                mingate_valid_ref = min(mingate_valid_ref,igate)
                maxgate_valid_ref = max(maxgate_valid_ref,igate)

              end if

            endif ! l_gate_ref(igate) = .true. and we process the reflectivity

          endif ! l_gate(igate) = .true. and we need to process this gate

  180   continue ! igate
  200 continue ! jray

      write(6,815,err=816)elevation_deg,n_obs_vel
     1                   ,mingate_valid_ref,maxgate_valid_ref
  815 format(' remap_process > end ray/gate loop: elev= ',f10.2
     :      ,'  n_obs_vel = ',i9,' min/max ref gates ',2i7)

  816 i4_elapsed = ishow_timer()

      if (i_last_scan .eq. 1) then
        write(6,820)
  820   format(
     :  ' remap_process > last sweep - dividing velocity & ref arrays')       


        n_vel_grids_prelim = 0
        n_vel_grids_final = 0
        n_ref_grids = 0
        n_ref_grids_qc_fail = 0

c
c     diagnostic print-out
c
        write(6,825)
  825   format(' remap_process > prepare reflectivity output')

        do 480 k = 1, k_low-1
        do 480 jo = 1, ny_r
        do 480 io = 1, nx_r
          grid_ref(io,jo,k)=r_missing_data
          grid_rvel(io,jo,k)=r_missing_data
          grid_nyq(io,jo,k)=r_missing_data
  480   continue

        do 500 k = k_low,nz_l

          write(6,826) k
  826     format(' remap_process > dividing: k = ',i2)

          do 400 jo = 1,ny_r
          do 400 io = 1,nx_r
c
c     note min_vel_samples must be greater than 1
c
            if(ngrids_vel(io,jo,k) .ge. min_vel_samples) then ! good gates
              vknt=float(ngrids_vel(io,jo,k))

              if (vknt .ge. float(n_pot_vel(io,jo,k))
     1                                          * coverage_min_vel) then       

                n_vel_grids_prelim = n_vel_grids_prelim + 1
                variance=( grid_rvel_sq(io,jo,k) - 
     :                    (grid_rvel(io,jo,k)*grid_rvel(io,jo,k)/vknt) )
     :                     /(vknt-1.)

                if (variance .lt. rv_var_lim) then ! increment good counter

                  n_vel_grids_final = n_vel_grids_final + 1
                  grid_rvel(io,jo,k) = grid_rvel(io,jo,k)/vknt
                  if(vel_nyquist .ne. r_missing_data)then
                      grid_nyq(io,jo,k) = grid_nyq(io,jo,k)/vknt
                  else
                      grid_nyq(io,jo,k) = r_missing_data
                  endif

                else ! failed vel qc test

                  grid_rvel(io,jo,k) = r_missing_data
                  grid_nyq(io,jo,k) = r_missing_data
    
                end if ! vel qc test

              else ! insufficient coverage

                grid_rvel(io,jo,k) = r_missing_data
                grid_nyq(io,jo,k) = r_missing_data

              end if ! velocity coverage check

            else ! insufficient velocity count

              grid_rvel(io,jo,k) = r_missing_data
              grid_nyq(io,jo,k) = r_missing_data

            end if ! first check of velocity count
c
c     reflectivity data
c
!        qc flags within the gridded reflectivities are as follows...
!
!            r_missing_data     insufficient number of "potential gates" in
!                               the grid volume
!
!            -101.              insufficient fractional coverage of "actual"
!                               gates in grid volume
!
!            -102.              reflectivity less than threshold value


            if(ngrids_ref(io,jo,k) .ge. min_ref_samples) then ! good gates
              rknt=float(ngrids_ref(io,jo,k))
              if (rknt .ge. float(n_pot_ref(io,jo,k)) 
     1                                          * coverage_min_ref) then       

!               calculate mean value of z
                grid_ref(io,jo,k) = grid_ref(io,jo,k)/rknt

!               convert from z to dbz
                grid_ref(io,jo,k) = alog10(grid_ref(io,jo,k)) * 10.

                if (grid_ref(io,jo,k) .ge. ref_min) then

                  n_ref_grids = n_ref_grids + 1
                  if(n_ref_grids .lt. 200 .and. idebug .ge. 1)
     :               write(6,835) io,jo,k,grid_ref(io,jo,k)
  835                format(' grid loc: ',3(i4,','),'  refl: ',f6.1)

                else        ! failed ref qc test
 
                  n_ref_grids_qc_fail = n_ref_grids_qc_fail + 1
                  grid_ref(io,jo,k) = -102.

                end if      ! passed ref qc test

              else       ! insufficent coverage

                grid_ref(io,jo,k) = -101.

              end if     ! coverage check of count

            else       ! insufficent data count

              grid_ref(io,jo,k) = r_missing_data

            end if     ! first check of count

  400     continue ! io,jo
  500   continue ! k

        i4_elapsed = ishow_timer()
c
c     call qc routine (now disabled)
c
        istatus_qc = 1
c       call radar_qc(nx_l,ny_l,nz_l,grid_rvel,istatus_qc)
        if (istatus_qc .ne. 1) then
          i_num_finished_products = 0
          write(6,840)
  840     format(' remap_process > bad data detected, no data written')       
          go to 998 ! abnormal return
        end if

        write(6,842) n_ref_grids_qc_fail,n_ref_grids
  842   format(' remap_process > n_ref_qc_fail/n_ref = ',2i12)

        if (n_ref_grids .lt. ref_grids_check) then
          i_num_finished_products = 0
          write(6,845) n_ref_grids,ref_grids_check
  845     format(' remap_process > ',i4,' ref grids < ',i4
     :                                 ,'no data file written...')
          go to 999 ! normal return
        end if

        write(6,851)n_ref_obs_old(1),n_ref_grids,i4time_old(1)
     1                                          ,i_product_i4time

  851   format(' remap_process > ref obs: old/new',2i6
     :        ,' i4time: old/new',2i11)

        i4time_old(1) = i_product_i4time
        n_ref_obs_old(1) = n_ref_grids
c
c     determine filename extension
        ext = laps_radar_ext
        write(6,*)' remap_process > laps_ext = ',laps_radar_ext
c
c     prepare to write out data
c
        i4_elapsed = ishow_timer()
c
        var_a(1) = 'ref'
        var_a(2) = 'vel'
        var_a(3) = 'nyq'
        units_a(1) = 'dbz'
        units_a(2) = 'm/s'
        units_a(3) = 'm/s'
        comment_a(1) = 'doppler reflectivity'
        comment_a(2) = 'doppler velocity'
        comment_a(3) = 'nyquist velocity'
        nf = 3

c       do 555 k=7,9
c       print *, 'sample data on level ',k
c       do 555 jo=1,ny_r
c       do 555 io=60,60
c         print *,io,jo,grid_ref(io,jo,k),grid_rvel(io,jo,k)
c 555   continue
c

        v_nyquist_vol = -999.
        write(6,875) i_tilt
  875   format(' determine v_nyquist for the ',i4,' tilt volume')
c
        do 600 i = 1,i_tilt
          write(6,880) i,v_nyquist_tilt(i)
  880     format(' i_tilt:',i6,'  v_nyquist_tilt:',e12.4)
          if (v_nyquist_tilt(i) .gt. 0.) then
            if (v_nyquist_vol .gt. 0.) then
              if (v_nyquist_tilt(i) .ne. v_nyquist_vol) then
                v_nyquist_vol = r_missing_data
                write(6,886)
  886           format(' nyquist has changed for the tilt',
     1                 ', set v_nyquist_vol to missing.')
                go to 601
              end if
            else
              v_nyquist_vol = v_nyquist_tilt(i)
            end if
          end if
  600   continue
  601   continue
c
c     write out header type info into the comment array
c
        l_unfold=.false.
        write(comment_a(1),888)rlat_radar,rlon_radar,rheight_radar
     1       ,n_ref_grids,c4_radarname
        write(comment_a(2),889)rlat_radar,rlon_radar,rheight_radar
     1       ,n_vel_grids_final,c4_radarname,v_nyquist_vol,l_unfold
        write(comment_a(3),888)rlat_radar,rlon_radar,rheight_radar
     1       ,n_vel_grids_final,c4_radarname

  888   format(2f9.3,f8.0,i7,a4,3x)
  889   format(2f9.3,f8.0,i7,a4,3x,e12.4,l2)

        write(6,890)comment_a(1)(1:80)
        write(6,890)comment_a(2)(1:80)
        write(6,890)comment_a(3)(1:80)
  890   format(a80)

        i4_elapsed = ishow_timer()

        if(laps_radar_ext(1:3) .ne. 'vrc')then ! vxx output

            i4_elapsed = ishow_timer()

            call ref_fill_horz(grid_ref,nx_l,ny_l,nz_l
     1                        ,lat,lon,dgr
     1                        ,nx_r,ny_r,ioffset,joffset
     1                        ,rlat_radar,rlon_radar,rheight_radar
     1                        ,istatus)       
            if(istatus .ne. 1)then
                write(6,*)' error calling ref_fill_horz'          
                return
            endif

            allocate( out_array_4d(nx_l,ny_l,nz_l,3), stat=istat_alloc )       
            if(istat_alloc .ne. 0)then
                write(6,*)' error: could not allocate out_array_4d'
                stop
            endif

            i4_elapsed = ishow_timer()

            write(6,*)' filling output arrays'

            if(l_offset_radar)then

                out_array_4d = r_missing_data ! initialize

                iomin = max((1-ioffset),1)
                iomax = min((nx_l-ioffset),nx_r)

                jomin = max((1-joffset),1)
                jomax = min((ny_l-joffset),ny_r)

                write(6,*)' io range: ',iomin,iomax,
     1                    ' jo range: ',jomin,jomax

                do k = 1,nz_l
                    do jo = jomin,jomax ! 1,ny_r
                        j = jo + joffset
                        do io = iomin,iomax ! 1,nx_r
                            i = io + ioffset
                            out_array_4d(i,j,k,1) = grid_ref(io,jo,k)
                            out_array_4d(i,j,k,2) = grid_rvel(io,jo,k)       
                            out_array_4d(i,j,k,3) = grid_nyq(io,jo,k)
                        enddo ! io
                    enddo ! jo
                enddo ! k

            else
                out_array_4d(:,:,:,1) = grid_ref(:,:,:)
                out_array_4d(:,:,:,2) = grid_rvel(:,:,:)
                out_array_4d(:,:,:,3) = grid_nyq(:,:,:)

            endif

            i4_elapsed = ishow_timer()

            call make_fnam_lp(i_product_i4time,a9time,istatus)
            if(istatus .ne. 1)return

            call s_len(ext,len_ext)

            if(l_compress_output)then
                write(6,865) c4_radarname,ext(1:len_ext),a9time
 865            format(
     1             ' remap_process > calling put_compressed_multi_3d'          
     1             ,1x,a4,2x,a,2x,a9)          

                call put_compressed_multi_3d(i_product_i4time,ext,var_a       
     1                                ,units_a,comment_a,out_array_4d
     1                                ,nx_l,ny_l,nz_l,nf,istatus)

            else
                write(6,866) c4_radarname,ext(1:len_ext),a9time
 866            format(' remap_process > calling put_laps_multi_3d'
     1                 ,1x,a4,2x,a,2x,a9)          

                call put_laps_multi_3d(i_product_i4time,ext,var_a
     1                                ,units_a,comment_a,out_array_4d
     1                                ,nx_l,ny_l,nz_l,nf,istatus)

            endif

            deallocate(out_array_4d)

        else ! single level of data (as per wfo)
            call put_remap_vrc(i_product_i4time,comment_a(1)
     1             ,rlat_radar,rlon_radar,rheight_radar
     1             ,dgr
     1             ,grid_ref,nx_l,ny_l,nz_l
     1             ,c3_radar_subdir,path_to_vrc,r_missing_data,istatus)       

        endif

        i4_elapsed = ishow_timer()

!       go to 900

900     continue

      end if ! i_last_scan

      go to 999 ! normal return

!     return section

998   i_status = 0
      write(6,*) ' warning: return from remap_process with 0 status'
      return

999   i_status = 1
      return

      end


        subroutine purge(ext,nfiles,ntime_min,i4time_now)

!       keeps number of files according to nfiles or time span according to
!       ntime_min, whichever is greater

        integer max_files
        parameter (max_files = 1000)

        character*9 asc_tim_9
        character*31 ext
        character*255 c_filespec
        character c_fnames(max_files)*80

        call s_len(ext,len_ext)

        c_filespec = '../lapsprd/'//ext(1:len_ext)//'/*.'
     1                            //ext(1:len_ext)//'*'         

        write(6,*)c_filespec

        call    get_file_names(  c_filespec,
     1			 i_nbr_files_ret,
     1			 c_fnames,max_files,
     1			 i_status )

        if(i_nbr_files_ret .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
            write(6,*)i_nbr_files_ret,' file(s) in directory'
        else ! error condition
            write(6,*)' no files in directory'
            istatus = 0
            return
        endif

        ntime_sec = ntime_min * 60

10      do i=1,i_nbr_files_ret-nfiles ! loop through excess versions
            asc_tim_9 = c_fnames(i)(lenf+1:lenf+9)
            call i4time_fname_lp(asc_tim_9,i4time_file,istatus)
            if(i4time_now - i4time_file .gt. ntime_sec)then ! file is too old

!               delete the file
!               call rm_file(c_fnames(i)(1:lenf+13),istatus)
                call rm_file(c_fnames(i),istatus)


            endif
        enddo

        return
        end


        subroutine rm_file(c_filename,istatus)

        character*(*) c_filename

        integer istatus

        lun = 151

        write(6,*)' rm_file ',c_filename

        open(lun,file=c_filename,status='unknown')

        close(lun,status='delete')

        istatus = 1
 
        return
        end



        subroutine put_remap_vrc(i4time,comment_2d 
     1                         ,rlat_radar,rlon_radar,rheight_radar
     1                         ,dgr
     1                         ,field_3d,imax,jmax,kmax,c3_radar_subdir        
     1                         ,path_to_vrc,r_missing_data,istatus)

!       stuff from 'put_laps_2d' except how do we handle radar subdir?

        character*150 directory
        character*150 directory1
        character*31 ext

        character*125 comment_2d
        character*125 comments_2d(2)
        character*10 units(2)
        character*3 vars(2)

        integer lvl_2d(2)

        real field_3d(imax,jmax,kmax)
        real fields_2d(imax,jmax,2)
        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)
        real dist(imax,jmax)

        character*9 a9time
        character*8 radar_subdir
        character*3 c3_radar_subdir
        character*(*) path_to_vrc
        character*4 lvl_coord_2d(2)

        call make_fnam_lp(i4time,a9time,istatus)
        if(istatus .ne. 1)return

        write(6,*)' subroutine put_remap_vrc for ',a9time

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)return

!       get column max reflectivity (now passing in r_missing_data)
        call get_max_reflect(field_3d,imax,jmax,kmax,r_missing_data
     1                      ,fields_2d(1,1,1) )

        call get_laps_domain(imax,jmax,'nest7grid',lat,lon,topo,istatus)       
        if(istatus .ne. 1)then
            write(6,*)' error calling get_laps_domain'
            return
        endif

!       calculate closest radar array
        write(6,*)' calculating closest radar array (dist to vrc radar)'       
        do i = 1,imax
        do j = 1,jmax
            call latlon_to_radar(lat(i,j),lon(i,j),topo(i,j)
     1                          ,azimuth,dist(i,j),elev
     1                          ,rlat_radar,rlon_radar,rheight_radar)       
        enddo ! j
        enddo ! i

!       call vrc_clutter_thresh(      fields_2d(1,1,1)                   ! i/o
!    1                               ,dist                               ! i
!    1                               ,imax,jmax,ref_base,r_missing_data) ! i

        call ref_fill_horz(fields_2d(1,1,1),imax,jmax,1,lat,lon,dgr
     1                    ,imax,jmax,0,0
     1                    ,rlat_radar,rlon_radar,rheight_radar,istatus)       
        if(istatus .ne. 1)then
            write(6,*)' error calling ref_fill_horz'          
            return
        endif

!       utilize closest radar array
        write(6,*)' utilizing closest radar array (dist to vrc radar)'       
        do i = 1,imax
        do j = 1,jmax
            if(fields_2d(i,j,1) .ne. r_missing_data)then
                fields_2d(i,j,2) = dist(i,j)
            else
                fields_2d(i,j,2) = r_missing_data
            endif
        enddo ! j
        enddo ! i

        ext = 'vrc'

        vars(1) = 'ref'
        vars(2) = 'dis'

        units(1) = 'dbz'
        units(2) = 'm'

        lvl_2d(1) = 0
        lvl_2d(2) = 0
        
        lvl_coord_2d(1) = 'msl'
        lvl_coord_2d(2) = 'msl'

        comments_2d(1) = comment_2d
        comments_2d(2) = comment_2d

        write(6,*)'path_to_vrc = ',path_to_vrc

        if(path_to_vrc .eq. 'rdr')then
            radar_subdir = c3_radar_subdir
            write(6,*)' radar_subdir = ',radar_subdir

            call get_directory('rdr',directory1,len_dir1)

            directory = directory1(1:len_dir1)//radar_subdir(1:3)
     1                                        //'/vrc/'  
            call s_len(directory,len_dir)

        else ! 'lapsprd'
            call get_directory('vrc',directory,len_dir)

        endif            

        write(6,11)directory(1:len_dir),ext(1:5),vars
11      format(' writing 2d ',a,1x,a5,2(1x,a3))

        call write_laps_data(i4time,directory,ext,imax,jmax,
     1                       2,2,vars,lvl_2d,lvl_coord_2d,units,
     1                       comments_2d,fields_2d,istatus)
        if(istatus .eq. 1)then
            write(6,*)' vrc successfully written'
        else
            write(6,*)' vrc not successfully written', istatus
        endif


        return
        end

        subroutine vrc_clutter_thresh(ref                             ! i/o
     1                               ,dist                            ! i
     1                               ,ni,nj,ref_base,r_missing_data)  ! i

!       apply a range dependent reflectivity threshold to filter ground clutter

        real ref(ni,nj), dist(ni,nj)

        do i = 1,ni
        do j = 1,nj
            if(dist(i,j) .lt. 80000.)then
                thresh = 10.
            elseif(dist(i,j) .gt. 100000.)then
                thresh = 0.
            else
                thresh = 10. * (100000. - dist(i,j)) / 20000.
            endif             

            if(      ref(i,j) .lt. thresh 
     1         .and. ref(i,j) .gt. ref_base
     1         .and. ref(i,j) .ne. r_missing_data )then       
                ref(i,j) = ref_base
            endif
        
        enddo ! j
        enddo ! i

        return
        end
