      subroutine config_satellite_lvd(istatus)
c
cdoc  Reads static/satellite_lvd.nl file.

      character nest7grid*150
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satellite_namelist_lvd.cmn'
      include 'grid_fname.cmn'

      include 'satdata_lvd.for'

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/satellite_lvd.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_lvd_nl,err=901)
      close(1)
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading satellite_nl in ',nest7grid
      write(*,satellite_lvd_nl)
      
      print*,'**************************************************'
      print*,'Now using repository version: satellite_lvd_rep.nl'
      print*,'**************************************************'

      nest7grid = nest7grid(1:len_dir)//'/satellite_lvd_rep.nl'
      open(1,file=nest7grid,status='old',err=902)
      read(1,satellite_lvd_nl,err=903)
      close(1)
      istatus = 1
      return

 902  print*,'error opening file ',nest7grid
      return
 903  print*,'error reading satellite_lvd_rep.nl in ',nest7grid
      write(*,satellite_lvd_nl)
      return

      end
c
c =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine get_sat_sounder_info(n_sndr,
     +c_sndr_id,n_sndr_channels,path_to_sat_sounder,
     +channel_wavelength_u,imsng_sndr_pix,pct_req_lsr,
     +istatus)
c
cdoc  Reads static/sat_sounder.nl file.

      include 'lsr_dims.inc'
      include 'grid_fname.cmn'       !grid_fnam_common

      integer len_dir
      integer n_sndr
      integer n_sndr_channels
      integer imsng_sndr_pix
      integer istatus
      character*6   c_sndr_id(max_sat)
      character*150 nest7grid
      character*200 path_to_sat_sounder(max_sat)

      real*4 channel_wavelength_u(max_ch,max_sat)
      real*4 pct_req_lsr

      namelist /satellite_sounder_nl/ n_sndr,c_sndr_id,path_to_sat_sound
     +er,n_sndr_channels,channel_wavelength_u,imsng_sndr
     +_pix,pct_req_lsr

      call get_directory(grid_fnam_common,nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/sat_sounder.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_sounder_nl,err=901)
      close(1)

      istatus = 1
      return
 900  print*,'error opening file ',nest7grid
      stop
 901  print*,'error reading satellite_sounder_nl in ',nest7grid
      write(*,satellite_sounder_nl)
      stop
      end
c
c-----------------------------------------------------------
c
      subroutine get_balance_nl(lrunbal,adv_anal_by_t_min
     .           ,cpads_type,incl_clom,istatus)
c
cdoc  Reads static/balance.nl file.

      implicit none

      integer    istatus
      integer    len_dir
      integer    adv_anal_by_t_min
      logical    lrunbal
      logical    incl_clom 
      character  nest7grid*150
      character  cpads_type*3

      include   'grid_fname.cmn'       !grid_fnam_common

      namelist /balance_nl/lrunbal,adv_anal_by_t_min,cpads_type
     1        ,incl_clom

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'balance.nl'

         incl_clom=.true.

      open(1,file=nest7grid,status='old',err=900)
      read(1,balance_nl,err=901)
      close(1)
      return
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading balance.nl in ',nest7grid
      write(*,balance_nl)
      stop
      end
c
c-----------------------------------------------------------
c
      subroutine get_sfcqc_nl(lrunqc,istatus)
c
cdoc  Reads static/balance.nl file.

      implicit none

      integer    istatus
      integer    len_dir
      logical    lrunqc
      character  nest7grid*150

      include   'grid_fname.cmn'       !grid_fnam_common

      namelist /sfc_qc_nl/lrunqc

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'sfc_qc.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,sfc_qc_nl,err=901)
      close(1)
      istatus = 1
      return
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading sfc_qc.nl in ',nest7grid
      write(*,sfc_qc_nl)
      stop
      end
c
c ---------------------------------------------------------
c
      subroutine mosaic_radar_nl(c_radar_mosaic_type,n_radars,
     & c_radar_ext,i_window,mosaic_cycle_time,imosaic_3d,
     & n_radars_wideband,n_radars_narrowband,istatus)
c
cdoc  Reads static/radar_mosaic.nl file

      include    'radar_mosaic_dim.inc'
      include    'grid_fname.cmn'              !grid_fnam_common

      integer    istatus
      integer    len_dir
      integer    n_radars
      integer    i_window
      integer    imosaic_3d
      character  c_radar_mosaic_type*3
      character  c_radar_ext(max_radars_mosaic)*3
      character  nest7grid*150

      namelist /radar_mosaic_nl/c_radar_mosaic_type,n_radars,
     & c_radar_ext,i_window,mosaic_cycle_time,imosaic_3d,
     & n_radars_wideband,n_radars_narrowband

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'radar_mosaic.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,radar_mosaic_nl,err=901)
      close(1)
      return
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading radar_mosaic.nl in ',nest7grid
      write(*,radar_mosaic_nl)
      stop
      end
c
c -------------------------------------------------------------
c
      subroutine get_background_info(bgpaths,bgmodels
     +,forecast_length
     +,use_analysis,cmodel,itime_inc,smooth_fields,luse_sfc_bkgd
     +,lgb_only)

cdoc reads static/background.nl

      implicit none
      include 'bgdata.inc'
      include 'grid_fname.cmn'             !grid_fnam_common

      character*150 nest7grid
      character*256 bgpaths(maxbgmodels)
      character*132 cmodel(maxbgmodels)
      integer bgmodels(maxbgmodels), len_dir
      integer forecast_length
      integer itime_inc
      logical luse_sfc_bkgd
      logical use_analysis
      logical smooth_fields
      logical lgb_only
      namelist /background_nl/bgpaths,bgmodels
     +,forecast_length
     +,use_analysis,cmodel,itime_inc,smooth_fields,luse_sfc_bkgd
     +,lgb_only

      smooth_fields = .false.
      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif
      nest7grid = nest7grid(1:len_dir)//'background.nl'

      open(1,file=nest7grid(1:len_dir+13),status='old',err=900)
      read(1,background_nl,err=901)
      close(1)
      return
 900  print*,'error opening file ',nest7grid
      stop
 901  print*,'error reading background_nl in ',nest7grid
      write(*,background_nl)
      stop
      end
c
c --- OSSE namelist reader ---
c
      subroutine get_osse_information(path_to_model,
     1                                cmodel,
     1                                c_obs_types,
     1                                n_sim_obs,
     1                                c_simob_fnames,
     1                                a9_time_init,
     1                                a4_time_fcst,
     1                                ifcst_intrvl,
     1                                isim_time_hr,
     1                                istatus)
c
cdoc  Reads static/osse.nl file.

      implicit none

      include 'grid_fname.cmn'                 !grid_fnam_common
      include 'osse.inc'

      integer        istatus
      integer        len_dir
      integer        i,n_sim_obs
      integer        ifcst_intrvl
      integer        isim_time_hr
      integer        n_simfiles

      integer        nsimobs
      data           nsimobs/0/
      save           nsimobs

      character      nest7grid*150
      character      path_to_model*150
      character      c_obs_types(maxobtype)*15
      character      c_simob_fnames(nsimfiles)*50
      character      a9_time_init*9
      character      a4_time_fcst*4
      character      cmodel*10

      namelist /osse_nl/path_to_model,cmodel,c_obs_types
     1,c_simob_fnames,a9_time_init,a4_time_fcst,ifcst_intrvl
     1,isim_time_hr

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif


      do i = 1,maxobtype
         c_obs_types(i) = ' '
      enddo
      do i = 1,nsimfiles
         c_simob_fnames(i) = ' '
      enddo


      nest7grid = nest7grid(1:len_dir)//'osse.nl'
      print*,'nest7grid = ',nest7grid(1:len_dir+7)

      open(1,file=nest7grid,status='old',err=900)
      read(1,osse_nl,err=901)
      close(1)

      if(nsimobs.eq.0)then
         do i = 1,maxobtype
            if(c_obs_types(i).ne.' ')then
               n_sim_obs = n_sim_obs + 1
            endif
         enddo
         nsimobs=n_sim_obs
      endif
      n_sim_obs=nsimobs
      print*,'n_sim_obs types from namelist= ',n_sim_obs
      istatus = 1
     
      return
       
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading osse.nl in ',nest7grid
      write(*,osse_nl)
      stop
      end
c
c ---------------------------------------------------------------
c
       subroutine get_wind_parms(l_use_raob,l_use_cdw,l_use_radial_vel
     1                          ,thresh_2_radarobs_lvl_unfltrd
     1                          ,thresh_4_radarobs_lvl_unfltrd
     1                          ,weight_bkg_const_wind
     1                          ,weight_radar
     1                          ,rms_thresh_wind
     1                          ,max_pr,max_pr_levels,i_3d
     1                          ,istatus)

       include 'grid_fname.cmn'                          !grid_fnam_common

       logical l_use_raob, l_use_cdw, l_use_radial_vel
       integer*4 thresh_2_radarobs_lvl_unfltrd
     1          ,thresh_4_radarobs_lvl_unfltrd

       namelist /wind_nl/ l_use_raob, l_use_cdw, l_use_radial_vel
     1                   ,thresh_2_radarobs_lvl_unfltrd
     1                   ,thresh_4_radarobs_lvl_unfltrd
     1                   ,weight_bkg_const_wind
     1                   ,weight_radar
     1                   ,rms_thresh_wind
     1                   ,max_pr,max_pr_levels,i_3d
 
       character*150 static_dir,filename
 
!      thresh_2_radarobs_lvl_unfltrd = 300 ! temporary assignment
!      thresh_4_radarobs_lvl_unfltrd = 600 ! temporary assignment

       call get_directory(grid_fnam_common,static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/wind.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,wind_nl,err=901)
       close(1)

       if(max_pr .le. 0)then
           write(6,*)' ERROR: invalid or uninitialized value of '
     1              ,'max_pr in wind.nl ',max_pr
           istatus = 0
           return
       endif

       if(max_pr_levels .le. 0)then
           write(6,*)' ERROR: invalid or uninitialized value of '
     1              ,'max_pr_levels in wind.nl ',max_pr_levels
           istatus = 0
           return
       endif

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading wind_nl in ',filename
       write(*,wind_nl)
       istatus = 0
       return

       end
c
c----------------------------------------
c
       subroutine get_gridnl(mode)
       implicit none
       integer mode
       integer len
       integer istatus
       character*256 directory
       character*256 fname
       character*200 cdataroot
       character*10  c10_grid_fname
       namelist /grid_nl/ mode

       mode = 0
       call find_domain_name(cdataroot,c10_grid_fname,istatus)
       if(istatus.ne.1)then
          print*,'Error returned from find_domain_name'
          return
       endif
       call get_directory(c10_grid_fname,directory,len)
       fname = directory(1:len)//'grid.nl'
       open(3,file=fname,status='old',err=101)
       read(3,grid_nl,err=101,end=101)
       close(3)

101    return
       end
c
c --------------------------------------------------------------
c
       subroutine read_verif_nl(type_obs,path_to_raw_profiler,
     1path_to_raw_sounding, raob_process_lag,  raob_process_lag_bal,
     1                          max_verif, verif_output_dir,
     1                          verif_missing_data, n_verif, istatus)

       implicit none

       include 'grid_fname.cmn'                          !grid_fnam_common

       integer          max_verif
       character*1	type_obs
       character*150	path_to_raw_profiler
       character*150	path_to_raw_sounding
       integer*4	raob_process_lag_Bal
       integer*4	raob_process_lag
       integer          n_verif
       integer          i,len
       character*150	verif_output_dir(4)
       real*4		verif_missing_data
       integer		istatus
       
       namelist /verif_nl/ type_obs,path_to_raw_profiler,
     1path_to_raw_sounding, raob_process_lag, raob_process_lag_bal,
     1                  verif_output_dir,
     1                  verif_missing_data

       character*150    static_dir,filename
       integer		len_dir

       istatus = 0
       call get_directory(grid_fnam_common,static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/verif.nl'

       open(1,file=filename,status='old',err=900)
       read(1,verif_nl,err=901)
       close(1)

       n_verif=0
       do i=1,max_verif
          call s_len(verif_output_dir(i),len)
          if(len.gt.0)then
             n_verif=n_verif+1
          endif
       enddo
       if(n_verif.le.0)then
          print*,'Error! Check namelist verif.nl variable'
          print*,'verif_output_dir. It appears to be empty'
          return
       endif
c      print*,'Number of verification types from verif.nl: ',n_verif
c      call s_len(verif_output_dir(1),len)
c      do i=1,n_verif
c         print*,i,' ',verif_output_dir(i)(1:len)
c      enddo

       istatus = 1
       return

  900  print*,'error opening file ',filename
       return

  901  print*,'error reading verif_nl in ',filename
       write(*,verif_nl)
       return

       end
c
c --------------------------------------------------------------
c
       subroutine read_sfc_nl(use_lso_qc,skip_internal_qc 
     1                       ,itheta, redp_lvl, del, gam, ak
     1                       ,l_require_lso
     1                       ,bad_t,bad_td,bad_u,bad_v,bad_p
     1                       ,bad_mp,bad_th,bad_the
     1                       ,bad_vis,bad_tb8
     1                       ,thresh_t,thresh_td,thresh_mslp
     1                       ,sfc_nl_parms,istatus)

       implicit none

       real    badflag
       include 'laps_sfc.inc'

       include 'grid_fname.cmn'                          !grid_fnam_common
       integer use_lso_qc, skip_internal_qc, itheta
       logical l_require_lso
       real    redp_lvl,del,gam,ak
       real    bad_t,bad_td,bad_u,bad_v,bad_p
       real    bad_mp,bad_th,bad_the
       real    bad_vis,bad_tb8
       real    thresh_t,thresh_td,thresh_mslp
       real    rms_wind, rms_temp, rms_dewpoint

       integer istatus
       
       namelist /surface_analysis/  use_lso_qc,skip_internal_qc,
     1                              itheta, redp_lvl, del, gam, ak,       
     1                              l_require_lso,
     1                              bad_t,bad_td,bad_u,bad_v,bad_p,
     1                              bad_mp,bad_th,bad_the,
     1                              bad_vis,bad_tb8,
     1                              thresh_t,thresh_td,thresh_mslp,
     1                              rms_wind, rms_temp, rms_dewpoint

       character*150    static_dir,filename
       integer len_dir

       istatus = 0
       call get_directory(grid_fnam_common,static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/surface_analysis.nl'

       open(1,file=filename,status='old',err=900)
       read(1,surface_analysis,err=901)
       close(1)

       sfc_nl_parms%rms_wind = rms_wind
       sfc_nl_parms%rms_temp = rms_temp
       sfc_nl_parms%rms_dewpoint = rms_dewpoint

       istatus = 1
       return

  900  print*,'error opening file ',filename
       return

  901  print*,'error reading verif_nl in ',filename
       write(*,surface_analysis)
       return

       end
c
c --------------------------------------------------------------
c
       subroutine get_laps_redp(redp_lvl,istatus)

       implicit none
       real    badflag
       include 'laps_sfc.inc'
       include 'grid_fname.cmn'                          !grid_fnam_common
       integer use_lso_qc, skip_internal_qc, itheta
       logical l_require_lso
       real    redp_lvl,del,gam,ak
       real    bad_t,bad_td,bad_u,bad_v,bad_p
       real    bad_mp,bad_th,bad_the
       real    bad_vis,bad_tb8
       real    thresh_t,thresh_td,thresh_mslp
       real    rms_wind, rms_temp, rms_dewpoint
       integer istatus
       namelist /surface_analysis/  use_lso_qc,skip_internal_qc,
     1                              itheta, redp_lvl, del, gam, ak,
     1                              l_require_lso,
     1                              bad_t,bad_td,bad_u,bad_v,bad_p,
     1                              bad_mp,bad_th,bad_the,
     1                              bad_vis,bad_tb8,
     1                              thresh_t,thresh_td,thresh_mslp,
     1                              rms_wind, rms_temp, rms_dewpoint
       character*150    static_dir,filename
       integer len_dir
       istatus = 0
       call get_directory(grid_fnam_common,static_dir,len_dir)
       filename = static_dir(1:len_dir)//'/surface_analysis.nl'
       open(1,file=filename,status='old',err=903)
       read(1,surface_analysis,err=904)
       close(1)
       sfc_nl_parms%rms_wind = rms_wind
       sfc_nl_parms%rms_temp = rms_temp
       sfc_nl_parms%rms_dewpoint = rms_dewpoint
       istatus = 1
       return
  903  print*,'error opening file ',filename
       return
  904  print*,'error reading sfc_nl in ',filename
       write(*,surface_analysis)
       return

       end
