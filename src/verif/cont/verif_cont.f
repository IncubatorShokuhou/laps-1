
        program verif_radar_main

        use mem_namelist, only: model_fcst_intvl,model_cycle_time
     1                                          ,model_fcst_len

        character*9 a9time

        character*150 verif_dir, n_plot_times_file              

        logical l_persist 

        l_persist = .true. ! add persistence forecast for evaluation

        call get_systime(i4time,a9time,istatus)
        if(istatus .ne. 1)go to 999

        write(6,*)' systime = ',a9time

        call get_grid_dim_xy(nx_l,ny_l,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting horizontal domain dimensions'
           go to 999
        endif

        call get_laps_dimensions(nz_l,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting vertical domain dimension'
           go to 999
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting r_missing_data'
           go to 999
        endif

        call get_laps_cycle_time(laps_cycle_time,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting laps_cycle_time'
           go to 999
        endif

        istat = init_timer()
        i4_elapsed = ishow_timer()

        call verif_radar(i4time,a9time,model_fcst_intvl,
     1                  model_fcst_len,
     1                  laps_cycle_time,
     1                  nx_l,ny_l,
     1                  nz_l,
     1                  r_missing_data,
     1                  l_persist,
     1                  j_status)

        i4_elapsed = ishow_timer()

!       read n_plot_times from file
        call get_directory('verif',verif_dir,len_verif)
        lun_plot_times = 42
        n_plot_times_file = verif_dir(1:len_verif)//'/n_fcst_times.dat'
        open(lun_plot_times,file=n_plot_times_file,status='old')
        read(lun_plot_times,*)n_plot_times
        close(lun_plot_times)

        write(6,*)                                           
        write(6,*)' calling verif_radar_composite...'
        write(6,*)' time is ',i4time,a9time

        call verif_radar_composite(i4time,a9time,model_fcst_intvl,
     1                  model_fcst_len,
     1                  model_cycle_time,
     1                  laps_cycle_time,
     1                  nx_l,ny_l,
     1                  nz_l,
     1                  r_missing_data,
     1                  n_plot_times,
     1                  l_persist,
     1                  j_status)

        i4_elapsed = ishow_timer()

999     continue

        write(6,*)' end of verif_radar_main program...'

        end
          
        subroutine verif_radar(i4time_sys,a9time,model_fcst_intvl,
     1                  model_fcst_len,
     1                  laps_cycle_time,
     1                  nx_l,ny_l,
     1                  nz_l,
     1                  r_missing_data,
     1                  l_persist,
     1                  j_status)

        use constants_laps, only: inch2m

        include 'lapsparms.for' ! maxbgmodels

        real var_anal_3d(nx_l,ny_l,nz_l)
        real var_fcst_3d(nx_l,ny_l,nz_l)
        real var_prst_3d(nx_l,ny_l,nz_l)
        real rqc(nx_l,ny_l)

        logical lmask_and_3d(nx_l,ny_l,nz_l)
        logical lmask_rqc_3d(nx_l,ny_l,nz_l)
        logical l_col /.true./

!       integer       maxbgmodels
!       parameter     (maxbgmodels=10)
        character*30  c_fdda_mdl_src(maxbgmodels)

        integer max_fcst_times
        parameter (max_fcst_times=200)

        integer max_regions
        parameter (max_regions=10)

        integer il(maxbgmodels,0:max_fcst_times,max_regions)
        integer ih(maxbgmodels,0:max_fcst_times,max_regions)
        integer jl(maxbgmodels,0:max_fcst_times,max_regions)
        integer jh(maxbgmodels,0:max_fcst_times,max_regions)

        character ext*31, directory*255, c_model*30

        character*10  units_2d
        character*125 comment_2d
        character*6 var_2d               
        character*3 var_2d_anal
        character*3 var_2d_fcst
        character*5 fcst_hh_mm
        character*9 a9time,a9time_valid,a9time_initial
        character*24 a24time_valid
        character*150 hist_dir, cont_dir, verif_dir, plot_dir
        character*150 hist_file, bias_file, ets_file, members_file

        integer n_fields
        parameter (n_fields=5)
        character*10 ext_anal_a(n_fields), ext_fcst_a(n_fields)
        character*10 var_a(n_fields)
        character*10 type_a(n_fields)
        integer nthr_a(n_fields)  ! number of thresholds for each field
        integer ndims_a(n_fields) ! number of dimensions for each field
        real dbz_an(7),dbz_fc(7),pcp_thr(7)
        character*2 c2_region
        character*10 c_thr

!       specify what is being verified
        data ext_fcst_a /'fua','fsf','fsf','fsf','fsf'/                        
        data ext_anal_a /'lps','lmr','st4','st4','st4'/                        
        data var_a      /'ref','lmr','pcp_01','pcp_03','pcp_06'/
        data type_a     /'rdr','rdr','pcp','pcp','pcp'/
        data nthr_a     /5,5,7,7,7/        
        data ndims_a    /3,2,2,2,2/        

        data pcp_thr /.01,.05,0.1,0.5,1.0,2.0,5.0/

        integer,parameter :: k12 = selected_int_kind(12)
        integer (kind=k12) :: contable(0:1,0:1)

        integer maxthr
        parameter (maxthr=7)

!       real cont_4d(nx_l,ny_l,nz_l,maxthr)
        real bias(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        real ets(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        real 
     1  frac_coverage(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        real frac_obs(maxbgmodels,0:max_fcst_times,max_regions,maxthr) 
        real frac_fcst(maxbgmodels,0:max_fcst_times,max_regions,maxthr) 
        integer 
     1  n(maxbgmodels,0:max_fcst_times,max_regions,maxthr,0:1,0:1)


        logical l_persist, l_good_persist

        rmiss = -999.
        imiss = -999

        thresh_var = 20. ! lowest threshold for this variable

        i4_initial = i4time_sys

!       determine verification timing
        model_verif_intvl = max(laps_cycle_time,model_fcst_intvl)
        n_fcst_times = (model_fcst_len*60) / model_verif_intvl

        write(6,*)' model_verif_intvl = ',model_verif_intvl
        write(6,*)' n_fcst_times = ',n_fcst_times

        lun_in = 21

!       get fdda_model_source and 'n_fdda_models' from static file
        call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

        if(l_persist .eqv. .true.)then
            n_fdda_models = n_fdda_models + 1
            c_fdda_mdl_src(n_fdda_models) = 'persistence'
            write(6,*)' adding persistence to fdda_models'
        endif

        write(6,*)' n_fdda_models = ',n_fdda_models
        write(6,*)' c_fdda_mdl_src = '
     1            ,(c_fdda_mdl_src(m),m=1,n_fdda_models)

!       read in data file with region points
        n_models = n_fdda_models
        call read_region_info(maxbgmodels,max_fcst_times,max_regions
     1                       ,n_models,n_fcst_times,n_regions
     1                       ,il,ih,jl,jh,lun_in)

!       in 'nest7grid.parms' model #1 is lga usually
!       in 'verif_regions.dat' model #1 represents the analysis

!       if(n_fdda_models .ne. n_models)then
!           write(6,*)' error n_models differs from n_fdda_models '
!    1                       ,n_models,n_fdda_models
!           stop
!       endif

        do ifield = 1,n_fields

!        initialize arrays
         bias = -999.
         ets = -999.
         frac_coverage = -999.
         n = -999

         var_prst_3d = r_missing_data
         l_good_persist = .false.

         var_2d = trim(var_a(ifield))
         call s_len(var_2d,lenvar)

         if(trim(type_a(ifield)) .eq. 'rdr')then
             nk_cont = nz_l
             istart = 0
         else
             nk_cont = 1
             istart = 1
         endif

         do imodel = 1,n_fdda_models

          c_model = c_fdda_mdl_src(imodel)

          if(c_model(1:3) .ne. 'lga')then

           write(6,*)
           write(6,*)' processing model ',c_model

           call s_len(c_model,len_model)

           call get_directory('verif',verif_dir,len_verif)

           hist_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                      //'/hist/'
     1                                      //c_model(1:len_model)
           len_hist = len_verif + 6 + lenvar + len_model

           cont_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                      //'/cont/'
     1                                      //c_model(1:len_model)
     1                                      //'/'
           len_cont = len_verif + 6 + lenvar + len_model

           do itime_fcst = istart,n_fcst_times

!            itime = itime_fcst + 1

             do iregion = 1,n_regions ! 1 for testing

                ilow  = il(imodel,itime_fcst,iregion)
                ihigh = ih(imodel,itime_fcst,iregion)
                jlow  = jl(imodel,itime_fcst,iregion)
                jhigh = jh(imodel,itime_fcst,iregion)

                i4_valid = i4_initial + itime_fcst * model_verif_intvl 

!               write(*,*)'beka1',i4_valid,l_good_persist,len_model

                call make_fnam_lp(i4_valid,a9time_valid,istatus)
                call make_fnam_lp(i4_initial,a9time_initial,istatus)

                write(6,*)
                write(6,*)' histograms for forecast time step '
     1                   ,itime_fcst,' region = ',iregion

                lun_out = 11
                lun_bias = 12
                lun_ets  = 13
                lun_mem  = 14

                write(c2_region,1)iregion
 1              format(i2.2)

!               perhaps a14time should be used (call make_fcst_time)?
                call make_fcst_time(i4_valid,i4_initial
     1                             ,fcst_hh_mm,istatus)

                hist_file = hist_dir(1:len_hist)//'/'//a9time_initial    
     1                                          //trim(fcst_hh_mm)       
     1                                          //'_'//c2_region     
     1                                          //'.hist'     

                write(6,*)'hist_file = ',hist_file

                open(lun_out,file=hist_file,status='unknown',err=401)
                goto 402
401             write(6,*)' error opening hist_file: ',trim(hist_file)
                goto 900
402             continue

                if(iregion .eq. 1)then

                  if(ndims_a(ifield) .eq. 3)then ! read analyzed 3d field
!                   write(*,*)'beka2a',i4_valid,l_good_persist,len_model
                    ext = ext_anal_a(ifield)
                    call get_laps_3d(i4_valid,nx_l,ny_l,nz_l
     1              ,ext,var_2d,units_2d,comment_2d,var_anal_3d,istatus)       
                    if(istatus .ne. 1)then
                        write(6,*)' error reading 3d analysis for '
     1                            ,var_2d
                        goto 900
                    endif

                  else ! read analyzed 2d field
!                   write(*,*)'beka2b',i4_valid,l_good_persist,len_model
                    ext = ext_anal_a(ifield)

                    intvl_pcp = -999

                    if(trim(var_2d) .eq. 'lmr')then
                      var_2d_anal = 'r'
                    elseif(trim(var_2d) .eq. 'pcp_01')then
                      var_2d_anal = 'ppt'
                      intvl_pcp = 3600 
                    elseif(trim(var_2d) .eq. 'pcp_03')then
                      var_2d_anal = 'nul'
                      intvl_pcp = 10800
                    elseif(trim(var_2d) .eq. 'pcp_06')then
                      var_2d_anal = 'nul'
                      intvl_pcp = 21600
                    elseif(trim(var_2d) .eq. 'pcp_24')then
                      var_2d_anal = 'nul'
                      intvl_pcp = 86400
                    else
                      var_2d_anal = trim(var_2d(1:3))
                    endif

                    i4_fcst = i4_valid - i4_initial

                    if(var_2d_anal .ne. 'nul')then
                        call get_laps_2d(i4_valid,ext,var_2d_anal
     1                  ,units_2d,comment_2d,nx_l,ny_l,var_anal_3d
     1                  ,istatus)       
                    elseif(i4_fcst .eq. (i4_fcst/intvl_pcp) 
     1                                         * intvl_pcp  )then  
                        write(6,*)' call get_interval_precip for st4...'
                        call get_interval_precip(' ',ext
     1                        ,i4_valid - intvl_pcp,i4_valid
     1                        ,laps_cycle_time,i4_dum
     1                        ,r_missing_data
     1                        ,nx_l,ny_l,1,var_anal_3d,istatus)
                    else
                        write(6,*)' no analysis for this time/field'    
                        goto 900
                    endif

                    if(istatus .ne. 1 .and. istatus .ne. -1)then
                        write(6,*)' error reading 2d analysis for '
     1                           ,trim(ext),' ',trim(var_2d_anal)
                        goto 900
                    endif

                    if(trim(type_a(ifield)) .eq. 'rdr')then 
                        do k = 2,nz_l
                            var_anal_3d(:,:,k) = var_anal_3d(:,:,1)
                        enddo ! k
                    endif

                  endif

                  if(itime_fcst .eq. 0)then
                      var_prst_3d = var_anal_3d
                      l_good_persist = .true.
                      write(6,*)
     1         ' setting persistence to analysis (and l_good_persist=t)'       
                  endif

!                 write(*,*)'beka3',i4_valid,l_good_persist,len_model

                  if(trim(type_a(ifield)) .eq. 'rdr')then
                      if(var_2d .eq. 'llr')then 
                          rqc_thresh = 2.0
                      else
                          rqc_thresh = 3.0
                      endif

                      ext = 'lcv'
                      call get_laps_2d(i4_valid,ext,'rqc',units_2d
     1                                ,comment_2d,nx_l,ny_l,rqc,istatus)
                      if(istatus .ne. 1)then
                          write(6,*)' error reading 2d rqc analysis'
                          goto 900
                      endif
                      write(6,*)' successful radar quality (rqc) read'
                  else ! determine precip quality
                      rqc_thresh = 0.
                      rqc(:,:) = var_anal_3d(:,:,1)
                      write(6,*)' range of rqc is ',minval(rqc)
     1                                             ,maxval(rqc)
                  endif

                  if(c_fdda_mdl_src(imodel) .ne. 'persistence')then

                     if(ndims_a(ifield) .eq. 3)then ! read forecast 3d field
                        ext = trim(ext_fcst_a(ifield))
                        call get_directory(ext,directory,len_dir)
                        directory=directory(1:len_dir)
     1                                      //c_model(1:len_model)    
     1                                      //'/'

                        call get_lapsdata_3d(i4_initial,i4_valid
     1                          ,nx_l,ny_l,nz_l       
     1                          ,directory,var_2d
     1                          ,units_2d,comment_2d,var_fcst_3d
     1                          ,istatus)
                        if(istatus .ne. 1)then
                             write(6,*)' error reading 3d forecast for '
     1                                 ,var_2d
                             goto 900
                        endif

                     else ! read forecast 2d field
                        ext = trim(ext_fcst_a(ifield))

                        call get_directory(ext,directory,len_dir)
!                       write(6,*)' len_model/len_dir/directory = '
!    1                      ,len_model,len_dir,trim(directory)
                        directory=directory(1:len_dir)
     1                                      //c_model(1:len_model)    
     1                                      //'/'

                        write(6,*)' directory = ',trim(directory)

                        if(trim(c_model) .eq. 'nam')then
                             model_pcp_intvl = 10800
                        elseif(trim(c_model) .eq. 'nam-nh')then
                             model_pcp_intvl = 21600
                        else
                             model_pcp_intvl = model_fcst_intvl
                        endif

                        if(trim(type_a(ifield)) .eq. 'rdr')then
                            write(6,*)' radar verif type block'
                            var_2d_fcst = trim(var_2d(1:3))
                        else ! find cases where 'r01' field is appropriate
                            write(6,*)' precip verif type block'
                            var_2d_fcst = 'nul'
                            if(trim(var_2d) .eq. 'pcp_01')then
                                if(model_fcst_intvl .eq. 3600 .and.
     1                             model_fcst_intvl .ge. model_pcp_intvl
     1                                                             )then
                                    var_2d_fcst = 'r01' 
                                endif
                            elseif(trim(var_2d) .eq. 'pcp_03')then
                                if(model_fcst_intvl .eq. 10800 .and.
     1                             model_fcst_intvl .ge. model_pcp_intvl
     1                                                             )then
                                    var_2d_fcst = 'r01' 
                                endif
                            elseif(trim(var_2d) .eq. 'pcp_06')then
                                if(model_fcst_intvl .eq. 21600 .and.
     1                             model_fcst_intvl .ge. model_pcp_intvl
     1                                                             )then
                                    var_2d_fcst = 'r01' 
                                endif
                            elseif(trim(var_2d) .eq. 'pcp_24')then
                                if(model_fcst_intvl .eq. 86400 .and.
     1                             model_fcst_intvl .ge. model_pcp_intvl
     1                                                             )then
                                    var_2d_fcst = 'r01' 
                                endif
                            endif
                        endif

                        write(6,*)' var_2d/var_2d_fcst = ',var_2d,' '
     1                                                    ,var_2d_fcst

                        if(trim(type_a(ifield)) .ne. 'pcp' .or.
     1                          var_2d_fcst .eq. 'r01'         )then 
                            call get_lapsdata_2d(i4_initial,i4_valid
     1                          ,directory,var_2d_fcst
     1                          ,units_2d,comment_2d
     1                          ,nx_l,ny_l
     1                          ,var_fcst_3d
     1                          ,istatus)
                            if(istatus .ne. 1)then
                                write(6,*)
     1                           ' error reading 2d forecast for '      
     1                           ,trim(ext),' ',trim(var_2d_fcst)
                                goto 900
                            endif
                        elseif(i4_fcst .eq. (i4_fcst/intvl_pcp) 
     1                                             * intvl_pcp  
     1                                  .and.
     1                         model_pcp_intvl .le. intvl_pcp)then ! fsf

                            write(6,*)
     1                      ' call get_interval_precip for ',trim(ext)
                            call get_interval_precip('r',c_model  
     1                        ,i4_valid - intvl_pcp,i4_valid
     1                        ,model_pcp_intvl,i4_initial
     1                        ,r_missing_data
     1                        ,nx_l,ny_l,1,var_fcst_3d,istatus)
                            if(istatus .ne. 1)then
                                write(6,*)
     1                          ' bad status from get_interval_precip'
                                goto 900
                            endif
                        else
                            write(6,*)' no precip fcst at this time...'
                            write(6,*)
     1                       ' i4_fcst/intvl_pcp/model_pcp_intvl ',
     1                         i4_fcst,intvl_pcp,model_pcp_intvl
                            goto 900
                        endif

                        if(trim(type_a(ifield)) .eq. 'rdr')then 
                            do k = 2,nz_l
                                var_fcst_3d(:,:,k) = var_fcst_3d(:,:,1)
                            enddo ! k
                         endif

                     endif

                  elseif(l_good_persist .eqv. .true.)then
                      write(6,*)' setting forecast to persistence'
                      var_fcst_3d = var_prst_3d
                  else
                      write(6,*)' persistence fcst unavailable'
                      goto 900
                  endif

!                 calculate "and" as well as "rqc" "and" mask
                  n_rqc_and = 0
                  n_rqc_and_pot = 0
                  do k = 1,nk_cont
                  do i = 1,nx_l
                  do j = 1,ny_l
                     lmask_and_3d(i,j,k) = .false.
                     if(var_anal_3d(i,j,k) .ne. r_missing_data .and. 
     1                  var_anal_3d(i,j,k) .ge. thresh_var .and.
     1                  var_fcst_3d(i,j,k) .ne. r_missing_data .and.
     1                  var_fcst_3d(i,j,k) .ge. thresh_var        )then
                         lmask_and_3d(i,j,k) = .true.
                     endif
                     if(lmask_and_3d(i,j,k))then
                         n_rqc_and_pot = n_rqc_and_pot + 1
                         if(rqc(i,j) .lt. rqc_thresh .or.
     1                      rqc(i,j) .eq. r_missing_data)then
                             n_rqc_and = n_rqc_and + 1
                             lmask_and_3d(i,j,k) = .false.
                         endif
                     endif
                  enddo ! j
                  enddo ! i
                  enddo ! k

!                 calculate "rqc" "all" mask
                  lmask_rqc_3d = .true.
                  n_rqc_all = 0
                  n_rqc_all_pot = 0
                  do k = 1,nk_cont
                  do i = 1,nx_l
                  do j = 1,ny_l
                    n_rqc_all_pot = n_rqc_all_pot + 1
                    if(rqc(i,j) .lt. rqc_thresh .or.
     1                 rqc(i,j) .eq. r_missing_data)then
                        n_rqc_all = n_rqc_all + 1
                        lmask_rqc_3d(i,j,k) = .false.
                    endif
                  enddo ! j
                  enddo ! i
                  enddo ! k

                  write(6,*)' time, # of qc (2d) points (and/all)'
     1                     ,itime_fcst,n_rqc_and,n_rqc_all
                  if(n_rqc_and_pot .eq. 0)then
                      write(6,*)
     1                ' n_rqc_and_pot = 0, no echoes in 20dbz and mask '       
!                     goto 900
                  else
                      write(6,*)' time, % of qc (2d) points (and/all)'
     1                     ,itime_fcst
     1                     ,float(n_rqc_and)/float(n_rqc_and_pot) * 100.
     1                     ,float(n_rqc_all)/float(n_rqc_all_pot) * 100.
                  endif

                  if(n_rqc_all .eq. n_rqc_all_pot)then
                      write(6,*)' warning, 100% failed 3d radar qc test'
                      write(6,*)' radar data is apparently all 2d'
                      write(6,*)' histograms will be zeroed out'
                  endif

                endif ! iregion .eq. 1

                write(lun_out,*)
                write(lun_out,*)' region/irange/jrange ',iregion
     1                         ,ilow,ihigh,jlow,jhigh

                write(lun_out,*)
                write(lun_out,*)' no mask is in place'

                write(lun_out,*)
                write(lun_out,*)
     1                ' calling radarhist for analysis at ',a9time_valid 
                call radarhist(nx_l,ny_l,nz_l,var_anal_3d
     1                        ,ilow,ihigh,jlow,jhigh
     1                        ,lmask_rqc_3d,lun_out)       

                write(lun_out,*)
                write(lun_out,*)' calling radarhist for',itime_fcst
     1                    ,' hr forecast valid at ',a9time_valid
                call radarhist(nx_l,ny_l,nz_l,var_fcst_3d
     1                        ,ilow,ihigh,jlow,jhigh
     1                        ,lmask_rqc_3d,lun_out)

                write(lun_out,*)
                write(lun_out,*)
     1    ' 3-d and mask is in place with dbz threshold of ',thresh_var       

                write(lun_out,*)
                write(lun_out,*)
     1                ' calling radarhist for analysis at ',a9time_valid       
                call radarhist(nx_l,ny_l,nz_l,var_anal_3d
     1                        ,ilow,ihigh,jlow,jhigh
     1                        ,lmask_and_3d,lun_out)

                write(lun_out,*)
                write(lun_out,*)' calling radarhist for',itime_fcst
     1               ,' hr forecast valid at ',a9time_valid
                call radarhist(nx_l,ny_l,nz_l,var_fcst_3d
     1                        ,ilow,ihigh,jlow,jhigh
     1                        ,lmask_and_3d,lun_out)

                nthr = nthr_a(ifield)

                do iter = 1,1

!                calculate contingency tables
!                radar qc (rqc) for "all" case is used via 'lmask_rqc_3d'
                 do idbz = 1,nthr
                  if(trim(type_a(ifield)) .eq. 'pcp')then
                      rdbz = pcp_thr(idbz) * inch2m
!                     write(6,*)' pcp_thr / inch2m / rdbz = '
!    1                         ,pcp_thr(idbz),inch2m,rdbz
                  else
                      rdbz = float(idbz*10) + 10
                  endif

                  if(iter .eq. 1)then
                   dbz_an(idbz) = rdbz
                   dbz_fc(idbz) = rdbz
                  endif

                  write(lun_out,*)
                  write(lun_out,*)' calculate contingency table for '
     1                           ,rdbz,' dbz'
                  write(6,*)' calculate contingency table for '
     1                           ,rdbz,' dbz'
                  write(lun_out,*)' region = ',iregion
     1                           ,ilow,ihigh,jlow,jhigh
                  call contingency_table(var_anal_3d,var_fcst_3d     ! i
     1                                  ,nx_l,ny_l,nk_cont           ! i
     1                                  ,dbz_an(idbz),dbz_fc(idbz)   ! i
     1                                  ,lun_out                     ! i
     1                                  ,ilow,ihigh,jlow,jhigh       ! i
     1                                  ,lmask_rqc_3d                ! i
     1                                  ,contable)                   ! o
                  if(trim(type_a(ifield)) .eq. 'pcp')then
                      write(6,801)rdbz / inch2m, contable
801                   format(' contingency table for ',f5.2
     1                                               ,' in is:',4i9)       
!                     write(6,*)' dbz_an = ',dbz_an(idbz)
                  else
                      write(6,802)rdbz,contable
802                   format(' contingency table for ',f5.0
     1                                               ,' dbz is:',4i9)       
                  endif


!                 calculate/write skill scores
                  call skill_scores(contable,lun_out                   ! i
     1                  ,frac_coverage(imodel,itime_fcst,iregion,idbz) ! o
     1                  ,frac_obs(imodel,itime_fcst,iregion,idbz)      ! o 
     1                  ,frac_fcst(imodel,itime_fcst,iregion,idbz)     ! o 
     1                  ,bias(imodel,itime_fcst,iregion,idbz)          ! o
     1                  , ets(imodel,itime_fcst,iregion,idbz))         ! o

                  n(imodel,itime_fcst,iregion,idbz,:,:) = contable 

                  if(iregion .eq. 1)then

!                     calculate contingency table (3-d)
!                     call calc_contable_3d(
!    1                        var_anal_3d,var_fcst_3d
!    1                       ,rdbz,nx_l,ny_l,nz_l                    ! i
!    1                       ,lmask_rqc_3d,r_missing_data            ! i
!    1                       ,cont_4d(1,1,1,idbz))                   ! o

                  endif ! iregion = 1

                 enddo ! idbz

!                correlate dbz values with bias
                 write(6,*)' dbz_an (before call to calc_thresholds) = '
     1                                                          ,dbz_an 
                 call calc_thresholds(itime_fcst,nthr,dbz_an,dbz_fc    ! i
     1                       ,frac_obs(imodel,itime_fcst,iregion,:)    ! i    
     1                       ,frac_fcst(imodel,itime_fcst,iregion,:)   ! i
     1                       ,bias(imodel,itime_fcst,iregion,:)        ! i
     1                       ,rmiss                                )   ! i

                 write(6,*)
                 write(6,*)' end of iteration ',iter

                enddo ! iter

             enddo ! iregion

             close (lun_out) 

 900       enddo ! itime_fcst

          endif ! c_model .ne. lga

         enddo ! model

         iregion = 1

         write(6,*)' nthr before assigning = ',nthr
         nthr = nthr_a(ifield) ! nthr may be unset or have earlier been stepped on
         do idbz = 1,nthr

           if(trim(type_a(ifield)) .eq. 'rdr')then
               rdbz = float(idbz*10) + 10
               write(c_thr,901)nint(rdbz)
 901           format(i2)
           else
               write(c_thr,902)nint(pcp_thr(idbz)*100.)
 902           format(i4.4)
           endif

           plot_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                      //'/plot'
!    1                                      //c_model(1:len_model)
!    1                                      //'/'
           len_plot = len_verif + 5 + lenvar ! + len_model

!          write gnuplot file for the time series of this model (region 1)
           bias_file = plot_dir(1:len_plot)//'/'
     1                                     //trim(c_thr)//'/'
     1                                     //a9time_initial     
     1                                     //'.bias'     

           write(6,*)'bias_file = ',bias_file

           ets_file  = plot_dir(1:len_plot)//'/'
     1                                     //trim(c_thr)//'/'
     1                                     //a9time_initial     
     1                                     //'.ets'      

           write(6,*)'ets_file = ',ets_file

           members_file  = verif_dir(1:len_verif)//'/'
     1                                     //'members.txt'      

           write(6,*)'members_file = ',members_file

           open(lun_bias,file=bias_file,status='unknown')
           open(lun_ets,file=ets_file,status='unknown')

!          write comment with model member names
           write(lun_bias,905)(trim(c_fdda_mdl_src(imodel))
     1                              ,imodel=2,n_fdda_models)  
           write(lun_ets,905)(trim(c_fdda_mdl_src(imodel))
     1                              ,imodel=2,n_fdda_models)  
905        format('# ',30(1x,a))  

           if(l_col)then
!              write bias and ets values
               do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                 ,istatus)

!                  if little 20dbz radar coverage set to missing
                   do imodel = 2,n_fdda_models
                       if(frac_coverage(imodel,itime_fcst,iregion,1)
     1                                                   .lt. 0.001)then
                           write(6,*)' set to missing model/itime/frac '
     1                       ,imodel,itime_fcst
     1                       ,frac_coverage(imodel,itime_fcst,iregion,1)       
                           bias(imodel,itime_fcst,iregion,idbz) = rmiss       
                           ets(imodel,itime_fcst,iregion,idbz)  = rmiss
                       endif
                   enddo ! imodel

                   write(lun_bias,911)a24time_valid,    
     1                 (bias(imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
                   write(lun_ets,911)a24time_valid,    
     1                 (ets(imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     

911                format(a24,3x,40f12.3)
               enddo ! itime_fcst

!              write n values in separate blocks                     
               do jn = 0,1
               do in = 0,1
                 write(lun_bias,*)
                 write(lun_bias,*)
                 do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                 ,istatus)
                   write(lun_bias,912)a24time_valid,    
     1                 (n(imodel,itime_fcst,iregion,idbz,in,jn)
     1                              ,imodel=2,n_fdda_models)     
912                format(a24,3x,40i12.3)
                 enddo ! itime_fcst
               enddo ! jn
               enddo ! in

!              write fractional coverage values
               write(lun_bias,*)
               write(lun_bias,*)
               do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                 ,istatus)
                   write(lun_bias,913)a24time_valid,    
     1                 (frac_coverage                
     1                   (imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
913                format(a24,3x,40f12.5)
               enddo ! itime_fcst

!              write to members.txt file
               open(lun_mem,file=members_file,status='unknown')
               do imodel = 2,n_fdda_models
                   c_model = c_fdda_mdl_src(imodel)
                   if(c_model(1:3) .ne. 'lga')then
                       write(lun_mem,*)trim(c_model)              
                   endif
               enddo ! imodel
               close(lun_mem)

           else ! l_col is false (older code)
             do imodel = 2,n_fdda_models
               do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                 ,istatus)
                   write(lun_bias,921)a24time_valid     
     1                             ,bias(imodel,itime_fcst,iregion,idbz)     
                   write(lun_ets,921)a24time_valid     
     1                              ,ets(imodel,itime_fcst,iregion,idbz)
921                format(a24,3x,2f12.3)
               enddo ! itime_fcst
               write(lun_bias,*)
               write(lun_ets,*)
             enddo ! imodel
           endif

           close(lun_bias)
           close(lun_ets)

         enddo ! idbz

        enddo ! fields

 999    write(6,*)' end of subroutine verif_radar'

        return

        end


        subroutine read_region_info(
     1                        maxbgmodels,max_fcst_times,max_regions
     1                       ,n_models,n_fcst_times,n_regions
     1                       ,il,ih,jl,jh,lun_in)

        integer il(maxbgmodels,0:max_fcst_times,max_regions)
        integer ih(maxbgmodels,0:max_fcst_times,max_regions)
        integer jl(maxbgmodels,0:max_fcst_times,max_regions)
        integer jh(maxbgmodels,0:max_fcst_times,max_regions)

        integer iperim_buff

        character*150 static_dir,static_file

        iperim_buff = 0

!       read in data file with region points
        call get_directory('static',static_dir,len_static)
        static_file = static_dir(1:len_static)//'/verif_regions.txt'
        open(lun_in,file=static_file,status='old')

        read(lun_in,*)n_regions
        read(lun_in,*)n_fcst_times_dum
        read(lun_in,*)n_models_dum

        if(n_regions .eq. 0)then
          write(6,*)' filling default region info automatically'
          n_regions = 1
          call get_grid_dim_xy(nx_l,ny_l,istatus)

          do ir = 1,n_regions
            do if = 0,n_fcst_times
                do im = 1,n_models
                    il(im,if,ir) = 1    + iperim_buff
                    ih(im,if,ir) = nx_l - iperim_buff
                    jl(im,if,ir) = 1    + iperim_buff
                    jh(im,if,ir) = ny_l - iperim_buff
                enddo ! im
            enddo ! if
          enddo ! ir

        else
          do ir = 1,n_regions
            do if = 1,n_fcst_times
                read(lun_in,*)
                read(lun_in,*)i_fcst_time

                if(if .ne. (i_fcst_time+1))then
                    write(6,*)' error in read_region_info '
                    write(6,*)' if differs from i_fcst_time+1: '
     1                        ,if,i_fcst_time,i_fcst_time+1
                    stop
                endif

                do im = 1,n_models
                    read(lun_in,*)il(im,if,ir),jl(im,if,ir)
                    read(lun_in,*)ih(im,if,ir),jh(im,if,ir)
                enddo ! im

            enddo ! if
          enddo ! ir
        endif

        close(lun_in)

        return
        end


        subroutine calc_thresholds(itime_fcst,nthr,dbz_an,dbz_fc
     1                            ,frac_obs,frac_fcst,bias,rmiss)

        real dbz_an(nthr)
        real dbz_fc(nthr)
        real frac_obs(nthr)
        real frac_fcst(nthr)
        real bias(nthr)
        real dbz_nw(nthr)
        integer jbr(nthr)

        jbr_min = nthr+1
        jbr_max = 0

        errmax = 0.

        do idbz = 1,nthr
            if(bias(idbz) .ne. rmiss)then
              if(bias(idbz) .ne. 0.)then
                errmax = max(errmax,abs(alog10(bias(idbz))))
              else
                errmax = 999.
              endif
            endif

!           find bracket values of obs that match fcst
            jbr(idbz) = 0
            do jdbz = 1,nthr-1
                if(frac_fcst(jdbz)   .ge. frac_obs(idbz) .and.
     1             frac_fcst(jdbz)   .ne. rmiss          .and.
     1             frac_fcst(jdbz+1) .le. frac_obs(idbz) .and.
     1             frac_fcst(jdbz+1) .ne. rmiss          .and.         
     1             bias(idbz)        .ne. rmiss                )then
                    jbr(idbz) = jdbz
                    jbr_min = min(jbr_min,idbz)
                    jbr_max = max(jbr_max,idbz)
                endif
            enddo ! jdbz

            dbz_nw(idbz) = rmiss

            if(bias(idbz) .ne. rmiss)then
              if(jbr(idbz) .gt. 0)then
                fracint = (frac_obs(idbz)        -frac_fcst(jbr(idbz)) ) 
     1                  / (frac_fcst(jbr(idbz)+1)-frac_fcst(jbr(idbz)) )       
!               write(6,*)' idbz/jbr/fracint = '
!    1                     ,idbz,jbr,fracint
                dbz_nw(idbz) = dbz_fc(jbr(idbz)) + 
     1              fracint * (dbz_fc(jbr(idbz)+1) - dbz_fc(jbr(idbz)))

              else ! we have to go beyond the end of the binning range
                if(bias(idbz) .lt. 1.0)then
                  dbz_nw(idbz) = dbz_fc(idbz) - 1.0
                else
                  dbz_nw(idbz) = dbz_fc(idbz) + 1.0
                endif

              endif ! jbr(idbz) .ne. 0

            endif ! valid bias value
        enddo ! idbz

        write(6,*)
        write(6,*)' subroutine calc_thresholds, errmax = ',errmax
        write(6,*)' dbz_an = ',dbz_an
        write(6,840)itime_fcst
 840    format('   bin  dbz/an  dbz/fc    bias      obs      ',
     1         'fcst  jbrckt dbz/nw   for itime_fcst: ',i4)       

!       additional passes to get tails
        do idbz = jbr_min+1,nthr
          if(bias(idbz) .ne. rmiss)then
             dbz_nw(idbz) = max(dbz_nw(idbz),dbz_nw(idbz-1)+1.0)
             write(6,*)' upper tail ',idbz,dbz_nw(idbz)
          else
             write(6,*)' upper tail ',idbz
          endif
        enddo ! idbz

        do idbz = jbr_max-1,1,-1
          if(bias(idbz) .ne. rmiss)then
             dbz_nw(idbz) = min(dbz_nw(idbz),dbz_nw(idbz+1)-1.0)
             write(6,*)' lower tail ',idbz,dbz_nw(idbz)
          else
             write(6,*)' lower tail ',idbz
          endif
        enddo ! idbz

        do idbz = 1,nthr
           write(6,850)idbz,dbz_an(idbz),dbz_fc(idbz)
     1                 ,bias(idbz)
     1                 ,frac_obs(idbz)
     1                 ,frac_fcst(idbz)
     1                 ,jbr(idbz)
     1                 ,dbz_nw(idbz)
 850       format(1x,i4,2f8.1,f10.3,2f10.6,i5,f8.1)
        enddo ! idbz

        dbz_fc = dbz_nw

        return
        end
