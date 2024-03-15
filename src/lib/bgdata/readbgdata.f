      subroutine read_bgdata(nx_bg,ny_bg,
     +	   nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,ctype
     +    ,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +    ,prbght,prbgsh,prbguv,prbgww
     +    ,htbg, tpbg,uwbg,vwbg,shbg,cwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc
     +    ,t_at_sfc,uwbg_sfc,vwbg_sfc,mslpbg,pcpbg,crefbg,tpw,cwat,swi
     +    ,istatus)

c kml: changes made april 2004
c tdbg_sfc (model 2m dew point) is now read in during subroutine read_eta_conusc
c tdbg_sfc array is checked for nan
c kml: end

      use mem_namelist, only: precip_cycle_time 

      implicit none
      include 'netcdf.inc'

      integer nx_bg
      integer ny_bg
      integer nx_l
      integer ny_l
      integer nzbg_ht
      integer nzbg_tp
      integer nzbg_sh
      integer nzbg_uv
      integer nzbg_ww
      integer bgmodel
      integer ntbg,nzbg
      integer i,j,k,l
      integer lencm,lend
      integer i4_initial
      integer i4_valid
      integer i4time
      integer i4hr
      integer nan_flag
      integer nf_fid,nf_vid,nf_status
      integer istatus,istat_cyc,laps_cycle_time,ncyc


c *** background model grid data.
c
      real :: prbght(nx_bg,ny_bg,nzbg_ht) !pressure (mb) ht and temp
      real :: prbgsh(nx_bg,ny_bg,nzbg_sh) !pressure (mb) q
      real :: prbguv(nx_bg,ny_bg,nzbg_uv) !pressure (mb) u- v-components
      real :: prbgww(nx_bg,ny_bg,nzbg_ww) !pressure (mb) omega
      real  pr(nzbg_ht)

      real :: htbg(nx_bg,ny_bg,nzbg_ht)   !height (m)
      real :: tpbg(nx_bg,ny_bg,nzbg_tp)   !temperature (k)
      real :: shbg(nx_bg,ny_bg,nzbg_sh)   !specific humidity (kg/kg)
      real :: uwbg(nx_bg,ny_bg,nzbg_uv)   !u-wind (m/s)
      real :: vwbg(nx_bg,ny_bg,nzbg_uv)   !v-wind (m/s)
      real :: cwbg(nx_bg,ny_bg,nzbg_ww)   !w-wind (pa/s)

      real :: mslpbg(nx_bg,ny_bg)         !mslp  (mb)
      real :: htbg_sfc(nx_bg,ny_bg)
      real :: prbg_sfc(nx_bg,ny_bg)
      real :: shbg_sfc(nx_bg,ny_bg)       !specific humidity (kg/kg)
      real :: uwbg_sfc(nx_bg,ny_bg)
      real :: vwbg_sfc(nx_bg,ny_bg)
      real :: tdbg_sfc(nx_bg,ny_bg)
      real :: tpbg_sfc(nx_bg,ny_bg)
      real :: t_at_sfc(nx_bg,ny_bg)       !skin/ground temperature
      real :: pcpbg(nx_bg,ny_bg)          !precip at surface (m)
      real :: crefbg(nx_bg,ny_bg)         !composite reflectivity
      real tpw(nx_bg,ny_bg)
      real cwat(nx_bg,ny_bg)
      real swi(nx_bg,ny_bg)

c     local variables for the time being
      real r01(nx_bg,ny_bg)
      real llr(nx_bg,ny_bg)
      real s8a(nx_bg,ny_bg)

      real      lon0,lat1,lat2
      real      ssh2
      real      r_missing_data
     
      real      argmin,argmax

      character*256   bgpath
      character*256   fname_bg
      character*200   fullname
      character*150   directory
      character*132   cmodel
      character*125   comment_2d
      character*10    units_2d
      character*13    fname13
      character*13    fname9_to_wfo_fname13
      character*6     c6_maproj
      character*5     ctype      !either "dprep" or "lapsb" depending on dprep or lga
      character*4     af_bg
      character*2     gproj

      write(6,*)' subroutine read_bgdata: dims are ',nx_bg,ny_bg,nzbg_ht

      call get_r_missing_data(r_missing_data,istatus)

!     initialize
      htbg_sfc = 0.
      prbg_sfc = r_missing_data
      tdbg_sfc = r_missing_data
      shbg_sfc = r_missing_data
      pcpbg = r_missing_data

      call s_len(cmodel,lencm)

      if(bgmodel.eq.0)then
        if(cmodel(1:lencm).eq.'laps_fua'.or.
     +     cmodel(1:lencm).eq.'laps'.or.
     +     cmodel(1:lencm).eq.'model_fua')then

           call get_directory_length(fullname,lend)
           directory=fullname(1:lend)
           call cv_asc_i4time(fullname(lend+1:lend+9),i4_initial)
           read(af_bg,100)i4hr
100     format(i4.4)
           i4_valid=i4_initial+i4hr/100*3600

c the following subroutine should also work for different
c domain fua/fsf but we'll try the get_lapsdata stuff first.

         if(cmodel(1:lencm).eq.'model_fua')then
            call read_fuafsf_cdf(fullname
     +                          ,nx_bg, ny_bg, nzbg_ht
     +                          ,htbg, pr, cwbg, shbg, tpbg, uwbg, vwbg       
     +                          ,uwbg_sfc, vwbg_sfc, tpbg_sfc, tdbg_sfc       
     +                          ,prbg_sfc, mslpbg, htbg_sfc, r01, pcpbg
     +                          ,crefbg, llr, s8a, swi, tpw
     +                          ,istatus)
            if(istatus.ne.1)then
               print*,'error returned: read_fuafsf_cdf'
               return
            endif
            do j=1,ny_bg
            do i=1,nx_bg
               prbght(i,j,:)=pr(:)
               prbgsh(i,j,:)=pr(:)
               prbguv(i,j,:)=pr(:)
               prbgww(i,j,:)=pr(:)
            enddo
            enddo

         elseif(cmodel.eq.'laps_fua')then

           call get_grid_dim_xy(nx_l,ny_l,istatus)

           if(nx_l .ne. nx_bg .or. ny_l .ne. ny_bg)then
              print*
              print*,' ***********************************'
              print*,' error: nx-laps ne nx_bg '
              print*,' background.nl var cmodel = ',cmodel
              print*,' ***********************************'
              print*
              return
           endif

           call s_len(fullname,lend)
           fullname=fullname(1:lend)//".fua"
           nf_status = nf_open(fullname,nf_nowrite,nf_fid)
           if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'nf_open ', fullname
              return
           endif

c        nf_status=nf_inq_varid(nf_fid,'z',nf_vid)
c        if(nf_status.ne.nf_noerr) then
c           print *, nf_strerror(nf_status)
c           print *,'in nf_get_var_ model '
c           return
c        endif
c        nf_status = nf_inq_dimlen(nf_fid,nf_vid,nzbg)
c        if(nf_status.ne.nf_noerr) then
c           print *, nf_strerror(nf_status)
c           print *,'dim n_valtimes'
c           return
c        endif

           nf_status=nf_inq_varid(nf_fid,'level',nf_vid)
           if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'in nf_get_var_ model '
              return
           endif
           nf_status=nf_get_var_real(nf_fid,nf_vid,pr)
           if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'in nf_get_var_ model '
              return
           endif

           do k = 1,nzbg_ht
           do j=1,ny_bg
           do i=1,nx_bg
              prbght(i,j,k)=pr(nzbg_ht-k+1)
              prbgsh(i,j,k)=pr(nzbg_ht-k+1)
              prbguv(i,j,k)=pr(nzbg_ht-k+1)
              prbgww(i,j,k)=pr(nzbg_ht-k+1)
           enddo
           enddo
           enddo

c upper air
           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_uv,directory,'u3 '
     1            ,units_2d,comment_2d,uwbg,istatus)
           if(istatus.ne.1)then
              print*,'error 3d bkgd file (u3): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_uv,directory,'v3 '
     1            ,units_2d,comment_2d,vwbg,istatus)
           if(istatus.ne.1)then
              print*,'error 3d bkgd file (v3): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ht,directory,'t3 '
     1            ,units_2d,comment_2d,tpbg,istatus)
           if(istatus.ne.1)then
              print*,'error 3d bkgd file (t3): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ht,directory,'ht '
     1            ,units_2d,comment_2d,htbg,istatus)
           if(istatus.ne.1)then
              print*,'error 3d bkgd file (ht): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_sh,directory,'sh '
     1            ,units_2d,comment_2d,shbg,istatus)
           if(istatus.ne.1)then
              print*,'error 3d bkgd file (sh): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ww,directory,'om '
     1            ,units_2d,comment_2d,cwbg,istatus)

           if(istatus.ne.1)then
              print*,'error 3d bkgd file (om): ',directory(1:lend)
              return
           endif
c sfc data
           search_dir: do l=lend,1,-1
              if(directory(l:l).eq.'f')then
                 if(directory(l:l+2).eq.'fua')then
                    exit search_dir
                 endif
              endif
           enddo search_dir

           if(l.le.1)then
              print*,'unable to determine location of fua in'
              print*,'directory string for fsf. return with no data'
              return
           endif

           directory(l:l+2)='fsf'

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'psf',units_2d,comment_2d,nx_bg,ny_bg,prbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'error 2d bkgd file (psf): ',directory(1:lend)
              return
           endif

c          prbg_sfc=prbg_sfc*100.

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'tsf',units_2d,comment_2d,nx_bg,ny_bg,tpbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'error 2d bkgd file (tsf): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'slp',units_2d,comment_2d,nx_bg,ny_bg,mslpbg
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'error 2d bkgd file (slp): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'dsf',units_2d,comment_2d,nx_bg,ny_bg,shbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'error 2d bkgd file (dsf): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'usf',units_2d,comment_2d,nx_bg,ny_bg,uwbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'error 2d bkgd file (usf): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'vsf',units_2d,comment_2d,nx_bg,ny_bg,vwbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'error 2d bkgd file (vsf): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'ter',units_2d,comment_2d,nx_bg,ny_bg,htbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'error 2d bkgd file (ter): ',directory(1:lend)
              return
           endif
c
c ---------------read laps analyses---------------------
c
         elseif(cmodel.eq.'laps')then

           call get_laps_3d_analysis_data(i4_initial,nx_bg,ny_bg
     +,nzbg_ht, htbg,tpbg,uwbg,vwbg,shbg,cwbg,istatus)

           call get_laps_2d(i4_initial,'lsx','ps ',units_2d,comment_2d
     +,nx_bg,ny_bg,prbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','u  ',units_2d,comment_2d
     +,nx_bg,ny_bg,uwbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','v  ',units_2d,comment_2d
     +,nx_bg,ny_bg,vwbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','t  ',units_2d,comment_2d
     +,nx_bg,ny_bg,tpbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','td ',units_2d,comment_2d
     +,nx_bg,ny_bg,shbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','msl',units_2d,comment_2d
     +,nx_bg,ny_bg,mslpbg,istatus)

           if(istatus.ne.1)then
              print*,'error returned: read_laps_analysis'
              return
           endif
           
         endif !model_fua?!

        endif

      elseif (bgmodel .eq. 1) then     ! process 60 km ruc data

          call read_ruc60_native(bgpath,fname_bg,af_bg,nx_bg,ny_bg
     .,nzbg_ht,prbght,htbg,tpbg,shbg,uwbg,vwbg,gproj,lon0,istatus)

      elseif (bgmodel .eq. 2) then ! process 48 km eta conus-c grid data
c
c for now all fields have 3d dimension of nzbg_ht
c
          call read_eta_conusc(fullname,nx_bg,ny_bg,nzbg_ht
     .,htbg,prbght,tpbg,uwbg,vwbg,shbg,cwbg
     .,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc
     .,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

          if(istatus.ne.0)goto 99

          if(ctype.eq."lapsb")then

c convert rh to sh.
             print*,'prepare grids for laps background-lga'
             call lprep_eta_conusc(nx_bg,ny_bg,nzbg_ht
     +,prbght,tpbg,shbg,tpbg_sfc,prbg_sfc,shbg_sfc,istatus)
 
          endif
          if(cmodel(1:lencm).eq.'orsm_hko')then
             print*,'in read_bgdata'
             print*,'convert ww to pa/sec for orsm_hko'
             cwbg=cwbg/36.
          endif

          prbgsh=prbght
          prbguv=prbght
          prbgww=prbght
c
      elseif (bgmodel .eq. 4) then

          call read_sbn_grids(fullname,af_bg,cmodel
     .,nx_bg,ny_bg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .,prbght,prbgsh,prbguv,prbgww
     .,htbg,tpbg,shbg,uwbg,vwbg,cwbg
     .,htbg_sfc,prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc
     .,tpbg_sfc,mslpbg,ctype,istatus)

          if(istatus.ne.0)goto 99

c original nz allows reading of sfc variables (ie., 2m and 10m)
c but now decrement nz since arrays only contain 3d info.

c         if(cmodel(1:lencm).eq.'ruc40_native')then
c            nzbg_sh=nzbg_sh-1
c            nzbg_uv=nzbg_uv-1
c         elseif(cmodel(1:lencm).eq.'eta48_conus')then
c            nzbg_sh=nzbg_sh-1
c            nzbg_uv=nzbg_uv-1
c         elseif(cmodel(1:lencm).eq.'avn_sbn_cyleq')then
c            nzbg_sh=nzbg_sh-2
c            nzbg_uv=nzbg_uv-2
c         endif
 
      elseif (bgmodel .eq. 5) then ! process 40km ruc public-netcdf data

          call read_ruc2_hybb(fullname,nx_bg,ny_bg,nzbg_ht
     +,mslpbg,htbg,prbght,shbg,uwbg,vwbg,tpbg,cwbg,istatus)

          if(istatus.ne.0)goto 99 

          print*,'read complete'

          if(ctype.eq.'lapsb')then

             print*,' entering lprep_ruc2_hybrid'

             call lprep_ruc2_hybrid(nx_bg,ny_bg,nzbg_ht
     +,htbg,prbght,shbg,uwbg,vwbg,tpbg,uwbg_sfc,vwbg_sfc
     +,tpbg_sfc,prbg_sfc,shbg_sfc,htbg_sfc,istatus)
             print*,'data prep complete'

          endif

          prbgsh=prbght
          prbguv=prbght
          prbgww=prbght
c
c eta grib ingest currently disabled (j. smart 9-4-98)
c also, nogaps 2.5 degree obsolete.
c bgmodel 3 = fa (taiwan). bgmodel 6 = nogaps1.0. bgmodel 8 = avn 1.0 deg
c
      elseif (bgmodel .eq. 3 .or.
     .        bgmodel .eq. 6 .or.
     .        bgmodel .eq. 8 .or.
     .        bgmodel .eq.12) then ! process avn, ecmwf or nogaps1.0 grib data

             call read_dgprep(bgmodel,cmodel,bgpath
     .,fname_bg,af_bg,nx_bg,ny_bg,nzbg_ht
     .,prbght,htbg,tpbg,shbg,uwbg,vwbg,cwbg
     .,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc,t_at_sfc
     .,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

             if(istatus.ne.0)goto 99

             if(bgmodel.eq.12.and.
     &          cmodel(1:lencm).eq.'fmi_netcdf_ll')then
                tdbg_sfc=shbg_sfc
             endif

             prbgsh=prbght
             prbguv=prbght
             prbgww=prbght

      elseif (bgmodel .eq. 9) then ! process nws conus data (ruc,eta,ngm,avn)

             call read_conus_nws(bgpath,fname_bg,af_bg
     .,nx_bg,ny_bg,nzbg_ht,prbght
     .,htbg,tpbg,shbg,uwbg,vwbg
     .,gproj,lon0,lat1,lat2,istatus)
c
c wni-bls
      elseif (bgmodel .eq. 10) then ! process unidata netcdf
       if ( (cmodel .eq. 'ruc_iso') .or.
     +      (cmodel .eq. 'gfs_iso') )then
          call read_unidata_iso(fullname,af_bg,cmodel
     .   ,nx_bg,ny_bg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .   ,prbght,prbgsh,prbguv,prbgww
     .   ,htbg,tpbg,shbg,uwbg,vwbg,cwbg
     .   ,htbg_sfc,prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc
     .   ,tpbg_sfc,mslpbg,ctype,istatus)
       elseif(cmodel .eq. 'ruc_hyb') then
           call read_unidata_ruc_hyb(fullname,af_bg,cmodel
     .   ,nx_bg,ny_bg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .   ,prbght,prbgsh,prbguv,prbgww
     .   ,htbg,tpbg,shbg,uwbg,vwbg,cwbg
     .   ,htbg_sfc,prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc
     .   ,tpbg_sfc,mslpbg,ctype,istatus)
         print *, "completed read of ruc_hyb"
       endif

      elseif (bgmodel .eq. 13) then ! process grib1/grib2

         write(*,*) 'call degrib_data: dims are ',nx_bg, ny_bg, nzbg_ht
         write(*,*) ' grib filename ',trim(fullname)

         call degrib_data(fullname, nx_bg, ny_bg, nzbg_ht, 
     &      prbght, htbg, tpbg, shbg, uwbg, vwbg, cwbg, 
     &      htbg_sfc, tpbg_sfc, shbg_sfc, uwbg_sfc, vwbg_sfc, 
     &      tdbg_sfc, t_at_sfc, prbg_sfc, mslpbg, pcpbg, crefbg, 
     &      tpw,cwat,istatus)

            prbgsh(:,:,:)=prbght(:,:,:) 
            prbguv(:,:,:)=prbght(:,:,:) 
            prbgww(:,:,:)=prbght(:,:,:) 

c           write(*, *) "readbgdata htbg(3,30,1)", htbg(3,30,1)
c           write(*, *) "readbgdata tpbg(3,30,1)", tpbg(3,30,1)
c           write(*, *) "readbgdata cwbg(3,30,1)", cwbg(3,30,1)
c           write(*, *) "readbgdata: shbg_sfc/tdbg_sfc is actually rh?"
            write(*, *) "readbgdata shbg_sfc(3,30)", shbg_sfc(3,30)
            write(*, *) "readbgdata tdbg_sfc(3,30)", tdbg_sfc(3,30)
            write(*, *) "readbgdata pcpbg(3,30)", pcpbg(3,30)

            if(precip_cycle_time .eq. 900)then
               write(6,*)' dividing precip based on cycle time'
               ncyc = 3600 / laps_cycle_time
               where (pcpbg(:,:) .ne. r_missing_data)
                  pcpbg(:,:) = pcpbg(:,:) / float(ncyc)
               endwhere
            endif

            write(6,*)' setting a precip floor value of zero'
            where (pcpbg(:,:) .ne. r_missing_data)
               pcpbg(:,:) = max(pcpbg(:,:),0.)
            endwhere
      endif
c      
c - 3d fields
c
      call check_nan3 (htbg,nx_bg,ny_bg,nzbg_ht,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in height array '
c        goto 99
      endif

      call check_nan3 (tpbg,nx_bg,ny_bg,nzbg_tp,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in temperature array '
c        goto 99
      endif

      call check_nan3 (shbg,nx_bg,ny_bg,nzbg_sh,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in specif hum array '
c        goto 99
      endif

      call check_nan3 (uwbg,nx_bg,ny_bg,nzbg_uv,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in u-comp array '
c        goto 99
      endif

      call check_nan3 (vwbg,nx_bg,ny_bg,nzbg_uv,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in v-comp array '
c        goto 99
      endif

      call check_nan3 (cwbg,nx_bg,ny_bg,nzbg_ww,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in w-comp array '
      endif
c
c - 2d fields
c
      call check_nan2(htbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc height array '
      endif

      call check_nan2(prbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc pressure array '
      endif

      call check_nan2(tdbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc temp array '
      endif

      call check_nan2(tpbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc temp array '
      endif

      call check_nan2(shbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc spec hum array '
      endif

      call check_nan2(uwbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc u-comp array '
      endif

      call check_nan2(vwbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc v-comp array '
      endif

      call check_nan2(mslpbg,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc mslp array '
      endif

      call check_nan2(pcpbg,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' error: nan found in sfc pcpbg array '
      endif

      argmin = minval(prbg_sfc)
      argmax = maxval(prbg_sfc)
      write(6,*)' prbg_sfc range = ',argmin,argmax
      if(argmax .lt. 100. .or. argmax .gt. 1e6)then
          write(6,*)' warning: prbg_sfc has questionable range'
      endif

      write(6,*)' tdbg_sfc range = ',minval(tdbg_sfc),maxval(tdbg_sfc)
      write(6,*)' shbg_sfc range = ',minval(shbg_sfc),maxval(shbg_sfc)
      write(6,*)' pcpbg range = ',minval(pcpbg),maxval(pcpbg)
      write(6,*)' returning from read_bgdata'

      istatus = 0

99    if(istatus.ne. 0)then
         print*,'error with background model data in read_bgdata'
      endif
      return
      end
