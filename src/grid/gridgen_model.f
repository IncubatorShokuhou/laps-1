cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

	program gridgen_model
c**********************************************************************
c portions of the following code were taken from the rams software and 
c   were used by permission of the *aster division/mission research 
c   corporation. 
c**********************************************************************
c*	this program will be used to map model grids based on the rams*
c*  version 2b which uses polar stereographic projections. other      *
c*  projections have since been added.                                *
c*                                                                    *
c*********************************************************************
c beginning july 2000, this software was updated to accomodate the    *
c wrfsi (weather research and forecasting standard initialization)    *
c grid requirements.  this included adding stagger calculations, and  *
c map projection factors (map factor and components of coriolis parm, *
c landuse, soil type, vegetation greenness fraction, mean annual soil *
c temperature, albedo, and sea surface temperature.  significant work *
c has occurred to reduce memory requirements.  many sections of code  *
c previously in this file now reside in file gridgen_utils.f          *
c for further information contact                                     *
c mr. john smart noaa/fsl, smart@fsl.noaa.gov (303-497-6590)          *
c*********************************************************************

        use horiz_interp

        integer istag
        integer n_staggers
        integer lf
        character c6_maproj*6
        character c8_maproj*8
        character c10_grid_f*10 ! type of domain (nest7grid, wrfsi)/
        character c_dataroot*200
        character cstaticdir*200
        logical   localize
        real,     allocatable  ::  mother_lat(:,:)
        real,     allocatable  ::  mother_lon(:,:)
        real,     allocatable  ::  dum2d(:,:)
        real      la1,lo1,la2,lo2
        real      lov,latin1,latin2
        real      dx,dy

        integer lli_orig
        integer llj_orig
        integer uri_orig
        integer urj_orig

        include 'grid_fname.cmn'
 
        write(6,*)
        write(6,*)' gridgen_model: start'

        call get_config(istatus)
        if(istatus.ne.1)then
           print*,'error returned: get_config'
           print*,'abort'
           stop
        endif

        call find_domain_name(c_dataroot,c10_grid_f,istatus)
        if(istatus .ne. 1)then
            write(6,*) 'error: find domain name - stop'
            stop
        endif
        call s_len(c10_grid_f,lf)

        if(c10_grid_f(1:lf).eq.'wrfsi')then
           call get_num_domains(num_domains,istatus)
           if(istatus.ne.1)then
              print *,'error returned: get_num_domains'
              print *,'abort'
              stop
           endif
        else
           num_domains=1
        endif
c -----------------------------------------------------------------
c variable "nest" is common in src/include/grid_fname.cmn
c and is used to communicate the domain # throughout the code.
c -----------------------------------------------------------------
        do nest=1,num_domains
c
c determine if this domain requires localizing or re-localizing
c
           call get_directory('static',cstaticdir,len_dir)

           call rd_static_attr(cstaticdir,nest
     .,c10_grid_f,nx_dom,ny_dom,dx,dy,la1,lo1,latin1
     .,latin2,lov,c8_maproj,istatus)

c -------------------------------------------------------------------------
c routine below accesses environment variable $force_localization which, if
c set, will return variable "localize" = .true.
c -------------------------------------------------------------------------
 
           call eval_localization(cstaticdir,nest,localize
     .,c10_grid_f,la1,lo1,istatus)

           call get_parent_id(iparent_id,istatus)

           if(localize)then

c get_mother_dims,  uses parent_id as index to determine mother domain specs
           if(nest.gt.1)then
              call get_mother_dims(nx_dom,ny_dom,istatus)
              allocate (mother_lat(nx_dom,ny_dom),
     +                  mother_lon(nx_dom,ny_dom),
     +                  dum2d(nx_dom,ny_dom),stat=istat)
              if(istat.ne.0)then
                 print*,'cannot deallocate mother arrays: ',istat
                 stop
              endif

              call get_mother_domain(nx_dom,ny_dom,mother_lat
     +,mother_lon,istatus)
              if(istatus.lt.1)then
                 print*
                 print *,'error reading mother lat, lon data.'
                 print *,'abort'
                 print*
                 stop
              endif
           endif   

           call get_grid_dim_xy(nx_dom,ny_dom,istatus)
           call get_ratio_to_parent(iratio_to_parent,istatus)
           call get_domain_origin(lli_orig,llj_orig,uri_orig,urj_orig
     +,istatus)
c
c use nx_dom and ny_dom generically regardless of a nest or not.
c
           call get_grid_spacing(dx,istatus)
           call get_c6_maproj(c6_maproj,istatus)
           call downcase(c6_maproj,c6_maproj)

           if(c6_maproj.ne.'rotlat')dy=dx

           print*
           print*,'---------------------------------------------'
           print*,'domain #         = ',nest
           print*,'parent of nest   = ',iparent_id
           print*,'ratio to parent  = ',iratio_to_parent
           print*,'grid_spacing(m)  = ',dx
           print*,'grid dims (nx/ny)= ',nx_dom,ny_dom
           print*,'lower left  i/j  = ',lli_orig,llj_orig
           print*,'upper right i/j  = ',uri_orig,urj_orig

           if(nest.gt.1)then
              print*,'calculate center lat/lon '
              la1=mother_lat(lli_orig,llj_orig)
              lo1=mother_lon(lli_orig,llj_orig)
              la2=mother_lat(uri_orig,urj_orig)
              lo2=mother_lon(uri_orig,urj_orig)
              print*,'*** sw lat/lon = ',la1,lo1,'***'
              print*,'*** ne lat/lon = ',la2,lo2,'***'
              call get_earth_radius(erad,istatus)
c get x/y for nested grid center
              call latlon_to_xy(la1,lo1,erad,xsw,ysw) 
              call latlon_to_xy(la2,lo2,erad,xne,yne) 
              xcen=(xne+xsw)/2.
              ycen=(yne+ysw)/2.
              call xy_to_latlon(xcen,ycen,erad,rlatcen,rloncen)
              print*,'*** center lat/lon= ',rlatcen,rloncen,' ***'
           endif

           call get_n_staggers(n_staggers,istatus)
           call get_stagger_index(istag,istatus)

           if(istatus.ne.1)then
              print*,'error getting namelist info for',
     +': ',c10_grid_f(1:lf)
              stop
           endif

           if(allocated(mother_lat))then
              deallocate (mother_lat, mother_lon, dum2d)
           endif

           print*
c
c -------------------------------------------------------------------
c
           call gridmap_sub(nx_dom,ny_dom,n_staggers,istag,nest
     +,rlatcen,rloncen,lli_orig,llj_orig,uri_orig,urj_orig
     +,iparent_id,iratio_to_parent,istatus)
c
c -------------------------------------------------------------------
c
           print*,' gridmap_sub finished '
           print*
           if(istatus.eq.1)then

              print*,' gridmap_sub finish successfully'

              if( allocated(mother_lat) )then
               deallocate (mother_lat,mother_lon,dum2d,stat=istat)
               if(istat.ne.0)then
                  print*
                  print*,'cannot deallocate mother arrays: ',istat
                  print*,'abort'
                  print*
                  stop
               endif
              endif

           else
              print*,'problem in gridmap_sub. terminating'
              stop
           endif

           else  !localize or not!

           print*
           print*,' no need to re-localize this domain '
           print*,' localization namelist parms havent changed'
           print*

           if( allocated(mother_lat) )then
              deallocate (mother_lat,mother_lon,dum2d,stat=istat)
              if(istat.ne.0)then
                 print*
                 print*,'cannot deallocate mother arrays: ',istat
                 print*,'abort'
                 print*
                 stop
              endif
           endif

           endif !localize case

        enddo

 999	print*,' gridgen_model finish: istatus = ',istatus
        print*

        stop
        end
c
c --------------------------------------------------------------------
c
        subroutine gridmap_sub(nnxp,nnyp,n_staggers,istag,nest
     +,mdlat,mdlon,lli,llj,iuri,iurj,iparent_id,iratio_to_parent
     +,istatus)

        include 'trigd.inc'

        use mem_namelist, only: laps_cycle_time, l_fsf_gridgen

        logical exist,new_dem
        logical lforce_ter_2_zero
!mp
	logical categorical, useland
!mp
        logical l_topo_wps, l_parse, l_topo_1s, l_albedo_bm /.false./

        integer nnxp,nnyp,mode
        integer ngrids
        integer n_staggers
        integer nest
        integer iter
        integer istag

        real  mdlat,mdlon
        real  xmn(nnxp),ymn(nnyp)
        real  xtn(nnxp,n_staggers)
        real  ytn(nnyp,n_staggers)

        real, allocatable ::  topt_10   (:,:)
        real, allocatable ::  topt_10_s (:,:)
        real, allocatable ::  topt_10_ln(:,:)
        real, allocatable ::  topt_10_lt(:,:)
        real, allocatable ::  topt_30   (:,:)
        real, allocatable ::  topt_30_s (:,:)
        real, allocatable ::  topt_30_ln(:,:)
        real, allocatable ::  topt_30_lt(:,:)

        real  topt_out   (nnxp,nnyp)
        real  topt_out_s (nnxp,nnyp)
        real  topt_out_ln(nnxp,nnyp)
        real  topt_out_lt(nnxp,nnyp)

        integer maxdatacat
        parameter (maxdatacat=24)

        real, allocatable ::  adum2d  (:,:)
        real, allocatable ::  adum3d  (:,:,:)

        real, allocatable ::  geodat2d(:,:)
        real, allocatable ::  geodat3d(:,:,:)

        real lats(nnxp,nnyp,n_staggers)
        real lons(nnxp,nnyp,n_staggers)

	real*8 r8lat(nnxp,nnyp), r8lon(nnxp,nnyp)

	real, allocatable, dimension(:,:):: hlat,hlon,vlat,vlon

        character (len=3),   allocatable :: var(:)
        character (len=125), allocatable :: comment(:)

        character*2   cnest
        character*30  ctopo,clatlon,corners,ctopog
        character*30  ctopos,clatlons,cornerss,ctopogs
        character*30  clatlon2d,clatlons2d
        character*131 model

        character*200 path_to_topt30s
        character*200 path_to_topt10m
        character*200 path_to_pctl10m
        character*200 path_to_soiltype_top_30s
        character*200 path_to_soiltype_bot_30s
        character*200 path_to_green_frac    !no reference to res since there is both
                                            !10m and 8.64m
        character*200 path_to_soiltemp_1deg
        character*200 path_to_luse_30s
        character*200 path_to_albedo        !albedo = 0.144 deg res or 8.64m 
        character*200 path_to_maxsnoalb
        character*200 path_to_islope

        character*255 filename
        character*255 filename_wps
        character*200 c_dataroot
        character*200 cdl_dir
        character*180 static_dir 
        character*10  c10_grid_f           ! type of domain (nest7grid, wrfsi)
        character*10  c10_grid_fname       ! actual filename, cdl (nest7grid for now)
        character*6   c6_maproj
        character*1   cdatatype
        integer len,lf,lfn,ns,avgelem,zinelem
        integer ishow_timer,init_timer
        integer itstatus

        real,  allocatable ::  rmapdata(:,:,:)
        real,  allocatable ::  data(:,:,:)     !primary output array storage unit.
 
        interface

          subroutine adjust_geog(nnxp,nnyp,ncat,ctype
     &,istat_dat,lat,topt_out,landmask,geog_data,istatus)

          integer nnxp,nnyp
          integer ncat
          integer istat_dat
          integer istatus

          character*(*) ctype
          real    lat(nnxp,nnyp)
          real    landmask(nnxp,nnyp)
          real    topt_out(nnxp,nnyp)
          real    geog_data(nnxp,nnyp,ncat)
          end subroutine

          subroutine proc_geodat(nx_dom,ny_dom,ncat
     1,path_to_tile_data,dom_lats_in,dom_lons_in,lmask_out
     1,geodat,istatus)

          integer nx_dom
          integer ny_dom
          integer ncat
          character*(*) path_to_tile_data
          real    dom_lats_in(nx_dom,ny_dom)
          real    dom_lons_in(nx_dom,ny_dom)
          real    lmask_out(nx_dom,ny_dom)
          real    geodat(nx_dom,ny_dom,ncat)
          integer istatus
          end subroutine

        end interface

c*********************************************************************

        itstatus = init_timer()

        call find_domain_name(c_dataroot,c10_grid_f,istatus)
        if(istatus .ne. 1)then
            write(6,*) 'error: returned fro find_domain_name '
            return
        endif

        call s_len(c10_grid_f,lf)

!mp  - need to know map projection information from the start
        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error calling laps routine'
            return
        endif
        write(6,*)' c6_maproj = ',c6_maproj
!mp

c add 12 for albedo month 1 - 12.
c add 1 more for islope js: 2-25-03
        if(c10_grid_f(1:lf).eq.'wrfsi')then
!mp
            if (c6_maproj .eq. 'rotlat') then
               write(6,*) 'ngrids set to 103'
               ngrids=103 !!!
            else !normal
               ngrids=112     !97+12   ! 2d grids (including %dist for landuse
                                       ! and two soiltype categories, green frac, albedo and others).
	    endif

        else
           ngrids=38      !26+12   doesn't include terrain slope index and other stuff.
        endif

        allocate (data(nnxp,nnyp,ngrids),stat=istat)
        
        if(istat.ne.0)then
           print*,'error allocating data object: array data ',istat
           print*,'terminating: maybe not enough memory'
           stop
        endif

        model = 'model 4 delta x smoothed filter\0'

        icount_10 = 0
        icount_30 = 0
        icount_ramp = 0

        call get_directory('static',static_dir,lens)

c iplttopo is 1 if you want to plot the topography
c   the 30s topo covers the continental us
cc        itoptfn_30=static_dir(1:len)//'model/topo_30s/u'
c   the 10m topo covers the world
cc        itoptfn_10=static_dir(1:len)//'model/topo_10m/h'

        call get_path_to_topo_10m(path_to_topt10m,istatus)
        if(istatus .ne. 1)then
           write(6,*) 'error getting path_to_topt10m'
           return
        endif

        call get_path_to_topo_30s(path_to_topt30s,istatus)
        if(istatus .ne. 1)then
           write(6,*) 'error getting path_to_topt30s'
           return
        endif

!       l_topo_wps = l_parse(path_to_topt30s,'wps')
        l_topo_1s = l_parse(path_to_topt30s,'topo_1s')

        call get_path_to_soiltype_top(path_to_soiltype_top_30s
     +,istatus)
        if(istatus .ne. 1)then
           write(6,*)'error getting path_to_soiltype_top_30s'
           return
        endif

        call get_path_to_soiltype_bot(path_to_soiltype_bot_30s
     +,istatus)
        if(istatus .ne. 1)then
           write(6,*)'error getting path_to_soiltype_bot_30s'
           return
        endif

        call get_path_to_landuse_30s(path_to_luse_30s,istatus)
        if(istatus .ne. 1)then
           write(6,*)'error getting path_to_luse_30s'
           return
        endif

        call get_path_to_green_frac(path_to_green_frac,istatus)
        if(istatus .ne. 1)then
           print*,'error getting path_to_green_frac'
           return
        endif

        call get_path_to_soiltemp_1deg(path_to_soiltemp_1deg
     &,istatus)
        if(istatus .ne. 1)then
           print*, 'error getting path_to_soiltemp_1deg'
           return
        endif

        call get_path_to_albedo(path_to_albedo,istatus)
        if(istatus .ne. 1)then
           print*, 'error getting path_to_albedo'
           return
        endif

        call get_path_to_maxsnoalb(path_to_maxsnoalb,istatus)
        if(istatus .ne. 1)then
           print*, 'error getting path_to_maxsnoalb'
           return
        endif

        call get_path_to_islope(path_to_islope,istatus)
        if(istatus .ne. 1)then
           print*, 'error getting path_to_islope'
           return
        endif

        call s_len(path_to_topt10m,len)
        if(len.gt.0)then
           print*,'path to toptl0m:        ',path_to_topt10m(1:len)
        endif
        path_to_topt10m(len+1:len+2)='/h'

        call s_len(path_to_soiltype_top_30s,len)
        print*,'path to soiltype_top:   '
     .,path_to_soiltype_top_30s(1:len)
        path_to_soiltype_top_30s(len+1:len+2)='/o'
        print*,'path to soiltype_bot:   '
     .,path_to_soiltype_bot_30s(1:len)
        path_to_soiltype_bot_30s(len+1:len+2)='/o'

        call s_len(path_to_luse_30s,len)
        print*,'path to landuse 30s: ',path_to_luse_30s(1:len)
        path_to_luse_30s(len+1:len+2)='/v'

        call s_len(path_to_green_frac,len)
        print*,'path to green frac:     ',path_to_green_frac(1:len)
        path_to_green_frac(len+1:len+2)='/g'

        call s_len(path_to_soiltemp_1deg,len)
        print*,'path to soiltemp 1deg:  ',path_to_soiltemp_1deg(1:len)
        path_to_soiltemp_1deg(len+1:len+2)='/t'

        call s_len(path_to_albedo,len)
        print*,'path to albedo:         ',path_to_albedo(1:len)
        path_to_albedo(len+1:len+2)='/a'

        call s_len(path_to_maxsnoalb,len)
        print*,'path to max snow albedo: ',path_to_maxsnoalb(1:len)
        path_to_maxsnoalb(len+1:len+2)='/m'

        call s_len(path_to_islope,len)
        print*,'path to islope categorical: ',path_to_islope(1:len)
        path_to_islope(len+1:len+2)='/i'

        call get_topo_parms(silavwt_parm,toptwvl_parm,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'error getting terrain smoothing parms'
	   return
	endif

        call get_r_missing_data(r_missing_data,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'error getting r_missing_data'
	   return
	endif

!       silhouette weighting parameter
        silavwt=silavwt_parm

!       terrain wavelength for filtering
        toptwvl=toptwvl_parm

        iplttopo=1

c*********************************************************************
        call get_gridnl(mode) 

c calculate delta x and delta y using grid and map projection parameters
        call get_standard_latitudes(std_lat,std_lat2,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error calling laps routine'
            return
        endif
        write(6,*)' standard_lats = ',std_lat,std_lat2

	if (c6_maproj .ne. 'rotlat') then

          call get_grid_spacing(grid_spacing_m,istatus)
          if(istatus .ne. 1)then
            write(6,*)' error calling laps routine'
            return
          endif

	else
	
          write(6,*)
     1  'error: get_nmm_grid_spacing not currently supported for rotlat'       
          stop
!         call get_nmm_grid_spacing(dlmd,dphd,grid_spacing_m)

	endif

        write(6,*)' grid_spacing = ',grid_spacing_m

        if(nest.eq.1)then
           print*,'get grid center from namelist'
           call get_grid_center(mdlat,mdlon,istatus)
           if(istatus .ne. 1)then
              write(6,*)' error calling laps routine'
              return
           endif
        endif
          
        print*,' gridmap_sub: grid_center = ',mdlat,mdlon

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error calling laps routine'
            return
        endif
        write(6,*)' c6_maproj = ',c6_maproj

        call get_standard_longitude(std_lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error calling laps routine'
            return
        endif
        write(6,*)' std_lon = ',std_lon

        call get_earth_radius(erad,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error calling laps routine'
            return
        endif
        write(6,*)' earth radius = ',erad

        if(c6_maproj .eq. 'plrstr')then
            call get_ps_parms(std_lat,std_lat2,grid_spacing_m,phi0
     1                       ,grid_spacing_proj_m)

            deltax = grid_spacing_proj_m

            if(phi0 .eq. 90.)then
                write(6,*)' projection is tangent to earths surface'
            else
                write(6,*)' projection is secant to earths surface'
            endif

            if(std_lat2 .eq. +90.)then
                write(6,*)' note, grid spacing will equal '
     1                    ,grid_spacing_m,' at a latitude of ',std_lat
!               write(6,*)' deltax, deltay ',deltax,deltay
!    1                   ,' at the north pole'

            elseif(std_lat2 .eq. -90.)then
                write(6,*)' note, grid spacing will equal '
     1                    ,grid_spacing_m,' at a latitude of ',std_lat       
!               write(6,*)' deltax, deltay ',deltax,deltay
!    1                   ,' at the south pole'

            else ! abs(std_lat2) .ne. 90. (local stereographic)
                write(6,*)' the standard latitude ',std_lat,' is'
     1                   ,' relative to where the pole'
                write(6,*)' of the map projection is: lat/lon '
     1                   ,std_lat2,std_lon
!               write(6,*)' deltax, deltay ',deltax,deltay
!    1                   ,' at the projection pole'
                if(std_lat .ne. +90.)then
                    write(6,*)' note: std_lat should usually be set'
     1                       ,' to +90. for local stereographic'
                endif

            endif

        else ! c6_maproj .ne. 'plrstr'
            deltax = grid_spacing_m

        endif ! c6_maproj .eq. 'plrstr'

        deltay = deltax
        write(6,*)' deltax, deltay (in projection plane) '
     1             ,deltax,deltay
 
c*********************************************************************
c in arrays lats/lons, the first stagger is actually not staggered
c but is the usual lat/lon values at non-staggered grid points. these
c are used for the analysis (laps-nest7grid). for wrfsi, we use the
c mass c-stagger (#4); staggered in both x and y. 

!mp


	if (c6_maproj .eq. 'rotlat') then
          write(6,*) 'dlmd, dphd, mdlat, mdlon: ',dlmd,dphd,mdlat,mdlon       
                                                                                
          allocate(hlat(nnxp,nnyp),hlon(nnxp,nnyp))
          allocate(vlat(nnxp,nnyp),vlon(nnxp,nnyp))
             
          call etall(nnxp,nnyp,mdlat,mdlon,dlmd,dphd,hlat,hlon
     &                  ,vlat,vlon)
             
	  write(6,*) 'egrid corners'
          write(6,*) 'hlat(1,1),hlon(1,1): ', hlat(1,1),hlon(1,1)
          write(6,*) 'hlat(1,jm),hlon(1,jm): ',hlat(1,nnyp),hlon(1,nnyp)       
          write(6,*) 'hlat(im,1),hlon(im,1): ', hlat(nnxp,1),
     &			hlon(nnxp,1)
          write(6,*) 'hlat(im,jm),hlon(im,jm): ', hlat(nnxp,nnyp),
     &			hlon(nnxp,nnyp)
             
             
          do j=1,nnyp
            do i=1,nnxp
                data(i,j,1)=hlat(i,j)
                data(i,j,2)=hlon(i,j)
                data(i,j,3)=vlat(i,j)
                data(i,j,4)=vlon(i,j)
            enddo
          enddo

	elseif(c6_maproj .eq. 'icshdr')then ! read icosahedral grid from a file
          write(6,*)' read icosahedral grid from a file'

          filename=static_dir(1:lens)//'ltln5c.txt'
          open(7,file=filename)
          read(7,*) lats, lons
          close(7) 

          rpd = 3.1415926535897932 / 180.

          do j=1,nnyp
          do i=1,nnxp
          do is = 1,n_staggers
              lats(i,j,is) = lats(i,j,is) / rpd
              lons(i,j,is) = lons(i,j,is) / rpd
              if(lons(i,j,is) .lt. -180.)then
                  lons(i,j,is) = lons(i,j,is) + 360.
              endif
              if(lons(i,j,is) .gt. +180.)then
                  lons(i,j,is) = lons(i,j,is) - 360.
              endif
          enddo
          enddo
          enddo

        else ! conformal grid ('plrstr','merctr',lambrt')
          call compute_latlon(nnxp,nnyp,n_staggers,mdlat,mdlon
     +                         ,deltax,xtn,ytn,lats,lons,istatus)
          if(istatus.ne.1)then
             print*,'error returned: compute_latlon'
             return
          endif
             
        endif

        if (c6_maproj .ne. 'rotlat' .and. c6_maproj .ne. 'icshdr') then       

        ns = istag

        print*,'ns [staggers array index] = ',ns
        if(ns.ne.1)
     &print*,' using c-stagger grid for lsm fields and other output'

c*****************************************************************
        write(6,*)
        write(6,*)'corner points.'

        call latlon_to_xy(lats(1,1,ns),lons(1,1,ns)      ,erad,xco,yco)
        write(6,701,err=702)1,1,   lats(1,1,ns)   ,lons(1,1,ns),xco,yco 

        call latlon_to_xy(lats(nnxp,1,ns),lons(nnxp,1,ns),erad,xco,yco)
        write(6,701,err=702)nnxp,1,lats(nnxp,1,ns),lons(nnxp,1,ns)
     +                                                         ,xco,yco   

        call latlon_to_xy(lats(1,nnyp,ns),lons(1,nnyp,ns),erad,xco,yco)
        write(6,701,err=702)1,nnyp,lats(1,nnyp,ns),lons(1,nnyp,ns)
     +                                                         ,xco,yco   

        call latlon_to_xy(lats(nnxp,nnyp,ns),lons(nnxp,nnyp,1)
     +                                                   ,erad,xco,yco)
        write(6,701,err=702)nnxp,nnyp,lats(nnxp,nnyp,ns)
     +                                       ,lons(nnxp,nnyp,1),xco,yco   

 701    format(' lat/lon/x/y at ',i5,',',i5,' =',2f12.5,2f12.0)
 702    continue

        call check_domain(lats(1,1,ns),lons(1,1,ns),nnxp,nnyp
     +,grid_spacing_m,1,istat_chk)  
        if(istat_chk.ne.1)then
           print*,'error returned from check_domain'
           istatus = istat_chk
           return
        endif

	endif ! note: totally avoided check_domain for rotlat & icshdr


       if(.true.)then       
           print*,'get perimeter of grid'
           call get_domain_perimeter_grid(nnxp,nnyp,c10_grid_f
     1                  ,lats(1,1,ns),lons(1,1,ns)
     1                  ,0.0,rnorth,south,east,west,istatus)
           print*,'static dir = ',static_dir(1:lens)
           open(10,file=static_dir(1:lens)//'/llbounds.dat'
     +         ,status='unknown')

           print*,'write llbounds.dat'
           write(10,*)rnorth                 
           write(10,*)south          
           write(10,*)east
           write(10,*)west
           close(10)
       endif

       write(6,*)'deltax = ',deltax
       write(6,*)'check_domain:status = ',istat_chk
c
c*****************************************************************
c calculate surface static fields
c
c ------------------------------------------------------
c ----------------- terrain usgs 30 sec ----------------
c type = u
c
       itstatus=ishow_timer()
       print*
       print*,' processing 30s topo data, l_fsf_gridgen = '
     1                                   ,l_fsf_gridgen      

       allocate (topt_30(nnxp,nnyp),
     +           topt_30_s(nnxp,nnyp),
     +           topt_30_ln(nnxp,nnyp),
     +           topt_30_lt(nnxp,nnyp))

       l_topo_wps = l_fsf_gridgen

 600   if(.not. l_topo_wps .and. .not. l_topo_1s)then

        call s_len(path_to_topt30s,len)
        print*,'path to topt30s:        ',path_to_topt30s(1:len)
        path_to_topt30s(len+1:len+2)='/u'

        if (c6_maproj .eq. 'rotlat') then
	 categorical=.false.

	 ncat=1
	 allocate(adum3d(nnxp,nnyp,ncat))
         allocate(adum2d(nnxp,nnyp))

         call alt_10by10_all(nnxp,nnyp,1200,grid_spacing_m/1000.,
     &   hlat,hlon,path_to_topt30s,topt_30,adum3d,ncat,topt_30_lt,
     &   topt_30_ln,topt_30_s,categorical,2,0)


	 write(6,*) 'topo before smooth'
         do j=nnyp,1,-nnyp/30
           write(6,633) (topt_30(i,j),i=1,nnxp,nnxp/15)
         enddo


         deallocate(adum3d)
                                                                                
  633    format(20(f5.0,1x))
                                                                                
!mp	 bogus status value to keep things smooth below
         istatus_30s=1
!mp

        elseif (c6_maproj .eq. 'icshdr') then
         allocate  (geodat3d(nnxp,nnyp,1))
         call proc_geodat(nnxp,nnyp,1,path_to_topt30s
     +       ,lats(1,1,1),lons(1,1,1),data(1,1,1)
     +       ,geodat3d,istatus_30s)

        else

         call geodat(nnxp,nnyp,erad,90.,std_lon,xtn(1,1)
     +   ,ytn(1,1),deltax,deltay,topt_30,topt_30_s,topt_30_ln
     +   ,topt_30_lt,path_to_topt30s,toptwvl,silavwt,new_dem,1
     +   ,istatus_30s)

        endif

        if(istatus_30s .ne. 1)then

          print*,'warning: file(s) missing for 30s terrain data'
          print*,' >>>> try processing 10m topo data <<<<'
          print*

          allocate (topt_10(nnxp,nnyp),
     +              topt_10_s(nnxp,nnyp),
     +              topt_10_ln(nnxp,nnyp),
     +              topt_10_lt(nnxp,nnyp))

          print*
          print*,' calling geodat:  10m terrain data ....'

          call geodat(nnxp,nnyp,erad,90.,std_lon,xtn(1,1)
     +,ytn(1,1),deltax,deltay,topt_10,topt_10_s,topt_10_ln
     +,topt_10_lt,path_to_topt10m,toptwvl,silavwt,new_dem,1
     +,istatus_10m)

          if(istatus_10m .ne. 1)then

             print*
             print*,'error: file(s) missing for both 10m and
     +30s terrain data'
             print*,'error: aborting --- no static file created'
             return

          else

             print *,'topt_10    =',topt_10(1,1),topt_10(nnxp,nnyp)
             print*
             print*,'try blending 10min and 30s terrain data'
             call blend_topo(nnxp,nnyp,lats(1,1,1),lons(1,1,1)
     1,topt_10,topt_10_s,topt_10_ln,topt_10_lt
     1,topt_30,topt_30_s,topt_30_ln,topt_30_lt
     1,topt_out,topt_out_s,topt_out_ln,topt_out_lt)

             deallocate (topt_10,
     1                   topt_10_s,
     1                   topt_10_ln,
     1                   topt_10_lt)

          endif

        else ! go with 30s topo data since all tiles available

          print*,'topo sw and ne corner values'
          print *,'topt_30    =',topt_30(1,1),topt_30(nnxp,nnyp)

          topt_out=topt_30
          topt_out_s=topt_30_s
          topt_out_ln=topt_30_ln
          topt_out_lt=topt_30_lt
          icount_30=nnyp*nnxp
          
        endif !istatus_30s ... data processed ok 

       elseif(l_topo_1s .eqv. .true.)then
          r8lat(:,:) = lats(:,:,1)
	  r8lon(:,:) = lons(:,:,1)
          call get_topo_1s(nnxp,nnyp,grid_spacing_m,r8lat,r8lon,topt_out
     1                    ,rnorth,south,east,west                        ! i
     1                    ,path_to_topt30s,istatus)
          if(istatus .ne. 1)then
              write(6,*)' error in get_topo_1s routine: returning'
              return
          endif

       elseif(.false.)then ! read topo data from wps output netcdf file
          write(6,*)' calling read_wrfstatic for wps topo'
          call s_len(path_to_topt30s,lenp)
          filename_wps = path_to_topt30s(1:lenp)//'/geo_em.d01.nc'
          call read_wrfstatic(nnxp,nnyp,lat,lon,filename_wps,topt_out
     1                       ,istatus)
          if(istatus .ne. 1)then
              write(6,*)' error - no wrf static data: returning'
              return
          endif

       else                ! read topo data from laps/wrf fsf file
          call get_systime_i4(i4time_sys,istatus)
          if(istatus .ne. 1)then
              write(6,*)' error - no systime: returning'
              return
          endif

          istatus = 0
          ntrys = 40
          itry = 0
          do while (itry .le. ntrys .and. istatus .eq. 0)
              i4time_try = i4time_sys - (itry * laps_cycle_time)
              write(6,*)' trying for fsf data at',itry,i4time_try
              call get_modelfg_2d(i4time_try,'ter',nnxp,nnyp,topt_out
     1                                                      ,istatus)       
              itry = itry + 1
          enddo

          if(itry .eq. ntrys)then
              write(6,*)' reached maximum number of tries'
              write(6,*)' use geog topo data instead'
              l_topo_wps = .false.
              goto 600
          endif

       endif

       deallocate (topt_30,
     1             topt_30_s,
     1             topt_30_ln,
     1             topt_30_lt)

       print*
c
c ---------------------------------------------------------
c ------------------- 30s landuse -------------------------
c type = v
c
       allocate  (geodat2d(nnxp,nnyp))
       allocate  (geodat3d(nnxp,nnyp,24))

        if (c6_maproj .ne. 'rotlat') then

       if(.not.allocated(adum2d))allocate (adum2d(nnxp,nnyp))

       itstatus=ishow_timer()
       print*
       print*,' calling geodat: processing 30s landuse data.'
       print*,' re-allocate geodat3d ',nnxp,nnyp,' 24'
       print*

       call geodat(nnxp,nnyp,erad,90.,std_lon,xtn(1,ns)
     +,ytn(1,ns),deltax,deltay,geodat2d,geodat3d
     +,adum2d,adum2d,path_to_luse_30s,2.0,0.0,new_dem,24
     +,istatus)

       if(istatus.ne.1)then
         print*,'********* error ***********'
         print*,'error:  file(s) missing for landuse data'
         print*,'error:  static file not created'
         print*
         return
       endif
	
	else

	categorical=.true.
        call alt_10by10_all(nnxp,nnyp,1200,grid_spacing_m/1000.,
     &		hlat,hlon,path_to_luse_30s,geodat2d,
     &          geodat3d,24,adum2d,adum2d,adum2d,categorical,1,9999)
 
        do j=nnyp,1,-nnyp/25
        write(6,733) (int(geodat2d(i,j)),i=1,nnxp,nnxp/15)
        enddo
        
  733   format(50(i2,1x))
 
	endif

       ilndmsk=12  !land mask is 12 for both wrfsi and laps
       if(c10_grid_f(1:lf).eq.'wrfsi')then

	 if (c6_maproj .ne. 'rotlat') then

c dominant category landuse
          ilnduse=11
          data(:,:,ilnduse)=geodat2d
c landmask for wrfsi
          data(:,:,ilndmsk) = 1.
          where(data(:,:,11) .eq. 16.)data(:,:,ilndmsk)=0.

c grids 15 thru 38 are percent distributions
          i=14
          do j=1,24
             data(:,:,i+j)=geodat3d(:,:,j)
          enddo
	else ! rotlat case
                                                                                
        ilndmsk=7
        ilnduse=6
                                                                                
c landmask for wrfsi
           data(:,:,ilnduse)=geodat2d
           data(:,:,ilndmsk) = 1.
           where(data(:,:,ilnduse) .eq. 16.)data(:,:,ilndmsk)=0.


!!!! require a large % of water to make a water point
!!	
	fraclk=0.7  !! require that 70% of points be water for output to be h2o

	do jj=1,nnyp
	 do ii=1,nnxp
	  if (geodat3d(ii,jj,16) .gt. fraclk) then
	    data(ii,jj,ilndmsk)=0.
          else	
	    if ( data(ii,jj,ilndmsk) .eq. 0.) then
!	      write(6,*) 'water water, becoming land:: ', ii,jj
	    endif
	    data(ii,jj,ilndmsk)=1.
          endif	
         enddo	
	enddo


!!	change dominant landuse  (otherwise will have land points with
!!	water as the dominant class)

	do jj=2,nnyp-1
        do ii=2,nnxp-1 

	if (data(ii,jj,ilndmsk).eq.1.and.data(ii,jj,6).eq.16.) then
	onepos=amax1(data(ii+1,jj,6),data(ii-1,jj,6),
     +   data(ii,jj+1,6),data(ii,jj-1,6))
	twopos=amin1(data(ii+1,jj,6),data(ii-1,jj,6),
     +   data(ii,jj+1,6),data(ii,jj-1,6))

	if (twopos .ne. 16) then
	 	data(ii,jj,6)=twopos
	elseif (onepos .ne. 16) then
		data(ii,jj,6)=onepos
	else
!	leave data alone, revert the land mask to water
	write(6,*) 'switching back to water at ii,jj: ', ii,jj
		data(ii,jj,ilndmsk)=0.
	endif

	endif

	enddo
	enddo

                                                                                
! putting land-use percentages in 19-42 for nmm
        i=18
           do j=1,24
              data(:,:,i+j)=geodat3d(:,:,j)
           enddo
                                                                                
        endif

       else
c dominant category landuse
          ilnduse=5
          data(:,:,ilnduse)=geodat2d   !landuse for laps
c landmask for laps
          data(:,:,ilndmsk)=1.
          where(data(:,:,5) .eq. 16.)data(:,:,ilndmsk)=0.
       endif
c
c **********************************************
c fix water areas to have 0.0m terrain elevation
c currently no fix going on!
c **********************************************
       lforce_ter_2_zero=.false.
       if(lforce_ter_2_zero)then
          where(data(:,:,ilndmsk).eq.0)topt_out=0.0
       endif


c for wrfsi we'll bilinearly interpolate to get the topo
c from the non-staggered topo

       if(c10_grid_f(1:lf).eq.'wrfsi')then

	if (c6_maproj .ne. 'rotlat') then

          data(:,:,9)=topt_out
          call move(topt_out_s,data(1,1,48),nnxp,nnyp)
          call move(topt_out_ln,data(1,1,49),nnxp,nnyp)
          call move(topt_out_lt,data(1,1,50),nnxp,nnyp)

          indxter=51
          data(:,:,indxter)=r_missing_data
          do j=2,nnyp
          do i=2,nnxp
             call bilinear_interp(i,j,nnxp,nnyp,topt_out,result)
             if(data(i,j,51).gt.30000. .and.
     & data(i,j,51).lt.r_missing_data)then
                print*,'here: large topo: i/j/topo ',i,j,data(i,j,51)
             endif
             data(i-1,j-1,indxter)=result
          enddo
          enddo

          where((data(:,:,indxter).lt. 0.01).and.
     1          (data(:,:,indxter).gt.-0.01))data(:,:,indxter)=0.0

	else

!       for rotlat, topt_out is already on mass point stagger.
!       no bilinear interp needed


!!	boundary smooth the topo data

	write(6,*) 'call smdhld'
	write(6,*) 'nnxp, nnyp: ', nnxp,nnyp
	call smdhld(nnxp,nnyp,topt_out,1.-data(:,:,ilndmsk),12,12)

c       write(6,*) 'topo after smoothing'
c       do j=nnyp,1,-nnyp/30
c       write(6,633) (topt_out(i,j),i=1,nnxp,nnxp/15)
c       enddo


!-------------4-point averaging of mountains along inner boundary-------

          do i=1,nnxp-1
      topt_out(i,2)=0.25*(topt_out(i,1)+topt_out(i+1,1)+ 
     +                    topt_out(i,3)+topt_out(i+1,3))
          enddo

          do i=1,nnxp-1
      topt_out(i,nnyp-1)=0.25*(topt_out(i,nnyp-2)+topt_out(i+1,nnyp-2)+ 
     +                         topt_out(i,nnyp)+topt_out(i+1,nnyp))
          enddo

          do j=4,nnyp-3,2
      topt_out(1,j)=0.25*(topt_out(1,j-1)+topt_out(2,j-1)+
     +                    topt_out(1,j+1)+topt_out(2,j+1))
          enddo

          do j=4,nnyp-3,2
      topt_out(nnxp,j)=0.25*(topt_out(nnxp-1,j-1)+topt_out(nnxp,j-1)+
     +                         topt_out(nnxp-1,j+1)+topt_out(nnxp,j+1))
          enddo


	indxter=18
        do j=1,nnyp
        do i=1,nnxp
          data(i,j,indxter)=topt_out(i,j)
        enddo
        enddo

          where((data(:,:,indxter).lt. 0.01).and.
     1          (data(:,:,indxter).gt.-0.01))data(:,:,indxter)=0.0

	endif


       else  !must be laps

          indxter=3
          call move(topt_out,data(1,1,indxter),nnxp,nnyp)         ! kwd
          call move(topt_out_s,data(1,1,7),nnxp,nnyp)             ! js
          call move(topt_out_ln,data(1,1,8),nnxp,nnyp)            ! js
          call move(topt_out_lt,data(1,1,9),nnxp,nnyp)            ! js

       endif 

       print*
       print *,'terrain  = ',data(1,1,indxter),data(nnxp,nnyp,indxter)

c      print *,'# of grid pts using 30 sec terrain data =  ',icount_30
c      print *,'# of grid pts using blended terrain data = ',icount_ramp
       print*
c
c output section for topo, lat-lons and corners
c -----------------------------------------------
       if(c10_grid_f(1:lf).eq.'wrfsi')then
         if(c6_maproj.ne.'rotlat')then
          write(cnest,'(i2.2)')nest
          clatlon='latlon.d'//cnest//'.dat'            !binary non-staggered lat-lon file
          clatlons='latlon-mass.d'//cnest//'.dat'      !binary staggered lat-lon file
          ctopo='topo.d'//cnest//'.dat'                !binary non-staggered topo file
          ctopos='topo-mass.d'//cnest//'.dat'          !binary staggered topo file
          ctopog='topography.d'//cnest//'.dat'         !ascii non-staggered topo file
          ctopogs='topography-mass.d'//cnest//'.dat'   !ascii staggered topo file
          clatlon2d='latlon2d.d'//cnest//'.dat'        !ascii lat-lon non-staggered file
          clatlons2d='latlon2d-mass.d'//cnest//'.dat'  !ascii lat-lon staggered file
          corners='corners.d'//cnest//'.dat'           !non-staggered corners
          cornerss='corners-mass.d'//cnest//'.dat'     !staggered corners
         else
          clatlon='latlon-rotlat.dat'        !binary lat/lon file
          clatlon2d='latlon2d-rotlat.dat'    !ascii lat/lon file
          ctopo='topo-rotlat.dat'            !binary topo file
          ctopog='topography-rotlat.dat'     !ascii topo file
          corners='corners-rotlat.dat'       !domain corners file
         endif
       else
          clatlon='latlon.dat'        !binary lat/lon file
          clatlon2d='latlon2d.dat'    !ascii lat/lon file
          ctopo='topo.dat'            !binary topo file
          ctopog='topography.dat'     !ascii topo file
          corners='corners.dat'       !domain corners file
       endif

       call s_len(static_dir,len)
       if(c10_grid_f(1:lf).eq.'wrfsi')then

         if(c6_maproj.ne.'rotlat')then
           in=0
c write both non-stagger and stagger data for wrfsi, non-rotlat
           do k=1,2
            if(k.eq.1)then
             is=1
             idxt=9
             filename=static_dir(1:len)//clatlon
             open(10,file=filename,status='unknown',form='unformatted')
             filename=static_dir(1:len)//ctopo
             open(11,file=filename,status='unknown',form='unformatted')
             filename=static_dir(1:len)//corners
             open(15,file=filename,status='unknown')
             filename=static_dir(1:len)//ctopog
             open(666,file=filename)
             filename=static_dir(1:len)//clatlon2d
             open(667,file=filename)
            else
             is = 4
             in = 1
             idxt=51
             filename=static_dir(1:len)//clatlons
             open(10,file=filename,status='unknown',form='unformatted')
             filename=static_dir(1:len)//ctopos
             open(11,file=filename,status='unknown',form='unformatted')
             filename=static_dir(1:len)//cornerss
             open(15,file=filename,status='unknown')
             filename=static_dir(1:len)//ctopogs
             open(666,file=filename)
             filename=static_dir(1:len)//clatlons2d
             open(667,file=filename)
            endif

            write(10)lats(1:nnxp-in,1:nnyp-in,is)
     +,lons(1:nnxp-in,1:nnyp-in,is)
            close(10)
            write(11)data(1:nnxp-in,1:nnyp-in,idxt)
            close(11)
            write(15,*)lats(1,1,is),lons(1,1,is)
            write(15,*)lats(1,nnyp-in,is),lons(1,nnyp-in,is)
            write(15,*)lats(nnxp-in,1,is),lons(nnxp-in,1,is)
            write(15,*)lats(nnxp-in,nnyp-in,is),lons(nnxp-in,nnyp-in,is)
            close(15)

            do j=1,nnyp-in
            do i=1,nnxp-in
               write(666,*) data(i,j,idxt)
               write(667,*) lats(i,j,is),lons(i,j,is)
            enddo
            write(666,'()')
            write(667,'()')
            enddo
            close(666)
            close(667)

           enddo

         else   !if(c6_maproj.eq.'rotlat')then

           is=1
           idxt=18
           filename=static_dir(1:len)//clatlon
           open(10,file=filename,status='unknown',form='unformatted')
           filename=static_dir(1:len)//ctopo
           open(11,file=filename,status='unknown',form='unformatted')
           filename=static_dir(1:len)//corners
           open(15,file=filename,status='unknown')
           filename=static_dir(1:len)//ctopog
           open(666,file=filename)
           filename=static_dir(1:len)//clatlon2d
           open(667,file=filename)

           write(10)lats(1:nnxp,1:nnyp,is),lons(1:nnxp,1:nnyp,is)
           close(10)
           write(11)data(:,:,idxt)
           close(11)
           write(15,*)lats(1,1,is),lons(1,1,is)
           write(15,*)lats(1,nnyp,is),lons(1,nnyp,is)
           write(15,*)lats(nnxp,1,is),lons(nnxp,1,is)
           write(15,*)lats(nnxp,nnyp,is),lons(nnxp,nnyp,is)
           close(15)

           do j=1,nnyp
           do i=1,nnxp
              write(666,*) data(i,j,idxt)
              write(667,*) lats(i,j,is),lons(i,j,is)
           enddo
           write(666,'()')
           write(667,'()')
           enddo
           close(666)
           close(667)

         endif

       else !laps case

         is=1
         idxt=3
!        filename=static_dir(1:len)//clatlon
!        open(10,file=filename,status='unknown',form='unformatted')
!        filename=static_dir(1:len)//ctopo
!        open(11,file=filename,status='unknown',form='unformatted')
         filename=static_dir(1:len)//corners
         open(15,file=filename,status='unknown')
         filename=static_dir(1:len)//ctopog
         open(666,file=filename)
         filename=static_dir(1:len)//clatlon2d
         open(667,file=filename)

!        write(10)lats(1:nnxp,1:nnyp,is),lons(1:nnxp,1:nnyp,is)
!        close(10)
!        write(11)data(:,:,idxt)
!        close(11)
         write(15,*)lats(1,1,is),lons(1,1,is)
         write(15,*)lats(1,nnyp,is),lons(1,nnyp,is)
         write(15,*)lats(nnxp,1,is),lons(nnxp,1,is)
         write(15,*)lats(nnxp,nnyp,is),lons(nnxp,nnyp,is)
         close(15)

         do j=1,nnyp
         do i=1,nnxp
            write(666,*) data(i,j,idxt)
            write(667,*) lats(i,j,is),lons(i,j,is)
         enddo
         write(666,'()')
         write(667,'()')
         enddo
         close(666)
         close(667)

       endif
c ----------------------------------------------------------
c
c ------ land fraction from 30s land use -------------------
c type = l (no longer used ... no land fraction data base)
c
        print*
        print*,' create from 30s land use fractional dist'
        print*
 
        geodat2d(:,:)=1.- geodat3d(:,:,16)
c
c this function allows variable amounts of smoothing of
c the land fraction data depending on grid spacing
c
!       iter=min(10,int(8000./deltax))
        iter = 0

	if (c6_maproj .ne. 'rotlat') then
        
          print*,'filter land fraction with 2dx ',iter
          do k=1,iter
             call filter_2dx(geodat2d,nnxp,nnyp,1, 0.5)
             call filter_2dx(geodat2d,nnxp,nnyp,1,-0.5)
             geodat2d = max(geodat2d,0.) 
          enddo

	endif

        if(c10_grid_f(1:lf).eq.'wrfsi')then
	
	    if (c6_maproj .ne. 'rotlat') then
		idx=10
	    else
		idx=5
	    endif

	else ! not wrfsi
		idx=4
	endif

        data(:,:,idx)=geodat2d
        print *,'pctlandfrac=',geodat2d(1,1),geodat2d(nnxp,nnyp)       
c
c -------------------------------------------------------------------
c potential fix of fictitous islands for certain resolution domains.
c story here is that west coast terrain can be "funky" possibly due
c to steepness and method to compute avg terrain.

	write(6,*) 'idx, indxter: ', idx, indxter
        where(data(:,:,idx).le.0.1 .and. data(:,:,indxter).lt.5.0)
     &data(:,:,indxter)=0.0

c
c ------------------------------------------------------------
c -------------- top layer soil type -------------------------
c type = o
c
        itstatus=ishow_timer()
        print*
        print*,' processing 30s soil type top layer data....'

        deallocate(geodat3d)
        allocate (geodat3d(nnxp,nnyp,16))

        if (c6_maproj .ne. 'rotlat') then

        call geodat(nnxp,nnyp,erad,90.,std_lon,xtn(1,ns),ytn(1,ns)
     +,deltax,deltay,geodat2d,geodat3d,adum2d,adum2d
     +,path_to_soiltype_top_30s,2.0,0.0,new_dem,16   !maxdatacat
     +,istatus_soil)

        if(istatus_soil.ne.1)then
           print*
           print*,' soil type data not processed completely'
           if(c10_grid_f(1:lf).eq.'wrfsi')then
             print*,' file(s) missing for soil type top layer data'
             print*,' error:  static file not created'
             print*
             return
           else
             print*,' file(s) missing for soil type top layer data'
             print*,'           *** warning ***'
             print*,' soil type top data not added to static file'
             print*
           endif
           geodat2d=r_missing_data
           geodat3d=r_missing_data
        endif

	else

	categorical=.true.
        call alt_10by10_all(nnxp,nnyp,1200,grid_spacing_m/1000.,
     &		hlat,hlon,path_to_soiltype_top_30s,geodat2d,
     &          geodat3d,16,adum2d,adum2d,adum2d,categorical,1,9999)
                                                                                
        write(6,*) 'soiltype_top'
        do j=nnyp,1,-nnyp/25
        write(6,733) (int(geodat2d(i,j)),i=1,nnxp,nnxp/15)
        enddo
                                                                                
        endif

        if(c10_grid_f(1:lf).eq.'wrfsi')then

         if (c6_maproj .ne. 'rotlat') then
c
c make water points = 0 for adjust_geog. later we'll put it back to original
c
           where(geodat2d.eq.14)geodat2d=0.0

           call adjust_geog(nnxp-1,nnyp-1,1,'soiltype',istatus_soil
     &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
     &,data(1:nnxp-1,1:nnyp-1,ilndmsk),geodat2d(1:nnxp-1,1:nnyp-1)
     &,istatus)

           where(geodat2d.eq.0.0)geodat2d=14.
           idxstl=13
           data(:,:,idxstl)=geodat2d
           i=51
           do j=1,16
              data(:,:,i+j)=geodat3d(:,:,j)
           enddo
	 else
           idxstl=8
           data(:,:,8)=geodat2d
           i=42
           do j=1,16
              data(:,:,i+j)=geodat3d(:,:,j)
           enddo
         endif

        else
c
c make water points = 0 for adjust_geog. after put it back to original
c
           where(geodat2d.eq.14)geodat2d=0.0

           call adjust_geog(nnxp,nnyp,1,'soiltype',istatus_soil
     &,lats(1,1,ns),data(1,1,indxter),data(1,1,ilndmsk),geodat2d
     &,istatus)

           where(geodat2d.eq.0.0)geodat2d=14.
           idxstl=10
           data(:,:,idxstl)=geodat2d

        endif
c
c --------------------------------------------------------------
c -------------- bottom layer soil type ------------------------
c type = o
c
        itstatus=ishow_timer()
        print*
        print*,' processing 30s soil type bottom layer data'

        if (c6_maproj .ne. 'rotlat') then

        call geodat(nnxp,nnyp,erad,90.,std_lon,xtn(1,ns)
     +,ytn(1,ns),deltax,deltay,geodat2d,geodat3d,adum2d,adum2d
     +,path_to_soiltype_bot_30s,2.0,0.0,new_dem,16   !maxdatacat
     +,istatus_soil)

        if(istatus_soil.ne.1)then
           print*
           print*,' bottom layer soil data not processed completely'
           if(c10_grid_f(1:lf).eq.'wrfsi')then
             print*,' file(s) missing for soil type bot layer data'
             print*,' error:  static file not created'
             print*
             return
           else
             print*,' file(s) missing for soil type bot layer data'
             print*,'           *** warning ***'
             print*,' soil type bot data not added to static file'
             print*
           endif
           geodat2d=r_missing_data
           geodat3d=r_missing_data
        endif

	else

	categorical=.true.
        call alt_10by10_all(nnxp,nnyp,1200,grid_spacing_m/1000.,
     &		hlat,hlon,path_to_soiltype_bot_30s,geodat2d,
     &          geodat3d,16,adum2d,adum2d,adum2d,categorical,1,9999)

        endif
c
c soiltype bottom layer ... 68 thru 83
c
        if(c10_grid_f(1:lf).eq.'wrfsi')then

          if (c6_maproj .ne. 'rotlat') then
c
c make water points = 0 for adjust_geog. after, put it back to original.
c
            where(geodat2d.eq.14)geodat2d=0.0

            call adjust_geog(nnxp-1,nnyp-1,1,'soiltype',istatus_soil
     &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
     &,data(1:nnxp-1,1:nnyp-1,ilndmsk),geodat2d(1:nnxp-1,1:nnyp-1)
     &,istatus)

            where(geodat2d.eq.0.)geodat2d=14.0
            idxsbl=14
            data(:,:,idxsbl)=geodat2d
            i=67
            do j=1,16
               data(:,:,i+j)=geodat3d(:,:,j)
            enddo
	  else
            idxsbl=9
	    data(:,:,9)=geodat2d
            i=58
            do j=1,16
               data(:,:,i+j)=geodat3d(:,:,j)
            enddo
          endif

        else
c
c make water points = 0 for adjust_geog. after, put it back to original.
c
            where(geodat2d.eq.14)geodat2d=0.0

            call adjust_geog(nnxp,nnyp,1,'soiltype',istatus_soil
     &,lats(1,1,ns),data(1,1,indxter),data(1,1,ilndmsk),geodat2d
     &,istatus)
            where(geodat2d.eq.0.)geodat2d=14.0
            idxsbl=11
            data(:,:,idxsbl)=geodat2d
        endif
c
c ---------------------------------------------------------
c ----------------- greenness fraction --------------------
c type = g
c
        itstatus=ishow_timer()
        print*
        print*,' calling proc_geodat: process green frac data.'
        print*,' re-allocate geodat3d ',nnxp,nnyp,' 12'
        print*

        deallocate(geodat3d)
        allocate  (geodat3d(nnxp,nnyp,12),stat=istat)
        if(istat.ne.0)then
           print*,'error allocating data object geodat3d ',istat
           if(c10_grid_f(1:lf).eq.'wrfsi')then
              print*,'error: aborting process'
              print*,'error: maybe not enough memory'
              print*,'error: cannot allocate geodat3d'
              stop
           else
              print*,'warning: continue without green fraction'         
           endif
        endif

	if (c6_maproj .ne. 'rotlat') then

        print*,'calling proc_geodat'
        call proc_geodat(nnxp,nnyp,12,path_to_green_frac
     +,lats(1,1,ns),lons(1,1,ns),data(1,1,ilndmsk)
     +,geodat3d,istatus_grn)

	else
                                                                                
	categorical=.false.
	useland=.true.
        call alt_hemi_all(nnxp,nnyp,hlat,hlon,
     &                       path_to_green_frac,
     &                       type,12,grid_spacing_m/1000.,geodat3d,
     &                       geodat2d,categorical,data(1,1,ilndmsk),
     &			     useland)


!	convert from 0-100% to 0-1 fraction

	do n=1,12
	do j=1,nnyp
	do i=1,nnxp
	    geodat3d(i,j,n)=geodat3d(i,j,n)/100.
	enddo
	enddo
	enddo
                                                                                
        istatus_grn=1 !bogus value
                                                                                
        endif

        if(istatus_grn.ne.1)then
         print*
         print*,'greenness fraction data not processed completely'
         if(c10_grid_f(1:lf).eq.'wrfsi')then
            print*,' error: file(s) missing for green frac data'
            print*,' error: static file not created'
            print*
            istatus=0
            return
         else
            print*,' warning: file(s) missing for green frac data'
            print*,' warning: green frac not added to static file'
            print*
         endif
         geodat3d=r_missing_data
        endif

        if(c10_grid_f(1:lf).eq.'wrfsi')then

           if (c6_maproj .ne. 'rotlat') then

              igrn=83

c             call adjust_geog(nnxp-1,nnyp-1,12,'greenfrac',istatus_grn
c    &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
c    &,data(1:nnxp-1,1:nnyp-1,ilndmsk),geodat3d(1:nnxp-1,1:nnyp-1,:)
c    &,istatus)

           else
               igrn=74
           endif

        else

           igrn=12
          
           call adjust_geog(nnxp,nnyp,igrn,'greenfrac',istatus_grn
     &,lats(:,:,ns),data(:,:,indxter),data(:,:,ilndmsk),geodat3d
     &,istatus)
        endif
        do j=1,12
           data(:,:,igrn+j)=geodat3d(:,:,j)
        enddo
c
c annual max/min greenness fraction in domain
c --------------------------------------------
        if(c10_grid_f(1:lf).eq.'wrfsi')then
           dommaxgf=0.0
           dommingf=99.0
           print*,' compute max/min greenness frac at grid points'
           if (c6_maproj .ne. 'rotlat') then
               mxgiwrf=110
               mngiwrf=111
           else
               mxgiwrf=101
               mngiwrf=102
           endif
           do j=1,nnyp
           do i=1,nnxp
              data(i,j,mxgiwrf)=maxval(geodat3d(i,j,:))
              data(i,j,mngiwrf)=minval(geodat3d(i,j,:))
              if(data(i,j,mxgiwrf).gt.dommaxgf)
     &dommaxgf=data(i,j,mxgiwrf)
              if(data(i,j,mngiwrf).lt.dommingf)
     &dommingf=data(i,j,mngiwrf)
           enddo
           enddo
           write(6,*)'domain annual max/min green fraction',
     &dommaxgf,dommingf
        endif
c
c ---------------------------------------------------------
c --------------- deep soil temperature -------------------
c type = t
c
	if (c6_maproj .ne. 'rotlat') then
        itstatus=ishow_timer()
        print*
        print*,' call proc_geodat: process 1 deg soiltemp data.'

        geodat2d=0.0

        call proc_geodat(nnxp,nnyp,1
     1,path_to_soiltemp_1deg,lats(1,1,ns),lons(1,1,ns)
     1,data(1,1,ilndmsk),geodat2d,istatus_tmp)

	else

	write(6,*) 'call alt_10by10 for soiltemp'
        geodat2d=0.0
	categorical=.false.

	ncat=1
        call alt_10by10_all(nnxp,nnyp,10,grid_spacing_m/1000.,
     &		hlat,hlon,path_to_soiltemp_1deg,geodat2d,
     &          adum3d,ncat,adum2d,adum2d,adum2d,categorical,2,0)

!       convert to degrees
                                                                                
	write(6,*) 'convert to degrees ', nnxp,nnyp
        do j=1,nnyp
        do i=1,nnxp
                geodat2d(i,j)=geodat2d(i,j)/100.
        enddo
        enddo
                                                                                
        istatus_tmp=1 !bogus
                                                                                
        endif
	
!       do j=nnyp,1,-nnyp/25
!         write(6,833) (geodat2d(i,j),i=1,nnxp,nnxp/15)
!       enddo

        if(istatus_tmp.ne.1)then
         print* 
         print*,'soiltemp data not processed completely' 
         if(c10_grid_f(1:lf).eq.'wrfsi')then 
            print*,' error:  file(s) missing for soiltemp data' 
            print*,' error:  static file not created'
            print*
            istatus=0
            return
         else
            print*,' file(s) missing for soiltemp data'
            print*,'           *** warning ***'
            print*,' soiltemp data not added to static file'
            print*
         endif
         geodat2d=r_missing_data
        endif

        idst=25
        if(c10_grid_f(1:lf).eq.'wrfsi')then 

            if (c6_maproj .ne. 'rotlat') then

                idst=96

                call adjust_geog(nnxp-1,nnyp-1,1,'soiltemp',istatus_tmp
     &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
     &,data(1:nnxp-1,1:nnyp-1,ilndmsk),geodat2d(1:nnxp-1,1:nnyp-1)
     &,istatus)

            else
                idst=87
            endif


        else

            call adjust_geog(nnxp,nnyp,1,'soiltemp',istatus_grn
     &,lats(:,:,ns),data(:,:,indxter),data(:,:,ilndmsk),geodat2d
     &,istatus)

        endif 
	write(6,*) 'idst= ', idst
        data(:,:,idst)=geodat2d
	write(6,*) 'idst written'
c
c -------------- terrain slope index categories -------------------------
c type = i
c
        itstatus=ishow_timer()
        print*
        print*,' processing 1 deg islope data'

        allocate (adum3d(nnxp,nnyp,9))  !temporarily holds %dist of islope cat

	if (c6_maproj .ne. 'rotlat') then

        call geodat(nnxp,nnyp,erad,90.,std_lon,xtn(1,ns)
     +,ytn(1,ns),deltax,deltay,geodat2d,adum3d,adum2d,adum2d
     +,path_to_islope,2.0,0.0,new_dem,9,istatus_slp)

	else

	categorical=.false.
!!	should be done as n.n. interp by alt_interp code
	useland=.false.
        call alt_hemi_all(nnxp,nnyp,hlat,hlon,
     &                       path_to_islope,
     &                       type,1,grid_spacing_m/1000.,adum3d,
     &                       geodat2d,categorical,data(:,:,ilndmsk),
     &			     useland)

	istatus_slp=1 ! bogus

	endif

        deallocate (adum3d)

        if(istatus_slp.ne.1)then
           print*
           print*,' ter slope category data not processed completely'
           if(c10_grid_f(1:lf).eq.'wrfsi')then
             print*,' error: file(s) missing for islope data'
             print*,' error: static file not created'
             print*
             return
           else
             print*,' warning: file(s) missing for islope data'
             print*,' warning: islope not added to static file'
             print*
           endif
           geodat2d=r_missing_data
           geodat3d=r_missing_data
        endif

c put the categories back to the original raw data. if it is a land
c point but islope indicates water, force islope = 1.
        where(geodat2d .eq. 8)geodat2d=13.0
        where(geodat2d .eq. 9)geodat2d=0.0

        if(c10_grid_f(1:lf).eq.'wrfsi')then
	
	   if (c6_maproj .ne. 'rotlat') then

               islp=109
               call adjust_geog(nnxp-1,nnyp-1,1,'islope',istatus_slp
     &,lats(1:nnxp-1,1:nnyp-1,ns),data(1:nnxp-1,1:nnyp-1,indxter)
     &,data(1:nnxp-1,1:nnyp-1,ilndmsk),geodat2d(1:nnxp-1,1:nnyp-1)
     &,istatus)

           else
               islp=100 ! will need more adjustments down the road for this
           endif

        else

           islp=38
           call adjust_geog(nnxp,nnyp,1,'islope',istatus_slp
     &,lats(1,1,ns),data(1,1,indxter),data(1,1,ilndmsk),geodat2d
     &,istatus)

        endif

        if(istatus.ne.1)then
           print*,' warning: processing incomplete: adjust_geog_data'
           if(c10_grid_f(1:lf).eq.'wrfsi')then
              print*,'error: wrfsi static file not generated'
              return
           endif
        endif

        data(:,:,islp)=geodat2d

c force land points to have the correct (default) islope value.
        where(data(:,:,islp)   .eq. 0.0 .and.
     .        data(:,:,ilndmsk).eq. 1.0)data(:,:,islp)=1.0
c force water points to have the correct category for islope
        where(data(:,:,ilndmsk).eq. 0.0)data(:,:,islp)=0.0

!write(6,*) 'islope: '
!      do j=nnyp,1,-nnyp/30
!       write(6,733) (int(data(i,j,islp)),i=1,nnxp,nnxp/20)
!       enddo
c
c    subroutine adjust_geog_data now named adjust_geog and called
c    after each data type is processed (as opposed to all geog data
c    getting adjusted with one call).
c
c!js: even though matt pyle suggests that adjust geog is fine for nmm
c     it currently isn't used.  if it should be used, then it would be
c     similar to the laps case which is not a staggered grid.
c
!mp	believe adjust_geog_data should be fine for nmm case
c
c -------------------------------------------------------------------
c ---------------monthly albedo climatology--------------------------
c type = a
c
        itstatus=ishow_timer()

        if(c10_grid_f(1:lf).eq.'wrfsi')then
          if (c6_maproj .ne. 'rotlat') then
            ialb=96
          else
            ialb=87
          endif
        else
          ialb=25
        endif

        print*
        print*,' processing albedo climo data into ialb/bm = '
     +        ,ialb,l_albedo_bm

        if(.not. l_albedo_bm)then ! obtain albedo from blue marble (bmng) data

	  if (c6_maproj .ne. 'rotlat') then

            call proc_geodat(nnxp,nnyp,12,path_to_albedo
     +         ,lats(1,1,ns),lons(1,1,ns),data(1,1,ilndmsk)
     +         ,geodat3d,istatus_alb)                            ! o

	  else

	    categorical=.false.
	    useland=.true.
            call alt_hemi_all(nnxp,nnyp,hlat,hlon,
     &                       path_to_albedo,
     &                       type,12,grid_spacing_m/1000.,geodat3d,
     &                       adum2d,categorical,data(1,1,ilndmsk),
     &			     useland)

            istatus_alb=1 !bogus

	    do n=1,12
	     do j=1,nnyp
	      do i=1,nnxp
	       geodat3d(i,j,n)=geodat3d(i,j,n)*.01
	      enddo
	     enddo
	    enddo

          endif

          print*,'done in proc_geodat: albedo'

          if(istatus_alb.ne.1)then
           print*
           print*,' error: albedo climo data not processed completely'
           print*,' error: file(s) missing for albedo data'
           print*,' error: static file not created: albedo missing'
           print*
           istatus=0
           return
          endif

c force water points to 0.08 albedo

          do k=1,12
           where(data(:,:,ilndmsk).eq.0.0)geodat3d(:,:,k)=0.08
          enddo

          do j=1,12
              data(:,:,ialb+j)=geodat3d(:,:,j)
          enddo

	else ! use bmng

	  call i4time_fname_lp('170010000',i4_yrstart,istatus)

          do j=1,12
            i4_midmonth = i4_yrstart + ((j*30)-15) * 86400 ! approximate
            call land_albedo_bm(lats(1,1,ns),lons(1,1,ns),nnxp,nnyp
     +                         ,i4_midmonth,data(:,:,ialb+j),istat_bm)
          enddo
	endif
c
c ---------------- max snow albedo ------------------
c type = m
c
	if (c6_maproj .ne. 'rotlat') then

        call proc_geodat(nnxp,nnyp,1,path_to_maxsnoalb
     +,lats(1,1,ns),lons(1,1,ns),data(1,1,ilndmsk)
     +,geodat2d,istatus_alb)

	else

	categorical=.false.

	useland=.true.
	allocate(adum3d(nnxp,nnyp,1))
        call alt_hemi_all(nnxp,nnyp,hlat,hlon,
     &                       path_to_maxsnoalb,
     &                       type,1,grid_spacing_m/1000.,adum3d,
     &                       geodat2d,categorical,data(1,1,ilndmsk),
     &			     useland)

	do j=1,nnyp
	do i=1,nnxp
	geodat2d(i,j)=geodat2d(i,j)*.01
	enddo
	enddo

	endif

        print*,'done in proc_geodat: max snow albedo'

        if(istatus_alb.ne.1)then
         print*
         if(c10_grid_f(1:lf).eq.'wrfsi')then
            print*,'--------------- wrfsi ------------------'
            print*,'error: max snow albedo not processed completely'
            print*,'error: static file not created '
            print*
            istatus=0
            return
         else
            print*,'warning: max snow albedo not processed completely'
         endif
        endif
c
c force max albedo = 0.08 over water. force max albedo = 0.7 over ice
c
        where(data(:,:,ilndmsk).eq.0.0)geodat2d=0.08
        where(data(:,:,ilnduse).eq.24.0)geodat2d=0.7

        if(c10_grid_f(1:lf).eq.'wrfsi')then

	if (c6_maproj .ne. 'rotlat') then
           mxsnalb=47
           data(:,:,mxsnalb)=geodat2d
	else
           mxsnalb=14
	   data(:,:,mxsnalb)=geodat2d	
	endif

        else
           mxsnalb=6
           data(:,:,mxsnalb)=geodat2d
        endif

        deallocate (geodat2d)
        deallocate (geodat3d)
c
c ---------------------------------------------------------------------------------
c let's compare the grids to landmask to assess their consistency (or lack thereof)
c ---------------------------------------------------------------------------------

       if(c10_grid_f(1:lf).eq.'wrfsi')then
        print*,'total number of gridpoints (nx*ny) = ',nnxp*nnyp
        print*
        do i=1,10
           if(i == 1)ii=indxter  !terrain
           if(i == 2)ii=idxstl   !soil top layer
           if(i == 3)ii=idxsbl   !soil bot layer
           if(i == 4)ii=igrn+6   !mxgiwrf  !max greenness
           if(i == 5)ii=mngiwrf  !min greenness
           if(i == 6)ii=idst     !deep soil temp
           if(i == 7)ii=islp     !terrain slope index
           if(i == 8)ii=ialb+6   !albedo; arbitrarily at month 6
           if(i == 9)ii=mxsnalb  !max snow albedo
           if(i == 10)ii=ilnduse !landuse dominant category
           call gridcompare(nnxp,nnyp,i,data(1,1,ii),data(1,1,ilndmsk)
     &,istatus)
        enddo
       endif
c
c ---------------------------------------------------------------------------------
c --------- this is were we prepare for writing the netcdf static file ------------
c ---------------------------------------------------------------------------------

        if(.not.allocated(var))allocate (var(ngrids)
     +,comment(ngrids))

        if(c10_grid_f(1:lf).eq.'wrfsi')then

	if (c6_maproj .ne. 'rotlat') then

           in1=0
           in2=0
           do j=1,ns
              in1=in2+1
              in2=in1+1
              data(:,:,in1)=lats(:,:,j)
              data(:,:,in2)=lons(:,:,j)
           enddo


c compute map/grid information. reallocate geodat3d array
           allocate (geodat3d(nnxp,nnyp,4))
           call get_projrot_grid(nnxp,nnyp,lats(1,1,ns)
     +,lons(1,1,ns),geodat3d,istatus)

           i=39

           data(:,:,i)  =geodat3d(:,:,1)
           data(:,:,i+1)=geodat3d(:,:,2)
 
           call get_map_factor_grid(nnxp,nnyp,n_staggers
     +,lats,lons ,geodat3d,istatus)
           if(istatus.ne.1)then
              print*,'error returned: get_maps_factor_grid'
              return
           endif

           do j=1,ns
              data(:,:,i+j+1)=geodat3d(:,:,j)
           enddo
c           
           call get_coriolis_components(nnxp,nnyp,lats(1,1,ns)
     +,geodat3d)
           data(:,:,i+6)=geodat3d(:,:,1)
           data(:,:,i+7)=geodat3d(:,:,2)

           deallocate(geodat3d)

c          call move(static_albedo,data(1,1,i+8),nnxp,nnyp)

           call move(topt_out_s,data(1,1,i+9),nnxp,nnyp)
           call move(topt_out_ln,data(1,1,i+10),nnxp,nnyp)
           call move(topt_out_lt,data(1,1,i+11),nnxp,nnyp)
c          call move(topt_stag_out,data(1,1,i+12),nnxp,nnyp) !51

           call get_gridgen_var(ngrids,ngrids,var,comment)

	else ! rotlat case
                                                                                

!!!!!	these doloops possibly redundant.  defined above????

!        do j=1,nnyp
!        do i=1,nnxp
!              data(i,j,1)=hlat(i,j)
!              data(i,j,2)=hlon(i,j)
!              data(i,j,3)=vlat(i,j)
!              data(i,j,4)=vlon(i,j)
!        enddo
!        enddo
                                                                                
           data(:,:,18)=topt_out
                                                                                
c compute map/grid information. reallocate geodat3d array
           allocate (geodat3d(nnxp,nnyp,4))
                                                                                
                                                                                
           call vecrot_rotlat(nnxp,nnyp,mdlat,mdlon
     +,vlat,vlon,geodat3d(:,:,1),geodat3d(:,:,2))

!        write(6,*) 'cos term'
        do j=nnyp,1,-nnyp/40
!        write(6,833) (geodat3d(i,j,1),i=1,nnxp,nnxp/15)
        enddo
!        write(6,*) 'sin term'
        do j=nnyp,1,-nnyp/40
!        write(6,833) (geodat3d(i,j,2),i=1,nnxp,nnxp/15)
        enddo

  833	format(30(f5.1,1x))

           i=10
                                                                                
!       sin, then cosine?
           data(:,:,i)  =geodat3d(:,:,2)
           data(:,:,i+1)=geodat3d(:,:,1)
c
                                                                                
!!!!    coriolis just dependent on lat (hopefully)
           call get_coriolis_components(nnxp,nnyp,vlat
     +,geodat3d)
           data(:,:,i+2)=geodat3d(:,:,1)
           data(:,:,i+3)=geodat3d(:,:,2)

           deallocate(geodat3d)

c          call move(static_albedo,data(1,1,i+4),nnxp,nnyp)
                                                                                
           call move(topt_out_s,data(1,1,i+5),nnxp,nnyp)
           call move(topt_out_ln,data(1,1,i+6),nnxp,nnyp)
           call move(topt_out_lt,data(1,1,i+7),nnxp,nnyp)
c          call move(topt_stag_out,data(1,1,i+8),nnxp,nnyp) !18
                                                                                
!030513  made modifications within the get_gridgen_var code
!	to handle nmm variable set
                                                                                
        write(6,*) 'calling get_gridgen_var with ngrids= ', ngrids
           call get_gridgen_var(ngrids,ngrids,var,comment)

        endif

        else ! non-wrfsi

           call move(lats(1,1,1),data(1,1,1),nnxp,nnyp)            ! kwd
           call move(lons(1,1,1),data(1,1,2),nnxp,nnyp)            ! kwd
           call move(topt_out,data(1,1,3),nnxp,nnyp)               ! kwd
           call move(topt_out_s,data(1,1,7),nnxp,nnyp)             ! js
           call move(topt_out_ln,data(1,1,8),nnxp,nnyp)            ! js
           call move(topt_out_lt,data(1,1,9),nnxp,nnyp)            ! js 

           call get_gridgen_var(ngrids,ngrids,var,comment)
 
        endif

!mp
        if (c6_maproj .ne. 'rotlat') then
          if(c10_grid_f(1:lf).eq.'wrfsi')then
             filename=c10_grid_f(1:lf)//'.d'//cnest//'.cdl'
          else
             filename = c10_grid_f(1:lf)//'.cdl'
          endif
        else
          filename = c10_grid_f(1:lf)//'.rotlat.cdl'
        endif

	write(6,*) 'using filename: ', filename
        call s_len(filename,lfn)
        call get_directory('cdl',cdl_dir,lcdl)

	write(6,*) 'file= ', cdl_dir(1:lcdl)//filename(1:lfn)
        inquire(file=cdl_dir(1:lcdl)//filename(1:lfn),exist=exist)

        if(.not.exist)then
           print*,'error: could not find file '
     +           ,cdl_dir(1:lcdl)//filename(1:lfn)
           print*,'c10_grid_f: ',c10_grid_f(1:lf)
           istatus = 0
           return
        endif
	
        if (c6_maproj .ne. 'rotlat') then

        call check_domain(lats(1,1,ns),lons(1,1,ns)
     +,nnxp,nnyp,grid_spacing_m,1,istat_chk)

        write(6,*)'deltax = ',deltax

        if(istat_chk .eq. 1)then
            write(6,*)'check_domain: status = ',istat_chk
        else
            write(6,*)'error in check_domain: status = ',istat_chk       
        endif

	endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	write(6,*) 'call put_laps_static ', ngrids

	write(6,*) 'model= ', model

        c10_grid_fname = 'nest7grid'
        call s_len(c10_grid_fname,len_fname)

!       do zin calc (note this is last [kmax] element in data array)
        if(c10_grid_fname(1:len_fname).eq.'nest7grid')then
          avgelem=3
          zinelem=ngrids
          write(6,*)' calculating zin field in array element ',zinelem
          do i = 1,nnxp
          do j = 1,nnyp
              psa = ztopsa(data(i,j,avgelem)) ! this is the avg data (3rd element)
              data(i,j,zinelem) = (20.0 - ((psa - 100.0) * 0.02))
          enddo ! j
          enddo ! i
        endif

        call put_laps_static(grid_spacing_m,model,comment,var
     1,data,nnxp,nnyp,ngrids,ngrids,std_lat,std_lat2,std_lon
     1,c6_maproj,deltax,deltay,mdlat,mdlon,lli,llj,iuri,iurj
     1,iparent_id,iratio_to_parent,nest,c10_grid_fname)

        istatus = istat_chk

        itstatus=ishow_timer()
        print*,' -------------------------------'
        print*,' total elapsed time: ', itstatus
        print*,' -------------------------------'

        deallocate (data)

	return
	end
c
c ================================================================
c
      subroutine geodat(n2,n3,erad,rlat,wlon1,xt,yt,deltax,deltay
     1 ,datr,dats,datln,datlt,ofn,wvln,silwt,which_data,maxdatacat
     1 ,istat_files)

      include 'trigd.inc'

      implicit none


      integer n2,n3
      integer lbf
      integer lb,mof,np,niq,njq,nx,ny,isbego,iwbego,
     1  iblksizo,no,iodim,istat_files,maxdatacat

      real erad,rlat,wlon1,deltax,deltay,wvln,silwt

c wvln = toptwvl_parm_wrf fror wrfsi.nl section &hgridspec
c siwt = silavwt_parm_wrf      "                

      real datr(n2,n3)
      real dats(n2,n3,maxdatacat)
      real datln(n2,n3)
      real datlt(n2,n3)

c this parameter can be increased if we ever read data with
c resolution finer than 30 sec, or if the tilesize for 30s
c data becomes greater than 10x10 deg.

      parameter(iodim=5800000)

      real xt(n2),yt(n3)
      real deltallo,deltaxq,deltayq,
     1  deltaxp,deltayp

      real rsoff,rwoff

      real std_lon
      integer istatus
      integer lcat

      character*(*) ofn
      character*180 title
      logical which_data
c
c *********************
      nx=n2-1
      ny=n3-1
c ****************************

      lcat = 1
      lbf=index(ofn,' ')-1
      title=ofn(1:lbf)//'header'


      lb=index(title,' ')-1

      call jclget(29,title(1:lb),'formatted',1,istat_files)
      if(istat_files .ne. 1)then
         write(6,*)' warning in gridgen_model opening header: '
     1,' geog path = ', title(1:lb)
         return
      endif

      read(29,*)iblksizo,no,isbego,iwbego,rsoff,rwoff
c 2    format(4i5,2(f10.8))
      print *,'title=',title
      print *,'rsoff,rwoff = ',rsoff,rwoff
      print *,'isbego,iwbego =',isbego,iwbego
      print *,'iblksizo,no =',iblksizo,no
      close(29)

      if(no .gt. 1)then
         deltallo=float(iblksizo)/float(no-1)
      elseif(no .eq. 1)then
         deltallo=float(iblksizo)/float(no)
      else
         print*,'header value no = 0'
         return
         istat_files = 0
      endif

      mof=iodim/(no*no)

      if(ofn(lbf:lbf).eq.'g'.or.ofn(lbf:lbf).eq.'a')then
         lcat=12
         if(no.eq.1250)then
            mof=1
         endif
      endif
c sg97 mof determines the number of files held in buffer while reading
c sg97 dem data; it saves some time when buffer data can be used instead
c sg97 of reading dem file again. originally mof was 4.
      if (mof.gt.10) mof=5
      
      deltaxq=0.5*wvln*deltax
      deltayq=0.5*wvln*deltay
      print *,'deltaxq,deltayq=',deltaxq,deltayq
      np=min(10,max(1,int(deltaxq/(deltallo*111000.))))
      print *,' np=',np
      deltaxp=deltaxq/float(np)
      deltayp=deltayq/float(np)
      niq=int(float(nx)*deltax/deltaxq)+4
      njq=int(float(ny)*deltay/deltayq)+4
c
      call get_standard_longitude(std_lon,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error calling laps routine'
          stop 
      endif

      call sfcopqr(no,mof,np,niq,njq,n2,n3,lcat
     +,xt,yt,90.,std_lon,erad,rwoff,rsoff
     +,deltallo,deltaxp,deltayp,deltaxq,deltayq
     +,iblksizo,isbego,iwbego,datr,dats,datln,datlt
     +,ofn,wvln,silwt,which_data,maxdatacat,istat_files)       
      return
      end


c
c     ******************************************************************
c
      subroutine sfcopqr(no,mof,np,niq,njq,n2,n3,lcat
     +          ,xt,yt,rlat,wlon1,erad,rwoff,rsoff
     +          ,deltallo,deltaxp,deltayp,deltaxq,deltayq
     +          ,iblksizo,isbego,iwbego,datr,dats,datln,datlt
     +          ,ofn,wvln,silwt,dem_data,maxdatacat,istat_files)

c js: removed dato array from subroutine argument list
c js: added rwoff/rsoff - west and south offset of tile data

      real,  allocatable ::  dato(:,:,:,:)    !dato(no,no,mof,lcat)

      real,  allocatable ::  datp(:,:,:)
      real,  allocatable ::  datq(:,:)
      real,  allocatable ::  datqs(:,:,:)
      real,  allocatable ::  datsm(:,:)
      real,  allocatable ::  datsmx(:,:)
      real,  allocatable ::  datsln(:,:)
      real,  allocatable ::  datslt(:,:)

      real,  intent(inout) ::  datln(n2,n3)
      real,  intent(out)   ::  datlt(n2,n3)
      real,  intent(out)   ::  datr(n2,n3)
      real,  intent(out)   ::  dats(n2,n3,maxdatacat)

      real iso(mof),iwo(mof),xt(n2),yt(n3),rlat,wlon1,
     +     erad,deltallo,deltaxp,deltayp,deltaxq,deltayq,
     +     wvln,silwt,xq,yq,xp,yp,xcentr,ycentr,glatp,               ! pla,plo,
     +     glonp,rio,rjo,wio2,wio1,wjo2,wjo1,xq1,yq1,
     +     rwoff,rsoff

      real r_missing_data
      real xr,yr,rval,sh,sha,rh,rha,rhn,rht,shn,sht
      real shln,shlt,rhln,rhlt
      real delta_ln(np,np),delta_lt(np,np)

c     real xpmn,xpmx,ypmn,ypmx
      real xp1,xp2,yp1,yp2
      real xpcentr,ypcentr

      real  pctcat(maxdatacat)

      integer lp
      integer ixr,iyr
      integer lent

      character*180 ofn
      character*180 title3,title3_last_read,title3_last_inquire
      character*3   title1
      character*4   title2
      character*10  cdatatype

      logical l1,l2,dem_data,l_string_contains

      data icnt/0/
      save icnt
c
      print *,'no,mof,np,niq,njq=',no,mof,np,niq,njq

      istat_files = 1

      nono=no*no
      xcentr=0.5*(xt(1)+xt(n2))
      ycentr=0.5*(yt(1)+yt(n3))
      print *,xt(1),xt(n2),xcentr
      print *,'deltaxp=',deltaxp
      nofr=0
      do 11 iof=1,mof
         iso(iof)=0
         iwo(iof)=0
  11  continue

      title3_last_read    = '/dev/null'
      title3_last_inquire = '/dev/null'

      lcat=1
      len=index(ofn,' ')
      if(ofn(len-1:len-1).eq.'v'.or.
     &   ofn(len-1:len-1).eq.'i')then
         icnt = 0
         cdatatype='landuse'
         if(ofn(len-1:len-1).eq.'i')cdatatype='islope'
      elseif(ofn(len-1:len-1).eq.'o')then
         icnt = 0
         cdatatype='soiltype'
      elseif(ofn(len-1:len-1).eq.'u' .or.
     &       ofn(len-1:len-1).eq.'h')then
         cdatatype='topography'
      endif

      print*,'sfcopqr: cdatatype = ',cdatatype

      call s_len(cdatatype,lent)

      allocate(dato(no,no,mof,lcat))

      allocate (datp(np,np,lcat),
     &          datq(niq,njq),
     &          datsm(niq,njq),
     &          datsmx(niq,njq),
     &          datsln(niq,njq),
     &          datslt(niq,njq),
     &          datqs(niq,njq,maxdatacat))

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'failed to get r_missing_data'
         return
      endif

      do 15 jq=1,njq
         print *,'jq,njq,niq,nofr=',jq,njq,niq,nofr
         do 16 iq=1,niq

            xq=(float(iq)-0.5*float(niq+1))*deltaxq+xcentr
            yq=(float(jq)-0.5*float(njq+1))*deltayq+ycentr

            xpmn=1.0e30
            ypmn=1.0e30
c           xpmx=-1.0e30
c           ypmx=-1.0e30

            do 17 jp=1,np
               do 18 ip=1,np

                  xp=xq+(float(ip)-0.5*float(np+1))*deltaxp
                  yp=yq+(float(jp)-0.5*float(np+1))*deltayp

                  call xy_to_latlon(xp,yp,erad,glatp,glonp) 

                  glatp = max(-89.9999,min(89.9999,glatp - rsoff))
                  glonp = glonp - rwoff
                  if(glonp.ge.180.) glonp = glonp - 360.
                  if(glonp.le.-180.) glonp = glonp + 360.

c                 print *,'rlat,wlon1=',rlat,wlon1

                  isoc=(int((glatp-float(isbego))/float(iblksizo)
     &          +200.)-200)*iblksizo+isbego
            iwoc=(int((glonp-float(iwbego))/float(iblksizo)
     &          +400.)-400)*iblksizo+iwbego

                  do 19 iofr=1,nofr
                     jofr=iofr
                     if(iso(iofr).eq.isoc.and.iwo(iofr).eq.iwoc)go to 10
 19                 continue
                  isocpt=abs(isoc)/10
                  isocpo=abs(isoc)-isocpt*10
                  iwocph=abs(iwoc)/100
                  iwocpt=(abs(iwoc)-iwocph*100)/10
                  iwocpo=abs(iwoc)-iwocph*100-iwocpt*10
                  if(isoc.ge.0)then
                     write(title1,'(2i1,a1)')isocpt,isocpo,'n'
                  else
                     write(title1,'(2i1,a1)')isocpt,isocpo,'s'
                  endif

                  if(iwoc.ge.0 
     1               .and. iwoc .ne. 180                    ! 1998 steve albers
     1                                      )then
                     write(title2,'(3i1,a1)')iwocph,iwocpt,iwocpo,'e'
                  else
                     write(title2,'(3i1,a1)')iwocph,iwocpt,iwocpo,'w'
                  endif

                  lb=index(ofn,' ')-1
                  title3=ofn(1:lb)//title1//title2
                  lb=index(title3,' ')-1

                  if(title3 .ne. title3_last_inquire)then
                     inquire(file=title3(1:lb),exist=l1,opened=l2)
                     title3_last_inquire = title3
                  endif

                  if(.not.l1)then
                     iwrite = 0

                     if(icnt .le. 100)then ! reduce the output
                         iwrite=1

                     elseif(icnt .le. 1000)then
                         if(icnt .eq. (icnt/100)*100)iwrite=1

                     elseif(icnt .le. 10000)then
                         if(icnt .eq. (icnt/1000)*1000)iwrite=1

                     elseif(icnt .le. 100000)then
                         if(icnt .eq. (icnt/10000)*10000)iwrite=1

                     else
                         if(icnt .eq. (icnt/100000)*100000)iwrite=1

                     endif

                     if(iwrite .eq. 1)then
                        if(l_string_contains(title3(1:lb),
     1                                       'world_topo_30s',
     1                                       istatus)             )then       
                           print*, ' error: ',title3(1:lb)
     1                            ,' does not exist ',icnt

                        elseif(l_string_contains(title3(1:lb),
     1                                           'topo_30s',
     1                                       istatus)             )then       
                           print*, ' topo_30s file ',title3(1:lb)
     1                            ,' does not exist, using topo_10m '
     1                            ,icnt

                        else ! generic warning message
                           print*, ' warning: ',title3(1:lb)
     1                            ,' does not exist ',icnt

                        endif

                     endif ! iwrite

                     icnt = icnt + 1

c initialize these arrays as they may have some garbage in them
c if we don't actually read in any data.
c
                     datp(ip,jp,:) = 0.
                     delta_ln(ip,jp) = 0.
                     delta_lt(ip,jp) = 0.
                     istat_files = 0
                     go to 20

                  endif

                  if(nofr.ge.mof)then
                     do 21 iof=1,mof
                        iso(iof)=0
                        iwo(iof)=0
21                    continue
                     nofr=0
                  endif
                  nofr=nofr+1
                  jofr=nofr

!                 read the tile
                  if(title3 .ne. title3_last_read)then
                    if( (ofn(len-1:len).eq.'u').and.(no.eq.1200).or.
     .                   no.eq.1201 )then
                         if(no.eq.1201)then
                            print*,'reading ', title3(1:lb)
                            call read_dem(29,title3(1:lb),no,no,4,4, ! topo_3s experimental
     .                              dato(1,1,nofr,1),istat)
                         else
                            call read_dem(29,title3(1:lb),no,no,2,2, ! world topo_30s
     .                              dato(1,1,nofr,1),istat)
                         endif
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'o') )then      ! soiltype top and bot layer
                      call read_dem(29,title3(1:lb),no,no,1,4,
     .                              dato(1,1,nofr,1),istat)
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'v') )then      ! world usgs 30s landuse
                      call read_dem(29,title3(1:lb),no,no,1,4,
     .                              dato(1,1,nofr,1),istat)
                      dem_data=.true.
                    elseif((ofn(len-1:len).eq.'i'))then      ! only islope in this code section
                      call read_dem_g(29,title3(1:lb),no,no,1,lcat
     .                     ,nofr, 1,4, dato,istat)
                      dem_data=.true.
                    elseif( (ofn(len-1:len).eq.'t') )then      ! soiltemp - obsolete in this code
                      call read_dem(29,title3(1:lb),no,no,2,2,
     .                              dato(1,1,nofr,1),istat)
                      dem_data=.true.
                    else                                       ! other
                      call jclget(29,title3(1:lb),'formatted',0,istatus)      
                      call vfirec(29,dato(1,1,nofr,1),nono,'lin')
                      if ((ofn(len-1:len).eq.'u').and.(no.eq.121)) then
                        dem_data=.false.                       ! topo_30s
                      endif
                    endif

                    if(istat.ne.0)then
                       print*,'error returned: sfcopqr: read_dem'
                       return
                    endif


                    title3_last_read = title3

c                   print *,'nofr,dato=',nofr,dato(1,1,nofr)
                    close(29)

                  else
                    write(6,*)' we have made the code more efficient'

                  endif ! is this a new file we haven't read yet?

                  iso(nofr)=isoc
                  iwo(nofr)=iwoc
10		  continue

                  rio=(glonp-float(iwoc))/deltallo+1.
                  rjo=(glatp-float(isoc))/deltallo+1.

!                 prevent bounds error (steve albers)
                  if(rio .lt. 1.0)then
                      if(rio .gt. 0.98)then
                          write(6,*)' reset rio for machine epsilon'      
                          rio = 1.0
                      elseif(rio .lt. 0.5)then
                          write(6,*)' error: rio out of bounds',rio
                          stop
                      endif
                  endif

                  if(rjo .lt. 1.0)then
                      if(rjo .gt. 0.98)then
                          write(6,*)' reset rjo for machine epsilon'      
                          write(6,*)jq,iq,
     1                          ip,jp,io1,jo1,jofr,rio,rjo,glatp,isoc
                          rjo = 1.0
                      elseif(rjo .lt. 0.5)then
                          write(6,*)' error: rjo out of bounds',rjo
                          write(6,*)jq,iq,
     1                          ip,jp,io1,jo1,jofr,rio,rjo,glatp,isoc
                          stop
                      endif
                  endif

c interp ok for continuous data such as topography

                  if(cdatatype.eq.'topography')then

                   io1=int(rio)
                   jo1=int(rjo)
                   io2=io1+1
                   jo2=jo1+1
                   wio2=rio-float(io1)
                   wjo2=rjo-float(jo1)
                   wio1=1.0-wio2
                   wjo1=1.0-wjo2

                   do lp = 1,lcat

                   datp(ip,jp,lp)=wio1*(wjo1*dato(io1,jo1,jofr,lp)
     +                                 +wjo2*dato(io1,jo2,jofr,lp))
     +                           +wio2*(wjo1*dato(io2,jo1,jofr,lp)
     +                                 +wjo2*dato(io2,jo2,jofr,lp))

!s & w-facing slopes > 0.
                   delta_ln(ip,jp)=
     .           ((dato(io2,jo1,jofr,lp)-dato(io1,jo1,jofr,lp))+
     .            (dato(io2,jo2,jofr,lp)-dato(io1,jo2,jofr,lp)))*.5

                   delta_lt(ip,jp)=
     .           ((dato(io1,jo2,jofr,lp)-dato(io1,jo1,jofr,lp))+
     .            (dato(io2,jo2,jofr,lp)-dato(io2,jo1,jofr,lp)))*.5

                   enddo !lp = 1,lcat

                  else

c nearest grid point for landuse and soiltype

                   io1=nint(rio)
                   jo1=nint(rjo)
                   do lp = 1,lcat
                    datp(ip,jp,lp)= dato(io1,jo1,jofr,lp)
                   enddo

                  endif ! cdatatype eq topography.
                   
20               continue
18             continue ! ip
17           continue ! jp


!           print*,'xpmx/xpmn//ypmx/ypmn/ ',xpmx,xpmn,ypmx,ypmn

! calculate average and silhouette terrain, then apply silwt weight

            if(cdatatype(1:lent).eq.'topography')then

             sha=0.
             rha=0.
             rhln=0.
             rhlt=0.
             shmax=0.

             do 22 jp=1,np
               sh=0.
               rh=0.
               rhn=0.
               rht=0.
               do 23 ip=1,np
!                 test for missing - then go to 16?
                  sh=max(sh,datp(ip,jp,1)) 
                  rh=rh+datp(ip,jp,1)
                  rhn=rhn+delta_ln(ip,jp)
                  rht=rht+delta_lt(ip,jp)
23             continue ! ip
               sha=sha+sh/(2.*float(np))
               rha=rha+rh
               rhln=rhln+rhn
               rhlt=rhlt+rht
               shmax=max(shmax,sh)
22           continue ! jp
 
             rha=rha/float(np*np)
             rms=0.0
             do 24 ip=1,np 
               sh=0.
               do 25 jp=1,np
                  sh=max(sh,datp(ip,jp,1))
                  rms=rms+((datp(ip,jp,1)-rha)*(datp(ip,jp,1)-rha))
25             continue ! jp
               sha=sha+sh/(2.*float(np))
24           continue ! ip

             datqs(iq,jq,1)=sqrt(rms/float(np*np))
             datq(iq,jq)=sha*silwt+rha*(1.-silwt)
             datsm(iq,jq)=rha                           !mean value of points used for iq,jq
             datsmx(iq,jq)=shmax                        !max value from points used for iq,jq
             datsln(iq,jq)=rhln/float(np*np)/deltaxp
             datslt(iq,jq)=rhlt/float(np*np)/deltayp

c            print *,'datq=',datq(iq,jq)

            elseif(cdatatype(1:lent).eq.'islope'    .or.
     &             cdatatype(1:lent).eq.'landuse'   .or.
     &             cdatatype(1:lent).eq.'soiltype'    )then

             call compute_categories(cdatatype,np*np,datp(1,1,1)
     &               ,maxdatacat,domcat,pctcat)
             datq(iq,jq)=domcat 
             datqs(iq,jq,:)=pctcat(:)

            elseif(cdatatype(1:lent).eq.'greenfrac'   )then

c dominant greenness fraction for each month

             do lp=1,lcat
              call compute_categories(cdatatype,np*np,datp(1,1,lp)
     &                               ,1,domcat,pctcat)
              datqs(iq,jq,lp)=domcat
             enddo

            endif

16       continue ! iq
15    continue ! jq

      print *,'after 15'
 
      xq1=(1.-0.5*float(niq+1))*deltaxq+xcentr
      yq1=(1.-0.5*float(njq+1))*deltayq+ycentr

      if(cdatatype(1:lent).eq.'topography')then

        print*
        print*,'before gdtost2'
        print*,'--------------'
        print*,'datq(1,1)/(niq,njq)= ',datq(1,1),datq(niq,njq)
        print*,'datqs(1,1,1)/(niq,njq)= ',datqs(1,1,1),datqs(niq,njq,1)
        print*,'datsln(1,1)/(niq,njq)= ',datsln(1,1),datsln(niq,njq)
        print*,'datslt(1,1)/(niq,njq)= ',datslt(1,1),datslt(niq,njq)
        print*,'mean/max topo at iq,jq (1,1)/(niq,njq): '
     +,datsm(1,1),datsmx(1,1),datsm(niq,njq),datsmx(niq,njq)

        do 28 jr=1,n3
         do 29 ir=1,n2
           xr=(xt(ir)-xq1)/deltaxq+1.
           yr=(yt(jr)-yq1)/deltayq+1.

           call gdtost2(datq,niq,njq,xr,yr,rval)
           datr(ir,jr)=max(0.,rval)
           if( datr(ir,jr).gt.30000. )then
               print*,'warning: value out of bounds'
           endif    

           call gdtost2(datqs,niq,njq,xr,yr,rval)
           dats(ir,jr,1)=max(0.,rval)
           call gdtost2(datsln,niq,njq,xr,yr,rval)
           datln(ir,jr)=rval
           call gdtost2(datslt,niq,njq,xr,yr,rval)
           datlt(ir,jr)=rval

 29      continue
 28     continue

        print*,'after gdtost2'
        print*,'-------------'
        print*,'datr(1,1)/(n2,n3)= ',datr(1,1),datr(n2,n3)
        print*,'dats(1,1,1)/(n2,n3)= ',dats(1,1,1),dats(n2,n3,1)
        print*,'datln(1,1)/(n2,n3)= ',datln(1,1),datln(n2,n3)
        print*,'datlt(1,1)/(n2,n3)= ',datlt(1,1),datlt(n2,n3)
 
      elseif(cdatatype(1:lent).eq.'landuse'.or.
     +       cdatatype(1:lent).eq.'islope' .or.
     +       cdatatype(1:lent).eq.'soiltype')then

        do 38 jr=1,n3
         do 39 ir=1,n2
            ixr=nint((xt(ir)-xq1)/deltaxq)+1.
            iyr=nint((yt(jr)-yq1)/deltayq)+1.
            if(ixr.lt.1)ixr=1
            if(iyr.lt.1)iyr=1
            if(ixr.gt.n2)ixr=niq
            if(iyr.gt.n3)iyr=njq

            datr(ir,jr)=  datq(ixr,iyr)     !dominant category
            dats(ir,jr,:)=datqs(ixr,iyr,:)  !percent dist for ea category

 39      continue
 38     continue

      endif

      deallocate(dato)
      deallocate(datp,
     &           datq,
     &           datqs,
     &           datsm, 
     &           datsmx,
     &           datsln,
     &           datslt)

      return
      end

      subroutine vfirec(iunit,a,n,type)
      character*1 vc
      character*(*) type
      common/vform/vc(0:63)
      character line*80, cs*1
      dimension a(*)

      if(vc(0).ne.'0') call vfinit

      ich0=ichar('0')
      ich9=ichar('9')
      ichcz=ichar('z')
      ichlz=ichar('z')
      ichca=ichar('a')
      ichla=ichar('a')
      
      read(iunit,10)nn,nbits,bias,fact
 10   format(2i8,2e20.10)
      if(nn.ne.n) then
         print*,' word count mismatch on vfirec record '
         print*,' words on record - ',nn
         print*,' words expected  - ',n
         stop 'vfirec'
      endif

      nvalline=(78*6)/nbits
      nchs=nbits/6
      do 20 i=1,n,nvalline
         read(iunit,'(a78)') line
         ic=0
         do 30 ii=i,i+nvalline-1
            isval=0
            if(ii.gt.n) go to 20
            do 40 iii=1,nchs
               ic=ic+1
               cs=line(ic:ic)
               ics=ichar(cs)
               if(ics.le.ich9)then
                  nc=ics-ich0
               elseif(ics.le.ichcz) then
                  nc=ics-ichca+10
               else
                  nc=ics-ichla+36
               endif
               isval=ior(ishft(nc,6*(nchs-iii)),isval)
 40         continue
            a(ii)=isval
 30      continue
 20   continue

      facti=1./fact
      if(type.eq.'lin') then
         do 48 i=1,n
            a(i)=a(i)*facti-bias
 48      continue
      elseif(type.eq.'log') then
         scfct=2.**(nbits-1)
         do 55 i=1,n
            a(i)=sign(1.,a(i)-scfct)
     +           *(10.**(abs(20.*(a(i)/scfct-1.))-10.))
 55      continue
      endif

      return
      end
c
cc ------------------------------------------------------------------
c
      subroutine vfinit                                                  
      character*1vc,vcscr(0:63)                                         
      common/vform/vc(0:63)                                             
      data vcscr/'0','1','2','3','4','5','6','7','8','9'                 
     +,'a','b','c','d','e','f','g','h','i','j'                          
     +,'k','l','m','n','o','p','q','r','s','t'                          
     +,'u','v','w','x','y','z','a','b','c','d'                          
     +,'e','f','g','h','i','j','k','l','m','n'                          
     +,'o','p','q','r','s','t','u','v','w','x'                          
     +,'y','z','{','|'/                                                 
                                                                        
      do10n=0,63                                                        
      vc(n)=vcscr(n)                                                    
  10  continue                                                          
                                                                        
      return                                                            
      end
c +------------------------------------------------------------------+
      function intlshft(iword,nshft)
c
c       this function shifts iword to the left nshft bits in a
c         circular manner.
c
      intlshft=ishft(iword,nshft)
      return
      end
c +------------------------------------------------------------------+
      function intor(iword1,iword2)
c
c       this function performs a bit-by-bit or between iword1 and
c         iword2.
c
      intor=ior(iword1,iword2)
      return
      end


      subroutine binom2(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
      implicit none
      real x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy,
     +   wt1,wt2,yz22,yz23,yz24,yz11,yz12,yz13,yoo
      integer istend
c      common/bin/itypp,i0x,i1x,i2x,yoo
       yyy=1e30
       if(x2.gt.1.e19.or.x3.gt.1.e19.or.
     +   y2.gt.1.e19.or.y3.gt.1.e19)return
      wt1=(xxx-x3)/(x2-x3)
      wt2=1.0-wt1
      istend=0
      if(y4.lt.1.e19.and.x4.lt.1.e19) go to 410
      yz22=wt1
      yz23=wt2
      yz24=0.0
      istend= 1
410   if(y1.lt.1.e19.and.x1.lt.1.e19) go to 430
      yz11=0.0
      yz12=wt1
      yz13=wt2
      if(istend.eq.1)go to 480
      go to 450
430   yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
      yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
      yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
      if(istend.eq.  1    ) go to 470
450   yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
      yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
      yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
470   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)+wt2*(yz22*y2+yz23*y3+yz24*y4)
       go to 490
480      yyy=wt1*y2+wt2*y3
490   yoo=yyy
      return
      end
c
c determine dominant category 1-05-01 js
c

      subroutine compute_categories(ctype,nnp,data,nlcat,domcat
     +,pctcat)

      implicit none

      character*(*) ctype

      integer nnp
      integer igc
      integer i,j,k
      integer nlcat
      integer maxcat
      integer lcat(nlcat)
      real    pctcat(nlcat)
      real    data(nnp)
      real    domcat
      real    sum_g

c categorical data types
      if(ctype.eq.'landuse'   .or.
     +   ctype.eq.'soiltype'  .or.
     +   ctype.eq.'islope'  )then

         do k=1,nlcat
            lcat(k)=0
            pctcat(k)=0.
         enddo

         do i=1,nnp
         do k=1,nlcat
            if(nint(data(i)).eq.k)then
               lcat(k)=lcat(k)+1
            endif
         enddo
         enddo

         maxcat=-1
         do k=1,nlcat
            pctcat(k)=lcat(k)/float(nnp)
            if(lcat(k).gt.maxcat)then
               maxcat=lcat(k)
               domcat=float(k)
            endif
         enddo
         if(ctype.eq.'landuse')then
          if(pctcat(16).ge.0.5) then                            !!jbresch
	    domcat = 16                                         !!jbresch
          else if(domcat.eq.16.and.pctcat(16).lt.0.5)then       !!jbresch
            maxcat=-1
            do k=1,nlcat
               if(k.ne.16)then
                  if(lcat(k).gt.maxcat)then
                     maxcat=lcat(k)
                     domcat=float(k)
                  endif
               endif
            enddo
          endif
         endif

c quantitative data types
      elseif(ctype.eq.'greenfrac'.or.ctype.eq.'soiltemp')then

         sum_g=0.
         igc=0
         do i=1,nnp
          if(data(i).gt.0.0)then
             sum_g=sum_g+data(i)
             igc = igc+1
          endif
         enddo
         if(igc.gt.0)then
            domcat=sum_g/float(igc)
         else
            domcat=0.0
         endif

      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccc
	subroutine zero2d(adum2d,nnxp,nnyp)

	integer i,j,nnxp,nnyp
	real adum2d(nnxp,nnyp)
	
	do j=1,nnyp
	do i=1,nnxp
	adum2d(i,j)=0.
	enddo
	enddo

	end subroutine zero2d

cccccccccccccccccccccccccccccccccccccccccccc

	subroutine dumtodata(adum2d,nnxp,nnyp,ival,data)

	integer nnxp,nnyp,ival
	real adum2d(:,:),data(:,:,:)

	do j=1,nnyp
	do i=1,nnxp
	data(i,j,ival)=adum2d(i,j)
	enddo
	enddo

	end subroutine dumtodata

ccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine smdhld(ime,jme,h,s,lines,nsmud)
!      parameter(ime=im+2,jme=jm+4)
      dimension ihw(jme),ihe(jme)
      dimension h(ime,jme),s(ime,jme)
     &         ,hbms(ime,jme),hne(ime,jme),hse(ime,jme)
!-----------------------------------------------------------------------
          do j=1,jme
      ihw(j)=-mod(j,2)
      ihe(j)=ihw(j)+1
          enddo
!-----------------------------------------------------------------------

              do j=1,jme
          do i=1,ime
      hbms(i,j)=1.-s(i,j)
          enddo
              enddo
!
      jmelin=jme-lines+1
      ibas=lines/2
      m2l=mod(lines,2)
!
              do j=lines,jmelin
          ihl=ibas+mod(j,2)+m2l*mod(j+1,2)
          ihh=ime-ibas-m2l*mod(j+1,2)
!
          do i=ihl,ihh
      hbms(i,j)=0.
          enddo
              enddo
!-----------------------------------------------------------------------
                  do ks=1,nsmud
!-----------------------------------------------------------------------
              do j=1,jme-1
          do i=1,ime-1
      hne(i,j)=h(i+ihe(j),j+1)-h(i,j)
          enddo
              enddo
              do j=2,jme
          do i=1,ime-1
      hse(i,j)=h(i+ihe(j),j-1)-h(i,j)
          enddo
              enddo
!
              do j=2,jme-1
          do i=1+mod(j,2),ime-1
      h(i,j)=(hne(i,j)-hne(i+ihw(j),j-1)
     &       +hse(i,j)-hse(i+ihw(j),j+1))*hbms(i,j)*0.125+h(i,j)
          enddo
              enddo
!-----------------------------------------------------------------------
                      enddo
!-----------------------------------------------------------------------
      return
      end
