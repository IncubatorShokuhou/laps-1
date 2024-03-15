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
  
        subroutine put_derived_wind_prods
     1    (imax,jmax,kmax                                         ! input
     1    ,nx_l,ny_l,nz_l                                         ! input
     1    ,max_radars,r_missing_data                              ! input
     1    ,i4time_sys                                             ! input
     1    ,dbz_max_2d,istat_lps                                   ! input
     1    ,lat,lon,topo,ldf                                       ! input
     1    ,tpw_2d                                                 ! input
     1    ,heights_3d                                             ! input
     1    ,uanl,vanl)                                             ! output

!       1997 jun     ken dritz      added nx_l, ny_l, nz_l, and max_radars
!                                   as dummy arguments, making arrays
!                                   (including those in 'windparms.inc')
!                                   declared with those dimensions automatic.
!                                   also added r_missing_data as dummy
!                                   argument.  (resizability change)

        include 'windparms.inc'

!       housekeeping
        integer ss_normal,rtsys_no_data
        parameter (ss_normal        =1, ! success
     1             rtsys_no_data    =3) ! no data

        integer j_status(20),i4time_array(20)
        character*3 exts(20)

!       array variables
        real lat(imax,jmax),lon(imax,jmax),topo(imax,jmax)
        real ldf(imax,jmax)      
        real dx(nx_l,ny_l)
        real dy(nx_l,ny_l)

        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)
     1        ,wanl_2d(imax,jmax)
        
        real grid_ra_ref(imax,jmax,kmax)
        real dbz_max_2d(imax,jmax)
        real heights_3d(imax,jmax,kmax)

!       stuff for sfc and mean winds
        real umean(nx_l,ny_l),vmean(nx_l,ny_l)
        real ushear(nx_l,ny_l),vshear(nx_l,ny_l)
        real ustorm(nx_l,ny_l),vstorm(nx_l,ny_l)
        real out_lhe_3d(nx_l,ny_l,5)

!       stuff for reading radar reflectivity
        character*4 radar_name
        character*255 c_filespec
        character*31 ext_radar

!       stuff for helicity
        real helicity(nx_l,ny_l)

        real tpw_2d(nx_l,ny_l)         ! units are m
        real upslope_flux(nx_l,ny_l)         

!       dummy arrays
        real dum1_2d(nx_l,ny_l)
        real dum2_2d(nx_l,ny_l)
        real dum3_2d(nx_l,ny_l)
        integer idum1_2d(nx_l,ny_l)

!       used for steer grid
!       real iiilut(-nx_l:nx_l,-ny_l:ny_l)

!       stuff for 2d fields
        real lifted(nx_l,ny_l),liw(nx_l,ny_l)
        real field_array(nx_l,ny_l,2)

!       radar stuff
        integer  n_fcst_radar
!       parameter (n_fcst_radar = 7200 / laps_cycle_time) ! out to 2 hours
        parameter (n_fcst_radar = 0) ! no forecasts for now

        real ref_max(nx_l,ny_l,0:n_fcst_radar)

        real ref_curr_2d(nx_l,ny_l)
        real ref_fcst_2d(nx_l,ny_l)

        character*125 comment_2d,comment_a(0:10)
        character*10 units_2d,units_a(0:10)
        character*3 var_2d,var_a(0:10)

        character*31 ext

        character*40 c_vars_req
        character*180 c_values_req

        write(6,*)
        write(6,*)' entering derived wind fields subroutine',i4time_sys       

!       housekeeping
        n_prods = 4
        do i = 1,n_prods
            j_status(i) = rtsys_no_data
            i4time_array(i) = i4time_sys ! default value
        enddo ! i

        istat_lhe = 0
        istat_lmr = 0
        istat_lf1 = 0
        istat_liw = 0

        n_lhe = 1
        n_liw = 2
        n_lmr = 3
        n_lf1 = 4
!       n_vrc = 5

!       exts(n_lwm) = 'lwm'
        exts(n_lhe) = 'lhe'
        exts(n_liw) = 'liw'
        exts(n_lmr) = 'lmr'
        exts(n_lf1) = 'lf1'
!       exts(n_vrc) = 'vrc'


        i4_elapsed = ishow_timer()

!       read in 3d u/v wind data
        ext = 'lw3'
        var_2d = 'u3'
        call get_laps_3d(i4time_sys,nx_l,ny_l,nz_l
     1          ,ext,var_2d,units_2d,comment_2d,uanl,istat_lw3u)

        var_2d = 'v3'
        call get_laps_3d(i4time_sys,nx_l,ny_l,nz_l
     1          ,ext,var_2d,units_2d,comment_2d,vanl,istat_lw3v)

        if(istat_lw3u .ne. 1 .or. istat_lw3v .ne. 1)then
            write(6,*)' error reading in lw3 u/v fields - abort'
            return
        endif


        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' error getting laps_cycle_time'
            return
        endif

        write(6,*)
        write(6,*)' reading reflectivity data'

        c_vars_req = 'radarext_3d'

        call get_static_info(c_vars_req,c_values_req,1,istatus)

        if(istatus .eq. 1)then
            ext_radar = c_values_req(1:3)
        else
            write(6,*)' error getting radarext_3d'
            return
        endif

        write(6,*)' ext_radar = ',ext_radar

!       i4_tol = max(ilaps_cycle_time / 2, iradar_cycle_time / 2)

        call get_filespec(ext_radar,2,c_filespec,istatus)

        write(6,*)' istat_lps = ',istat_lps

        if(istat_lps .eq. 1)then
            write(6,*)' passing in "dbz_max_2d" for "ref_max" array'
            ref_max(:,:,0) = dbz_max_2d
            istat_radar_2dref = 1 ! since we can now use the data

        elseif(ext_radar .ne. 'vrc')then  ! original way to get column max ref
            call get_ref_base(ref_base,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error getting ref_base'
                return
            endif

            i4_tol = 1200

            call read_radar_3dref(i4time_sys,
     1                 i4_tol,i4_ret,                                   ! i/o
     1                 .true.,ref_base,
     1                 nx_l,ny_l,nz_l,ext_radar,
     1                 lat,lon,topo,.true.,.true.,
     1                 heights_3d,
     1                 grid_ra_ref,
     1                 rlat_radar,rlon_radar,rheight_radar,radar_name,     
     1                 n_ref_grids,istat_radar_2dref,istat_radar_3dref)       

            if(istat_radar_2dref .eq.  1 .or. 
     1         istat_radar_2dref .eq. -1          )then
                write(6,*)' call get_max_reflect'

                call get_max_reflect(grid_ra_ref,nx_l,ny_l,nz_l
     1                              ,ref_base,ref_max(1,1,0))

                istat_radar_2dref = 1 ! since we can now use the data

            endif

        else ! ext_radar = 'vrc' (and not lps)
            call read_radar_2dref(i4time_sys,radar_name,
     1                  nx_l,ny_l,
     1                  ref_max(1,1,0),istat_radar_2dref)

            if(istat_radar_2dref .eq. 1)then
                write(6,*)' radar 2d ref data successfully read in'
            elseif(istat_radar_2dref .eq. -1)then
                write(6,*)' radar 2d ref: fill missing data'
                do i = 1,nx_l
                do j = 1,ny_l
                    if(ref_max(i,j,0) .eq. r_missing_data)then
                        ref_max(i,j,0) = ref_base
                    endif
                enddo ! j
                enddo ! i
                istat_radar_2dref = 1 ! since we can now use the data
            else
                write(6,*)' radar 2d ref data not successfully'
     1                   ,' read in'
            endif

        endif ! ext_radar

!       calculate and write out storm steering wind field (ustorm, vstorm)

!       get layer mean wind
        if(.false.)then
            call mean_wind(uanl,vanl,topo,imax,jmax,kmax
     1                    ,umean,vmean,ustorm,vstorm,istatus)
        else
            call mean_wind_bunkers(uanl,vanl,topo,imax,jmax,kmax       ! i
     1                    ,heights_3d                                  ! i
     1                    ,umean,vmean,ushear,vshear,ustorm,vstorm     ! o
     1                    ,istatus)                                    ! o
        endif

        if(istatus .ne. 1)then
            write(6,*)' fatal error in mean_wind routine'
            return
        endif

        i4_elapsed = ishow_timer()

!       write(6,*)
!       write(6,*)' calculating storm motion grid'

!       get storm motion from mean wind and tracking info
!       call steer_grid(i4time_lapswind,imax,jmax,kmax
!    1    ,dum1_2d,dum2_2d,dum3_2d,dum4_2d,grid_ra_veldum,grid_ra_ref
!    1    ,dum5_2d,dum6_2d
!    1    ,lat,lon,standard_latdum,standard_londum
!    1                  ,iiilut,umean,vmean
!    1                  ,ustorm,vstorm,istatus)

        i4_elapsed = ishow_timer()

        if(istat_radar_2dref.eq.1)then ! get and advect the reflectivity stuff (lmr)
                                       ! note: ustorm and vstorm need to be determined

            i4_elapsed = ishow_timer()

            if(n_fcst_radar .gt. 0)then
              index = n_fcst_radar
              do ifcst = 1,index
                write(6,*)
     1            ' generating advected max reflectivity for pd ',ifcst       

                write(6,*)' grid spacing (m) = ',grid_spacing_m

                call advect(ustorm,vstorm,ref_max(1,1,0)
     1                  ,dum1_2d,grid_spacing_m,imax,jmax
     1                  ,ref_max(1,1,ifcst)
     1                  ,float(ilaps_cycle_time*ifcst),1.0,lon
     1                  ,r_missing_data)

              enddo ! i
            endif

            i4_elapsed = ishow_timer()

!           move writing of lmr product to here with appropriate changes
!           write out lmr file
            write(6,*)' writing out lmr file for max reflectivity'
            do ifcst = 0,n_fcst_radar
                minutes_10 = (ilaps_cycle_time * ifcst) / 600
                write(var_a(ifcst),101) 'r',minutes_10 
 101            format(a,i2.2)
                if(minutes_10 .lt. 10)then
!                   write(var_a(ifcst),101)minutes_10
!101                format('r0',i1)
cc                    var_a(ifcst) = 'r'
                    write(comment_a(ifcst),111)minutes_10
111                 format('laps max reflectivity  ',i1,'0 min fcst')
                else
cc                    write(var_a(ifcst),102)minutes_10
cc102                 format('r',i2)
                    write(comment_a(ifcst),112)minutes_10
112                 format('laps max reflectivity ' ,i2,'0 min fcst')
                endif
                units_a(ifcst) = 'dbz'
            enddo

            ext = 'lmr'

            call put_laps_multi_2d(i4time_sys,ext,var_a
     1      ,units_a,comment_a,ref_max,nx_l,ny_l,n_fcst_radar+1
     1                                                     ,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error writing out lmr file'
            endif

            istat_lmr = istatus

        else
            write(6,*)
     1       ' not writing out max reflectivity lmr file'

        endif

!       get and advect the reflectivity history stuff (lf1)
        if(.false.)then ! 'get_radar_max_pd' not yet modernized for radar input
            write(6,*)
            write(6,*)' generating max reflectivity history analysis'
            call get_radar_max_pd(i4time_sys-ilaps_cycle_time
     1          ,i4time_sys,imax,jmax,kmax,heights_3d,ext_radar
     1          ,max_radar_files                                         ! i
     1          ,lat,lon,topo
     1          ,ref_max(1,1,0),frac_sum,istat_radar_hist)

        else
            istat_radar_hist = 0

        endif

        if(istat_radar_hist.eq.1)then

            i4_elapsed = ishow_timer()

            if(n_fcst_radar .gt. 0)then
              index = n_fcst_radar
              do ifcst = 1,index
                write(6,*)
     1          ' generating advected max reflectivity history for pd'       
     1          ,ifcst

                call advect(ustorm,vstorm,ref_max(1,1,0)
     1                  ,dum1_2d,grid_spacing_m,imax,jmax
     1                  ,ref_max(1,1,ifcst)
     1         ,float(ilaps_cycle_time*ifcst),1.0,lon,r_missing_data)

              enddo ! ifcst
            endif

            i4_elapsed = ishow_timer()

!           write out lf1 file
            write(6,*)
     1      ' writing out lf1 file for max reflectivity history'

            do ifcst = 0,n_fcst_radar
                minutes_10 = (ilaps_cycle_time * ifcst) / 600
                write(var_a(ifcst),101) 'h',minutes_10
                if(minutes_10 .lt. 10)then
!                   write(var_a(ifcst),201)minutes_10
!201                format('h0',i1)
ccc                 var_a(ifcst) = 'h'
                    write(comment_a(ifcst),211)minutes_10
211                 format('laps max reflectivity history  ',i1
     1                    ,'0 min fcst')
                else
ccc                 write(var_a(ifcst),202)minutes_10
ccc202              format('h',i2)
                    write(comment_a(ifcst),212)minutes_10
212                 format('laps max reflectivity history ' ,i2
     1                    ,'0 min fcst')
                endif
                units_a(ifcst) = 'dbz'
            enddo

            ext = 'lf1'

            call put_laps_multi_2d(i4time_sys,ext,var_a,units_a
     1          ,comment_a,ref_max,nx_l,ny_l,n_fcst_radar+1,istatus)

            if(istatus .ne. 1)then
                write(6,*)' error writing out lf1 file'
            endif

            istat_lf1 = istatus

        else
            write(6,*)
     1      ' not writing out max reflectivity history lf1 file'

        endif


        if(.true.)then ! calculate helicity
                       ! note: (u+vstorm must be calculated first)

!           write out helicity field
            write(6,*)' calculating helicity'
            call helicity_laps(uanl,vanl,ustorm,vstorm
     1                        ,heights_3d,topo        
     1                        ,imax,jmax,kmax,helicity,istatus)

            if(istatus .eq. 1)then
                ext = 'lhe'

                var_a(0) = 'lhe'
                var_a(1) = 'mu'
                var_a(2) = 'mv'
                var_a(3) = 'shu'
                var_a(4) = 'shv'

                units_a(0) = 'm/s**2'
                units_a(1) = 'm/s'
                units_a(2) = 'm/s'
                units_a(3) = 'm/s'
                units_a(4) = 'm/s'

                comment_a(0) = 'laps helicity'
                comment_a(1) = 'laps mean wind 0-6km agl'
                comment_a(2) = 'laps mean wind 0-6km agl'
                comment_a(3) = 'u shear component 0-6km agl'
                comment_a(4) = 'v shear component 0-6km agl'

                call move(helicity,out_lhe_3d(1,1,1),nx_l,ny_l)
                call move(umean   ,out_lhe_3d(1,1,2),nx_l,ny_l)
                call move(vmean   ,out_lhe_3d(1,1,3),nx_l,ny_l)
                call move(ushear  ,out_lhe_3d(1,1,4),nx_l,ny_l)
                call move(vshear  ,out_lhe_3d(1,1,5),nx_l,ny_l)

                call put_laps_multi_2d(i4time_sys,ext,var_a
     1          ,units_a,comment_a,out_lhe_3d,nx_l,ny_l,5,istatus)
                 istat_lhe = istatus

            else
                write(6,*)' warning: lhe not written out'

            endif ! istatus

        endif

        i4_elapsed = ishow_timer()

!       calculate upslope component of moisture flux (psd conventions)
        call get_grid_spacing_array(lat,lon,nx_l,ny_l,dx,dy)
	call up_mflux(nx_l,ny_l,nz_l,topo,ldf,dx,dy
     1                     ,uanl,vanl,tpw_2d,upslope_flux
     1                     ,heights_3d,r_missing_data)

        if(.true.)then ! liw

!           calculate li * 600mb omega
            write(6,*)' generating li * omega file'

!           read in li data
            var_2d = 'li'
            ext = 'lst'
            call get_laps_2dgrid(i4time_sys,0,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                                          ,lifted,0,istatus)

!           read in omega data

!           determine k coordinate for 600mb for passing in omega (wanl array)
            lvl_liw = nint(zcoord_of_pressure(60000.))
            ipres = nint( pressure_of_level(lvl_liw) / 100. )
            write(6,*)' lvl for lifted index * omega = ',lvl_liw,ipres

            var_2d = 'om'
            ext = 'lw3'
            call get_laps_2dgrid(i4time_sys,0,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                ,wanl_2d,ipres,istat_omega)

            if(istatus .ne. 1 .or. istat_omega .ne. 1)then
                write(6,*)' error reading lifted index / om data'
                liw = r_missing_data

            else ! good li data
                call cpt_liw(lifted,wanl_2d,imax,jmax,liw)

            endif

            call move(liw    ,field_array(1,1,1),nx_l,ny_l)
!           call move(wanl_2d,field_array(1,1,2),nx_l,ny_l)
            call move(upslope_flux,field_array(1,1,2),nx_l,ny_l)        

!           write out liw field
!           note that these arrays start off with 0 as the first index
            var_a(0) = 'liw'
            var_a(1) = 'umf'
            ext = 'liw'
            units_a(0) = 'k-pa/s'
!           units_a(1) = 'pa/s  '
            units_a(1) = 'm**2/s'
            comment_a(0) = 'log laps li * 600mb omega'
!           comment_a(1) = 'laps 600mb omega'
            comment_a(1) = 'upslope component of moisture flux'

            call put_laps_multi_2d(i4time_sys,ext,var_a,units_a
     1                ,comment_a,field_array,imax,jmax,2,istatus)

            istat_liw = istatus

        endif

        if(istat_lhe .eq. 1)then
            i4time_array(n_lhe) = i4time_sys
            j_status(n_lhe) = ss_normal
        endif

        if(istat_lmr .eq. 1)then
            i4time_array(n_lmr) = i4time_sys
            j_status(n_lmr) = ss_normal
        endif

        if(istat_lf1 .eq. 1)then
            i4time_array(n_lf1) = i4time_sys
            j_status(n_lf1) = ss_normal
        endif

        if(istat_liw .eq. 1)then
            i4time_array(n_liw) = i4time_sys
            j_status(n_liw) = ss_normal
        endif


        write(6,*)' status of output products'
        do i = 1,n_prods
            write(6,*)i,' ',exts(i),' '
     1         ,i4time_array(i),j_status(i)
        enddo ! i

        return
        end
