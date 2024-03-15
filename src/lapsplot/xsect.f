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

        subroutine xsect(c_display,i4time_ref,lun,l_atms
     1                  ,standard_longitude
     1                  ,nx_l,ny_l,nz_l,nx_c,nz_c,nx_p,nx_t
     1                  ,r_missing_data,laps_cycle_time,maxstns      
     1                  ,dyn_low,dyn_high,dx,dy
     1                  ,density,plot_parms,namelist_parms,ifield_found)       

!       97-aug-14     ken dritz     added nx_l, ny_l, nz_l as dummy arguments
!       97-aug-14     ken dritz     added nx_c, nz_c as dummy arguments
!       97-aug-14     ken dritz     removed parameter declarations for
!                                   nx_c, nz_c (commented them out)
!       97-aug-14     ken dritz     added r_missing_data, laps_cycle_time
!                                   as dummy arguments
!       97-aug-14     ken dritz     added maxstns as dummy argument
!       97-aug-14     ken dritz     removed include of lapsparms.for
!       97-aug-14     ken dritz     pass nx_l, ny_l, nz_l, nx_c, nz_c
!                                   to interp_3d
!       97-aug-14     ken dritz     pass r_missing_data to interp_3d
!       97-aug-14     ken dritz     pass nx_l, ny_l, nz_l, nx_c, nz_c
!                                   to interp_3dn
!       97-aug-14     ken dritz     pass nx_l, ny_l, nz_l, nx_c, nz_c
!                                   to interp_3d_spread
!       97-aug-14     ken dritz     pass r_missing_data to interp_3d_spread
!       97-aug-14     ken dritz     pass nx_l, ny_l, nx_c, r_missing_data
!                                   to interp_2d
!       97-aug-14     ken dritz     pass maxstns to label_other_stations
!       97-aug-25     steve albers  removed /lapsplot_cmn1/
!       97-aug-25     steve albers  removed /lapsplot_cmn2/
!       97-aug-25     steve albers  removed equivalence of pcp_type_3d,rh_3d.

        include 'lapsplot.inc'

        real lat(nx_l,ny_l),lon(nx_l,ny_l),topo(nx_l,ny_l)
        real rlaps_land_frac(nx_l,ny_l)
        real dx(nx_l,ny_l)
        real dy(nx_l,ny_l)

        integer nx_c,nz_c,nz_b
!       parameter (nx_c = 61)   ! nx_l ! number of horizontal points in x-sect
!       parameter (nz_c = nz_l) ! number of vertical levels in laps

        parameter (nz_b = 5)    ! bottom level of atms x-sect

        include 'laps_cloud.inc'

        real clouds_3d(nx_l,ny_l,kcloud)

        common/lapsplot_omega/l_convert

        logical l_sta,l_convert,lapsplot_pregen,l_atms,l_pregen 
        logical l_radar_read, l_wind_read,l_arrival_gate
        logical iflag_mvd,iflag_icing_index,iflag_cloud_type
        logical iflag_snow_potential,iflag_bogus_w
        logical l_low_fill, l_high_fill
        logical l_latlon, l_parse

        data lapsplot_pregen /.true./

        real cld_pres(kcloud)

        character*2 c2_cloud_type,c2_cloud_types(0:10)
        data c2_cloud_types
     1  /'  ','st','sc','cu','ns','ac','as','cs','ci','cc','cb'/

        character*2 c2_precip_type,c2_precip_types(0:10)
        data c2_precip_types
     1  /'  ','rn','sn','zr','sl','ha','  ','  ','  ','  ','  '/

        character*1 c1_precip_types(0:10)
        data c1_precip_types
     1       /' ','r','*','z','i','h',' ',' ',' ',' ',' '/


        character*1 c_prodtype
        character*5 fcst_hhmm

        character*4 c4_log
        character*4 radar_name

        character*10 colortable

        real dum1_array(nx_l,1)
        real dum2_array(nx_l,1)
        real dum3_array(nx_l,1)
        real dum4_array(nx_l,1)

        data mode_lwc/2/

!       character*255 c_filespec_wd/'*.lw3'/
!       character*255 c_filespec_wc/'*.lco'/
!       character*255 c_filespec_wb/'*.lba'/
        character*255 c_filespec_qg
        character*255 c_filespec

        data c_filespec_qg/'user_data:*.lqo'/

        character*31 ext_wind, ext_radar

        character*3 var_2d
        character*150  directory
        character*31  ext
        character*10  units_2d
        character*125 comment_2d
        character*40 c_model
        character*200 new_dataroot

      ! used for "potential" precip type
        logical l_mask_pcptype(nx_c,1)
        integer ibase_array(nx_c,1)
        integer itop_array(nx_c,1)

        real u_3d(nx_l,ny_l,nz_l)
        real v_3d(nx_l,ny_l,nz_l)
        real field_3d(nx_l,ny_l,nz_l)
        real temp_3d(nx_l,ny_l,nz_l)
        real rh_3d(nx_l,ny_l,nz_l)
        real q_3d(nx_l,ny_l,nz_l)
        real slwc_3d(nx_l,ny_l,nz_l)
        real cice_3d(nx_l,ny_l,nz_l)
        real grid_ra_ref(nx_l,ny_l,nz_l)
        real grid_ra_vel(nx_l,ny_l,nz_l)
!       real grid_ra_rfill(nx_l,ny_l,nz_l)

        real pcp_type_3d(nx_l,ny_l,nz_l)
!       equivalence(pcp_type_3d,rh_3d)

!       real pcp_type_2d(nx_c,nz_c)
!       equivalence(pcp_type_2d,rh_2d)

        real field_2d(nx_c,nz_c)
        real fieldt_2d(nx_t,nz_c)

!       common/labs/ia(2),nc,nrep,ncrt,ilab,nulbll,sizel,sizem,sizep
        character c9_radarage*9

        real lifted(nx_l,ny_l)

!       real cloud_cvr_2d(nx_l,ny_l)
        real cloud_ceil_2d(nx_l,ny_l)
        real vis_2d(nx_l,ny_l)
        real cloud_top_2d(nx_l,ny_l)

        real clouds_vert(nx_c,kcloud)

        real cloud_ceil_1d(nx_c)
        real cloud_top_1d(nx_c)

        real k_to_c, make_rh

        common /mcolor/mini,maxi

        real xcoord(nx_c),ycoord(nx_c)

!       common /conre1/ioffp,spval,epsval,cntmin,cntmax,cntint,ioffm

        include 'icolors.inc'

        integer n_contours
        parameter (n_contours = 29)
        real factor(n_contours)
        data factor/
     1  .00001,
     1  .00002,
     1  .00005,
     1  .0001,
     1  .0002,
     1  .0005,
     1  .001,
     1  .002,
     1  .005,
     1  .01,
     1  .02,
     1  .05,
     1  .1,
     1  .2,
     1  .5,
     1  1.,
     1  2.,
     1  5.,
     1  10.,
     1  20.,
     1  50.,
     1  100.,
     1  200.,
     1  500.,
     1  1000.,
     1  2000.,
     1  5000.,
     1  10000.,
     1  20000.
     1                  /

        real lat_1d(nx_c)
        real lon_1d(nx_c)
        real snow_1d(nx_c)

        real pres_3d(nx_l,ny_l,nz_l)
        real pres_2d(nx_l,ny_l)
        real pres_1d(nx_c)

        real u_vert(nx_c,nz_c)
        real v_vert(nx_c,nz_c)
        real pres_vert(nx_c,nz_c)
        real temp_2d(nx_c,nz_c)
        real rh_2d(nx_c,nz_c)
        real heights_2d(nx_c,nz_c)
        real radar_2d(nx_c,nz_c)
        real slwc_2d(nx_c,nz_c)
        real cice_2d(nx_c,nz_c)
        real field_vert(nx_c,nz_c)
        real field_vert_buf(nx_c,nz_c)
        real field_vert_diff(nx_c,nz_c)
        real field_vert2(nx_c,nz_c)
        real field_vert3(nx_p,nx_p)
        real w_2d(nx_c,nz_c)
        real mvd_2d(nx_c,nz_c)
        integer icing_index_2d(nx_c,nz_c)
        real terrain_vert(nx_c,nz_c)
        real terrain_vert1d(nx_c)
        real lon_vert(nx_c)

        integer iarg

        real mspkt 
        data mspkt /.518/

        character*100 c_label
        character*1 c_display
        character*1 c1_string,c_wind
        character*2 c_metacode 
        character*3 c3_string,c3_sta,c3_type
        character*3 c3_ylow,c3_xlow
        character*5 c5_arrival_gate
        character*4 c4_string,c_field
        character*7 c7_string
        character*24 asc_tim_24
        character*9 a9time,c9_string
        character*20 c20_sta
        integer ity

        integer n_stations
        parameter (n_stations = 32)

        character*3 c3_sta_array(n_stations)
        data c3_sta_array
     1  /'wig','ftc','lov','elb','flg','ptv','stp',
     1         'elb','bjc','den','apa','cos','cys','lar',
     1         'lic','ako','gld','lhx','bou','kio','gxy',
     1   'eri','mhr','cp3','chl','und','okc','mkc',
     1   'ict','dsm','gri','ase'/

        real sta_lat(n_stations)
        data sta_lat
     1    /  40.29,  40.59,  40.59,  39.20,  39.36,  40.26,  39.75,
     1       39.23,  39.90,  39.75,  39.57,  38.82,  41.15,  41.32,
     1       39.18,  40.17,  39.37,  38.05,  40.01,  39.35,  40.42,
     1       40.10,  39.87,  39.95,  40.44,  40.10,  35.40,  39.12,
     1       37.65,  41.53,  40.97,  39.22/

        real sta_lon(n_stations)
        data sta_lon
     1    /-103.05,-105.14,-105.14,-104.50,-103.04,-104.87,-104.87,
     1     -104.63,-105.12,-104.87,-104.85,-104.72,-104.82,-105.68,
     1     -103.70,-103.22,-101.70,-103.52,-105.25,-104.42,-104.63,
     1     -105.03,-104.76,-105.19,-104.64,-104.34, -97.60, -94.60,
     1      -97.43, -93.65, -98.32,-106.87/

        o_k(t_k,p_pa)   =   o( t_k-273.15 , p_pa/100. )  + 273.15

        ifield_found = 0

        zoom = 1.0
        denslogthr = sqrt(10.) ** (-(log(density) / log(2.)))

!       sizem = 1.0
        sizel = 2.0

        vxmin = .10
        vxmax = .90
        vymin = .20
        vymax = .80

        vymin2 = .50 - (.50-vymin) * 1.25
        vymax2 = .50 + (vymax-.50) * 1.25

!       vymin3 = .50 - (.50-vymin) * 1.20
!       vymax3 = .50 + (vymax-.50) * 1.20

        if(vymin .eq. .10)then
!           iyl_remap = 13
!           iyh_remap = 8
            iyl_remap = 10
            iyh_remap = 10
        elseif(vymin .eq. .20)then
            ixl_remap = nint(float(nx_p-1) * .0630)
            ixh_remap = nint(float(nx_p-1) * .0630)

            if(nz_l .eq. 41)then ! e.g. rsa
                iyl_remap = nint(float(nx_p-1) * .175)
                iyh_remap = nint(float(nx_p-1) * .166)
            else
                iyl_remap = nint(float(nx_p-1) * .175)
                iyh_remap = nint(float(nx_p-1) * .166)
            endif

        else
            write(6,*)' error, invalid vymin ',vymin

        endif

        ioffm = 1 ! don't plot label stuff in conrec

        igrid=0
        idiff=0

        plot_parms%iraster = -1 ! ensure raster image plots are off
                                ! there are "background" effects preventing
                                ! the use of iraster
        plot_parms%l_discrete = namelist_parms%l_discrete

        chigh_3dwind = namelist_parms%chigh_3dwind

        lapsplot_pregen = .true.

c read in laps lat/lon and topo
        call get_laps_domain_95(nx_l,ny_l,lat,lon,topo,rlaps_land_frac
     1                      ,grid_spacing_dum,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps domain'
            return
        endif

        call get_grid_spacing_cen(grid_spacing_m,istatus)

        i_graphics_overlay = 0
        i_label_overlay = 0
        i_map = 0
        i_initialize = 0

        l_wind_read = .false.
        l_radar_read = .false.

!       if(c_display .eq. 't')call setusv_dum(2hin,16)

        rleft = 1
        right = nx_c

!       decide whether the bottom of the x-sect is at the bottom of the laps domain

!       if(l_atms)then
            topo_min = 1e10
            do j = 1,ny_l
            do i = 1,nx_l
                topo_min = min(topo_min,topo(i,j))
            enddo
            enddo

            ibottom_terrain = height_to_zcoord(topo_min,istatus)

            ibottom_terrain = max(ibottom_terrain,1)

!           ibottom_terrain = 1

            write(6,*)'    lowest displayed level = ',ibottom_terrain
            write(6,*)'    nx_p/nz_c/nz_l',nx_p,nz_c,nz_l
            write(6,*)'    iyl_remap,iyh_remap',iyl_remap,iyh_remap

            bottom = ibottom_terrain
!       else
!           bottom = 1
!       endif

        ibottom = bottom

        top = nz_c
        width = right - rleft
        r_height = top - bottom

!       this lets up plot outside the main box
!       call set(.00, 1.0, .00, 1.0, rleft - width/8., right + width/8.,
!       1                            bottom - r_height/8., top + r_height/8.,1)

        if(i_initialize .eq. 1)goto100

        i_initialize = 1

!       define segment for cross section on laps grid
80      continue

        if(.true.)then
            write(6,102)
102         format(/
     1  '    type of xsect ',
     1  ' [we, sn, xxx (azimuth-true), arr (arrival gate)]   ? ',$)
        else
            write(6,103)
103         format(/
     1 '    type of xsect ',
     1  ' [we, sn, xxx (azimuth-true)]                       ? ',$)
        endif

        if(l_atms)write(6,*)' reading x-sect type from lun = ',lun

        read(lun,1201)c3_type
1201    format(a3)

        if(l_atms)write(6,*)' just read x-sect type  = ',c3_type,lun

        l_sta = .true.
        l_arrival_gate = .false.

        if(c3_type(1:2) .eq. 'we')then
            xlow = 1.
            xhigh = nx_l

            if(.true.)then
                write(6,111)ny_l,ny_l/2+1,ny_l
111             format(/'     n-s location ',
     1        '[1 to ',i3,'; 1 = s edge, '
     1        ,i3,' = center, '
     1        ,i3,' = n edge]   or '//
     1        6x,' enter latitude of way point              or'/
     1        6x,' class:       wig,ftc,lov,elb,flg'/
     1        6x,' radiometers: ptv,stp,elb,eri'/
     1        6x,' radars:      mhr,cp3,chl,und'/
     1        6x,' saos:        bjc,den,apa,cos,cys,lar,lic,ako,gld,lhx,
     1gxy'/
     1        6x,' vors:        kiw'/
     1        ' ',6x,'mesonet:     ',
     1        'bou                                               ? ',$)

            else ! stormfest
                write(6,1110)ny_l,ny_l/2+1,ny_l
1110            format(/'     n-s location ',
     1        '[1 to ',i3,'; 1 = s edge, '
     1        ,i3,' = center, '
     1        ,i3,' = n edge]   or '//
     1        '$',6x,'saos   :     ',
     1        'okc,mkc,ict,dsm,gri                               ? ')

            endif

            read(lun,120)c3_ylow
120         format(a3)

            call upcase(c3_ylow,c3_ylow)

            l_sta = .false.

            do i = 1,n_stations
              if(c3_ylow .eq. c3_sta_array(i))then
                call latlon_to_rlapsgrid(sta_lat(i),sta_lon(i),lat,lon
     1                          ,nx_l,ny_l,xsta,ysta,istatus)

                if(xsta .lt. 1 .or. xsta .gt. nx_l .or.
     1             ysta .lt. 1 .or. ysta .gt. ny_l)then
                    write(6,*)' station is outside domain - try again...
     1'
                    goto80
                endif

                ylow = ysta
                l_sta = .true.
                i_sta = i

                pos_sta = 1. + (nx_c-1.) * (xsta-xlow)/(xhigh-xlow)

!               pos_sta = xsta

                c3_sta = c3_ylow
              endif

            enddo

            if(.not. l_sta)then
                read(c3_ylow,*,err=80)ylow
                write(6,*)
                write(6,*)'      j = ',nint(ylow)

                ylow = max(min(nint(ylow),ny_l),1)

!75              if(ylow .lt. 1 .or. ylow .gt. ny_l)then
!                    write(6,*)' grid point is outside domain - try again
!     1...'
!                    goto80
!                endif

            else
                write(6,72)sta_lat(i_sta),sta_lon(i_sta)
     1                          ,nint(ylow),nint(pos_sta)
72              format(/7x,'lat/lon ',2f8.2,' j/i =',2i4)
            endif

            yhigh = ylow

        elseif(c3_type(1:2) .eq. 'sn')then
            ylow = 1.
            yhigh = ny_l

            if(.true.)then
                write(6,112)nx_l,nx_l/2+1,nx_l
112             format(/'     e-w location ',
     1        '[1 to ',i3,'; 1 = w edge, '
     1        ,i3,' = center, '
     1        ,i3,' = e edge]   or '//
     1        6x,' class:       wig,ftc,lov,elb,flg'/
     1        6x,' radiometers: ptv,stp,elb,eri'/
     1        6x,' radars:      mhr,cp3,chl,und'/
     1        6x,' saos:        bjc,den,apa,cos,cys,lar,lic,ako,gld,lhx,
     1gxy'/
     1        6x,' vors:        kiw'/
     1        ' ',6x,'mesonet:     ',
     1        'bou                                               ? ',$)

            else ! stormfest
                write(6,1120)nx_l,nx_l/2+1,nx_l
1120            format(/'     e-w location ',
     1        '[1 to ',i3,'; 1 = w edge, '
     1        ,i3,' = center, '
     1        ,i3,' = e edge]   or '//
     1        '$',6x,'saos   :     ',
     1        'okc,mkc,ict,dsm,gri                               ? ')

            endif

            read(lun,120)c3_xlow

            call upcase(c3_xlow,c3_xlow)

            l_sta = .false.

            do i = 1,n_stations
              if(c3_xlow .eq. c3_sta_array(i))then
!               call latlon_grid(sta_lat(i),sta_lon(i),igrid,jgrid)
                call latlon_to_rlapsgrid(sta_lat(i),sta_lon(i),lat,lon
     1                          ,nx_l,ny_l,xsta,ysta,istatus)
                if(xsta .lt. 1 .or. xsta .gt. nx_l .or.
     1           ysta .lt. 1 .or. ysta .gt. ny_l)then
                    write(6,*)' station is outside domain - try again...
     1'
                    goto80
                endif

                xlow = xsta
                l_sta = .true.
                i_sta = i

                pos_sta = 1. + (nx_c-1.) * (ysta-ylow)/(yhigh-ylow)

!               pos_sta = ysta

                c3_sta = c3_xlow
              endif

            enddo

            if(.not. l_sta)then
                read(c3_xlow,*,err=80)xlow
                write(6,*)
                write(6,*)'      i = ',nint(xlow)

                xlow = max(min(nint(xlow),nx_l),1)

! 76             if(xlow .lt. 1 .or. xlow .gt. nx_l)then
!                    write(6,*)' grid point is outside domain - try again
!     1...'
!                    goto80
!                endif

            else
                write(6,82)sta_lat(i_sta),sta_lon(i_sta)
     1                          ,nint(xlow),nint(pos_sta)
82              format(/7x,'lat/lon ',2f8.2,' i/j =',2i4)
            endif

            xhigh = xlow


        else ! try to get an azimuth for the x-sect
            read(c3_type,*,err=85)azi_xsect

            write(6,113)
113         format(/'     waypoint for x-sect: '/
     1      6x,' class:       wig,ftc,lov,elb,flg'/
     1      6x,' radiometers: ptv,stp,elb,eri'/
     1      6x,' radars:      mhr,cp3,chl,und'/
     1      6x,' saos:        bjc,den,apa,cos,cys,lar,lic,ako,gld,lhx,
     1gxy'/
     1      6x,' vors:        kiw'/
     1      6x,' mesonet:     bou'/
     1      7x,'lat/lon (e.g. 40.0l,-105.l) or i,j location:'
     1                                               ,4x,'? ',$)       
            read(lun,130)c20_sta
130         format(a20)
            c3_sta = c20_sta(1:3)
            call upcase(c3_sta,c3_sta)

            l_sta = .false.

            do i = 1,n_stations
              if(c3_sta .eq. c3_sta_array(i))then
                call latlon_to_rlapsgrid(sta_lat(i),sta_lon(i),lat,lon
     1                          ,nx_l,ny_l,xsta,ysta,istatus)
                l_sta = .true.
!               i_sta = i
              endif
            enddo

            if(.not. l_sta)then ! get i,j of waypoint
              call s_len2(c20_sta,lenxy)
              l_latlon = l_parse(c20_sta(1:lenxy),'l')

              if(l_latlon)then ! xy values flagged as with "l" at the end 
                  do is = 1,lenxy
                      if(c20_sta(is:is) .eq. 'l')then
                          c20_sta(is:is) = ' '
                      endif
                  enddo 

                  read(c20_sta,*,err=85)waylat,waylon

                  call latlon_to_rlapsgrid(waylat,waylon,lat,lon
     1                              ,nx_l,ny_l,xsta,ysta,istatus)       

                  if(istatus .ne. 1)then
                      return
                  endif

                  if(xsta .lt. 1. .or. xsta .gt. float(nx_l) .or.
     1               ysta .lt. 1. .or. ysta .gt. float(ny_l)   )then      
                      write(6,*)
     1                       ' station is outside domain - try again...'
                      return
                  endif

              else ! i/j mode
                read(c20_sta,*,err=85)xsta,ysta

                write(6,86)xsta,ysta

!               convert from (0-1) space to gridpoint space
                if(xsta .ge. 0.0 .and. xsta .lt. 1.0)then
                    xsta = 1.0 + xsta * float(nx_l-1)
                endif

                if(ysta .ge. 0.0 .and. ysta .lt. 1.0)then
                    ysta = 1.0 + ysta * float(ny_l-1)
                endif

                write(6,86)xsta,ysta
86              format(/2x,'waypt x,y = ',2f7.2)

                nxsta = nint(xsta)
                nysta = nint(ysta)
                write(6,96)lat(nxsta,nysta),lon(nxsta,nysta)
96              format('  waypt lat/lon = ',2f9.3)

              endif

            else
                write(6,87)c3_sta,xsta,ysta
87              format(/2x,a3,' x,y =   ',2f7.2)

            endif

!           calculate endpoints of x-sect from waypoint and azimuth
            call xsect_endpoints
     1  (xsta,ysta,azi_xsect,xlow,ylow,xhigh,yhigh,pos_sta,istatus,
     1   nx_l,ny_l,nx_c)

            goto90

85          write(6,*)' try again...'
            goto80

90      endif ! type of x-sect


100    write(6,95)
95     format(
     1  /'  field - (add "i" for image):'
     1  /
     1  /'          [wind: di,sp,u,v,om,dv,vo,pv,va,vc (barbs)'
     1  /
     1  /'           temp: [t,pt,tb,pb] (t, theta, t blnc, theta blnc)'       
     1  /
     1  /'           ht: [ht] (height)'       
     1  /
     1  /'           humid: sh,rh,rl (specific/relative humidity)'
     1  /
     1  /'           ts (thetae sat), tw (wetbulb), td (dewpoint)'
     1  /
     1  /'           cg/cf (3d cloud image),  tc (cloud type),  '
     1  ,'tp (precip type)'
!       1 /'           la (lwc - adiabatic),         lj (lwc - adjusted)'
!       1 /'                                         sj (slwc - adjusted)'
     1  /'           ls (cloud liquid)' 
                                          ! ,        ss (slwc - smith-feddes)'
     1  /'           ci (cloud ice)'
     1  /
     1  /'           ix (icing index)    mv (mean vol diam)'
     1  /'           pc/rn/sn/ic (precip/rain/snow/ice concentration)'
     1  /
     1  /'           cv (cloud cover contours)'
     1  /'           rf (ref-graphic), ri (ref-image), rv (ref-obs)]'
     1  //'     difference field: [df]  '
     1  /' ',49x,'q (quit/display)]   ? ',$)

        nulbll = 3 ! for conrec (number of lines between labels)

        read(lun,1301)c_field
1301    format(a3)

        istatus = 1

        if(c_field(1:2) .eq. 'q ')goto9999

        plot_parms%color_power = 1.0

!       c4_log = 'x '//c_field
!       if(lun .eq. 5)call logit(c4_log)

        write(6,*)' generating cross section'

        i_image = 0

        scale = 1.0                                      ! default

        call s_len(c_field,len_field)
        
        if(c_field(len_field:len_field) .eq. 'i' 
     1                    .and. c_field .ne. 'ri'
     1                    .and. c_field .ne. 'ci'
     1                    .and. c_field .ne. 'dii')then
            i_image = 1
            c_field = c_field(1:len_field-1)
            colortable = 'hues'                          ! default
        endif

        call interp_2d(lat,lat_1d,xlow,xhigh,ylow,yhigh,
     1                 nx_l,ny_l,nx_c,r_missing_data)
        call interp_2d(lon,lon_1d,xlow,xhigh,ylow,yhigh,
     1                 nx_l,ny_l,nx_c,r_missing_data)

        call interp_2d(topo,terrain_vert1d,xlow,xhigh,ylow,yhigh,
     1                 nx_l,ny_l,nx_c,r_missing_data)

        if(c_field(1:2) .eq. 'fc')then ! force config with new dataroot
          write(6,*)' enter new dataroot:'
          read(lun,17)new_dataroot
 17       format(a)
          call s_len(new_dataroot,lenroot)
          call force_get_laps_config(new_dataroot(1:lenroot),istatus)
          if(istatus .ne. 1)then
            write(6,*)' bad status returned from force_laps_config'
            return
          endif
          call get_lapsplot_parms(namelist_parms,istatus)
          if(c_field(1:3) .eq. 'fcf')then ! montage case generally

!-------------splice

!           contour in the terrain surface
            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)
            n_div = 20
            call setusv_dum(2hin,3)

!           read in sfc pressure
            i4time_tol = 43200
            var_2d = 'ps'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,i4time_tol,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,pres_2d,0,istatus)
            if(istatus .ne. 1)then
                write(6,*)' error reading surface pressure analysis'
                write(6,*)
     1        ' converting terrain to sfc pressure with std atmosphere'            
                istat_sfc_pres = 0
            else
                call interp_2d(pres_2d,pres_1d,xlow,xhigh,ylow,yhigh,
     1                         nx_l,ny_l,nx_c,r_missing_data)
                istat_sfc_pres = 1
            endif

            do i = 1,nx_c
                xcoord(i) = i
                if(istat_sfc_pres .eq. 1)then
                    ycoord(i) = max(zcoord_of_pressure(pres_1d(i)),1.0)
                else
                    ycoord(i) = 
     1              max(height_to_zcoord(terrain_vert1d(i),istatus),1.0)      
                endif

                if(i .gt. 1)then
                    do j = 0,n_div
                        frac = float(j) / float(n_div)
                        ybottom = bottom
                        ytop = ycoord(i-1) * (1.-frac) + ycoord(i)*frac
                        xval = float(i-1) + frac
                        call line(xval,ybottom,xval,ytop)
                    enddo ! j
                endif
            enddo ! i

            call setusv_dum(2hin,7)

!-----------end of splice
!           call make_xsect_labels(  vxmin, vxmax, vymin, vymax, 
!    1                               rleft, right, bottom, top,
!    1                               width, vymin2, vymax2, i_map, 
!    1                               ibottom, r_height,
!    1                               xlow,xhigh,ylow,yhigh,
!    1                               lat_1d, lon_1d, nx_c, nz_c,
!    1                               nx_l, ny_l, 1)

!           splice start for labeling other stations

!              this lets up plot outside the main box
               call set(.00, 1.0, vymin2 , vymax2
     1                , rleft - width/8., right + width/8.
     1                , bottom - r_height/8., top + r_height/8.,1)

               xsta = -10000.

               i4time_label = i4time_ref/laps_cycle_time*laps_cycle_time
!    1                                          -laps_cycle_time

               y = bottom - .015 * r_height

               call label_other_stations(i4time_label,standard_longitude       
     1                                  ,y,xsta,lat,lon,nx_l,ny_l
     1                                  ,grid_spacing_m
     1                                  ,xlow,xhigh,ylow,yhigh,nx_c
     1                                  ,bottom,r_height,maxstns)

!           splice end for labeling other stations

            call frame
            i_map = 0
            i_graphics_overlay = 0
            i_label_overlay = 0 ! i_label_overlay - 1
            n_image = 0
          endif ! c_field = 'fcf'
          goto 100
        endif

        if(c_field(1:2) .eq. 'df' .and. idiff .eq. 0)then
            if(ifield_found .eq. 0)then
                write(6,*)' skip difference plot - no field was found'
                goto100
            endif

            write(6,*)' plotting difference field of last two entries'       
            call diff(field_vert,field_vert_buf,field_vert ! _diff
     1               ,nx_c,nz_c)       

            c_label = 'difference field (b - a)'

            scale = 1.
            cint = 0.

            if(dyn_low  .eq. r_missing_data .or. 
     1         dyn_high .eq. r_missing_data)then ! initialize range
                call contour_settings(field_vert,nx_c,nz_c
     1                   ,clow,chigh,cint,zoom,density,scale)      
                dyn_low  = clow  ! save for subsequent times
                dyn_high = chigh ! save for subsequent times
            else ! use previously initialized values
                clow  = dyn_low
                chigh = dyn_high
            endif

            i_contour = 1
            idiff = 1

        else
            if(igrid .eq. 1)then
                write(6,*)' copying field_vert to field_buf'
                call move(field_vert,field_vert_buf,nx_c,nz_c)       
            endif
            idiff = 0
        endif

        igrid = 1

        if(    c_field .eq. 'di' .or. c_field      .eq. 'sp'
     1    .or. c_field .eq. 'u ' .or. c_field      .eq. 'v '
     1    .or. c_field .eq. 'w ' .or. c_field      .eq. 'dv'
     1    .or. c_field .eq. 'vo' .or. c_field      .eq. 'va' 
     1    .or. c_field .eq. 'pv'
     1    .or. c_field .eq. 'vc' .or. c_field(1:2) .eq. 'om' )then

            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            if(istatus .ne. 1)then
                write(6,*)' error: could not determine product time'
                goto100
            endif

            i4time_3dw = i4_valid

            if(c_field .ne. 'w ' .and. c_field(1:2) .ne. 'om')then
                if(c_prodtype .eq. 'n')then
                    c_wind = 'b'
                else
                    c_wind = 'k'
                endif

            elseif(c_field(1:2) .eq. 'om' .and. 
     1             (c_prodtype .eq. 'a' .or. c_prodtype .eq. 'n') )then       
                write(6,105)
105             format('  omega field: kinematic (lw3), '
     1                ,'cloud bogused'
     1                  ,' [k,c]  ',7x,'? ',$)

                read(lun,106)c_wind
106             format(a)
                call downcase(c_wind,c_wind)

            endif

            if  (c_prodtype .eq. 'n')then
                ext_wind = 'balance'
                call get_directory(ext_wind,directory,len_dir)
                c_filespec = directory(1:len_dir)//'lw3/*.lw3'
                ext = 'lw3'

            elseif(c_wind .eq. 'c')then
                ext_wind = 'lco'
                call get_directory(ext_wind,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext_wind(1:3)

            elseif(c_prodtype .eq. 'a')then
                ext_wind = 'lw3'
                call get_directory(ext_wind,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext_wind(1:3)

            else ! background or forecast
                ext_wind = ext
!               call get_directory(ext_wind,directory,len_dir)
                call s_len(directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext_wind(1:3)
                call directory_to_cmodel(directory,c_model)

            endif

!           if(.not. l_wind_read)then
            if(.true.)then
                write(6,*)
                write(6,*)' looking for 3d wind data: ',ext_wind(1:10)
     1                   ,' ',trim(ext),trim(c_field),trim(c_wind)

                if(c_field .ne. 'w ' .and. c_field(1:2) .ne. 'om')then ! non-omega
                    if(c_prodtype .eq. 'n')then
                        write(6,*)' reading u/v via get_3dgrid_dname'
                        directory = directory(1:len_dir)//'lw3'
                        ext = 'lw3'

                        var_2d = 'u3'

                        call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_3dw
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,nx_l,ny_l,nz_l,u_3d,istatus)       

                        var_2d = 'v3'

                        call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_3dw       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,nx_l,ny_l,nz_l,v_3d,istatus)       

                    elseif(c_prodtype .eq. 'a')then
                        write(6,*)' reading u/v via get_uv_3d'
                        call get_file_time(c_filespec,i4time_ref
     1                                               ,i4time_3dw)
                        call get_uv_3d(i4time_3dw,nx_l,ny_l,nz_l
     1                                  ,u_3d,v_3d,ext_wind,istatus)

                    else ! background or forecast
                        write(6,*)' reading u/v via get_lapsdata_3d'
                        var_2d = 'u3'
                        call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,u_3d
     1                              ,istatus)
                        write(6,*)' istatus from u3 read is ',istatus

                        var_2d = 'v3'
                        call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,v_3d
     1                              ,istatus)
                        write(6,*)' istatus from v3 read is ',istatus

                    endif

                    call make_fnam_lp(i4time_3dw,a9time,istat_a9)
                    write(6,*)' a9time = ',a9time

                elseif(c_field .eq. 'w ' .or. 
     1                 c_field(1:2) .eq. 'om')then ! omega
                    write(6,*)' reading omega/w ',c_wind,c_prodtype

                    if(c_prodtype .eq. 'b' .or. 
     1                 c_prodtype .eq. 'f')then ! c_prodtype
                        var_2d = 'om'
                        write(6,*)'call get_lapsdata_3d',trim(directory)      
                        call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istatus)

                        call make_fnam_lp(i4_valid,a9time,istat_a9)

                    else 
                        call get_file_time(c_filespec
     1                                    ,i4time_ref,i4time_3dw)

                        if(c_wind .eq. 'c')then
                            write(6,*)'call get_w_3d ',trim(ext_wind)
                            call get_w_3d(i4time_3dw,nx_l,ny_l,nz_l
     1                                      ,field_3d,ext_wind,istatus)  

                        elseif(c_prodtype .eq. 'a')then
                            write(6,*)'call get_w_3d ',trim(ext_wind)
                            call get_w_3d(i4time_3dw,nx_l,ny_l,nz_l
     1                                      ,field_3d,ext_wind,istatus)

                        elseif(c_prodtype .eq. 'n')then
                            directory = directory(1:len_dir)//'lw3'
                            ext = 'lw3'
                            var_2d = 'om'

                            write(6,*)'call get_3dgrid_dname ',trim(ext)

                            call get_3dgrid_dname(directory
     1                      ,i4time_ref,laps_cycle_time*10000,i4time_3dw       
     1                      ,ext,var_2d,units_2d
     1                      ,comment_2d,nx_l,ny_l,nz_l,field_3d,istatus)       
                        endif

                        call make_fnam_lp(i4time_3dw,a9time,istat_a9)

                    endif ! c_prodtype
                endif
            endif
!           l_wind_read = .true.

            if(istatus .ne. 1)then
                write(6,*)' error reading in wind field'
                goto100
            else
                write(6,*)' istatus from wind read is ',istatus
            endif

        endif

        if(c_field .eq. 'vc')then

!           remap from 3d grid to vert xsect grid
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 nx_l,ny_l,nx_c,r_missing_data)

            i_contour = 2

            if       (c_prodtype .eq. 'n')then
                c_label = 'wind (balanced)             (kt) '
            elseif   (c_prodtype .eq. 'a')then
                c_label = 'wind (analyzed)             (kt) '
            elseif(c_prodtype .eq. 'b' .or. c_prodtype .eq. 'f')then   
                comment_2d = 'wind'
                units_2d = 'kt'

                len_label = len(c_label)                          ! debug
                write(6,*)' label length = ',len_label            ! debug
 
                call mk_fcst_xlabel(comment_2d,fcst_hhmm
     1                             ,ext(1:3),units_2d,c_model,c_label)       
            else
                c_label = 'laps wind (??????????)     knots '
            endif

        elseif(c_field(1:2) .eq. 'om' )then

            if(ext_wind .eq. 'lco')then ! cloud omega

!               take out missing data values to ensure better interpolation
                do k = 1,nz_l
                do j = 1,ny_l
                do i = 1,nx_l
                    if(field_3d(i,j,k) .eq. r_missing_data)then
                        field_3d(i,j,k) = -1e-30
                    endif
                enddo ! i
                enddo ! j
                enddo ! k

                call interp_3d(field_3d,field_vert,xlow,xhigh,ylow
     1                        ,yhigh,nx_l,ny_l,nz_l,nx_c,nz_c
     1                        ,r_missing_data)

                call get_pres_3d(i4time_ref,nx_l,ny_l,nz_l,pres_3d
     1                                                    ,istatus)
                call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow
     1                            ,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

                do k = 1,nz_c
                do i = 1,nx_c
                    if(field_vert(i,k) .ne. r_missing_data)then
!                       field_vert(i,k) = field_vert(i,k)*10. ! pa/s to ubar/s
                        field_vert(i,k) = field_vert(i,k) * 100.! .1 ubar/s
                    else
                        field_vert(i,k) = -1e-30
                    endif
                enddo ! i
                enddo ! k

            else ! not lco field
                call interp_3d(field_3d,field_vert
     1                        ,xlow,xhigh,ylow,yhigh
     1                        ,nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)       


                call get_pres_3d(i4time_ref,nx_l,ny_l,nz_l,pres_3d
     1                                                    ,istatus)
                call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,       
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)


                do k = nz_c,1,-1
                do i = 1,nx_c
                    if(field_vert(i,k) .ne. r_missing_data)then
!                       field_vert(i,k) = field_vert(i,k) * 10. ! pa/s to ubar/s
                        field_vert(i,k) = field_vert(i,k) * 100.! .1 ubar/s
!                   else
!                       field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                    endif
                enddo ! i
                enddo ! k

            endif ! lco field

            colortable = 'omega'

            if(i_image .eq. 0)then ! contour plot
                cint = -1. * 2. ** (-density)
                i_contour = 1
            else                   ! image plot
!               cint = 0.
                cint = 10. 

                if(grid_spacing_m .ge. 5500.)then
                    chigh = 80.
                    clow = -80.
                elseif(grid_spacing_m .ge. 3500.)then
                    chigh = 800.
                    clow = -800.
                elseif(grid_spacing_m .ge. 2000.)then
                    chigh = 800.
                    clow = -800.
                else
                    chigh = 800.
                    clow = -800.
                endif
                i_contour = 1

            endif

            write(6,*)' omega i_image/i_contour is ',i_image,i_contour

            if       (c_prodtype .eq. 'n')then
                c_label = 'laps omega (balanced)  0.1ubar/s'
            else   if(c_wind .eq. 'c')then
                c_label = 'laps omega (cloud)     0.1ubar/s'
            else   if(c_prodtype .eq. 'a')then
                c_label = 'laps omega (analyzed)  0.1ubar/s'
            else   if(c_prodtype .eq. 'b')then
                c_label = 'laps  bkgnd   omega  '//fcst_hhmm
     1                                          //'  0.1ubar/s'
            else   if(c_prodtype .eq. 'f')then
!               c_label = 'laps  fua     omega  '//fcst_hhmm
!    1                                             //'  ubar/s'
!               c_label = 'laps  '//c_model//' omega  '//fcst_hhmm
!    1                                             //'  ubar/s'
                call mk_fcst_xlabel('omega',fcst_hhmm
     1                       ,ext(1:3),'0.1ubar/s',c_model,c_label)
            else
                c_label = '                                 '
            endif

        elseif(c_field .eq. 'w ' )then

            i_contour = 1

            if(ext_wind .eq. 'lco')then ! cloud omega

!               take out missing data values to insure better interpolation
                do k = 1,nz_l
                do j = 1,ny_l
                do i = 1,nx_l
                    if(field_3d(i,j,k) .eq. r_missing_data)then
                        field_3d(i,j,k) = -1e-30
                    endif
                enddo ! i
                enddo ! j
                enddo ! k

                call interp_3d(field_3d,field_vert,xlow,xhigh,ylow
     1                        ,yhigh,nx_l,ny_l,nz_l,nx_c,nz_c
     1                        ,r_missing_data)

                call get_pres_3d(i4time_ref,nx_l,ny_l,nz_l,pres_3d
     1                                                    ,istatus)
                call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow
     1                            ,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

                do k = 1,nz_c
                do i = 1,nx_c
                    if(field_vert(i,k) .ne. r_missing_data)then
                        if(l_convert)then
                            field_vert(i,k) = ! always should be .true.
     1           omega_to_w(field_vert(i,k),pres_vert(i,k)) * 100.
                        endif
                    else
                        field_vert(i,k) = -1e-30
                    endif
                enddo ! i
                enddo ! k

                cint = -1. * 2. ** (-density)

            else ! not lco field
                call interp_3d(field_3d,field_vert
     1                        ,xlow,xhigh,ylow,yhigh
     1                        ,nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)       


                call get_pres_3d(i4time_ref,nx_l,ny_l,nz_l,pres_3d
     1                                                    ,istatus)
                call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,       
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)


                do k = nz_c,1,-1
                do i = 1,nx_c
                    if(field_vert(i,k) .ne. r_missing_data)then
!                       l_convert is .false. if old 'w' data is read in (not om)
                        if(l_convert)field_vert(i,k) =
     1             omega_to_w(field_vert(i,k),pres_vert(i,k)) * 100.       
                    else
                        field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                    endif
                enddo ! i
                enddo ! k

                cint = -1. * 2. ** (-density)

            endif ! lco field

            if       (c_prodtype .eq. 'n')then
                c_label = 'laps w (bal)   vert x-sect (cm/s)'
            else   if(c_wind .eq. 'c')then
                c_label = 'laps w (cloud) vert x-sect (cm/s)'
            else ! if(c_prodtype .eq. 'a')then
                c_label = 'laps w (kinem) vert x-sect (cm/s)'
            endif


        elseif(c_field .eq. 'u ' )then
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 nx_l,ny_l,nx_c,r_missing_data)
            do k = nz_c,1,-1
            do i = 1,nx_c
                if(u_vert(i,k) .ne. r_missing_data)then
                    field_vert(i,k) = u_vert(i,k)/mspkt
                else
                    field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                endif
            enddo ! i
            enddo ! k

            clow =  -chigh_3dwind
            chigh = +chigh_3dwind
            if(i_image .eq. 0)then
                cint = 10. / density
            else
                colortable = 'spectral'
                cint = 20.
            endif

            i_contour = 1

            if    (c_prodtype .eq. 'n')then
                c_label = 'u component (balanced)      (kt) '
            elseif(c_prodtype .eq. 'a')then   
                c_label = 'u component (analyzed)      (kt) '
            elseif(c_prodtype .eq. 'b' .or. c_prodtype .eq. 'f')then   
                comment_2d = 'u component'
                units_2d = 'kt'

                len_label = len(c_label)                          ! debug
                write(6,*)' label length = ',len_label            ! debug

                call mk_fcst_xlabel(comment_2d,fcst_hhmm
     1                             ,ext(1:3),units_2d,c_model,c_label)       
            endif

        elseif(c_field .eq. 'v ' )then
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 nx_l,ny_l,nx_c,r_missing_data)
            do k = nz_c,1,-1
            do i = 1,nx_c
                if(v_vert(i,k) .ne. r_missing_data)then
                    field_vert(i,k) = v_vert(i,k)/mspkt
                else
                    field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                endif
            enddo ! i
            enddo ! k

            clow =  -chigh_3dwind
            chigh = +chigh_3dwind
            if(i_image .eq. 0)then
                cint = 10. / density
            else
                colortable = 'spectral'
                cint = 20.
            endif

            i_contour = 1

            if    (c_prodtype .eq. 'n')then
                c_label = 'v component (balanced)      (kt) '
            elseif(c_prodtype .eq. 'a')then   
                c_label = 'v component (analyzed)      (kt) '
            elseif(c_prodtype .eq. 'b' .or. c_prodtype .eq. 'f')then   
                comment_2d = 'v component'
                units_2d = 'kt'
                call mk_fcst_xlabel(comment_2d,fcst_hhmm
     1                             ,ext(1:3),units_2d,c_model,c_label)       
            endif

        elseif(c_field .eq. 'sp' )then
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 nx_l,ny_l,nx_c,r_missing_data)
            do k = nz_c,1,-1
            do i = 1,nx_c
                if(v_vert(i,k) .ne. r_missing_data)then
                    call uv_to_disp(u_vert(i,k),
     1                              v_vert(i,k),
     1                              di_dum,
     1                              speed_ms)
                    field_vert(i,k) = speed_ms/mspkt
                else
                    field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                endif
            enddo ! i
            enddo ! k

            clow =  0.
            chigh = chigh_3dwind
            if(i_image .eq. 0)then
                cint = 10. / density
            else
                colortable = 'spectral'
                cint = 20.
            endif

            i_contour = 1

            if       (c_prodtype .eq. 'n')then
                c_label = 'laps isotachs (balanced)   knots '
            elseif   (c_prodtype .eq. 'a')then
                c_label = 'laps isotachs (analyzed)   knots '
            elseif   (c_prodtype .eq. 'b')then
                c_label = 'laps isotachs (background) knots '
            elseif   (c_prodtype .eq. 'f')then
                c_label = 'laps isotachs (forecast)   knots '
            endif

        elseif(c_field .eq. 'di' )then
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 nx_l,ny_l,nx_c,r_missing_data)
            do k = nz_c,1,-1
            do i = 1,nx_c
                if(u_vert(i,k) .ne. r_missing_data)then
                    call uv_to_disp(u_vert(i,k),
     1                              v_vert(i,k),
     1                              field_vert(i,k),
     1                              speed_dum)
                else
                    field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                endif
            enddo ! i
            enddo ! k
            clow = -100.
            chigh = +1000.
            cint = 10. / density
            i_contour = 1
            if       (c_prodtype .eq. 'n')then
                c_label = 'laps isogons (balanced)    knots '
            elseif   (c_prodtype .eq. 'a')then
                c_label = 'laps isogons (analysis)    knots '
            elseif   (c_prodtype .eq. 'b')then
                c_label = 'laps isogons (background)  knots '
            elseif   (c_prodtype .eq. 'f')then
                c_label = 'laps isogons (forecast)    knots '
            endif

        elseif(c_field .eq. 'dv' )then
            do k = 1,nz_l
                call divergence(u_3d(1,1,k),v_3d(1,1,k),field_3d(1,1,k)       
     1                         ,lat,lon,nx_l,ny_l,.true.,r_missing_data)       
            enddo

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            do k = nz_c,1,-1
              do i = 1,nx_c
                if(field_vert(i,k) .eq. r_missing_data)then
                    field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                else
                    field_vert(i,k) = field_vert(i,k) * 1e5
                endif
              enddo ! i
            enddo ! k

            if(i_image .eq. 0)then
!               clow = -100.
!               chigh = +1000.
!               cint = 2. / density
                cint = -1. * 2. ** (-density)
                i_contour = 1
                i_contour = 1
            else
                cint = 10. 
                clow = -40.
                chigh = +40.
                i_contour = -1
            endif

            if       (c_prodtype .eq. 'n')then
                c_label = 'laps divergence  (bal)   [1e-5/s]'
            elseif   (c_prodtype .eq. 'a')then
                c_label = 'laps divergence  (anal)  [1e-5/s]'
            elseif   (c_prodtype .eq. 'b')then
                c_label = 'laps divergence  (bkgnd) [1e-5/s]'
            elseif   (c_prodtype .eq. 'f')then
!               c_label = 'laps divergence  (fcst)  [1e-5/s]'
                call mk_fcst_xlabel('divergence',fcst_hhmm
     1                             ,ext(1:3),'1e-5/s',c_model,c_label)       
            endif

            plot_parms%iraster = 1

        elseif(c_field .eq. 'vo' )then
            do k = 1,nz_l
                call vorticity_abs(u_3d(1,1,k),v_3d(1,1,k)
     1                            ,field_3d(1,1,k),lat,lon
     1                            ,nx_l,ny_l,dx,dy
     1                            ,.true.,r_missing_data)       
            enddo

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            do k = nz_c,1,-1
              do i = 1,nx_c
                if(field_vert(i,k) .eq. r_missing_data)then
                    field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                else
                    field_vert(i,k) = field_vert(i,k) * 1e5
                endif
              enddo ! i
            enddo ! k

            clow = -100.
            chigh = +100.
            cint = 2. / density

            i_contour = 1

            if       (c_prodtype .eq. 'n')then
                c_label = 'laps abs vort  (bal)     [1e-5/s]'
            elseif   (c_prodtype .eq. 'a')then
                c_label = 'laps abs vort  (anal)    [1e-5/s]'
            elseif   (c_prodtype .eq. 'b')then
                c_label = 'laps abs vort  (bkgnd)   [1e-5/s]'
            elseif   (c_prodtype .eq. 'f')then
                c_label = 'laps abs vort  (fcst)    [1e-5/s]'
            endif

            plot_parms%iraster = 1

        elseif(c_field .eq. 'pv' )then
!           read 3d temperature
            write(6,*)' obtaining 3d temperature'
            if(c_prodtype .eq. 'a')then
                iflag_temp = 1 ! returns ambient temp (k)
                call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,nx_l,ny_l,nz_l,temp_3d,istatus)

            elseif(c_prodtype .eq. 'n')then
                iflag_temp = 4 ! returns balanced temp (k)
                call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,nx_l,ny_l,nz_l,temp_3d,istatus)

            elseif(c_prodtype .eq. 'b' .or. 
     1             c_prodtype .eq. 'f')then
                var_2d = 't3'
                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,temp_3d
     1                              ,istatus)
                if(istatus .ne. 1)goto100

            else
                write(6,*)' sorry, temperature source not yet supported'       
                goto100

            endif

            call calc_potvort(i4time_ref,u_3d,v_3d,temp_3d,field_3d
     1                       ,lat,lon,nx_l,ny_l,nz_l,nz_l,0,.true.       
     1                       ,dx,dy,r_missing_data,istatus)       

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            do k = nz_c,1,-1
              do i = 1,nx_c
                if(field_vert(i,k) .eq. r_missing_data)then
                    field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                else
                    field_vert(i,k) = field_vert(i,k) * 1e6
                endif
              enddo ! i
            enddo ! k

            clow = -10.
            chigh = +10.
            cint = 2. / density

            i_contour = 1

            if       (c_prodtype .eq. 'n')then
                c_label = 'laps pvort  (bal)    pvu '
            elseif   (c_prodtype .eq. 'a')then
                c_label = 'laps pvort  (anal)   pvu '
            elseif   (c_prodtype .eq. 'b')then
                c_label = 'laps pvort  (bkgnd)  pvu '
            elseif   (c_prodtype .eq. 'f')then
                c_label = 'laps pvort  (fcst)   pvu '
            endif

            plot_parms%iraster = 1

        elseif(c_field .eq. 'va' )then
            do k = 1,nz_l
                call vorticity_abs(u_3d(1,1,k),v_3d(1,1,k)
     1                            ,field_2d,lat,lon
     1                            ,nx_l,ny_l
     1                            ,dx,dy,.true.,r_missing_data)       

                call cpt_advection(field_2d,u_3d(1,1,k),v_3d(1,1,k)
     1                            ,dx,dy,nx_l,ny_l
     1                            ,field_3d(1,1,k))

            enddo

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            do k = nz_c,1,-1
              do i = 1,nx_c
                if(field_vert(i,k) .eq. r_missing_data)then
                    field_vert(i,k) = field_vert(i,min(k+1,nz_c))
                else
                    field_vert(i,k) = field_vert(i,k) * 1e8
                endif
              enddo ! i
            enddo ! k

            clow = -100.
            chigh = +1000.
            cint = 5. / density

            i_contour = 1

            if       (c_prodtype .eq. 'n')then
                c_label = 'laps vort adv  (bal)   [1e-8/s^2]'
            elseif   (c_prodtype .eq. 'a')then
                c_label = 'laps vort adv  (anal)  [1e-8/s^2]'
            elseif   (c_prodtype .eq. 'b')then
                c_label = 'laps vort adv  (bkgnd) [1e-8/s^2]'
            elseif   (c_prodtype .eq. 'f')then
                c_label = 'laps vort adv  (fcst)  [1e-8/s^2]'
            endif

            plot_parms%iraster = 1

        elseif(c_field .eq. 'rf' .or. c_field .eq. 'rg'
     1                           .or. c_field .eq. 'rk')then
            if(c_field .ne. 'rg')then
                i4time_get = i4time_ref/laps_cycle_time 
     1                     * laps_cycle_time
                goto1300
            endif

1300        write(6,*)' getting radar data via get_laps_3dgrid'
!           var_2d = 'ref'
!           ext = 'lps'

!           call get_laps_3dgrid(i4time_ref,86400,i4time_radar,
!    1          nx_l,ny_l,nz_l,ext,var_2d
!    1                  ,units_2d,comment_2d,grid_ra_ref,istatus)

            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            var_2d = 'ref'

            if(c_prodtype .eq. 'a')then
                write(6,*)' getting radar data via get_laps_3dgrid'
                ext = 'lps'

                call get_laps_3dgrid(i4time_get,86400,i4time_radar
     1              ,nx_l,ny_l,nz_l,ext,var_2d
     1              ,units_2d,comment_2d,grid_ra_ref,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' could not read lps via get_laps_3dgrid'       
                    goto100
                endif

                call make_fnam_lp(i4time_radar,a9time,istatus)
                c_label = 'laps  reflectivity  vert x-sect  '

            elseif(c_prodtype .eq. 'f')then
                call get_lapsdata_3d(i4_initial,i4_valid,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,grid_ra_ref
     1                              ,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' could not read forecast ref'       
                    goto100
                endif
                c_label = 'laps  fua reflectivity '//fcst_hhmm
     1                    //'   dbz'

            else
                goto100

            endif

            call interp_3d(grid_ra_ref,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)
            clow = 0.
            chigh = +100.
            cint = 10. / density
            i_contour = 1
!           c_label = 'laps  reflectivity  vert x-sect  '
!           call make_fnam_lp(i4time_radar,a9time,istatus)

        elseif(c_field .eq. 'ri' .or. c_field .eq. 'rj' 
     1    .or. c_field .eq. 'rs' .or. c_field .eq. 'rv')then ! reflectivity image
            i_image = 1
            if(c_field .eq. 'ri')then
                i4time_get = i4time_ref/laps_cycle_time 
     1                                * laps_cycle_time
            else
                i4time_get = i4time_ref
            endif

            if(c_field .ne. 'ri')then
                if(c_field .eq. 'rv')then
!                   ask which radar number (extension)
                    write(6,*)
                    write(6,2026)
2026                format('  enter radar extension (for reflectivity)'      
     1                                                     ,28x,'? ',$)       
                    read(lun,*)ext_radar

                    write(6,*)' reading reflectivity data from radar '
     1                                          ,ext_radar

                    write(6,*)

                    call get_ref_base(ref_base,istatus)

                    call get_filespec(ext_radar,2,c_filespec,istatus)
                    call get_file_time(c_filespec,i4time_get
     1                                ,i4time_radar)    

                    if(ext_radar .ne. 'vrz')then ! vxx

                        write(6,*)' read height field for rfill purpose'

                        l_low_fill = namelist_parms%l_low_fill
                        l_high_fill = namelist_parms%l_high_fill

                        if(l_low_fill .or. l_high_fill)then
!                           obtain height field
                            ext = 'lt1'
                            var_2d = 'ht'
                            call get_laps_3dgrid(
     1                         i4time_radar,10000000,i4time_ht,
     1                         nx_l,ny_l,nz_l,ext,var_2d,
     1                         units_2d,comment_2d,field_3d,istatus)
                            if(istatus .ne. 1)then
                                write(6,*)' error locating height field'
                                call get_pres_3d(i4time_radar
     1                                ,nx_l,ny_l,nz_l,pres_3d,istatus)
                                if(istatus .ne. 1)then
                                    write(6,*)
     1                                 ' error getting pressure field'      
                                    goto 100
                                else
                                    write(6,*)' convert pres to ht'
     1                              ,' - using standard atmosphere'     
                                    do k = 1,nz_l
                                    do j = 1,ny_l
                                    do i = 1,nx_l
                                        field_3d(i,j,k) = 
     1                                  psatoz(pres_3d(i,j,k)*.01) 
                                    enddo ! i
                                    enddo ! j
                                    enddo ! k
                                endif
                            endif

                        endif ! height field may be necessary

                        i4_tol = 1200

                        call read_radar_3dref(i4time_radar,
     1                   i4_tol,i4_ret,                                  ! i/o
!    1                   .true.,ref_base,
     1                   .true.,r_missing_data,
     1                   nx_l,ny_l,nz_l,ext_radar,
     1                   lat,lon,topo,l_low_fill,l_high_fill,
     1                   field_3d,         
     1                   grid_ra_ref,
     1                   rlat_radar,rlon_radar,rheight_radar,radar_name,       
     1                   n_ref_grids,istat_radar_2dref,
     1                   istat_radar_3dref)      

                        call make_fnam_lp(i4time_radar,a9time,istatus)
                        call filter_string(radar_name)
                        c_label = 'reflectivity dbz          '
     1                            //ext_radar(1:3)//'    '
     1                            //radar_name(1:4)       

                    else
                        write(6,*)
     1                  ' getting radar vrz mosaic via get_laps_3dgrid'       
                        ext = 'vrz'
                        var_2d = 'ref'

                        call get_laps_3dgrid(i4time_get,86400
     1                    ,i4time_radar
     1                    ,nx_l,ny_l,nz_l,ext,var_2d
     1                    ,units_2d,comment_2d,grid_ra_ref,istatus)

                        if(istatus .ne. 1)then
                            write(6,*)
     1                      ' could not read lps via get_laps_3dgrid'       
                            goto100
                        endif

                        call make_fnam_lp(i4time_radar,a9time,istatus)
                        c_label = 'reflectivity dbz             '
     1                              //ext_radar(1:3)

                    endif

!                   i4time_radar = i4time_get

                elseif(.not. l_radar_read)then
!                   obtain height field
                    ext = 'lt1'
                    var_2d = 'ht'
                    call get_laps_3dgrid(
     1                   i4time_get,10000000,i4time_ht,
     1                   nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,heights_3d,istatus)
                    if(istatus .ne. 1)then
                        write(6,*)' error locating height field'
                        go to 100
                    endif

                    call get_radar_ref(i4time_get,2000,i4time_radar,1
     1                ,1,nx_l,ny_l,nz_l,lat,lon,topo,.true.,.true.
     1                ,heights_3d
     1                ,grid_ra_ref,n_ref
     1                ,rlat_radar,rlon_radar,rheight_radar,istat_2dref
     1                ,istat_3dref)

                    if(istat_2dref .le. 0)goto 100

                    l_radar_read = .true.

                    if(istat_3dref .le. 0)then
                        if(istat_2dref .eq. 1)then
                            write(6,*)
     1                      ' radar xsect unavailable, try earlier time'
                        endif
                        goto 100
                    endif

                endif

            else ! 'ri'
                call input_product_info(i4time_ref              ! i
     1                                 ,laps_cycle_time         ! i
     1                                 ,3                       ! i
     1                                 ,c_prodtype              ! o
     1                                 ,ext                     ! o
     1                                 ,directory               ! o
     1                                 ,a9time                  ! o
     1                                 ,fcst_hhmm               ! o
     1                                 ,i4_initial              ! o
     1                                 ,i4_valid                ! o
     1                                 ,istatus)                ! o

                var_2d = 'ref'

                if(c_prodtype .eq. 'a')then
                    write(6,*)' getting radar data via get_laps_3dgrid'
                    ext = 'lps'

                    call get_laps_3dgrid(i4time_get,86400,i4time_radar
     1                  ,nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,grid_ra_ref,istatus)
                    if(istatus .ne. 1)then
                        write(6,*)
     1                  ' could not read lps via get_laps_3dgrid'       
                        goto100
                    endif

                    call make_fnam_lp(i4time_radar,a9time,istatus)
                    c_label = 'laps  reflectivity  vert x-sect  '

                elseif(c_prodtype .eq. 'f')then
                    call get_lapsdata_3d(i4_initial,i4_valid
     1                                  ,nx_l,ny_l,nz_l       
     1                                  ,directory,var_2d
     1                                  ,units_2d,comment_2d,grid_ra_ref       
     1                                  ,istatus)
                    if(istatus .ne. 1)then
                        write(6,*)' could not read forecast ref'       
                        goto100
                    endif
                    c_label = 'laps  fua reflectivity '//fcst_hhmm
     1                        //'   dbz'

                else
                    goto100

                endif

            endif

            call get_ref_base(ref_base,istatus)

            do i=1,nx_l
            do j=1,ny_l
            do k=1,nz_l
                if(grid_ra_ref(i,j,k) .eq. r_missing_data)then
                    grid_ra_ref(i,j,k) = ref_base
                endif  
            enddo ! k
            enddo ! j
            enddo ! i

            if(c_field .ne. 'rs')then
                call interp_3d(grid_ra_ref,field_vert
     1                        ,xlow,xhigh,ylow,yhigh
     1                        ,nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)       

            else ! get spread out data
                call interp_3d_spread
     1          (grid_ra_ref,field_vert,xlow,xhigh,ylow,yhigh,
     1           nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            endif

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

            write(6,*)' generating reflectivity image'
            do i = 1,nx_c-1
            do k = ibottom,nz_l-1
                x1 = i
                y1 = k

!               perform a bi-linear interpolation to provide an image
!               this image consists of a set of line segments which will look
!               smoother than blocks the size of the grid.

                z1=field_vert(i  , k  )
                z2=field_vert(i+1, k  )
                z3=field_vert(i+1, k+1)
                z4=field_vert(i  , k+1)

                nii = 15
                njj = 60

                do ii = 0,nii
                  fraci = float(ii) / float(nii)
                  xx = x1 + fraci

                  do jj = 0,njj
                    fracj = float(jj) / float(njj)
                    yy = y1 + fracj
                    r_dbz =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                  - (z2+z4-z3-z1)*fraci*fracj

                    i_dbz = nint(r_dbz)
                    i_dbz = i_dbz/5 * 5 ! this reduces the resolution and saves on graphics

!                   plot line segments of the same color
                    if(i_dbz .ne. i_dbz_ref
     1                        .and. fracj .gt. 0.)then

                        y_last = yy - 1./float(njj)

                        if(i_dbz_ref .ge. 0)then
                            icol = 180 + i_dbz_ref / 5

!                           write(6,432)i_dbz,i_dbz_ref,xx,y_ref,y_last
 432                        format(2i5,3f9.4)

                            call setusv_dum(2hin,icol)

                            call line(xx,y_ref,xx,y_last)
                        endif

                        i_dbz_ref = i_dbz
                        y_ref = yy

                    elseif(fracj .eq. 0.0)then
                        i_dbz_ref = i_dbz
                        y_ref = yy

                    endif

                    if(fracj .eq. 1.0)then
                        if(i_dbz .ge. 0)then
                           icol = 180 + i_dbz_ref / 5

!                          write(6,432)i_dbz,i_dbz_ref,xx,y_ref,y_last        

                           call setusv_dum(2hin,icol)

                           call line(xx,y_ref,xx,yy)

                        endif
                    endif

                  enddo
                enddo
            enddo ! k
            enddo ! i

            i_contour = 0

        elseif(c_field .eq. 'cf' )then ! cloud gridded image
            i_image = 1
            i_contour = -1

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

            call setusv_dum(2hin,2)

            write(6,*)' plotting cloud gridded image'

            ext = 'lcp'
            var_2d = 'lcp'
            call get_laps_3dgrid(
     1                   i4time_ref,10000000,i4time_nearest,
     1                   nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,field_3d,istatus)
            if(istatus .ne. 1)then
                write(6,*)' no cloud grid available'
            endif

            call make_fnam_lp(i4time_nearest,a9time,istatus)

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            c_label = 'gridded cloud cover        x-sect'

            call remap_field_2d(
     1                            nx_c,1,nx_c
     1                           ,nz_c,ibottom,nz_c
     1                           ,nx_p, ixl_remap, nx_p-ixh_remap+1
     1                           ,nx_p, 1+iyl_remap, nx_p-iyh_remap
     1                           ,field_vert,field_vert3,r_missing_data)

!           blank out the edges external to the x-section
!           write(6,*)' blackening the edges'
!           do i = 1,nx_p
!           do j = 1,nx_p
!               if(field_vert3(i,j) .eq. r_missing_data)then
!                   field_vert3(i,j) = 0.
!               endif
!           enddo ! j
!           enddo ! i
!           where(field_vert3(:,:).eq.r_missing_data)field_vert3(:,:)=0.

            write(6,*)' calling solid fill cloud plot - commented out'       
            colortable = 'linear'
            clow = 0.0
            chigh = 1.0
            cint = 0.1
            scale = 1e0

        elseif(c_field .eq. 'cg' )then ! cloud gridded image
            i_image = 1

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

            call setusv_dum(2hin,2)

            write(6,*)' plotting cloud gridded image'

            ext = 'lc3'
            call get_clouds_3dgrid(i4time_ref,i4time_nearest
     1               ,nx_l,ny_l,kcloud
     1               ,ext,clouds_3d,cld_hts,cld_pres,istatus)

            if(istatus .ne. 1)then
                write(6,*)' no cloud grid available'
            endif

            call make_fnam_lp(i4time_nearest,a9time,istatus)

            call interp_3dc(clouds_3d,clouds_vert,xlow,xhigh,ylow,yhigh,
     1               nx_l,ny_l,nx_c,r_missing_data)

            niii = 12 ! horizontal resolution of cloud plot

            do i = 1,nx_c

              i_eighths_ref = 0
              k_ref = ibottom+1

              do k = ibottom+1,kcloud-1

                if(clouds_vert(i,k) .ne. r_missing_data)then

                    if(l_atms)then ! turn into a binary (bipolar) field
                        if(clouds_vert(i,k) .gt. 0.65)then
                            clouds_vert(i,k) = 1.0
                        else
                            clouds_vert(i,k) = 0.0
                        endif
                    endif

                    i_eighths = nint(clouds_vert(i,k)*8.)

                else
                    i_eighths = 0
                endif

c               if(clouds_vert(i,k) .gt. 0.01)
c       1       write(6,1100)i,k,nint(cld_hts(k)),clouds_vert(i,k)
c       1                                               ,i_eighths
1100            format(2i3,i6,f6.2,i3)

                if(i_eighths .ne. i_eighths_ref)then

!                 remap clouds using standard atmosphere
                  chigh = (cld_hts(k) + cld_hts(k-1))/2.
                  clow  = (cld_hts(k_ref) + cld_hts(k_ref-1))/2.

!                 remap clouds using ambient pressure in center of domain
                  phigh = (cld_pres(k) + cld_pres(k-1))/2.
                  plow  = (cld_pres(k_ref) + cld_pres(k_ref-1))/2.

c                 write(6,1101)i_eighths_ref,nint(clow),nint(chigh)
1101              format(1x,i2,2i6)

                  if(i_eighths_ref .ge. 1)then
                    do ii = -niii,+niii

                      fraci = float(ii) / float(2*niii)

                      if((fraci+0.5) .lt. i_eighths_ref/8.0)then

!                       remap clouds using standard atmosphere
!                       r_low = height_to_zcoord(clow,istatus)
!                       r_high = height_to_zcoord(chigh,istatus)

!                       remap clouds using ambient pressure in center of domain
                        r_low = zcoord_of_pressure(plow)
                        r_high = zcoord_of_pressure(phigh)

                        uu = i + fraci

                        if(r_low .lt. nz_l)then ! keep below the top of the domain
                          call line(uu, r_low, uu, r_high)
                        endif

                      endif ! this line is covered
                    enddo ! ii
                  endif ! finite cloud cover in this grid segment

                  i_eighths_ref = i_eighths
                  k_ref = k

                endif ! change in cloud cover

              enddo ! k

            enddo ! i

!           generate cloud key
            if(.not. l_atms)then
              do i_eighths = 1,8
                i = 7 + 5 * i_eighths

                do ii = -niii,+niii
                      x = i-2
                      y = bottom - r_height * .055
                      write(c3_string,2013)i_eighths
                      call pwrity (x, y, c3_string, 3, 0, 0, 0)
2013                  format(i1,'/8')

                      fraci = float(ii) / float(2*niii)

                      if((fraci+0.5) .lt. i_eighths/8.0)then
                        r_low  = ibottom - 1. + 0.05
                        r_high = ibottom - 1. + 0.45
                        uu = i + fraci

                        if(r_low .lt. nz_l)then ! keep below the top of the domain
                          call line(uu, r_low, uu, r_high)
                        endif

                      endif ! this line is covered
                enddo ! ii
              enddo ! i

              c_label = 'laps gridded cloud cover   x-sect'

            endif ! l_atms

        elseif(c_field .eq. 'pt')then
            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            if(c_prodtype .eq. 'a')then ! original code
                iflag_temp = 0 ! returns potential temperature
                call get_temp_3d(i4time_ref,i4_valid,iflag_temp
     1                      ,nx_l,ny_l,nz_l,field_3d,istatus)
!               if(istatus .ne. 1)goto100
                c_label = 'laps potl temp vert x-sect    k  '

            elseif(c_prodtype .eq. 'b' .or. 
     1             c_prodtype .eq. 'f')then
                var_2d = 't3'
                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,temp_3d
     1                              ,istatus)
                if(istatus .ne. 1)goto100

!               convert from t to theta
                do i = 1,nx_l
                do j = 1,ny_l

                    do k = 1,nz_l
                        theta = o_k(temp_3d(i,j,k),zcoord_of_level(k))         
                        field_3d(i,j,k) = theta
                    enddo ! k

                enddo ! j
                enddo ! i

                if(c_prodtype .eq. 'b')then
                    c_label = 'laps  bkgnd   theta  '//fcst_hhmm
     1                                                 //'  deg k '
                elseif(c_prodtype .eq. 'f')then
                    call directory_to_cmodel(directory,c_model)
                    call mk_fcst_xlabel('theta',fcst_hhmm
     1                          ,ext(1:3),'deg k',c_model,c_label)       

                endif

            endif

            call make_fnam_lp(i4_valid,a9time,istatus)
            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            clow = 200.
            chigh = +500.
            cint = 5. / density
            i_contour = 1

        elseif(c_field .eq. 'te')then ! under construction
            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            if(c_prodtype .eq. 'a')then ! original code
                iflag_temp = 0 ! returns potential temperature
                call get_temp_3d(i4time_ref,i4_valid,iflag_temp
     1                      ,nx_l,ny_l,nz_l,field_3d,istatus)
!               if(istatus .ne. 1)goto100
                c_label = 'laps theta(e) vert x-sect    k  '

            elseif(c_prodtype .eq. 'b' .or. 
     1             c_prodtype .eq. 'f')then
                var_2d = 't3'
                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,temp_3d
     1                              ,istatus)
                if(istatus .ne. 1)goto100

!               convert from t to theta
                do i = 1,nx_l
                do j = 1,ny_l

                    do k = 1,nz_l
                        theta = o_k(temp_3d(i,j,k),zcoord_of_level(k))         
                        field_3d(i,j,k) = theta
                    enddo ! k

                enddo ! j
                enddo ! i

                if(c_prodtype .eq. 'b')then
                    c_label = 'laps  bkgnd  theta(e) '//fcst_hhmm
     1                                                //'  deg k '
                elseif(c_prodtype .eq. 'f')then
                    call directory_to_cmodel(directory,c_model)
                    call mk_fcst_xlabel('theta(e)',fcst_hhmm
     1                          ,ext(1:3),'deg k',c_model,c_label)       

                endif

            endif

            call make_fnam_lp(i4_valid,a9time,istatus)
            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            clow = 200.
            chigh = +500.
            cint = 5. / density
            i_contour = 1

        elseif(c_field .eq. 'pb')then
            iflag_temp = 3 ! returns balanced potential temperature
            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,nx_l,ny_l,nz_l,temp_3d,istatus)
!           if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            clow = 200.
            chigh = +500.
            cint = 5. / density
            i_contour = 1
            c_label = 'laps potl temp (balanced)     k  '

        elseif(c_field .eq. 'ht')then ! height field
            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            if(c_prodtype .eq. 'a')then
!               obtain height field
                ext = 'lt1'
                var_2d = 'ht'
                call get_laps_3dgrid(
     1                   i4time_ref,10000000,i4time_ht,
     1                   nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,field_3d,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' error locating height field'
                    go to 100
                endif

                c_label = 'laps height    vert x-sect  dm'

                call make_fnam_lp(i4time_ht,a9time,istatus)

            elseif(c_prodtype .eq. 'b' .or. 
     1             c_prodtype .eq. 'f')then
                var_2d = 'ht'
                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d       
     1                              ,istatus)
                if(istatus .ne. 1)goto100

                if(c_prodtype .eq. 'b')then
                    c_label = 'laps  bkgnd   ht     '//fcst_hhmm
     1                                                 //'  dm'
                elseif(c_prodtype .eq. 'f')then
                    c_label = 'laps  fua     ht     '//fcst_hhmm
     1                                                 //'  dm'
                endif

            else
                write(6,*)' sorry, not yet supported'
                goto100

            endif

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            clow = -100.
            chigh = +2000.
            cint = 100. / density
            i_contour = 1
            scale = 10.

        elseif(c_field .eq. 't ')then
            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            if(c_prodtype .eq. 'a')then
                iflag_temp = 1 ! returns ambient temp (k)

                call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,nx_l,ny_l,nz_l,temp_3d,istatus)

                c_label = 'laps temp      vert x-sect  deg c'

                call make_fnam_lp(i4time_nearest,a9time,istatus)

            elseif(c_prodtype .eq. 'b' .or. 
     1             c_prodtype .eq. 'f')then
                var_2d = 't3'
                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,temp_3d
     1                              ,istatus)
                if(istatus .ne. 1)goto100

                if(c_prodtype .eq. 'b')then
                    c_label = 'laps  bkgnd   t      '//fcst_hhmm
     1                                                 //'  deg c '
                elseif(c_prodtype .eq. 'f')then
!                   c_label = 'laps  fua     t      '//fcst_hhmm
!                                                      //'  deg c '

                    call directory_to_cmodel(directory,c_model)
                    call mk_fcst_xlabel('temperature',fcst_hhmm
     1                          ,ext(1:3),'deg c',c_model,c_label)       

                endif

            else
                write(6,*)' sorry, not yet supported'
                goto100

            endif

            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            do k = 1,nz_l
            do i = 1,nx_c
                field_vert(i,k) = k_to_c(field_vert(i,k))
            enddo ! i
            enddo ! k

            clow = -100.
            chigh = +100.
            cint = 5. / density
            i_contour = 1

        elseif(c_field .eq. 'tb')then
            var_2d = 't3'
            call make_fnam_lp(i4time_ref,a9time,istatus)
            ext='lt1'

            call get_directory('balance',directory,lend)
            directory=directory(1:lend)//'lt1/'

            call get_3dgrid_dname(directory
     1           ,i4time_ref,laps_cycle_time*10000,i4time_nearest
     1           ,ext,var_2d,units_2d
     1           ,comment_2d,nx_l,ny_l,nz_l,temp_3d,istatus)       

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            do k = 1,nz_l
            do i = 1,nx_c
                field_vert(i,k) = k_to_c(field_vert(i,k))
            enddo ! i
            enddo ! k

            clow = -100.
            chigh = +100.
            cint = 5. / density
            i_contour = 1
            c_label = 'laps temp (balanced)        deg c'

        elseif(c_field .eq. 'hb')then
            var_2d = 'ht'
            call make_fnam_lp(i4time_ref,a9time,istatus)
            ext='lt1'

            call get_directory('balance',directory,lend)
            directory=directory(1:lend)//'lt1/'

            call get_3dgrid_dname(directory
     1           ,i4time_ref,laps_cycle_time*10000,i4time_nearest
     1           ,ext,var_2d,units_2d
     1           ,comment_2d,nx_l,ny_l,nz_l,temp_3d,istatus)       

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            clow = -100.
            chigh = +2000.
            cint = 100. / density
            i_contour = 1
            scale = 10.

            c_label = 'laps ht (balanced)     (dm)'

        elseif(c_field .eq. 'ts')then
            iflag_temp = 1 ! returns ambient temp (k)
            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,nx_l,ny_l,nz_l,temp_3d,istatus)
!           if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)


            call get_pres_3d(i4time_nearest,nx_l,ny_l,nz_l,pres_3d
     1                                                    ,istatus)
            call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)


            do k = 1,nz_l
            do i = 1,nx_c
                field_vert(i,k) =
     1           os(field_vert(i,k)-273.15,pres_vert(i,k)/100.) + 273.15       
            enddo ! i
            enddo ! k

            clow = +250.
            chigh = +450.
            cint = 5. / density
            i_contour = 1
            c_label = 'laps theta(e) sat   x-sect  deg k'

        elseif(c_field .eq. 'be')then
            var_2d = 't3'
            call make_fnam_lp(i4time_ref,a9time,istatus)
            ext='lt1'

            call get_directory('balance',directory,lend)
            directory=directory(1:lend)//'lt1/'

            call get_3dgrid_dname(directory
     1           ,i4time_ref,laps_cycle_time*10000,i4time_nearest
     1           ,ext,var_2d,units_2d
     1           ,comment_2d,nx_l,ny_l,nz_l,temp_3d,istatus)       

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)


            call get_pres_3d(i4time_nearest,nx_l,ny_l,nz_l,pres_3d
     1                                                    ,istatus)
            call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)


            do k = 1,nz_l
            do i = 1,nx_c
                field_vert(i,k) =
     1           os(field_vert(i,k)-273.15,pres_vert(i,k)/100.) + 273.15       
            enddo ! i
            enddo ! k

            clow = +250.
            chigh = +450.
            cint = 5. / density
            i_contour = 1
            c_label = 'laps theta(e) sat (balanced) deg k'

        elseif(c_field .eq. 'tw' .or. c_field .eq. 'td')then
            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            iflag_temp = 1 ! returns ambient temp (k)
            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,nx_l,ny_l,nz_l,temp_3d,istatus)
!           if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            var_2d = 'rhl'
            ext = 'lh3'
            call get_laps_3dgrid
     1          (i4time_nearest,1000000,i4time_nearest,nx_l,ny_l,nz_l       
     1          ,ext,var_2d,units_2d,comment_2d,rh_3d,istatus)
            if(istatus .ne. 1)goto100

            call interp_3d(rh_3d,field_vert2,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            call get_pres_3d(i4time_nearest,nx_l,ny_l,nz_l,pres_3d
     1                                                    ,istatus)
            call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)


            do k = 1,nz_l
            do i = 1,nx_c
                t_c         = field_vert(i,k) - 273.15
                td_c        = dwpt(t_c,field_vert2(i,k))
                pressure_mb = pres_vert(i,k)/100.

!               this function call here is fast but returns a t_wb_c
!               equal to t_c if pres < 500mb. this approximation should
!               not hurt the algorithm.

                if(c_field .eq. 'tw')then
                    t_wb_c = twet_fast(t_c,td_c,pressure_mb)
                    field_vert(i,k) = t_wb_c
                else 
                    field_vert(i,k) = td_c
                endif
            enddo ! i
            enddo ! k

            clow = -100.
            chigh = +100.
            cint = 5. / density
            i_contour = 1
            if(c_field .eq. 'tw')then
                c_label = 'wet bulb temp       x-sect  deg c'
            else
                c_label = 'dew point           x-sect  deg c'
            endif

        elseif(c_field .eq. 'sh')then
            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            var_2d = 'sh '

            if(c_prodtype .eq. 'a')then
                ext = 'lq3'
                call get_laps_3dgrid
     1          (i4time_ref,1000000,i4time_nearest,nx_l,ny_l,nz_l
     1              ,ext,var_2d,units_2d,comment_2d
     1              ,q_3d,istatus)
                if(istatus .ne. 1)goto100

                call make_fnam_lp(i4time_nearest,a9time,istatus)

                c_label = 'laps specific humidity    (g/kg) '

            elseif(c_prodtype .eq. 'b' .or. 
     1             c_prodtype .eq. 'f')then
                var_2d = 'sh'
                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,q_3d
     1                              ,istatus)
                if(istatus .ne. 1)goto100

                call directory_to_cmodel(directory,c_model)
                call mk_fcst_xlabel('specific humidity',fcst_hhmm
     1                          ,ext(1:3),'g/kg',c_model,c_label)       

            elseif(c_prodtype .eq. 'n')then
                ext = 'lq3'
                call get_directory('balance',directory,len_dir)
                directory = directory(1:len_dir)//'lq3/'

                call get_3dgrid_dname(directory
     1             ,i4time_ref,1000000000,i4time_nearest
     1             ,ext,var_2d,units_2d
     1             ,comment_2d,nx_l,ny_l,nz_l,q_3d,istatus)       

                call make_fnam_lp(i4time_nearest,a9time,istatus)

                c_label = 'balanced specific humidity (g/kg) '

            endif

            call interp_3d(q_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            clow = 0.
            chigh = +25.
            cint = 1.0 / density
!           cint = -1.
            i_contour = 1
            scale = 1e-3

!           colortable = 'moist'
            colortable = 'tpw'

        elseif(c_field(1:2) .eq. 'rh' .or. c_field(1:2) .eq. 'rl')then
            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            if(c_prodtype .eq. 'a')then
                if(c_field(1:2) .eq. 'rh')then
                    var_2d = 'rh3'
                    c_label = 'laps relative humidity (anl)    %'
                elseif(c_field(1:2) .eq. 'rl')then
                    var_2d = 'rhl'
                    c_label = 'laps relative humidity (anl-liq) %'
                endif

                write(6,*)' reading lh3 / ',var_2d

                ext = 'lh3'
                call get_laps_3dgrid(i4time_ref,1000000,i4time_nearest
     1                              ,nx_l,ny_l,nz_l,ext,var_2d,units_2d
     1                              ,comment_2d,rh_3d,istatus)
                if(istatus .ne. 1)goto100

                call make_fnam_lp(i4time_nearest,a9time,istatus)

                call interp_3d(rh_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                         nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)      

            elseif(c_prodtype .eq. 'b' .or. 
     1             c_prodtype .eq. 'f')then
                var_2d = 'sh'
                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,q_3d
     1                              ,istatus)
                if(istatus .ne. 1)goto100

                var_2d = 't3'
                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,temp_3d
     1                              ,istatus)
                if(istatus .ne. 1)goto100

                call interp_3d(rh_3d,field_vert
     1                        ,xlow,xhigh,ylow,yhigh
     1                        ,nx_l,ny_l,nz_l,nx_c,nz_c
     1                        ,r_missing_data)      

                call directory_to_cmodel(directory,c_model)

                call get_pres_3d(i4time_nearest,nx_l,ny_l,nz_l,pres_3d
     1                          ,istatus)

                if(c_field(1:2) .eq. 'rh')then
                    t_ref = -10. 
                    call mk_fcst_xlabel('rh (tref = -10c)',fcst_hhmm
     1                              ,ext(1:3),'%',c_model,c_label)       
                elseif(c_field(1:2) .eq. 'rl')then
                    t_ref = -192.
                    call mk_fcst_xlabel('rh (liq)',fcst_hhmm
     1                              ,ext(1:3),'%',c_model,c_label)       
                endif

                write(6,*)' t_ref = ',t_ref

                do i = 1,nx_l
                do j = 1,ny_l
                do k = 1,nz_l
                    rh_3d(i,j,k)=make_rh(pres_3d(i,j,k)/100.
     1                                  ,k_to_c(temp_3d(i,j,k))
     1                                  ,q_3d(i,j,k)*1000.
     !                                  ,t_ref                )  * 100. 
                enddo ! k
                enddo ! j
                enddo ! i

            elseif(c_prodtype .eq. 'n')then
                call get_directory('balance',directory,len_dir)
                directory = directory(1:len_dir)//'lh3/'

                ext = 'lh3'

                if(c_field(1:2) .eq. 'rh')then
                    var_2d = 'rh3'
                    c_label = 'laps relative humidity (bal)    %'
                elseif(c_field(1:2) .eq. 'rl')then
                    var_2d = 'rhl'
                    c_label = 'laps relative humidity (bal-liq) %'
                endif

                write(6,*)' reading balanced lh3 / ',var_2d

                call get_3dgrid_dname(directory
     1                  ,i4time_ref,10800,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,nx_l,ny_l,nz_l,rh_3d,istatus)       
                if(istatus .ne. 1)then
                    write(6,*)' data not found'
                    goto100
                endif

                call make_fnam_lp(i4time_nearest,a9time,istatus)

            endif ! c_prodtype

            call interp_3d(rh_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)      

            do k = nz_c,1,-1
            do i = 1,nx_c
                if(field_vert(i,k) .ne. r_missing_data)then
                    field_vert(i,k) = field_vert(i,k)
                endif
            enddo ! i
            enddo ! k

            if(c_field(3:3) .ne. 'i')then
                clow = 0.
                chigh = +100.
                cint = 10. / density
                i_contour = +1
            else ! image
                clow = 220.
                chigh = -40.
                cint = 10.
                i_contour = -1
                scale = 1.
            endif

            nulbll = 1 ! for conrec (number of lines between labels)

            colortable = 'moist'

        elseif(c_field .eq. 'cv')then
            var_2d = 'lcp'
            ext = 'lcp'
            call get_laps_3dgrid(i4time_ref,1000000,i4time_nearest
     1                          ,nx_l,ny_l,nz_l
     1                          ,ext,var_2d,units_2d,comment_2d
     1                          ,rh_3d,istatus)
            if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(rh_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            clow = 0.
            chigh = +1.
            cint = 0.2 / density
            i_contour = 1
            c_label = 'laps cloud fraction              '

            nulbll = 1 ! for conrec (number of lines between labels)

        elseif(c_field .eq. 'la' .or. c_field .eq. 'lj'
     1                         .or. c_field .eq. 'sj'
     1                         .or. c_field .eq. 'ls'
     1                         .or. c_field .eq. 'ss'
     1                         .or. c_field .eq. 'ci'
     1                         .or. c_field .eq. 'cn'
     1                         .or. c_field .eq. 'pc'
     1                         .or. c_field .eq. 'rn'
     1                         .or. c_field .eq. 'sn'
     1                         .or. c_field .eq. 'ic'
     1                                          )then

            call input_product_info(i4time_ref              ! i
     1                             ,laps_cycle_time         ! i
     1                             ,3                       ! i
     1                             ,c_prodtype              ! o
     1                             ,ext                     ! o
     1                             ,directory               ! o
     1                             ,a9time                  ! o
     1                             ,fcst_hhmm               ! o
     1                             ,i4_initial              ! o
     1                             ,i4_valid                ! o
     1                             ,istatus)                ! o

            i4time_lwc = i4_valid

            if(c_field .eq. 'la')then
                c_label = 'laps adiabatic lwc      g/m^3    '
            elseif(c_field .eq. 'lj')then
                c_label = 'laps adjusted  lwc      g/m^3    '
            elseif(c_field .eq. 'sj')then
                c_label = 'laps adjusted  slwc     g/m^3    '
            elseif(c_field .eq. 'ls')then
                c_label = 'cloud liquid            g/m^3    '
            elseif(c_field .eq. 'ci')then
                c_label = 'cloud ice               g/m^3    '
            elseif(c_field .eq. 'cn')then
                c_label = 'cloud condensate        g/m^3    '
            elseif(c_field .eq. 'ss')then
                c_label = 'laps smith-feddes slwc  g/m^3    '
            elseif(c_field .eq. 'pc')then
                c_label = 'precip concentration    g/m^3    '
            elseif(c_field .eq. 'rn')then
                c_label = 'rain concentration      g/m^3    '
            elseif(c_field .eq. 'sn')then
                c_label = 'snow concentration      g/m^3    '
            elseif(c_field .eq. 'ic')then
                c_label = 'ice concentration       g/m^3    '
            endif

            l_pregen = lapsplot_pregen

            if(c_field .eq. 'ls')then
                var_2d = 'lwc'
            elseif(c_field .eq. 'ci')then
                var_2d = 'ice'
            elseif(c_field .eq. 'cn')then
                var_2d = 'lwc'
            elseif(c_field .eq. 'pc')then
                var_2d = 'pcn'
            elseif(c_field .eq. 'rn')then
                var_2d = 'rai'
            elseif(c_field .eq. 'sn')then
                var_2d = 'sno'
            elseif(c_field .eq. 'ic')then
                var_2d = 'pic'
            endif

            if(c_prodtype .eq. 'a')then   
                write(6,*)' getting pregenerated lwc file'
                ext = 'lwc'
                call get_laps_3dgrid(i4time_lwc,86400,i4time_valid,
     1              nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,slwc_3d,istatus)

                if(c_field .eq. 'cn')then
                    write(6,*)' adding ice to get total condensate'
                    var_2d = 'ice'
                    call get_laps_3dgrid(i4time_lwc,86400,i4time_valid,
     1                  nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,cice_3d,istatus)
                    slwc_3d = slwc_3d + cice_3d
                    comment_2d = 'cloud condensate'
                endif

                call make_fnam_lp(i4time_valid,a9time,istatus)

            elseif(c_prodtype .eq. 'b' .or. 
     1             c_prodtype .eq. 'f')then

                call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,slwc_3d
     1                              ,istatus)
                if(istatus .ne. 1)goto100

                if(c_field .eq. 'cn')then
                    write(6,*)' adding ice to get total condensate'
                    var_2d = 'ice'
                    call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,nx_l,ny_l,nz_l       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,cice_3d
     1                              ,istatus)
                    slwc_3d = slwc_3d + cice_3d
                    comment_2d = 'cloud condensate'
                endif

                call downcase(units_2d,units_2d)
                if(units_2d .eq. 'kg/m**3')then
                    units_2d = 'g/m**3'
                endif

                call directory_to_cmodel(directory,c_model)

                write(6,*)' comment_2d = ',trim(comment_2d)

                call mk_fcst_xlabel(comment_2d,fcst_hhmm
     1                                 ,ext(1:3),units_2d
     1                                 ,c_model,c_label)


            endif

            call interp_3d(slwc_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

            scale = 1e-3

            i_contour = 1

            if(i_image .eq. 1)then
                clow = 0.
                if(c_field .eq. 'sn')then
                    chigh = 2.
                elseif(c_field .eq. 'rn')then
                    chigh = 2.
                else
                    chigh = 1.
                endif
                cint = 0.1 
!               cint = 0.0 ! this should force a dynamic setting of clow/chigh
                colortable = 'linear' 
                plot_parms%color_power = 0.3
            else
                clow = 0.
                chigh = 0.
                if(c_field .eq. 'pc')then
                    cint = -0.005 * denslogthr            
                else
                    if(c_prodtype .eq. 'f')then
                        cint = -0.0005 * denslogthr            
		    else
                        cint = -0.0005 * denslogthr                
                    endif
                endif
            endif

        elseif(c_field .eq. 'ix')then

            c_label = '    laps icing index             '

            i4time_lrp = i4time_ref/laps_cycle_time * laps_cycle_time

            if(lapsplot_pregen)then
!           if(.false.)then
                write(6,*)' getting pregenerated lrp file'
                var_2d = 'lrp'
                ext = 'lrp'
                call get_laps_3dgrid(i4time_lrp,86400,i4time_cloud,
     1          nx_l,ny_l,nz_l,ext,var_2d
     1                  ,units_2d,comment_2d,slwc_3d,istatus)
                call interp_3dn(slwc_3d,field_vert,xlow,xhigh,ylow,yhigh
     1                         ,nx_l,ny_l,nz_l,nx_c,nz_c)

                do k = 1,nz_c
                do i = 1,nx_c
                    iarg = nint(field_vert(i,k))
                    icing_index_2d(i,k) = iarg
                enddo ! i
                enddo ! k

            else ! calculate on the fly

            endif ! read pregenerated file

            clow = 0.
            chigh = 0.
            cint = -0.1 * 2. ** (-density)
            i_contour = 4
            call make_fnam_lp(i4time_cloud,a9time,istatus)

        elseif(c_field .eq. 'mv')then
            c_label = 'laps mean volume diameter  m^-6  '

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

            if(.false.)then ! calculate on the fly

            else
                write(6,*)' reading pregenerated mvd'
                var_2d = 'lmd'
                ext =    'lmd'
                call get_laps_3dgrid(i4time_lwc,86400,i4time_cloud,
     1          nx_l,ny_l,nz_l,ext,var_2d
     1          ,units_2d,comment_2d,pcp_type_3d,istatus)

                call interp_3dn(pcp_type_3d,field_vert,xlow,xhigh
     1                         ,ylow,yhigh,nx_l,ny_l,nz_l,nx_c,nz_c)

            endif ! (read pregenerated file)

            do k = 1,nz_c
            do i = 1,nx_c
                field_vert(i,k) = field_vert(i,k) * 1e6 + .01
            enddo ! i
            enddo ! k

            clow = 10.
            chigh = 26.
            cint = 2. / density
            i_contour = 1
            call make_fnam_lp(i4time_cloud,a9time,istatus)

        elseif(c_field .eq. 'tc')then
            c_label = '        laps cloud type          '

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

            write(6,*)' reading pregenerated cloud type'
            var_2d = 'cty'
            ext =    'cty'
            call get_laps_3dgrid(i4time_lwc,86400,i4time_nearest,
     1          nx_l,ny_l,nz_l,ext,var_2d
     1          ,units_2d,comment_2d,pcp_type_3d,istatus)

            call interp_3dn(pcp_type_3d,fieldt_2d,xlow,xhigh
     1                     ,ylow,yhigh,nx_l,ny_l,nz_l,nx_t,nz_c)

            i_contour = 3
            call make_fnam_lp(i4time_nearest,a9time,istatus)

        elseif(c_field .eq. 'tp')then
            c_label = 'laps precip type                  '

            i4time_pcp = i4time_ref

            if(.true.)then

                write(6,*)' reading pregenerated precip type'
                var_2d = 'pty'
                ext =    'pty'
                call get_laps_3dgrid(i4time_pcp,86400,i4time_radar,
     1          nx_l,ny_l,nz_l,ext,var_2d
     1          ,units_2d,comment_2d,pcp_type_3d,istatus)

                call interp_3dn(pcp_type_3d,fieldt_2d,xlow,xhigh
     1                         ,ylow,yhigh,nx_l,ny_l,nz_l,nx_t,nz_c)

            else ! calculate precip type on the fly


            endif ! pregen

            i_contour = 5
            call make_fnam_lp(i4time_radar,a9time,istatus)

        elseif(c_field .eq. 'ix')then
            c_label = '     laps icing severity index   '

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

            if(lapsplot_pregen)then

            else

            endif ! get pregenerated file

        endif ! c_field

        write(6,1605)c_label,a9time
1605    format(2x,a33,2x,a9)

        c_metacode = 'c '

        call i4time_fname_lp(a9time,i4time_dum,istatus)
        call cv_i4tim_asc_lp(i4time_dum,asc_tim_24,istatus)
        asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

        i_label_overlay = i_label_overlay + 1

        if(i_image .eq. 0)then
            i_graphics_overlay = i_graphics_overlay + 1
            call setusv_dum(2hin,icolors(i_graphics_overlay))
        else
            call setusv_dum(2hin,7) ! yellow
        endif

        write(6,*)' plotting, overlay = ',i_graphics_overlay
     1                                   ,i_label_overlay
     1                                   ,' i_image ',i_image
     1                                   ,' i_contour ',i_contour

        ifield_found = 1

        if(abs(i_contour) .eq. 1)then
            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

!           mini = icolors(i_graphics_overlay)   ! is this effective?
!           maxi = icolors(i_graphics_overlay)   ! is this effective?

!           if(i_contour .eq. 1)then
            if(i_image .eq. 0)then
                call setusv_dum(2hin,icolors(i_graphics_overlay)) 
            endif

            vmax = -1e30
            vmin = 1e30

            do k = ibottom,nz_c
            do i = 1,nx_c
                vmax = max(vmax,field_vert(i,k))
                vmin = min(vmin,field_vert(i,k))
            enddo ! k
            enddo ! i

            write(6,*)' clow,high,cint ',clow,chigh,cint
            write(6,*)' max/min = ',vmax,vmin

            if(i_image .eq. 0)then
                write(6,*)' scaling data by ',scale
                do k = nz_c,1,-1
                do i = 1,nx_c
                    if(field_vert(i,k) .ne. r_missing_data)then
                        field_vert(i,k) = field_vert(i,k) / scale       
                    endif
                enddo ! i
                enddo ! k
            endif ! contour plot

            if(cint .ge. 0. .or. i_image .eq. 1)then
                if(.true.)then 
!                   can we make this work (anamorphically) for color plots?
!                   call get_border(ni,nj,x_1,x_2,y_1,y_2)
!                   call set(x_1,x_2,y_1,y_2,0.05,0.95,0.05,0.95,1)

                    call get_border(nx_c,nz_c-ibottom+1,x_1,x_2,y_1,y_2)       

                    call remap_field_2d(
     1                            nx_c,1,nx_c
     1                           ,nz_c,ibottom,nz_c
     1                           ,nx_p, ixl_remap, nx_p-ixh_remap+1
     1                           ,nx_p, 1+iyl_remap, nx_p-iyh_remap
     1                           ,field_vert,field_vert3,r_missing_data)       


!                   if(i_contour .eq. 1)then
                    if(i_image .eq. 0)then ! more generic test hopefully
                        write(6,*)' calling set for conrec_line'
                        call set(0.10,0.90,0.05,0.95
     1                          ,0.10,0.90,0.05,0.95,1) ! orig

                        write(6,*)' arguments for conrec_line, '
     1                           ,'clow/chigh/cint:'       
     1                            ,clow,chigh,cint

                        call conrec_line(field_vert3,nx_p,nx_p,nx_p
     1                             ,clow,chigh,cint,plot_parms
     1                             ,-1,0,-1848,0)

                    else ! image plot
                        call set(0.10,0.90,0.05,0.95
     1                          ,0.10,0.90,0.05,0.95,1) ! orig

!                       call set(x_1,x_2,y_1,y_2,0.15,0.85,0.15,0.85,1) ! new

                        if(cint .le. 0.)then
                            write(6,*)' cint <= 0, calling array_range'
                            call array_range(field_vert3,nx_p,nx_p
     1                                      ,rmin,rmax,r_missing_data)
                            clow = rmin
                            chigh = rmax
                        endif

                        write(6,*)
     1                      ' calling solid fill plot (clow/chigh/cint)'
     1                      ,clow,chigh,cint       

!                       blank out the edges external to the x-section
!                       write(6,*)' blackening the edges'
!                       do i = 1,nx_p
!                       do j = 1,nx_p
!                           if(field_vert3(i,j) .eq. r_missing_data)then       
!                               field_vert3(i,j) = clow
!                           endif
!                       enddo ! j
!                       enddo ! i

!                       if(scale .ne. 0.)then
!                           scale_inv = 1. / scale
!                       endif

!                       this method remaps well - with large border artifacts
!                       if(plot_parms%iraster .eq. 1)then
!                           where(field_vert3(:,:) .eq. r_missing_data)
!    1                            field_vert3(:,:) = -1. ! dark hues color
!                       endif
                        call ccpfil(field_vert3,nx_p,nx_p
     1                             ,clow,chigh
     1                             ,colortable,n_image,scale,'xsect'
     1                             ,plot_parms,namelist_parms)       

!                       this method may avoid the artifacts
!                       call ccpfil(field_vert(1,ibottom),nx_c
!    1                             ,(nz_c-ibottom+1)
!    1                             ,clow,chigh
!    1                             ,colortable,n_image,scale,'xsect'
!    1                             ,plot_parms,namelist_parms)       

                    endif

                endif

            elseif(i_image .eq. 0)then ! logarithmic contour plot
              if(.false.)then
                call conrec(field_vert(1,ibottom)
     1                     ,nx_c,nx_c,(nz_c-ibottom+1)
     1                     ,0.,1e8,1e8,-1,0,-1848,0)

                do i = 1,n_contours
                    cvalue = factor(i)
                    if(cvalue .ge. abs(cint))then
                        call conrec(field_vert(1,ibottom)
     1                             ,nx_c,nx_c,(nz_c-ibottom+1)
     1                             ,cvalue,cvalue,1e-6,-1,0,-1848,0)
                        call conrec(field_vert(1,ibottom)
     1                             ,nx_c,nx_c,(nz_c-ibottom+1)
     1                             ,-cvalue,-cvalue,1e-6,-1,0,-1848,0)
                    endif
                enddo ! i

              else
                write(6,*)' logarithmic contour line plot, cint = ',cint
                write(6,*)' density, denslogthr = ',density,denslogthr
                call remap_field_2d(
     1                            nx_c,1,nx_c
     1                           ,nz_c,ibottom,nz_c
     1                           ,nx_p, ixl_remap, nx_p-ixh_remap+1
     1                           ,nx_p, 1+iyl_remap, nx_p-iyh_remap
     1                           ,field_vert,field_vert3,r_missing_data)       


                call array_range(field_vert3,nx_p,nx_p,rmin,rmax
     1                          ,r_missing_data)

                cmax = max(abs(rmin),abs(rmax))

                call conrec_line(field_vert3,nx_p,nx_p,nx_p
     1                             ,0.,0.,cint,plot_parms,-1,0,-1848,0)       

                do i = 1,n_contours
                    cvalue = factor(i)
                    if(cvalue .ge. abs(cint) .and. cvalue .le. cmax)then       
                        cint_in = 2 * cvalue
                        call conrec_line(field_vert3,nx_p,nx_p,nx_p
     1                             ,-cvalue,cvalue,cint_in,plot_parms
     1                             ,-1,0,-1848,0)       
                    endif
                enddo ! i
 
             endif

           endif ! cint > 0
        endif ! i_contour = 1

!       call upcase(c_label,c_label)
        call set(0., 1., vymin2, vymax2, 0.,1.,0.,1.,1)

!       write bottom label
!       if(i_label_overlay .le. 1)then
!           ity = 35
!           call pwrity(cpux(320),cpux(ity),c_label,33,1,0,0)
!           call pwrity(cpux(800),cpux(ity),asc_tim_24(1:17),17,1,0,0)
!       endif

        call write_label_lplot(100,94,c_label,a9time,plot_parms
     1                        ,namelist_parms,i_label_overlay,'xsect')       

        call make_xsect_labels(      vxmin, vxmax, vymin, vymax, 
     1                               rleft, right, bottom, top,
     1                               width, vymin2, vymax2, i_map, 
     1                               ibottom, r_height,
     1                               xlow,xhigh,ylow,yhigh,
     1                               lat_1d, lon_1d, nx_c, nz_c,
     1                               nx_l, ny_l, 2)

        if(i_contour .eq. 2)then ! plot wind barbs
            call setusv_dum(2hin,icolors(i_graphics_overlay))

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, rleft, right,1)

            barb_factor = float(nx_c) / 249.

            du=(0.4 * barb_factor) / density
            rot = 0.

            iskip_barbs = max(nint( (3.0*barb_factor) / density),1)

            do i = nx_c,1,-iskip_barbs
                rk_terrain = 
     1          max(height_to_zcoord(terrain_vert1d(i),istatus),1.0)
                do k = ibottom,nz_c
                    if(u_vert(i,k) .ne. r_missing_data .and.
     1                 v_vert(i,k) .ne. r_missing_data .and.
     1                abs(u_vert(i,k)) .lt. 1e6      )then
                        x1 = i
                        y1 = (k-ibottom) * float(nx_c-1)
     1                                   / float(nz_c-ibottom) + 1.
                        call uv_to_disp(u_vert(i,k),
     1                                  v_vert(i,k),
     1                                  dir,
     1                                  spd_ms)
                        if(dir .gt. -400. .and. 
     1                       k .ge. int(rk_terrain))then
                            call barbs(spd_ms/mspkt,dir,x1,y1,du,rot
     1                             ,-1e10,+1e10,-1e10,+1e10,1.0,1.0)
                        endif
                    endif
                enddo ! k
            enddo ! i
        endif ! i_contour = 2

        if(i_contour .eq. 3)then ! plot cloud type
            call set(vxmin, vxmax, vymin, vymax
     1             ,    1., float(nx_t), 1., float(nx_t), 1)
            du=.0063 * (121./float(nx_t))
            rot = 0.
            do i = nx_t,1,-2
                do k = ibottom,nz_c
                    i_cloud_type = fieldt_2d(i,k)

                    if(i_cloud_type .ne. 0)then
                        x = i
                        y = (k-ibottom) * float(nx_t-1)
     1                                  / float(nz_c-ibottom) + 1.

                        c2_cloud_type = c2_cloud_types(i_cloud_type)
!                       call upcase(c2_cloud_type,c2_cloud_type)
                        write(6,*)i_cloud_type,c2_cloud_type,x,y
!                       call pwrity (x, y, c2_cloud_type, 2, 0, 0, 0)
                        call pcmequ(x, y, c2_cloud_type,du,0,0)       
                    endif

                enddo ! k
            enddo ! i
        endif ! i_contour = 3

        if(i_contour .eq. 4)then ! plot icing index
            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, rleft, right,1)
            du=0.4
            rot = 0.
            do i = nx_c,1,-1
                do k = ibottom,nz_c
                    if(icing_index_2d(i,k) .ne. 0)then
                        x = i
                        y = (k-ibottom) * float(nx_c-1)
     1                                  / float(nz_c-ibottom) + 1.
!                       write(c1_string,2011)icing_index_2d(i,k)
                        c1_string = '*'
2011                    format(i1)

!                       normal overlay colors
                        call setusv_dum(2hin,icolors(i_graphics_over
     1lay))

!                       multicolored isi display
                        if(.false.)then
                        if(icing_index_2d(i,k) .gt. 3)then
                            i_color_value = icing_index_2d(i,k) - 3
                        else
                            i_color_value = icing_index_2d(i,k)
                        endif

                        if(i_color_value .eq. 1)call setusv_dum(2hin
     1,7) ! y
                        if(i_color_value .eq. 2)call setusv_dum(2hin
     1,245) ! o
                        if(i_color_value .eq. 3)call setusv_dum(2hin
     1,111) ! r
                        endif

                        call pwrity (x, y, c1_string, 1, 0, 0, 0)
                    endif
                enddo ! k
            enddo ! i
        endif ! i_contour = 4


        if(i_contour .eq. 5)then ! plot precip type
            call setusv_dum(2hin,icolors(i_graphics_overlay))

            call set(vxmin, vxmax, vymin, vymax
     1             ,    1., float(nx_t), 1., float(nx_t), 1)
            do i = nx_t,1,-1
                do k = ibottom+1,nz_c
                    i_precip_type = fieldt_2d(i,k)
                    if(i_precip_type .ne. 0)then
                        x = i
                        y = (k-ibottom) * float(nx_t-1)
     1                                  / float(nz_c-ibottom) + 1.

                        c2_precip_type = c1_precip_types(i_precip_type)
                        call pwrity (x, y, c2_precip_type, 2, 0, 0, 0)
                    endif
                enddo ! k
            enddo ! i
        endif ! i_contour = 5


        goto100

9999    continue

!       contour in the terrain surface
        call set(vxmin, vxmax, vymin, vymax
     1         , rleft, right, bottom, top,1)
        n_div = 20
        call setusv_dum(2hin,3)

!       read in sfc pressure
        i4time_tol = 43200
        var_2d = 'ps'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_ref,i4time_tol,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l
     1                      ,pres_2d,0,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading surface pressure analysis'
            write(6,*)
     1        ' converting terrain to sfc pressure with std atmosphere'       
            istat_sfc_pres = 0
        else
            call interp_2d(pres_2d,pres_1d,xlow,xhigh,ylow,yhigh,
     1                     nx_l,ny_l,nx_c,r_missing_data)
            istat_sfc_pres = 1
        endif

        do i = 1,nx_c
            xcoord(i) = i
            if(istat_sfc_pres .eq. 1)then
                ycoord(i) = max(zcoord_of_pressure(pres_1d(i)),1.0)
            else
                ycoord(i) = 
     1            max(height_to_zcoord(terrain_vert1d(i),istatus),1.0)
            endif

            if(i .gt. 1)then
                do j = 0,n_div
                    frac = float(j) / float(n_div)
                    ybottom = bottom
                    ytop = ycoord(i-1) * (1.-frac) + ycoord(i) * frac
                    xval = float(i-1) + frac
                    call line(xval,ybottom,xval,ytop)
                enddo ! j
            endif
        enddo ! i

        call setusv_dum(2hin,7)

        if(l_sta)then ! label location of station

!          this lets up plot outside the main box
           call set(.00, 1.0, vymin2 , vymax2
     1             , rleft  - width/8.   , right + width/8.
     1             , bottom - r_height/8., top   + r_height/8. ,1)

           x = pos_sta
           write(6,*)'     labeling ',c3_sta,pos_sta

           call line(x,bottom,x,bottom - .004 * r_height)

           y = bottom - .010 * r_height
           call upcase(c3_type,c3_type)

           if(l_arrival_gate)then
               write(c9_string,2039)c3_sta,c5_arrival_gate(1:3)
           else
               write(c9_string,2039)c3_sta,c3_type
           endif

2039       format(1x,a3,'-',a3)

!          call pwrity ((x+0.15), y, c9_string, 9, 0, 0, 0)
           call pcmequ ((x+0.15), y, c9_string, .0050, 0, 0)

           if(.not. l_atms)then
               i4time_label = i4time_ref/laps_cycle_time*laps_cycle_time
!    1                                          -laps_cycle_time

               call label_other_stations(i4time_label,standard_longitude       
     1                                  ,y,xsta,lat,lon,nx_l,ny_l
     1                                  ,grid_spacing_m
     1                                  ,xlow,xhigh,ylow,yhigh
     1                                  ,nx_c,bottom,r_height,maxstns)

           endif

        else ! l_sta = .false.

           if(.not. l_atms)then

!              this lets up plot outside the main box
               call set(.00, 1.0, vymin2 , vymax2
     1                , rleft - width/8., right + width/8.
     1                , bottom - r_height/8., top + r_height/8.,1)

               xsta = -10000.

               i4time_label = i4time_ref/laps_cycle_time*laps_cycle_time
!    1                                          -laps_cycle_time

               y = bottom - .015 * r_height

               call label_other_stations(i4time_label,standard_longitude       
     1                                  ,y,xsta,lat,lon,nx_l,ny_l
     1                                  ,grid_spacing_m
     1                                  ,xlow,xhigh,ylow,yhigh,nx_c
     1                                  ,bottom,r_height,maxstns)

           endif

        endif ! l_sta = .true.

!       call frame

        return
        end


        subroutine interp_3d(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                       nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

!       97-aug-14     ken dritz     added nx_l, ny_l, nz_l, nx_c, nz_c as
!                                   dummy arguments
!       97-aug-14     ken dritz     added r_missing_data as dummy argument
!       97-aug-14     ken dritz     removed (commented out) parameter
!                                   declarations for nx_c, nz_c
!       97-aug-14     ken dritz     removed include of lapsparms.for
!       97-aug-14     ken dritz     pass nx_l, ny_l, nx_c, r_missing_data
!                                   to interp_2d

        integer nx_c,nz_c
!       parameter (nx_c = 61) ! nx_l
!       parameter (nz_c = nz_l)

        real array_in(nx_l,ny_l,nz_l)
        real array_out(nx_c,nz_c)

!       array_out = r_missing_data

        do k = 1,nz_l
             call interp_2d(array_in(1,1,k),array_out(1,k)
     1                     ,xlow,xhigh,ylow,yhigh
     1                     ,nx_l,ny_l,nx_c,r_missing_data)
        enddo ! i

        return
        end

        subroutine interp_3dn(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                        nx_l,ny_l,nz_l,nx_c,nz_c)

!       97-aug-14     ken dritz     added nx_l, ny_l, nz_l, nx_c, nz_c as
!                                   dummy arguments
!       97-aug-14     ken dritz     removed (commented out) parameter
!                                   declarations for nx_c, nz_c
!       97-aug-14     ken dritz     removed include of lapsparms.for
!       97-aug-14     ken dritz     pass nx_l, ny_l, nx_c to interp_2dn

!       nearest neighbor interpolation of a vertical x-sect

        integer nx_c,nz_c
!       parameter (nx_c = 61) ! nx_l
!       parameter (nz_c = nz_l)

        real array_in(nx_l,ny_l,nz_l)
        real array_out(nx_c,nz_c)

        do k = 1,nz_l
             call interp_2dn(array_in(1,1,k),array_out(1,k)
     1                                 ,xlow,xhigh,ylow,yhigh,
     1                       nx_l,ny_l,nx_c)
        enddo ! i

        return
        end

        subroutine interp_3d_spread(array_in,array_out,xlow,xhigh,ylow,y
     1high,nx_l,ny_l,nz_l,nx_c,nz_c,r_missing_data)

!       97-aug-14     ken dritz     added nx_l, ny_l, nz_l, nx_c, nz_c as
!                                   dummy arguments
!       97-aug-14     ken dritz     added r_missing_data as dummy argument
!       97-aug-14     ken dritz     removed (commented out) parameter
!                                   declarations for nx_c, nz_c
!       97-aug-14     ken dritz     removed include of lapsparms.for
!       97-aug-14     ken dritz     pass nx_l, ny_l, nx_c, r_missing_data
!                                   to interp_2d_spread

        integer nx_c,nz_c
!       parameter (nx_c = 61) ! nx_l
!       parameter (nz_c = nz_l)

        real array_in(nx_l,ny_l,nz_l)
        real array_out(nx_c,nz_c)

        do k = 1,nz_l
             call interp_2d_spread(array_in(1,1,k),array_out(1,k)
     1                                 ,xlow,xhigh,ylow,yhigh,
     1                             nx_l,ny_l,nx_c,r_missing_data)
        enddo ! i

        return
        end


        subroutine interp_2d(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                       nx_l,ny_l,nx_c,r_missing_data)

!       97-aug-14     ken dritz     added nx_l, ny_l, nx_c, and r_missing_data
!                                   as dummy arguments
!       97-aug-14     ken dritz     removed (commented out) parameter
!                                   declarations for nx_c, nz_c (the latter
!                                   was unused)
!       97-aug-14     ken dritz     removed include of lapsparms.for

        integer nx_c,nz_c
!       parameter (nx_c = 61) ! nx_l
!       parameter (nz_c = nz_l)

        real array_in(nx_l,ny_l)
        real array_out(nx_c)

        if(nx_c .gt. 1)then
            deltax = (xhigh - xlow) / (float(nx_c) - 1.)
            deltay = (yhigh - ylow) / (float(nx_c) - 1.)
        else
            deltax = 0.
            deltay = 0.
        endif

!       bilinearly interpolate from 2d grid to 1d array
        do ii = 1,nx_c
            ri = xlow + (float(ii-1) * deltax)
            rj = ylow + (float(ii-1) * deltay)

            if(ri .lt. 1.0)ri = 1.0
            if(rj .lt. 1.0)rj = 1.0

            i = int(ri)
            if(i .eq. nx_l)i=i-1

            j = int(rj)
            if(j .eq. ny_l)j=j-1

            fraci = ri - i
            fracj = rj - j

            z1=array_in(i  , j  )
            z2=array_in(i+1, j  )
            z3=array_in(i+1, j+1)
            z4=array_in(i  , j+1)

            if(  z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                array_out(ii) =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                - (z2+z4-z3-z1)*fraci*fracj

            else
                array_out(ii) = r_missing_data

            endif

!            if(array_out(ii) .ne. r_missing_data
!       1       .and. abs(array_out(ii) .lt. 1e-20)then
!                array_out(ii) = r_missing_data
!            endif

        enddo ! ii

        return
        end


        subroutine interp_2dn(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                        nx_l,ny_l,nx_c)

!       97-aug-14     ken dritz     added nx_l, ny_l, nx_c as dummy arguments
!       97-aug-14     ken dritz     removed (commented out) parameter
!                                   declarations for nx_c, nz_c (the latter
!                                   was unused)
!       97-aug-14     ken dritz     removed include of lapsparms.for

!       nearest neighbor interpolation of a 2d array to a line.

        integer nx_c,nz_c
!       parameter (nx_c = 61) ! nx_l
!       parameter (nz_c = nz_l)

        real array_in(nx_l,ny_l)
        real array_out(nx_c)

        deltax = (xhigh - xlow) / (float(nx_c) - 1.)
        deltay = (yhigh - ylow) / (float(nx_c) - 1.)

!       interpolate from 2d grid to 1d array using the nearest neighbor
        do ii = 1,nx_c
            ri = xlow + (float(ii-1) * deltax)
            rj = ylow + (float(ii-1) * deltay)

            i = nint(ri)
            j = nint(rj)

            array_out(ii) = array_in(i,j)

        enddo ! ii

        return
        end


        subroutine interp_2d_spread(array_in,array_out,xlow,xhigh
     1                             ,ylow,yhigh,nx_l,ny_l,nx_c
     1                             ,r_missing_data)

!       97-aug-14     ken dritz     added nx_l, ny_l, nx_c, r_missing_data
!                                   as dummy arguments
!       97-aug-14     ken dritz     removed (commented out) parameter
!                                   declarations for nx_c and nz_c (the
!                                   latter was unused)
!       97-aug-14     ken dritz     removed include of lapsparms.for

        integer nx_c,nz_c
!       parameter (nx_c = 61) ! nx_l
!       parameter (nz_c = nz_l)

        real array_in(nx_l,ny_l)
        real array_out(nx_c)

        deltax = (xhigh - xlow) / (float(nx_c) - 1.)
        deltay = (yhigh - ylow) / (float(nx_c) - 1.)

!       bilinearly interpolate from 2d grid to 1d array
        do ii = 1,nx_c
            ri = xlow + (float(ii-1) * deltax)
            rj = ylow + (float(ii-1) * deltay)

            if(ri .lt. 1.0)ri = 1.0
            if(rj .lt. 1.0)rj = 1.0

            i = int(ri)
            if(i .eq. nx_l)i=i-1

            j = int(rj)
            if(j .eq. ny_l)j=j-1

            fraci = ri - i
            fracj = rj - j

            z1=array_in(i  , j  )
            z2=array_in(i+1, j  )
            z3=array_in(i+1, j+1)
            z4=array_in(i  , j+1)

            if(  z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                array_out(ii) =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                - (z2+z4-z3-z1)*fraci*fracj

            else
                array_out(ii) = r_missing_data

            endif

!            if(array_out(ii) .ne. r_missing_data
!       1       .and. abs(array_out(ii) .lt. 1e-20)then
!                array_out(ii) = r_missing_data
!            endif

        enddo ! ii

        return
        end

        subroutine interp_3dc(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                        nx_l,ny_l,nx_c,r_missing_data)

!       97-aug-14     ken dritz     added nx_l, ny_l, nx_c, r_missing_data
!                                   as dummy arguments
!       97-aug-14     ken dritz     removed (commented out) parameter
!                                   declarations for nx_c, nz_c (the
!                                   latter was unused)
!       97-aug-14     ken dritz     removed include of lapsparms.for
!       97-aug-14     ken dritz     pass nx_l, ny_l, nx_c, r_missing_data
!                                   to interp_2d

        include 'laps_cloud.inc'

        integer nx_c,nz_c
!       parameter (nx_c = 61) ! nx_l
!       parameter (nz_c = kcloud)

        real array_in(nx_l,ny_l,kcloud)
        real array_out(nx_c,kcloud)

        do k = 1,kcloud
             call interp_2d(array_in(1,1,k),array_out(1,k)
     1                                 ,xlow,xhigh,ylow,yhigh,
     1                      nx_l,ny_l,nx_c,r_missing_data)
        enddo ! i

        return
        end


        subroutine xsect_endpoints(xsta,ysta,azi_xsect,
     1                  xlow,ylow,xhigh,yhigh,pos_sta,istatus,
     1                  nx_l,ny_l,nx_c)

!       97-aug-14     ken dritz     added nx_l, ny_l, nx_c as dummy arguments
!       97-aug-14     ken dritz     removed (commented out) parameter
!                                   declarations for nx_c, nz_c (the latter
!                                   was unused)
!       97-aug-14     ken dritz     removed include of lapsparms.for

        include 'trigd.inc'

        logical l_left, l_right, l_top, l_bottom

        integer nx_c,nz_c
!       parameter (nx_c = 61) ! nx_l
!       parameter (nz_c = nz_l)

        angdif(x,y)=mod(x-y+540.,360.)-180.
        cotand(x) = tand(90.-x)

!           calculate endpoints of x-sect from waypoint and azimuth
            azi_met = 90. - azi_xsect

!           intersection with right edge
            l_right = .false.
            if(cotand(azi_met) .ne. 0.)then
                yright = ysta + (nx_l - xsta) * tand(azi_met)
                if(yright .le. float(ny_l) .and. yright .ge. 1.)then
                    l_right = .true.
                    write(6,311)float(nx_l),yright
311                 format('  right edge  ',2f7.2)
                endif
            endif

!           intersection with left edge
            l_left = .false.
            if(cotand(azi_met) .ne. 0.)then
                yleft = ysta - (xsta - 1.) * tand(azi_met)
                if(yleft .le. float(ny_l) .and. yleft .ge. 1.)then
                    l_left = .true.
                    write(6,312)1.,yleft
312                 format('  left edge   ',2f7.2)
                endif
            endif

!           intersection with top edge
            l_top = .false.
            if(tand(azi_met) .ne. 0.)then
                xtop = xsta + (ny_l - ysta) * cotand(azi_met)
                if(xtop .le. float(nx_l) .and. xtop .ge. 1.)then
                    l_top = .true.
                    write(6,313)xtop,float(ny_l)
313                 format('  top edge    ',2f7.2)
                endif
            endif

!           intersection with bottom edge
            l_bottom = .false.
            if(tand(azi_met) .ne. 0.)then
                xbottom = xsta - (ysta - 1.) * cotand(azi_met)
                if(xbottom .le. float(nx_l) .and. xbottom .ge. 1.)then
                    l_bottom = .true.
                    write(6,314)xbottom,1.
314                 format('  bottom edge ',2f7.2)
                endif
            endif

!           now we can get the endpoints
            if(abs(angdif(azi_xsect,45.)) .le. 45.)then
                if(l_left)then
                    xlow = 1.
                    ylow = yleft
                elseif(l_bottom)then
                    xlow = xbottom
                    ylow = 1.
                endif

                if(l_right)then
                    xhigh = nx_l
                    yhigh = yright
                elseif(l_top)then
                    xhigh = xtop
                    yhigh = ny_l
                endif

            elseif(abs(angdif(azi_xsect,135.)) .le. 45.)then
                if(l_left)then
                    xlow = 1.
                    ylow = yleft
                elseif(l_top)then
                    xlow = xtop
                    ylow = ny_l
                endif

                if(l_right)then
                    xhigh = nx_l
                    yhigh = yright
                elseif(l_bottom)then
                    xhigh = xbottom
                    yhigh = 1
                endif

            elseif(abs(angdif(azi_xsect,225.)) .le. 45.)then
                if(l_left)then
                    xhigh = 1.
                    yhigh = yleft
                elseif(l_bottom)then
                    xhigh = xbottom
                    yhigh = 1.
                endif

                if(l_right)then
                    xlow = nx_l
                    ylow = yright
                elseif(l_top)then
                    xlow = xtop
                    ylow = ny_l
                endif

            elseif(abs(angdif(azi_xsect,315.)) .le. 45.)then
                if(l_left)then
                    xhigh = 1.
                    yhigh = yleft
                elseif(l_top)then
                    xhigh = xtop
                    yhigh = ny_l
                endif

                if(l_right)then
                    xlow = nx_l
                    ylow = yright
                elseif(l_bottom)then
                    xlow = xbottom
                    ylow = 1
                endif

            endif

            if(xlow .ne. xhigh)then
                pos_sta = 1. + (nx_c-1.) * (xsta-xlow)/(xhigh-xlow)
            else
                pos_sta = 1. + (nx_c-1.) * (ysta-ylow)/(yhigh-ylow)
            endif

            write(6,91)xlow,ylow,xhigh,yhigh,pos_sta
91          format(1x,' end points  ',2f7.2,2x,2f7.2,'  pos_sta',f7.2)


        return
        end

        subroutine label_other_stations(i4time,standard_longitude,y,xsta
     1          ,lat,lon,ni,nj
     1          ,grid_spacing_m
     1          ,xlow,xhigh,ylow,yhigh,nx_c,bottom,r_height,maxstns)

!       97-aug-14     ken dritz     added maxstns as dummy argument
!       97-aug-14     ken dritz     removed include of lapsparms.for
!       97-aug-25     steve albers  removed /read_sfc_cmn/.

!       this routine labels stations on the x-sect in a logical manner

        real stapos_a(maxstns+1)

        real lat(ni,nj),lon(ni,nj)

        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cover_s(maxstns), hgt_ceil(maxstns), hgt_low(maxstns)
        real t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns)
     1       , ffg_s(maxstns)
        real vis_s(maxstns)
        character stations(maxstns)*3, wx_s(maxstns)*8    ! c5_stamus

!       declarations for new read_surface routine
!       new arrays for reading in the sao data from the lso files
        real   pstn(maxstns),pmsl(maxstns),alt(maxstns)
     1          ,store_hgt(maxstns,5)
        real   ceil(maxstns),lowcld(maxstns),cover_a(maxstns)
     1          ,vis(maxstns)
     1                                          ,rad(maxstns)

        integer   obstime(maxstns),kloud(maxstns),idp3(maxstns)

        character   obstype(maxstns)*8
     1             ,store_emv(maxstns,5)*1,store_amt(maxstns,5)*4

!       common /read_sfc_cmn/ lat_s,lon_s,elev_s,cover_s,hgt_ceil,hgt_lo
!    1w
!    1                ,t_s,td_s,pr_s,sr_s,dd_s,ff_s,ddg_s,ffg_s,vis_s
c
        character atime*24, infile*255
        character directory*150,ext*31

        character*9 c9_string
        character*13 filename13

        character*2 icompass(8)
        data icompass/'n ','ne','e ','se','s ','sw','w ','nw'/

        write(6,*)
     1  ' reading station locations from read_sfc for labelling '
        ext = 'lso'
        call get_directory(ext,directory,len_dir) ! returns top level directory
        infile = directory(1:len_dir)//filename13(i4time,ext(1:3))

        call read_surface_old(infile,maxstns,atime,n_meso_g,n_meso_pos,
     &      n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g,
     &      n_obs_b,n_obs_pos_b,stations,obstype,lat_s,lon_s,
     &      elev_s,wx_s,t_s,td_s,dd_s,ff_s,ddg_s,
     &      ffg_s,pstn,pmsl,alt,kloud,ceil,lowcld,cover_a,rad,idp3,
     &      store_emv,store_amt,store_hgt,vis,obstime,istatus)

100     write(6,*)'     n_obs_b',n_obs_b

        if(n_obs_b .gt. maxstns .or. istatus .ne. 1)then
            write(6,*)' warning: too many stations, or no file present'
            istatus = 0
            return
        endif

        stapos_a(1) = xsta

        sect_length = sqrt((xhigh-xlow)**2 + (yhigh-ylow)**2)

        write(6,*)'i,xsta,ysta,frac,rmin,stations(i),'
     1                            ,'stapos,dist_min,ran,azi'

        iplot = 1

        do isweep = 0,7
          do i = 1,n_obs_b
            call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon
     1                          ,ni,nj,xsta,ysta,istatus)

!           calculate distance from station to x-sect
!           rmin is in terms of laps grid points
            call closest(xlow-xsta,ylow-ysta,xhigh-xlow,yhigh-ylow,frac,
     1rmin)

            stapos = 1. + frac * float(nx_c-1)

!           rmin = rmin * float(nx_c-1) / sect_length

            if(frac .ge. 0. .and. frac .le. 1.)then
              if(int(abs(rmin*2.)) .eq. isweep)then

!               find distance from previously written out locations
                dist_min = min(abs(nx_c-stapos),abs(stapos-1.))
                do ii = 1,iplot
                    dist = abs(stapos-stapos_a(ii))
                    if(dist .lt. dist_min)then
                        dist_min = dist
                    endif
                enddo ! ii


                if(dist_min .gt. 5.0)then ! no other stations in the way
                    xclo = xlow + frac * (xhigh - xlow)
                    yclo = ylow + frac * (yhigh - ylow)

                    iplot = iplot + 1
                    stapos_a(iplot) = stapos

                    xdelt = xclo - xsta
                    ydelt = yclo - ysta

                    call xy_to_met_xm(xdelt,ydelt,rang,azi,istatus)

                    ran = rang * (grid_spacing_m / 1000.) ! station dist in km

                    azi = azi + (lon_s(i) - standard_longitude)

                    i_cardinal_pt = mod(int((azi+22.5)/45.),8) + 1

                    if(icompass(i_cardinal_pt)(2:2) .eq. ' ')then
                        write(c9_string,2038)nint(ran)
     1                    ,icompass(i_cardinal_pt),stations(i)(1:3)
2038                    format(i2,a2,a3)
!                       call pwrity ((stapos+0.6), y, c9_string
!    1                                                  , 9, 0, 0, 0)
                        call pcmequ ((stapos+0.8), y, c9_string
     1                               , .0060, 0, 0)
                    else
                        write(c9_string,2039)nint(ran)
     1                    ,icompass(i_cardinal_pt),stations(i)(1:3)
2039                    format(i2,a2,' ',a3)

!                       call pwrity ((stapos+.15), y, c9_string
!    1                                                  , 9, 0, 0, 0)
                        call pcmequ ((stapos+.35), y, c9_string
     1                               , .0060, 0, 0)
                    endif


                    write(6,11)i,xsta,ysta,frac,rmin,stations(i)(1:3)
     1                  ,stapos,dist_min
     1                  ,ran,azi,c9_string,isweep
11                  format(i4,2f6.1,f6.3,f6.1,1x,a3,4f6.1,1x,a9,i2)

                    call line(stapos,bottom,stapos
     1                       ,bottom - .012 * r_height)
!                   write(6,*) ' call line',stapos,bottom,r_height

                endif

              endif
            endif
          enddo ! stations
        enddo ! isweep

        return
        end

        subroutine closest(r1,r2,rd1,rd2,tclo,rmin)
        implicit real (a-z)

c       write(6,*)

        rnum =  -(      r1 * rd1
     1  +       r2 * rd2)
        denom =  (      rd1 * rd1
     1  +       rd2 * rd2)

        tclo = rnum/denom ! + t0

        rdot = sqrt(rd1**2 + rd2**2)
        rnum = (r1*rd2) - (r2*rd1)
        rmin = rnum/rdot
c       write(6,*)' rnum,rdot = ',rnum,rdot

c       write(6,1)rmin,r1,r2,rd1,rd2,tclo
1       format(f12.10,5f12.8)
        return
        end

c

        subroutine remap_field_2d(nx_in,ixlow_in,ixhigh_in
     1                           ,ny_in,iylow_in,iyhigh_in
     1                           ,nx_out,ixlow_out,ixhigh_out
     1                           ,ny_out,iylow_out,iyhigh_out
     1                           ,field_in,field_out,r_missing_data)

        real field_in(nx_in,ny_in)
        real field_out(nx_out,ny_out)

        call constant(field_out,r_missing_data,nx_out,ny_out)

        rxlow_out  = ixlow_out
        rxhigh_out = ixhigh_out
        rylow_out  = iylow_out
        ryhigh_out = iyhigh_out

        rxlow_in   = ixlow_in
        rxhigh_in  = ixhigh_in 
        rylow_in   = iylow_in 
        ryhigh_in  = iyhigh_in

        write(6,*)' remap_field_2d x:'
     1            ,rxlow_out,rxhigh_out,rxlow_in,rxhigh_in
        write(6,*)' remap_field_2d y:'
     1            ,rylow_out,ryhigh_out,rylow_in,ryhigh_in

!       loop over subset of output (large square) array
        iytest = ny_out/2

        do ixout = ixlow_out, ixhigh_out
        do iyout = iylow_out, iyhigh_out

!           determine indices of input (smaller) array
            rxout = ixout
            ryout = iyout

            arg = rxout
            call stretch2(rxlow_out,rxhigh_out,rxlow_in,rxhigh_in,arg)      
            rxin = arg

            arg = ryout
            call stretch2(rylow_out,ryhigh_out,rylow_in,ryhigh_in,arg)
            ryin = arg

!           interpolate input array to find value of output field array element
            call bilinear_laps(rxin,ryin,nx_in,ny_in,field_in,result)       
            field_out(ixout,iyout) = result

            if(iyout .eq. iytest .and. .false.)then
                write(6,*)ixout,iyout,rxin,nx_in,nx_out
     1                   ,field_in(nint(rxin),nint(ryin))
     1                   ,result
            endif

        enddo
        enddo

        return
        end

        subroutine mk_fcst_xlabel(comment_2d,fcst_hhmm_in
     1                                ,ext
     1                                ,units_2d
     1                                ,c_model_in
     1                                ,c_label)

        character*(*) comment_2d,ext,units_2d,c_model_in
        character*100 c_label, c_model

        character*5 fcst_hhmm_in,fcst_hhmm

        len_model_max = 7

!       initialize with blanks
        do i = 1,100
            c_label(i:i) = ' '
        enddo ! i

!       call downcase(units_2d,units_2d)

        call s_len2(comment_2d,len_fcst)
        call s_len2(units_2d,len_units)

        if(ext .eq. 'lga')then
            c_model = 'lga'
        endif

        idash = index(c_model_in,'-')
        write(6,*)' c_model_in/idash = ',c_model_in,idash
        if(idash .gt. 0)then
            c_model = c_model_in(idash+1:len(c_model_in))
        else
            c_model = c_model_in
        endif

        call s_len2(c_model,len_model)
        call upcase(c_model,c_model)

        call s_len(fcst_hhmm_in,length_fcst_in)

!        write(c_label,102)k_mb
!102     format(i5,' hpa ')

        if(fcst_hhmm_in(length_fcst_in-1:length_fcst_in) .eq. '00')then      
            fcst_hhmm = fcst_hhmm_in(1:length_fcst_in-2)//'hr '
        else
            fcst_hhmm = fcst_hhmm_in
        endif

        ist = 27

!       field and units info
        if(len_units .gt. 0)then
            c_label(1:ist-1) = comment_2d(1:len_fcst)
     1                         //' ('//units_2d(1:len_units)//')'
        else
            c_label(1:ist-1) = comment_2d(1:len_fcst)
        endif

!       fcst time info
        c_label(ist:ist+length_fcst_in+1) = 
     1                fcst_hhmm(1:length_fcst_in)//' '

!       model info
        if(len_model .gt. 0)then
            len_model = min(len_model,len_model_max) 
            c_label(ist+5:ist+11+len_model) = 
     1                            c_model(1:len_model)//' fcst'       
        else
            c_label(ist+5:ist+8) = 'fcst'       
        endif

        len_label = len(c_label)
        if(len_label .le. 200)then
            write(6,*)' mk_fcst_xlabel:',c_label
        else
            write(6,*)' mk_fcst_xlabel: label length = ',len_label
        endif

        return
        end

        subroutine make_xsect_labels(vxmin, vxmax, vymin, vymax, 
     1                               rleft, right, bottom, top,
     1                               width, vymin2, vymax2, i_map, 
     1                               ibottom, r_height,
     1                               xlow,xhigh,ylow,yhigh,
     1                               lat_1d, lon_1d, nx_c, nz_c,
     1                               nx_l, ny_l, i_pos)

        real xcoord(nx_c),ycoord(nx_c)

        character*4 c4_string
        character*7 c7_string

        real lat_1d(nx_c)
        real lon_1d(nx_c)

        write(6,*)' subroutine make_xsect_labels: i_pos, i_map = '
     1           ,i_pos,i_map

        if(i_map .eq. 0)then

            i_map = 1

            call setusv_dum(2hin,7)

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

!           this lets us plot outside the main box
            call set(.00, 1.0, vymin2 , vymax2
     1             , rleft - width/8., right + width/8.,
     1               bottom - r_height/8., top + r_height/8.,1)

!           draw box enclosing graph
            xcoord(1) = rleft
            ycoord(1) = bottom
            xcoord(2) = right
            ycoord(2) = bottom
            xcoord(3) = right
            ycoord(3) = top
            xcoord(4) = rleft
            ycoord(4) = top
            xcoord(5) = xcoord(1)
            ycoord(5) = ycoord(1)
            npts = 5
            call curve (xcoord, ycoord, npts)

!           label left axis
            if(.true.)then ! label pressure on left axis
                if(nz_c .gt. 60)then
                    iskip = 2
                else
                    iskip = 1
                endif

                do i = ibottom,nz_c
                    y = i
                    call line (rleft, y, rleft + width * .015, y )

!                   pressure
                    if(i-1 .eq. ((i-1)/iskip)*iskip)then
                        x = rleft - width * .030
                        ipres_mb = nint(zcoord_of_level(i)/100.)
                        write(c4_string,2014)ipres_mb
                        call pwrity (x, y, c4_string, 4, 0, 0, 0)
2014                    format(i4)
                    endif

                end do
                call pwrity (rleft - .070 * width,bottom + r_height*0.5,
     1          ' pressure (hpa) ',16,1,90,0)
!           endif


!           if(.not. l_atms)then
!               label height (km - msl) on right axis
                do i = 0,16
                    y = height_to_zcoord(i*1000.,istatus)

                    if(y .ge. bottom)then
                        call line (right, y, right - width * .015, y )

!                       height
                        x = right + width * .015
                        iht_km = i
                        write(c4_string,2014)iht_km
                        call pwrity (x, y, c4_string, 4, 0, 0, 0)
                    endif

                end do
                call pwrity (right + .070 * width,bottom + r_height*0.5,
     1          'height  (km msl)',16,1,270,0)
            endif


!           put in lat/lon of endpoints
            x = rleft-1.0
            y = bottom - r_height * .03
            write(c7_string,2017)lat_1d(1)
            call pwrity (x, y, c7_string, 7, 0, 0, 0)
            y = bottom - r_height * .05
            write(c7_string,2017)lon_1d(1)
            call pwrity (x, y, c7_string, 7, 0, 0, 0)

            x = right+1.0
            y = bottom - r_height * .03
            write(c7_string,2017)lat_1d(nx_c)
            call pwrity (x, y, c7_string, 7, 0, 0, 0)
            y = bottom - r_height * .05
            write(c7_string,2017)lon_1d(nx_c)
            call pwrity (x, y, c7_string, 7, 0, 0, 0)

2017        format(f7.2)

!           put in "minibox"
            call get_border(nx_l,ny_l,x_1,x_2,y_1,y_2)
            size_mini = .05
            size_mini_x = size_mini * width * 0.8 ! aspect ratio of xsect
            size_mini_y = size_mini * r_height
            rleft_mini  = rleft + width    * .28 ! .26
            bottom_mini = top   + r_height * .025
            xl = rleft_mini  + x_1 * size_mini_x
            xh = rleft_mini  + x_2 * size_mini_x
            yl = bottom_mini + y_1 * size_mini_y
            yh = bottom_mini + y_2 * size_mini_y
            xcoord(1) = xl  ! ul
            ycoord(1) = yh
            xcoord(2) = xh  ! ur
            ycoord(2) = yh
            xcoord(3) = xh  ! lr
            ycoord(3) = yl
            xcoord(4) = xl  ! ll
            ycoord(4) = yl
            xcoord(5) = xcoord(1)
            ycoord(5) = ycoord(1)
            npts = 5
            call curve (xcoord, ycoord, npts)

!           get x-section line segment
            xfracl = (xlow  - 1.0) / float(nx_l)
            xfrach = (xhigh - 1.0) / float(nx_l)
            yfracl = (ylow  - 1.0) / float(ny_l)
            yfrach = (yhigh - 1.0) / float(ny_l)

            x1 = xl + xfracl*(xh-xl)
            x2 = xl + xfrach*(xh-xl)
            y1 = yl + yfracl*(yh-yl)
            y2 = yl + yfrach*(yh-yl)

            call line(x1,y1,x2,y2)

        endif ! i_map .eq. 0

        return
        end
