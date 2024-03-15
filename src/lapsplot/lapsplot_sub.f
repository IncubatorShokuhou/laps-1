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

        subroutine lapsplot_3d()

        include 'lapsplot.inc'

        character*1 c_display
        integer aski4t
        character*9 asc9_tim
        character*4 rm
        character*35 time
        character*8 czoom
        logical l_atms, l_plotobs
        logical l_parse, l_polar, l_cyl

        character*2 c_section

        common/ramtdims/nxmin,nxmax,nymin,nymax

        common /zoom/ zoom,density

!       integer nk
!       common /get_packed_data2/ nk

        real, allocatable :: dx(:,:)            
        real, allocatable :: dy(:,:)            

!       nk = 17 ! a temporary fix as the blockdata won't work

        call get_grid_dim_xy(nx_l,ny_l,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting horizontal domain dimensions'
           stop
        endif

        call get_laps_dimensions(nz_l,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting vertical domain dimension'
           stop
        endif

        call get_max_radars(max_radars,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting value of max_radars'
           stop
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting value of r_missing_data'
           stop
        endif

        call get_standard_longitude(standard_longitude,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting value of standard_longitude'
           stop
        endif

        call get_maxstns(maxstns,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting value of maxstns'
           stop
        endif

        dyn_low = r_missing_data
        dyn_high = r_missing_data

        allocate (dx(nx_l,ny_l),dy(nx_l,ny_l))

!       read(5,2)c_display
2       format(a1)
        c_display = 'r'

        if(c_display .eq. 'r')then
            write(6,*)
            write(6,*)
            write(6,*)' this program will produce a gmeta file for displ
     1ay'
            write(6,*)
            write(6,*)
     1    ' it will be located in ./gmeta'
        endif

!       read namelist parameters
        call get_lapsplot_parms(namelist_parms,istatus)       

!       set default plot parameters
        plot_parms%rimage_intensity = 1.0
        plot_parms%zoom = 1.0
        plot_parms%xcen = 0.5
        plot_parms%ycen = 0.5
        plot_parms%zoom_wdw = 1.0
        plot_parms%obs_size = 1.0
        plot_parms%ncols = 0.
        plot_parms%icol_barbs = namelist_parms%icol_barbs
        plot_parms%l_hinterp_zoom = .false.

        if(.true.)then       ! testcode
!           call gopks(6,10000)
!           call gopwk(1,51,3)
!           call gacwk(1)

            call gopks (6,idum)
            call gopwk (1,2,1)
            call gacwk (1)

            write(6,*)' calling color'
            iwhite = namelist_parms%i_background_color
            call color(iwhite)

            write(6,*)' calling gopwk'
            call gopwk (9,1,3)

        else
            call opngks
        endif

!       supply reference time for matching the data
        lun = 5
1090    write(6,*)
        write(6,*)
        i4time_ref = aski4t()
!       i4time_ref = (i4time_ref)/60*60

1100    write(6,1110)
1110    format(/////'     [h/hz]  horizontal plan view '
     1        ,'   (const pressure level or sfc)'
     1        /'     [x]  vertical cross section'
     1        /'     [s]  sounding'
     1        /'     [nt]  new time'
     1        /' ',60x,'[q] quit ? ',$)
        read(lun,1111)c_section
1111    format(a2)

        write(6,*)' c_section is :',c_section

        if(c_section .eq. 'nt')then
!           write(6,*)' setting namelist_parms%iraster to -1'
!           write(6,*)' raster may not work with multiple frames'
!           namelist_parms%iraster = -1
            goto1090
        endif

        if(c_section .eq. 'q' .or. c_section .eq. 'q')goto999

        call get_laps_cycle_time(laps_cycle_time, istatus)
        if(istatus .ne. 1)then
            write(6,*)' bad call to get_laps_cycle_time'
        endif

        if(c_section .eq. 'in')then
            write(6,*)' input image intensity...'
            read(lun,*)plot_parms%rimage_intensity
            go to 1100
        endif

        ifield_found = 1
        i_overlay = 0

        if(     c_section .eq. 'h' .or. c_section .eq. 'h' 
     1     .or. c_section .eq. '1'
     1     .or. c_section .eq. 'hz'
     1                                                          )then

            if(c_section .eq. 'hz')then
                write(6,101)
 101            format(
     1       '    zoom,density (contours/sfc wind barbs), contour width'       
     1                ,15x,'     ? ',$)
 
                read(lun,*)zoom,density
     1                    ,plot_parms%contour_line_width       
     1                    ,plot_parms%xcen
     1                    ,plot_parms%ycen
     1                    ,plot_parms%zoom_wdw
     1                    ,plot_parms%obs_size

            else
                zoom = 1.0
                density = 1.0
                plot_parms%contour_line_width = 1
            endif

            plot_parms%zoom = zoom

            if(max_radars .ge. 1)then
                l_radars = 1
            else
                l_radars = 0
            endif

            call lapswind_plot(c_display,i4time_ref,lun,nx_l,ny_l,nz_l,
     1                         max_radars,l_radars,r_missing_data,
     1                         laps_cycle_time,zoom,density,
     1                         dyn_low,dyn_high,dx,dy,
     1                         plot_parms,namelist_parms,ifield_found)       

            if(ifield_found .eq. 1)then
                write(6,*)' field found - calling frame'
                call frame
            else
                write(6,*)' field not found - not calling frame'
            endif

        elseif(c_section .eq. 'x' .or. c_section .eq. 'x'
     1                            .or. c_section .eq. 'xz'
     1                            .or. c_section .eq. '2')then
            zoom = 1.0
            plot_parms%zoom = zoom

            if(c_section .eq. 'xz')then
                write(6,102)
 102            format('    density (contours), contour line width'       
     1                ,13x,'     ? ',$)
                read(lun,*)density,plot_parms%contour_line_width
!               plot_parms%contour_line_width = 1
            else
                density = 1.0
                plot_parms%contour_line_width = 1
            endif

            l_atms = .false.

            nxsect = nint(sqrt(float(nx_l**2 + ny_l**2))) ! 541

            write(6,*)' nxsect (nx_p) = ',nxsect

            nx_c = max(nx_l,ny_l)

            nx_t = nint(121. * density)                      

            write(6,*)' nx_c/nx_t= ',nx_c,nx_t

            call xsect(c_display,i4time_ref,lun,l_atms
     1                ,standard_longitude,nx_l,ny_l,nz_l
     1                ,nx_c,nz_l,nxsect,nx_t      
     1                ,r_missing_data,laps_cycle_time,maxstns
     1                ,dyn_low,dyn_high,dx,dy
     1                ,density,plot_parms,namelist_parms,ifield_found)       

            if(ifield_found .eq. 1)then
                write(6,*)' field found - calling frame'
                call frame
            else
                write(6,*)' field not found - not calling frame'
            endif

        elseif(c_section .eq. 's' .or. c_section .eq. 'sz')then
            if(c_section .eq. 'sz')then
                read(lun,*)plot_parms%obs_size           
            endif
            
            call plot_sounding(i4time_ref,lun,nx_l,ny_l,nz_l
     1                        ,r_missing_data,laps_cycle_time,maxstns
     1                        ,i_overlay,plot_parms,namelist_parms
     1                        ,l_plotobs)       

            if(ifield_found .eq. 1)then
                write(6,*)' field found - calling frame'
                call frame
            endif

        elseif(c_section .eq. 'a' .or. c_section .eq. 'az')then
            if(c_section .eq. 'az')then
                write(6,*)' enter plot info (e.g. 180p, 90c)'
                read(lun,*)czoom                   
                if(l_parse(trim(czoom),'b') .eqv. .true.)then
                    l_polar = .true.
                    l_cyl   = .true.
                else
                    l_polar = l_parse(trim(czoom),'p')
                    l_cyl   = l_parse(trim(czoom),'c')
                endif
                if((l_polar .or. l_cyl) .eqv. .false.)then
                    write(6,*)' set l_polar to true'
                    l_polar = .true.
                endif
                write(6,*)' czoom = ',trim(czoom),' l_polar = ',l_polar
     1                   ,' l_cyl = ',l_cyl
!               l_polar = .true.
!               l_cyl = .true.
!               exposure will be controlled by density
                read(lun,*)ipolar_sizeparm
            endif

!           ni_polar = (256 * 2**(ipolar_sizeparm)) - 1
!           nj_polar = (256 * 2**(ipolar_sizeparm)) - 1
            ni_polar = (256 * 2 * (ipolar_sizeparm)) - 1
            nj_polar = (256 * 2 * (ipolar_sizeparm)) - 1

            if(ni_polar .le. 8192)then
                iplo = 1; iphi = ni_polar; jplo = 1; jphi = nj_polar
!               iplo = ni_polar-6000; iphi = ni_polar-200
!               iplo = 3026         ; iphi = ni_polar-3026
!               jplo = nj_polar-1400; jphi = nj_polar-200
            else ! crop polar image
                iplo = 7232         ; iphi = ni_polar-7232
!               jplo = nj_polar-1600; jphi = nj_polar-400
                jplo = 400          ; jphi = 1479
            endif

            write(6,*)' call plot_allsky',ni_polar,nj_polar
     1                                   ,ipolar_sizeparm

            call plot_allsky(i4time_ref,lun,nx_l,ny_l,nz_l
     1                        ,ni_polar,nj_polar,ipolar_sizeparm,density
     1                        ,iplo,iphi,jplo,jphi
     1                        ,r_missing_data,laps_cycle_time
     1                        ,l_polar,l_cyl)       

            if(ifield_found .eq. 1)then
                write(6,*)' field found - calling frame'
                call frame
            endif

        elseif(c_section .eq. 'hi')then ! hsect images
            call hsect_img(i4time_ref,lun,nx_l,ny_l,nz_l
     1                        ,r_missing_data,laps_cycle_time
     1                        ,plot_parms,namelist_parms)

        endif ! c_section

        write(6,901)
901     format(////20x,'   ----   add new frame??    ----')

        goto 1100 ! loop back for user input

999     continue ! quit/end program

        call sflush

        write(6,*)
        write(6,*)' closing gks, look for your ./gmeta file'
!       call clsgks

        if(c_display .ne. 'v')then
!           write(6,*)' saving metacode file in for008.dat'
        endif

        deallocate(dx,dy)

        return

        end
