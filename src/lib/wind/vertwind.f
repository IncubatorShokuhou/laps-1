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

        subroutine vert_wind(uanl,vanl                              ! i ! /o
     1          ,u_sfc,v_sfc,ni,nj,nk                               ! i
     1          ,wanl                                               ! o
     1          ,topo,lat,lon,grid_spacing_m                        ! i
     1          ,rk_terrain,r_missing_data,l_grid_north             ! i
     1          ,istatus)                                           ! o

!      ~1990        steve albers  orig version
!       1997 jun    ken dritz     added ni and nj as dummy
!                                 arguments, thus making wsum and one and
!                                 several other arrays automatic.
!       1997 jun    ken dritz     initialize wsum and one dynamically, instead
!                                 of by data statements.
!       1997 jun    ken dritz     added r_missing_data as dummy argument.
!       1997 jun    ken dritz     removed include of 'lapsparms.for'.
!       1997 oct    steve albers  pass lon to fflxc. misc cleanup.
!       1997 dec    steve albers  changed nx_l_max/ny_l_max to ni/nj
!       1998 nov    steve albers  change to m, hopefully cancels out change to
!                                 fflxc/sigma

        include 'trigd.inc'
        real m ! grid points per meter

        real wsum(ni,nj)
        real  one(ni,nj)

        logical l_grid_north ! sfc & 3d winds

        real rk_terrain(ni,nj)
        integer k_terrain(ni,nj)

        data scale/1./

        real uanl(ni,nj,nk),vanl(ni,nj,nk)
        real wanl(ni,nj,nk) ! omega (pascals/second)
        real terrain_w(ni,nj),conv(ni,nj)
        real u_sfc(ni,nj),v_sfc(ni,nj)
        real topo_pa(ni,nj)

        real lat(ni,nj)
        real lon(ni,nj)
        real topo(ni,nj)

        real flu(ni,nj)
        real flv(ni,nj)
        real sigma(ni,nj)

        real u_sfc_grid(ni,nj),
     1         v_sfc_grid(ni,nj)
        real beta_factor(ni,nj)

        real radius_earth
        parameter (radius_earth = 6371e3)

        character*6 c6_maproj

        do i=1,ni
           do j=1,nj
              wsum(i,j) = 0.0
              one(i,j) = 1.0
           enddo
        enddo

        ierrcnt = 0
        istatus = 1

!       grid_spacing_m = sqrt(
!       1                      (  lat(1,2) - lat(1,1)                  )**2
!       1                    + ( (lon(1,2) - lon(1,1))*cosd(lat(1,1))  )**2
!       1                                       )    * 111317. ! grid spacing m

        call get_standard_latitudes(slat1,slat2,istatus)
        if(istatus .ne. 1)then
            return
        endif

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            return
        endif

        if(c6_maproj .eq. 'plrstr')then
            call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                 ,grid_spacing_proj_m)
        else
            grid_spacing_proj_m = grid_spacing_m
        endif

        write(6,*)' grid spacings (m) = ',grid_spacing_m
     1                                   ,grid_spacing_proj_m

        m = 1.0 / grid_spacing_proj_m

        imaxm1 = ni - 1
        jmaxm1 = nj - 1

!       calculate surface pressure array
        do j=1,nj
        do i=1,ni
            topo_pa(i,j) = ztopsa(topo(i,j)) * 100.
        enddo
        enddo

!       rotate sfc winds
        if(.not. l_grid_north)then
          write(6,*)' rotating sfc winds to calculate terrain forcing'
          write(6,*)' calculating beta factor'
          do j=1,nj
          do i=1,ni
            if(u_sfc(i,j) .ne. r_missing_data
     1                   .and. v_sfc(i,j) .ne. r_missing_data)then

                call uvtrue_to_uvgrid(u_sfc(i,j),
     1                              v_sfc(i,j),
     1                              u_sfc_grid(i,j),
     1                              v_sfc_grid(i,j),
     1                              lon(i,j))

            endif

            beta_factor(i,j) = sind(lat(i,j)) / radius_earth

          enddo
          enddo

        else ! winds are already rotated
          do j=1,nj
          do i=1,ni
              u_sfc_grid(i,j) = u_sfc(i,j)
              v_sfc_grid(i,j) = v_sfc(i,j)
              beta_factor(i,j) = 0.
          enddo
          enddo

        endif

        if(c6_maproj .eq. 'latlon')then
            write(6,*)' latlon grid, skipping omega...'
            wanl = r_missing_data
            return
        endif

!       calculate terrain induced omega
        do j=2,jmaxm1
        do i=2,imaxm1

!           units of dterdx are pascals/meter
            dterdx=(topo_pa(i+1,j) - topo_pa(i-1,j)) * .5/grid_spacing_m       
            dterdy=(topo_pa(i,j+1) - topo_pa(i,j-1)) * .5/grid_spacing_m

!           units of ubar,vbar are m/s
            ubar=(u_sfc_grid(i-1,j-1) + u_sfc_grid(i,j-1)) *.5
            vbar=(v_sfc_grid(i-1,j-1) + v_sfc_grid(i-1,j)) *.5

!           units of terrain_w are pascals/second
            terrain_w(i,j) = ubar*dterdx+vbar*dterdy

        enddo ! i
        enddo ! j

!       fill in edges
        do i = 1,ni
            terrain_w(i,   1) = terrain_w(i,     2)
            terrain_w(i,nj)   = terrain_w(i,jmaxm1)
        enddo ! i

        do j = 1,nj
            terrain_w(1,   j) = terrain_w(2,     j)
            terrain_w(ni,j)   = terrain_w(imaxm1,j)
        enddo ! j

        do j = 1,nj
        do i = 1,ni
            k_terrain(i,j) = max( nint(rk_terrain(i,j)) ,1 )
        enddo ! i
        enddo ! j


        write(6,*)'            terrain       conv       omega      beta 
     1corr'

        do k = 1,nk
            call fflxc(ni,nj,m,scale
     1                ,uanl(1,1,k),vanl(1,1,k),one,conv,lat,lon
     1                ,flu,flv,sigma,r_missing_data)

            if(k.gt.1)then
                z_interval=abs(zcoord_of_level(k)-zcoord_of_level(k-1))    
            else
                z_interval=abs(zcoord_of_level(2)-zcoord_of_level(2-1))    
            endif

            do j = 1,nj
            do i = 1,ni

                k_terr = k_terrain(i,j)

                if(k .lt. k_terr)then !
                    wanl(i,j,k) = r_missing_data
!                   uanl(i,j,k) = 0.
!                   vanl(i,j,k) = 0.

                elseif(k .eq. k_terr)then !
                    wsum(i,j)   = terrain_w(i,j)
!    1       - (conv(i,j) + vanl(i,j,k) * beta_factor(i,j)) * z_interval
                    wanl(i,j,k) = wsum(i,j)


                else ! k .gt. k_terr
                    wsum(i,j)   = wsum(i,j)
     1       - (conv(i,j) + vanl(i,j,k) * beta_factor(i,j)) * z_interval       
                    wanl(i,j,k) = wsum(i,j)

                endif

                if(abs(conv(i,j)) .gt. 1.0)then
                    ierrcnt = ierrcnt + 1
                    if(ierrcnt .lt. 20)then
                        write(6,*)' error: large convergence'
     1                          ,i,j,k,k_terrain(i,j),conv(i,j)
                    endif

                    istatus = 0
                endif

                if(j .eq. 29 .and. i .eq. 29)then
                    write(6,111)i,j,k,terrain_w(i,j)
     1              ,conv(i,j),wanl(i,j,k),beta_factor(i,j)*vanl(i,j,k)       
111                 format(3i3,4e12.3)
                endif


            enddo ! i
            enddo ! j

        enddo ! k

        return
        end
