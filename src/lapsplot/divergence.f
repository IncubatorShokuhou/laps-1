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

        subroutine divergence(uanl,vanl,div,lat,lon,ni,nj
     1                       ,l_grid_north,r_missing_data)

!      ~90            steve albers  original version
!       97-aug-17     ken dritz     added r_missing_data as dummy argument
!       97-aug-17     ken dritz     removed include of lapsparms.for
!       97-oct        steve albers  add lon in call to fflxc.

        include 'trigd.inc'

        real m ! grid points per meter

        data scale/1./

        real lat(ni,nj),lon(ni,nj)
        real uanl(ni,nj),vanl(ni,nj)
        real uanl_grid(ni,nj),vanl_grid(ni,nj)
        real div(ni,nj)
        real one(ni,nj)

        real dum1(ni,nj)
        real dum2(ni,nj)
        real dum3(ni,nj)

        character*6 c6_maproj

        logical l_grid_north

!       grid_spacing_m = sqrt(
!    1                 (  lat(1,2) - lat(1,1)                  )**2
!    1               + ( (lon(1,2) - lon(1,1))*cosd(lat(1,1))  )**2
!    1                                  )    * 111317. ! grid spacing m

        call get_standard_latitudes(slat1,slat2,istatus)
        if(istatus .ne. 1)then
            return
        endif

        call get_grid_spacing(grid_spacing_m,istatus)
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

        do j = 1,nj
        do i = 1,ni
!           uanl(i,j) = 0.
!           vanl(i,j) = 500.0
            one(i,j) = 1.0

            if(l_grid_north)then
                uanl_grid(i,j) = uanl(i,j)
                vanl_grid(i,j) = vanl(i,j)
            else
                call uvtrue_to_uvgrid(uanl(i,j),vanl(i,j)
     1                               ,uanl_grid(i,j),vanl_grid(i,j)
     1                               ,lon(i,j))
            endif

        enddo ! i
        enddo ! j

        call fflxc(ni,nj,m,scale
     1  ,uanl_grid,vanl_grid,one,div,lat,lon
     1  ,dum1,dum2,dum3,r_missing_data)

        do j = 1,nj
        do i = 1,ni

            if(abs(div(i,j)) .gt. 1e10)then
                div(i,j) = 0.
            else
                div(i,j) = -div(i,j)
            endif

        enddo ! j
        enddo ! i

        return
        end


        subroutine vorticity_abs(uanl,vanl,vort,lat,lon,ni,nj,dx,dy
     1                          ,l_grid_north,r_missing_data)

        real lat(ni,nj),lon(ni,nj)                           ! i
        real uanl(ni,nj),vanl(ni,nj)                         ! i
        real coriolis(ni,nj)                                 ! l
        real vort(ni,nj)                                     ! o
        real div(ni,nj)                                      ! l
        real dx(ni,nj), dy(ni,nj)                            ! o

        logical l_grid_north                                 ! i

        data init/0/
        save init

        if(init .eq. 0)then ! call this only once to save time
            write(6,*)'     call get_grid_spacing_array'
            call get_grid_spacing_array(lat,lon,ni,nj,dx,dy)
            init = 1
        endif

        write(6,*)'     call vortdiv'
        call vortdiv(uanl,vanl,vort,div,ni,nj,dx,dy)

        write(6,*)'     call get_coriolis_rotation'
        call get_coriolis_rotation(ni,nj,lat,coriolis)

!       call add(vort,coriolis,vort,ni,nj)
        vort = vort + coriolis 

        return
        end


        subroutine get_coriolis_rotation(nx,ny,lat,coriolis_rotation)
c
        include 'trigd.inc'
        implicit none

        integer nx,ny
        integer i,j

        real  lat(nx,ny)
        real  coriolis_rotation(nx,ny)

        real  omega_ear
        data    omega_ear/7.292e-5/

        do j=1,ny
        do i=1,nx
           coriolis_rotation(i,j)=2*omega_ear*sind(lat(i,j))
        enddo
        enddo

        return
        end



        subroutine calc_potvort(i4time,uanl,vanl,temp_3d,potvort,lat,lon       
     1                         ,ni,nj,nk,nkuv,k_in,l_grid_north,dx,dy
     1                         ,r_missing_data,istatus)       

        include 'constants.inc'

        real lat(ni,nj),lon(ni,nj)                           ! i
        real uanl(ni,nj,nkuv),vanl(ni,nj,nkuv)               ! i
        real temp_3d(ni,nj,nk)                               ! i
        real vort_2d(ni,nj)                                  ! l
        real dx(ni,nj)                                       ! o
        real dy(ni,nj)                                       ! o
        real theta(ni,nj,nk)                                 ! l
        real pres_3d(ni,nj,nk)                               ! l
        real dtheta_dp(ni,nj,nk)                             ! l
        real potvort(ni,nj,nk)                               ! o

        logical l_grid_north                                 ! i

        write(6,*)' subroutine calc_potvort'

!       obtain 3d pressure field
        call get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)
        if(istatus .ne. 1)return

!       calculate theta
        do k = 1,nk
        do i = 1,ni
        do j = 1,nj
            theta(i,j,k) = o(temp_3d(i,j,k),pres_3d(i,j,k))
        enddo ! j
        enddo ! i
        enddo ! k

!       calculate dtheta / dp
        do k = 2,nk-1
            kp1 = k+1
            km1 = k-1
            do i = 1,ni
            do j = 1,nj
                dtheta = theta(i,j,kp1) - theta(i,j,km1)              
                dp     = pres_3d(i,j,kp1) - pres_3d(i,j,km1)
                dtheta_dp(i,j,k) = dtheta / dp
            enddo ! j
            enddo ! i
        enddo ! k
        dtheta_dp(:,:,1)  = dtheta_dp(:,:,2)
        dtheta_dp(:,:,nk) = dtheta_dp(:,:,nk-1)


!       calculate absolute vorticity
        if(k_in .gt. 0)then ! 2d
            kstart = k_in
            kend = k_in
        else                ! 3d where k_in = 0
            kstart = 1
            kend = nk
        endif

        do k = kstart,kend
            if(k_in .gt. 0)then ! 2d
                kuv = 1
            else                ! 3d where k_in = 0
                kuv = k
            endif

            write(6,*)' calling vorticity_abs for kuv = ',kuv
            call vorticity_abs(uanl(1,1,kuv),vanl(1,1,kuv),vort_2d
     1                   ,lat,lon,ni,nj,dx,dy
     1                   ,l_grid_north,r_missing_data)

!           see http://www-das.uwyo.edu/~geerts/cwx/notes/chap12/pot_vort.html

            if(k_in .gt. 0)then     ! just returning 2d field (for efficiency)
                write(6,*)'     calculating 2d potvort - with debugging'       
                do i = 1,ni
                do j = 1,nj
                    potvort(i,j,1) = vort_2d(i,j) * dtheta_dp(i,j,k) 
     1                                            * grav       
                    if(i .eq. ni/2)then ! debug
                        write(6,*)i,j,vort_2d(i,j)
     1                           ,uanl(i,j,1),vanl(i,j,1)
     1                           ,dtheta_dp(i,j,k)
     1                           ,potvort(i,j,1)
                    endif
                enddo ! j
                enddo ! i

            elseif(k_in .eq. 0)then ! return 3d field
                write(6,*)'     calculating 3d potvort'
                potvort(:,:,k) = vort_2d(:,:) * dtheta_dp(:,:,k) 
     1                                        * grav       
            endif
        enddo ! k

        return
        end
