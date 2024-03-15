
        subroutine get_uprad_lyr(ni,nj,gnd_radc,ht,uprad_3d,xcos,ycos)

        use mem_namelist, only: r_missing_data,earth_radius,grid_spacing_m 
        use mem_allsky, only: nc

        include 'trigd.inc'
        include 'rad_nodata.inc'

        real gnd_radc(nc,ni,nj) ! spectral radiance (from sfc lights - wm2srnm)
        real sumrad(nc,ni,nj)
        real uprad_3d(ni,nj,nc) ! spectral upward irradiance for layer (wm2nm)
        real xcos(ni,nj)
        real ycos(ni,nj)
        real, allocatable :: drad(:,:)
        real, allocatable :: aef(:,:)
        real, allocatable :: disti_a(:,:)
        real, allocatable :: distj_a(:,:)

        if(.false.)then
          radius = 120000.
          iradius = nint(radius / grid_spacing_m)
          write(6,*)' ht/radius/iradius = ',ht,radius,iradius
        else
          iradius = max(ni,nj)
        endif

        allocate(drad(-iradius:+iradius,-iradius:+iradius))
        allocate(aef(-iradius:+iradius,-iradius:+iradius))
        allocate(disti_a(-iradius:+iradius,-iradius:+iradius))
        allocate(distj_a(-iradius:+iradius,-iradius:+iradius))

!       determine radiation weighting function array
        do ii = -iradius,+iradius
        do jj = -iradius,+iradius
            disti = float(ii) * grid_spacing_m
            distj = float(jj) * grid_spacing_m
            distr = sqrt(disti**2+distj**2+ht**2)
            sin_theta_r = ht/distr
            stearadians = sin_theta_r * (grid_spacing_m / distr)**2
            drad(ii,jj) = stearadians
            aef(ii,jj) = aef_f(asind(sin_theta_r))
            disti_a(ii,jj) = disti
            distj_a(ii,jj) = distj
!           if(ii .eq. 0)then
!               write(6,*)'jj|disti|distj|ht|distr|stearadians',jj,disti,distj,ht,distr,stearadians
!           endif
        enddo ! ii
        enddo ! jj

        gndmax = maxval(gnd_radc(1,:,:))

        write(6,*)' range of gnd_radc (red: wm2nm) = ',minval(gnd_radc(1,:,:)),gndmax
        write(6,*)' range of drad = ',minval(drad),maxval(drad)
        write(6,*)' sum of drad = ',sum(drad)
        sumrad = 0. ! initialize
        uprad_3d(:,:,:) = 0.

        if(gndmax .eq. 0.)then
           write(6,*)' skip rest of get_uprad'
           goto 999
        endif

        if(iradius .ge. 100)then
          iskip = 10
        elseif(iradius .ge. 50)then
          iskip = 3
        elseif(iradius .ge. 20)then
          iskip = 2
        else
          iskip = 1
        endif

        write(6,*)' iskip = ',iskip

        if(.false.)then

!        loop through each gridpoint on the layer array
         do i = 1,ni,iskip
         do j = 1,nj,iskip

!         index limits of layer array
          imin = max(i-iradius,1)
          imax = min(i+iradius,ni)
          jmin = max(j-iradius,1)
          jmax = min(j+iradius,nj)

!         index limits on the weighting array
          iimin = imin-i
          iimax = imax-i
          jjmin = jmin-j
          jjmax = jmax-j
          do ic = 1,nc
            uprad_3d(i,j,ic) = sum( drad(iimin:iimax,jjmin:jjmax) &
                                  * gnd_radc(ic,imin:imax,jmin:jmax) &
                                  * aef(iimin:iimax,jjmin:jjmax) )
          enddo ! ic

          if(.true.)then
         
            ic = 2
            xsum             = sum( drad(iimin:iimax,jjmin:jjmax) &
                                  * gnd_radc(ic,imin:imax,jmin:jmax) &
                                  * aef(iimin:iimax,jjmin:jjmax) &
!                                 * 2000.0)
                                  * disti_a(iimin:iimax,jjmin:jjmax))

            ysum             = sum( drad(iimin:iimax,jjmin:jjmax) &
                                  * gnd_radc(ic,imin:imax,jmin:jmax) &
                                  * aef(iimin:iimax,jjmin:jjmax) &
!                                 * 2000.0)
                                  * distj_a(iimin:iimax,jjmin:jjmax))

            npts = (iimax - iimin + 1) * (jjmax - jjmin + 1)

            if(uprad_3d(i,j,ic) .gt. 0.)then
              xave = xsum / uprad_3d(i,j,ic) 
              yave = ysum / uprad_3d(i,j,ic) 
            else
              xave = 0.
              yave = 0.             
            endif
            
            vecmag = sqrt(xave**2 + yave**2 + ht**2)

            if(vecmag .le. 0.)then
               xcos(i,j) = r_missing_data
               ycos(i,j) = r_missing_data
            else
               xcos(i,j) = xave / vecmag
               ycos(i,j) = yave / vecmag
            endif

            if(i .eq. 1 .and. j .eq. 1)then
                write(6,*)'xyh ',i,j,xsum,ysum,xave,yave,ht
                if(xcos(i,j) .eq. r_missing_data)then
                   write(6,*)' warning in get_uprad_lyr: xcos/ycos has missing value'
                endif
            endif

          endif
         enddo ! j
         enddo ! i

!        interpolate to fill in the horizontal
!        note that the integrated radiance over the hemisphere gives
!        irradiance. we can multiply that by the extinction
!        coefficient to obtain the emission density "s", equivalent to
!        radiant flux per unit volume.
         do ic = 1,nc
           write(6,*)' fill horizontal layer to yield uprad (wm2nm) for color',ic
           write(6,*)' range of uprad ',minval(uprad_3d(:,:,ic)),maxval(uprad_3d(:,:,ic))
           call bilinear_fill(uprad_3d(:,:,ic),ni,nj,iskip,r_missing_data)
           write(6,*)' range of uprad ',minval(uprad_3d(:,:,ic)),maxval(uprad_3d(:,:,ic))
         enddo ! ic
         call bilinear_fill(xcos(:,:),ni,nj,iskip,r_missing_data)
         call bilinear_fill(ycos(:,:),ni,nj,iskip,r_missing_data)

        else ! alternate looping strategy (under construction)

!        loop through each gridpoint on the layer array
         radius_max = 0.
         do i = 1,ni
         do j = 1,nj
           if(gnd_radc(1,i,j) .gt. 0.)then
             radius = (gnd_radc(1,i,j)**0.4) * 3e7
             radius_max = max(radius,radius_max) 
             iradius = nint(radius / grid_spacing_m)

             imin = max(i-iradius,1)
             imax = min(i+iradius,ni)
             jmin = max(j-iradius,1)
             jmax = min(j+iradius,nj)

!            surrounding area loop of layer array
             do ii = imin,imax
             do jj = jmin,jmax
               idif = ii-i
               jdif = jj-j 
!              do ic = 1,nc
                 uprad_3d(ii,jj,:) = uprad_3d(ii,jj,:) + drad(idif,jdif) * gnd_radc(:,i,j) * aef(idif,jdif)
!              enddo ! ic
             enddo ! jj
             enddo ! ii

           endif
         enddo ! j
         enddo ! i

         write(6,*)' radius_max = ',radius_max

        endif

999     deallocate(drad)
        deallocate(aef)
        deallocate(disti_a)
        deallocate(distj_a)

        return
        end
