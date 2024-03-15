
      subroutine mosaic_ref_multi(i_ra_count,maxradars,l_low_level,       ! i
     &                            radar_name,lat,lon,nx,ny,nz,            ! i
     &                            rlat_radar,rlon_radar,rheight_radar,    ! i
     &                            topo,rheight_laps,grid_ra_ref,          ! i
     &                            grid_ra_ref_offset,ioffset,joffset,     ! i
     &                            nx_r,ny_r,                              ! i
     &                            maxradarsg,maxradarso,                  ! i
     &                            imosaic_3d,                             ! i
     &                            dist_multiradar_2d,                     ! i  
     &                            l_offset,                               ! i
     &                            grid_mosaic_2dref,grid_mosaic_3dref,    ! i/o
     &                            closest_radar_m,istatus)                ! o
c
c routine mosaics vxx radar files. uses radar closest to the
c laps grid point. depending on the value of l_low_level, uses 
c either the maximum reflectivity in the column (l_low_level = false)
c or uses the reflectivity on the first level above the laps elevation
c (l_low_level = true).
c
      real    lat(nx,ny)
      real    lon(nx,ny)
      real    grid_ra_ref(nx,ny,nz,maxradarsg)
      real    grid_ra_ref_offset(nx_r,ny_r,nz,maxradarso)
      real    grid_mosaic_2dref(nx,ny)
      real    dist_multiradar_2d(nx,ny,maxradars)
      real    closest_radar_m(nx,ny)
      real    grid_mosaic_3dref(nx,ny,nz)
      real    topo(nx,ny)
      real    rheight_laps(nx,ny,nz)
      real    rlat_radar(maxradars)
      real    rlon_radar(maxradars)
      real    rheight_radar(maxradars)
      real    ri(maxradars)
      real    rj(maxradars)
      real    dist_radar_m(maxradars)

      logical   l_low_level
      logical   found_height
      logical   l_valid
      logical   l_valid_latlon(maxradars)
      logical   l_offset 
!     parameter (l_offset = .true.)

      integer   lr_2d(nx,ny)                     ! closest radar
      integer   ioffset(maxradars),joffset(maxradars)

      character*4 radar_name(maxradars)

      write(6,*)
      write(6,*)' subroutine mosaic_ref_multi: imosaic_3d =', imosaic_3d       

!     if(l_offset)then
!         write(6,*)' grid_ra_ref_offset range: '
!    1              ,minval(grid_ra_ref_offset)
!    1              ,maxval(grid_ra_ref_offset)
!     else
!         write(6,*)' grid_ra_ref range: '
!    1              ,minval(grid_ra_ref)
!    1              ,maxval(grid_ra_ref)
!     endif

      call get_ref_base(ref_base, istatus)
      call get_r_missing_data(r_missing_data, istatus)
      call get_grid_spacing_cen(grid_spacing_cen_m,istatus)

!     initialize
      grid_mosaic_2dref = r_missing_data
      grid_mosaic_3dref = r_missing_data
      lr_2d = 0
c
c first find the radar location in ri/rj laps-space.
c
      istatus = 1
      if(i_ra_count .ge. 1)then ! essentially all the time

!        determine which radars have valid lat/lon
         do k = 1,i_ra_count 
            if(rlat_radar(k) .eq. r_missing_data .or.
     1         rlon_radar(k) .eq. r_missing_data      )then
                write(6,*)' no valid or single lat/lon for radar ',k
     1                   ,' ',radar_name(k)
                l_valid_latlon(k) = .false.

            else
                call latlon_to_rlapsgrid(rlat_radar(k),
     &                                   rlon_radar(k),
     &                                   lat,lon,
     &                                   nx,ny,
     &                                   ri(k),rj(k),
     &                                   jstatus)
                if(jstatus.ne.1)then
                    write(6,*)
     1               'error computing ri/rj for radar (outside domain)?'    
                endif
                write(6,*)radar_name(k),k,ri(k),rj(k)
     1                   ,rlat_radar(k),rlon_radar(k)
51              format('valid lat/lon: ',a,i5,2f8.1,2f8.2)
                l_valid_latlon(k) = .true.

                if(l_offset)then
                    write(6,*)'      offsets ',ioffset(k),joffset(k)               
                endif

            endif

         enddo ! radars

!        loop through all horizontal gridpoints to define best radar array
         icntn=0
         icntp=0
         icntz=0
         icntb=0
         do j = 1,ny
         do i = 1,nx

            r_min_dist_m = sqrt(float(nx*nx+ny*ny))*grid_spacing_cen_m       
c
c find the valid radar with the minimum distance to the grid point in question.
c
            lr = 0 

!           loop through all radars
            do l = 1,i_ra_count
               l_valid = .false.

               if(.not. l_offset)then
                 do k = 1,nz ! look for non-missing reflectivity in column
                   if(grid_ra_ref(i,j,k,l) .ne. r_missing_data)then     
                     l_valid = .true.
                   endif
                 enddo ! k

               else ! use offset array
                 io = i - ioffset(l)
                 jo = j - joffset(l)
                 if(io .lt. 1 .or. io .gt. nx_r .or. 
     1              jo .lt. 1 .or. jo .gt. ny_r)then
!                  write(6,*)' offset out of domain ',i,j,io,jo
!    1                      ,ioffset(l),joffset(l),radar_name(l)
!                  stop
                   continue

                 else   
                   do k = 1,nz ! look for non-missing reflectivity in column
                     if(grid_ra_ref_offset(io,jo,k,l) .ne. 
     1                                  r_missing_data)then     
                       l_valid = .true.
                       goto 50
                     endif
                   enddo ! k

                 endif ! in bounds of offset array                  

               endif ! use full array (not offset)

 50            if(l_valid_latlon(l))then
                   ridist = float(i)-ri(l)
                   rjdist = float(j)-rj(l)
                   rijdist=sqrt(ridist*ridist + rjdist*rjdist)
                   dist_radar_m(l) = rijdist * grid_spacing_cen_m

               else ! no valid or single lat/lon, use distance array
                   dist_radar_m(l) = dist_multiradar_2d(i,j,l)

               endif

               if(dist_radar_m(l) .lt. r_min_dist_m .and. l_valid)then       
                  lr=l
                  r_min_dist_m = dist_radar_m(l)
                  closest_radar_m(i,j) = r_min_dist_m       
               endif

            enddo ! l

            lr_2d(i,j) = lr

         enddo ! i
         enddo ! j

         do l = 1,i_ra_count ! mosaic in one radar at a time
            do j = 1,ny
            do i = 1,nx

            if(l .eq. lr_2d(i,j))then ! closest radar at this grid point 
               r_dbzmax=ref_base

               if(.not. l_offset)then
!                 get max ref in column
                  do k=1,nz
                     if(grid_ra_ref(i,j,k,l).ne.ref_base .and.
     1                  grid_ra_ref(i,j,k,l).ne.r_missing_data)then
                        r_dbzmax=max(r_dbzmax,grid_ra_ref(i,j,k,l))
                     endif
                  enddo

                  if(imosaic_3d .eq. 0)then      ! vrc output only
                     do k=1,nz 
                        grid_mosaic_2dref(i,j)=r_dbzmax 
                        grid_mosaic_3dref(i,j,k)=r_dbzmax 
                     enddo ! k 

                  elseif(imosaic_3d .eq. 1)then  ! vrz output only
                     grid_mosaic_3dref(i,j,:)=grid_ra_ref(i,j,:,l)   

                  elseif(imosaic_3d .eq. 2)then  ! both vrc & vrz
                     do k=1,nz 
                        grid_mosaic_2dref(i,j)=r_dbzmax 
                        grid_mosaic_3dref(i,j,k)=grid_ra_ref(i,j,k,l)   
                     enddo ! k 

                  endif ! imosaic_3d

               else ! l_offset
                 io = i - ioffset(l)
                 jo = j - joffset(l)
                 
                 if(io .ge. 1 .and. io .le. nx_r .and.
     1              jo .ge. 1 .and. jo .le. ny_r)then

!                  get max ref in column
                   do k=1,nz
                     if(grid_ra_ref_offset(io,jo,k,l).ne.ref_base .and.
     1                  grid_ra_ref_offset(io,jo,k,l).ne.r_missing_data
     1                                                             )then
                        r_dbzmax=max(r_dbzmax,
     1                               grid_ra_ref_offset(io,jo,k,l))
                     endif
                   enddo

                   if(imosaic_3d .eq. 0)then      ! vrc output only
                     grid_mosaic_2dref(i,j)=r_dbzmax 
                     grid_mosaic_3dref(i,j,:)=r_dbzmax 

                   elseif(imosaic_3d .eq. 1)then  ! vrz output only
                     grid_mosaic_3dref(i,j,:)=
     1               grid_ra_ref_offset(io,jo,:,l)       

                   elseif(imosaic_3d .eq. 2)then  ! both vrc & vrz
                     grid_mosaic_2dref(i,j)=r_dbzmax 
                     grid_mosaic_3dref(i,j,:)=
     1               grid_ra_ref_offset(io,jo,:,l)       

                   endif ! imosaic_3d

                 endif ! in bounds of offset array

               endif ! .true.

!              increment stats
               if(r_dbzmax.ne.ref_base)then
                  if(r_dbzmax.lt.0.0)then
                     icntn=icntn+1
                  elseif(r_dbzmax.eq.0.0)then
                     icntz=icntz+1
                  else
                     icntp=icntp+1
                  endif
               else
                  icntb = icntb + 1
               endif

            endif ! .true.

            enddo ! i
            enddo ! j
         enddo ! l

      endif ! i_ra_count > 1

      print*,'statistics for this mosaic'
      print*,'--------------------------'
      print*,'num points > 0.0  ',icntp
      print*,'num points = 0.0  ',icntz
      print*,'num points < 0.0  ',icntn
      print*,'num points = base ',icntb

      intvl = int(nx/80) + 1

      write(6,*)
      write(6,*)' map of ',i_ra_count,' valid radars used:'

      do j = ny,1,-intvl
          write(6,101)(lr_2d(i,j),i=1,nx,intvl) 
 101      format(75i2)
      enddo ! j

      return
      end
