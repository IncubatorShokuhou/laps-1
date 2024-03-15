       subroutine getlapsvxx(imax,jmax,kmax,maxradar,c_radar_id,         ! i
     &      n_radars,c_extension_proc,i4timefile_proc,i4_tol,rheight_3d, ! i   
     &      lat,lon,topo,i4_file_closest,                                ! i
     &      nx_r,ny_r,igrid_r,                                           ! i
     &      rlat_radar,rlon_radar,rheight_radar,n_valid_radars,          ! o
     &      grid_ra_ref,maxradarg,                                       ! o
     &      grid_ra_ref_offset,ioffset,joffset,maxradaro,                ! o
     &      l_offset,                                                    ! i
     &      istatus)                                                     ! o
c
       integer       imax,jmax,kmax  !same as imax,jmax,kmax in lapsparms.for
       integer       maxradar
       integer       maxfiles
       integer       i,j,k
       integer       lvl_3d(kmax)
       integer       i4time
       integer       i4_tol
       integer       n_radars
       integer       n_ref_grids
       integer       istatus
       integer       istatus_2dref
       integer       istatus_3dref
       integer       level

       integer       i4timefile_proc
       integer       i4_file_closest(n_radars)
       integer       ioffset(maxradar)
       integer       joffset(maxradar)

!      note that only one of maxradarg, maxradaro is non-zero
       real          grid_ra_ref(imax,jmax,kmax,maxradarg)
       real          grid_ra_ref_offset(nx_r,ny_r,kmax,maxradaro)
       real          grid_ra_ref_3d(imax,jmax,kmax)             
       real          rheight_3d (imax,jmax,kmax)
       real          lat(imax,jmax)
       real          lon(imax,jmax)
       real          topo(imax,jmax)
       real          rlat_radar(maxradar)
       real          rlon_radar(maxradar)
       real          rheight_radar(maxradar)
       real          r_missing_data

       real          zcoord_of_level
c
c readlaps stuff
c
       character ext*31, var_3d(kmax)*3,
     &lvl_coord_3d(kmax)*4,units_3d(kmax)*10,
     &comment_3d(kmax)*125

       character       c_extension_proc(maxradar)*4
       character*31    ext_a(maxradar)
       character*4     c_radar_id(maxradar)
       character*1     c_again

       logical l_apply_map
       logical l_low_fill
       logical l_high_fill
       logical l_offset

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)return

c ========================================================
c read intermediate radar files. this will return the array only with those
c valid radars within the internal 20min time window of 'read_radar_3dref'.
c
      l_apply_map=.true.
      l_low_fill = .true.
      l_high_fill= .true.

      write(6,*)
      write(6,*)'get_laps_vxx: reading v-file reflectivity, ',
     1          '# of potential radars = ',n_radars

!     number of valid radars - change this index to 'l'?
      k = 0

      do kcount=1,n_radars

         write(6,*)

         i4_elapsed = ishow_timer()

         ext = trim(c_extension_proc(kcount))

         i4_diff = i4timefile_proc - i4_file_closest(kcount)

         k = k+1      ! used for output arrays
         write(6,*)' radar #, i4_diff = ',kcount,k,i4_diff

         i4_tol_radar = 0 

         call read_radar_3dref(i4_file_closest(kcount),               ! i
     1   i4_tol_radar,i4_ret,                                         ! i/o
     1   l_apply_map,r_missing_data,
     1   imax,jmax,kmax,ext,                                          ! i
     1   lat,lon,topo,l_low_fill,l_high_fill,
     1   rheight_3d,
     1   grid_ra_ref_3d,                                              ! o
     1   rlat_radar(k),rlon_radar(k),rheight_radar(k),c_radar_id(k),  ! o
     1   n_ref_grids,istatus_2dref,istatus_3dref)                     ! o

c check laps analysis values
         if(istatus_3dref.ne.1 .and. istatus_3dref.ne.-1)then
            write(6,*)'error: unsuccessful reading radar ',kcount,k,ext       
            k=k-1
         else
            write(6,*)'successful reading radar ',kcount,k,ext
            i4_elapsed = ishow_timer()
            write(6,*)'radar lat/lon/elev: ',rlat_radar(k),
     &                         rlon_radar(k),rheight_radar(k)
            level=9
            write(6,*)
            write(6,*)'level ',level,' analysis output'
            write(6,*)'------------------------------'
29          format(1x,'  i  j    lat   lon     topo    ref ')
            write(6,29)

!           move radar data to second array via offset
            if(l_offset)then
              if(rlat_radar(k) .eq. r_missing_data .or.
     1           rlon_radar(k) .eq. r_missing_data      )then
                write(6,*)' no valid or single lat/lon for radar ',k
!    1                   ,' ',radar_name(k)
!               l_valid_latlon(k) = .false.

                ioffset(k) = 0
                joffset(k) = 0

              else
                call latlon_to_rlapsgrid(rlat_radar(k),
     &                                   rlon_radar(k),
     &                                   lat,lon,
     &                                   imax,jmax,
     &                                   ri,rj,
     &                                   jstatus)
                if(jstatus.ne.1)then
                    write(6,*)
     1               'computing ri/rj for radar (outside domain)'    
                endif
!               write(6,*)'name: ',radar_name(k),ri(k),rj(k),k
!               l_valid_latlon(k) = .true.

!               offset is location of lower left corner of small array in the large array
                ioffset(k) = (nint(ri) - igrid_r) - 1
                joffset(k) = (nint(rj) - igrid_r) - 1

              endif

              i4_elapsed = ishow_timer()

              write(6,*)' getlapsvxx - offset info '
     1               ,'ri,rj,ioffset(k),joffset(k),igrid_r : '
     1                ,ri,rj,ioffset(k),joffset(k),igrid_r

              nfill = 0  

              do jo = 1,ny_r
                j = jo + joffset(k)
                if(j .ge. 1 .and. j .le. jmax)then
                  do io = 1,nx_r
                    i = io + ioffset(k)
                    if(i .ge. 1 .and. i .le. imax)then
                      grid_ra_ref_offset(io,jo,:,k) = 
     1                grid_ra_ref_3d(i,j,:)
                      nfill = nfill + 1
                    endif ! in i bounds
                  enddo ! io
                endif ! in j bounds
              enddo ! j

              write(6,*)' nfill = ',nfill

            else ! fill 4d array with 3d array contents
              grid_ra_ref(:,:,:,k) = grid_ra_ref_3d(:,:,:)

            endif ! l_offset

!           write sample of radar data
            do j=1,jmax,20
            do i=1,imax,20
               write(6,30)i,j,lat(i,j),lon(i,j),topo(i,j)
     &                   ,grid_ra_ref_3d(i,j,level)
            end do
            end do

30          format(1x,2i5,1x,3f7.1,1x,2(f8.1,1x))

            write(6,*)'midpoint column ',grid_ra_ref_3d(imax/2,jmax/2,:)      

            i4_elapsed = ishow_timer()

         endif

      enddo

      n_valid_radars = k

      return
      end
