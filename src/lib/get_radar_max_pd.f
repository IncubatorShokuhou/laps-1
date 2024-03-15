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

        subroutine get_radar_max_pd(i4time_beg,i4time_end,imax,jmax,kmax
     1                             ,heights_3d,ext_radar
     1                             ,max_radar_files                      ! i
     1                             ,lat,lon,topo
     1                             ,radar_max,frac_sum,istatus)

!       1992         steve albers
!       1996 feb     steve albers  call read_radar_2dref for radar data

!       input
        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)
        real heights_3d(imax,jmax,kmax)

!       output
        real radar_max(imax,jmax) ! m

        real dbz_2d(imax,jmax)

        character*9 asc_tim_9,asc_tim_9_beg,asc_tim_9_end
        integer i4time_file(max_radar_files)
        real frac(max_radar_files)
        character c_fnames(max_radar_files)*80

        character*255 c255_radar_filename
        character*255 c_filespec

        real grid_ra_ref(imax,jmax,kmax)

        character*3 ext_radar
        character*150  directory

        character*4 radar_name

        integer max_radar_gap
        parameter (max_radar_gap = 4800)

        data mode_radar/1/

        istatus = 1

c       rmax_so_far = 0.

        call make_fnam_lp(i4time_beg,asc_tim_9_beg,istatus)
        call make_fnam_lp(i4time_end,asc_tim_9_end,istatus)

        write(6,*)' max radar reflectivity from ',asc_tim_9_beg,
     1                         ' to ',asc_tim_9_end

!       get file times
        call get_directory(ext_radar,directory,len_dir)
        c_filespec = directory(1:len_dir)//'*.'//ext_radar

        call get_file_names(c_filespec,
     1                      i_nbr_files_ret,
     1                      c_fnames,
     1                      max_radar_files,
     1                      i_status )

        min_diff = 1999999999

        if(i_nbr_files_ret .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
        else ! error condition
            write(6,*)' no radar data available'
            istatus = 0
            return
        endif

10      do i=1,i_nbr_files_ret
            asc_tim_9 = c_fnames(i)(lenf+1:lenf+9)
            call i4time_fname_lp(asc_tim_9,i4time_file(i),istatus)
        enddo

        call get_fracs_radar(i4time_beg,i4time_end,max_radar_gap
     1                      ,max_radar_files
     1                      ,i_nbr_files_ret
     1                      ,i4time_file,frac,frac_sum,istatus)

        if(istatus .ne. 1)then
            return
        endif

!       initialize
        do j = 1,jmax
        do i = 1,imax
            radar_max(i,j) = 0.
        enddo ! i
        enddo ! j

        i4_interval = (i4time_end - i4time_beg)
        i4time_mid = i4time_beg + i4_interval / 2
        i4time_temp = 0

        do ifile = 1,i_nbr_files_ret
          if(frac(ifile) .gt. 0.)then

            i4time_radar = i4time_file(ifile)

            call make_fnam_lp(i4time_radar,asc_tim_9,istatus)

!           write(6,101)asc_tim_9,frac(ifile)
!101        format(' time, frac = ',a9,2x,f6.3)

            c255_radar_filename
     1      = c_fnames(ifile)(1:lenf)//asc_tim_9//'.'//ext_radar

            write(6,*)

            if(ext_radar .ne. 'vrc')then

              call get_ref_base(ref_base,istatus)
              if(istatus .ne. 1)return

              i4_tol = 1200

              call read_radar_3dref(i4time_radar,
     1                 i4_tol,i4_ret,                                   ! i/o
     1                 .true.,ref_base,
     1                 imax,jmax,kmax,ext_radar,
     1                 lat,lon,topo,.false.,.false.,
     1                 heights_3d,
     1                 grid_ra_ref,
     1                 rlat_radar,rlon_radar,rheight_radar,radar_name,     
     1                 n_ref_grids,istatus_2dref,istatus_3dref)       


              if(istatus_2dref .eq. 0)then
                write(6,*)' error in reading radar data'
                istatus = 0
                return
              endif

              call get_max_ref(grid_ra_ref,imax,jmax,kmax,dbz_2d)

            else

              call read_radar_2dref(i4time_radar,radar_name,
     1                  imax,jmax,
     1                  dbz_2d,istatus_2dref)

              if(istatus_2dref .eq. 0)then
                write(6,*)' error in reading radar data'
                istatus = 0
                return
              endif

            endif

            write(6,*)' increment radar data'

            do j = 1,jmax
            do i = 1,imax
                radar_max(i,j) = max(radar_max(i,j),dbz_2d(i,j))
c               rmax_so_far = max(rmax_so_far,dbz_2d(i,j))
            enddo ! i
            enddo ! j

          endif ! frac > 0

c         write(6,*)' cycle to next file',rmax_so_far

        enddo ! i

        istatus = 1

        return
        end

        subroutine get_fracs_radar(i4time_beg,i4time_end,max_radar_gap
     1          ,max_radar_files
     1          ,i_nbr_files_ret,i4time_file,frac,frac_sum,istatus)

        integer max_radar_files

        character*9 asc_tim_9
        real frac(max_radar_files)
        integer i4time_file(max_radar_files)

        i4_interval = i4time_end - i4time_beg

        do i = 1,i_nbr_files_ret
!           write(6,301)i4time_beg,i4time_file(i),i4time_end
!301        format(1x,3i12)
            frac(i) = 0.
        enddo

        do i = 1,i_nbr_files_ret
            if(i4time_file(i) .ge. i4time_beg .and.
     1       i4time_file(i) .le. i4time_end        )then ! within pd
                frac(i) = 1.0

            endif
        enddo ! i

        frac_sum = 0.
        do ifile = 1,i_nbr_files_ret
          if(frac(ifile) .gt. 0.)then
            call make_fnam_lp(i4time_file(ifile),asc_tim_9,istat_fnam)
            write(6,101)asc_tim_9,frac(ifile)
 101        format('        filetime, frac = ',a9,2x,f8.5)
            frac_sum = frac_sum + frac(ifile)
          endif
        enddo

        if(frac_sum .eq. 0.)then
            write(6,*)' no radar data in the time window'
            istatus = 0
        endif

999     if(istatus .ne. 1)then
            write(6,*)' insufficient files within time window'
            return
        endif

        return
        end
