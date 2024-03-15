
        subroutine tower_sfc_driver(maxsta,i4time_sys                     ! i
     1                             ,path_to_tower_data                    ! i
     1                             ,lat,lon,ni,nj,grid_spacing            ! i
     1                             ,laps_cycle_time                       ! i
     1                             ,itime_before,itime_after              ! i
     1                             ,nn,n_local_g,n_local_b,stations       ! i/o
     1                             ,store_1,store_2,store_2ea             ! o
     1                             ,store_3,store_3ea,store_4,store_4ea   ! o    
     1                             ,store_5,store_5ea,store_6,store_6ea   ! o
     1                             ,store_7,store_cldht,store_cldamt      ! o
     1                             ,provider,istatus)                     ! o
c        
        integer ni, nj, maxsta, maxobs 
        integer maxlvls ! raw/processed stations for snd file
        parameter (maxlvls=10)
c
	real    lat(ni,nj), lon(ni,nj), topo(ni,nj)
c
        integer    wmoid(maxsta), istatus
        integer    dpchar(maxsta), narg, iargc
c
        character  stations(maxsta)*20, provider(maxsta)*11
        character  reptype(maxsta)*6, atype(maxsta)*6
	character  atime*24
	character  dir_s*256, ext_s*31, units*10, comment*125,var_s*3
	character  filename9*9, filename13*13, a9time_metar_file*9
        character  fname9_to_wfo_fname13*13
	character  data_file_l*150
c
        character*200 path_to_metar
        character*200 path_to_tower_data
        character*8   metar_format
        character*8   a9_to_a8, a8_time

c       declared then used in 'get_local_towerobs' for snd purposes
        real     stalat_s(maxsta,maxlvls),stalon_s(maxsta,maxlvls)
        real     staelev_s(maxsta)
        real     soilmoist_p(maxsta)       
	character  stname_s(maxsta)*5
c
c.....  output arrays.
c
	real  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)
c
c.....	start here.  
c
        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

c.....  get the time from the scheduler or from the user if interactive.
c
        call make_fnam_lp(i4time_sys,filename9,istatus)
c
        write(6,*)' systime = ',filename9
c
c.....  call the routine that reads the mesonet data files, then get the data.
c
        write(6,*)
	write(6,*)'getting mesonet tower data...'
c
        ext_s = 'lso'

        call get_local_towerobs(maxsta,maxlvls,                          ! i
     &                      i4time_sys,lun_out,
     &                      path_to_tower_data,metar_format,
     &                      ext_s,                                       ! i
     &                      itime_before,itime_after,
!    &                      grid_east,grid_west,grid_north,grid_south,
     &                      lat,lon,ni,nj,                               ! i
     &                      nobs,                                        ! o
!    &                      stations,
!    &                      reptype,atype,wmoid,
!    &                      laps_cycle_time, 
     &                      stalat_s,stalon_s,staelev_s,                 ! o
     &                      stname_s,                                    ! o
     &                      soilmoist_p,                                 ! o
     &                      istatus)

	if(istatus .ne. 1) then
	   print *, ' istatus=0 returned from get_local_towerobs...'
	   return
	endif

        if(nobs .gt. maxsta)then
           write(6,*)' error: nobs > maxsta ',nobs,maxsta
           return
        endif
c
!       final qc check 
        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

!       check for no obs
        if(nobs .eq. 0)then
            write(6,*)' note: no snd appended due to no tower obs'
            return
        endif
c
c
c.....	that's about it...let's go home.
c
	write(6,*)' normal completion of tower_driver, nobs = ',nobs

        return
	end


 


