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

      program gen_verif_raob

      implicit none

      include 'grid_fname.cmn'

      integer           max_verif
      parameter        (max_verif=4)

      character*256     model_dir	!location of model data directories
                                        !lapsprd, or location of fua, fsf
      integer		i4time		!i4time of laps/model file to read
      character*9	a9_time
      character*150     verif_output_dir(max_verif) ! -->
                                        !type of verif  desired ... builds output directory
                                        !names (nobal, bal, bkgd, grid)
      character*256	nl_dir		!directory where verify_raob.nl located
      integer           ni, nj, nk	!i, j and k grid dimensions
      real 		stdlon		!standard longitude
      integer		istatus		!return value from subroutines
      integer		len, balance,i
      real		r_missing_data
      character*2       cur_hr
      integer		i4_raob_window 
      real		max_ht_m_proc     !maximum height(m) to process up to

      character*1       type_obs
      character*150     path_to_raw_profiler
      character*150     path_to_raw_sounding
      integer           raob_process_lag_bal
      integer           raob_process_lag
      integer           n_verif
      real              verif_missing_data
c
c     begin
c

c     max height to process up to in meters
      max_ht_m_proc = 25000.0

c     do laps housekeeping

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
         write(6,*)' error in get_laps_config'
         goto 999
      endif

      call get_grid_dim_xy(ni,nj,istatus)

      if (istatus .ne. 1) then
         write (6,*) 'error getting horizontal domain dimensions'
         go to 999
      endif

      call get_laps_dimensions(nk,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'error getting vertical domain dimension'
         go to 999
      endif

c     get stdlon
      call get_standard_longitude(stdlon,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'error getting standard longitude'
         go to 999
      endif

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1) then
        write(6,*) 'error gettin r_missing_data value'
        goto 999
      endif

c     get current time from systime.dat
      call get_systime(i4time,a9_time,istatus)

      if(istatus .ne. 1) then
        write(6,*) 'error from get_systime'
        go to 999
      endif

c debug lw
c     i4time = 1346508000
c     a9_time = '022441400'

      call read_verif_nl(type_obs,path_to_raw_profiler,
     1 path_to_raw_sounding, raob_process_lag,raob_process_lag_bal,
     1                   max_verif, verif_output_dir,
     1                   verif_missing_data, n_verif, istatus)
      if (istatus .ne. 1) then
        write(6,*)' error in read_verif_nl '
        istatus = 0
        stop
      endif

      

c     get the directory for model_dir (lapsprd directory)
c     get the directory for nl_dir (static directory)

      call get_directory('lw3',model_dir,len)
      model_dir = model_dir(1:len-4)
      call get_directory(grid_fnam_common,nl_dir,len)

      do i=1,n_verif

         call downcase(verif_output_dir(i),verif_output_dir(i))
         if(verif_output_dir(i)(1:len).eq.'bal')then
            balance=1
         elseif(verif_output_dir(i)(1:len).eq.'nobal')then
            balance=0
         elseif(verif_output_dir(i)(1:len).eq.'bkgd')then
            balance=2
         endif
         

         call gen_verif_raob_sub(model_dir, i4time, 
     1                        nl_dir, ni, nj, nk, stdlon,
     1                        i, balance, r_missing_data, 
     1                        max_ht_m_proc, max_verif,
     1                        istatus)

         if(istatus .ne. 1) then
             write (6,*) 'error in gen_verif_raob_sub'
         endif

      enddo

999   continue

      end
