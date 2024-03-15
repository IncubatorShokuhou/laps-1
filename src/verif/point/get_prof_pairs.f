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

      subroutine get_prof_pairs(prof_fname, model_dir, i4time, 
     1                          output_fname, nl_dir, ni, nj,
     1                          nk, lats, lons, stdlon, 
     1                          laps_levels_mb,laps_levels_pa,
     1                          balance, r_missing_data, 
     1                          verif_missing_data, istatus)

      implicit none

      character*(*)     prof_fname	!path and name of profiler file to read
      character*(*)     model_dir	!location of model data directories
                                        !lapsprd, or location of fua, fsf
      integer		i4time		!i4time of laps/model file to read
      character*(*)     output_fname	!path and name of output file
      character*(*)	nl_dir		!directory where verify_prof.nl located
      integer           ni, nj, nk	!i, j and k grid dimensions
      real		lats(ni,nj) 	!domain lats
      real		lons(ni,nj)	!domain lons
      real 		stdlon		!standard longitude
      real            laps_levels_mb(nk) !laps pressure levels
      real            laps_levels_pa(nk) !laps pressure levels
      integer		balance
      real		r_missing_data
      real		verif_missing_data
      integer		istatus		!return value from subroutine

      integer		max_profs
      parameter         (max_profs=50)
      integer           numhts
      parameter         (numhts=64)
      integer		max_var
      parameter         (max_var=150)

      integer 		nf_status, ncid, varid
      integer           all_profs	! 1=use all profilers in domain
      integer           use_prof(max_profs), n_profs_use
      integer		kmax, kdim
      integer           lun, i, j, lend
      character*256     filename
      logical           l_eof
      integer 		wmonum_use(max_profs)
      character*6       staname_use(max_profs)

      integer           nprofs, nhts(max_profs)
      integer		wmostanum(max_profs) 
      character*6       staname(max_profs) 
      real            stalat(max_profs), stalon(max_profs), 
     1                  staelev(max_profs), 
     1                  ri(numhts,max_profs),
     1                  rj(numhts,max_profs), 
     1                  rk(numhts,max_profs)
      character*9	a9_time(max_profs)
      real            htmsl(numhts, max_profs), 
     1                  up(numhts,max_profs), 
     1                  vp(numhts,max_profs),
     1                  wp(numhts,max_profs), 
     1                  baselevels(numhts)
      real            tsp(max_profs), rhsp(max_profs),
     1                  usp(max_profs), vsp(max_profs)

      real            ui(numhts,max_profs), 
     1                  vi(numhts,max_profs),
     1                  wi(numhts,max_profs), 
     1                  tsi(max_profs), rhsi(max_profs),
     1                  usi(max_profs), vsi(max_profs)

      real            htm(ni,nj,nk) 
      real            tm(ni,nj,nk)
      real            um(ni,nj,nk)
      real            vm(ni,nj,nk)
      real            omm(ni,nj,nk)
      real		tsm(ni,nj)
      real		rhsm(ni,nj)
      real		usm(ni,nj)
      real		vsm(ni,nj)
c
c     begin
c

      istatus = 1   !assume good return

c     set up nhts(max_profs)
      do i = 1, max_profs
        nhts(i) = numhts
      enddo

c     read verify_prof.nl

      lun = 10
      call s_len(nl_dir,lend)
      filename = nl_dir(1:lend)//'/verif_prof.nl'
      open(lun,file=filename,status='old',err=900)

      l_eof = .false.
      n_profs_use = 0
      i = 1

50    format(i5,1x,a6)

      do while (.not.l_eof)
        read(lun,50,end=55,err=901)wmonum_use(i), staname_use(i)
        i = i + 1
        goto 56
55      l_eof = .true.
56      continue
      enddo

      close(lun)
      n_profs_use = i - 1

      all_profs = 0	!set to 1 if use all profs in file

c     see if wmonum_use(1) .eq. -1...if so, use all profilers in file
      if (wmonum_use(1) .eq. -1) then
        all_profs = 1
      endif

c     read profiler file "prof_fname" and return requested profiler data
      call get_prof_data(prof_fname, max_profs, numhts, 
     1                   nprofs, a9_time,
     1                   stdlon, wmostanum, staname, 
     1                   stalat, stalon, staelev,
     1                   htmsl, tsp, rhsp, usp,
     1                   vsp, up, vp, wp, 
     1                   verif_missing_data, istatus)
       
      if (istatus .ne. 1) then
        write(6,*) ' unable to read profiler file: ',prof_fname
        istatus = 0
        return
      endif

c     set which profilers to pull laps/model data from
      if (all_profs .ne. 1) then   !set use_prof(i) to 1 if profiler is on list
        do i = 1, n_profs_use
          do j = 1, nprofs
            if ((wmonum_use(i) .eq. wmostanum(j)) .or.
     1           (staname_use(i) .eq. staname(j))) then
              use_prof(j) = 1    !set to 1 if current profiler is on list to use
            else
              use_prof(j) = 0  
            endif
          enddo
        enddo
      else
        do i = 1, nprofs
          use_prof(i) = 1    !set to 1 if current profiler is on list to use
        enddo
      endif

c     get laps/model data
      call get_model_data(i4time, model_dir, ni, nj, nk, balance, 
     1                    laps_levels_mb, htm, um, vm, omm, tm, 
     1                    tsm, rhsm, usm, vsm, istatus)

      if (istatus .ne. 1) then
        write(6,*) ' unable to read model data '
        istatus = 0
        return
      endif

c     write(6,*) 'usm(121,1)=',usm(121,1)
c     write(6,*) 'vsm(121,1)=',vsm(121,1)
c     write(6,*) 'rhsm(121,1)=',rhsm(121,1)
c     write(6,*) 'tsm(121,1)=',tsm(121,1)
c     write(6,*) 
c     write(6,*) 'usm(70,44)=',usm(70,44)
c     write(6,*) 'vsm(70,44)=',vsm(70,44)
c     write(6,*) 'rhsm(70,44)=',rhsm(70,44)
c     write(6,*) 'tsm(70,44)=',tsm(70,44)
c     write(6,*) 
c     write(6,*) 'usm(41,7)=',usm(41,7)
c     write(6,*) 'vsm(41,7)=',vsm(41,7)
c     write(6,*) 'rhsm(41,7)=',rhsm(41,7)
c     write(6,*) 'tsm(41,7)=',tsm(41,7)
c     write(6,*) 

c     loop through use_prof and pull model data if use_prof(i) .eq. 1
      do i = 1, nprofs
        if (use_prof(i) .eq. 1) then

c         determine ri, rj, rk for ob in the model domain
          call get_rijk_ob_ht(i, max_profs, numhts, 
     1                     laps_levels_pa, nhts,
     1                     stalat(i), stalon(i), lats,
     1                     lons, ni, nj, nk, htmsl, htm, 
     1                     ri, rj, rk, istatus)

          if (istatus .ne. 1) then
            use_prof(i) = 0
            if (istatus .eq. -1) then
            else
              write(6,*) ' unable to determine ri,rj,rk for profiler ',
     1                 wmostanum(i), ' ',staname(i)
            endif
          else
c           pull values interpolated from grid to ob
            call interp_to_ob(i, max_profs, numhts,nhts,
     1                        ni, nj, nk, ri, rj, rk, um, 
     1                        vm, omm, tm, tsm, rhsm, usm,
     1                        vsm, ui, vi, wi, tsi, rhsi,
     1                        usi, vsi, r_missing_data,
     1                        verif_missing_data, istatus) 

            if (istatus .ne. 1) then
              use_prof(i) = 0
              write(6,*) ' unable to interp obs for profiler ',
     1                   wmostanum(i), ' ',staname(i)
            endif

          endif
        endif
      enddo

c     re-calc n_profs_use based on interpolation errors
      n_profs_use = 0
      do i = 1, nprofs
        if (use_prof(i) .eq. 1) n_profs_use = n_profs_use + 1
      enddo
       

c     write output_file

      call write_verif_prof(output_fname, max_profs, numhts,
     1                      nhts, n_profs_use, use_prof, 
     1                      nprofs, wmostanum,
     1                      staname, stalat, stalon, staelev,
     1                      htmsl, ri, rj, rk, up, vp, wp, tsp,
     1                      rhsp, usp, vsp, ui, vi, wi, tsi,
     1                      rhsi, usi, vsi, a9_time, istatus)

      if (istatus .ne. 1) then
        write(6,*) ' error writing profiler verif file: ',
     1              output_fname
      endif

      goto 999

900   print*,'error opening namelist file ', filename
      goto 999

901   print*,'error reading namelist file ', filename
      goto 999


999   continue

      return
      end
!1.........................................................................
      subroutine get_prof_data(infile, max_profs, numhts, 
     1                         nprofs, a9_time,
     1                         stdlon, wmostanum, staname, 
     1                         stalat, stalon, staelev,
     1                         htmsl, tsfc, rhsfc, usfc,
     1                         vsfc, u, v, w, 
     1                         verif_missing_data, istatus)

      implicit none
      include 'netcdf.inc'

      character*256	infile
      integer  		max_profs, numhts, nprofs
      character*9       a9_time(max_profs)
      real		stdlon
      integer  		wmostanum(max_profs) 
      character*6       staname(max_profs) 
      real            stalat(max_profs), stalon(max_profs), 
     1                  staelev(max_profs), pi
      real            htmsl(numhts, max_profs), 
     1                  u(numhts,max_profs), v(numhts,max_profs),
     1                  w(numhts,max_profs), 
     1 			verif_missing_data,
     1                  baselevels(numhts)
      real            tsfc(max_profs), rhsfc(max_profs),
     1                  wssfc(max_profs), wdsfc(max_profs),
     1                  usfc(max_profs), vsfc(max_profs)
      real*8		timeobs(max_profs), prev_timeobs
      real		prof_missing

      integer           istatus, i, j, i4time
      integer           lenf
      character*9	prev_a9time
      integer           nf_fid, nf_vid, start_1(1), count_1(1),
     1                  start_2(2), count_2(2), nf_status

!..........................................................................
c     begin

      istatus = 1   !assume good return

c  open input file for reading
c
      call s_len(infile,lenf)
      print*,'reading ',infile(1:lenf)
      nf_status = nf_open(infile,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',infile
        istatus = 0
        return
      endif
c
c  read number of profilers in file and number of levels
c
c get size of recnum
c
      nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
        istatus = 0
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nprofs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
        istatus = 0
        return
      endif

      if (nprofs .gt. max_profs) then
        print *,' max_profs too small to hold profilers in file.'
        istatus = 0
        return
      endif
c
c  fill baselevels
c
      start_1(1) = 1
      count_1(1) = 28

      nf_status = nf_inq_varid(nf_fid, 'levels', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid levels'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                            count_1, baselevels(1))
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading levels 1:28'
        istatus = 0
        return
      endif

      start_1(1) = 37
      count_1(1) = 36

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                            count_1, baselevels(29))
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading levels 37:72'
        istatus = 0
        return
      endif
c
c  set start_1 and count_1 for reading 1d variables 
      start_1(1) = 1
      count_1(1) = nprofs
c
c  read wmostanum 
c
      nf_status = nf_inq_varid(nf_fid, 'wmostanum', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid wmostanum'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_int(nf_fid, nf_vid, start_1,
     1                            count_1, wmostanum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading wmostanum'
        istatus = 0
        return
      endif
c
c  read staname
c
      start_2(2) = 1
      count_2(2) = nprofs
      start_2(1) = 1
      count_2(1) = 6
c
      nf_status = nf_inq_varid(nf_fid, 'staname', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid stalat'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_text(nf_fid, nf_vid, start_2, 
     1                             count_2, staname)

      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading staname'
        istatus = 0
        return
      endif

      do i = 1, nprofs
        staname(i)(6:6) = ' '
      enddo
c
c  read stalat
c
      nf_status = nf_inq_varid(nf_fid, 'stalat', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid stalat'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                             count_1, stalat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading stalat'
        istatus = 0
        return
      endif
c
c  read stalon
c
      nf_status = nf_inq_varid(nf_fid, 'stalon', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid stalon'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                             count_1, stalon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading stalon'
        istatus = 0
        return
      endif
c
c  read staelev
c
      nf_status = nf_inq_varid(nf_fid, 'staelev', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid staelev'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                             count_1, staelev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading staelev'
        istatus = 0
        return
      endif
c
c  read timeobs
c
      nf_status = nf_inq_varid(nf_fid, 'timeobs', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid timeobs'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_double(nf_fid, nf_vid, start_1,
     1                             count_1, timeobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading timeobs'
        istatus = 0
        return
      endif

c     convert timeobs to a9_time string
      prev_timeobs = 0.
      prev_a9time = ' '
      do i = 1, nprofs
        if (timeobs(i) .eq. prev_timeobs) then
          a9_time(i) = prev_a9time
        else
          prev_timeobs = timeobs(i)
          i4time = nint(timeobs(i))  + 315619200
          call make_fnam_lp(i4time,prev_a9time,istatus)
          if(istatus .ne. 1)then
            prev_a9time = '         '
          endif
          a9_time(i) = prev_a9time
        endif
      enddo
c
c  read surface temperature
c
      nf_status = nf_inq_varid(nf_fid, 'temperature', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid temperature'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                             count_1, tsfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading temperature'
        istatus = 0
        return
      endif
c 
c  read fill value for numerical variables
c
      nf_status = nf_get_att_real(nf_fid, nf_vid, 
     1                            '_fillvalue', 
     1                            prof_missing)
      
      do j = 1, nprofs
        if (tsfc(j) .eq. prof_missing) 
     1     tsfc(j) = verif_missing_data
      enddo
c
c  read surface rh
c
      nf_status = nf_inq_varid(nf_fid, 'relhumidity', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid relhumidity'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                             count_1, rhsfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading relhumidity'
        istatus = 0
        return
      else
        do j = 1, nprofs
          if (rhsfc(j) .gt. 100.0)
     1       rhsfc(j) = verif_missing_data
        enddo
      endif
c
c  read surface wind speed
c
      nf_status = nf_inq_varid(nf_fid, 'windspeedsfc', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid windspeedsfc'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                             count_1, wssfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading windspeedsfc'
        istatus = 0
        return
      else
        do j = 1, nprofs
          if (wssfc(j) .eq. prof_missing) 
     1       wssfc(j) = verif_missing_data
        enddo
      endif
c
c  read surface wind direction 
c
      nf_status = nf_inq_varid(nf_fid, 'winddirsfc', nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'varid winddirsfc'
        istatus = 0
        return
      endif

      nf_status = nf_get_vara_real(nf_fid, nf_vid, start_1,
     1                             count_1, wdsfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading winddirsfc'
        istatus = 0
        return
      else
        do j = 1, nprofs
          if (wdsfc(j) .eq. prof_missing) 
     1       wdsfc(j) = verif_missing_data
        enddo
      endif
c
c  loop through nprofs to read 1:28 and 37:72 of u, v, w
c
      do i = 1, nprofs
        start_2(2) = i
        count_2(2) = 1
        start_2(1) = 1
        count_2(1) = 28
c
c  read u values
c
        nf_status = nf_inq_varid(nf_fid, 'ucomponent', nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'varid ucomponent'
          istatus = 0
          return
        endif

        nf_status = nf_get_vara_real(nf_fid, nf_vid, start_2,
     1                               count_2, u(1,i))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading u 1:28'
          istatus = 0
          return
        endif

        do j = 1, 28
          if (u(j,i) .eq. prof_missing) then
            u(j,i) = verif_missing_data
          else
            if (abs(u(j,i)) .gt. 200) then
              u(j,i) = verif_missing_data
              write(6,*) '*** ',nprofs, i, j, u(j,i) 
            endif
          endif
        enddo
        
        start_2(1) = 37
        count_2(1) = 36

        nf_status = nf_get_vara_real(nf_fid, nf_vid, start_2,
     1                               count_2, u(29,i))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading u 37:72'
          istatus = 0
          return
        endif
 
c       if ((i .eq. 3) .or. (i .eq. nprofs)) then
c         write(6,*) 'u bad? 54,3:',u(54,3),' 55,3:',u(55,3)
c       endif

        do j = 29, 64
          if (u(j,i) .eq. prof_missing) then
            u(j,i) = verif_missing_data
          else
            if (abs(u(j,i)) .gt. 200) then
              write(6,*) '*** ',nprofs, i, j, u(j,i) 
              u(j,i) = verif_missing_data
            endif
          endif
        enddo
c
c  read v values
c
        nf_status = nf_inq_varid(nf_fid, 'vcomponent', nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'varid vcomponent'
          istatus = 0
          return
        endif

        start_2(1) = 1
        count_2(1) = 28

        nf_status = nf_get_vara_real(nf_fid, nf_vid, start_2,
     1                               count_2, v(1,i))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading v 1:28'
          istatus = 0
          return
        endif

        do j = 1, 28
          if (v(j,i) .eq. prof_missing) then
            v(j,i) = verif_missing_data
          else
            if (abs(v(j,i)) .gt. 200) then
              write(6,*) '*** ',nprofs, i, j, v(j,i) 
              v(j,i) = verif_missing_data
            endif
          endif
        enddo

        start_2(1) = 37
        count_2(1) = 36

        nf_status = nf_get_vara_real(nf_fid, nf_vid, start_2,
     1                               count_2, v(29,i))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading v 37:72'
          istatus = 0
          return
        endif

c       if ((i .eq. 3) .or. (i .eq. nprofs)) then
c         write(6,*) 'v bad? 54,3:',v(54,3),' 55,3:',v(55,3)
c       endif

        do j = 29, 64
          if (v(j,i) .eq. prof_missing) then
            v(j,i) = verif_missing_data
          else
            if (abs(v(j,i)) .gt. 200) then
              write(6,*) '*** ',nprofs, i, j, v(j,i) 
              v(j,i) = verif_missing_data
            endif
          endif
        enddo
c
c  read w values
c
        nf_status = nf_inq_varid(nf_fid, 'wcomponent', nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'varid wcomponent'
          istatus = 0
          return
        endif

        nf_status = nf_get_vara_real(nf_fid, nf_vid, start_2,
     1                               count_2, w(1,i))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading w 1:28'
          istatus = 0
          return
        endif

        do j = 1, 28
          if (w(j,i) .eq. prof_missing) then
            w(j,i) = verif_missing_data
          else
            if (abs(w(j,i)) .gt. 20000.) then
              write(6,*) '*** ',nprofs, i, j, w(j,i) 
              w(j,i) = verif_missing_data
            endif
          endif
        enddo

        start_2(1) = 37
        count_2(1) = 36

        nf_status = nf_get_vara_real(nf_fid, nf_vid, start_2,
     1                               count_2, w(29,i))
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading w 37:72'
          istatus = 0
          return
        endif

c     if ((i .eq. 15) .or. (i .eq. nprofs)) then
c       write(6,*) 'w bad? 29,15:',v(29,15)
c       write(6,*) 'w bad? 39,15:',v(39,15)
c       write(6,*) 'w bad? 45,15:',v(45,15)
c       write(6,*) 'w bad? 47,15:',v(47,15)
c       write(6,*) 'w bad? 49,15:',v(49,15)
c       write(6,*) 'w bad? 55,15:',v(55,15)
c       write(6,*) 'w bad? 57,15:',v(57,15)
c     endif

        do j = 29, 64
          if (w(j,i) .eq. prof_missing) then
            w(j,i) = verif_missing_data
          else
            if (abs(w(j,i)) .gt. 20000.) then
              write(6,*) '*** ',nprofs, i, j, w(j,i) 
              w(j,i) = verif_missing_data
            endif
          endif
        enddo

      enddo
c
c  close input file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close ',infile
        istatus = 0
        return
      endif
c
c  fill htmsl for each profiler
c
      do i = 1, nprofs
        do j = 1, numhts
          htmsl(j,i) = staelev(i) + baselevels(j)
        enddo
      enddo
c
c  convert wssfc and wdsfc to usfc and vsfc
c
      do i = 1, nprofs
        if ((wdsfc(i) .eq. verif_missing_data) .or.
     1      (wdsfc(i) .eq. verif_missing_data)) then
          usfc(i) = verif_missing_data
          vsfc(i) = verif_missing_data
        else
          call disptrue_to_uvgrid(wdsfc(i), 
     1                            wssfc(i),
     1                            usfc(i),
     1                            vsfc(i),
     1                            stdlon)
          if (usfc(i) .gt. 100) then
            usfc(i) = verif_missing_data
            write(6,*) '*** ',nprofs, i, usfc(i) 
          endif
          if (vsfc(i) .gt. 100) then
            write(6,*) '*** ',nprofs, i, vsfc(i) 
            vsfc(i) = verif_missing_data
          endif
        endif
      enddo

c     if ((i .eq. 15) .or. (i .eq. nprofs)) then
c       write(6,*) 'tsfc bad? 15:',tsfc(15)
c       write(6,*) 'rhsfc bad? 15:',rhsfc(15)
c       write(6,*) 'usfc bad? 15:',usfc(15)
c       write(6,*) 'vsfc bad? 15:',vsfc(15)
c     endif

      return
      end
!2.........................................................................
      subroutine get_model_data(i4time, model_dir, ni, nj, nk, 
     1                          balance, laps_levels,
     1                          ht, u, v, om, t, 
     1                          ts, rhs, us, vs, istatus)

      implicit none

      integer		max_var
      parameter         (max_var=150)

      integer		i4time		!i4time of laps/model file to read
      character*(*)     model_dir	!location of model data directories
      integer           ni, nj, nk	!i, j and k grid dimensions
      integer		balance
      real            ht(ni,nj,nk)     !data to be read
      real            t(ni,nj,nk)     !data to be read
      real            u(ni,nj,nk)     !data to be read
      real            v(ni,nj,nk)     !data to be read
      real            om(ni,nj,nk)     !data to be read
      real		ts(ni,nj)
      real		rhs(ni,nj)
      real		us(ni,nj)
      real		vs(ni,nj)
      real		ps(ni,nj)
      real            make_rh
      integer		istatus		!return value from subroutine

      character*256     dir_in
      real            laps_levels(nk) !laps pressure levels
      character*3       var(max_var)
      character*4       lvl_coord(max_var)      !vertical coordinate for each field
      character*10      units(max_var)    	!units of each field
      character*125     comment(max_var)  	!comments for each field
      integer           lvl(nk)
      integer           lend, i,j, nks
      integer           dir_len
      character*31      ext
      character*3       csubdir
      integer           lt
      integer           i4time_fcst
      integer           i4time_init
      character*20      cmdltype
      character*9       a9_time
c
c     begin
c
      istatus = 1   !assume good return

      call make_fnam_lp(i4time,a9_time,istatus)
      call bkgd_wgi(a9_time,i4time_init,i4time_fcst
     +,csubdir,cmdltype,balance,istatus)
      if(istatus.ne.1)then
         print*,'failure in bkgd_wgi to get model bkgd time'
         return
      endif

      call s_len(model_dir, dir_len)
      lend=dir_len
      dir_in=model_dir
      if(dir_in(lend:lend).ne.'/')then
         lend=lend+1
         dir_in(lend:lend)='/'
      endif

      if(balance .eq. 1)then
        dir_in = model_dir(1:lend)//'balance/'
      elseif(balance .eq. 2)then      !this is the background;
c                                      must be that used in the analysis.
          if(csubdir.eq.'fua')then
             call s_len(cmdltype,lt)
      dir_in=model_dir(1:lend)//csubdir//'/'//cmdltype(1:lt)//'/'
          else
             dir_in=model_dir(1:lend)//'lga/'
          endif

      endif

c     read laps/model grid  ht
      ext = 'lt1'
      var = 'ht'
      lvl = laps_levels

      if(balance .ne. 2)then

         dir_in=dir_in(1:lend)//'lt1/'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,ht,istatus)
         if (istatus .ne. 1) then
             print*,' error reading ht variable '
             return
         endif
         var = 't3'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,t,istatus)

         if (istatus .ne. 1) then
             print*,' error reading t3 variable '
             return
         endif

      else

         ext = csubdir
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,ht,istatus)
         if (istatus .ne. 1) then
             print*,' error reading ht variable '
             return
         endif
         var = 't3'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,t,istatus)
         if (istatus .ne. 1) then
             print*,' error reading t3 variable '
             return
         endif

      endif

c     read laps/model grid  u
      ext = 'lw3'
      var = 'u3'
      if(balance.ne.2)then

         dir_in=dir_in(1:lend)//'lw3/'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,u,istatus)
         if (istatus .ne. 1) then
             print*,' error reading u3 variable '
             return
         endif
         var = 'v3'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,v,istatus)
         if (istatus .ne. 1) then
             print*,' error reading v3 variable '
             return
         endif
         var = 'om'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nk,
     1                     nk,var,lvl,lvl_coord,
     1                     units,comment,om,istatus)
         if (istatus .ne. 1) then
             print*,' error reading om variable '
             return
         endif

      else

         ext = csubdir
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,u,istatus)
         if (istatus .ne. 1) then
             print*,' error reading u3 variable '
             return
         endif
         var = 'v3'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,v,istatus)
         if (istatus .ne. 1) then
             print*,' error reading v3 variable '
             return
         endif
         var = 'om'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nk,nk,var,lvl,lvl_coord,
     1                  units,comment,om,istatus)
         if (istatus .ne. 1) then
             print*,' error reading om variable '
             return
         endif

      endif

      print*,'laps 3d: t, ht, u, v, om successfully read'

c     read laps/model grid  sfc u, v, t, rh

      if(balance.eq.0.or.balance.eq.1)then
         ext = 'lsx'
         dir_in = model_dir(1:dir_len)//'/lsx/'
         nks = 1
         lvl(1) = 0
         var(1) = 'u'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nks,
     1                     nks,var,lvl,lvl_coord,
     1                     units,comment,us,istatus)
         if (istatus .ne. 1) then
             print*,' unable to read model surface variable: u '
             return
         endif

         var(1) = 'v'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nks,
     1                     nks,var,lvl,lvl_coord,
     1                     units,comment,vs,istatus)
         if (istatus .ne. 1) then
             print*,' unable to read model surface variable: v '
             return
         endif

         var(1) = 't'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nks,
     1                     nks,var,lvl,lvl_coord,
     1                     units,comment,ts,istatus)
         if (istatus .ne. 1) then
             print*,' unable to read model surface variable: t '
             return
         endif

         var(1) = 'rh'
         call read_laps_data(i4time,dir_in,ext,ni,nj,nks,
     1                     nks,var,lvl,lvl_coord,
     1                     units,comment,rhs,istatus)
         if (istatus .ne. 1) then
             print*,' unable to read model surface variable: rh '
             return
         endif

      elseif(balance.eq.2)then

         if(csubdir.eq.'fua')then
      dir_in=model_dir(1:dir_len)//'fsf/'//cmdltype(1:lt)//'/'
            ext = 'fsf'
         else
            dir_in = model_dir(1:dir_len)//'/lgb/'
            ext = 'lgb'
         endif
         nks = 1
         lvl(1) = 0
         var(1) = 'usf'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nks,nks,var,lvl,lvl_coord,
     1                  units,comment,us,istatus)
         if (istatus .ne. 1) then
             print*,' unable to read model surface variable: usf'
             return
         endif

         var(1) = 'vsf'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nks,nks,var,lvl,lvl_coord,
     1                  units,comment,vs,istatus)
         if (istatus .ne. 1) then
             print*,' unable to read model surface variable: vsf'
             return
         endif

         var(1) = 'tsf'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nks,nks,var,lvl,lvl_coord,
     1                  units,comment,ts,istatus)
         if (istatus .ne. 1) then
             print*,' unable to read model surface variable: tsf'
             return
         endif

         var(1) = 'psf'
         call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nks,nks,var,lvl,lvl_coord,
     1                  units,comment,ps,istatus)
         if (istatus .ne. 1) then
             print*,' unable to read model surface variable: psf '
             return
         endif

         if(ext.eq.'fsf')then
            var(1) = 'rh'
            call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nks,nks,var,lvl,lvl_coord,
     1                  units,comment,rhs,istatus)
            if (istatus .ne. 1) then
                print*,' unable to read model sfc var: rh '
             return
            endif

         else

            var(1) = 'rsf'
            call read_laps(i4time_init,i4time_fcst,dir_in,ext,
     1                  ni,nj,nks,nks,var,lvl,lvl_coord,
     1                  units,comment,rhs,istatus)
            if (istatus .ne. 1) then
                print*,' unable to read model sfc var: rsf '
                return
            endif

c note: ps is in pa so no need to multiply result by 100 (rh 0-100%)
            do j=1,nj
            do i=1,ni
               rhs(i,j)=make_rh(ps(i,j),ts(i,j)-273.15,rhs(i,j)*1000.
     +,-132.)
            enddo
            enddo

         endif

      endif

      write(6,*) 'laps sfc t, rh, u, v successfully read'

      return
      end

!3.........................................................................
      subroutine get_rijk_ob_ht(whichindex, max_profs, numhts,
     1                          laps_levels, nhts,
     1                          lat_ob, lon_ob, lats,
     1                          lons, ni, nj, nk, htmsl, 
     1                          ht, ri, rj, rk, istatus)

      implicit none

      integer		whichindex
      integer		max_profs, numhts
      real            laps_levels(nk) !laps pressure levels
      integer		nhts(max_profs)
      real		lat_ob, lon_ob
      real		lats(ni,nj) 	!domain lats
      real		lons(ni,nj)	!domain lons
      integer           ni, nj, nk	!i, j and k grid dimensions
      real            htmsl(numhts,max_profs), 
     1                  ht(ni,nj,nk),
     1                  ri(numhts,max_profs),
     1                  rj(numhts,max_profs), 
     1                  rk(numhts,max_profs)
      integer		istatus

      integer		i, j, int_ri, int_rj
      real		height_to_zcoord3

c
c     begin
c
      istatus = 1   !assume good return

c     get ri, rj of ob 
      call latlon_to_rlapsgrid(lat_ob,lon_ob,lats,lons,ni,
     1                         nj,ri(1,whichindex),
     1                         rj(1,whichindex), istatus)
 
      if (istatus .eq. 1) then
        write(6,*)
        write(6,*) ' profiler in domain: ', whichindex
      else
c       write(6,*) ' profiler out of domain: ',whichindex
        istatus = -1
        return
      endif

c     fill rest of ri,rj for height column of ob
      do i = 2, numhts
        ri(i,whichindex) = ri(1,whichindex)
        rj(i,whichindex) = rj(1,whichindex)
      enddo

c     get integral ri, rj for height_to_zcoord3  
      int_ri = nint(ri(1,whichindex))
      int_rj = nint(rj(1,whichindex))

c     calc rk for each level
      do i = 1, numhts
        rk(i,whichindex) = height_to_zcoord3(htmsl(i,whichindex),
     1                                       ht,laps_levels,
     1                                       ni,nj,nk,int_ri,
     1                                       int_rj,istatus)

        if (istatus .ne. 1) then
          if (i .gt. 1) then
            nhts(whichindex) = i - 1
            istatus = 1
            return
          else
            write(6,*) ' unable to get rk for profiler ',whichindex,
     1               '  height index= ',i,'maxht=',
     1               ht(int_ri,int_rj,nk),'pht=',htmsl(i,whichindex)
            istatus = 0
            return
          endif
        endif
      enddo

      return
      end
!4.........................................................................
      subroutine write_verif_prof(output_fname, max_profs, numhts,
     1                            nhts, n_profs_use, use_prof, 
     1                            nprofs, wmostanum,
     1                            staname, stalat, stalon, staelev,
     1                            htmsl, ri, rj, rk, up, vp, wp, tsp,
     1                            rhsp, usp, vsp, ui, vi, wi, tsi,
     1                            rhsi, usi, vsi, a9_time, istatus)

      implicit none

      character*(*)     output_fname	!path and name of output file
      integer		max_profs, numhts, nhts(max_profs)
      integer           n_profs_use, use_prof(max_profs), nprofs
      integer		wmostanum(max_profs) 
      character*6       staname(max_profs) 
      real            stalat(max_profs), stalon(max_profs), 
     1                  staelev(max_profs), 
     1                  htmsl(numhts, max_profs), 
     1                  ri(numhts,max_profs),
     1                  rj(numhts,max_profs), 
     1                  rk(numhts,max_profs),
     1                  up(numhts,max_profs), 
     1                  vp(numhts,max_profs),
     1                  wp(numhts,max_profs), 
     1                  tsp(max_profs), rhsp(max_profs),
     1                  usp(max_profs), vsp(max_profs)
      real            ui(numhts,max_profs), 
     1                  vi(numhts,max_profs),
     1                  wi(numhts,max_profs), 
     1                  tsi(max_profs), rhsi(max_profs),
     1                  usi(max_profs), vsi(max_profs)
      character*9	a9_time(max_profs)
      integer		istatus

      integer		i, j

c
c     begin
c
      istatus = 1   !assume good return

      write(6,*) output_fname

c     open output_fname
      open(1,file=output_fname,status='unknown',err=98)
      go to 99

98    write(6,*)' error opening verif file: ',output_fname
      istatus = 0
      return

99    continue

c 100 writes: n_profs_use
100   format(i3)
c 101 writes: nhts+sfc,a9_time,wmostanum,staname,stalat,stalon,staelev
101   format(i3,1x,a9,1x,i7,1x,a6,1x,f7.3,1x,f8.3,1x,f7.0)
c 102 writes: staelev,ri,rj,usp,usi,vsp,vsi,tsp,tsi,rhsp,rhsi
102   format(f7.0,1x,2(f8.3,1x),9x,8(f6.1,1x))
c 103 writes: htmsl,ri,rj,rk,up,ui,vp,vi,wp,wi
103   format(f7.0,1x,3(f8.3,1x),6(f6.1,1x))
c104  write: header for sfc variables
104   format(1x,'stnhgt',4x,'ri',7x,'rj',16x,'usp',4x,'usi',
     +4x,'vsp',4x,'vsi',4x,'tsp',4x,'tsi',3x,'rhsp',3x,'rhsi')
c 1010writes: header for ua variables
105   format(4x,'hgt',4x,'ri',7x,'rj',7x,'rk',8x,'up',5x,
     +'ui',5x,'vp',5x,'vi',5x,'wp',5x,'wi')
c     write number of profilers into file
      write(1,100) n_profs_use
      
c     write each profiler out where use_prof(i) .eq. 1
      do i = 1, nprofs
        if (use_prof(i) .eq. 1) then
c         write header info
          write(1,101) nhts(i)+1,a9_time(i),wmostanum(i),
     1                 staname(i)(1:6),
     1                 stalat(i),stalon(i),staelev(i) 
c         write surface data
          write(1,104)
          write(1,102) staelev(i),ri(1,i),rj(1,i),usp(i),usi(i),
     1                 vsp(i),vsi(i),tsp(i),tsi(i),rhsp(i),rhsi(i)

          write(1,105)
          do j = 1, nhts(i)
            write(1,103) htmsl(j,i),ri(j,i),rj(j,i),rk(j,i),up(j,i),
     1                   ui(j,i),vp(j,i),vi(j,i),wp(j,i),wi(j,i)
          enddo
        endif
      enddo

      close(1)

      return
      end
!5.........................................................................
      subroutine interp_to_ob(whichindex, max_profs, numhts,
     1                        nhts,
     1                        ni, nj, nk, ri, rj, rk, um, 
     1                        vm, omm, tm, tsm, rhsm, usm,
     1                        vsm, ui, vi, wi, tsi, rhsi,
     1                        usi, vsi, r_missing_data,
     1                        verif_missing_data, istatus) 

      implicit none

      integer		whichindex, max_profs, numhts
      integer           nhts(max_profs)
      integer           ni, nj, nk	!i, j and k grid dimensions
      real            ri(numhts,max_profs),
     1                  rj(numhts,max_profs), 
     1                  rk(numhts,max_profs),
     1                  ui(numhts,max_profs), 
     1                  vi(numhts,max_profs),
     1                  wi(numhts,max_profs), 
     1                  tsi(max_profs), rhsi(max_profs),
     1                  usi(max_profs), vsi(max_profs),
     1                  tm(ni,nj,nk), um(ni,nj,nk),
     1                  vm(ni,nj,nk), omm(ni,nj,nk),
     1      		tsm(ni,nj), rhsm(ni,nj),
     1      		usm(ni,nj), vsm(ni,nj)
      real		r_missing_data, verif_missing_data
      integer		istatus

      integer           int_ri, int_rj, i, extrap
      real		om, temp, pres
      real		pressure_of_rlevel
c
c     begin
c
      istatus = 1   !assume good return

c     write(6,*) 'ri=',ri(1,whichindex),' rj=',rj(1,whichindex)
c     int_ri = nint(ri(1,whichindex))
c     int_rj = nint(rj(1,whichindex))
c     write(6,*) 'int_ri=',int_ri,' int_rj=',int_rj
c     write(6,*)

c     determine if extrapolation is needed 
      if (((ri(1,whichindex) .lt. 1.0) .and. 
     1     (ri(1,whichindex) .gt. 0.5)) 
     1    .or.
     1    ((ri(1,whichindex) .lt. ni+0.5) .and. 
     1     (ri(1,whichindex) .gt. ni))
     1    .or.
     1    ((rj(1,whichindex) .lt. 1.0) .and. 
     1     (ri(1,whichindex) .gt. 0.5))
     1    .or.
     1    ((rj(1,whichindex) .lt. nj+0.5) .and. 
     1     (ri(1,whichindex) .gt. nj))) then
        extrap = 1
      else
        extrap = 0
      endif

c     interpolate surface variable: usi from usm
      if (extrap .eq. 1) then
        call bilinear_interp_extrap(ri(1,whichindex),
     1                              rj(1,whichindex),ni,nj,
     1                              usm,usi(whichindex),
     1                              istatus)
      else
        call bilinear_laps(ri(1,whichindex),rj(1,whichindex),
     1                     ni,nj,usm,usi(whichindex))
      endif

      if (usi(whichindex) .eq. r_missing_data)
     1  usi(whichindex) = verif_missing_data

c     write(6,*) 'usm(int_ri,int_rj)=',usm(int_ri,int_rj)
c     write(6,*) 'usi(',whichindex,')=',usi(whichindex)
c     write(6,*)

c     interpolate surface variable: vsi from vsm
      if (extrap .eq. 1) then
        call bilinear_interp_extrap(ri(1,whichindex),
     1                              rj(1,whichindex),ni,nj,
     1                              vsm,vsi(whichindex),
     1                              istatus)
      else
        call bilinear_laps(ri(1,whichindex),rj(1,whichindex),
     1                     ni,nj,vsm,vsi(whichindex))
      endif

      if (vsi(whichindex) .eq. r_missing_data)
     1  vsi(whichindex) = verif_missing_data

c     write(6,*) 'vsm(int_ri,int_rj)=',vsm(int_ri,int_rj)
c     write(6,*) 'vsi(',whichindex,')=',vsi(whichindex)
c     write(6,*)

c     interpolate surface variable: rhsi from rhsm
      if (extrap .eq. 1) then
        call bilinear_interp_extrap(ri(1,whichindex),
     1                              rj(1,whichindex),ni,nj,
     1                              rhsm,rhsi(whichindex),
     1                              istatus)
      else

        call bilinear_laps(ri(1,whichindex),rj(1,whichindex),
     1                     ni,nj,rhsm,rhsi(whichindex))
      endif

      if (rhsi(whichindex) .eq. r_missing_data)
     1  rhsi(whichindex) = verif_missing_data

c     write(6,*) 'rhsm(int_ri,int_rj)=',rhsm(int_ri,int_rj)
c     write(6,*) 'rhsi(',whichindex,')=',rhsi(whichindex)
c     write(6,*)

c     interpolate surface variable: tsi from tsm
      if (extrap .eq. 1) then
        call bilinear_interp_extrap(ri(1,whichindex),
     1                              rj(1,whichindex),ni,nj,
     1                              tsm,tsi(whichindex),
     1                              istatus)
      else
        call bilinear_laps(ri(1,whichindex),rj(1,whichindex),
     1                     ni,nj,tsm,tsi(whichindex))
      endif

      if (tsi(whichindex) .eq. r_missing_data)
     1  tsi(whichindex) = verif_missing_data

c     write(6,*) 'tsm(int_ri,int_rj)=',tsm(int_ri,int_rj)
c     write(6,*) 'tsi(',whichindex,')=',tsi(whichindex)
c     write(6,*)

c     loop through levels for u,v,w
      do i = 1, nhts(whichindex)

c       interpolate 3d variable: ui from um
        if (extrap .eq. 1) then
          call trilinear_interp_extrap(ri(i,whichindex),
     1                                 rj(i,whichindex),
     1                                 rk(i,whichindex),
     1                                 ni,nj,nk, um,
     1                                 ui(i,whichindex),
     1                                 istatus)
        else
          call trilinear_laps(ri(i,whichindex),rj(i,whichindex),
     1                        rk(i,whichindex),ni,nj,nk, um,
     1                        ui(i,whichindex))
        endif

        if (ui(i,whichindex) .eq. r_missing_data)
     1    ui(i,whichindex) = verif_missing_data

c       interpolate 3d variable: vi from vm
        if (extrap .eq. 1) then
          call trilinear_interp_extrap(ri(i,whichindex),
     1                                 rj(i,whichindex),
     1                                 rk(i,whichindex),
     1                                 ni,nj,nk, vm,
     1                                 vi(i,whichindex),
     1                                 istatus)
        else
          call trilinear_laps(ri(i,whichindex),rj(i,whichindex),
     1                        rk(i,whichindex),ni,nj,nk, vm,
     1                        vi(i,whichindex))
        endif

        if (vi(i,whichindex) .eq. r_missing_data)
     1    vi(i,whichindex) = verif_missing_data

c       interpolate 3d variable: om from omm
        if (extrap .eq. 1) then
          call trilinear_interp_extrap(ri(i,whichindex),
     1                                 rj(i,whichindex),
     1                                 rk(i,whichindex),
     1                                 ni,nj,nk, omm,
     1                                 om,
     1                                 istatus)
        else
          call trilinear_laps(ri(i,whichindex),rj(i,whichindex),
     1                        rk(i,whichindex),ni,nj,nk, omm,
     1                        om)
        endif

c       interpolate 3d variable: temp from tm
        if (extrap .eq. 1) then
          call trilinear_interp_extrap(ri(i,whichindex),
     1                                 rj(i,whichindex),
     1                                 rk(i,whichindex),
     1                                 ni,nj,nk, tm,
     1                                 temp,
     1                                 istatus)
        else
          call trilinear_laps(ri(i,whichindex),rj(i,whichindex),
     1                        rk(i,whichindex),ni,nj,nk, tm,
     1                        temp)
        endif

        if ((om .eq. r_missing_data) .or.
     1      (temp .eq. r_missing_data))  then
          wi(i,whichindex) = verif_missing_data
        else
c         get pressure of rk
          pres = pressure_of_rlevel(rk(i,whichindex))

c         calclulate wi(i,whichindex) from om   w = -om*(r*temp/g*pres)
          wi(i,whichindex) = -1*om*((287.04*temp)/(9.80815*pres))
        endif

      enddo

      return
      end
!6.........................................................................

