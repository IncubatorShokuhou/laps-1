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
      subroutine get_raob_pairs(raob_fname,model_dir,i4time_sys,i4time,
     1                          i4time_raob_earliest,i4time_raob_latest,
     1                          output_fname, nl_dir, ni, nj,
     1                          nk, lats, lons, stdlon, 
     1                          laps_levels_mb,laps_levels_pa,
     1                          max_ht_m_proc, min_pres_mb_proc,
     1                          balance, r_missing_data, 
     1                          verif_missing_data, istatus)

      implicit none

      character*(*)     raob_fname	!path and name of raob file to read
      character*(*)     model_dir	!location of model data directories
                                        !lapsprd, or location of fua, fsf
      integer		i4time		!i4time of laps/model file to read
      integer           i4time_sys
      integer         i4time_raob_latest, i4time_raob_earliest
      character*256     output_fname	!path and name of output file
      character*(*)	nl_dir		!directory where verify_raob.nl located
      integer           ni, nj, nk	!i, j and k grid dimensions
      real		lats(ni,nj) 	!domain lats
      real		lons(ni,nj)	!domain lons
      real 		stdlon		!standard longitude
      real            laps_levels_mb(nk) !laps pressure levels
      real            laps_levels_pa(nk) !laps pressure levels
      real            max_ht_m_proc   !maximum height(m) to process up to
      real            min_pres_mb_proc  !minimum pressure(mb) to process up to

      integer		balance
      real		r_missing_data, raob_missing_data
      real		verif_missing_data
      integer		istatus		!return value from subroutine

      integer		max_raobs
      parameter         (max_raobs=50)
      integer           maxm
      parameter         (maxm=22)
      integer           maxw
      parameter         (maxw=97)  !75+22, for mingle mand and sigw
      integer		maxt
      parameter         (maxt=177) !150+22 for mingle mand and sigt
      integer		max_hts
      parameter		(max_hts=275)

      integer 		nf_status, ncid, varid
      integer           all_raobs	! 1=use all raobs in domain
      integer           use_raob(max_raobs), n_raobs_use
      integer           lun, i, j, dir_len
      integer 	statusl(2,6)
      character*256     filename
      logical           l_eof

c     1d raob data
      integer		nraobs, n_raobs_avail
      integer 		wmonum_use(max_raobs)
      character*6       staname_use(max_raobs)
      integer 	timesyn(max_raobs),
     1			wmostanum(max_raobs),
     1			timerel(max_raobs),
     1			numsigt(max_raobs),
     1			numsigw(max_raobs) 
      character*6       staname(max_raobs) 
      real            stalat(max_raobs), 
     1			stalon(max_raobs), 
     1                  staelev(max_raobs)
      character*9	a9_time(max_raobs)

c     temp variables
      integer 	numt(max_raobs),
     1                  timelapst(maxt,max_raobs),
     1			fileavailthr
      character*1	typet(maxt,max_raobs)
      real            tsigt(maxt,max_raobs),
     1                  tdsigt(maxt,max_raobs),
     1                  htsigt(maxt,max_raobs),
     1                  prsigt(maxt,max_raobs),
     1			prit(maxt,max_raobs), 
     1                  tit(maxt,max_raobs),
     1                  tdit(maxt,max_raobs),
     1                  prpt(maxt,max_raobs), 
     1                  tpt(maxt,max_raobs),
     1                  tdpt(maxt,max_raobs),
     1                  rit(maxt,max_raobs),
     1                  rjt(maxt,max_raobs), 
     1                  rkt(maxt,max_raobs),
     1                  latt(maxt,max_raobs), 
     1                  lont(maxt,max_raobs), 
     1                  htt(maxt,max_raobs) 

c     wind variables
      integer 	numw(max_raobs),
     1                  timelapsw(maxw,max_raobs),
     1                  fileavailuv
      character*1	typew(maxw,max_raobs)
      real            htsigw(maxw,max_raobs), 
     1                  wssigw(maxw,max_raobs), 
     1                  wdsigw(maxw,max_raobs), 
     1                  uiw(maxw,max_raobs), 
     1                  viw(maxw,max_raobs),
     1                  upw(maxw,max_raobs), 
     1                  vpw(maxw,max_raobs),
     1                  riw(maxw,max_raobs),
     1                  rjw(maxw,max_raobs), 
     1                  rkw(maxw,max_raobs),
     1                  latw(maxw,max_raobs), 
     1                  lonw(maxw,max_raobs), 
     1                  htw(maxw,max_raobs) 

c     data written out
      integer         nhts(max_raobs),
     1                  timelaps(max_hts,max_raobs)
      character*1	type(max_hts,max_raobs)
      real            ri(max_hts,max_raobs),
     1                  rj(max_hts,max_raobs), 
     1                  rk(max_hts,max_raobs),
     1                  lat(max_hts,max_raobs), 
     1                  lon(max_hts,max_raobs), 
     1                  hts(max_hts,max_raobs),
     1			lapstime(max_hts,max_raobs),
     1                  up(max_hts,max_raobs), 
     1                  ui(max_hts,max_raobs), 
     1                  vp(max_hts,max_raobs), 
     1                  vi(max_hts,max_raobs),
     1                  tp(max_hts,max_raobs),
     1                  ti(max_hts,max_raobs), 
     1                  tdp(max_hts,max_raobs),
     1                  tdi(max_hts,max_raobs),
     1                  prp(max_hts,max_raobs), 
     1                  pri(max_hts,max_raobs)

c     laps data read in
      real            ulapsgs(ni,nj,nk), vlapsgs(ni,nj,nk),
     1                  tlapsgs(ni,nj,nk), rhlapsgs(ni,nj,nk),
     1                  htlapsgs(ni,nj,nk),ulapsgp(ni,nj,nk),
     1                  vlapsgp(ni,nj,nk),tlapsgp(ni,nj,nk),
     1                  rhlapsgp(ni,nj,nk), htlapsgp(ni,nj,nk),
     1                  htlgags(ni,nj,nk), htlgagp(ni,nj,nk)
       
      integer		writet,writew
c
c     begin
c
      istatus = 1   !assume good return
      writew = 0
      writet = 0

c     set up nhts(max_raobs)
      do i = 1, max_raobs
        nhts(i) = max_hts
      enddo

c     read verify_raob.nl

      lun = 10
      call s_len(nl_dir,dir_len)
      filename = nl_dir(1:dir_len)//'/verif_raob.txt'
      open(lun,file=filename,status='old',err=900)

      l_eof = .false.
      n_raobs_use = 0
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
      n_raobs_use = i - 1

      all_raobs = 0	!set to 1 if use all raobs in file

c     see if wmonum_use(1) .eq. -1...if so, use all raobs in file
      if (wmonum_use(1) .eq. -1) then
        all_raobs = 1
      endif

c     read raob file "raob_fname" and return requested raob data
c     with man/sigt mingled by pressure and man/sigw mingled by height
      call get_raob_data_a(ni, nj, i4time,i4time_raob_earliest,
     1                     i4time_raob_latest,raob_fname,max_raobs,
     1                     maxm,maxt, maxw,lats,lons,timesyn,timerel, 
     1                     numsigt, numsigw, wmostanum,staname,typew,
     1                     typet,prsigt, tsigt, tdsigt, htsigt, htsigw,
     1                     wdsigw,wssigw, stalat, stalon, staelev,
     1                     max_ht_m_proc, min_pres_mb_proc,
     1                     nraobs, n_raobs_avail,verif_missing_data,
     1                     raob_missing_data, istatus)

      if (istatus .ne. 1) then
        write(6,*) ' unable to read raob file: ',raob_fname
        istatus = 0
        return
      endif

      write(6,*)
      write(6,*) 'number of raobs available: ',n_raobs_avail
      write(6,*)

      if (n_raobs_avail .gt. 0) then
   
c       set which raobs to pull laps/model data from
        if (all_raobs .ne. 1) then   !set use_raob(i) to 1 if raob is on list
          do i = 1, n_raobs_use
            do j = 1, nraobs
              if ((wmonum_use(i) .eq. wmostanum(j)) .or.
     1            (staname_use(i) .eq. staname(j))) then
                use_raob(j) = 1    !set to 1 if current raob is on list to use
              else
                use_raob(j) = 0  
              endif
            enddo
          enddo
        else
          do i = 1, nraobs
            use_raob(i) = 1    !set to 1 if current raob is on list to use
          enddo
        endif

c       for i4time read laps lw3 file for u and v/ lt1 file for t and ht
c       and lh3 file for rh
c       statusl(1,x) = syn  (2,x) = prev
c       statusl(x,1)=u (x,2)=v (x,3)=t (x,4)=ht (x,5)=rh (x,6)=lgaht

        call make_fnam_lp(i4time,a9_time,istatus)

        call get_laps(ni,nj,nk,model_dir,i4time_sys,a9_time,
     1              ulapsgs,vlapsgs,tlapsgs,rhlapsgs,htlapsgs,
     1                ulapsgp,vlapsgp,tlapsgp,rhlapsgp,
     1                htlapsgp,htlgags,htlgagp,balance,
     1                laps_levels_mb,statusl)

        fileavailuv = 0
        if ((statusl(1,1) .eq. 1) .and. (statusl(1,2) .eq. 1)) then
          fileavailuv = fileavailuv + 1
        elseif ((statusl(2,1) .eq. 1) .and.
     1          (statusl(2,2) .eq. 1)) then
          fileavailuv = fileavailuv + 2
        else
          istatus = 0
          write(6,*) ' missing u/v data ',statusl(1,1),
     1    statusl(1,2),' 11z ',statusl(2,1), statusl(2,2)
        endif

        fileavailthr = 0
        if ((statusl(1,3) .eq. 1) .and. (statusl(1,4) .eq. 1) .and.
     1      (statusl(1,5) .eq. 1)) then
          fileavailthr = fileavailthr + 1
        elseif ((statusl(2,3) .eq. 1) .and.
     1          (statusl(2,4) .eq. 1) .and.
     1          (statusl(2,5) .eq. 1)) then
          fileavailthr = fileavailthr + 2
        else
          istatus = 0
          write(6,*) ' missing t/ht/td data ',statusl(1,1),
     1    statusl(1,2),' 11z ',statusl(2,1), statusl(2,2)
        endif

        if (fileavailuv .gt. 0) then !need uv for both..if not there, can't run

c         loop through use_raob and pull model data if use_raob(i) .eq. 1
          do i = 1, nraobs
            if (use_raob(i) .eq. 1) then

c             ascend balloon for w and return ri,rj,rk,latw,lonw,htw, uiw,viw,
c                                           timelapsw,istatus
              if (fileavailuv .gt. 0) then

                call ascend_w(timesyn(i), timerel(i), numsigw(i),
     1                      numw(i),fileavailuv,wmostanum(i),
     1                      stalat(i), stalon(i), staelev(i), maxw,
     1                      max_ht_m_proc, typew,
     1                      max_raobs, ni, nj, nk, i, statusl,
     1                      htsigw, wdsigw, wssigw,  !remember (maxw,max_raobs)
     1                      ulapsgs,vlapsgs,ulapsgp,vlapsgp,  !(ni,nj,nk)
     1                      htlapsgs, htlapsgp, raob_missing_data,
     1                      verif_missing_data,
     1                      lats,lons,laps_levels_pa,    !(ni,nj,nk)
!................................variables below returned...................
     1                      riw,rjw,rkw,latw,lonw,htw, uiw,viw,
     1                      upw,vpw,timelapsw,istatus)

                if (istatus .ne. 1) then
                  write(6,*) 'wind data for raob ',wmostanum(i),
     1                     ' not available.'
           
                  numw(i) = 0
                endif
              else
                numw(i) = 0
              endif

              if (fileavailthr .gt. 0) then

                call ascend_t(timesyn(i), timerel(i), numsigt(i),
     1                        numt(i),fileavailuv,fileavailthr, 
     1                        wmostanum(i),stalat(i), stalon(i), 
     1                        staelev(i), typet,
     1                        maxt, maxw,max_ht_m_proc,
     1                        max_raobs, ni, nj, nk, i, statusl,
     1                        prsigt, tsigt, tdsigt,htsigt,  !remember (maxt,max_raobs)
     1                        tlapsgs,rhlapsgs,htlapsgs,tlapsgp,  !(ni,nj,nk)
     1                        rhlapsgp, htlapsgp, 
     1                        htlgags, htlgagp,raob_missing_data,
     1                        verif_missing_data, lats,lons, !(ni,nj,nk)
     1                        laps_levels_pa, numsigw(i),
     1                        htsigw, wdsigw, wssigw, !remember (maxw,max_raobs)
!................................variables below returned...................
     1                        rit,rjt,rkt,latt,lont,htt,prit,tit,
     1                        tdit,prpt,tpt,tdpt,timelapst,istatus)

                if (istatus .ne. 1) then
                  numt(i) = 0
                  write(6,*) 'temp data for raob ',wmostanum(i),
     1                       ' not available.'
                endif
              else
                numt(i) = 0
              endif

              if ((numt(i) .eq. 0) .and. (numw(i) .eq. 0)) then
                use_raob(i) = 0
              else
                if (numt(i) .ne. 0) writet = writet + 1
                if (numw(i) .ne. 0) writew = writew + 1
              endif

            endif
          enddo

c         re-calc n_raobs_use based on interpolation errors
          n_raobs_use = 0
          do i = 1, nraobs
            if (use_raob(i) .eq. 1) n_raobs_use = n_raobs_use + 1
          enddo

c         merge t and w data
c         call mergetw(max_hts,max_raobs,maxw,maxt,numw,numt,
c    1               nraobs,use_raob,n_raobs_use,
c    1               ri,rj,rk,lat,lon,hts,type,up,ui,vp,vi,
c    1               tp,ti,tdp,tdi,prp,pri,timelaps,nhts,
c    1               rit,rjt,rkt,latt,lont,htt,typet,
c    1               riw,rjw,rkw,latw,lonw,htw,typew, 
c    1               uiw,viw,upw,vpw,timelapsw,
c    1               prit,tit,tdit,prpt,tpt,tdpt,timelapst,
c    1               verif_missing_data, raob_missing_data,
c    1               istatus)


c         write output files

          if (n_raobs_use .gt. 0) then
            call write_verif_raob(output_fname, max_raobs,maxw,
     1                          maxt,numw,numt,nraobs,
     1                          writet, writew,
     1                          n_raobs_use,use_raob, wmostanum,
     1                          staname,stalat,stalon, staelev,
     1                          a9_time,timesyn,timerel,
     1                          rit,rjt,rkt,latt,lont,htt,
     1                          riw,rjw,rkw,latw,lonw,htw, 
     1                          uiw,viw,upw,vpw,timelapsw,
     1                          prit,tit,tdit,prpt,tpt,tdpt,
     1                          timelapst, istatus)

            if (istatus .ne. 1) then
              write(6,*) ' error writing raob verif file: ',
     1                    output_fname
            endif
          else
            write(6,*) 'no raobs to output'
          endif

        else !no model data files to process
          write(6,*) 'no uv data available...cannot verify raobs'
        endif
        goto 999

      else
        write(6,*)'no raobs to output' 
        goto 999
      endif

900   print*,'error opening namelist file ', filename
      goto 999

901   print*,'error reading namelist file ', filename
      goto 999


999   continue

      return
      end
!1.........................................................................
      subroutine write_verif_raob(output_fname, max_raobs, maxw,
     1                          maxt,numw,numt,nraobs,
     1                          writet,writew,
     1                          n_raobs_use,use_raob, wmostanum,
     1                          staname,stalat,stalon,staelev,
     1                          a9_time,timesyn,timerel,
     1                          rit,rjt,rkt,latt,lont,htt,
     1                          riw,rjw,rkw,latw,lonw,htw, 
     1                          uiw,viw,upw,vpw,timelapsw,
     1                          prit,tit,tdit,prpt,tpt,tdpt,
     1                          timelapst,istatus)

      implicit none

      character*256     output_fname	!path and name of output file
      integer		max_raobs, maxw,maxt
      integer           nraobs, n_raobs_use, use_raob(max_raobs)
      integer		writet,writew,wmostanum(max_raobs) 
      character*6       staname(max_raobs) 
      real            stalat(max_raobs), stalon(max_raobs), 
     1                  staelev(max_raobs)
      character*9	a9_time(max_raobs)
      integer		timesyn(max_raobs),
     1         		timerel(max_raobs)

c     temp variables
      integer         numt(max_raobs),
     1			timelapst(maxt,max_raobs)
      real            rit(maxt,max_raobs),
     1                  rjt(maxt,max_raobs),
     1                  rkt(maxt,max_raobs),
     1                  latt(maxt,max_raobs),
     1                  lont(maxt,max_raobs),
     1                  htt(maxt,max_raobs),
     1                  prit(maxt,max_raobs),
     1                  tit(maxt,max_raobs),
     1                  tdit(maxt,max_raobs),
     1                  prpt(maxt,max_raobs),
     1                  tpt(maxt,max_raobs),
     1                  tdpt(maxt,max_raobs)

c     wind variables
      integer         numw(max_raobs),
     1 			timelapsw(maxw,max_raobs)
      real            riw(maxw,max_raobs),
     1                  rjw(maxw,max_raobs),
     1                  rkw(maxw,max_raobs),
     1                  latw(maxw,max_raobs),
     1                  lonw(maxw,max_raobs),
     1                  htw(maxw,max_raobs),
     1                  uiw(maxw,max_raobs),
     1                  viw(maxw,max_raobs),
     1                  upw(maxw,max_raobs),
     1                  vpw(maxw,max_raobs)

      integer		istatus

c     local variables
      integer		i, j
      character*256     output_fnamew   !path/name of wind output file
      character*256     output_fnamet   !path/name of temp output file
      integer		output_len

c
c     begin
c
      istatus = 1   !assume good return

      call s_len(output_fname,output_len)
      output_fnamet = output_fname(1:output_len)//'t'

      if (writet .ne. 0) then
        write(6,*) output_fnamet

c       open output_fname
        open(1,file=output_fnamet,status='unknown',err=98)
        go to 99

98      write(6,*)' error opening temp verif file: ',output_fnamet
        istatus = 0
        return

99      continue

c 100   writes: n_raobs_use
100     format(i3)
c 101   writes: nhts,a9_time,wmostanum,staname,stalat,stalon,staelev,timesyn,timerel
101     format(i3,1x,a9,1x,i7,1x,a6,1x,f7.3,1x,f8.3,1x,f7.0,2(1x,i12))
c 102   writes: htt,latt,lont,rit,rjt,rkt,tpt,tit,tdpt,tdit,prpt,prit,timelapst
102     format(f7.1,1x,f6.2,1x,f8.2,1x,3(f8.3,1x),6(f6.1,1x),i12)


c       write number of raob into file
        write(1,100) writet
      
c       write temp data out for each raob where use_raob(i) .eq. 1
        do i = 1, nraobs
          if ((use_raob(i) .eq. 1).and.(numt(i).gt.0)) then
c           write header info
            write(1,101)numt(i),a9_time(i),wmostanum(i),
     1                 staname(i)(1:6),
     1                 stalat(i),stalon(i),staelev(i),
     1                 timesyn(i),timerel(i)


            write(1,104)
            do j = 1, numt(i)
              write(1,102) htt(j,i),latt(j,i),lont(j,i),
     1                   rit(j,i),rjt(j,i),rkt(j,i),tpt(j,i),
     1                   tit(j,i),tdpt(j,i),tdit(j,i),
     1                   prpt(j,i),prit(j,i),timelapst(j,i)
            enddo
          endif
        enddo

        close(1)
      else
        write(6,*) 'no temp data written: ',output_fnamet
      endif
104     format(2x,'hgt',5x,'lat',5x,'lon',6x,'ri',7x,'rj',7x,'rk',
     +7x,'tpt',4x,'tit',3x,'tdpt',3x,'tdit',3x,'prpt',3x,'prit'
     +,3x,'i4timelapst')

      output_fnamew = output_fname(1:output_len)//'w'
      if (writew .gt. 0) then
        write(6,*) output_fnamew

c       open output_fname
        open(1,file=output_fnamew,status='unknown',err=198)
        go to 199

198     write(6,*)' error opening wind verif file: ',output_fnamew
        istatus = 0
        return

199     continue

c 103   writes: htw,latw,lonw,riw,rjw,rkw,upw,uiw,vpw,viw,timelapsw
103     format(f7.1,1x,f6.2,1x,f8.2,1x,3(f8.3,1x),4(f6.1,1x),i12)

c       write number of raob into file
        write(1,100) n_raobs_use
      
c       write wind data out for each raob where use_raob(i) .eq. 1
        do i = 1, nraobs
          if ((use_raob(i) .eq. 1).and.(numw(i).gt.0)) then
c           write header info
            write(1,101)numw(i),a9_time(i),wmostanum(i),
     1                 staname(i)(1:6),
     1                 stalat(i),stalon(i),staelev(i),
     1                 timesyn(i),timerel(i)

            write(1,105)
            do j = 1, numw(i)
            write(1,103) htw(j,i),latw(j,i),lonw(j,i),riw(j,i),
     1                   rjw(j,i),rkw(j,i),upw(j,i),
     1                   uiw(j,i),vpw(j,i),viw(j,i),timelapsw(j,i)
            enddo
          endif
        enddo

        close(1)
      else
        write(6,*) 'no wind data to write out: ',output_fnamew
      endif

105   format(2x,'hgt',5x,'lat',5x,'lon',6x,'ri',7x,'rj',7x,'rk',
     +7x,'upw',4x,'uiw',4x,'vpw',4x,'viw',3x,'i4timelapsw')

      return
      end
!2.........................................................................
