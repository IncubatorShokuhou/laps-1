
      subroutine get_raob_data(
     1                         i4time_sys,ilaps_cycle_time,nx_l,ny_l
     1                        ,i4time_raob_earliest,i4time_raob_latest       
     1                        ,filename,lun_out,l_fill_ht
     1                        ,istatus)

      character*170 filename
      logical l_fill_ht

!...........................................................................

      include 'netcdf.inc'
      integer mtropnum, mwndnum, manlevel, recnum, sigtlevel,
     +     sigwlevel,stanamelen, nf_fid, nf_vid, nf_status,
     +     nlvl_out

c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',filename
      endif
c
c  fill all dimension values
c
c
c get size of mtropnum
c
      nf_status = nf_inq_dimid(nf_fid,'mtropnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim mtropnum'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,mtropnum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim mtropnum'
      endif
c
c get size of mwndnum
c
      nf_status = nf_inq_dimid(nf_fid,'mwndnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim mwndnum'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,mwndnum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim mwndnum'
      endif
c
c get size of manlevel
c
      nf_status = nf_inq_dimid(nf_fid,'manlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim manlevel'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,manlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim manlevel'
      endif
c
c get size of recnum
c
      nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,recnum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
      endif
c
c get size of sigtlevel
c
      nf_status = nf_inq_dimid(nf_fid,'sigtlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigtlevel'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,sigtlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigtlevel'
      endif
c
c get size of sigwlevel
c
      nf_status = nf_inq_dimid(nf_fid,'sigwlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigwlevel'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,sigwlevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim sigwlevel'
      endif
c
c get size of stanamelen
c
      nf_status = nf_inq_dimid(nf_fid,'stanamelen',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim stanamelen'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,stanamelen)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim stanamelen'
      endif

      nlvl_out = manlevel + sigtlevel + sigwlevel
      call main_sub(nf_fid, mtropnum, mwndnum, manlevel, recnum,
     +     sigtlevel, sigwlevel, stanamelen, nlvl_out
!.............................................................................
     1                        ,i4time_sys,ilaps_cycle_time,nx_l,ny_l
     1                        ,i4time_raob_earliest,i4time_raob_latest  
     1                        ,l_fill_ht,lun_out     
     1                        ,istatus)


      return
!.............................................................................
      end
c
c
      subroutine main_sub(nf_fid, mtropnum, mwndnum, manlevel, recnum,
     +     sigtlevel, sigwlevel, stanamelen, nlvl_out
!.............................................................................
     1                        ,i4time_sys,ilaps_cycle_time,nx_l,ny_l
     1                        ,i4time_raob_earliest,i4time_raob_latest       
     1                        ,l_fill_ht,lun_out     
     1                        ,istatus)
!.............................................................................

      include 'netcdf.inc'
      integer mtropnum, mwndnum, manlevel, recnum, sigtlevel,
     +     sigwlevel,stanamelen, nlvl_out, nf_fid, nf_vid, nf_status
      integer nummand(recnum), nummwnd(recnum), numsigt(recnum),
     +     numsigw(recnum), numtrop(recnum), sondtyp(recnum),
     +     wmostanum(recnum)
      real htman( manlevel, recnum), htsigw( sigwlevel, recnum),
     +     prman( manlevel, recnum), prmaxw( mwndnum, recnum),
     +     prsigt( sigtlevel, recnum), prtrop( mtropnum, recnum),
     +     staelev(recnum), stalat(recnum), stalon(recnum), tdman(
     +     manlevel, recnum), tdsigt( sigtlevel, recnum), tdtrop(
     +     mtropnum, recnum), tpman( manlevel, recnum), tpsigt(
     +     sigtlevel, recnum), tptrop( mtropnum, recnum), wdman(
     +     manlevel, recnum), wdmaxw( mwndnum, recnum), wdsigw(
     +     sigwlevel, recnum), wdtrop( mtropnum, recnum), wsman(
     +     manlevel, recnum), wsmaxw( mwndnum, recnum), wssigw(
     +     sigwlevel, recnum), wstrop( mtropnum, recnum)
      double precision reltime(recnum), syntime(recnum)
      character*6 staname(recnum)
      character   stanamefile(stanamelen,recnum)

      real      prsigw                         (sigwlevel,recnum)
!..............................................................................

      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

      character*9 a9time_syn, a9time_release, a9time_raob, a9time_sys
      character*8 c8_obstype
      logical l_fill_ht

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting r_missing_data'
          return
      endif
      call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat_a,lon_a,  
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_domain_perimeter'
          return
      else
          write(6,*)' nsew perimeter is ',rnorth,south,east,west
      endif

!.............................................................................


      call read_netcdf(nf_fid, mtropnum, mwndnum, manlevel, recnum, 
     +     sigtlevel, sigwlevel, stanamelen, nummand, nummwnd, numsigt, 
     +     numsigw, 
     +     numtrop, sondtyp, wmostanum, htman, htsigw, prman, prmaxw, 
     +     prsigt, prtrop, staelev, stalat, stalon, tdman, tdsigt, 
     +     tdtrop, tpman, tpsigt, tptrop, wdman, wdmaxw, wdsigw, 
     +     wdtrop, wsman, wsmaxw, wssigw, wstrop, reltime, syntime, 
     +     stanamefile, staname)
c
c the netcdf variables are filled - your code goes here
c
!     write all raobs to laps snd file

      prsigw = r_missing_data

      r_nc_missing_data = 1e20

      n_snd = recnum

      do isnd = 1,n_snd

!         qc and write out the sounding
          i4time_raob = 0

          if(abs(syntime(isnd)) .lt. 1e10 .and.
     1       abs(syntime(isnd)) .gt.    0.      )then
              i4time_syn  = idint(syntime(isnd))+315619200
              i4time_raob = i4time_syn
          else
              i4time_syn = 0
          endif

          if(abs(reltime(isnd)) .lt. 1e10 .and.
     1       abs(reltime(isnd)) .gt.    0.      )then
              i4time_release = idint(reltime(isnd))+315619200

              i4time_diff = i4time_release - i4time_sys

              if(abs(i4time_diff) .gt. 20000)then
                  write(6,*)' warning: i4time_release is not '
     1                     ,'consistent with i4time_diff'
     1                     ,i4time_release,i4time_sys
              endif

!             correction for balloon rise time to mid-troposphere
              i4time_raob = i4time_release + 1800

          else
              i4time_release = 0
              i4time_diff = 0

          endif

          write(6,*)
          write(6,*)' raob #',isnd,i4time_sys,i4time_release,i4time_diff       
     1                         ,i4time_syn

          call make_fnam_lp(i4time_sys    , a9time_sys    , istatus)
          call make_fnam_lp(i4time_release, a9time_release, istatus)
          call make_fnam_lp(i4time_syn    , a9time_syn    , istatus)
          call make_fnam_lp(i4time_raob   , a9time_raob   , istatus)

          write(6,*)' times - sys/release/syn/raob: '
     1             ,a9time_sys,' ',a9time_release,' '
     1             ,a9time_syn,' ',a9time_raob

          if(stalat(isnd) .ge. r_nc_missing_data)then
              write(6,*)' missing first latitude',isnd
              goto 999
          endif

          if(stalon(isnd) .ge. r_nc_missing_data)then
              write(6,*)' missing first longitude',isnd
              goto 999
          endif

          if(stalat(isnd) .le. rnorth .and. stalat(isnd) .ge. south 
     1                                .and.      
     1       stalon(isnd) .ge. west   .and. stalon(isnd) .le. east      
     1                                                            )then       

!         if(.true.)then      ! for testing

              write(6,*)' raob is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' outside domain lat/lon perimeter - reject'
     1           ,stalat(isnd),stalon(isnd)
              goto 999
          endif

          if(i4time_raob .ne. 0)then ! test window
              if(i4time_raob .ge. i4time_raob_earliest .and.
     1           i4time_raob .le. i4time_raob_latest)then
                  write(6,*)' inside time window'
              else
                  write(6,*)' outside time window - reject'
                  goto 999
              endif
          endif

          c8_obstype = 'raob'

          call sort_and_write(i4time_sys,lun_out,l_fill_ht              ! i
     1                       ,recnum,isnd,r_missing_data,a9time_raob
     1                       ,c8_obstype
     1                       ,wmostanum,staname,stalat,stalon,staelev
     1                       ,nummand,htman,prman,tpman,tdman      
     1                       ,wdman,wsman
     1                       ,numsigt,prsigt,tpsigt,tdsigt
     1                       ,numsigw,prsigw,htsigw,wdsigw,wssigw
     1                       ,nlvl_out 
     1                       ,manlevel,sigtlevel,sigwlevel,istatus)

          go to 999

 998      write(6,*)' error writing out raob'

 999      continue

      enddo ! i

      return
      end


      subroutine sort_and_write(i4time_sys,lun_out,l_fill_ht            ! i
     1                       ,nrec,isnd,r_missing_data,a9time_raob
     1                       ,c8_obstype
     1                       ,wmostanum,staname,stalat,stalon,staelev
     1                       ,nummand,htman,prman,tpman,tdman      
     1                       ,wdman,wsman
     1                       ,numsigt,prsigt,tpsigt,tdsigt
     1                       ,numsigw,prsigw,htsigw,wdsigw,wssigw
     1                       ,nlvl_out
     1                       ,manlevel,sigtlevel,sigwlevel,istatus)

      integer     nlvl_out
      integer     manlevel
      integer     sigtlevel
      integer     sigwlevel
      integer   wmostanum                      (nrec)
      character*1 staname                        (6,nrec)
      real      stalat                         (nrec)
      real      stalon                         (nrec)
      real      staelev                        (nrec)

      integer   nummand                        (nrec)
      real      prman                          (manlevel ,nrec)
      real      htman                          (manlevel ,nrec)
      real      tpman                          (manlevel ,nrec)
      real      tdman                          (manlevel ,nrec) ! dwpt dprs
      real      wdman                          (manlevel ,nrec)
      real      wsman                          (manlevel ,nrec)

      real      prman_good                     (manlevel)
      real      htman_good                     (manlevel)
      real      tpman_good                     (manlevel)
      real      tdman_good                     (manlevel)       ! dwpt dprs
      real      wdman_good                     (manlevel)
      real      wsman_good                     (manlevel)

      integer   numsigt                        (nrec)
      real      prsigt                         (sigtlevel,nrec)
      real      tpsigt                         (sigtlevel,nrec)
      real      tdsigt                         (sigtlevel,nrec) ! dwpt dprs

      integer   numsigw                        (nrec)
      real      htsigw                         (sigwlevel,nrec)
      real      prsigw                         (sigwlevel,nrec)
      real      wdsigw                         (sigwlevel,nrec)
      real      wssigw                         (sigwlevel,nrec)

      integer   indx(nlvl_out)  
      character*9 a9time_out_sort                (nlvl_out)
      real      latout_sort                    (nlvl_out)
      real      lonout_sort                    (nlvl_out)
      real      prout                          (nlvl_out)
      real      prout_sort                     (nlvl_out)
      real      htout                          (nlvl_out)
      real      htout_sort                     (nlvl_out)
      real      tpout                          (nlvl_out)
      real      tpout_sort_c                   (nlvl_out)
      real      tpout_c_zman                   (nlvl_out)
      real      tdout                          (nlvl_out)  ! dewpoint depress
      real      tdout_sort_c                   (nlvl_out)  ! dewpoint deg c
      real      tdout_c_zman                   (nlvl_out)  ! dewpoint deg c
      real      wdout                          (nlvl_out)
      real      wdout_sort                     (nlvl_out)
      real      wsout                          (nlvl_out)
      real      wsout_sort                     (nlvl_out)

      real k_to_c

      character*9 a9time_raob
      character*8 c8_obstype
      character*132 c_line

      logical l_fill_ht, l_fill_ht_a(nlvl_out)

!     generate info for sorting/qc, write original mandatory data to log file
      write(6,*)
      n_good_levels = 0

      l_fill_ht_a = .false.

      write(6,*)' subroutine sort_and_write - initial mandatory data'       
      if(nummand(isnd) .le. manlevel)then
        do ilevel = 1,nummand(isnd)
          call check_nan(htman(ilevel,isnd),istat_nan)
          if(htman(ilevel,isnd) .lt. 90000. .and.
     1       htman(ilevel,isnd) .ge. staelev(isnd) .and.
     1       istat_nan .eq. 1                         )then ! valid height agl

            if(prman(ilevel,isnd) .le. prman(1,isnd))then ! pres is <= sfcp 
              n_good_levels = n_good_levels + 1
              write(6,*) htman(ilevel,isnd),prman(ilevel,isnd)
     1                  ,tpman(ilevel,isnd),tdman(ilevel,isnd)
     1                  ,wdman(ilevel,isnd),wsman(ilevel,isnd)

              htman_good(n_good_levels) = htman(ilevel,isnd)
              prman_good(n_good_levels) = prman(ilevel,isnd)
              tpman_good(n_good_levels) = tpman(ilevel,isnd)
              tdman_good(n_good_levels) = tdman(ilevel,isnd)
              wdman_good(n_good_levels) = wdman(ilevel,isnd)
              wsman_good(n_good_levels) = wsman(ilevel,isnd)
              indx(n_good_levels) = n_good_levels

            else
              write(6,*)' reject pres > sfcp ',ilevel,prman(ilevel,isnd)
     1                                               ,prman(1,isnd)
            endif

          endif
        enddo

      else
        write(6,*)' note: nummand(isnd) > manlevel'
     1                   ,nummand(isnd),manlevel      
      endif

      n_good_man = n_good_levels

!     bubble sort the mandatory levels by height 
!     (this may not be all that essential given that underground 
!     levels have been filtered out)

 300  iswitch = 0
      do i = 2,n_good_man
          if(htman_good(indx(i)) .lt. htman_good(indx(i-1)))then
              izz = indx(i-1)
              indx(i-1) = indx(i)
              indx(i) = izz
              iswitch = 1
          endif
      enddo

      if(iswitch .eq. 1)go to 300

      do i = 1,n_good_man
          ilevel = indx(i)
          htout(i) = htman_good(ilevel)
          prout(i) = prman_good(ilevel)
          tpout(i) = tpman_good(ilevel)
          tdout(i) = tdman_good(ilevel)
          wdout(i) = wdman_good(ilevel)
          wsout(i) = wsman_good(ilevel)
      enddo

!     generate sounding from man lvls for subsequent use in height integration
      n_good_zman = 0
      tdout_ref = 10.
      do i = 1,n_good_man
          if(abs(tpout(i)) .le. 500.)then ! good t value
              tpout_c_zman(i) = k_to_c(tpout(i))
              n_good_zman = n_good_zman + 1

              if(abs(tdout(i)) .le. 500.)then ! good td value
                  tdout_c_zman(i) = tpout_c_zman(i) - tdout(i)
!                 tdout_c_zman(i) = k_to_c(tdout(i))

                  tdout_ref = tdout(i)

              else ! generate approximate moisture value for height integration
                  tdout_c_zman(i) = tpout_c_zman(i) - 10.
!                 tdout_c_zman(i) = tpout_c_zman(i) - tdout_ref

              endif ! good td value
          endif ! good t value
      enddo ! i

      write(6,*)' subroutine sort_and_write - sig wind data'       
      if(numsigw(isnd) .le. sigwlevel)then
        do ilevel = 1,numsigw(isnd)
          call check_nan(htsigw(ilevel,isnd),istat_nan)
          if(htsigw(ilevel,isnd) .lt. 90000. .and.
     1       htsigw(ilevel,isnd) .ne. 0.     .and.
     1       istat_nan .eq. 1                      )then ! height is valid
              n_good_levels = n_good_levels + 1
              write(6,*) htsigw(ilevel,isnd),r_missing_data
     1                  ,r_missing_data,r_missing_data
     1                  ,wdsigw(ilevel,isnd),wssigw(ilevel,isnd)

              indx(n_good_levels) = n_good_levels
              htout(n_good_levels) = htsigw(ilevel,isnd)
              prout(n_good_levels) = r_missing_data
              tpout(n_good_levels) = r_missing_data
              tdout(n_good_levels) = r_missing_data
              wdout(n_good_levels) = wdsigw(ilevel,isnd)
              wsout(n_good_levels) = wssigw(ilevel,isnd)

          elseif(prsigw(ilevel,isnd) .lt. 2000. .and.
     1           prsigw(ilevel,isnd) .gt. 0.            )then ! pres is valid

!           attempt to calculate height based on good mandatory level data
            if(n_good_zman .gt. 0)then
                ht_calc = z(prsigw(ilevel,isnd),prman_good
     1                     ,tpout_c_zman,tdout_c_zman,n_good_zman)  

                ht_ref = htout(1)     

                if(nanf(ht_calc) .eq. 1 
     1         .or. ht_calc .gt. 99999. .or. ht_calc .lt. -1000.)then
                    ht_calc = -1.0        ! flag value for invalid height
                endif

                if(ht_calc .ne. -1.0)then ! valid height returned
                    ht_calc = ht_calc + ht_ref

                    n_good_levels = n_good_levels + 1
                    l_fill_ht_a(n_good_levels) = .true.

                    write(6,*) ht_calc,prsigw(ilevel,isnd)
     1                        ,tpsigt(ilevel,isnd),tdsigt(ilevel,isnd)
     1                        ,r_missing_data,r_missing_data

                    indx(n_good_levels) = n_good_levels
                    htout(n_good_levels) = ht_calc
                    prout(n_good_levels) = prsigw(ilevel,isnd)
                    tpout(n_good_levels) = r_missing_data
                    tdout(n_good_levels) = r_missing_data
                    wdout(n_good_levels) = wdsigw(ilevel,isnd)
                    wsout(n_good_levels) = wssigw(ilevel,isnd)
                endif ! valid height
            endif ! n_good_zman > 0
          endif ! htsigw/prsigw in bounds
        enddo
      else
        write(6,*)' note: numsigw(isnd) > sigwlevel'
     1                   ,numsigw(isnd),sigwlevel      
      endif

      write(6,*)' subroutine sort_and_write - sig t data'       
      if(numsigt(isnd) .le. sigtlevel)then
        do ilevel = 1,numsigt(isnd)
          if(prsigt(ilevel,isnd) .lt. 2000. .and.
     1       prsigt(ilevel,isnd) .gt. 0.            )then

!           attempt to calculate height based on good mandatory level data
            if(n_good_zman .gt. 0)then
                ht_calc = z(prsigt(ilevel,isnd),prman_good
     1                     ,tpout_c_zman,tdout_c_zman,n_good_zman)

                ht_ref = htout(1)

                if(nanf(ht_calc) .eq. 1 
     1         .or. ht_calc .gt. 99999. .or. ht_calc .lt. -1000.)then
                    ht_calc = -1.0        ! flag value for invalid height
                endif

                if(ht_calc .eq. -1.0)then ! invalid height returned

!                   we may be above the highest mandatory level, so calculate 
!                   layer thickness/ht with reference to previous sigt level

                    if(ilevel .gt. 1)then
                        p1 = prout(n_good_levels)
                        p2 = prsigt(ilevel,isnd)
                        t1 = tpout(n_good_levels)
                        t2 = tpsigt(ilevel,isnd)
                        td1 = tpout(n_good_levels)-tdout(n_good_levels)
                        td2 = tpsigt(ilevel,isnd) -tdsigt(ilevel,isnd)
                        h1 = htout(n_good_levels)
                        call calc_new_ht(p1,p2,t1,t2,td1,td2       ! i
     1                                  ,r_missing_data            ! i
     1                                  ,h1,h2)                    ! i/o
                        if(h2 .ne. r_missing_data)then
                            ht_ref = 0.
                            ht_calc = h2
                        endif
                    endif

                endif

                if(ht_calc .ne. -1.0)then ! valid height returned
                    ht_calc = ht_calc + ht_ref

                    n_good_levels = n_good_levels + 1
                    l_fill_ht_a(n_good_levels) = .true.

                    write(6,*) ht_calc,prsigt(ilevel,isnd)
     1                        ,tpsigt(ilevel,isnd),tdsigt(ilevel,isnd)
     1                        ,r_missing_data,r_missing_data

                    indx(n_good_levels) = n_good_levels
                    htout(n_good_levels) = ht_calc
                    prout(n_good_levels) = prsigt(ilevel,isnd)
                    tpout(n_good_levels) = tpsigt(ilevel,isnd)
                    tdout(n_good_levels) = tdsigt(ilevel,isnd)
                    wdout(n_good_levels) = r_missing_data
                    wsout(n_good_levels) = r_missing_data
                endif ! valid height

            endif ! n_good_zman > 0
          endif ! prsigt in bounds
        enddo ! ilevel
      else
        write(6,*)' note: numsigt(isnd) > sigtlevel'
     1                   ,numsigt(isnd),sigtlevel      
      endif

!     bubble sort all the levels by height
      do i = 1,n_good_levels
          indx(i) = i
      enddo ! i

 400  iswitch = 0
      do i = 2,n_good_levels
          if(htout(indx(i)) .lt. htout(indx(i-1)))then
              izz = indx(i-1)
              indx(i-1) = indx(i)
              indx(i) = izz
              iswitch = 1
          endif
      enddo

      if(iswitch .eq. 1)go to 400
      
!     detect and remove duplicate levels
      do ipass = 1,2
        i = 2
        do while(i .le. n_good_levels)
          idupe = 0
          if(htout(indx(i)) .eq. htout(indx(i-1)))then          
              idupe = i
              write(6,*)' remove duplicate ht level '
     1                  ,idupe,htout(indx(idupe))

          elseif( prout(indx(i)) .eq. prout(indx(i-1))
     1     .and.  prout(indx(i)) .ne. r_missing_data    )then          
              idupe_pr = 0

!             we would like to retain the mandatory report that also has winds
              if(wdout(indx(i)) .ne. r_missing_data .and.
     1           tpout(indx(i)) .ne. r_missing_data)then
                  idupe_pr = i-1 ! previous sig t level can be removed
              elseif(wdout(indx(i-1)) .ne. r_missing_data .and.
     1               tpout(indx(i-1)) .ne. r_missing_data)then
                  idupe_pr = i   ! current sig t level can be removed
              elseif(wdout(indx(i)) .ne. r_missing_data)then
                  idupe_pr = i-1 ! previous sig t level can be removed
              elseif(wdout(indx(i-1)) .ne. r_missing_data)then
                  idupe_pr = i   ! current sig t level can be removed
              endif
            
              write(6,*)' remove duplicate pr level'
     1                  ,i,idupe_pr,prout(indx(i))

              idupe = idupe_pr

          endif ! duplicate level

          if(idupe .gt. 0)then
              do j = idupe,n_good_levels-1
                l_fill_ht_a(indx(j)) = l_fill_ht_a(indx(j+1))
                htout(indx(j)) = htout(indx(j+1))
                prout(indx(j)) = prout(indx(j+1))
                tpout(indx(j)) = tpout(indx(j+1))
                tdout(indx(j)) = tdout(indx(j+1))
                wdout(indx(j)) = wdout(indx(j+1))
                wsout(indx(j)) = wsout(indx(j+1))
              enddo ! j
              n_good_levels = n_good_levels - 1
          endif
          i = i+1
        enddo ! i
      enddo ! ipass

!     qc and convert units, t and td are converted to deg c
      do i = 1,n_good_levels
          ilevel = indx(i)

          htout_sort(i) = htout(ilevel)
          if((.not. l_fill_ht) .and. l_fill_ht_a(ilevel))then
              htout_sort(i) = r_missing_data
          endif 

          prout_sort(i) = prout(ilevel)

          if(tpout(ilevel) .eq. 99999. .or.
     1       tpout(ilevel) .eq. r_missing_data     )then
              tpout_sort_c(i) = r_missing_data
          else
              tpout_sort_c(i) = k_to_c(tpout(ilevel))
          endif

          if(tdout(ilevel)   .eq. 99999. .or. 
     1       tpout_sort_c(i) .eq. r_missing_data)then       
              tdout_sort_c(i) = r_missing_data
          else
              tdout_sort_c(i) = k_to_c(tpout(ilevel)) 
     1                        - tdout(ilevel)      
          endif

          if(abs(wdout(ilevel)) .ge. 99999. .or.
     1       abs(wsout(ilevel)) .ge. 99999.)then
              wdout(ilevel) = r_missing_data
              wsout(ilevel) = r_missing_data
          endif

          wdout_sort(i) = wdout(ilevel)
          wsout_sort(i) = wsout(ilevel)

      enddo ! i

      latout_sort = stalat(isnd)      ! assign entire array for this sounding
      lonout_sort = stalon(isnd)      ! assign entire array for this sounding
      a9time_out_sort = a9time_raob   ! assign entire array for this sounding

      call open_ext(lun_out,i4time_sys,'snd',istatus)

      maxlvl = nlvl_out

      call write_snd  (lun_out                                    ! i
     1                ,1,maxlvl,1                                 ! i
     1                ,wmostanum(isnd)                            ! i
     1                ,latout_sort,lonout_sort,staelev(isnd)      ! i
     1                ,staname(1,isnd),a9time_out_sort,c8_obstype ! i
     1                ,n_good_levels                              ! i
     1                ,htout_sort                                 ! i
     1                ,prout_sort                                 ! i
     1                ,tpout_sort_c                               ! i
     1                ,tdout_sort_c                               ! i
     1                ,wdout_sort                                 ! i
     1                ,wsout_sort                                 ! i
     1                ,istatus)                                   ! o

      return
      end

!.............................................................................
c
c  subroutine to read the file "raob data : selected by ob time : time range from 887191200 to 887202000" 
c
      subroutine read_netcdf(nf_fid, mtropnum, mwndnum, manlevel, 
     +     recnum, sigtlevel, sigwlevel, stanamelen, nummand, nummwnd, 
     +     numsigt, numsigw, 
     +     numtrop, sondtyp, wmostanum, htman, htsigw, 
     +     prman, prmaxw, prsigt, prtrop, staelev, stalat, stalon, 
     +     tdman, tdsigt, tdtrop, tpman, tpsigt, tptrop, wdman, 
     +     wdmaxw, wdsigw, wdtrop, wsman, wsmaxw, wssigw, wstrop, 
     +     reltime, syntime, stanamefile, staname)
c
      include 'netcdf.inc'
      integer mtropnum, mwndnum, manlevel, recnum, sigtlevel, 
     +     sigwlevel,stanamelen, nf_fid, nf_vid, nf_status
      integer nummand(recnum), nummwnd(recnum), numsigt(recnum),
     +     numsigw(recnum), numtrop(recnum), sondtyp(recnum),
     +     wmostanum(recnum)
      real htman( manlevel, recnum), htsigw( sigwlevel, recnum),
     +     prman( manlevel, recnum), prmaxw( mwndnum, recnum),
     +     prsigt( sigtlevel, recnum), prtrop( mtropnum, recnum),
     +     staelev(recnum), stalat(recnum), stalon(recnum), tdman(
     +     manlevel, recnum), tdsigt( sigtlevel, recnum), tdtrop(
     +     mtropnum, recnum), tpman( manlevel, recnum), tpsigt(
     +     sigtlevel, recnum), tptrop( mtropnum, recnum), wdman(
     +     manlevel, recnum), wdmaxw( mwndnum, recnum), wdsigw(
     +     sigwlevel, recnum), wdtrop( mtropnum, recnum), wsman(
     +     manlevel, recnum), wsmaxw( mwndnum, recnum), wssigw(
     +     sigwlevel, recnum), wstrop( mtropnum, recnum)
      double precision reltime(recnum), syntime(recnum)
      character*6 staname(recnum), name
      character   stanamefile(stanamelen,recnum)


c   variables of type real
c
c     variable        netcdf long name
c      htman        "geopotential - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'htman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var htman'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,htman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var htman'
      endif
c
c     variable        netcdf long name
c      htsigw       "geopotential - significant level wrt w"
c
        nf_status = nf_inq_varid(nf_fid,'htsigw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var htsigw'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,htsigw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var htsigw'
      endif
c
c     variable        netcdf long name
c      prman        "pressure - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'prman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prman'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,prman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prman'
      endif
c
c     variable        netcdf long name
c      prmaxw       "pressure - maximum wind level"
c
        nf_status = nf_inq_varid(nf_fid,'prmaxw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prmaxw'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,prmaxw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prmaxw'
      endif
c
c     variable        netcdf long name
c      prsigt       "pressure - significant level wrt t"
c
        nf_status = nf_inq_varid(nf_fid,'prsigt',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prsigt'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,prsigt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prsigt'
      endif
c
c     variable        netcdf long name
c      prtrop       "pressure - tropopause level"
c
        nf_status = nf_inq_varid(nf_fid,'prtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prtrop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,prtrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var prtrop'
      endif
c
c     variable        netcdf long name
c      staelev      "station elevation"
c
        nf_status = nf_inq_varid(nf_fid,'staelev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staelev'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,staelev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staelev'
      endif
c
c     variable        netcdf long name
c      stalat       "station latitude"
c
        nf_status = nf_inq_varid(nf_fid,'stalat',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalat'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stalat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalat'
      endif
c
c     variable        netcdf long name
c      stalon       "station longitude"
c
        nf_status = nf_inq_varid(nf_fid,'stalon',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalon'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stalon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalon'
      endif
c
c     variable        netcdf long name
c      tdman        "dew point depression - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'tdman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tdman'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,tdman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tdman'
      endif
c
c     variable        netcdf long name
c      tdsigt       "dew point depression - significant level wrt t"
c
        nf_status = nf_inq_varid(nf_fid,'tdsigt',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tdsigt'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,tdsigt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tdsigt'
      endif
c
c     variable        netcdf long name
c      tdtrop       "dew point depression - tropopause level"
c
        nf_status = nf_inq_varid(nf_fid,'tdtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tdtrop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,tdtrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tdtrop'
      endif
c
c     variable        netcdf long name
c      tpman        "temperature - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'tpman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpman'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,tpman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpman'
      endif
c
c     variable        netcdf long name
c      tpsigt       "temperature - significant level wrt t"
c
        nf_status = nf_inq_varid(nf_fid,'tpsigt',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpsigt'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,tpsigt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpsigt'
      endif
c
c     variable        netcdf long name
c      tptrop       "temperature - tropopause level"
c
        nf_status = nf_inq_varid(nf_fid,'tptrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tptrop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,tptrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tptrop'
      endif
c
c     variable        netcdf long name
c      wdman        "wind direction - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'wdman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdman'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wdman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdman'
      endif
c
c     variable        netcdf long name
c      wdmaxw       "wind direction - maximum wind level"
c
        nf_status = nf_inq_varid(nf_fid,'wdmaxw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdmaxw'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wdmaxw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdmaxw'
      endif
c
c     variable        netcdf long name
c      wdsigw       "wind direction - significant level wrt w"
c
        nf_status = nf_inq_varid(nf_fid,'wdsigw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdsigw'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wdsigw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdsigw'
      endif
c
c     variable        netcdf long name
c      wdtrop       "wind direction - tropopause level"
c
        nf_status = nf_inq_varid(nf_fid,'wdtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdtrop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wdtrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wdtrop'
      endif
c
c     variable        netcdf long name
c      wsman        "wind speed - mandatory level"
c
        nf_status = nf_inq_varid(nf_fid,'wsman',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wsman'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wsman)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wsman'
      endif
c
c     variable        netcdf long name
c      wsmaxw       "wind speed - maximum wind level"
c
        nf_status = nf_inq_varid(nf_fid,'wsmaxw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wsmaxw'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wsmaxw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wsmaxw'
      endif
c
c     variable        netcdf long name
c      wssigw       "wind speed - significant level wrt w"
c
        nf_status = nf_inq_varid(nf_fid,'wssigw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wssigw'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wssigw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wssigw'
      endif
c
c     variable        netcdf long name
c      wstrop       "wind speed - tropopause level"
c
        nf_status = nf_inq_varid(nf_fid,'wstrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wstrop'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wstrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wstrop'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      nummand      "number of mandatory levels"
c
        nf_status = nf_inq_varid(nf_fid,'nummand',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nummand'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,nummand)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nummand'
      endif
c
c     variable        netcdf long name
c      nummwnd      "number of maximum wind levels"
c
        nf_status = nf_inq_varid(nf_fid,'nummwnd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nummwnd'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,nummwnd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nummwnd'
      endif
c
c     variable        netcdf long name
c      numsigt      "number of significant levels wrt t"
c
        nf_status = nf_inq_varid(nf_fid,'numsigt',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numsigt'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,numsigt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numsigt'
      endif
c
c     variable        netcdf long name
c      numsigw      "number of significant levels wrt w"
c
        nf_status = nf_inq_varid(nf_fid,'numsigw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numsigw'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,numsigw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numsigw'
      endif
c
c     variable        netcdf long name
c      numtrop      "number of tropopause levels"
c
        nf_status = nf_inq_varid(nf_fid,'numtrop',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numtrop'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,numtrop)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numtrop'
      endif
c
c     variable        netcdf long name
c      sondtyp      "instrument type"
c
        nf_status = nf_inq_varid(nf_fid,'sondtyp',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sondtyp'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,sondtyp)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sondtyp'
      endif
c
c     variable        netcdf long name
c      wmostanum    "wmo station number"
c
        nf_status = nf_inq_varid(nf_fid,'wmostanum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wmostanum'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,wmostanum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wmostanum'
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      reltime      "sounding release time"
c
        nf_status = nf_inq_varid(nf_fid,'reltime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reltime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,reltime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reltime'
      endif
c
c     variable        netcdf long name
c      syntime      "synoptic time"
c
        nf_status = nf_inq_varid(nf_fid,'syntime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var syntime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,syntime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var syntime'
      endif


c   variables of type char
c
c     variable        netcdf long name
c      staname      "station identifier"
c
      nf_status = nf_inq_varid(nf_fid,'staname',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staname'
      endif

      if (stanamelen .eq. 6) then  !read directly into staname variable
        nf_status = nf_get_var_text(nf_fid,nf_vid,staname)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var staname'
        endif
      else
        nf_status = nf_get_var_text(nf_fid,nf_vid,stanamefile)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var staname'
        endif

        do i = 1, recnum
          do j = 1, 6
            name(j:j) = stanamefile(j,i)
          enddo
          call filter_string(name)
          staname(i) = name
        enddo
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end

      subroutine calc_new_ht(p1,p2,t1,t2,td1,td2       ! i
     1                      ,r_missing_data            ! i
     1                      ,h1,h2)                    ! i/o

!     use the hypsometric equation to calculate the height at the top of
!     a layer

      real pr_z(2)
      real tp_z(2)
      real td_z(2)

      h2 = r_missing_data

!     apply qc and fill arrays
      if(p1 .lt. 2000. .and. p1 .gt. 0.)then
          pr_z(1) = p1
      else
          return
      endif

      if(p2 .lt. 2000. .and. p2 .gt. 0.)then
          pr_z(2) = p2
      else
          return
      endif

      if(t1 .gt. 100. .and. t1 .lt. 400.)then
          tp_z(1) = t1
      else
          return
      endif

      if(t2 .gt. 100. .and. t2 .lt. 400.)then
          tp_z(2) = t2
      else
          return
      endif

      if(td1 .gt. 100. .and. td1 .lt. 400.)then
          td_z(1) = td1
      else
          return
      endif

      if(td2 .gt. 100. .and. td2 .lt. 400.)then
          td_z(2) = td2
      else
          return
      endif

      if(h1 .gt. 90000. .or. h1 .le. -1000.)then
          return
      endif

      thk = z(p2,pr_z,tp_z,td_z,2)

      if(thk .ne. -1.0)then
          h2 = h1 + thk
      endif

      return
      end

