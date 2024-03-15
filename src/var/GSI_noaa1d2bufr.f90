       subroutine gsi_noaa1d2bufr
!***********************************************************************
! program name : gsi_noaa1d2bufr(gsi_noaa1d2bufr.f90)
!
! description: to read in noaa-nn satellite 1d data (amsu-a,amsu-b,airs3,airs4,or mhs)
!              and transform to ncep bufr format.
!
! environmental variable:
!      polarsat : directory of polar satellite input data.
!
! files:
!      satellite.bufrtable : satellite bufr table (input file)
!      amsual1d_noaann_yyyymmdd_hhmm*.l1d: noaa-nn amsu-a 1d format file. (input file)
!      amsuabufrnn_yyyymmdd_hhmm or amsuabufr: noaa-nn amsu-a ncep bufr data file. (output file)
!
!  or  amsubl1d_noaann_yyyymmdd_hhmm*.l1d: noaa-nn amsu-b 1d format file. (input file)
!      amsubbufrnn_yyyymmdd_hhmm or amsubbufr: noaa-nn amsu-b ncep bufr data file. (output file)
!
!  or  mhsl1d_noaa18_yyyymmdd_hhmm*.l1d: noaa-18 mhs 1d format file. (input file)
!      amsubbufr18_yyyymmdd_hhmm or amsubbufr: noaa-18 mhs ncep bufr data file. (output file)
!
!  or  hirsl1d_noaann_yyyymmdd_hhmm*.l1d: noaa-nn hirs 1d format file. (input file)
!      hirs3bufrnn_yyyymmdd_hhmm or hirs3bufr: noaa-nn hirs ncep bufr data file. (output file)
!      (yyyymmdd_hhmm: year, month, day, hour and minute; nn is satellite id.)
!
! 1-d format structure:
!          general information:        4*24 bytes
!          1-d data structure: depending on different instruments
! amsua-1d data structure:
!          header record: 4*606 bytes
!          spare record: 2088 bytes
!          data record scan line 1: 4*1152 bytes
!          data record scan line 2: 4*1152 bytes
!          .....
!          data record scan line n: 4*1152 bytes
! amsub-1d or mhs-d data structure:
!          header record: 4*2278 bytes
!          spare record: 3080 bytes
!          data record scan line 1: 4*3072 bytes
!          data record scan line 2: 4*3072 bytes
!          .....
!          data record scan line n: 4*3072 bytes
! hirs-1d data structure:
!          header record: 4*3823 bytes
!          spare record: 484 bytes
!          data record scan line 1: 4*3968 bytes
!          data record scan line 2: 4*3968 bytes
!          .....
!          data record scan line n: 4*3968 bytes
!
! ncep noaa satellite bufr table-a:
!     nc021021: msg type 021-021 processed hirs-2 1b data (noaa 14)
!     nc021022: msg type 021-022 processed msu-2  1b data (noaa 14)
!     nc021023: msg type 021-023 processed amsu-a 1b data (noaa 15-18)
!     nc021024: msg type 021-024 processed amsu-b 1b data (noaa 15-17)
!     nc021025: msg type 021-025 processed hirs-3 1b data (noaa 15-17)
!     nc021027: msg type 021-027 processed mhs    1b data (noaa 18)
!     nc021028: msg type 021-028 processed hirs-4 1b data (noaa 18)
!
! ncep bufr format:
!   amsu-a, amsu-b, mhs:
!-----------------------------------------------------------------------
!  nc021sss  | year  mnth  days  hour  minu  seco  clat  clon  said
!  nc021sss  | siid  fovn  lsql  saza  soza  hols  hmsl  solazi  bearaz
!  nc021sss  | "brit"xx
!  brit      | chnm  tmbr
!-----------------------------------------------------------------------
!  where xx=15 for amsu-a,  =5 for amsu-b/mhs
!
!   hirs-3, hirs4:
!-----------------------------------------------------------------------
!  nc021sss  | year  mnth  days  hour  minu  seco  clat  clon  said
!  nc021sss  | siid  fovn  lsql  saza  soza  hols  hmsl  solazi  bearaz
!  nc021sss  | "brit"20
!  brit      | chnm  tmbr
!-----------------------------------------------------------------------
!
!   hirs-2, msu:
!-----------------------------------------------------------------------
!  nc021sss  | year  mnth  days  hour  minu  seco  clat  clon  said
!  nc021sss  | siid  fovn  lsql  saza  soza  hols  hmsl
!  nc021sss  | "brit"xx
!  brit      | chnm  tmbr
!-----------------------------------------------------------------------
!  where xx=20 for hirs-2,  =4 for msu
!  for siid:
!      =605, for hirs-2;    =606, for hirs-3;    =607, for hirs-4;
!      =623, for msu;       =570, for amsu-a:    =574, for amsu-b;
!      =203, for mhs
!
! called function:
!    rioreadbyte(getiofile.c),
!    getnoaafile,  dc_gen_inf,  dc_hirs,  dc_amsua,  dc_mhs
!    openbf, openmb, ufbseq, writsb, closbf (ncep bufr routine)
!
! date :
!   original     -- may  02, 2007 (shiow-ming deng)
!***********************************************************************

          implicit none

          character*80 path_satdat
          integer ipath, numfile, numcha(900)
          character*180 file(900)

          character*180 dir, sat_table
          integer len, len_sat_table
          logical sattab

          character filein*80, fileout*80
          integer isize, machine
          parameter(isize=90000000, machine=1)
          integer ifile, ier, isize1
          byte buf(isize), genbyte(96)
          byte scanbyte(15872)

          integer i, j, k, kk, l, ll
          integer sateid, instrument, nline
          integer year, month, day, hour, minute, icheck
          integer siid, nbyte, xtrack, nchannel

          integer iscan, isqc
          real second, hmsl

          integer ifqc(90), isty(90)
          real rlat(90), rlon(90), ht(90), stzn(90), sozn(90)
          real staz(90), soaz(90)
          real brtmp_hirs(20, 56), brtmp_amsua(15, 30), brtmp_mhs(5, 90)

          integer lnbufr, idate, iret, n
          character subset*8
          integer ndat
          parameter(ndat=100)
          real*8 bufrf(ndat)

!-----------------------------------------------------------------------
!c  to get directory of satellite input data.
!c  environmental variable: polarsat

          call getenv('polarsat', path_satdat)
          ipath = index(path_satdat, ' ')
          if (ipath .le. 1) then
             print *, 'no setting environmental variable: polarsat'
             return
          end if
          print *, 'directory of input data: ', path_satdat(1:ipath - 1)
          path_satdat(ipath:ipath) = char(0)

!-----------------------------------------------------------------------
!c  to get noaa satellite 1d format input files.

          call getnoaafile(path_satdat, numfile, numcha, file, ier)
          if ((ier .ne. 0) .or. (numfile .lt. 1)) then
             print *, 'no any noaa satellite hdf format file.'
             print *, 'directory of input data: ', path_satdat(1:ipath - 1)
             return
          end if

!-----------------------------------------------------------------------
!c  to get satellite bufr table: satellite.bufrtable

          call getenv('laps_data_root', dir)
          len = index(dir, ' ')
          if (len .le. 1) then
             print *, 'cannot get laps_data_root directory.'
             return
          end if
          sat_table(1:len - 1) = dir(1:len - 1)
          sat_table(len:len + 23) = '/log/satellite.bufrtable'
          len_sat_table = len + 23
          inquire (file=sat_table(1:len_sat_table), exist=sattab)
          if (.not. sattab) then
             print *, 'satellite bufr table is not exist.'
             print *, 'satellite bufr table: ', sat_table(1:len_sat_table)
             return
          end if
          sat_table(len_sat_table + 1:len_sat_table + 1) = char(0)
          print *, 'satellite bufr table: ', sat_table(1:len_sat_table)

          do 100 l = 1, numfile

             filein(1:ipath - 1) = path_satdat(1:ipath - 1)
             filein(ipath:ipath) = '/'
             ll = numcha(l)
             filein(ipath + 1:ipath + ll) = file(l) (1:ll)
             filein(ipath + ll + 1:ipath + ll + 1) = char(0)
             print *, ' '
             print *, 'process file: ', filein(1:ipath + ll)
             print *, ' '

!-----------------------------------------------------------------------
!c  to read in noaa-nn 1d data bytes.

             isize1 = isize
             call rioreadfile(filein, buf, isize1, ier)
             if (ier .ne. 0) then
                print *, 'cannot read in byte data.'
                print *, 'file: ', filein(1:ipath + ll)
                go to 100
             end if
             if (isize1 .lt. 4608) then
                print *, 'this is not noaa-nn level 1d file.'
                print *, 'isize1: ', isize1
                print *, 'file: ', filein(1:ipath + ll)
                go to 100
             end if

             do k = 1, ndat
                bufrf(k) = 10.0e10
             end do

!-----------------------------------------------------------------------
!c  to decode header record.

             do i = 1, 96
                genbyte(i) = buf(i)
             end do
             call dc_gen_inf(genbyte, machine, sateid, instrument, nline &
                             , year, month, day, hour, minute, ier)
             if (ier .ne. 0) then
                print *, 'cannot decode 1d general information.'
                print *, 'input file may be not a noaa-nn 1d data.'
                print *, 'input file: ', filein(1:ipath + ll)
                go to 100
             end if
             if (instrument .eq. 10) icheck = 4608*(nline + 1) !amsu-a
             if (instrument .eq. 11) icheck = 12288*(nline + 1) !amsu-b
             if (instrument .eq. 12) icheck = 12288*(nline + 1) !mhs
             if (instrument .eq. 5) icheck = 15872*(nline + 1) !hirs
             if (icheck .ne. isize1) then
                print *, 'the size of input 1d data has question.'
                print *, 'nline: ', nline, ' icheck: ', icheck
                print *, 'size of read in data: ', isize1, ' (bytes)'
                print *, 'input file: ', filein(1:ipath + ll)
                go to 100
             end if

!-----------------------------------------------------------------------
!c  to get output bufr format file name.

             if (instrument .eq. 5) then  !hirs
                write (fileout, '(8x,i2.2,1x,i4.4,2i2.2,1x,2i2.2)') sateid, year &
                   , month, day, hour, minute
                fileout(1:8) = 'hir3bufr'
                fileout(11:11) = '_'
                subset(1:8) = 'nc021025'
                siid = 606
                if (sateid .ge. 18) then
                   subset(1:8) = 'nc021028'
                   siid = 607
                end if
                if (sateid .le. 15) then
                   subset(1:8) = 'nc021021'
                   siid = 605
                   fileout(4:4) = '2'
                end if
                fileout(20:20) = '_'
                fileout(25:25) = char(0)
                nbyte = 15872
                xtrack = 56
                nchannel = 20
             end if
             if (instrument .eq. 10) then  !amsu-a
                write (fileout, '(9x,i2.2,1x,i4.4,2i2.2,1x,2i2.2)') sateid, year &
                   , month, day, hour, minute
                fileout(1:9) = 'amsuabufr'
                fileout(12:12) = '_'
                fileout(21:21) = '_'
                fileout(26:26) = char(0)
                siid = 570
                subset(1:8) = 'nc021023'
                nbyte = 4608
                xtrack = 30
                nchannel = 15
             end if
             if (instrument .eq. 11) then  !amsu-b
                write (fileout, '(9x,i2.2,1x,i4.4,2i2.2,1x,2i2.2)') sateid, year &
                   , month, day, hour, minute
                fileout(1:9) = 'amsubbufr'
                fileout(12:12) = '_'
                fileout(21:21) = '_'
                fileout(26:26) = char(0)
                siid = 574
                subset(1:8) = 'nc021024'
                nbyte = 12288
                xtrack = 90
                nchannel = 5
             end if
             if (instrument .eq. 12) then  !mhs
                write (fileout, '(9x,i2.2,1x,i4.4,2i2.2,1x,2i2.2)') sateid, year &
                   , month, day, hour, minute
                fileout(1:9) = 'amsubbufr'
                fileout(12:12) = '_'
                fileout(21:21) = '_'
                fileout(26:26) = char(0)
                siid = 203
                subset(1:8) = 'nc021027'
                nbyte = 12288
                xtrack = 90
                nchannel = 5
             end if

!-----------------------------------------------------------------------
!c  to write out bufr.

             print *, 'output bufr file name: ', fileout
             print *, ' '

             open (11, file=fileout, form='unformatted')
             open (12, file=sat_table)
             lnbufr = 11
             call openbf(lnbufr, 'out', 12)

             n = 0
             do k = 1, nline

!-----------------------------------------------------------------------
!c  to decode scan line data.

                do i = 1, nbyte
                   j = nbyte + nbyte*(k - 1) + i
                   scanbyte(i) = buf(j)
                end do

!-----------------------------------------------------------------------
!c  to decode hirs line data.

                if (instrument .eq. 5) then
                   call dc_hirs(nbyte, scanbyte, machine, xtrack, iscan, isqc &
                                , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                                , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp_hirs, ier)
                end if

!-----------------------------------------------------------------------
!c  to decode amsu-a line data.

                if (instrument .eq. 10) then
                   call dc_amsua(nbyte, scanbyte, machine, xtrack, iscan, isqc &
                                 , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                                 , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp_amsua, ier)
                end if

!-----------------------------------------------------------------------
!c  to decode amsu-b or mhs line data.

                if ((instrument .eq. 11) .or. (instrument .eq. 12)) then
                   call dc_mhs(nbyte, scanbyte, machine, xtrack, iscan, isqc &
                               , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                               , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp_mhs, ier)
                end if

!-----------------------------------------------------------------------
!c  to write data into bufr

                bufrf(1) = year
                bufrf(2) = month
                bufrf(3) = day
                bufrf(4) = hour
                bufrf(5) = minute
                bufrf(6) = second
                bufrf(9) = sateid + 191
                bufrf(10) = siid
                bufrf(16) = hmsl
                idate = year*1000000 + month*10000 + day*100 + hour
                do i = 1, xtrack
                   if (ifqc(i) .eq. 0) then
                      bufrf(7) = rlat(i)
                      bufrf(8) = rlon(i)
                      bufrf(11) = i
                      bufrf(12) = isty(i)
                      bufrf(13) = stzn(i)
                      bufrf(14) = sozn(i)
                      bufrf(15) = ht(i)
                      bufrf(17) = soaz(i)
                      bufrf(18) = staz(i)
                      kk = 18
                      if (instrument .eq. 5) then
                         do j = 1, nchannel
                            kk = kk + 1
                            bufrf(kk) = j
                            kk = kk + 1
                            bufrf(kk) = brtmp_hirs(j, i)
                         end do
                      end if
                      if (instrument .eq. 10) then
                         do j = 1, nchannel
                            kk = kk + 1
                            bufrf(kk) = j
                            kk = kk + 1
                            bufrf(kk) = brtmp_amsua(j, i)
                         end do
                      end if
                      if ((instrument .eq. 11) .or. (instrument .eq. 12)) then
                         do j = 1, nchannel
                            kk = kk + 1
                            bufrf(kk) = j
                            kk = kk + 1
                            bufrf(kk) = brtmp_mhs(j, i)
                         end do
                      end if
                      call openmb(lnbufr, subset, idate)
                      call ufbseq(lnbufr, bufrf, ndat, 1, iret, subset)
                      call writsb(lnbufr)
                      n = n + 1
                   end if
                end do
             end do

             call closbf(lnbufr)
             close (lnbufr)
             close (12)
             print *, 'number of output bufr: ', n
             print *, ' '

!-----------------------------------------------------------------------
!c  to write out bufr again.

             if (instrument .eq. 5) then  !hirs
                fileout(1:8) = 'hir3bufr'
                fileout(9:9) = char(0)
             end if
             if (instrument .eq. 10) then  !amsu-a
                fileout(1:9) = 'amsuabufr'
                fileout(10:10) = char(0)
             end if
             if ((instrument .eq. 11) .or. (instrument .eq. 12)) then  !amsu-b
                fileout(1:9) = 'amsubbufr'
                fileout(10:10) = char(0)
             end if
             print *, 'output bufr file name: ', fileout
             print *, ' '

             open (11, file=fileout, form='unformatted')
             open (12, file=sat_table)
             lnbufr = 11
             call openbf(lnbufr, 'out', 12)

             n = 0
             do k = 1, nline

!-----------------------------------------------------------------------
!c  to decode scan line data.

                do i = 1, nbyte
                   j = nbyte + nbyte*(k - 1) + i
                   scanbyte(i) = buf(j)
                end do

!-----------------------------------------------------------------------
!c  to decode hirs line data.

                if (instrument .eq. 5) then
                   call dc_hirs(nbyte, scanbyte, machine, xtrack, iscan, isqc &
                                , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                                , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp_hirs, ier)
                end if

!-----------------------------------------------------------------------
!c  to decode amsu-a line data.

                if (instrument .eq. 10) then
                   call dc_amsua(nbyte, scanbyte, machine, xtrack, iscan, isqc &
                                 , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                                 , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp_amsua, ier)
                end if

!-----------------------------------------------------------------------
!c  to decode amsu-b or mhs line data.

                if ((instrument .eq. 11) .or. (instrument .eq. 12)) then
                   call dc_mhs(nbyte, scanbyte, machine, xtrack, iscan, isqc &
                               , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                               , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp_mhs, ier)
                end if

!-----------------------------------------------------------------------
!c  to write data into bufr

                bufrf(1) = year
                bufrf(2) = month
                bufrf(3) = day
                bufrf(4) = hour
                bufrf(5) = minute
                bufrf(6) = second
                bufrf(9) = sateid + 191
                bufrf(10) = siid
                bufrf(16) = hmsl
                idate = year*1000000 + month*10000 + day*100 + hour
                do i = 1, xtrack
                   if (ifqc(i) .eq. 0) then
                      bufrf(7) = rlat(i)
                      bufrf(8) = rlon(i)
                      bufrf(11) = i
                      bufrf(12) = isty(i)
                      bufrf(13) = stzn(i)
                      bufrf(14) = sozn(i)
                      bufrf(15) = ht(i)
                      bufrf(17) = soaz(i)
                      bufrf(18) = staz(i)
                      kk = 18
                      if (instrument .eq. 5) then
                         do j = 1, nchannel
                            kk = kk + 1
                            bufrf(kk) = j
                            kk = kk + 1
                            bufrf(kk) = brtmp_hirs(j, i)
                         end do
                      end if
                      if (instrument .eq. 10) then
                         do j = 1, nchannel
                            kk = kk + 1
                            bufrf(kk) = j
                            kk = kk + 1
                            bufrf(kk) = brtmp_amsua(j, i)
                         end do
                      end if
                      if ((instrument .eq. 11) .or. (instrument .eq. 12)) then
                         do j = 1, nchannel
                            kk = kk + 1
                            bufrf(kk) = j
                            kk = kk + 1
                            bufrf(kk) = brtmp_mhs(j, i)
                         end do
                      end if
                      call openmb(lnbufr, subset, idate)
                      call ufbseq(lnbufr, bufrf, ndat, 1, iret, subset)
                      call writsb(lnbufr)
                      n = n + 1
                   end if
                end do
             end do

             call closbf(lnbufr)
             close (lnbufr)
             close (12)
             print *, 'number of output bufr: ', n
             print *, ' '

100          continue

             return
          end

          subroutine getnoaafile(path_satdat, numfile, numcha, file, ier)
!***********************************************************************
! subroutine/function : getnoaafile
!
! usage :
!    call getnoaafile(path_satdat,numfile,numcha,file,ier)
!
! description      : to get noaa satellinte 1d format input file names.
!
! arguments :
!  i/o/w   name,      type,       description
!    i    path_satdat c*80        directory of noaa satellite input data.
!    o     numfile    integer     total number of input file names.
!    o     file(900)  c*180       noaa satellite 1-d file names.
!    o    numcha(900) int array
!    o     ier        integer     error message.
!                                 =0, success;  =1, failure.
!
! modules called :
!   getfilenames (getiofile.c)
!***********************************************************************

             implicit none
             character*80 path_satdat
             character*180 file(900)
             integer numfile, numcha(900), ier

             character*900000 cdat
             integer i, isize, ii, j, i1, i2, jj

             ier = 0
             numfile = 0
             isize = 900000
             do i = 1, isize
                cdat(i:i) = char(0)
             end do
             call getfilenames(path_satdat, isize, cdat, ier)
             if (ier .ne. 0) return

             ii = 1
             do i = isize, 1, -1
                if (cdat(i:i) .ne. char(0)) then
                   ii = i
                   go to 10
                end if
             end do
             return
10           continue
             if (ii .lt. 5) then
                ier = 1
                return
             end if

             j = 1
20           continue
             j = j + 1
             if (j .gt. ii) go to 30
             if (cdat(j:j) .eq. ' ') then
                i2 = j - 1
                jj = i2 - i1 + 1
                if (cdat(i2 - 3:i2) .eq. '.l1d') then
                   numfile = numfile + 1
                   if (numfile .gt. 900) then
                      print *, 'error in subroutine: getnoaafile'
                      print *, 'please modify parameter file(900) and numcha(900).'
                      print *, 'numfile: ', numfile
                      ier = 1
                      return
                   end if
                   file(numfile) (1:jj) = cdat(i1:i2)
                   numcha(numfile) = jj
                end if
                i1 = j + 1
             else
                if (j .eq. ii) then
                   i2 = ii
                   jj = i2 - i1 + 1
                   if (cdat(i2 - 3:i2) .eq. '.l1d') then
                      numfile = numfile + 1
                      if (numfile .gt. 900) then
                         print *, 'error in subroutine: getnoaafile'
                         print *, 'please modify parameter file(900) and numcha(900).'
                         print *, 'numfile: ', numfile
                         ier = 1
                         return
                      end if
                      file(numfile) (1:jj) = cdat(i1:i2)
                      numcha(numfile) = jj
                   end if
                   go to 30
                end if
             end if
             go to 20
30           continue

             return
          end

          subroutine dc_mhs(nbt, scan, machine, xtrack, iscan, isqc &
                            , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                            , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp, ier)
!***********************************************************************
! subroutine/function : dc_mhs
!
! usage
!    call dc_mhs(nbt,head,machine,xtrack,iscan,isqc
!   1       ,year,month,day,hour,minute,second,hmsl,ifqc,isty
!   2       ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!
! description      : to decode scan line data of amsub-1d or mhs data.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     nbt        integer     total number of input byte. (=12288)
!    i     scan(nbt)  byte        scan line data array.
!    i     machine    integer     machine index,
!                                 =0, for hp, sgi machine.
!                                 .not. 0, for dec, pc machine.
!    i     xtrack     integer     number of fov. (=90)
!    o     iscan      integer     scan line number.
!    o     isqc       integer     scan line quality flags.
!                                 =0, good
!                                 =1, bad
!    o     year       integer     scan line year.
!    o     month      integer     scan line month.
!    o     day        integer     scan line day.
!    o     hour       integer     scan line hour.
!    o     minute     integer     scan line minute.
!    o     second     real        scan line second.
!    o     hmsl       real        satellite altitude above refernce ellipsoid. (meters)
!    o   ifqc(xtrack) int array   fov quality flag.
!                                 =0, good
!                                 =1, bad
!    o   isty(xtrack) int array   surface type of fov.
!                                 =0, land
!                                 =1, sea
!    o   rlat(xtrack) real array  latitude in degrees of fov.
!    o   rlon(xtrack) real array  longitude in degrees of fov.
!    o     ht(xtrack) real array  surface height of fov in metres.
!    o   stzn(xtrack) real array  satellite zenith angle in degrees of fov.
!    o   sozn(xtrack) real array  solar zenith angle in degrees of fov.
!    o   staz(xtrack) real array  satellite azimuth angle in degrees of fov.
!    o   soaz(xtrack) real array  solar azimuth angle in degrees of fov.
!    o    brtmp(5,90) real array  scene brightness temperature in k ch.16 - ch.20
!    o     ier        integer     error message.
!                                 =0, success.
!                                 =1, failure.
!
! modules called : by2int4,  jumonday1
!
! note: the line scan data structure
!
!       scan line information:  4*9 bytes
!       navigation:             4*723 bytes
!       earth observations:     4*1800 bytes
!       pre-processing output:  4*540 bytes
!
!  data    field
!  type    description
!  ....    .................................................
!                     scan line information (4*9 bytes)
!  i  4    scan line number
!  i  4    scan line year
!  i  4    scan line day of year
!  i  4    scan line utc time of day in milliseconds
!  i  4    quality indicator bit field -- in all of the follwong,
!          if the bit is on (i.e., if it is set to 1) then the statement is true.
!          otherwise it is false.
!          bit 31:  do not use scan for product generation
!          bit 30:  time sequence error detected with this scan (see below)
!          bit 29:  data gap precedes this scan
!          bit 28:  no calibration (see below)
!          bit 27:  no earth location (see below)
!          bit 26:  first good time following a clock update
!          bit 25:  instrument status changed with this scan
!          bit 24-0: spare <zero fill>
!  i  4    scan line quality flags -- if bit is on (=1) then true
!          time problem code: (all bits off implies the scan time is as expected.)
!          bit 31-24: spare <zero fill>
!          bit 23:  time field is bad but can probaly be inferred from the previous
!                   good time.
!          bit 22:  time field is bad and cannot be inferred from the previous good time.
!          bit 21:  this record starts a sequence that is inconsistent with previous times
!                   (i.e., there is a time discontinuity).  this may or may not be
!                   associated with a spacecraft clock update.
!          bit 20:  start of a sequence that apparently repeates scan times that have
!                   been prevously accepted.
!          bit 19-16: spare <zero fill>
!          calibration problem code: (note these bits compliment the channel indicators;
!            all bits set to 0 indicates normal calibration.)
!          bit 15:  scan line was not calibrated because of bad time.
!          bit 14:  scan line was calibrated using fewer than the preferred number of scan
!                   lines because of proximity to start or end of data set or to a data gap.
!          bit 13:  scan line was not calibrated because of bad or insufficient prt data.
!          bit 12:  scan line was calibrated but with marginal prt data.
!          bit 11:  some uncalibrated channels on this scan. (see channel indicators.)
!          bit 10:  uncalibrated due to instrument mode.
!          bit 09:  questionable calibration because of antenna position error of space view.
!          bit 08:  questionable calibration because of antenna position error of blackbody.
!          earth location problem code: (all bits set to 0 implies the earth location
!            was normal)
!          bit 07:  not earth located because of bad time.
!          bit 06:  earth location questionable because of questionable time code.
!                   (see time problem flags above.)
!          bit 05:  earth location questionable -- only marginal agreement with
!                   reasonableness check.
!          bit 04:  earth location questionable -- fails reasonableness check.
!          bit 03:  earth location questionable because of antenna position check.
!          bit 02-0: spare <zero fill>
!  i  4    mixer chan 18-20 instrument temperature (k*100)
!  i  4*2  spare
!  ....    .................................................
!                       navigation (4*723 bytes)
!  i  4    10000*(latitude in degrees of position 1)
!  i  4    10000*(longitude in degrees of position 1)
!  i  4    10000*(latitude in degrees of position 2)
!  i  4    10000*(longitude in degrees of position 2)
!          ....
!  i  4    10000*(latitude in degrees of position 90)
!  i  4    10000*(longitude in degrees of position 90)
!  i  4    surface height of position 1 in metres
!  i  4    surface type of position 1 (0=sea, 1=mixed, 2=land)
!  i  4    surface height of position 2 in metres
!  i  4    surface type of position 2 (0=sea, 1=mixed, 2=land)
!          ....
!  i  4    surface height of position 90 in metres
!  i  4    surface type of position 90 (0=sea, 1=mixed, 2=land)
!  i  4    100*(local zenith angle in degrees of position 1)
!  i  4    100*(local azimuth angle in degrees of position 1)
!  i  4    100*(solar zenith angle in degrees of position 1)
!  i  4    100*(solar azimuth angle in degrees of position 1)
!  i  4    100*(local zenith angle in degrees of position 2)
!  i  4    100*(local azimuth angle in degrees of position 2)
!  i  4    100*(solar zenith angle in degrees of position 2)
!  i  4    100*(solar azimuth angle in degrees of position 2)
!          ....
!  i  4    100*(local zenith angle in degrees of position 90)
!  i  4    100*(local azimuth angle in degrees of position 90)
!  i  4    100*(solar zenith angle in degrees of position 90)
!  i  4    100*(solar azimuth angle in degrees of position 90)
!  i  4    10*(satellite altitude above refernce ellipsoid. km)
!  i  4*2  spare
!  ....    .................................................
!                    earth observations (4*1800 bytes)
!  i  4    100* scene brightness temp. in k. fov 1, amsu-a ch.1
!          (missing data indicator is -999999)
!  i  4    100* scene brightness temp. in k. fov 1, amsu-a ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 1, amsu-a ch.15
!  i  4    100* scene brightness temp. in k. fov 1, amsu-b ch.16
!  i  4    100* scene brightness temp. in k. fov 1, amsu-b ch.17
!          ....
!  i  4    100* scene brightness temp. in k. fov 1, amsu-b ch.20
!  i  4    100* scene brightness temp. in k. fov 2, amsu-a ch.1
!  i  4    100* scene brightness temp. in k. fov 2, amsu-a ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 2, amsu-a ch.15
!  i  4    100* scene brightness temp. in k. fov 2, amsu-b ch.16
!  i  4    100* scene brightness temp. in k. fov 2, amsu-b ch.17
!          ....
!  i  4    100* scene brightness temp. in k. fov 2, amsu-b ch.20
!          ........
!  i  4    100* scene brightness temp. in k. fov 90, amsu-a ch.1
!  i  4    100* scene brightness temp. in k. fov 90, amsu-a ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 90, amsu-a ch.15
!  i  4    100* scene brightness temp. in k. fov 90, amsu-b ch.16
!  i  4    100* scene brightness temp. in k. fov 90, amsu-b ch.17
!          ....
!  i  4    100* scene brightness temp. in k. fov 90, amsu-b ch.20
!  ....    .................................................
!                   pre-processing output (4*540 bytes)
!  i  4    pre-processing quality flags for fov 1: (all bits off implies acceptable data)
!          bit 31:  spare <zero fill>
!          bit 30:  set if amsu-b used secondary calibration from nearest neighbour amsu-a:
!          bit 29:  set if amsu-a used secondary calibration
!          bit 28:  set if amsu-a data missing
!          bit 27:  maximum probability scheme cloud flag
!          bit 26:  scattering test (only set over the sea)
!          bit 25:  logistic precipitation probability test
!          bit 24:  grody light rainfall test
!          bit 23:  mismatch between amsu-a/b 89ghz values
!          bit 22:  mismatch between surface type from topography dataset and from
!                   pre-processing
!          bit 21:  spare
!          bit 20:  set if 89ghz channel greatly different from that in surrounding fovs
!          bit 19:  scattering test (only set over the sea) - using amsu-b 89ghz channel
!          bit 18:  mismatch between amsu-a/b 89ghz values
!          bit 17-1: spare <zero fill>
!          bit 00:  set if amsu-b data missing
!  i  4    nearest-neightbour amsu-a estimated surface type for fov 1:
!          1 = bare young ice (i.e. new ice, no snow)
!          2 = dry land (i.e. dry. with or without vegetation)
!          3 = dry snow (i.e. snow with water less than 2%, over land)
!          4 = multi-year ice (i.e. old ice with dry snow cover)
!          5 = sea (i.e. open water, no islands, ice-free, ws < 14m/s)
!          6 = wet forest (i.e. established forest with wet canopy)
!          7 = wet land (i.e. non-forestedland with a wet suface)
!          8 = wet snow (i.e. water content > 2%, over land or ice)
!          9 = desert
!  i  4    cost function from ppasurf surface identification for fov 1
!  i  4    scattering index (recalculated with amsu-b 89ghz) for fov 1
!  i  4*2  spare
!  i  4    pre-processing quality flags for fov 2
!  i  4    nearest-neightbour amsu-a estimated surface type for fov 2
!  i  4    cost function from ppasurf surface identification for fov 2
!  i  4    scattering index (recalculated with amsu-b 89ghz) for fov 2
!  i  4*2  spare
!          ....
!  i  4    pre-processing quality flags for fov 90
!  i  4    nearest-neightbour amsu-a estimated surface type for fov 90
!  i  4    cost function from ppasurf surface identification for fov 90
!  i  4    scattering index (recalculated with amsu-b 89ghz) for fov 90
!  i  4*2  spare
!
!***********************************************************************
             implicit none
             integer nbt, xtrack
             byte scan(nbt), by4(4)
             integer machine
             integer iscan, isqc, year, month, day, hour, minute, ier
             real second, hmsl
             integer ifqc(xtrack), isty(xtrack)
             real rlat(xtrack), rlon(xtrack), ht(xtrack)
             real stzn(xtrack), sozn(xtrack), staz(xtrack), soaz(xtrack)
             real brtmp(5, 90)
             integer i, j, k, kk, ij, jday, mn, ihour
             integer in4

             ier = 0
             if (nbt .ne. 12288) then
                print *, 'error in subroutine: dc_mhs.'
                print *, 'input total byte number is error.'
                print *, 'nbt(=12288): ', nbt
                ier = 1
                return
             end if
             if (xtrack .ne. 90) then
                print *, 'error in subroutine: dc_mhs.'
                print *, 'input total xtrack number is error.'
                print *, 'xtrack(=90): ', xtrack
                ier = 1
                return
             end if

!-----------------------------------------------------------------------
!c  to decode scan line information.

             do i = 1, 4
                by4(i) = scan(i)
             end do
             call by2int4(machine, by4, in4)
             iscan = in4

             do i = 1, 4
                j = i + 4
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             year = in4

             do i = 1, 4
                j = i + 8
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             jday = in4
             call jumonday1(year, jday, month, day)

             do i = 1, 4
                j = i + 12
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             mn = in4/60000
             second = 0.001*(in4 - 60000*mn)
             hour = mn/60
             minute = mn - 60*hour

             do i = 1, 4
                j = i + 16
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             i = in4/2
             isqc = in4 - 2*i

!-----------------------------------------------------------------------
!c  to decode navigation.

             do k = 1, 90
                do i = 1, 4
                   j = 36 + 8*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                rlat(k) = 0.0001*in4
                do i = 1, 4
                   j = 36 + 8*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                rlon(k) = 0.0001*in4
             end do
             do k = 1, 90
                do i = 1, 4
                   j = 756 + 8*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                ht(k) = in4
                do i = 1, 4
                   j = 756 + 8*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                isty(k) = -9
                if (in4 .eq. 2) isty(k) = 0
                if (in4 .eq. 0) isty(k) = 1
                if (in4 .eq. 1) isty(k) = 2
                if ((isty(k) .ne. 0) .and. (isty(k) .ne. 1)) then
                   isty(k) = 1
                   if (ht(k) .gt. 0) isty(k) = 0
                end if
             end do
             do k = 1, 90
                do i = 1, 4
                   j = 1476 + 16*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                stzn(k) = 0.01*in4
                do i = 1, 4
                   j = 1476 + 16*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                staz(k) = 0.01*in4
                do i = 1, 4
                   j = 1476 + 16*(k - 1) + i + 8
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                sozn(k) = 0.01*in4
                do i = 1, 4
                   j = 1476 + 16*(k - 1) + i + 12
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                soaz(k) = 0.01*in4
             end do
             do i = 1, 4
                j = 2916 + i
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             hmsl = 100.*in4

!-----------------------------------------------------------------------
!c  to decode earth observations.

             do k = 1, 90
                ij = 2928 + 80*(k - 1)
                ifqc(k) = 0
                do kk = 1, 5
                   do i = 1, 4
                      j = ij + 60 + 4*(kk - 1) + i
                      by4(i) = scan(j)
                   end do
                   call by2int4(machine, by4, in4)
                   brtmp(kk, k) = 0.01*in4
                   if (in4 .lt. 0) ifqc(k) = 1
                end do
             end do
             if (isqc .eq. 0) then
                ij = 0
                do k = 1, 90
                   ij = ij + ifqc(k)
                end do
                if (ij .eq. 90) isqc = 1
             end if

             return
          end

          subroutine dc_amsua(nbt, scan, machine, xtrack, iscan, isqc &
                              , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                              , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp, ier)
!***********************************************************************
! subroutine/function : dc_amsua
!
! usage
!    call dc_amsua(nbt,scan,machine,xtrack,iscan,isqc
!   1       ,year,month,day,hour,minute,second,hmsl,ifqc,isty
!   2       ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!
! description      : to decode scan line data of amsua-1d data.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     nbt        integer     total number of input byte. (=4608)
!    i     scan(nbt)  byte        scan line data array.
!    i     machine    integer     machine index,
!                                 =0, for hp, sgi machine.
!                                 .not. 0, for dec, pc machine.
!    i     xtrack     integer     number of fov. (=30)
!    o     iscan      integer     scan line number.
!    o     isqc       integer     scan line quality flags.
!                                 =0, good
!                                 =1, bad
!    o     year       integer     scan line year.
!    o     month      integer     scan line month.
!    o     day        integer     scan line day.
!    o     hour       integer     scan line hour.
!    o     minute     integer     scan line minute.
!    o     second     real        scan line second.
!    o     hmsl       real        satellite altitude above refernce ellipsoid. (meters)
!    o   ifqc(xtrack) int array   fov quality flag.
!                                 =0, good
!                                 =1, bad
!    o   isty(xtrack) int array   surface type of fov.
!                                 =0, land
!                                 =1, sea
!                                 =2, mixed
!    o   rlat(xtrack) real array  latitude in degrees of fov.
!    o   rlon(xtrack) real array  longitude in degrees of fov.
!    o     ht(xtrack) real array  surface height of fov in metres.
!    o   stzn(xtrack) real array  satellite zenith angle in degrees of fov.
!    o   sozn(xtrack) real array  solar zenith angle in degrees of fov.
!    o   staz(xtrack) real array  satellite azimuth angle in degrees of fov.
!    o   soaz(xtrack) real array  solar azimuth angle in degrees of fov.
!    o   brtmp(15,30) real array  scene brightness temperature in k ch.1 - ch.15
!    o     ier        integer     error message.
!                                 =0, success.
!                                 =1, failure.
!
! modules called : by2int4,  jumonday1
!
! note: the scan line data structure
!
!       scan line information:  4*11 bytes
!       navigation:             4*243 bytes
!       earth observations:     4*600 bytes
!       pre-processing output:  4*210 bytes
!       spare bytes:            4*88 bytes
!
!  data    field
!  type    description
!  ....    .................................................
!                      scan line information (4*11 bytes)
!  i  4    scan line number
!  i  4    scan line year
!  i  4    scan line day of year
!  i  4    scan line utc time of day in milliseconds
!  i  4    quality indicator bit field -- in all of the follwong,
!          if the bit is on (i.e., if it is set to 1) then the statement is true.
!          otherwise it is false.
!          bit 31:  do not use scan for product generation
!          bit 30:  time sequence error detected with this scan (see below)
!          bit 29:  data gap precedes this scan
!          bit 28:  no calibration (see below)
!          bit 27:  no earth location (see below)
!          bit 26:  first good time following a clock update
!          bit 25:  instrument status changed with this scan
!          bit 24-0: spare <zero fill>
!  i  4    scan line quality flags -- if bit is on (=1) then true
!          time problem code: (all bits off implies the scan time is as expected.)
!          bit 31-24: spare <zero fill>
!          bit 23:  time field is bad but can probaly be inferred from the previous
!                   good time.
!          bit 22:  time field is bad and cannot be inferred from the previous good time.
!          bit 21:  this record starts a sequence that is inconsistent with previous times
!                   (i.e., there is a time discontinuity).  this may or may not be
!                   associated with a spacecraft clock update.
!          bit 20:  start of a sequence that apparently repeates scan times that have
!                   been prevously accepted.
!          bit 19-16: spare <zero fill>
!          calibration problem code: (note these bits compliment the channel indicators;
!            all bits set to 0 indicates normal calibration.)
!          bit 15:  scan line was not calibrated because of bad time.
!          bit 14:  scan line was calibrated using fewer than the preferred number of scan
!                   lines because of proximity to start or end of data set or to a data gap.
!          bit 13:  scan line was not calibrated because of bad or insufficient prt data.
!          bit 12:  scan line was calibrated but with marginal prt data.
!          bit 11:  some uncalibrated channels on this scan. (see channel indicators.)
!          bit 10:  uncalibrated due to instrument mode.
!          bit 09:  questionable calibration because of antenna position error of space view.
!          bit 08:  questionable calibration because of antenna position error of blackbody.
!          earth location problem code: (all bits set to 0 implies the earth location
!            was normal)
!          bit 07:  not earth located because of bad time.
!          bit 06:  earth location questionable because of questionable time code.
!                   (see time problem flags above.)
!          bit 05:  earth location questionable -- only marginal agreement with
!                   reasonableness check.
!          bit 04:  earth location questionable -- fails reasonableness check.
!          bit 03:  earth location questionable because of antenna position check.
!          bit 02-0: spare <zero fill>
!  i  4    amsu-a1 instrument rf shelf temperature (k*100)
!  i  4    amsu-a2 instrument rf shelf temperature (k*100)
!  i  4    amsu-b instrument mixer scan 18-20 temperature (k*100)
!  i  4*2  spare
!  ....    .................................................
!                       navigation (4*243 bytes)
!  i  4    10000*(latitude in degrees of position 1)
!  i  4    10000*(longitude in degrees of position 1)
!  i  4    10000*(latitude in degrees of position 2)
!  i  4    10000*(longitude in degrees of position 2)
!          ....
!  i  4    10000*(latitude in degrees of position 30)
!  i  4    10000*(longitude in degrees of position 30)
!  i  4    surface height of position 1 in metres
!  i  4    surface type of position 1 (0=sea, 1=mixed, 2=land)
!  i  4    surface height of position 2 in metres
!  i  4    surface type of position 2 (0=sea, 1=mixed, 2=land)
!          ....
!  i  4    surface height of position 30 in metres
!  i  4    surface type of position 30 (0=sea, 1=mixed, 2=land)
!  i  4    100*(local zenith angle in degrees of position 1)
!  i  4    100*(local azimuth angle in degrees of position 1)
!  i  4    100*(solar zenith angle in degrees of position 1)
!  i  4    100*(solar azimuth angle in degrees of position 1)
!  i  4    100*(local zenith angle in degrees of position 2)
!  i  4    100*(local azimuth angle in degrees of position 2)
!  i  4    100*(solar zenith angle in degrees of position 2)
!  i  4    100*(solar azimuth angle in degrees of position 2)
!          ....
!  i  4    100*(local zenith angle in degrees of position 30)
!  i  4    100*(local azimuth angle in degrees of position 30)
!  i  4    100*(solar zenith angle in degrees of position 30)
!  i  4    100*(solar azimuth angle in degrees of position 30)
!  i  4    10*(satellite altitude above refernce ellipsoid. km)
!  i  4*2  spare
!  ....    .................................................
!                   earth observations (4*600 bytes)
!  i  4    100* scene brightness temp. in k. fov 1, ch.1
!          (missing data indicator is -999999)
!  i  4    100* scene brightness temp. in k. fov 1, ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 1, ch.15
!  i  4    100* scene brightness temp. in k. fov 1. amsu-b ch.16
!  i  4    100* scene brightness temp. in k. fov 1. amsu-b ch.17
!          ....
!  i  4    100* scene brightness temp. in k. fov 1. amsu-b ch.20
!  i  4    100* scene brightness temp. in k. fov 2, ch.1
!  i  4    100* scene brightness temp. in k. fov 2, ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 2. amsu-b ch.20
!          ........
!  i  4    100* scene brightness temp. in k. fov 30, ch.1
!  i  4    100* scene brightness temp. in k. fov 30, ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 30. amsu-b ch.20
!  ....    .................................................
!                   pre-processing output (4*210 bytes)
!  i  4    pre-processing quality flags for fov 1: (all bits off implies acceptable data)
!          bit 31:  amsu bts considered contaminated. due e.g. to precip or surface type
!          bit 30:  set if amsu-a used secondary calibration
!          bit 29:  set if amsu-b used secondary calibration
!          bit 28:  set if amsu-b missing due to insufficient good data
!          bit 27:  maximum probability scheme cloud flag
!          bit 26:  scattering test (only set over the sea)
!          bit 25:  logistic precipitation probability test
!          bit 24:  grady light rainfall test
!          bit 23:  mismatch between amsu-a/b 89ghz values
!          bit 22:  mismatch between surface type from topography data and from pre-processing
!          bit 21-1: space <zero fill>
!          bit 00:  set if amsu-a data missing
!  i  4    estimated amsu-a surface type, for fov 1:
!          1 = bare young ice (i.e. new ice, no snow)
!          2 = dry land (i.e. dry. with or without vegetation)
!          3 = dry snow (i.e. snow with water less than 2%, over land)
!          4 = multi-year ice (i.e. old ice with dry snow cover)
!          5 = sea (i.e. open water, no islands, ice-free, ws < 14m/s)
!          6 = wet forest (i.e. established forest with wet canopy)
!          7 = wet land (i.e. non-forestedland with a wet suface)
!          8 = wet snow (i.e. water content > 2%, over land or ice)
!          9 = desert
!  i  4    cost function from surface identification (ppasurf), for fov 1
!  i  4    scattering index (ppascat), for fov1
!  i  4    logistic precipitation probability (ppcrosby), for fov 1
!  i  4*2  spare
!  i  4    pre-processing quality flags for fov 2
!  i  4    estimated amsu-a surface type, for fov 2
!  i  4    cost function from surface identification for fov 2
!  i  4    scattering index for fov 2
!  i  4    logistic precipitation probability for fov 2
!  i  4*2  spare
!          ....
!  i  4    pre-processing quality flags for fov 30
!  i  4    estimated amsu-a surface type, for fov 30
!  i  4    cost function from surface identification for fov 30
!  i  4    scattering index for fov 30
!  i  4    logistic precipitation probability for fov 30
!  i  4*2  spare
!  ....    .................................................
!                           spare (4*88 bytes)
!
!***********************************************************************
             implicit none
             integer nbt, xtrack
             byte scan(nbt), by4(4)
             integer machine
             integer iscan, isqc, year, month, day, hour, minute, ier
             real second, hmsl
             integer ifqc(xtrack), isty(xtrack)
             real rlat(xtrack), rlon(xtrack), ht(xtrack)
             real stzn(xtrack), sozn(xtrack), staz(xtrack), soaz(xtrack)
             real brtmp(15, 30)
             integer i, j, k, kk, ij, jday, mn, ihour
             integer in4

             ier = 0
             if (nbt .ne. 4608) then
                print *, 'error in subroutine: dc_amsua.'
                print *, 'input total byte number is error.'
                print *, 'nbt(=4608): ', nbt
                ier = 1
                return
             end if
             if (xtrack .ne. 30) then
                print *, 'error in subroutine: dc_amsua.'
                print *, 'input total xtrack number is error.'
                print *, 'xtrack(=30): ', xtrack
                ier = 1
                return
             end if

!-----------------------------------------------------------------------
!c  to decode scan line information.

             do i = 1, 4
                by4(i) = scan(i)
             end do
             call by2int4(machine, by4, in4)
             iscan = in4

             do i = 1, 4
                j = i + 4
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             year = in4

             do i = 1, 4
                j = i + 8
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             jday = in4
             call jumonday1(year, jday, month, day)

             do i = 1, 4
                j = i + 12
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             mn = in4/60000
             second = 0.001*(in4 - 60000*mn)
             hour = mn/60
             minute = mn - 60*hour

             do i = 1, 4
                j = i + 16
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             i = in4/2
             isqc = in4 - 2*i

!-----------------------------------------------------------------------
!c  to decode navigation.

             do k = 1, 30
                do i = 1, 4
                   j = 44 + 8*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                rlat(k) = 0.0001*in4
                do i = 1, 4
                   j = 44 + 8*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                rlon(k) = 0.0001*in4
             end do
             do k = 1, 30
                do i = 1, 4
                   j = 284 + 8*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                ht(k) = in4
                do i = 1, 4
                   j = 284 + 8*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                isty(k) = -9
                if (in4 .eq. 2) isty(k) = 0
                if (in4 .eq. 0) isty(k) = 1
                if (in4 .eq. 1) isty(k) = 2
                if ((isty(k) .ne. 0) .and. (isty(k) .ne. 1)) then
                   isty(k) = 1
                   if (ht(k) .gt. 0) isty(k) = 0
                end if
             end do
             do k = 1, 30
                do i = 1, 4
                   j = 524 + 16*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                stzn(k) = 0.01*in4
                do i = 1, 4
                   j = 524 + 16*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                staz(k) = 0.01*in4
                do i = 1, 4
                   j = 524 + 16*(k - 1) + i + 8
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                sozn(k) = 0.01*in4
                do i = 1, 4
                   j = 524 + 16*(k - 1) + i + 12
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                soaz(k) = 0.01*in4
             end do
             do i = 1, 4
                j = 1004 + i
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             hmsl = 100.*in4

!-----------------------------------------------------------------------
!c  to decode earth observations.

             do k = 1, 30
                ij = 1016 + 80*(k - 1)
                ifqc(k) = 0
                do kk = 1, 15
                   do i = 1, 4
                      j = ij + 4*(kk - 1) + i
                      by4(i) = scan(j)
                   end do
                   call by2int4(machine, by4, in4)
                   brtmp(kk, k) = 0.01*in4
                   if (in4 .lt. 0) ifqc(k) = 1
                end do
             end do
             if (isqc .eq. 0) then
                ij = 0
                do k = 1, 30
                   ij = ij + ifqc(k)
                end do
                if (ij .eq. 30) isqc = 1
             end if

             return
          end

          subroutine dc_hirs(nbt, scan, machine, xtrack, iscan, isqc &
                             , year, month, day, hour, minute, second, hmsl, ifqc, isty &
                             , rlat, rlon, ht, stzn, sozn, staz, soaz, brtmp, ier)
!***********************************************************************
! subroutine/function : dc_hirs
!
! usage
!    call dc_hirs(nbt,scan,machine,xtrack,iscan,isqc
!   1       ,year,month,day,hour,minute,second,hmsl,ifqc,isty
!   2       ,rlat,rlon,ht,stzn,sozn,staz,soaz,brtmp,ier)
!
! description      : to decode scan line data of hirs-1d data.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     nbt        integer     total number of input byte. (=15872)
!    i     scan(nbt)  byte        scan line data array.
!    i     machine    integer     machine index,
!                                 =0, for hp, sgi machine.
!                                 .not. 0, for dec, pc machine.
!    i     xtrack     integer     number of fov. (=56)
!    o     iscan      integer     scan line number.
!    o     isqc       integer     scan line quality flags.
!                                 =0, good
!                                 =1, bad
!    o     year       integer     scan line year.
!    o     month      integer     scan line month.
!    o     day        integer     scan line day.
!    o     hour       integer     scan line hour.
!    o     minute     integer     scan line minute.
!    o     second     real        scan line second.
!    o     hmsl       real        satellite altitude above refernce ellipsoid. (meters)
!    o   ifqc(xtrack) int array   fov quality flag.
!                                 =0, good
!                                 =1, bad
!    o   isty(xtrack) int array   surface type of fov.
!                                 =0, land
!                                 =1, sea
!                                 =2, mixed
!    o   rlat(xtrack) real array  latitude in degrees of fov.
!    o   rlon(xtrack) real array  longitude in degrees of fov.
!    o     ht(xtrack) real array  surface height of fov in metres.
!    o   stzn(xtrack) real array  satellite zenith angle in degrees of fov.
!    o   sozn(xtrack) real array  solar zenith angle in degrees of fov.
!    o   staz(xtrack) real array  satellite azimuth angle in degrees of fov.
!    o   soaz(xtrack) real array  solar azimuth angle in degrees of fov.
!    o   brtmp(20,56) real array  scene brightness temperature in k for ch.1 - ch.19
!                                 scene radiance in wm-2sr-1(cm-1)-1 for ch.20
!    o     ier        integer     error message.
!                                 =0, success.
!                                 =1, failure.
!
! modules called : by2int4,  jumonday1
!
! note: the line scan data structure
!
!       scan line information:  4*12 bytes
!       navigation:             4*451 bytes
!       earth observations:     4*2240 bytes
!       avhrr:                  4*728 bytes
!       pre-processing output:  4*448 bytes
!       spare:                  4*89 bytes
!
!  data    field
!  type    description
!  ....    .................................................
!                    scan line information (4*12 bytes)
!  i  4    scan line number
!  i  4    scan line year
!  i  4    scan line day of year
!  i  4    scan line utc time of day in milliseconds
!  i  4    quality indicator bit field -- in all of the follwong,
!          if the bit is on (i.e., if it is set to 1) then the statement is true.
!          otherwise it is false.
!          bit 31:  do not use scan for product generation
!          bit 30:  time sequence error detected with this scan (see below)
!          bit 29:  data gap precedes this scan
!          bit 28:  no calibration (see below)
!          bit 27:  no earth location (see below)
!          bit 26:  first good time following a clock update
!          bit 25:  instrument status changed with this scan
!          bit 24-0: spare <zero fill>
!  i  4    scan line quality flags -- if bit is on (=1) then true
!          time problem code: (all bits off implies the scan time is as expected.)
!          bit 31-24: spare <zero fill>
!          bit 23:  time field is bad but can probaly be inferred from the previous
!                   good time.
!          bit 22:  time field is bad and cannot be inferred from the previous good time.
!          bit 21:  this record starts a sequence that is inconsistent with previous times
!                   (i.e., there is a time discontinuity).  this may or may not be
!                   associated with a spacecraft clock update.
!          bit 20:  start of a sequence that apparently repeates scan times that have
!                   been prevously accepted.
!          bit 19-16: spare <zero fill>
!          calibration problem code: (note these bits compliment the channel indicators;
!            all bits set to 0 indicates normal calibration.)
!          bit 15:  scan line was not calibrated because of bad time.
!          bit 14:  scan line was calibrated using fewer than the preferred number of scan
!                   lines because of proximity to start or end of data set or to a data gap.
!          bit 13:  scan line was not calibrated because of bad or insufficient prt data.
!          bit 12:  scan line was calibrated but with marginal prt data.
!          bit 11:  some uncalibrated channels on this scan. (see channel indicators.)
!          bit 10:  uncalibrated due to instrument mode.
!          bit 09 and 08:  spare <zero fill>
!          earth location problem code: (all bits set to 0 implies the earth location
!            was normal)
!          bit 07:  not earth located because of bad time.
!          bit 06:  earth location questionable because of questionable time code.
!                   (see time problem flags above.)
!          bit 05:  earth location questionable -- only marginal agreement with
!                   reasonableness check.
!          bit 04:  earth location questionable -- fails reasonableness check.
!          bit 03-0:  spare <zero fill>
!  i  4    hirs instrument baseplate temperature (k*100)
!  i  4    amsu-a1 rf shelf instrument temperature (k*100)
!  i  4    amsu-a2 rf shelf instrument temperature (k*100)
!  i  4    amsu-b mixer chan 18-20 instrument temperature (k*100)
!  i  4*2  spare
!  ....    .................................................
!                     navigation (4*451 bytes)
!  i  4    10000*(latitude in degrees of position 1)
!  i  4    10000*(longitude in degrees of position 1)
!  i  4    10000*(latitude in degrees of position 2)
!  i  4    10000*(longitude in degrees of position 2)
!          ....
!  i  4    10000*(latitude in degrees of position 56)
!  i  4    10000*(longitude in degrees of position 56)
!  i  4    surface height of position 1 in metres
!  i  4    surface type of position 1 (0=sea, 1=mixed, 2=land)
!  i  4    surface height of position 2 in metres
!  i  4    surface type of position 2 (0=sea, 1=mixed, 2=land)
!          ....
!  i  4    surface height of position 56 in metres
!  i  4    surface type of position 56 (0=sea, 1=mixed, 2=land)
!  i  4    100*(local zenith angle in degrees of position 1)
!  i  4    100*(local azimuth angle in degrees of position 1)
!  i  4    100*(solar zenith angle in degrees of position 1)
!  i  4    100*(solar azimuth angle in degrees of position 1)
!  i  4    100*(local zenith angle in degrees of position 2)
!  i  4    100*(local azimuth angle in degrees of position 2)
!  i  4    100*(solar zenith angle in degrees of position 2)
!  i  4    100*(solar azimuth angle in degrees of position 2)
!          ....
!  i  4    100*(local zenith angle in degrees of position 56)
!  i  4    100*(local azimuth angle in degrees of position 56)
!  i  4    100*(solar zenith angle in degrees of position 56)
!  i  4    100*(solar azimuth angle in degrees of position 56)
!  i  4    10*(satellite altitude above refernce ellipsoid. km)
!  i  4*2  spare
!  ....    .................................................
!                    earth observations (4*2240 bytes)
!  i  4    100* scene brightness temp. in k. fov 1, hirs ch.1
!          (missing data indicator is -999999)
!  i  4    100* scene brightness temp. in k. fov 1, hirs ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 1, hirs ch.19
!  i  4    1000* scene radiance in w*m**2sr**-1*cm**-1 for fov1, ch.20
!  i  4    100* scene brightness temp. in k. fov 1, amsu-a ch.1
!  i  4    100* scene brightness temp. in k. fov 1, amsu-a ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 1, amsu-a ch.15
!  i  4    100* scene brightness temp. in k. fov 1, amsu-b ch.16
!  i  4    100* scene brightness temp. in k. fov 1, amsu-b ch.17
!          ....
!  i  4    100* scene brightness temp. in k. fov 1, amsu-b ch.20
!  i  4    100* scene brightness temp. in k. fov 2, hirs ch.1
!  i  4    100* scene brightness temp. in k. fov 2, hirs ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 2, hirs ch.19
!  i  4    1000* scene radiance in w*m**2sr**-1*cm**-1 for fov2, ch.20
!  i  4    100* scene brightness temp. in k. fov 2, amsu-a ch.1
!  i  4    100* scene brightness temp. in k. fov 2, amsu-a ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 2, amsu-a ch.15
!  i  4    100* scene brightness temp. in k. fov 2, amsu-b ch.16
!  i  4    100* scene brightness temp. in k. fov 2, amsu-b ch.17
!          ....
!  i  4    100* scene brightness temp. in k. fov 2, amsu-b ch.20
!          ........
!  i  4    100* scene brightness temp. in k. fov 56, hirs ch.1
!  i  4    100* scene brightness temp. in k. fov 56, hirs ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 56, hirs ch.19
!  i  4    1000* scene radiance in w*m**2sr**-1*cm**-1 for fov56, ch.20
!  i  4    100* scene brightness temp. in k. fov 56, amsu-a ch.1
!  i  4    100* scene brightness temp. in k. fov 56, amsu-a ch.2
!          ....
!  i  4    100* scene brightness temp. in k. fov 56, amsu-a ch.15
!  i  4    100* scene brightness temp. in k. fov 56, amsu-b ch.16
!  i  4    100* scene brightness temp. in k. fov 56, amsu-b ch.17
!          ....
!  i  4    100* scene brightness temp. in k. fov 56, amsu-b ch.20
!  ....    .................................................
!                      avhrr (4*728 bytes)
!  i  4*12*56  13 words for each hirs fov, 56 fovs
!  ....    .................................................
!                  pre-processing output (4*448 bytes)
!  i  4    fov quality flags (hirs) for fov 1: (all bits off implies acceptable data)
!          bit 31:  spare <zero fill>
!          bit 30:  set if secondary calibration used
!          bit 29-22:  spare <zero fill>
!          bit 21:  hirs cloud test (tbd)
!          bit 20-1:  bit n set to 1 if brightness temperature in hirs channel n is
!                     missing or unreasonable
!          bit 0:  bad or missing data (in any or all channels)
!  i  4    fov quality flags (hirs) for fov 2
!          ....
!  i  4    fov quality flags (hirs) for fov 56
!  i  4    pre-processing quality flags for fov 1: (all bits off implies acceptable data)
!          bit 31:  set if amsu-a surface types not all the same
!          bit 30:  set if amsu-a used secondary calibration
!          bit 29:  set if amsu-b used secondary calibration
!          bit 28:  set if amsu-b missing
!          bit 27:  flag for cloud cost set for any amsu-a
!          bit 26:  scattering tflag set for any amsu-a (only set over the set)
!          bit 25:  logistic precipitation probability test
!                   calculated from amsu-a data mapped to hirs grid
!          bit 24:  grody light rainfall test calculated on hirs grid
!          bit 23:  mismatch between amsu-a/b 89ghz values for any amsu-a
!          bit 22:  mismatch between surface type from topography dataset and from
!                   pre-processing (any amsu-a)
!          bit 21-4:  spare <zero fill>
!          bit 03:  set when avhrr chan 3 is albedo, not bright. temp.
!          bit 02:  cloud cost flag (recalculated on hirs grid)
!          bit 01:  scattering flag (recalculated on hirs grid)
!          bit 00:  set if amsu-a and amsu-b data missing
!  i  4    estimated nearest amsu surface type for fov 1:
!          1 = bare young ice (i.e. new ice, no snow)
!          2 = dry land (i.e. dry, with or without vegetation)
!          3 = dry snow (i.e. snow with water less than 2%, over land)
!          4 = multi-year ice (i.e. old ice with dry snow cover)
!          5 = sea (i.e. open water, no islands, ice-free, ws < 14 m/s)
!          6 = wet forest (i.e. established forest with wt canopy)
!          7 = wet land (i.e. non-forested land with a wet surface)
!          8 = wet snow (i.e. water content > 2%, over land or ice)
!          9 = desert
!  i  4    cost function from ppasurf surface identification for fov 1
!  i  4    scattering index for fov 1
!  i  4    logistic precipitation probability for fov 1
!  i  4*2  spare
!  i  4    pre-processing quality flags for fov 2
!  i  4    estimated nearest amsu surface type for fov 2
!  i  4    cost function from ppasurf surface identification for fov 2
!  i  4    scattering index for fov 2
!  i  4    logistic precipitation probability for fov 2
!  i  4*2  spare
!          ....
!  i  4    pre-processing quality flags for fov 56
!  i  4    estimated nearest amsu surface type for fov 56
!  i  4    cost function from ppasurf surface identification for fov 56
!  i  4    scattering index for fov 56
!  i  4    logistic precipitation probability for fov 56
!  i  4*2  spare
!  ....    .................................................
!                      spare (4*89 bytes)
!
!***********************************************************************

             implicit none
             integer nbt, xtrack
             byte scan(nbt), by4(4)
             integer machine
             integer iscan, isqc, year, month, day, hour, minute, ier
             real second, hmsl
             integer ifqc(xtrack), isty(xtrack)
             real rlat(xtrack), rlon(xtrack), ht(xtrack)
             real stzn(xtrack), sozn(xtrack), staz(xtrack), soaz(xtrack)
             real brtmp(20, 56), radiate
             integer i, j, k, kk, ij, jday, mn, ihour
             integer in4

             ier = 0
             if (nbt .ne. 15872) then
                print *, 'error in subroutine: dc_hirs.'
                print *, 'input total byte number is error.'
                print *, 'nbt(=15872): ', nbt
                ier = 1
                return
             end if
             if (xtrack .ne. 56) then
                print *, 'error in subroutine: dc_hirs.'
                print *, 'input total xtrack number is error.'
                print *, 'xtrack(=56): ', xtrack
                ier = 1
                return
             end if

!-----------------------------------------------------------------------
!c  to decode scan line information.

             do i = 1, 4
                by4(i) = scan(i)
             end do
             call by2int4(machine, by4, in4)
             iscan = in4

             do i = 1, 4
                j = i + 4
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             year = in4

             do i = 1, 4
                j = i + 8
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             jday = in4
             call jumonday1(year, jday, month, day)

             do i = 1, 4
                j = i + 12
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             mn = in4/60000
             second = 0.001*(in4 - 60000*mn)
             hour = mn/60
             minute = mn - 60*hour

             do i = 1, 4
                j = i + 16
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             i = in4/2
             isqc = in4 - 2*i

!-----------------------------------------------------------------------
!c  to decode navigation.

             do k = 1, 56
                do i = 1, 4
                   j = 48 + 8*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                rlat(k) = 0.0001*in4
                do i = 1, 4
                   j = 48 + 8*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                rlon(k) = 0.0001*in4
             end do
             do k = 1, 56
                do i = 1, 4
                   j = 496 + 8*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                ht(k) = in4
                do i = 1, 4
                   j = 496 + 8*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                isty(k) = -9
                if (in4 .eq. 2) isty(k) = 0
                if (in4 .eq. 0) isty(k) = 1
                if (in4 .eq. 1) isty(k) = 2
                if ((isty(k) .ne. 0) .and. (isty(k) .ne. 1)) then
                   isty(k) = 1
                   if (ht(k) .gt. 0) isty(k) = 0
                end if
             end do
             do k = 1, 56
                do i = 1, 4
                   j = 944 + 16*(k - 1) + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                stzn(k) = 0.01*in4
                do i = 1, 4
                   j = 944 + 16*(k - 1) + i + 4
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                staz(k) = 0.01*in4
                do i = 1, 4
                   j = 944 + 16*(k - 1) + i + 8
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                sozn(k) = 0.01*in4
                do i = 1, 4
                   j = 944 + 16*(k - 1) + i + 12
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                soaz(k) = 0.01*in4
             end do

             do i = 1, 4
                j = 1840 + i
                by4(i) = scan(j)
             end do
             call by2int4(machine, by4, in4)
             hmsl = 100.*in4

!-----------------------------------------------------------------------
!c  to decode earth observation.

             do k = 1, 56
                ij = 1852 + 160*(k - 1)
                ifqc(k) = 0
                do kk = 1, 19
                   do i = 1, 4
                      j = ij + 4*(kk - 1) + i
                      by4(i) = scan(j)
                   end do
                   call by2int4(machine, by4, in4)
                   brtmp(kk, k) = 0.01*in4
                   if (in4 .lt. 0) ifqc(k) = 1
                end do
                do i = 1, 4
                   j = ij + 76 + i
                   by4(i) = scan(j)
                end do
                call by2int4(machine, by4, in4)
                brtmp(20, k) = 10.e10
!         brtmp(20,k)=0.0001*in4
!         radiate=0.0001*in4
!         brtmp(20,k)=20862.6/log(1.+36310465./radiate)
             end do
             if (isqc .eq. 0) then
                ij = 0
                do k = 1, 56
                   ij = ij + ifqc(k)
                end do
                if (ij .eq. 56) isqc = 1
             end if

             return
          end

          subroutine dc_gen_inf(genbyte, machine, sateid, instrument, nline &
                                , year, month, day, hour, minute, ier)
!***********************************************************************
! subroutine/function : dc_gen_inf
!
! usage
!    call dc_gen_inf(genbyte,machine,sateid,instrument,nline
!   1               ,year,month,day,hour,minute,ier)
!
! description      : to decode general information  of noaa-nn 1d data.
!
! arguments :
!  i/o/w   name,      type,       description
!    i    genbyte(96) byte        general information byte data array.
!    i     machine    integer     machine index,
!                                 =0, for hp, sgi machine.
!                                 .not. 0, for dec, pc machine.
!    o     sateid     integer     satellite id.
!    o     instrument integer     instrument code.
!                                 =5, hirs;  =6, msu;  =10, amsu-a;
!                                 =11, amsu-b;  =12, mhs.
!    o     nline      integer     count of scan lines in this data set.
!    o     year       integer     start of data set year.
!    o     month      integer     start of data set month.
!    o     day        integer     start of data set day.
!    o     hour       integer     start of data set hour.
!    o     minute     integer     start of data set minute.
!    o     ier        integer     error message.
!                                 =0, success.
!                                 =1, failure.
!
! modules called : by2int4,  jumonday1
!
! general information structure: 4*24 bytes
!  data    field
!  type    description
!  ....    .................................................
!  c  3    1d dataset creation id (nss for nesdis, etc)
!  c  1    filler
!  c  3    site of originating centre for 1b data
!          (msc for cwb satellite center)
!  c  1    filler
!  i  4    level 1d format version number
!  i  4    level 1d format version year
!  i  4    level 1d format version day of year
!  i  4    count of header records in this data set (=1)
!  i  4    satellite id (wmo code)
!  i  4    instrument code (5=hirs; 6=msu; 10=amsu-a; 11=amsu-b; 12=mhs)
!  i  4    10 * (nominal satellite altitude, km) (=8790)
!  i  4    nominal orbit period (seconds)
!  i  4    orbit number (at start of dataset)
!  i  4    start of data set year
!  i  4    start of data set day of year
!  i  4    start of data set utc time of day in milliseconds
!  i  4    orbit number (at end of dataset)
!  i  4    end of data set year
!  i  4    end of data set day of year
!  i  4    end of data set utc time of day in milliseconds
!  i  4    count of scan lines in this data set
!  i  4    count of missing scan lines
!  i  4    atovpp version number (values above 9000 indicate test vm)
!  i  4    instruments present (bit0=hirs, bit1=msu, bit2=amusu-a,
!          bit3=amusu-b, bit4=avhrr)
!  i  4    version number of data set for antenna corrections
!          (=0 if data not corrected)
!  i  4    spare
!***********************************************************************
             implicit none
             byte genbyte(96), by4(4)
             integer machine, sateid, instrument, nline, ier
             integer year, month, day, hour, minute
             integer in4
             integer i, j, jday, mn

             ier = 0

!-----------------------------------------------------------------------
!c  to decode satellite id

             do i = 1, 4
                j = i + 24
                by4(i) = genbyte(j)
             end do
             call by2int4(machine, by4, in4)
             sateid = in4

!-----------------------------------------------------------------------
!c  to decode general information.

             do i = 1, 4
                j = i + 28
                by4(i) = genbyte(j)
             end do
             call by2int4(machine, by4, in4)
             if ((in4 .ne. 5) .and. (in4 .ne. 10) .and. (in4 .ne. 11) .and. &
                 (in4 .ne. 12)) then
                ier = 1
                print *, 'error in subroutine: dc_amsua.'
                print *, 'this data header is not amsu-a data header.'
                print *, 'instrument code 5=hirs, 10=amsu-a, 11=amsu-b, 12=mhs'
                print *, 'instrument code: ', in4
                return
             end if
             instrument = in4

             do i = 1, 4
                j = i + 44
                by4(i) = genbyte(j)
             end do
             call by2int4(machine, by4, in4)
             year = in4
             do i = 1, 4
                j = i + 48
                by4(i) = genbyte(j)
             end do
             call by2int4(machine, by4, in4)
             jday = in4
             call jumonday1(year, jday, month, day)
             do i = 1, 4
                j = i + 52
                by4(i) = genbyte(j)
             end do
             call by2int4(machine, by4, in4)
             mn = in4/60000
             hour = mn/60
             minute = mn - 60*hour

             do i = 1, 4
                j = i + 72
                by4(i) = genbyte(j)
             end do
             call by2int4(machine, by4, in4)
             nline = in4

             return
          end

          subroutine jumonday1(iyear, jday, month, iday)
!***********************************************************************
! subroutine/function : jumonday1
!
! usage :
!    call jumonday1(iyear,jday,month,iday)
!
! description      : to get month and day giving the day of year.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     iyear      integer     year.
!    i     jday       integer     the day of year.
!    o     month      integer     month.
!    o     iday       integer     day.
!
! modules called : none
!***********************************************************************

             implicit none
             integer iyear, jday, month, iday
             integer mon(12), mon1(12)
             data mon/31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/
             data mon1/31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366/
             integer id, i1, i2, i3, i4, i

             id = 0
             i1 = iyear/4
             i4 = iyear - 4*i1
             if (i4 .eq. 0) id = 1
             i1 = iyear/400
             i2 = iyear - 400*i1
             if (i2 .eq. 0) id = 1
             i1 = iyear/100
             i3 = iyear - 100*i1
             if ((i3 .eq. 0) .and. (i2 .ne. 0)) id = 0
             if (id .eq. 0) then
                do i = 1, 12
                   if (jday .le. mon(i)) then
                      month = i
                      if (month .eq. 1) then
                         iday = jday
                      else
                         iday = jday - mon(i - 1)
                      end if
                      return
                   end if
                end do
             else
                do i = 1, 12
                   if (jday .le. mon1(i)) then
                      month = i
                      if (month .eq. 1) then
                         iday = jday
                      else
                         iday = jday - mon1(i - 1)
                      end if
                      return
                   end if
                end do
             end if
             return
          end

          subroutine by2int4(machine, by4, in4)
!***********************************************************************
! subroutine/function : by2int4
!
! usage :
!    call by2int4(by4,in4)
!
! description      : to transform 4-bytes to integer for hp, sgi machine.
!
! arguments :
!  i/o/w   name,      type,       description
!    i     machine    integer     machine index,
!                                 =0, for hp, sgi machine.
!                                 .not. 0, for dec, pc machine.
!    i     by4(4)     byte array  input 4 bytes.
!    o     in4        integer   output integer.
!
! modules called : none
!***********************************************************************

             integer in4, ib
             byte by4(4), by_tem(4)
             equivalence(ib, by_tem(1))
             if (machine .eq. 0) then
                by_tem(1) = by4(4)
                by_tem(2) = by4(3)
                by_tem(3) = by4(2)
                by_tem(4) = by4(1)
             else
                by_tem(1) = by4(1)
                by_tem(2) = by4(2)
                by_tem(3) = by4(3)
                by_tem(4) = by4(4)
             end if
             in4 = ib
             return
          end
