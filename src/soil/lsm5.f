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
      program lsm5

      implicit none

      integer    istatus
      integer    nx_l,ny_l
       
      call get_grid_dim_xy(nx_l,ny_l,istatus)
      if(istatus.eq.1)then
        write(6,*)'laps parameters obtained'
      else
        write(6,*)'istatus = ',istatus,'error - get_laps_config'
        write(6,*)'terminating laps-lsm5. soil moisture'
        stop
      end if
      call lsm5_sub(nx_l,ny_l)
      stop
      end

c
c ************************************************************************** 
c
      subroutine lsm5_sub(imax,jmax)

c     program laps soilmoisture analysis
c     changed to enable gridded analysis and real time
c     chandran subramaniam
c     2/8/93
      
c     john smart 12/1/93: adapting to run in real time on the
c     			  unix platform.  set up laps standard i/o.
c                         enhanced modularity.
c     j smart    9/22/97  adapt for dynamic array memory allocation
c
      integer   imax,jmax

      include 'soilm.inc'
c
c**** model12 is a soil moisture content model developed in june l986, and
c     based upon the remote sensing based watershed model developed at
c     the university of maryland from 1982 to 1985.  this is version 12.4,
c     as submitted to water resources bulletin in january 1989.
c     created by groves.

      real 	ksat,lamda,in
      data 	day,sumr,in,oldwea/4*0./
      integer   istatus, i, j, ii
      integer   istatus_precip, istatus_m
      integer   istatus_w,istatus_n,istatus_e
      integer   soiltype(imax,jmax)
      real      laps_u(imax,jmax)     !u-component
      real      laps_v(imax,jmax)     !v-component
      real      laps_t(imax,jmax)     !temperature
      real      laps_td(imax,jmax)    !dewpoint temperature
      real      laps_in(imax,jmax)    !infiltration
      real      laps_wx(imax,jmax)    !weather data; grid pt wet or dry.
      real      laps_wfz(imax,jmax)   !wetting front depth, z
      real      laps_mwf(imax,jmax)   !wetting front moisture content
      real      laps_mwf_pre(imax,jmax)!previous wetting front moist content.
      real      laps_evap(imax,jmax)  !evaporation
      real      laps_rain(imax, jmax) !radar-estimated 1hr liq. precip
      real      laps_sc(imax,jmax)    !snow cover; fractional  0.0 to 1.0.
      real      laps_sc_pre(imax,jmax)!previous hour snow cover   " .
c     real      laps_sm(imax,jmax)    !snow melt, undefined currently
      real      laps_smc_3d(imax,jmax,3)!three layer soil moisture content
      real      soilm_field_cap(imax,jmax)!field cap soil moistr (m**3/m**3)
      real      soilm_sat(imax,jmax) !saturated soil moistr (m**3/m**3)
      real      data(imax,jmax,7)    !holds lm2 output data - current time.
      real      data_s(imax,jmax,4)  !holds lsx input data - current time

      logical   griddry,filefound
      integer   loop_bound
      integer   i4time_smcur, i4time_smpre(25),
     &          lvl_s(4),lvl_1(3),lvl_2(7),lvl_l
      character ftime_smcur*9
c laps precip
      character ext_l*31, dir_l*150, var_l*3, lvl_coord_l*4,
     &          units_l*10, comment_l*125
c laps surface
      character ext_s*31, dir_s*150, var_s(4)*3, lvl_coord_s(4)*4,
     &          units_s(4)*10, comment_s(4)*125
c background
      character var_bkg(4)*3
c laps lm1...3-layer % soil saturation
      character dir_1*150, ext1*31, var_1(3)*3, lvl_coord1(3)*4,
     &          units1(3)*10, comment1(3)*125
c laps lm2...variety of soil characteristic variables
      character dir_2*150, ext2*31, var_2(7)*3, lvl_coord2(7)*4, 
     &          units2(7)*10, comment2(7)*125
      character atime_smpre(25)*24
      character fname*200

      data	ext_s/'lsx'/
      data	ext_l/'l1s'/
      data	ext1/'lm1'/
      data	ext2/'lm2'/
      data      var_s/'u  ','v  ','t  ','td '/
      data      var_bkg/'usf','vsf','tsf','dsf'/
      data      var_l/'r01'/
      data      var_1/'lsm', 'lsm', 'lsm'/
      data      var_2/'civ', 'dwf', 'mwf', 'wx ', 'evp', 'sc ', 'sm'/
      data      units1/'%','%','%'/
      data      units2/'m**3','m','m**3/m**3',' ','m/s','%',' '/
      data      comment1/
     &          'layer 1 (0-6in [0-0.152m]) % soil saturation',
     &          'layer 2 (6-12in [.152-.305m]) % soil saturation',
     &          'layer 3 (12-36in [.305-0.914m]) % soil saturation'/
      data comment2(1)/'cumulative infiltration volume (m)'/
      data comment2(2)/'depth of wetting front (m)'/
      data comment2(3)/'moisture content of wetting front (m**3/m**3)'/
      data comment2(4)/'weather indicator - pos = wet, neg = dry'/
      data comment2(5)/'evaporation (m/s)'/
      data comment2(6)/'snow cover value (fractional, 0 to 1.0)'/
      data comment2(7)/'snow melt value, (gm/s)'/
      data	lvl_1/-1,-2,-3/
c 
c first read in current time
c
      call get_directory('etc',fname,len)
!     open(11, file=fname(1:len)//'systime.dat',status='unknown')
!     read(11,*,err=999)i4time_smcur
!     read(11,22,err=999)ftime_smcur
!22   format(1x,a9)
!     close(11)
      call get_systime(i4time_smcur,ftime_smcur,istatus)
c -------------------------
      print*,'systime: ',ftime_smcur,' i4time: ',i4time_smcur

c**** read soil description and simulation time step 

      call get_directory('lm1',dir_1,len)

      call soil_in5(imax,jmax,soiltype,istatus)   	! get soil texture group
      if(istatus.ne.1)then
         write(6,*)'soil textures obtained'
         write(6,*)'soil type = ',soiltype(1,1)
      else
         write(6,*)'soil textures not obtained - terminating'
         stop
      end if
c
c  get the current and previous i4time and ascii times. allow for 4 previous
c  times for the previous soil moisture product.  mandatory to have 1 previous
c  time so if the 4th previous time is the current time then the 5th is prev.

      loop_bound=5
      icnt=0
      do i=1,loop_bound
         i4time_smpre(i) = i4time_smcur - 3600*icnt
         call cv_i4tim_asc_lp(i4time_smpre(i),
     &                        atime_smpre(i),istatus)
         icnt = icnt + 1
      end do
      print*,'current time: ',atime_smpre(1)

c  get current surface data: laps lsx; u- v-component, temp and dew point.
       call get_directory('lsx',dir_s,len)
       do i=1,4
          lvl_s(i)=0
       end do

       filefound = .false.
       ii = 0
       do while (.not.filefound)
          ii = ii + 1
          call read_laps_data(i4time_smpre(ii), dir_s, ext_s,
     &              imax, jmax, 4, 4, var_s, lvl_s, lvl_coord_s, 
     &              units_s, comment_s, data_s, istatus)

         if(istatus.eq.1)then
            write(6,*)'laps surface data retrieved'
            write(6,*)'sfc directory [dir_s] = ',dir_s
            write(6,*)'time retrieved:',atime_smpre(ii)
            filefound = .true.

            laps_u=data_s(:,:,1)
            laps_v=data_s(:,:,2)
            laps_t=data_s(:,:,3)
            laps_td=data_s(:,:,4)

         elseif(ii .ge. loop_bound)then
            filefound = .true.
            write(6,*)'laps surface data not available'
            print*,'lets try model bkg for surface u/v/t/td'
            call get_modelfg_2d(i4time_smcur,var_bkg(1)
     1,imax,jmax,laps_u,istatus)
            if(istatus.ne.1)goto 59
            call get_modelfg_2d(i4time_smcur,var_bkg(2)
     1,imax,jmax,laps_v,istatus)
            if(istatus.ne.1)goto 59
            call get_modelfg_2d(i4time_smcur,var_bkg(3)
     1,imax,jmax,laps_t,istatus)
            if(istatus.ne.1)goto 59
            call get_modelfg_2d(i4time_smcur,var_bkg(4)
     1,imax,jmax,laps_td,istatus)

59          if(istatus.ne.1)then
               print*,'failed to get background in get_modelfg_2d'
               print*,'terminating laps soil moisture'
               return
            endif
         end if
       end do
c
c get laps snow cover
c -------------------
       write(6,*)'compute snow cover'
       write(6,*)

       call readcsc(i4time_smcur
     &              ,imax,jmax
     &              ,laps_sc)

       write(6,*)'snow cover computed'
       write(6,*)'***************************************'
       write(6,*)'getting field cap and saturated soil moisture'

       do j = 1, jmax
          do i = 1, imax
             isoil = soiltype(i,j)
c   get default soil hydraulic parameters and initial soil moisture
             call soils(isoil,ksat,ths,thr,psif,psiae,lamda)
             call amc(th_field_cap,isoil,ths)
             soilm_sat(i,j) = ths
             soilm_field_cap(i,j) = th_field_cap
          enddo
       enddo
c
c get previous soil moisture data: laps lm2; infiltration, depth of wetting
c front (wf), previous weather (ie, wet or dry), and moisture content of wf.
c
       write(6,*) 'obtaining soil moisture data'

       call get_directory('lm2',dir_2,len)
c       dir_2 = '../lapsprd/lm2/'
       do i=1,7
          lvl_2(i) = 0
       end do
c
       filefound = .false.
       ii=2
       do while(.not.filefound)
c
          call read_laps_data(i4time_smpre(ii), dir_2, ext2,
     &       imax, jmax, 7, 7, var_2, lvl_2, lvl_coord2,
     &       units2, comment2, data, istatus)
c
          if(istatus.eq.1)then
            write(6,*)'laps lm2 data retrieved'
            write(6,*)'lm2 directory [dir_2] = ',dir_2
            write(6,*)'time retrieved:',i4time_smpre(ii)
             do j = 1,jmax
             do i = 1,imax
                laps_in(i,j) = data(i,j,1)
                laps_wfz(i,j) = data(i,j,2)
                laps_mwf(i,j) = data(i,j,3)
                laps_wx(i,j) = data(i,j,4)
                laps_evap(i,j) = data(i,j,5)
                laps_sc_pre(i,j) = data(i,j,6)
                laps_mwf_pre(i,j) = data(i,j,7)
             end do
             end do
             filefound = .true.
c
c             write(6,*)'current soil moisture product retrieved'
c
          else
c
             if(ii.gt.14)then
c
c   use initial sm from field cap
c
             write(6,*) '***************************************'
             write(6,*) 'using default (field cap) soil moisture'
             do j = 1, jmax
             do i = 1, imax
                laps_mwf(i,j) = soilm_field_cap(i,j)
                laps_mwf_pre(i,j) = soilm_field_cap(i,j)
             enddo
             enddo
             filefound = .true.
c
             end if
c
             ii = ii + 1
c
          end if
       end do
c
c   get radar estimated precipitation.  laps l1s. use "data_s" to hold.
c
       write(6,*)'get laps precip and snow total'
       call get_directory('l1s',dir_l,len)
c       dir_l = '../lapsprd/l1s/'
       lvl_l = 0
       call read_laps_data(i4time_smcur, dir_l, ext_l,
     &           imax, jmax, 1, 1, var_l, lvl_l, lvl_coord_l, 
     &           units_l, comment_l, laps_rain, istatus_precip)
       if(istatus_precip.eq.1)then
c
          write(6,*)
          write(6,*)'retrieved laps rain data'
          write(6,*)'l1s directory [dir_l] = ',dir_l
          write(6,*)'time retrieved:',i4time_smcur
          do j = 1,jmax
             do i = 1,imax
                laps_rain(i,j)=laps_rain(i,j)*100.    !cm/hr
             end do
          end do
          griddry = .false.
c
       else
c
          write(6,*)'laps precip not available - assume dry'
          write(6,*)'l1s directory [dir_l] = ',dir_l
          write(6,*)'time attempted:',i4time_smcur
          griddry = .true.

       end if
c
c this concludes getting the initial data for the soil moisture model
c ****************************************************************************
c ********************** soil moisture subroutine ****************************
c
       call soil_moisture(imax,jmax,laps_u,laps_v,
     &      laps_t,laps_td,laps_rain,laps_sc,
     &      laps_in,laps_wfz,laps_mwf,laps_mwf_pre,
     &      laps_wx,soiltype,griddry,laps_evap,laps_smc_3d,
     &                  istatus)
c
c arrays laps_in, _wfz, _mwf, _wx return the new values of in, wfz, mwf and wx
c from the above subroutine, thus they are both input and output.  laps_smc_3d
c is returned with current three layer soil moisture content.
c
       if(istatus.eq.1)then
          write(6,*)' soil moisture computation complete'
          write(6,*)' check results for bad values'
c
          bad_lower_bndry = 0.0
          bad_upper_bndry = 100.0
          kmax=1
c
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_mwf,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_m)

          if(istatus_m .lt. 0) write(6,910) istatus_m
c
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_wfz,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_w)

          if(istatus_w .lt. 0) write(6,911) istatus_w
c
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_in,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_n)

          if(istatus_n .lt. 0) write(6,912) istatus_n
c
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_evap,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_e)

          if(istatus_e .lt. 0) write(6,913) istatus_e
c
          kmax=3
          call lsm_qc_check(imax,jmax,kmax,
     &                      laps_smc_3d,
     &                      bad_lower_bndry,
     &                      bad_upper_bndry,
     &                      istatus_3d)
c
          if(istatus_3d .lt. 0) write(6,914) istatus_3d
c
910       format(' +++ warning.  mwf:  istatus_m = ',i8,'+++')
911       format(' +++ warning.  wfz:  istatus_w = ',i8,'+++')
912       format(' +++ warning.  civ:  istatus_n = ',i8,'+++')
913       format(' +++ warning.  evap:  istatus_e = ',i8,'+++')
914       format(' +++ warning.  smc-3d:  istatus_3d = ',i8,'+++')
c
          write(6,*)' writing current hour soil moisture for grid'
          do j=1,jmax
          do i=1,imax
             data(i,j,1)=laps_in(i,j)
             data(i,j,2)=laps_wfz(i,j)
             data(i,j,3)=laps_mwf(i,j)
             data(i,j,4)=laps_wx(i,j)
             data(i,j,5)=laps_evap(i,j)
             data(i,j,6)=laps_sc(i,j)
             data(i,j,7)=laps_mwf_pre(i,j)
          end do
          end do
c
          call write_laps_data(i4time_smcur, dir_1, ext1,
     &             imax, jmax,3, 3, var_1, lvl_1, lvl_coord1,
     &             units1, comment1, laps_smc_3d, istatus)
          if(istatus.eq.1)then
             write(6,*)'3d soil moisture successfully written'
          else
             write(6,*)'error writing 3d soil moisture - laps_smc_3d'
          end if
c
          call write_laps_data(i4time_smcur, dir_2, ext2,
     &             imax, jmax,7, 7, var_2, lvl_2, lvl_coord2,
     &             units2, comment2, data, istatus)
c
          if(istatus.eq.1)then
             write(6,*)'2d soil moisture fields successfully written'
          else
             write(6,*)'error writing 2d soil moisture fields'
          end if
c
       else
c
          write(6,*)'error during soil moisture computation'
c
       end if
       goto 100
c
 999   write(6,*)'error reading systime.dat file'
       go to 100
 998   write(6,*)'error opening log file at:',ftime_smcur
c
100    write(6,*)'end of soil moisture simulation'
       close(6)
c
       stop
       end
