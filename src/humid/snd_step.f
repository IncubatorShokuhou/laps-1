cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
c
c
c
      subroutine snd_step (
     1     i4time,              ! current time i/
     1     p_3d,                !laps pressure grid i/
     1     radiometer_switch,   !hight to use data meters i/
     1     snd_lookback,        !paramter indicating time window i/
     1     lat,                 !laps lat grid i/
     1     lon,                 !laps lon grid i/
     1     laps_t,              !laps lt1 field i/
     1     ii,jj,kk,            !laps grid dimensions i/
     1     q_snd,               !data output (q or sh) g/kg /o
     1     weight_field,        !weight field for this data for variational /o
     1     raob_radius,         !radius of influence (meters)
     1     abort                !abort flag, 0= abort, 1= go i/o
     1     )



      implicit none


c     parameters for sounding acquisition

      integer ll, nnn
      parameter (ll = 300)      ! max levels in a sounding
      parameter (nnn = 10000)   ! max numbers of soundings handled

c     input parameters (defined above)

      integer i4time, ii,jj,kk, snd_lookback, abort
      integer radiometer_switch
      real q_snd(ii,jj,kk), p_3d (ii,jj,kk)
      real raob_radius
      real weight_field(ii,jj,kk)
      real lat(ii, jj), lon(ii, jj)
      real laps_t (ii,jj,kk)  ! laps temp field

c  dynamic dependent parameters  

      real td(ll,nnn)           !sounding dewpoint temp (c)
      integer ll_n (nnn)        !levels in given sounding n
      real lat_s(ll,nnn)        !sounding lat by level (account for drift)
      real lon_s(ll,nnn)        !sounding lon by level (account for drift)
      real t(ll,nnn)            !sounding temp by level (c)
      real rh(ll,nnn)           !sounding rh by level
      real meter_ht(ll,nnn)     !msl ht in meters (added for stick ware data)
      real p_s(ll,nnn)          !pressure at sounding level l (hpa)
      character*8 instrument(nnn) !created to track radiometer data 
      real rh_fields(ii,jj,kk)  !analyzed rh field
      
      



      integer mask(ii,jj,kk)    !mask for determining weights
      real climo_w(kk)          !climo for qc checking (unknown use now)
      
c     externals
      
      real make_rh              ! external function
      external make_rh
      real ssh2                 ! external function
      external ssh2
      real make_ssh             ! external function
      external make_ssh
      
      
c     normal internal parameters
      
      integer i,j,k,ks,is,l,n,nn
      integer look_back_time
      integer raob_i4time       ! actual time of raob data used
      integer raob_i4time_chk   ! used to validate actual raob
      integer numoffiles, istatus, maxfiles
      character*256 c_filenames(200),pathname_in
      character*9 filename
      character*200 fname
      integer idx, len
      real lat_r,lon_r
      integer  dummy, isound
      real rdummy
      character*5 cdummy
      integer idummy
      character*9 r_filename
      character*8 snd_type
      real temt, temtd
      real r, rmin
      real d2r,pi
      real tot_weight,tot_scale_weight,scale_used,weight_avg,counter
      real wt_ck(ii,jj)         ! weight check to test for ample raob data
      real rmd
      integer x_sum
      real ave(kk),adev(kk),sdev(kk),var(kk),skew(kk),curt(kk)
      integer satsnd_counter, raob_counter
      
c     climate model (for qc)
      
      real center_lat
      integer jday
      real 
     1     standard_press(40),
     1     tempertur_guess(40),
     1     mixratio_guess(40)


      write (6,*)
      write (6,*)' begin snd step for raob insertion'
      write (6,*)
      
      
     
c     *** begin routine
      
      if (abort .eq. 0) return  ! abort = 1 means run with data
      
      pi = acos(-1.0)
      d2r = pi/180.
      look_back_time = snd_lookback
      maxfiles = nnn
      
      do j = 1,jj
         do i = 1,ii
            do k = 1,kk
            weight_field(i,j,k) = 0.0
            rh_fields(i,j,k) = 0.5 !assign arbitrary 50%rh
            mask(i,j,k) = 0
            enddo
         enddo
      enddo
      
      call get_r_missing_data(rmd, istatus)   
      
c     *** read in the raob data
      
c     +++ get latest raob file name
      
      call get_directory('snd',pathname_in,len)
      pathname_in = pathname_in(1:len)//'*'
c     pathname_in = '../lapsprd/snd/*'
      
      call get_file_names (pathname_in,numoffiles,c_filenames,
     1     maxfiles,istatus)
      
      if (istatus.ne.1) then    ! no data to read
         return
         
      elseif(numoffiles.eq.0) then ! no data to read
         return
         
      endif

c     
c     new code to set up for reruns
c     
      do i =  numoffiles, 1, -1
         
         idx = index(c_filenames(i), ' ')
         filename = c_filenames(i)(idx-13:idx-4)
         call i4time_fname_lp (filename, raob_i4time, istatus)
         if (
     1        (i4time-look_back_time) .lt. raob_i4time
     1        .and.
     1        raob_i4time .le. i4time
     1        ) go to 27

      enddo
 27   continue

c     +++ validate the raob time
      
      if (raob_i4time + look_back_time .lt. i4time) then ! too old
         write (6,*) 'raob data found is too old to be used'
         abort = 0              !do not include in variational system
         return
      endif
      
      if (raob_i4time .gt. i4time) then ! too new
         write (6,*) 'warning raob data found is too new to be used'
         abort = 0              !do not include in variational system
         return
      endif

      write (6,*) 'succeeded in finding raob file, ',filename


c     +++ read snd file  ++++++++++++++  selectsnd modify n section++++++++
      call get_directory('snd',fname,len)

      open (12, file = fname(1:len)//filename//'.snd',
     1     form='formatted',status='old',err=18)

      n  = 0  ! initialize n 

      do isound = 1, nnn         ! isound mere counter here, index is n

c     read the header record for main data

 15      read (12,511,end=16,err=18) idummy,idx, 
     1        lat_r,lon_r,
     1        rdummy, cdummy,r_filename,snd_type
 511     format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)
         write(6,*) idx, 'advertized number of levels in sounding'


c     prepare to read individual levels in snd
         if (n.eq.nnn .or. n.lt.0) then
            write (6,*) 'n is going to be larger than nnn'
            write (6,*) 'or n is going to be zero index'
            write (6,*) n, ' is the value of n'
            abort = 0
            return
         endif


         n = n + 1               !count as sounding valid header and file

         ll_n(n) = 0            !initialize incremental level counter

         do i = 1, idx          !iterate through each level (idx)

            ll_n(n) = ll_n(n) + 1 !increment level counter

            read(12,*) meter_ht(ll_n(n),n),
     1           p_s (ll_n(n) , n),
     1           t (ll_n(n), n),
     1           td(ll_n(n), n),rdummy,rdummy

            if (
     1           ((p_s (ll_n(n),n).ne.rmd).or.(meter_ht(ll_n(n),n)
     1           .ne.rmd))
     1           .and.
     1           (t (ll_n(n),n).ne.rmd)
     1           .and.
     1           (td (ll_n(n),n).ne.rmd)
     1           ) then

c     now fill arrays lat_s and lon_s to enable future raob drift
c     placeholders
c     increment counter for next record


               lat_s(ll_n(n),n) = lat_r
               lon_s(ll_n(n),n) = lon_r
               instrument (n) = snd_type

            else

               ll_n(n) = ll_n(n) - 1 !decrement level counter invalid level
               
            endif


         enddo !  i !finished with all levels





c     levels now known
c     reject on time condition (one hour lookback)
         if(r_filename .ne. filename) then ! examine closer
            call i4time_fname_lp (r_filename, raob_i4time_chk, 
     1           istatus)
            if (abs(raob_i4time-raob_i4time_chk).gt.snd_lookback)then ! over limit
               write(6,*) 'rejecting on time bounds', r_filename
               n = n - 1 !reject -- out of time bounds
            else
               write(6,*) 'accepting.. ', r_filename,' ',
     1              snd_type
            endif

         else                   ! accept implicitly
            write(6,*) 'accepting.. time exact ', snd_type
         endif

         if (ll_n(n) .eq. 0) then
            n =n-1              ! no data in vertical test for blank radiometer data
         endif

      enddo  ! dummy isound loop (n determined)

 16   close (12)                ! end of file read, close and move on

c     invoke vaiable nn which is total number of soundings read 
      nn = n

      if (nn .le. 0 ) then ! no data found
         write(6,*) 'no usable raob data in database'
         abort = 0              ! 0 = skip in variational step
         return
      else
         write(6,*) nn, ' number of raobs considered in analysis'
      endif

c     -------end select snd modify n section  (see snd interface document 1.0)
c
c------------------------- modify radiometer soundings ------------------------

      do n = 1, nn
         if (instrument(n) .eq. 'radiomtr' ) then
            write (6,*) 'n is radiometer', n
            write (6,*) 'calling hypsometric equation for pressure'
            call hypso (meter_ht(1,n),t(1,n),rmd,ll_n(n),p_s(1,n),abort)

c     now perform second correction for radiometrics data.... data
c     is only good up until 2km or 2000 meters.  starting with the surface
c     loop throught the data and re-assign the "top" layer ll_n(n) of
c     each radiometer sounding to be under 2km.

            do i = 1, ll_n(n)
               if ((meter_ht(i,n)-meter_ht(1,n)).gt.radiometer_switch) 
     1              then        ! above limit of use
                  ll_n(n) = i ! this should maybe terminate loop
                  go to 212 ! preserve assignment
               endif
            enddo               ! hight limit correction
 212        continue

         endif

      enddo                     !radiometer correction for p and ht



c     ++++++++++++++++++++placeholder for qc step+++++++++++++++++++++++++++++
c     if there is eventually a qc step as there was for the old raob step, it would
c     go here.  it should be noted that the qc was deemed most useful for satellite
c     sounding data.  



c     ----------------------end qc section -------------------------------

c     +++++++++++++++++++++ start comput rh step +++++++++++++++++++++

      do n = 1, nn              ! for all soundings
         do l = 1, ll_n(n)      ! for all layers per given sounding

c     convert the dewpoint temp to rh for interpolation later in routine
c     note that rh is assumed to be wrt liquid at all temps per the standard
c     raob convention. temps are in c and pres is in hpa. ssh is g/kg as is the
c     appropriate input to make_rh.

c     rh is a fraction not percent at this point.

            rh(l,n) = make_rh(p_s(l,n),t(l,n),ssh2(p_s(l,n),t(l,n),
     1           td(l,n),-132.),-132.)

         enddo                  ! l levels done
      enddo                     ! n soundings done


c     an rh value is now available for all soundings and all levels

c     ------------------------end compute rh step -------------------------



c     +++++++++++++++++++++++start analz_snd+++++++++++++++++++++++++++++++++++

c     this routine is an analog of analz_gvap 
c     and analz_gps that produces analyzed fields of the 
c     rh with capability to follow raob or dropsonde drift.
c     this is by far the guts and most complex routine called.
c     the basic function of this routine are:
c     1) interpolate rh in 3-space to the laps grid
c     2) generate a snoothed field representing that space
c     3) generate a corresponding weight field for the variational step

      call analz_snd (
     1     rh,
     1     ll_n,
     1     lat_s,
     1     lon_s,
     1     p_s,
     1     ll,
     1     nnn,
     1     nn,
     1     raob_radius,
     1     lat,
     1     lon,
     1     p_3d,
     1     rh_fields,
     1     weight_field,
     1     ii,jj,kk)


c     ----------------------end analz_snd----------------------------------

c     ++++++++++++++++++++++convert to sh++++++++++++++++++++++++++++++++++
c
c     this module converts the analyzed rh field that is now in the laps system
c     to q or sh that is ready for 
c     utilization in the variational analysis system


      do k = 1,kk
         do j = 1,jj
            do i = 1,ii

               if (weight_field(i,j,k) .ne. 0.0
     1              .and. laps_t(i,j,k).ne. rmd) then      ! convert to q

                  q_snd(i,j,k) = make_ssh(p_3d(i,j,k), 
     1                 laps_t(i,j,k) -273.15, rh_fields(i,j,k), -132.)

               else             !pass rmd flag to data
                  
                  q_snd(i,j,k) = rmd

               endif

            enddo               ! end i
         enddo                  ! end j
      enddo                     ! end k


c     data(i,j,k) is now in g/kg sh for introduction to variational scheme

      abort = 1                 ! indicates ready for variational scheme
      write (6,*) 'success finishing snd_step.f'
      return
 18   write (6,*) 'error reading or operning file, abort'
      write (6,*) 'snd step routine has failed'
      abort = 0
      return
      end

