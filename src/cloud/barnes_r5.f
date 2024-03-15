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

      subroutine barnes_r5(t,imax,jmax,kmax,to,wt_p,cf_modelfg
     1  ,l_perimeter,cld_snd_in,wt_snd_in,r_missing_data
     1  ,grid_spacing_m,i_snd,j_snd,n_cld_snd,max_cld_snd
     1  ,max_obs,weight_modelfg,nx_dim_lut,ny_dim_lut
     1  ,ix_low,ix_high,iy_low,iy_high,n_fnorm,istatus)

!     1997 aug 01  k. dritz  - added nx_dim_lut, ny_dim_lut, ix_low,
!                              ix_high, iy_low, iy_high, and n_fnorm as
!                              dummy arguments.
!     1997 aug 01  k. dritz  - removed parameter statements for the above.
!     1997 aug 01  k. dritz  - changed nx_l_max to imax and ny_l_max to jmax.
!     1997 aug 01  k. dritz  - added r_missing_data as dummy argument.
!     1997 aug 01  k. dritz  - removed include of lapsparms.for.

      include 'laps_cloud.inc'

      integer nx_dim_lut,ny_dim_lut,ix_low,ix_high,iy_low,iy_high

      integer nskip,n_fnorm
!     parameter (nskip = 2)

      integer lowi_lut(imax)
      integer lowj_lut(jmax)

      dimension to(imax,jmax,kmax),t(imax,jmax,kmax),fnorm(n_fnorm)
     1  ,iob(max_obs),job(max_obs),kob(max_obs),nob(max_obs)
     1                          ,cf_modelfg(imax,jmax,kmax)
      dimension wt_p(imax,jmax,kmax)

      logical l_perimeter, l_use_snd, l_bterm
      parameter (l_bterm = .true.) ! limit search radius for soundings

      real cld_snd_in(max_cld_snd,kmax)
      real wt_snd_in(max_cld_snd,kmax)

      real*8, allocatable, dimension(:,:) :: cld_snd
      real*8, allocatable, dimension(:,:) :: wt_snd
      real*8, allocatable, dimension(:,:) :: cld_snd_diff
      real*8, allocatable, dimension(:,:) :: wt_snd_diff
!     real*8 cld_snd(max_cld_snd,kmax)
!     real*8 wt_snd(max_cld_snd,kmax)
!     real*8 cld_snd_diff(max_cld_snd,kmax)
!     real*8 wt_snd_diff(max_cld_snd,kmax)

      integer i_snd(max_cld_snd)
      integer j_snd(max_cld_snd)

      real*8 sum_a(imax,jmax)
      real*8 sumwt_a(imax,jmax)

      real*8 weight,sum,sumwt,fraci,fracj,z1,z2,z3,z4

      real*8 iiilut(-nx_dim_lut:nx_dim_lut,-ny_dim_lut:ny_dim_lut)
      integer nlast(kcloud)
      logical l_analyze(kcloud), l_diff_snd, l_debug

      write(6,*)' subroutine barnes_r5...'

      allocate( cld_snd(max_cld_snd,kmax), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate cldsnd'
      endif

      allocate( wt_snd(max_cld_snd,kcloud), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate wt_snd'
      endif

      allocate( cld_snd_diff(max_cld_snd,kmax), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate cldsnd_diff'
      endif

      allocate( wt_snd_diff(max_cld_snd,kcloud), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate wt_snd_diff'
      endif

!     convert from real to real*8
      cld_snd = cld_snd_in 
      wt_snd = wt_snd_in

      l_diff_snd = .false.  ! allow more efficiency

      l_debug = .false.     ! debug more efficient calculations

!     n_cld_snd = 9 ! enable for debug testing (use 'n_cld_snd_in' in arg list)

      if(n_cld_snd .gt. max_cld_snd)then
          write(6,*)' barnes_r5: error, too many cloud soundings'
          istatus = 0
          return
      endif 

!     note that this routine is considerably longer than it might need to
!     be to allow maximum efficiency. 'nskip' controls how many points
!     to skip in the i & j directions during the looping/analysis. the
!     remaining points are filled in later by a faster bilinear interpolation
!     step. also, l_analyze(k) controls which levels are analyzed. if the
!     obs at a given level are a repeat of the obs at the next lower level,
!     the looping again is skipped and the relevant summations are copied
!     up from the lower level.

!     obtain and/or iterate for value of nskip
      nskip = nint(20000. / grid_spacing_m)
      nskip = min(max(nskip,1),4)

100   rden_ratioi = ((float(imax)-1.)/float(nskip))
      rden_ratioj = ((float(jmax)-1.)/float(nskip))

      if(abs(rden_ratioi - nint(rden_ratioi))  .gt. .001
     1  .or.(rden_ratioj - nint(rden_ratioj))  .gt. .001)then
          write(6,*)' bad value of nskip'
          if(nskip .gt. 1)then
              nskip = nskip - 1
              goto 100
          elseif(nskip .eq. 1)then
              write(6,*)' code error - stop'
              stop
          endif
      else
          write(6,*)' good value of nskip = ',nskip
      endif

      if(l_diff_snd .and. l_debug)then
          nskip = 1
          write(6,*)' resetting nskip to ',nskip
      endif

      if(l_debug)then
          n_debug = 1000

!         extra output for test sounding (aia in roc domain)
          it1 = 169
          jt1 = 129

          it2 = it1-5
          jt2 = jt1

      else
          n_debug = 10

      endif

      write(6,*)' number of cloud soundings = ',n_cld_snd

      istatus = 1
      exponent_distance_wt = 5.0

      fnrmscl = 1.533 ! simplify later by setting to one?
      iiizero = 0
      iiione = 0
      iiifg = 0
      iiib = 0
      grid_spacing_mc = grid_spacing_m * 10000. / 10080.8 ! for testing

      dist_norm = 100000. ! at 100km distance, ob weight = 1.
      dist_norm = dist_norm * sqrt(fnrmscl/1.533)
      dist_norm_grid = dist_norm / grid_spacing_mc ! normalized dist in grid pts
      dist_norm_grid_sq = dist_norm_grid**2

      iiicorners = ((imax-1)**2 * (jmax-1)**2) * fnrmscl + 1

!     more efficient way to say weight = (d / 100km) ** exponent_distance_wt
      do iii = 1,n_fnorm ! iii is loosely the dist in grid points squared
        fnorm(iii) = (dist_norm_grid_sq / float(iii)) 
     1                                 ** (exponent_distance_wt / 2.0)       

        rnorm = sqrt(float(iii) / dist_norm_grid_sq) ! normalized radius
        range = (dist_norm * rnorm) / sqrt(fnrmscl)

!       add ramp so we reach zero earlier
        if(l_bterm)then
            b = max( min((6.-rnorm),1.) ,0.)
            fnorm_orig = fnorm(iii)
            fnorm(iii) = fnorm(iii) * b
            if(iii .eq. iiicorners)then
                write(6,*)'iiicorners/rnorm/b/fnorms/range ',  
     1          iiicorners,rnorm,b,fnorm_orig,fnorm(iii),range
            endif

            if(fnorm(iii) .eq. 0. .and. iiizero .eq. 0)then
                iiizero = iii
                write(6,*)
     1          ' fnorm array reached zero, iii/rnorm/b/range/orig=',
     1            iiizero,rnorm,b,range,fnorm_orig
            endif

            if(b .lt. 1. .and. iiib .eq. 0)then
                iiib = iii
                write(6,*)
     1          ' b slipped below 1.0, iii/rnorm/b/range/orig=',
     1            iii,rnorm,b,range,fnorm_orig
            endif

        endif

        if(fnorm(iii) .lt. 1. .and. iiione .eq. 0)then
            iiione = iii
            write(6,*)
     1      ' fnorm array reached one, iii/rnorm/range=',
     1        iii,rnorm,range
        endif

        if(fnorm(iii) .lt. weight_modelfg .and. iiifg .eq. 0)then
            iiifg = iii
            write(6,*)
     1      ' fnorm array reached wtmodelfg, iii/rnorm/range/fnorm=',
     1        iii,rnorm,range,fnorm(iii)
        endif

      enddo

      epsilon = 1e-5
      weight_epsilon = 1e-5 * (weight_modelfg / float(n_cld_snd))

      radius_of_influence_km =  ( weight_modelfg ** 
     1                           (-1./exponent_distance_wt) )
     1                                     * (dist_norm/1000.)
     1                                     * (1./sqrt(fnrmscl))

      radius_of_influence_eps =  ( weight_epsilon ** 
     1                           (-1./exponent_distance_wt) )
     1                                     * (dist_norm/1000.)
     1                                     * (1./sqrt(fnrmscl))

! lsw comment added for awips logging 7/1/04
      write(6,*)' optimizing model background weight ...'
      write(6,*)' model background weight  = ',weight_modelfg
      write(6,*)' dist_norm = ',dist_norm
      write(6,*)' approx radius of influence (km) = '
     1         ,radius_of_influence_km
      write(6,*)' min possible ob weight   = ',fnorm(n_fnorm)
      write(6,*)' max possible ob weight   = ',fnorm(1)
      write(6,*)' weight epsilon   = ',weight_epsilon
      write(6,*)' radius of influence (epsilon) = '
     1                               ,radius_of_influence_eps
      write(6,*)' l_bterm = ',l_bterm

      ncnt=0
      do k=1,kmax

        if(.true.)then

            nlast(k) = ncnt
            do n = 1,n_cld_snd
              if(cld_snd(n,k) .eq. r_missing_data) go to 233

!             test if out of bounds of established perimeter around laps domain
              if(i_snd(n) .lt. ix_low .or. 
     1           i_snd(n) .gt. ix_high) go to 233
              if(j_snd(n) .lt. iy_low .or. 
     1           j_snd(n) .gt. iy_high) go to 233

              ncnt=ncnt+1
              iob(ncnt)=i_snd(n)
              job(ncnt)=j_snd(n)
              kob(ncnt)=k
              nob(ncnt)=n
              nlast(k) = ncnt
233           continue
            enddo ! n

        endif 

        if(k .gt. 1)then
            km1 = k - 1
            l_analyze(k) = .false.

            if(.true.)then

              do n=1,n_cld_snd
                if(cld_snd(n,k) .ne. cld_snd(n,km1))then
                    l_analyze(k) = .true.
                    goto250
                endif
              enddo ! n

            endif 


250         continue

        else ! k .eq. 1
            l_analyze(k) = .true.
        endif

      enddo ! k

!     calculate difference soundings (soundings relative to level below)
!     this will allow more efficient summing
!     it is assumed for now that 'wt_snd' is always 1.00 or 'r_missing_data'
      if(l_diff_snd)then

          write(6,*)' calculating difference soundings'

!         initialize arrays
          cld_snd_diff = r_missing_data
          wt_snd_diff = r_missing_data

          do n = 1,n_cld_snd

              k = 1
              cld_snd_diff(n,1) = cld_snd(n,1)
              wt_snd_diff(n,1) = wt_snd(n,1)
              if(n .le. n_debug)then 
                  write(6,*)
                  write(6,*)' cloud sounding # ',n,' at ',i_snd(n)
     1                                                   ,j_snd(n)

                  if(l_debug)then
                      if(i_snd(n) .eq. it1 .and. j_snd(n) .eq. jt1)then
                          write(6,*)' ******test sounding ******'
                      endif
                  endif

                  write(6,901,err=902)k
     1                       ,cld_snd(n,k),cld_snd_diff(n,k)
     1                       ,wt_snd(n,k),wt_snd_diff(n,k)
 902              continue
              endif
              
              do k = 2,kmax              
                  if(cld_snd(n,k) .ne. cld_snd(n,k-1))then ! cloud 
                      if(cld_snd(n,k) .eq. r_missing_data)then
                          wt_snd_diff(n,k) = -wt_snd(n,k-1)
                          cld_snd_diff(n,k) = cld_snd(n,k-1)
                      elseif(cld_snd(n,k-1) .eq. r_missing_data)then
                          wt_snd_diff(n,k) = wt_snd(n,k)
                          cld_snd_diff(n,k) = cld_snd(n,k)
                      else ! diff the two soundings (both levels are present)
                          cld_snd_diff(n,k) = cld_snd(n,k)
     1                                      - cld_snd(n,k-1)      
!                         if(cld_snd_diff(n,k) .eq. 0.)then ! identical
!                             wt_snd_diff(n,k) = r_missing_data
!                         else                              ! different
                              wt_snd_diff(n,k) = 0d0
!                         endif
                      endif
                  endif

                  if(n .le. n_debug)then
                      write(6,901,err=904)k
     1                           ,cld_snd(n,k),cld_snd_diff(n,k)
     1                           ,wt_snd(n,k),wt_snd_diff(n,k)
 901                  format(i4,4(1x,f19.13))
 904                  continue
                  endif

              enddo ! k

          enddo ! n      
      endif

      if(ncnt.eq.0) then
         write(6,1002)
 1002    format(1x,'no data for barnes: result = first guess')
         t = cf_modelfg
         istatus = 1
         return
      else
         write(6,*)' ncnt/l_diff_snd = ',ncnt,l_diff_snd
      endif

      rr_max = (nx_dim_lut-1)**2 + (ny_dim_lut-1)**2
      iii_max = fnrmscl * rr_max + 1.

      write(6,*)' fnrmscl/iii_max/n_fnorm = '
     1           ,fnrmscl,iii_max,n_fnorm

      if(iii_max .gt. n_fnorm)then
          write(6,*)' iii_max is too large, increase n_fnorm'
     1             ,iii_max,n_fnorm
          istatus = 0
          return
      endif

!     create a lookup table for fnorm(iii)
      do i = -nx_dim_lut,nx_dim_lut
      do j = -ny_dim_lut,ny_dim_lut
          rr=i*i+j*j
          iii=fnrmscl*rr+1.
          if(iii .gt. n_fnorm)iii=n_fnorm
          iiilut(i,j) = fnorm(iii)
      enddo
      enddo

      if(l_diff_snd)then
          sum_a=0d0  
          sumwt_a=0d0
      endif

      do k=1,kmax

        if(k .eq. 1)then
          nstart = 1
        else
          nstart = nlast(k-1) + 1
        endif
        nstop = nlast(k)

        nobs = nstop-nstart+1
        nanl = 0

        if(l_debug)then
            write(6,*)k
     1                 ,sum_a(it1,jt1),sumwt_a(it1,jt1)
!    1                 ,sum_a(it2,jt2),sumwt_a(it2,jt2)
     1                 ,' sum-a'
        endif

        if((l_analyze(k) .and. nobs .ge. 1) .or. k .eq. 1 
     1                                      .or. l_diff_snd)then

          if( (.not. l_diff_snd) .or. k .eq. 1)then
              sum_a=0d0  
              sumwt_a=0d0
          endif

          if(l_diff_snd)then
              n1 = 1
              n2 = n_cld_snd
          else
              n1 = nstart
              n2 = nstop
          endif

          do n=n1,n2
              if(.not. l_diff_snd)then
                ii=iob(n)
                jj=job(n)
                nn=nob(n)
                l_use_snd = .true.
              else
                l_use_snd = .true.

!               test if out of bounds of established perimeter around laps domain
                if(i_snd(n) .lt. ix_low .or. 
     1             i_snd(n) .gt. ix_high) l_use_snd = .false.
                if(j_snd(n) .lt. iy_low .or. 
     1             j_snd(n) .gt. iy_high) l_use_snd = .false.

                ii = i_snd(n)
                jj = j_snd(n)
                nn = n
              endif

              if(l_use_snd)then

                if( (.not. l_diff_snd) .or. k .eq. 1)then

                  nanl = nanl + 1

                  if(.not. l_bterm)then
                      ilow = 1
                      ihigh = imax
                      jlow = 1
                      jhigh = jmax
                  else ! l_bterm
                      if(nskip .eq. 1)then
                          rdist_lim = dist_norm_grid * 6.
                          ilow  = max(int(float(ii)-rdist_lim)  ,1)
                          ihigh = min(int(float(ii)+rdist_lim)+1,imax)
                          jlow  = max(int(float(jj)-rdist_lim)  ,1)
                          jhigh = min(int(float(jj)+rdist_lim)+1,jmax)
                      else ! nskip .gt. 1 
                          ilow = 1
                          ihigh = imax
                          jlow = 1
                          jhigh = jmax
                      endif
                  endif

!                 analyze every few grid points
                  do j=jlow,jhigh,nskip
                  jmjj = j-jj
                  do i=ilow,ihigh,nskip
                      weight = iiilut(i-ii,jmjj) * wt_snd(nn,k) 

!                     obs are being weighted
                      sum_a(i,j)=weight*cld_snd(nn,k)+sum_a(i,j)
                      sumwt_a(i,j)=sumwt_a(i,j)+weight

                  enddo ! i
                  enddo ! j

                  if(l_debug)then
                      weight_dbg=iiilut(it1-ii,jt1-jj)*wt_snd(nn,k)      

                      write(6,*)k,nn
     1                              ,weight_dbg
     1                              ,sum_a(it1,jt1),sumwt_a(it1,jt1)
!    1                              ,sum_a(it2,jt2),sumwt_a(it2,jt2)
     1                              ,' a'

                  endif

                else ! process difference soundings

!                 check if both sounding levels have valid values          
                  if(wt_snd_diff(nn,k) .eq. 0.)then 
                      nanl = nanl + 1

!                     analyze every few grid points
                      do j=1,jmax,nskip
                      jmjj = j-jj
                      do i=1,imax,nskip
                          weight = iiilut(i-ii,jmjj) * wt_snd(nn,k) 

!                         obs are being weighted
                          sum_a(i,j)=weight*cld_snd_diff(nn,k)
     1                                     +sum_a(i,j)      

                      enddo ! i
                      enddo ! j

                      if(l_debug)then
                          write(6,*)k,nn
     1                              ,' b'

                      endif

!                 check if we are changing between valid and missing values
                  elseif(wt_snd_diff(nn,k) .ne. r_missing_data)then

                      nanl = nanl + 1

!                     analyze every few grid points
                      do j=1,jmax,nskip
                      jmjj = j-jj
                      do i=1,imax,nskip
                          weight = iiilut(i-ii,jmjj) * wt_snd_diff(nn,k)

!                         obs are being weighted
                          sum_a(i,j)=weight*cld_snd_diff(nn,k)
     1                                     +sum_a(i,j)      
                          sumwt_a(i,j)=sumwt_a(i,j)+weight

                      enddo ! i
                      enddo ! j

                      if(l_debug)then
                          weight_dbg=
     1                       iiilut(it1-ii,jt1-jj)*wt_snd_diff(nn,k)      

                          write(6,*)k,nn
     1                              ,weight_dbg
     1                              ,sum_a(it1,jt1),sumwt_a(it1,jt1)
!    1                              ,sum_a(it2,jt2),sumwt_a(it2,jt2)
     1                              ,' c'
                      endif

                  endif ! wt_snd_diff

                endif ! l_diff_snd

              endif ! l_use_snd

          enddo ! n

          if(l_debug .and. .false.)then
                write(6,*)k
     1                 ,sum_a(it1,jt1),sumwt_a(it1,jt1)
!    1                 ,sum_a(it2,jt2),sumwt_a(it2,jt2)
     1                 ,' sum-b'
          endif

          write(6,50)k,nstart,nstop,nobs,nanl
50        format(' lvl,nstart,nstop,nobs,nanl=',5i6)

          do j=1,jmax,nskip
          do i=1,imax,nskip

            if(l_diff_snd)then ! qc summations

!             prevent below zero values                 
              sum_a(i,j) = dmax1(sum_a(i,j),0d0) 
              sumwt_a(i,j) = dmax1(sumwt_a(i,j),0d0)

!             prevent sum from exceeding weight (cover > 1)
              sum_a(i,j) = dmin1(sum_a(i,j),sumwt_a(i,j))

            endif

!           add in model first guess as an ob
            sum = sum_a(i,j) + weight_modelfg * cf_modelfg(i,j,k)
            sumwt = sumwt_a(i,j) + weight_modelfg

!           divide weights to get analysis = f(obs + background)
            if (sumwt.eq.0.)then
              t(i,j,k) = r_missing_data
              istatus = 0
            else
              t(i,j,k)=sum/sumwt
            end if

            if(l_debug)then
                if(i .eq. it1 .and. j .eq. jt1)then
                  write(6,*)' total sum',i,j,k,cf_modelfg(i,j,k)
     1                        ,weight_modelfg,sum,sumwt,t(i,j,k)
                endif
                if(i .eq. it2 .and. j .eq. jt2)then
                  write(6,*)' total sum',i,j,k,cf_modelfg(i,j,k)
     1                        ,weight_modelfg,sum,sumwt,t(i,j,k)
                endif
            endif ! l_debug

          enddo ! i
          enddo ! j

          if(l_debug .and. .false.)then
                write(6,*)k
     1                 ,sum_a(it1,jt1),sumwt_a(it1,jt1)
!    1                 ,sum_a(it2,jt2),sumwt_a(it2,jt2)
     1                 ,' sum-c'
          endif

!         bilinearly interpolate to fill in rest of domain
!         fills in final analysis value and weights from obs alone
!         we may have to extrapolate at the n and e edges
          do i = 1,imax
              lowi_lut(i) = (i-1)/nskip*nskip + 1
              il = lowi_lut(i)
              ih = il + nskip
              if(ih .gt. imax)lowi_lut(i) = lowi_lut(i) - nskip
          enddo ! i
          do j = 1,jmax
              lowj_lut(j) = (j-1)/nskip*nskip + 1
              jl = lowj_lut(j)
              jh = jl + nskip
              if(jh .gt. jmax)lowj_lut(j) = lowj_lut(j) - nskip
          enddo ! i

          do j=1,jmax
              jl = lowj_lut(j)
              jh = jl + nskip
              fracj = dble(j-jl)/dble(nskip)

              do i=1,imax

                  if(l_debug .and. .false.)then
                    if(i .eq. it1 .and. j .eq. jt1)then
                        write(6,*)k
     1                         ,sum_a(it1,jt1),sumwt_a(it1,jt1)
     1                         ,' sum-c1'
                    endif
                  endif

                  il = lowi_lut(i)
                  ih = il + nskip
                  fraci = dble(i-il)/dble(nskip)

!                 calculate interpolated cloud cover
                  z1=t(il,jl,k)
                  z2=t(ih,jl,k)
                  z3=t(ih,jh,k)
                  z4=t(il,jh,k)

                  t(i,j,k) =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                - (z2+z4-z3-z1)*fraci*fracj

!                 calculate interpolated ob summation
                  z1=sum_a(il,jl)
                  z2=sum_a(ih,jl)
                  z3=sum_a(ih,jh)
                  z4=sum_a(il,jh)

                  sum_a(i,j) =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                        - (z2+z4-z3-z1)*fraci*fracj

!                 calculate interpolated ob weight summation
                  z1=sumwt_a(il,jl)
                  z2=sumwt_a(ih,jl)
                  z3=sumwt_a(ih,jh)
                  z4=sumwt_a(il,jh)

                  sumwt_a(i,j) =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                        - (z2+z4-z3-z1)*fraci*fracj

              enddo ! i
          enddo ! j

          if(l_debug .and. .false.)then
                write(6,*)k
     1                 ,sum_a(it1,jt1),sumwt_a(it1,jt1)
!    1                 ,sum_a(it2,jt2),sumwt_a(it2,jt2)
     1                 ,' sum-d'
          endif

        elseif(nobs .gt. 0)then ! obs are identical to lvl below; 
                                ! use analysis weights from last analyzed level
          write(6,51)k,nstart,nstop,nobs,nanl
51        format(' lvl,nstart,nstop,nobs,nanl=',5i6,
     1           ' identical obs; copy wts from last analyzed lvl')       

          km1 = k - 1

          do j=1,jmax
          do i=1,imax

!             recover weight summations from last analyzed level below
              sum = sum_a(i,j)
              sumwt = sumwt_a(i,j)

!             add in model first guess as an ob to 'sum'
!             note that sumwt does not need to be modified

              sum   = sum   + weight_modelfg * cf_modelfg(i,j,k)
              sumwt = sumwt + weight_modelfg

!             divide weights to get analysis = f(obs + background)
              if (sumwt.eq.0.)then
                t(i,j,k) = r_missing_data
                istatus = 0
              else
                t(i,j,k)=sum/sumwt
              end if

          enddo ! i
          enddo ! j

        else ! no obs; set level to model first guess
          write(6,52)k,nstart,nstop,nobs,nanl
52        format(' lvl,nstart,nstop,nobs,nanl=',5i6,
     1                  ' no obs; set level to model fg')

          do j=1,jmax
          do i=1,imax
              if(weight_modelfg .gt. 0.)then
                  t(i,j,k) = cf_modelfg(i,j,k)
              else
                  t(i,j,k) = cf_modelfg(i,j,k)
!                 t(i,j,k) = .01
              endif
          enddo ! i
          enddo ! j

        endif

      enddo ! k

      if(istatus .eq. 0)then
          write(6,*)' warning: distant obs given zero weight in barnes_r
     15'
          write(6,*)' try increasing bias_iii'
          write(6,*)' weight_modelfg = ',weight_modelfg
      endif

      deallocate(cld_snd)
      deallocate(wt_snd)
      deallocate(cld_snd_diff)
      deallocate(wt_snd_diff)

      return
      end

