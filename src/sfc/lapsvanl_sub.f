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

c
c=====  here are john's subroutines...(abandon hope, ye who enter)
c
	subroutine channel(u,v,topo,imax,jmax,top,pblht,dx,dy,
     &                     z,div)
c
c=====================================================================
c
c	routine to channel winds around terrain features.
c       includes option to conserve surface convergence from 
c       raw wind data since channeling routine acts to eliminate
c       convergence totally.
c
c       original:   j. mcginley, 2nd half of 20th century
c       changes:    p. stamus  26 aug 1997  changes for dynamic laps
c
c=====================================================================
c
	real u(imax,jmax),v(imax,jmax),z(imax,jmax)
	real dx(imax,jmax),dy(imax,jmax),top(imax,jmax),topo(imax,jmax)
	real div(imax,jmax)
c
	real phi(imax,jmax),ter(imax,jmax),du(imax,jmax),dv(imax,jmax)
	real dpbl(imax,jmax),b(imax,jmax),c(imax,jmax)
c
	call zero(phi,imax,jmax)	! zero the work arrays
	call zero(ter,imax,jmax)
	call zero(du,imax,jmax)
	call zero(dv,imax,jmax)
	call zero(dpbl,imax,jmax)
	call zero(b,imax,jmax)
	call zero(c,imax,jmax)
c
	do j=1,jmax
	do i=1,imax
	  dpbl(i,j) = top(i,j) - topo(i,j)
	enddo !i
	enddo !j
c
	do j=2,jmax-1
	do i=2,imax-1
	  dzx2 = (dpbl(i+1,j)+dpbl(i+1,j+1))*.5
	  dzx1 = (dpbl(i-1,j)+dpbl(i-1,j-1))*.5
	  dzy2 = (dpbl(i-1,j)+dpbl(i,j))*.5
	  dzy1 = (dpbl(i-1,j-1)+dpbl(i,j-1))*.5
	  u2 = u(i,j-1)
	  u1 = u(i-1,j-1)
	  v2 = v(i-1,j)
	  v1 = v(i-1,j-1)
          zbar=(dzx2+dzy2+dzx1+dzy1)*.25
          zbars=zbar*zbar
          b(i,j)=2.*(dzx2-dzx1)/dx(i,j)/zbar
          c(i,j)=2.*(dzy2-dzy1)/dy(i,j)/zbar
	  ter(i,j) = -((u2*dzx2-u1*dzx1)/dx(i,j)
     1                +(v2*dzy2-v1*dzy1)/dy(i,j)) /zbars
     2                -div(i,j)/zbar
	enddo !i
	enddo !j
c
c	print *,' calculating the solution for streamfunction'
	call zero(phi,imax,jmax)
	call leib_2d(phi,ter,100,.1,imax,jmax,z,b,c,z,z,dx,dy)
c
c.....	adjust the winds.
c
	do j=1,jmax-1
	do i=1,imax-1
	  dzx2 = (dpbl(i,j)+dpbl(i,j+1))*.5
	  dzy2 = (dpbl(i,j)+dpbl(i+1,j+1))*.5
	  du(i,j) = (phi(i+1,j+1) - phi(i,j+1)) / dx(i,j) * dzx2
	  u(i,j) = du(i,j) + u(i,j)
	  dv(i,j) = (phi(i+1,j+1) - phi(i+1,j)) / dy(i,j) * dzy2
	  v(i,j) = dv(i,j) + v(i,j)
	enddo !i
	enddo !j
c
	return
	end
c
c
	subroutine frict_2d(fu,fv,u,v,uo,vo,imax,jmax,ak,akk)
c
	real fu(imax,jmax),fv(imax,jmax),u(imax,jmax),v(imax,jmax)
	real vo(imax,jmax),uo(imax,jmax),akk(imax,jmax)
c
	do j=1,jmax
	do i=1,imax
	    uu = u(i,j) * .75 + uo(i,j) * .25
	    vv = v(i,j) * .75 + vo(i,j) * .25
	    fu(i,j) = ak * uu * abs(uu) * akk(i,j)
	    fv(i,j) = ak * vv * abs(vv) * akk(i,j)
	enddo !i
	enddo !j
c
	return
	end
c
c
	subroutine nonlin_2d(nu,nv,u,v,uo,vo,imax,jmax,dx,dy)
c
	real nu(imax,jmax),nv(imax,jmax),u(imax,jmax),v(imax,jmax)
	real uo(imax,jmax),vo(imax,jmax),dx(imax,jmax),dy(imax,jmax)
c
	do j=2,jmax-1
	do i=2,imax-1
	    dudx = (u(i+1,j)-u(i-1,j)) / dx(i,j) * .375 +
     &           (uo(i+1,j)-uo(i-1,j)) / dx(i,j) * .125
	    dvdy = (v(i,j+1)-v(i,j-1)) / dy(i,j) * .375 +
     &           (vo(i,j+1)-vo(i,j-1)) / dy(i,j) * .125
	    dudy = (u(i,j+1)-u(i,j-1)) / dy(i,j) * .375 +
     &           (uo(i,j+1)-uo(i,j-1)) / dy(i,j) * .125
	    dvdx = (v(i+1,j)-v(i-1,j)) / dx(i,j) * .375 +
     &           (vo(i+1,j)-vo(i-1,j)) / dx(i,j) * .125
	    uu = (u(i,j)+u(i,j-1)+u(i+1,j)+u(i+1,j-1)) * .1875 +
     &         (uo(i,j)+uo(i,j-1)+uo(i+1,j)+uo(i+1,j-1)) * .0675
	    vv = (v(i,j)+v(i-1,j)+v(i,j+1)+v(i-1,j+1)) * .1875 +
     &         (vo(i,j)+vo(i-1,j)+vo(i,j+1)+vo(i-1,j+1)) * .0675
	    utt = u(i,j) * .75 + uo(i,j) * .25
	    vtt = v(i,j) * .75 + vo(i,j) * .25
	    nu(i,j) = utt * dudx + vv * dudy
	    nv(i,j) = uu * dvdx + vtt * dvdy
	enddo !i
	enddo !j
c
	return
	end
c
c
	subroutine leib_2d(sol,force,itmax,erf,imax,jmax,a,b,c,d,e,
     &                   dx,dy)
c
c.....  relaxation routine.
c.....  changes:  p. stamus, noaa/fsl  13 aug 1999
c.....                 cleaned up, added tagit routine.
c       
	real sol(imax,jmax),force(imax,jmax),a(imax,jmax)
	real b(imax,jmax),c(imax,jmax),d(imax,jmax),e(imax,jmax)
	real dx(imax,jmax), dy(imax,jmax)
c
	call tagit('leib_2d', 19990813)
	ovr = 1.
	reslmm = 0.
	erb = 0.
c  first guess here
	ittr = 0
	do 1 it=1,itmax
	  ertm = 0.
	  ermm = 0.
	  ia = 0
	  corlm = 0.
	  do 2 j=2,jmax-1
	  do 2 i=2,imax-1
	    dx2 = dx(i,j) * 2.
	    dxs = dx(i,j) * dx(i,j)
	    dy2 = dy(i,j) * 2.
	    dys = dy(i,j) * dy(i,j)
	    aa = a(i,j)
	    bb = b(i,j)
	    cc = c(i,j)
	    dd = d(i,j)
	    cortm = (-2. / dxs) - (2. / dys) + e(i,j)
20	    res = (sol(i+1,j) + sol(i-1,j)) / dxs +
     &            (sol(i,j+1) + sol(i,j-1)) / dys +
     &            (cortm * sol(i,j)) - force(i,j) + bb * 
     &            ((sol(i+1,j) - sol(i-1,j)) / dx2) + cc * 
     &            (sol(i,j+1) - sol(i,j-1)) / dy2
	    cor = res / cortm
	    if(abs(cor) .gt. erf) ia = 1
	    if(abs(cor) .gt. corlm) corlm = abs(cor)
	    sol(i,j) = sol(i,j) - cor * ovr
2	  continue
5	  ittr = ittr + 1
	  cor5 = corlm
	  if(ittr .ne. 5) go to 15
	  ittr = 0
	  rho = (cor5 / cor0) ** .2
	  if(rho .gt. 1) go to 16
	  ovr = 2. / (1. + sqrt(1. - rho))
16	  continue
	  cor0 = cor5
15	  if(ia .ne. 1) go to 4
	  if(it .ne. 1) go to 1
	  corlmm = corlm
	  cor0 = corlmm
1	continue
4	continue
	reslm = corlm * cortm
	write(6,1001) it,reslm,corlm,corlmm,erb
	write(6,1002) ovr
1002	format(1x,'ovr rlxtn const at fnl ittr = ',e10.4)
1001	format(1x,'iterations= ',i4,' max residual= ',e10.3,
     & ' max correction= ',e10.3, ' first iter max cor= ',e10.3,
     & 'max bndry resid= ',e10.3)
c
	return
	end
c
c
	subroutine spline(t,to,tb,alf_in,alf2a_in,beta_in,a_in,s_in,
     &                    cormax,err,imax,jmax,rms_thresh_norm,bad_mult,
     &                    imiss,mxstn,obs_error,name,topo,ldf,wt_bkg_a)
c
c*******************************************************************************
c	laps spline routine...based on one by j. mcginley.
c
c	changes:  
c	  p. stamus	10-18-90  started to clean code. made alf2/alf2o arrays.
c			11-11-91  pass in dummy work arrays.
c			07-27-93  changes for new barnes2 routine. 
c                       07-20-95  put wt calcs here...call to dynamic_wts.
c                       08-26-97  changes for dynamic laps. pass in obs_error.
c         j. mcginley   09-22-98  changes to fully exploit background info.
c          and p.stamus           routine modified to always use a background
c                                 field for qc and 1st guess.  spline solves for
c                                 a solution difference from 1st guess.  removed
c                                 2 barnes calls. when no satellite data for hsm,
c                                 'a' weight set to zero.  bkg added back into
c                                 t, to, and s arrays on exit.
c         p. stamus     09-29-98  calc std dev from just obs (not boundaries+obs)
c                       01-28-99  temp. replace spline section with barnes. fix
c                                   boundary normalization.
c                       07-24-99  turn spline back on, rm barnes.  adj weights.
c                                   turn satellite back on.
c                       08-13-99  change call to allow diff alf/alf2a/beta/a for
c                                   diff variables.  rm alf2 as array.
c                       11-23-99  put background (tb) into output (t) array if
c                                   no obs or all obs bad in data array (to).
c
c         s. albers     2001      adding in data structures towards the goal of
c                                 directly analyzing obs external to the grid
c*******************************************************************************
c
        logical l_barnes_wide, l_struct
        data l_struct /.false./        ! using data structures?
        data l_barnes_wide /.true./    ! using barnes_wide routine on boundary?

        integer max_obs
        parameter (max_obs = 40000)       

        include 'barnesob.inc'
        type (barnesob_qc) obs_barnes(max_obs)

	real t(imax,jmax), to(imax,jmax), s(imax,jmax), s_in(imax,jmax)
	real ress(1000), tb(imax,jmax) !, alf2(imax,jmax)
        real topo(imax,jmax),ldf(imax,jmax)
        real wt_bkg_a(imax,jmax)                         
c
	real fnorm(0:imax-1,0:jmax-1)
c	real alf2o(imax,jmax)  !work array
c
	character name*10
	logical iteration
        logical l_boundary(imax,jmax)

        integer analysis_mode

        lun_s = 6
c
!       analysis mode controls the strategy of doing the analysis
!
!       1) method in place early 2001. start spline with constant field based
!          on mean of observation increments.
!
!       2) start spline with incremental field as analyzed by 
!          'barnes_multivariate' routine. 
!
!       3) use barnes_multivariate routine as complete substitute to spline.

        analysis_mode = 3

 	write(6,910)analysis_mode
 910	format(' subroutine spline: analysis_mode = ',i3)

!       fill data structure
        rinst_err = 1.5
        n_obs = 0

        if(l_struct)then ! define "to" array from data structure
            to = 0.0     ! f90 assignment
            do iob = 1,n_obs
                i = obs_barnes(iob)%i
                j = obs_barnes(iob)%j
                to(i,j) = obs_barnes(n_obs)%value(1)
            enddo ! iob

        else             ! define data structure from input "to" array
            do i = 1,imax
            do j = 1,jmax
                if(to(i,j) .ne. 0.0)then
                    n_obs = n_obs + 1
                    obs_barnes(n_obs)%i = i
                    obs_barnes(n_obs)%j = j
                    obs_barnes(n_obs)%k = 1
                    obs_barnes(n_obs)%weight   = 1./rinst_err**2       
                    obs_barnes(n_obs)%vert_rad_rat = 1.0
                    obs_barnes(n_obs)%value(1) = to(i,j)
                    obs_barnes(n_obs)%qc = .true.
c
c th: 29 november 2002 begin hack.
c
                    if (name .eq. 'pressure') then
                       obs_barnes(n_obs)%mask_sea = 0
                    else
                       obs_barnes(n_obs)%mask_sea = 1
                    endif
c
c th: 29 november 2002 end hack.
c
                endif
            enddo ! j
            enddo ! i
        endif

	imiss = 0
	ovr = 1.4
	itmax = 1000	! max number of iterations
	zeros = 1.e-30
	smsng = 1.e37
	cormax = 1.
	call move(s_in, s, imax, jmax)
	call zero(t, imax,jmax)
c
c.....	first guess use barnes
c
	npass = 1
c
c.....  count the number of observations in the field (not counting the 
c.....  boundaries.

c
        if(l_barnes_wide)then
            i_boundary = 0 ! 2
        else
            i_boundary = 0
        endif

        ilow = i_boundary + 1
        ihigh = imax-i_boundary 
        jlow = i_boundary + 1
        jhigh = jmax-i_boundary

        n_obs_var = 0
        if(l_barnes_wide)l_boundary = .true. ! f90 assignment

        if(l_struct)then
            if(l_barnes_wide)then
                do iob = 1,n_obs
                    i = obs_barnes(iob)%i
                    j = obs_barnes(iob)%j
                    if(i.ge.ilow .and. i.le.ihigh .and. 
     1                 j.ge.jlow .and. j.le.jhigh      )then ! interior point
                        l_boundary(i,j) = .false.            
                        n_obs_var = n_obs_var + 1
                    endif
                enddo ! iob
            else
                n_obs_var = n_obs
            endif
        else
	    do j=jlow,jhigh
	    do i=ilow,ihigh
                if(l_barnes_wide)l_boundary(i,j) = .false.  ! interior point
  	        if(to(i,j) .ne. 0.) n_obs_var = n_obs_var + 1
	    enddo !i
	    enddo !j
        endif

        n_valid_obs = n_obs

	if(n_obs_var .eq. 0) then
	  print *,'  warning.  no observations found in data array. '
	  imiss = 1
	  go to 950
	else
	   print *,'  observations in data array: ', n_obs_var
	endif
c
	if(name.ne.'noplot' .and. name(1:3).ne.'tb8') then
	  write(lun_s,912)
 912	  format('  data passed into spline:')
	endif
c
c.....	data check algorithm and computation of difference 
c.....  from background.
c
	sum = 0.
	cnt = 0.
	sum1 = 0.
	icnt = 0
c
c.....	compute standard deviation of the obs
c
	do j=jlow,jhigh
	do i=ilow,ihigh
	  if(to(i,j) .eq. 0.) go to 99
	  sum = sum + ((to(i,j) - tb(i,j)) ** 2)
	  cnt = cnt + 1.
99	continue
        enddo !i
	enddo !j
c
	if(cnt .eq. 0.) then
	  print *,'  warning.  zero observations found in data array. '
	  go to 950
	else
	  std = sqrt(sum / cnt)
	endif
	if(std .eq. 0.) then
	  write(6,927) 
 927	  format(1x,'  warning. standard deviation is zero.',
     &           ' observations equal backgroud at all locations.')
	  std = zeros
	endif
c
c.....  bad data defined as deviating 'bad_mult' sigmas from 1st guess
c
	bad = bad_mult * std
	print *,' std dev: ',std,', bad value: ',bad
     1                          ,', ratio: ',bad_mult      
c
c.....  normalize the obs with respect to the bkg.
c
	sumdif = 0.
	numdif = 0

        if(l_struct)then
            do iob = 1,n_obs
                i = obs_barnes(iob)%i
                j = obs_barnes(iob)%j
                if(i.ge.ilow .and. i.le.ihigh .and. 
     1             j.ge.jlow .and. j.le.jhigh      )then ! interior point

!                   eliminate bad data from the interior while normalizing
!                   to the background
                   

                else ! normalize boundary ob
 	            to(i,j) = to(i,j) - tb(i,j)
 
                endif

            enddo ! iob

        else ! .not. l_struct
          do j = 1,jmax
          do i = 1,imax
            if(l_boundary(i,j) .and. l_barnes_wide)then ! normalize boundary ob

!              here it seems to be important that the stations are mapped onto
!              the grid with rounding up allowed to get the best possible 
!              departures

	       to(i,j) = to(i,j) - tb(i,j)

            else ! eliminate bad data from the interior while normalizing to
                 ! the background

	       if(to(i,j) .ne. 0.) then
	           diff = to(i,j) - tb(i,j)
	           if(abs(diff) .lt. bad) then
	               to(i,j) = diff
	               sumdif = sumdif + diff
	               numdif = numdif + 1
	           else
	               write(6,1099) i,j,to(i,j), diff
	               to(i,j) = 0.
	               n_obs_var = n_obs_var - 1
	           endif

               endif
 
            endif ! boundary point

          enddo ! j
          enddo ! i

        endif ! l_struct       
 
1099	format(1x,'bad data at i,j ',2i5,': value ',e12.4
     1        ,', diff ',e12.4)
c
	print *,' '
	if(numdif .ne. n_obs_var) then
	   print *,' hmmmm...numdif= ',numdif,' ; n_obs_var= ',n_obs_var
	endif
	if(n_obs_var.gt.0 .and. numdif.gt.0) then
	   print *,
     &       ' observations in data array after spline qc: ',n_obs_var
	else
	   print *,
     &       '  warning. no observations in data array after qc check.'
	   go to 950
	endif


c.....  subtract background from satellite obs and calculate sat stats (deg f)
c
	isat_flag = 0
        num_sat_bkg = 0
        sum_sat_bkg = 0.
        sumsq_sat_bkg = 0.
        num_sat_obs = 0
        sum_sat_obs = 0.
        sumsq_sat_obs = 0.

	do j=1,jmax
	do i=1,imax
	   if(s(i,j) .ne. 0.) then
	      isat_flag = 1
	      s(i,j) = s(i,j) - tb(i,j)           ! sat diff from background
              num_sat_bkg = num_sat_bkg + 1
              sum_sat_bkg = sum_sat_bkg + s(i,j)
              sumsq_sat_bkg = sumsq_sat_bkg + s(i,j)**2

              if(to(i,j) .ne. 0)then
                 diff_sat_obs = s(i,j) - to(i,j)  ! sat diff from obs
                 num_sat_obs = num_sat_obs + 1
                 sum_sat_obs = sum_sat_obs + diff_sat_obs
                 sumsq_sat_obs = sumsq_sat_obs + diff_sat_obs**2
              endif
	   endif
	enddo !i
	enddo !j

        write(6,*)' num sat tb8 points = ',num_sat_bkg

        if(num_sat_bkg .gt. 1)then
            rmean = sum_sat_bkg   / float(num_sat_bkg)
            rms   = sqrt(sumsq_sat_bkg / float(num_sat_bkg))
            rms_abm = sqrt(rms**2 - rmean**2)
            write(6,*)
     1          ' mean/rms/rms about mean: sat vs. background = '
     1          ,rmean,rms,rms_abm       

            write(6,*)' num of sat/obs comparisons = ',num_sat_obs

            if(num_sat_obs .gt. 1)then
                rmean = sum_sat_obs   / float(num_sat_obs)
                rms   = sqrt(sumsq_sat_obs / float(num_sat_obs))
                rms_abm = sqrt(rms**2 - rmean**2)
                write(6,*)
     1              ' mean/rms/rms about mean: sat vs. obs = '
     1              ,rmean,rms,rms_abm       
            endif

        endif
c
c.....  have obs, so set starting field so spline converges faster.
c
        if(analysis_mode .eq. 1)then
	    amean_start = sumdif / numdif
	    print *,' using mean of '
     1             , amean_start, ' to start analysis.'
	    call constant(t, amean_start, imax,jmax)

        elseif(analysis_mode .ge. 2)then
            write(6,*)' calling barnes_multivariate_sfc to start spline'       

            call get_fnorm_max(imax,jmax,r0_norm,r0_value_min,fnorm_max)   
            n_fnorm = int(fnorm_max) + 1

            call barnes_multivariate_sfc(to,imax,jmax               ! inputs
     1                                ,smsng                        ! input
     1                                ,rms_thresh_norm              ! input
     1                                ,rinst_err                    ! input
     1                                ,bad                          ! input
     1                                ,wt_bkg_a                     ! input
     1                                ,n_fnorm                      ! input
     1                                ,l_boundary,.true.,.false.    ! input
     1                                ,n_valid_obs,obs_barnes       ! input
     1                                ,topo,ldf                     ! input
     1                                ,t                            ! output
     1                                ,istatus)                     ! output
        endif
c
cc	print *,' using smooth barnes to start the analysis.'
cc	rom2 = 0.005
cc	npass = 1
cc	idum = 0
cc	call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
cc	call barnes2(t,imax,jmax,to,smsng,idum,npass,fnorm)
cc	print *,' done.'
c       
        if(analysis_mode .le. 2)then
c
c.....      set the weights for the spline.
c
c           we need parity with respect to obs and filtering. since number
c           of obs can change depending on variable adjust alf accordingly
c           so that it is representative of the inverse of observation
c           error**2 times the number of working gridpoints divided by the 
c           number of ob points so that the term is roughly comprable to the 
c           beta term
	    alf = beta_in*(imax-4)*(jmax-4)/n_obs_var
	    alf2a = alf2a_in
	    beta = beta_in
	    a = a_in
	    if(isat_flag .eq. 0) a = 0.
c
	    write(6,9995) alf, beta, alf2a, a
 9995	    format(5x,'using spline wts: alf, beta, alf2a, a = ',4f12.4)       
c
c.....      now do the spline.
c
	    iteration = .true.

	    do it=1,itmax
	      cormax = 0.
	      if(iteration) then
c	         print *,' it = ', it
	         do j=3,jmax-2
	         do i=3,imax-2
		    alfo = alf
		    alf2o = alf2a
		    ao = a
c
		    if(to(i,j) .eq. 0.) alfo = 0.
		    if(s(i,j).eq.0. .or. s(i-1,j).eq.0. 
     &                              .or. s(i+1,j).eq.0.
     &                 .or. s(i,j-1).eq.0. .or. s(i,j+1).eq.0.) then       
		       ao = 0.
		       sxx = 0.
		       syy = 0.
		    else
		       sxx = (s(i+1,j) + s(i-1,j) - 2. * s(i,j))
		       syy = (s(i,j+1) + s(i,j-1) - 2. * s(i,j))
		    endif
c
		    dtxx = t(i+1,j) + t(i-1,j) - 2. * t(i,j)
		    dtyy = t(i,j+1) + t(i,j-1) - 2. * t(i,j)
		    d4t = 20. * t(i,j) 
     &                - 8. * (t(i+1,j) + t(i,j+1) + t(i-1,j) + t(i,j-1))
     &                + 2.*(t(i+1,j+1)+t(i+1,j-1) 
     &                + t(i-1,j+1) + t(i-1,j-1))
     &                + (t(i+2,j) + t(i-2,j) + t(i,j+2) + t(i,j-2))
		    d2t = dtxx + dtyy
		    dtx = (t(i+1,j) - t(i-1,j)) * .5
		    dty = (t(i,j+1) - t(i,j-1)) * .5
c       
		    res = d4t - ao * (d2t - sxx - syy) / beta
     &                   + alfo/beta * (t(i,j) - to(i,j)) ! stations
     &                   + alf2o/beta * t(i,j)            ! background
		    cortm = 20. + ao*4./beta + alfo/beta + alf2o/beta
		    tcor = abs(res / cortm)
		    t(i,j) = t(i,j) - res / cortm * ovr
		    if(tcor .le. cormax) go to 5
		    cormax = tcor
		    ress(it) = tcor
c	            write(6,1010) i,j,res,cortm,tcor
c1010	            format(1x,2i5,3e12.4)
c	            write(6,1009)beta,d4t,d2t,dtxy,dtx,dty,gam,sxx
c    1                          ,syy,sxy,sx,sy
 5		    continue
	         enddo !i
	         enddo !j
c
c	         write(6,1000) it,cormax
	         if(cormax .lt. err) iteration = .false.
	         corhold = cormax
	         ithold = it
	      endif ! iteration
	    enddo !it
c
cc	    do j=1,jmax
cc	       do i=1,imax
cc	          if(to(i,j) .ne. 0.) then
cc		     write(6,7119) i,j,to(i,j)
cc	          endif
cc	       enddo
cc	    enddo

 7119	    format(2i5,f10.2)

  	    write(6,1000) ithold ,corhold !it, cormax
 1000	    format(1x,' it/cormax= ',i4,e12.4)

        endif ! do the spline

cc	return

 876	continue
c
c.....  add backgrounds back to t, to, and s
c
	do j=1,jmax
	do i=1,imax
	   t(i,j) = t(i,j) + tb(i,j)
	   if(s(i,j)  .ne. 0.)  s(i,j) =  s(i,j) + tb(i,j)
	   if(to(i,j) .ne. 0.) to(i,j) = to(i,j) + tb(i,j)
	enddo !i
	enddo !j
c
	if(name.ne.'noplot' .and. name(1:3).ne.'tb8') then
	   write(lun_s,923)
 923	   format(1x,' solution after spline')
	endif
	if(cormax.lt.err .and. it.eq.1) return
 3	continue
c       
c.....  that's all.
c
!	print *,' leaving spline'
	return
c
 950	continue
	print *,
     &  '    no observations available. setting analysis to background.'
	call move(tb, t, imax, jmax)
        write(6,*)'    bkgnd analysis range is ',minval(t),maxval(t)
	return
c
	end
c
c
	subroutine barnes2(t,imax,jmax,to,smsng,mxstn,npass,fnorm)
c
	real to(imax,jmax), t(imax,jmax), val(imax*jmax)
	real fnorm(0:imax-1,0:jmax-1)
	real h1(imax,jmax), h2(imax,jmax)  !work arrays
c
	integer iob(imax*jmax), job(imax*jmax), dx, dy
c
c
	call zero(h1,imax,jmax)
	call zero(h2,imax,jmax)
c
	badd = 1.e6 - 2.	! bad data value
c
c.....	loop over field npass times 
c
!	print *,' *** in barnes2 ***'
        ncnt = 0
        do j=1,jmax
        do i=1,imax
          if (to(i,j).ne.0. .and. to(i,j).lt.badd) then
            ncnt = ncnt + 1
            iob(ncnt) = i
            job(ncnt) = j 
            val(ncnt) = to(i,j)
cc	     write(6,999) ncnt, i, j, to(i,j)
          endif
        enddo !i
        enddo !j
	if(ncnt .eq. 0) then
	  print *,' *** ncnt = 0 in barnes2. ***'
	  return
	endif
 999	format('   ncnt: ',i4,' at ',2i5,f10.3)
c
	do ipass=1,npass
c
	  do j=1,jmax
	  do i=1,imax
	    sum = 0.
	    sumwt = 0.
	    sum2 = 0.
	    sumwt2 = 0.
	    do n=1,ncnt
	      dy = abs(j - job(n))
	      dx = abs(i - iob(n))
	      sum2 = fnorm(dx,dy) * val(n) + sum2
	      sumwt2 = sumwt2 + fnorm(dx,dy)
	    enddo !n
c
	    if(sumwt2 .eq. 0.) then
	      if(ipass .eq. 1) then
c		   print *,' got into wierd loop.............'
	            sum2 = 0.
	            sumwt2 = 0.
	            do n=1,ncnt
	              dx = abs(i - iob(n))
	              dy = abs(j - job(n))
	              sum2 = (fnorm(dx,dy) + .01) * val(n) + sum2
	              sumwt2 = sumwt2 + (fnorm(dx,dy) + .01)
	            enddo !n
	            if(sumwt2 .ne. 0.) go to 490
	      else
	        go to 500
	      endif
	    endif
c
490	    continue
	    t(i,j) = sum2 / sumwt2
c
500	  continue
          enddo !i
	  enddo !j
c
	  if(ipass .eq. 2) then
	    call move(h2,to,imax,jmax)
	    call diff(h1,t,t,imax,jmax)
!	    write(lun_s,915)
!915	    format(' after 2nd barnes pass')
!	    write(lun_s,909) rom2
!909	    format(' radm2 = ',f8.4)
	    go to 550
	  endif
!	  write(lun_s,912)
912	  format(' after 1st barnes pass')
	  if(npass .eq. 1) go to 550
	  call move(t,h1,imax,jmax)
 	  call move(to,h2,imax,jmax)
	  do n=1,ncnt
	    val(n) = t(iob(n),job(n)) - val(n)  
	  enddo !n
550	continue
        enddo !ipass
c
!	print *,' *** barnes2 done. ***'
	return
	end
c
c
	subroutine barnes_wide(t,imax,jmax,ii,jj,t_ob,numsta,smsng,
     &                         mxstn,npass,fnorm,istatus)
c
c.....	routine to do a barnes analysis that will consider stations in
c.....	the 't_ob' array that are outside the boundaries of the 't' array.
c
c       changes:  p.stamus noaa/fsl  7 jan 1999  add status flag.
c

	real t(imax,jmax), t_ob(mxstn) 
	real fnorm(0:imax-1,0:jmax-1)
	real h1(imax,jmax), val(mxstn)
c	
	integer iob(mxstn), job(mxstn), ii(mxstn), jj(mxstn)
	integer dx, dy 
c

        lun_s = 6

!	print *,' *** in barnes_wide ***'
	istatus = -1
	call zero(h1,imax,jmax)
	im1 = imax - 1
	jm1 = jmax - 1
c
c.....	loop over field npass times 
c
	ncnt = 0
	do n=1,numsta
          if (t_ob(n).ne.0. .and. t_ob(n).ne.smsng) then
	    ncnt = ncnt + 1
	    iob(ncnt) = ii(n)
	    job(ncnt) = jj(n)
	    val(ncnt) = t_ob(n)
          endif
	enddo !n 
c
	if(ncnt .eq. 0) then
	   print *,' **warning. no obs for analysis in barnes_wide. **'
	   istatus = 0
	   return
	endif
c
	write(6,900) ncnt, numsta
900	format('   selected ',i5,' obs out of ',i5
     1        ,' total (barnes_wide)')
c
	do ipass=1,npass
c
	  do j=1,jmax
	  do i=1,imax
	    sum2 = 0.
	    sumwt2 = 0.
	    do n=1,ncnt
	      dy = min(abs(j - job(n)), jm1) 
	      dx = min(abs(i - iob(n)), im1) 
	      sum2 = fnorm(dx,dy) * val(n) + sum2
	      sumwt2 = sumwt2 + fnorm(dx,dy)
	    enddo !n
	    if(sumwt2 .eq. 0.) then
	      if(ipass .eq. 1) then
!               print *,' barneswide weird loop....'

                print *,' barneswide weird loop....skip rest of routine'
                istatus = 0
                if(.true.)return

	        sum2 = 0.
	        sumwt2 = 0.
	        do n=1,ncnt
	          dx = min(abs(i - iob(n)), im1) 
	          dy = min(abs(j - job(n)), jm1) 
	          sum2 = (fnorm(dx,dy)+.01) * val(n) + sum2
	          sumwt2 = sumwt2 + (fnorm(dx,dy) + .01) 
	        enddo !n
	      else
	        go to 500
	      endif 
	    endif 
c
	    if(sumwt2 .ne. 0.) t(i,j) = sum2 / sumwt2
c       
 500	    continue
	 enddo !i
  	 enddo !j
c
	  if(ipass .eq. 2) then
	    call diff(h1,t,t,imax,jmax)
!	    write(lun_s,915)
!915	    format(' after 2nd barnes pass')
!	    write(lun_s,909) rom2
!909	    format(' radm2 = ',f8.4)
	    go to 550
	  endif
!	  write(lun_s,912)
912	  format(' after 1st barnes pass')
	  if(npass .eq. 1) go to 550
	  call move(t,h1,imax,jmax)
c
	  do n=1,ncnt
	    if(iob(n).lt.1 .or. iob(n).gt.imax  .or.
     &         job(n).lt.1 .or. job(n).gt.jmax) then
	      val(n) = 0. ! which is 't_ob(n)-val(n)' at stns outside the grid
	    else
	      val(n) = t(iob(n),job(n)) - val(n)
	    endif
	  enddo !n
550	continue
        enddo !ipass
c
!	print *,'   leaving barnes_wide'
	istatus = 1
	return
	end
c
	subroutine make_cssi(t,td,pmsl,u,v,cssi,ni,nj,badflag)
c
c======================================================================
c
c	routine to calculate the cssi (rodgers and maddox 81) at each 
c       laps gridpoint.  the temp and dewpt enter in deg f, the msl 
c       pressure in mb, and the wind components in m/s, which have to be 
c       converted to speed and direction in kts.
c
c	original version: 05-03-91  peter a. stamus noaa/fsl
c	changes:          11-11-91  pass in dummy arrays.
c                         08-26-97  changes for dynamic laps
c
c======================================================================
c
	real t(ni,nj), td(ni,nj), pmsl(ni,nj), u(ni,nj), v(ni,nj)
	real cssi(ni,nj)
c
	real spd(ni,nj), dir(ni,nj)  !work arrays
c
c.....	start.  convert u,v in m/s to spd/dir in kts.
c
	call windconvert(u,v,dir,spd,ni,nj,badflag)
	call conv_ms2kt(spd,spd,ni,nj)	
c
c.....	calculate each of the 4 terms involved, then combine.
c
	do j=1,nj
	do i=1,ni
	  term1 = t(i,j) - 60.		! temperature
	  term2 = 2. * (td(i,j) - 45.)	! moisture
	  term3 = abs(1010.0 - pmsl(i,j))  ! pressure: abs of diff 
	  if(dir(i,j).gt.180. .and. dir(i,j).lt.360.) then	! west wind
	    term4 = -2. * spd(i,j)
	  else							! east wind
	    if(td(i,j) .ge. 45.) then				! that's moist
	      term4 = 2. * spd(i,j)
	    else						! that's not...
	      term4 = spd(i,j)
	    endif
	  endif
c
	  cssi(i,j) = term1 + term2 - term3 + term4
c
	enddo !i
	enddo !j
c
	return
	end
c
c
	subroutine windconvert(uwind,vwind,direction,speed,
     &                         ni,nj,badflag)
c
c======================================================================
c
c       given wind components, calculate the corresponding speed and 
c       direction.  hacked up from the windcnvrt_gm program.
c
c
c       argument     i/o   type       description
c      --------	     ---   ----   -----------------------------------
c       uwind         i    r*4a    u-component of wind
c       vwind         i	   r*4a    v-component of wind
c       direction     o    r*4a    wind direction (meteoro. degrees)
c       speed         o    r*4a    wind speed (same units as input)
c       ni,nj         i    i       grid dimensions
c       badflag       i    r*4     bad flag value
c
c       notes:
c       1.  if magnitude of uwind or vwind > 500, set the speed and 
c           direction set to the badflag value.
c
c       2.  units are not changed in this routine.
c
c======================================================================
c
	real  uwind(ni,nj), vwind(ni,nj)
	real  direction(ni,nj), speed(ni,nj)
c
	do j=1,nj
	do i=1,ni
	   if(abs(uwind(i,j)).gt.500. .or. 
     &                           abs(vwind(i,j)).gt.500.) then
	      speed(i,j) = badflag
	      direction(i,j) = badflag
c
	   elseif(uwind(i,j).eq.0.0 .and. vwind(i,j).eq.0.0) then
	      speed(i,j) = 0.0
	      direction(i,j) = 0.0			!undefined
c
	   else
	      speed(i,j) = 
     &          sqrt(uwind(i,j)*uwind(i,j) + vwind(i,j)*vwind(i,j))  !speed
	      direction(i,j) = 
     &          57.2957795 * (atan2(uwind(i,j),vwind(i,j))) + 180.   !dir
	   endif
	enddo !i
	enddo !j
c
	return
	end
c
c
c
	subroutine dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
c
c=====================================================================
c
c     routine to calculate the weights to be used in the barnes
c     analysis.  the data density is used to set the cutoff for
c     the response function.  then that cutoff is used to calculate
c     the exp, based on differences so that no additional distance
c     calculations are required in the barnes routine.  all of this
c     is done in gridpoint space.
c
c     original:  07-14-95  p. stamus, noaa/fsl
c     changes:   p. stamus  08-28-97  declare dx,dy integers.
c                           09-10-97  bag include. pass in fnorm.
c
c     notes:
c
c       1.  if variable 'rom2' is passed in as zero, it is calculated
c           from the data density.  otherwise, the value passed in is
c           used in the weight calculation.
c
c       2.  the response for 2d waves is hard-wired in this routine.
c           this is the 'con' variable, and comes from the eqn:
c                     d = exp -(pi**2 r**2)/lamba**2
c           if we set d (the response) to our desired cutoff, set 
c           lamba to the desired wavelength in gridpt space (2d),
c           then solve for r in terms of d, we get the 'con' value
c           (i.e.,  r = (con)*d).  here are some values for different
c           cutoffs:
c                     d = 0.01     r = 1.36616d
c                         0.10     r = 0.96602d
c                         0.25     r = 0.74956d
c                         0.50     r = 0.53002d
c
c=====================================================================
c
	real fnorm(0:imax-1,0:jmax-1)
	integer dx, dy
c
c.... first, find the area that each ob covers in gridpt space (this
c.... of course assumes a uniform coverage).
c
cc	con = 0.96602     ! resp of 0.10
cc	con = 0.74956     ! resp of 0.25
	con = 0.53002     ! resp of 0.50
	if(rom2 .eq. 0.) then
	   area = float(imax * jmax) / n_obs_var
	   d = sqrt( area )
	   rom2 = 1. / ((con * con) * (d * d))
	   write(6,900) n_obs_var, area, d, rom2
 900	   format(1x,'num obs: ',i5,'  area: ',f8.2,'  d: ',f8.2,
     &       '  rom2: ',f8.5)
	else
	   d = sqrt(1./(con * con * rom2))
	   write(6,902) rom2, d
 902	   format(' using preset rom2 of: ',f8.5,'  calc d: ',f8.2)
	endif
c
c.... now calculate the weights for all the possible distances.
c
	pi = 4. * atan(1.)
	fno = 1. / (sqrt(2. * pi))
c
	do dy=0,jmax-1
	do dx=0,imax-1
	   rr = float(dx*dx + dy*dy)
	   fnorm(dx,dy) = fno * (exp( -(rr * rom2)))
	enddo !dx
	enddo !dy
c
c.... that's it.
c
	return
	end
c
c
	subroutine calc_beta(d,obs_error,beta)
c
c=======================================================================
c
c       routine to calculate the 'beta' weight for the spline.  'beta'
c       is calculated based on the 'd' from the gridpt to data distance
c       and an expected observation error for the partictular variable.
c
c
c       original:  07-19-95  p. stamus, noaa/fsl
c       changes:
c
c=======================================================================
c
	pi = 4. * atan( 1. )
	pi4 = pi * pi * pi * pi
	d4 = d * d * d * d
	alpha = -99.9
	beta = -99.9
c
	if(obs_error .ne. 0.) then
	   alpha = 1. / (obs_error * obs_error)
	else
	   print *,' **error. obs_error = 0 in calc_beta.**'
	   go to 100
	endif
c
	beta = alpha * d4 / pi4
c
 100	continue
	write(6,900) obs_error, d, alpha, beta
 900	format(1x,'obs error: ',f9.4,'  d: ',f9.4,
     &                      '  alpha: ',f9.4,'  beta: ',f15.4)
	if(beta .eq. 0.) then
	   print *,' **error. beta = 0 in calc_beta.**'
	   beta = -99.9
	endif
c
	return
	end

