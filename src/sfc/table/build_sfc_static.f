cdis
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
c
        program locpost
        integer istatus

        real, allocatable, dimension(:,:) :: lat
        real, allocatable, dimension(:,:) :: lon
        real, allocatable, dimension(:,:) :: topo
        real, allocatable, dimension(:,:) :: ldf

        call get_grid_dim_xy(nx_l, ny_l, istatus)
        if (istatus .ne. 1) then
            write(6,*) 'return get_grid_dim_xy, status: ', istatus
            stop
        endif

!       allocate static arrays (lat, lon, topo, ldf)
        allocate( lat(nx_l,ny_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate lat'
            stop
        endif

        allocate( lon(nx_l,ny_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate lon'
            stop
        endif

        allocate( topo(nx_l,ny_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate topo'
            stop
        endif

        allocate( ldf(nx_l,ny_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate ldf'
            stop
        endif

!       read static arrays (lat, lon, topo, ldf)
        call read_static_grid(nx_l,ny_l,'lat',lat,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps lat'
            stop
        endif

        call read_static_grid(nx_l,ny_l,'lon',lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps lon'
            stop
        endif

        call read_static_grid(nx_l,ny_l,'avg',topo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps avg'
            stop
        endif

        call read_static_grid(nx_l,ny_l,'ldf',ldf,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps ldf'
            stop
        endif

        call build_sfc_static_sub(nx_l,ny_l,lat,lon,topo)

        lun = 31
        lun_out = 19

        call locpost_radar(nx_l,ny_l,lat,lon,topo,ldf
     1                    ,lun,lun_out,istatus)

        end

	subroutine build_sfc_static_sub(ni,nj,lat,lon,topo)
c
c*****************************************************************************
c
c	program to build the necessary static files for the laps surface
c       analysis.  files created include:
c
c           removed jul 95 --> 'fnorm2.lut' (the barnes look-up table), 
c           removed jan 97 --> 'bkgwts.wts' (background weights), 
c                 'drag_coef.dat' (surface roughness), 
c             and 'pbl_top.dat' (boundary layer depth).
c
c	changes:
c 	p.a. stamus	04-12-94  original from the individual programs.
c                       09-01-94  porting changes.
c                       07-26-95  remove fnorm2 calcs.
c                       01-08-97  porting changes. remove 'stations.in'
c                                  requirement, general clean up.
c                       04-09-97  remove equivalences.
c
c****************************************************************************
c
        implicit none
ccccc	include 'laps_sfc.inc'
c
c.....	grids for the outputs, weights, and stuff 
c
        integer ni,nj
	real akk(ni,nj)
	real top(ni,nj), d1(ni,nj), d2(ni,nj)
        real fnorm(0:ni-1,0:nj-1)
c
        character*80 grid_fnam_common
        common/grid_fnam_cmn/ grid_fnam_common
c
c..... laps lat/lon grids.
c
	real lat(ni,nj), lon(ni,nj), topo(ni,nj)
	integer grid_spacing, len
	character dir_s*150,ext_s*31,units*10,comment*125,var_s*3
        integer imx,jmx, imn, jmn, icnt, npass, i_tmn, j_tmn
        integer i_tmx, j_tmx, n_obs_var, i, j, imaxm1, jmaxm1
        integer imax, jmax, istatus
        real slope1, slope2, slp, std, ave, smin, smax, con, scale
        real co_max_scl,co_topo_range , termx, x, pi, fno, smsng, rom2
        real topo_range, d, pbldz, top_adj, topmx, topo_mx, topo_mn
c
c
c.....  start:
c.....	get the laps lat/lon and topo data here so we can pass them to the 
c.....	routines that need them.
c

!        call get_laps_config(laps_domain,istatus)
!	dir_s = '../static/' 
!	dir_s = './' 
        call get_directory('static',dir_s,len)
        ext_s = 'nest7grid' 
c
	imax = ni
	jmax = nj
	imaxm1 = imax - 1
	jmaxm1 = jmax - 1
c
	topo_mx = -1.e30
	topo_mn =  1.e30
	do j=1,jmax
	do i=1,imax
	   if(topo(i,j) .gt. topo_mx) then
	      topo_mx = topo(i,j)
	      i_tmx = i
	      j_tmx = j
	   endif
	   if(topo(i,j) .lt. topo_mn) then
	      topo_mn = topo(i,j)
	      i_tmn = i
	      j_tmn = j
	   endif
	enddo !i
	enddo !j
	topo_range = topo_mx - topo_mn
c
	print *,' '
	print *,' topography: '
	write(6,211) topo_mx, i_tmx, j_tmx
 211	format(1x,'max of ',f10.1,' at ',i3,',',i3)
	write(6,212) topo_mn, i_tmn, j_tmn
 212	format(1x,'min of ',f10.1,' at ',i3,',',i3)
	write(6,213) topo_range
 213	format(1x,'range of ',f10.1)
	print *,' '
c
	pi = 4. * atan(1.0) 
	x = 1.e-3
	fno = 1. / (sqrt(2. * pi)) 
c
c..........................................................................
c
c	this section calculates the drag coefficient for the laps domain.
c
c	original by john mcginley.  date unknown.
c	changes:
c		p. stamus	11-12-19  changes for new laps grids.
c				11-20-91  changed scaling to reduce akk.
c                               03-17-95  auto scaling code for other domains.
c
c..........................................................................
c
	co_max_scl = 51.
	co_topo_range = 3099.8
	call constant(akk,1.,imax,jmax)
c
	smax = -1.e30
	smin = 1.e30
	icnt = 0
	std = 0.
	ave = 0.
	do j=2,jmaxm1
	do i=2,imaxm1
           slope1 = topo(i,j+1) - topo(i,j-1)
           slope2 = topo(i+1,j) - topo(i-1,j)
           if(slope1 .eq. 0.) slope1 = 0.001
           if(slope2 .eq. 0.) slope2 = 0.001
           slp = sqrt((slope1)**2 + (slope2)**2)
           icnt = icnt + 1
           akk(i,j) = slp
           ave = ave + slp
           std = std + slp ** 2
           if(slp .gt. smax) then
	      smax = slp
	      imx = i
	      jmx = j
	   endif
           if(slp .lt. smin) then
	      smin = slp
	      imn = i
	      jmn = j
	   endif
	enddo !i
	enddo !j
c
	ave = ave / float(icnt)
	std = sqrt(std / float(icnt))
c
c.....	we want this factor to serve as a  multiplier on the drag coef 
c.....  constant in lapsanal.  the scale factor was arbitrarilly set so
c.....  'akk' ranged between 51 and 1.1.  here, the range of the topography
c.....  in the domain is scaled against the colorado version, so the scaling
c.....  will adjust to smoother (or rougher) domains.
c
	con = co_max_scl * (topo_range / co_topo_range)
	scale = (con / smax) * (smin + 1.)
c
	do j=2,jmaxm1
	do i=2,imaxm1
	   akk(i,j) = scale * akk(i,j) / (smin + 1.)  !akk goes from 51 to 1.1
	   if(akk(i,j) .lt. 1.) akk(i,j) = 1.
	enddo !i
	enddo !j
c
         
	open(19,file=dir_s(1:len)//'drag_coef.dat',
     &                    form='unformatted',status='unknown')
	print *,' drag coef.: '
	write(6,1000)ave,std,smax,imx,jmx,smin,imn,jmn,scale
1000	format(1x,'ave: ',f12.4,' std dev: ',f12.4,/,' max: ',f12.4,
     &         ' at ',2i3,/,' min: ',f12.4,' at ',2i3,/,
     &         ' calc scale factor: ',f12.4)
	write(19) akk
        close(19)
	print *,' '
	print *,' normal completion of drag_coef.'



        goto 9999 ! bypass pbl top code

c
c
c...........................................................................
c
c	this section calculates the top of the pbl across the laps domain.
c
c	original by john mcginley (somewhere in the depths of time).
c	changes:
c		p. stamus	11-12-91  changes for new laps grids.
c				11-15-91  incr plbtop to prevent hi wnds
c
c...........................................................................
c
	termx = 0.
	do j=1,nj
	do i=1,ni
	  if(topo(i,j) .gt. termx) termx = topo(i,j)
        enddo !i
        enddo !j
c
        npass = 1
	smsng = 1.e37
c
	rom2 = 0.01
	call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	call barnes2(top,imax,jmax,topo,smsng,npass,d1,d2,fnorm)
c
	topmx = 0.
	do j=1,nj
	do i=1,ni
	  if(top(i,j) .gt. topmx) topmx = top(i,j)
	enddo !i
	enddo !j
c
        write(6,302) termx, topmx
 302   format(1x,' max. terrain ',f8.0,'; max. smooth terrain ',f8.0)
c
!	pbldz=475.	! pbl depth...old value 1000 m
!	pbldz=775.	! pbl depth...old value 1000 m
!	pbldz=925.	! pbl depth...old value 1000 m
	pbldz = 1225.
c
	top_adj = termx - topmx + pbldz
	do j=1,nj
	do i=1,ni
	  top(i,j) = top(i,j) + top_adj
	  if((top(i,j) - topo(i,j)) .lt. 100.) then
	    write (6,300) top(i,j),topo(i,j),i,j
	    top(i,j) = topo(i,j) + 100.
	  endif
	enddo !i
	enddo !j
 300    format(1x,'top pbl ',f8.0,' top terrain ',f8.0,' at ',2i6)
c
	open(10,file=dir_s(1:len)//'/pbl_top.dat',
     &          form='unformatted',status='unknown')
	write(10) top
	close(10)
	print *,' '
	print *,' normal completion of pbl_top'
c
 9999   print *,' '
	print *,' normal completion of build_sfc_static'
c
	return
	end
c
c
	subroutine barnes2(t,imax,jmax,to,smsng,npass,h1,h2, fnorm)
c
ccc	include '../../source/sfc/laps_sfc.inc'
        implicit none
        integer imax, jmax
	real to(imax,jmax), t(imax,jmax), val(imax*jmax)
	real h1(imax,jmax), h2(imax,jmax)
	integer iob(imax*jmax), job(imax*jmax), dx, dy
        real badd
	parameter(badd = 1.e6 - 2.)	! bad data value
        real sum, sumwt, sum2, sumwt2, smsng
        integer i,j,n, ipass, ncnt, npass 
	real fnorm(0:imax-1,0:jmax-1)

	call zero(h1,imax,jmax)
	call zero(h2,imax,jmax)
c
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
          endif
        enddo !i
        enddo !j
	if(ncnt .eq. 0) then
	  print *,' *** ncnt = 0 in barnes2. ***'
	  return
	endif
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
cc	      sum = fnorm2(dx,dy,kdim) * val(n) + sum
	      sum2 = fnorm(dx,dy) * val(n) + sum2
cc	      sumwt = sumwt + fnorm2(dx,dy,kdim)
	      sumwt2 = sumwt2 + fnorm(dx,dy)
	    enddo !n
c
	    if(sumwt2 .eq. 0.) then
	      if(ipass .eq. 1) then
cc	        if(kdim .eq. 1) then
cc	          print *,
cc     &              ' *** warning. sumwt = 0, kdim = 1 in barnes2 ***'
cc	          go to 500
cc	        else
		   print *,' got into wierd loop.............'
cc	          do kk=kdim,1,-2
	            sum2 = 0.
	            sumwt2 = 0.
	            do n=1,ncnt
	              dx = abs(i - iob(n))
	              dy = abs(j - job(n))
cc	              sum = fnorm2(dx,dy,kk) * val(n) + sum
	              sum2 = (fnorm(dx,dy) + .01) * val(n) + sum2
	              sumwt2 = sumwt2 + (fnorm(dx,dy) + .01)
	            enddo !n
	            if(sumwt2 .ne. 0.) go to 490
cc	          enddo !kk
cc	        endif
	      else
	        go to 500
	      endif
	    endif
c
490	    continue
cc	    t(i,j) = sum / sumwt
	    t(i,j) = sum2 / sumwt2
c
500	  continue
          enddo !i
	  enddo !j
c
	  if(ipass .eq. 2) then
	    call move(h2,to,imax,jmax)
	    call diff(h1,t,t,imax,jmax)
!	    write(9,915)
!915	    format(' after 2nd barnes pass')
!	    write(9,909) rom2
!909	    format(' radm2 = ',f8.4)
	    go to 550
	  endif
!	  write(9,912)
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
	subroutine dynamic_wts(imax,jmax,n_obs_var,rom2,d, fnorm)
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
c     changes:
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
	implicit none
        integer imax,jmax 
        real fnorm(0:imax-1,0:jmax-1)
        real pi, con, area, fno, rr, d, rom2
        integer n_obs_var,  dx, dy
          
ccc	include '../../source/sfc/laps_sfc.inc'
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
	   rr = dx*dx + dy*dy
	   fnorm(dx,dy) = fno * (exp( -(rr * rom2)))
	enddo !dx
	enddo !dy
c
c.... that's it.
c
	return
	end
