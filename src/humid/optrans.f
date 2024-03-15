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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
               subroutine optrans(
     &		nlev,
     &		p,
     &		t,
     &		q,
     &		pw,
     &		o3,
     &		angle,
     &		nchannels,
     &		channels,
     &		ptrans
     &                     )

*	version 4.0  25 oct 96
*	
*	name:	subroutine optrans
*
*
*	purpose: computes the transmittance using the optical path
*		 regression algorithm.
*	
*	input:  
*		i*4 nlev      number of levels
*		r*4 p(nlev)   pressure levels
*		r*4 t(nlev)   temperature (k)
*		r*4 q(nlev)   mixing ratio 
*		r*4 pw(nlev)  precipitable water (cm) 
*		r*4 o3(nlev)  integrated ozone profile
*		r*4 angle   zenith angle to use (degrees)
*		i*4 nchannels number of channels to use
*		i*4 channels(nchannels) array of channels for which to 
*		    compute transmittance.
*
*		from for042 wet coefficients
*		from for043 dry coefficients
*		from for044 ozo coefficients
*		from for061 wet control file
*		from for062 dry control file
*		from for063 ozo control file
*
*
*	output: r*4 ptrans(0:nlev,nchannels) transmittances at nlev levels
*
*	author:	thomas j. kleespies
*		physics branch
*		satellite research laboratory
*		office of research and applications
*		noaa/nesdis
*		301-763-8136
*		301-763-8108 fax
*
*	mailing address:
*		810 nsc e/ra-14
*		noaa/nesdis
*		washington, d.c. 20233
*
*	email:
*		kleespies@nzms.wwb.noaa.gov
*
*

	implicit none
        save

	include 'constants_optran.inc'

*	input
	integer nlev  	! number of levels
	real p(nlev)  	! pressure (mb)
	real t(nlev)	! temperature (k)
	real q(nlev)	! mixing ratio
	real pw(nlev) 	! precipitable water profile (cm)
	real o3(nlev)   ! integrated ozone profile
	real angle	! zenith angle (degrees)
	integer nchannels ! number of channels to process
	integer channels(nchannels) ! channel list

*	output
	! estimated transmittance in pressure space
	real ptrans(0:nlevel,nchannels)

*	local
		
	real*8 bw(nwet+1, nw ,nchan)  ! vector of wet regression coefficients
	real*8 bd(ndry+1, nw ,nchan)  ! vector of dry regression coefficients
	real*8 bo(nozo+1, nw ,nchan)  ! vector of dry regression coefficients

	real ww(0:nw)		! wet absorber amount space
	real aa(0:nw)		! dry absorber amount space
	real oo(0:nw)		! ozo absorber amount space
	real pwmax 
        data pwmax /13.87981/   ! max wet absorber at ww(nw)
	real pmax  
        data pmax /2100./       ! max dry absorber at aa(nw)
	real ozmax  
        data ozmax /1.1403912/  ! max ozo absorber at oo(nw)
c        real pwmax /12.33761/   ! maximum value for wet absorber space
c	real pmax  /2000./	! max dry absorber at aa(nw)
c        real ozmax /153.4007/   ! maximum value for ozo absorber space
	real slantpw(nlevel)	! wet absorber on a generalized slant path
	real slantp(nlevel)	! dry absorber on a generalized slant path
	real slantoz(nlevel)	! ozo absorber on a generalized slant path

	real odwet(nlevel,nchan)    ! wet optical depth on pressure levels
	real oddry(nlevel,nchan)    ! dry optical depth on pressure levels
	real odozo(nlevel,nchan)    ! ozo optical depth on pressure levels

	real dtr 
        data dtr /1.74532952e-2/ ! degrees to radians conversion

	real secant		! secant of angle


	integer level

	integer i ! local index
	integer ichan 	! channel index for read
	integer ipred	! predictor index for read

	integer pred_index_wet(0:nwet,nchan) ! wet predictors
	integer pred_index_dry(0:ndry,nchan) ! dry predictors
	integer pred_index_ozo(0:nozo,nchan) ! dry predictors

c	integer ios ! io return from open


*  this may need a save statement on some systems
	logical first 
        data first /.true./     ! flag to initialize coefficients and such
				! first time routine is called.
c     laps apps
        character*200 fname
        integer len

	if(first) then ! first time called, initialize variables and parameters


           call get_directory ('static',fname,len)

	
	open(42,file=fname(1:len)//'optranlib/wet_coeff.dat',
     &          form='unformatted')
	open(43,file=fname(1:len)//'optranlib/dry_coeff.dat',
     &          form='unformatted')
	open(44,file=fname(1:len)//'optranlib/ozo_coeff.dat',
     &          form='unformatted')
	

	 read(42) bw	! read wet coefficients
	 close(42)
	 read(43) bd	! read dry coefficients
	 close(43)
	 read(44) bo	! read ozo coefficients
	 close(44)

*	read the control files

 3100     format(14i3)

         open(61,
     &		file=fname(1:len)//'optranlib/wet_control_file.dat',
     &          form='formatted',
     &          status='old',
     &          err=310)

         do ichan = 1 , nchan
          read(61,3100)  (pred_index_wet(ipred,ichan), ipred=0,nwet)
         enddo

         close(61)

         open(62,
     &		file=fname(1:len)//'optranlib/dry_control_file.dat',
     &          form='formatted',
     &          status='old',
     &          err=310)

         do ichan = 1 , nchan
          read(62,3100)  (pred_index_dry(ipred,ichan), ipred=0,ndry)
         enddo

         close(62)

         open(63,
     &		file=fname(1:len)//'optranlib/ozo_control_file.dat',
     &          form='formatted',
     &          status='old',
     &          err=310)

         do ichan = 1 , nchan
          read(63,3100)  (pred_index_ozo(ipred,ichan), ipred=0,nozo)
         enddo

         close(63)

*	construct the absorber amount spaces for wet and dry cases

	 call make_wlevels(1, pwmax, 1.0e-7, ww)
         call make_wlevels(1, pmax ,   0.01, aa)
         call make_wlevels(1, ozmax ,2.0e-5, oo)

	 first = .false.  ! assure that we don't come here again

	endif ! first


*	this is where the actual work starts

*	it may be desirable to put in a check to assure that nlev <= nlevel

	secant = 1.0 / cos(angle*dtr)
	do level= 1 , nlev	! compute slant path absorber amounts
         slantp(level) = p(level)*secant - 0.005
	 slantpw(level) = pw(level)*secant
	 slantoz(level) = o3(level)*secant
	enddo

*	get wet optical depth profile
           call optrans_species     (
     &		nlev,
     &		nwet+1,
     &		nwet,
     &          bw,
     &          ww,
     &          p,
     &          t,
     &		q,
     &          slantpw,
     &          nchannels,
     &          channels,
     &		pred_index_wet,
     &          odwet
     &                          )
*	get dry optical depth profile
           call optrans_species     (
     &		nlev,
     &		ndry+1,
     &		ndry,
     &          bd,
     &          aa,
     &          p,
     &          t,
     &		q,
     &          slantp,
     &          nchannels,
     &          channels,
     &		pred_index_dry,
     &          oddry
     &                          )
*	get ozone optical depth profile
           call optrans_species     (
     &		nlev,
     &		nozo+1,
     &		nozo,
     &          bo,
     &          oo,
     &          p,
     &          t,
     &		q,
     &          slantoz,
     &          nchannels,
     &          channels,
     &		pred_index_ozo,
     &          odozo
     &                          )

*	compute total transmittance profile

	do i = 1 ,  nchannels ! now process only requested channels
	 ichan = channels(i)
	 ptrans(0,i) = 1.0
	 do level = 1 , nlev
	  ptrans(level,i) = ptrans(level-1,i)
     &		    *exp(-(
     &			   odwet(level,i)
     &			  +oddry(level,i)
     &			  +odozo(level,i)
     &			) )
	 enddo ! nlev

	enddo  ! channel

*	that's all folks!

	return

  310   continue
	write(6,*) 'error on control file open '
c	write(6,*) ios
	stop
	end
	subroutine optrans_species	(
     &		nlev,
     &		ncoeff,
     &		npred,
     &		b,
     &		ww,
     &		p,
     &		t,
     &		q,
     &		absorber,
     &		nchannels,
     &		channels,
     &		pred_index,
     &		od
     &					)

*
*	version 4: 25 october 1995
*	
*	name: optrans_species
*
*
*	purpose:
*	this is the routine that determines the optical depth profile
*	given the pressure profile, temperature profile, and absorber
*	profile.  this is generalized such that it will work with any
*	species.
*	
*	input:
*	integer nlev		! number of levels in input atmosphere
*	integer ncoeff		! number of coefficients
*	integer npred		! number of predictors
*	real b(ncoeff,nw,nchan)	! regression coefficients for this species
*	real ww(0:nw)		! absorber space for this species
*	real p(nlev)		! atmospheric pressure profile
*	real t(nlev)		! atmospheric temperature profile
*	real q(nlev)		! atmospheric mixing ratio profile
*	real absorber(nlev)	! slant path absorber profile
*	integer nchannels	! number of channels 
*	integer channels(nchannels) ! vector of channels to use
*       integer pred_index(0:npred,nchan)  ! index set of predictors to use
*					! for each channel.
*
*	output variables
*	real od(nlev,nchannels) ! optical depth in pressure space
*
*	author:	thomas j. kleespies
*		physics branch
*		satellite research laboratory
*		office of research and applications
*		noaa/nesdis
*		301-763-8136
*		301-763-8108 fax
*
*	mailing address:
*		810 nsc e/ra-14
*		noaa/nesdis
*		washington, d.c. 20233
*
*	email:
*		kleespies@nzms.wwb.noaa.gov
*
	implicit none

	include 'constants_optran.inc'

*	input variables
	integer nlev		! number of levels in input atmosphere
	integer ncoeff		! number of coefficients
	integer npred		! number of predictors
	real*8 b(ncoeff,nw,nchan)! regression coefficients
	real ww(0:nw)		! absorber space
	real p(nlev)		! atmospheric pressure profile
	real t(nlev)		! atmospheric temperature profile
	real q(nlev)		! atmospheric mixing ratio profile
	real absorber(nlev)	! atmospheric absorber profile
	integer nchannels	! number of channels 
	integer channels(nchannels) ! vector of channels to use
        integer pred_index(0:npred,nchan) ! index set of predictors to use
     &                                    ! for each channel.

*	output variables
	real od(nlevel,nchannels)	! optical depth

*	local variables

	real*8 sum		! summation 

	real average_absorber	! average absorber amount

        real kw(nw)             ! k = predictand in absorber space
	
	real kp(nlevel)		! absorption coeff in pressure space

	real dp(nlevel)		! absorber within layer
	
	real factor(nlevel)	! linear interpolation factor

	integer i,level,kt	! utility variables
	integer ichan,jchan
	integer k2,k1
	integer ipred

	integer maxwlevels	! actual number of levels used in w space

	integer wlevels(nw)	! index array pointing to the w space
				! levels actually used

        integer m		! index for interpolation

        logical search		! .true. for first call to search_plevel_linear

        real*8 xx(nw,15) ! presently hard coded for 15 potential predictors 

*	find values in w space that bracket the absorber levels
	search = .true.

	do level = 1 , nlev
	 kt = 2*(level-1) + 1
	 
	 if(level .eq. 1) then
	  average_absorber = absorber(level) / 2.
	 else
	  average_absorber = (absorber(level)+absorber(level-1))/2.
	 endif

	 call search_plevel_linear(search,ww,average_absorber,nw,m)
	 wlevels(kt) = m
	 wlevels(kt+1) = m+1
	 factor(level) = (ww(m+1)-average_absorber)/(ww(m+1)-ww(m))
	 if(level .eq.1) then
	  dp(level) = absorber(level)
	 else
	  dp(level) = (absorber(level)-absorber(level-1))
	 endif
	enddo

	maxwlevels = kt + 1

*	compute the predictors for this atmosphere

        call get_predictors_all(
     &          p,
     &          absorber,
     &          t,
     &		q,
     &		nlev,
     &          ww,
     &          maxwlevels,
     &		wlevels,
     &		xx
     &                            )

	do  300 jchan = 1 , nchannels

	 ichan = channels(jchan)

	 if(pred_index(0,ichan) .lt. 1) then
	  do level = 1 , nlev ! no absorption for this channel
	   od(level,jchan) = 0.0
	  enddo
	  go to 300
	 endif

*	zero kw array

	 do i = 1 , nw
	  kw(i) = 0.0
	 enddo

*	compute k in w space

	 do 200 i = 1 , maxwlevels 

	  k1 = wlevels(i)
		
	  if(kw(k1) .ne. 0.0) go to 200 ! already calculated

           sum = b(1,k1,ichan)

           do ipred = 2 , ncoeff
            sum = sum + b(ipred,k1,ichan)
     &                      * xx(k1,pred_index(ipred-1,ichan))
           enddo

	   kw(k1) = sum

  200    continue

*       interpolate back to p space

	 do level = 1 , nlev
	  kt = 2*(level-1) + 1
	  k1 = wlevels(kt)
	  k2 = wlevels(kt+1)
	  kp(level) = kw(k2) - (kw(k2)-kw(k1))*factor(level)
	  if(kp(level) .lt. 0.0) kp(level) = 0.0
	 enddo

*	compute optical depth

         od(1,jchan) = 0.0

         do level = 1 , nlev

	  od(level,jchan) = dp(level)*kp(level) 

	 enddo

  300   continue ! jchan

	return
	end
	subroutine get_predictors_all(
     &		press,
     &		absorb,
     &		temp,
     &		mixr,
     &		nlev,
     &		ww,
     &		maxwlevels,
     &		wlevels,
     &		x )
*	
*	name: subroutine get_predictors_all
*
*	version 4: 25 october 95
*
*	purpose:
*	gets all predictors out of set of 14
*	
*	input:
*	press(nlev)	atmospheric pressure levels
*	absorb(nlev)	atmospheric absorber profile on slant path
*	temp(nlev)	atmospheric temperature profile
*	mixr(nlev)	atmospheric mixing ratio profile
*	nlev		number of atmospheric levels
*	ww(0:nw)	standard absorber space
*	maxwlevels	number of levels to be filled in absorber space
*	wlevels(maxwlevels)index of absorber space levels to be interpolated to 
*
*	output:
*
*	x - various predictors
*
*	note in the following description 'x' is multiplication,
*					  '^' is exponention
*					  '*' is a "star'd" quantity
*	x(1)	t	
*	x(2)	p	
*	x(3)	t^2
*	x(4)	p^2
*	x(5)    txp
*	x(6)	t^2xp
*	x(7)	txp^2
*	x(8)	t^2xp^2
*	x(9)	t*
*	x(10)	p*
*	x(11)	t**
*	x(12)	p**
*	x(13)	t***
*	x(14)	p***	
*	x(15)   q
*
*	author:	thomas j. kleespies
*		physics branch
*		satellite research laboratory
*		office of research and applications
*		noaa/nesdis
*               301-763-8136
*               301-763-8108 fax
*
*	mailing address:
*		810 nsc e/ra-14
*		noaa/nesdis
*		washington, d.c. 20233
*
*	email:
*		kleespies@nzms.wwb.noaa.gov
*
*

	implicit none

	include 'constants_optran.inc'

*	input
	integer nlev		! number of atmospheric levels
	real press(nlev)	! atmospheric pressure levels
	real absorb(nlev)		! atmospheric absorber profile
	real temp(nlev)		! atmospheric temperature profile
	real mixr(nlev)		! atmospheric mixing ratio profile
	real ww(0:nw)		! standard absorber profile

	integer maxwlevels	! number of levels to be filled in absorber space
	integer wlevels(maxwlevels) ! index of absorber space levels to be used

*	output

	real*8 x(nw,15)

	real*8 dela

	real*8 s1,s2,s3,s4,s5,s6

*	local

	integer i,j,k,kt,level

	logical search_absorber
	integer m1

	real*8 linear_interpolation_d ! function to do dp linear interpolation

	real*8 xx1(nlevel)	! these are the star'd quantities
	real*8 xx2(nlevel)	! in pressure space
	real*8 xx3(nlevel)
	real*8 xx4(nlevel)
	real*8 xx5(nlevel)
	real*8 xx6(nlevel)

	integer maxlevel
	parameter (maxlevel=100)
	real*8 p(maxlevel)
	real*8 t(maxlevel)
	real*8 absorber(maxlevel)
	real*8 q(maxlevel)

	search_absorber = .true.

*	zero arrays

	do i = 1 , nw
	 do j = 1 , 14
	   x(i,j) = 0.0d0
	 enddo
	enddo	 

*	make atmospheric variables double precision

	if(nlev .gt. maxlevel) then
	 write(6,*) 'too many levels in input profile'
	 write(6,*) maxlevel , ' is maximum allowed'
	 write(6,*) 'change parameter maxwlevels in get_predictors_all'
	 return
	endif

	do level = 1 , nlev
	 t(level) = temp(level)
	 p(level) = press(level)
	 absorber(level) = absorb(level)
	 q(level) = mixr(level)
	enddo

*	first compute star'd terms for the specified atmosphere
*	i.e. on pressure levels

	if(absorber(1) .gt. 0.0) then
	 xx1(1)  = t(1) / 2. ! if we treat t(0) = 0, absorber(0) = 0
	 xx2(1)  = p(1) / 2. ! if we treat p(0) = 0, absorber(0) = 0
	 xx3(1) = t(1) / 2. ! if we treat t(0) = 0, absorber(0) = 0
	 xx4(1) = p(1) / 2. ! if we treat p(0) = 0, absorber(0) = 0
	 xx5(1) = t(1) / 2. ! if we treat t(0) = 0, absorber(0) = 0
	 xx6(1) = p(1) / 2. ! if we treat p(0) = 0, absorber(0) = 0
 	else
	 xx1(1) = 0.0
	 xx2(1) = 0.0
	 xx3(1) = 0.0
	 xx4(1) = 0.0
	 xx5(1) = 0.0
	 xx6(1) = 0.0
	endif

	s1 = 0.0
	s2 = 0.0
	s3 = 0.0
	s4 = 0.0
	s5 = 0.0
	s6 = 0.0

	do k = 2 , nlev
	 dela = absorber(k) - absorber(k-1)
	 s1 = s1 + (t(k)+t(k-1))*dela	! t*
	 s2 = s2 + (p(k)+p(k-1))*dela	! p*
	 s3 = s3 + (absorber(k)*t(k) + absorber(k-1)*t(k-1))*dela !t**
	 s4 = s4 + (absorber(k)*p(k) + absorber(k-1)*p(k-1))*dela !p**
	 s5 = s5 + (absorber(k)*absorber(k)*t(k) + 
     &		    absorber(k-1)*absorber(k-1)*t(k-1))*dela !t***
	 s6 = s6 + (absorber(k)*absorber(k)*p(k) + 
     &		    absorber(k-1)*absorber(k-1)*p(k-1))*dela !p***

	 if(s1 .eq. 0.0) then 
	   xx1(k) = 0.0
	 else
	   xx1(k) = s1/(2.0*absorber(k))
	 endif
	 if(s2 .eq. 0.0) then 
	   xx2(k) = 0.0
	 else
	   xx2(k) = s2/(2.0*absorber(k))
	 endif
	 if(s3 .eq. 0.0) then 
	   xx3(k) = 0.0
	 else
	   xx3(k) = s3/(absorber(k)*absorber(k))
	 endif
	 if(s4 .eq. 0.0) then 
	   xx4(k) = 0.0
	 else
	   xx4(k) = s4/(absorber(k)*absorber(k))
	 endif
	 if(s5 .eq. 0.0) then 
	   xx5(k) = 0.0
	 else
	   xx5(k) = 1.5*s5/(absorber(k)*absorber(k)*absorber(k))
	 endif
	 if(s6 .eq. 0.0) then 
	   xx6(k) = 0.0
	 else
	   xx6(k) = 1.5*s6/(absorber(k)*absorber(k)*absorber(k))
	 endif

	enddo

*	next compute for the standard absorber levels
	 
	do 200 k = 1 , maxwlevels	

	  kt = wlevels(k)
	  if(x(kt,1) .ne. 0) go to 200 ! already did this level


          call search_wlevel_linear(
     &		search_absorber,absorb,nlev,ww(kt),m1)

          x(kt,1)    = linear_interpolation_d(absorber(m1),t(m1), 	! t
     &                                  absorber(m1+1),t(m1+1),
     &                                  dble(ww(kt)))

          x(kt,2)    = linear_interpolation_d(absorber(m1),p(m1), 	! p
     &                                  absorber(m1+1),p(m1+1),
     &                                  dble(ww(kt)))

	  x(kt,3)    =  x(kt,1)*x(kt,1) 			! t**2
	  x(kt,4)    =  x(kt,2)*x(kt,2) 			! p**2
	  x(kt,5)    =  x(kt,1)*x(kt,2) 			! t*p
	  x(kt,6)    =  x(kt,3)*x(kt,2) 			! t**2 p
	  x(kt,7)    =  x(kt,1)*x(kt,4) 			! t p**2
	  x(kt,8)    =  x(kt,3)*x(kt,4) 			! t**2 p**2

          x(kt,9)= 	linear_interpolation_d(absorber(m1),xx1(m1),
     &                                  absorber(m1+1),xx1(m1+1),
     &                                  dble(ww(kt)))

          x(kt,10)= 	linear_interpolation_d(absorber(m1),xx2(m1),
     &                                  absorber(m1+1),xx2(m1+1),
     &                                  dble(ww(kt)))

          x(kt,11)= 	linear_interpolation_d(absorber(m1),xx3(m1),
     &                                  absorber(m1+1),xx3(m1+1),
     &                                  dble(ww(kt)))

          x(kt,12)= 	linear_interpolation_d(absorber(m1),xx4(m1),
     &                                  absorber(m1+1),xx4(m1+1),
     &                                  dble(ww(kt)))

          x(kt,13)= 	linear_interpolation_d(absorber(m1),xx5(m1),
     &                                  absorber(m1+1),xx5(m1+1),
     &                                  dble(ww(kt)))

          x(kt,14)= 	linear_interpolation_d(absorber(m1),xx6(m1),
     &                                  absorber(m1+1),xx6(m1+1),
     &                                  dble(ww(kt)))

          x(kt,15)= 	linear_interpolation_d(absorber(m1),q(m1),
     &                                  absorber(m1+1),q(m1+1),
     &                                  dble(ww(kt)))

  200   continue	! maxwlevels loop

	return
	end
	real*8 function linear_interpolation_d(x1,y1,x2,y2,x)
*
*	performs simple linear interpolation
*	
******* warning: don't forget to declare real in calling routine **********
*
*	input-
*
*	x1,y1,x2,y2 - pairs bounding the interval
*	x - independent variable - we want to find the y value associated to x
*
*	returned value: dependent variable associated with x
*
*       3 dec 90
*
*
*	author:	thomas j. kleespies
*		physics branch
*		satellite research laboratory
*		office of research and applications
*		noaa/nesdis
*               301-763-8136
*               301-763-8108 fax
*
*	mailing address:
*		810 nsc e/ra-14
*		noaa/nesdis
*		washington, d.c. 20233
*
*	email:
*		kleespies@nzms.wwb.noaa.gov

	implicit none

	real*8 x1,y1,x2,y2,x
	real*8 denominator

	denominator = x2-x1

	if(denominator .ne. 0.0) then
	   linear_interpolation_d = y2 - (y2-y1)*(x2-x)/denominator
	else
	   linear_interpolation_d = y2
	endif

	return
	end
	subroutine search_wlevel_linear(first,absorber,nlev,ww,m)
*	
*	name:	subroutine search_wlevel_linear
*
*	version 4: 25 october 1995
*
*	purpose:
*	searches the atmospheric absorber profile for
*	the 'best' bracketing values for subsequent linear interpolation.
*	
*	input:
*	first		set .true. for first call for a given profile
*	absorber(nlev)  absorber amount profile
*	nlev		number of levels in absorber
*	ww		value of absorber to bracket
*
*	output:
*	m	 	returned value of m which represents
*			the first point in the interpolation array.
*
*	author:	thomas j. kleespies
*		physics branch
*		satellite research laboratory
*		office of research and applications
*		noaa/nesdis
*		301-763-8136
*		301-763-8108 fax
*
*	mailing address:
*		810 nsc e/ra-14
*		noaa/nesdis
*		washington, d.c. 20233
*
*	email:
*		kleespies@nzms.wwb.noaa.gov
*

	implicit none

	logical first
	integer nlev		! number of levels
	real absorber(nlev)	! absorber amount profile
	real ww			! value of absorber to bracket
	integer m		! returned value of m which represents
				! the first point in the interpolation array.

	integer i,j		! utility variables
	integer start_level	! starting level for search

	if(first) then
	  start_level = 1
	  first = .false.
	else
	  start_level = m  ! m should be retained by calling routine
	endif

	do i = start_level , nlev
	 if(ww .lt. absorber(i)) then
	  j = i
	  go to 110
	 endif
	enddo

	j = nlev

  110   continue

	if(j .le. 1) then 	! extrapolate above atmosphere
	  m = 1 
	else if (j .eq. nlev) then ! extrapolate below atmosphere
	  m = nlev-1
	else
	  m = j-1		! bracket within atmosphere
	endif

	return
	end
	subroutine search_plevel_linear(first,ww,absorber,wlevels,m)
*	
*	name: 	subroutine search_plevel_linear
*
*	version 4: 25 october 1995
*		
*
*	purpose:
*	searches the standard absorber profile for
*	the 'best' bracketing values for subsequent linear interpolation.
*	
*	
*	input:
*	first		set .true. for first call for a given profile
*	ww(0:nw)	standard absorber profile
*	absorber	absorber amount to bracket
*	wlevels		actual number of levels used in ww space
*
*	output:
*	m		returned value of m which represents
*			the first point in the interpolation array.
*
*	author:	thomas j. kleespies
*		physics branch
*		satellite research laboratory
*		office of research and applications
*		noaa/nesdis
*		301-763-8136
*		301-763-8108 fax
*
*	mailing address:
*		810 nsc e/ra-14
*		noaa/nesdis
*		washington, d.c. 20233
*
*	email:
*		kleespies@nzms.wwb.noaa.gov
*
*

	implicit none

	include 'constants_optran.inc'

	logical first 		! set .true. for initial invocation
	real ww(0:nw)		! standard absorber profile
	real absorber		! absorber amount to bracket
	integer wlevels		! actual number of levels used in ww space
	integer m		! returned value of m which represents
				! the first point in the interpolation array.

	integer i,j		! utility variables

	integer start_level	! starting level for search

	if(first) then
	  start_level = 0
	  first = .false.
	else
	  start_level = m  ! m should be retained by calling routine
	endif

	do i = start_level , wlevels
	 if(absorber .lt. ww(i)) then
	  j = i
	  go to 110
	 endif
	enddo

	j = wlevels

  110   continue

	if (j .gt. 1) then
	 m = j-1 
	else
	 m = 1
	endif

	return
	end
	subroutine make_wlevels(nprofiles,maxabsorber,minabsorber,w)
*	
*	name: subroutine make_wlevels
*
*	version 4: 25 october 1995
*
*	purpose:
*	makes the absorber space of nw levels,
*	identified as w(0:nw).  the most important value for
*	defining the parameter 'alpha' is w(200), which must be
*	at least the maximum possible absorber amount. the code 
*	is structured to produce different standard profiles for each channel.
*	
*	input:
*	nprofiles	number of profiles to use
*			this is the dimension of maxabsorber
*			and the second dimension of w.
*	maxabsorber(nprofiles)	! absorber amount at w(200)
*	minabsorber	! absorber amount at w(1), w(0) == 0.0
*
*	output:
*	w(0:nw,nprofiles) ! standard absorber profiles
*
*	author:	thomas j. kleespies
*		physics branch
*		satellite research laboratory
*		office of research and applications
*		noaa/nesdis
*		301-763-8136
*		301-763-8108 fax
*
*	mailing address:
*		810 nsc e/ra-14
*		noaa/nesdis
*		washington, d.c. 20233
*
*	email:
*		kleespies@nzms.wwb.noaa.gov
*
*

	implicit none
        save

	include 'constants_optran.inc'

	integer nprofiles	! number of profiles to use
				! this is the dimension of maxabsorber
				! and the second dimension of w.
	real maxabsorber(nprofiles)	! absorber amount at w(200)
	real minabsorber	! absorber amount at w(1), w(0) == 0.0

	real w(0:nw,nprofiles) ! standard absorber profiles

	real*8 alpha    ! exponential parameter

	integer k     ! utility variables

	real*8 tolerance 
        data tolerance /1.0e-10/ ! iteration-to-iteration convergence criterion

	integer maxiter 
        data maxiter /100/      ! arbitrary value to test non-convergence

	real*8 x1,x2,fx,fp ! newton's method variables
			   ! fx is the function
			   ! fp is the function primed (as in it's derivative)
			   ! x1 and x2 are previous and new approximations

	integer ichan,iter ! local variables

	do ichan = 1 , nprofiles  ! profile loop

	 w(0,ichan) = 0		  ! top of atmosphere always has zero absorber
	 w(1,ichan) = minabsorber
	 w(nw,ichan) = maxabsorber(ichan)
	
*	solve for alpha by newton's method

* these are somewhat arbitrary initial guesses
* the lower one is important for the wet atmospheres, and the
* upper one is important for the dry atmospheres.

	 x2 = 0.3	
	 x1 = 3000.
	 iter = 0
	 do while (abs(x2-x1) .ge. tolerance)

	  iter = iter + 1
	  if(iter .gt. maxiter) then
	   write(6,*) 'make_levels failed to converge '
	   stop
	  endif

	  x1 = x2
! 	the function
	  fx = (exp(nw*x1)-1.) - (w(nw,ichan)/w(1,ichan))*(exp(x1)-1.)
! 	function primed 
	  fp = nw*exp(nw*x1) - (w(nw,ichan)/w(1,ichan))* exp(x1)   
	  x2 = x1 - fx/fp

	 enddo

	 alpha = x2

*	now that we have alpha, construct the profile 

	 do k = 2 , nw
	  w(k,ichan) = w(1,ichan) * (exp(k*alpha)-1.)/(exp(alpha)-1.)
	 enddo

	enddo ! profile loop

	return
	end





	real function compbright(vnu,t,tau,tskin,n)
c	computes brightness temperature for a temperature and
c	transmittance profile.

*	input
*
*	r*4 vnu		wavenumber
*	r*4 t(n)	temperature profile
*	r*4 tau(n)	transmittance profile
*	r*4 tskin	skin temperature
*	i*4 n		number of levels

	implicit none
	
	real vnu
	integer n
	real t(n)
	real tau(n)
	real tskin

	real sum,b1,b2,bs,tau1,tau2
	integer i

	real c1 
        data c1 /1.1905e-5/
	real c2 
        data c2 /1.4385/

	real planck,bright,v,temp,radiance
	planck(v,temp) = (c1*v*v*v)/(exp(c2*v/temp) - 1)
        bright(v,radiance) = c2*v / log(c1*v*v*v/radiance + 1.0)

	sum=0.

	b1=planck(vnu,t(1))
	tau1=tau(1)

	do i=2,n
	 b2=planck(vnu,t(i))
	 tau2 = tau(i)
	 sum=sum+.5*(b1+b2)*(tau1-tau2)
	 b1=b2
	 tau1=tau2
	enddo

	bs=planck(vnu,tskin)
	sum=sum+bs*tau(n)

	compbright = 0.0
	if(sum.gt.0.) compbright=bright(vnu,sum)

	return
	end

