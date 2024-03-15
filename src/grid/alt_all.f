cccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine alt_10by10_all(im,jm,id2,testrad,
     &          glatr,glonr,path,out2d,out3d,ncat,
     &		slope_lt,slope_ln,stdev,categorical,wrdlen,
     &		ignoreval)

        parameter(r2d=57.29577951,dtr=0.017453293)

        integer gds(200),wrdlen
        character*2 uname
        character*8 filename
        character*(*) path
        character*256 fullname

        real out2d(im,jm),out3d(im,jm,ncat),stdev(im,jm)
	real slope_lt(im,jm),slope_ln(im,jm)
	logical pure_nn,categorical,topogen
        real glatr(im,jm),glonr(im,jm)
        integer*2, allocatable:: topoin(:,:)
        integer*2  topotmp(id2,id2), ignoreval
	integer*1 topotmp_1(id2,id2)
        real dum2d(im,jm)
	real, allocatable:: topolat(:,:),topolon(:,:)
	real, allocatable:: xpts(:,:),ypts(:,:),dum(:,:)



!!	id2=1200 is 30" topo, landuse, soil
!!	id2=10 is 1 degree soil temp


	if (id2 .eq. 1200) then
	iperdeg=id2/10
	ispan=10
	pure_nn=.false.
	elseif (id2 .eq. 10) then
	iperdeg=1
	ispan=10
	pure_nn=.true.
	else
	write(6,*) 'not appropriate...quitting'
	stop
	endif
!
	out2d=-999.
!
        apmx=maxval(glatr)
        apmn=minval(glatr)

	if (glatr(im/2,jm/2) .gt. 0) then
        aloneast=glonr(im,jm)
	else
	aloneast=glonr(im,1)
	endif

        alonwest1=glonr(1,1)
        alonwest2=glonr(1,jm)

	if (apmn .gt. 0) then
!nh
        if ( (alonwest1 .lt. 0. .and. alonwest2 .lt. 0.) .or. 
     &	     (alonwest1 .gt. 0. .and. alonwest2 .gt. 0.) .or. 
     &	     (alonwest1 .gt. 0. .and. alonwest2 .lt. 0.) ) then
          alonwest=amin1(alonwest1,alonwest2)
	elseif (alonwest1 .lt. 0 .and. alonwest2 .gt. 0) then
! nh western boundary straddles dl
	  alonwest=amax1(alonwest1,alonwest2)
	endif

	else
!sh
        if ( (alonwest1 .lt. 0 .and. alonwest2 .lt. 0) .or. 
     &	     (alonwest1 .gt. 0 .and. alonwest2 .gt. 0) .or. 
     &	     (alonwest1 .lt. 0. .and. alonwest2 .gt. 0.) ) then
          alonwest=amin1(alonwest1,alonwest2)
	elseif (alonwest1 .gt. 0 .and. alonwest2 .lt. 0) then
! nh western boundary straddles dl
	  alonwest=amax1(alonwest1,alonwest2)
	endif

	endif



        istrad=0

        if (alonwest .gt. 0 .and. aloneast .lt. 0) then
        write(6,*) 'straddling the dl!!!'
        istrad=1
        endif

        if (alonwest .lt. -180.) then
        write(6,*) 'straddling the dl!!!'
        istrad=1
        endif

        if (alonwest .lt. 0 .and. aloneast .gt. 0) then
        write(6,*) 'straddling the pm!!!'
        istrad=2
        endif

        rnlat=apmx+0.1
        slat=apmn-0.1
        west=almn-0.1
        west=alonwest-0.1
        east=aloneast+0.1

        itemp=int(west/ispan)
        if (west .lt. 0) then
        west=float(itemp-1)*ispan
        else
        west=float(itemp)*ispan
        endif

        itemp=int(east/ispan)
        if (east .gt. 0) then
        east=float(itemp+1)*ispan
        else
        east=float(itemp)*ispan
        endif

        itemp=int(rnlat/ispan)
        rnlat=float(itemp+1)*ispan
        itemp=int(slat/ispan)

	if (slat .gt. 0) then
	slat=float(itemp)*ispan
	else
	slat=float(itemp)*ispan-ispan
	endif


        write(6,*) 'limits on discretized grid'
        write(6,*) 'west= ', west
        write(6,*) 'east= ', east
        write(6,*) 'north= ', rnlat
        write(6,*) 'south= ', slat


        if (istrad .eq. 0 .or. istrad .eq. 2) then
        ilim=iperdeg*(east-west)
        jlim=iperdeg*(rnlat-slat)
        endif
        if (istrad .eq. 1) then
        ilim=iperdeg*((east+360.)-west)
        jlim=iperdeg*(rnlat-slat)
        endif

!
!	create 10 degree longitude strip
!
	ilim=id2

        allocate(topoin(ilim,jlim))
        allocate(topolat(ilim,jlim))
        allocate(topolon(ilim,jlim))
        allocate(xpts(ilim,jlim))
        allocate(ypts(ilim,jlim))
        allocate(dum(ilim,jlim))

        do j=1,jlim
        do i=1,ilim
        xpts(i,j)=i
        ypts(i,j)=j
        enddo
        enddo

        if (istrad .eq. 1) then
        easttmp=east+360
        else
        easttmp=east
        endif

         write(6,*) 'loop i from: ', int(west), int(easttmp)-ispan

!!!!!!!!
!!!!!!!!
        do iiii=int(west),int(easttmp)-ispan,ispan
!!!!!!!!
!!!!!!!!

        if (east .ne. easttmp) then
          if (iiii .ge. 180) then
           itmp=iiii-360
          else
           itmp=iiii
          endif
        else
          itmp=iiii
        endif

!!!!!!!!
!!!!!!!!
        do j=int(slat),int(rnlat)-ispan,ispan
!!!!!!!!
!!!!!!!!

        call gen_filename(j,itmp,filename)
        write(6,*) 'want to process ', filename


	call s_len(path,len)

        fullname=path(1:len)//filename(1:7)

	if (wrdlen .eq. 2) then

!endianissue
        open(unit=1,file=trim(fullname),status='old'
     .          ,form='unformatted',access='direct',recl=2*id2)
        do nrrd=1,id2
        read(1,rec=nrrd) (topotmp(ii,nrrd),ii=1,id2)

!!! add an #ifdef structure, swap if needed?

	if (nrrd .eq. id2/2) then
	write(6,*) 'read vals: ', (topotmp(ii,nrrd),ii=1,id2,id2/10)
	endif

        enddo
        close(1)

	elseif (wrdlen .eq. 1) then
!endianissue
        open(unit=1,file=trim(fullname),status='old'
     .          ,form='unformatted',access='direct',recl=1*id2)
        do nrrd=1,id2
        read(1,rec=nrrd) (topotmp_1(ii,nrrd),ii=1,id2)
        enddo
        close(1)
	
	endif

cc      put this tile into longitudinal strip

        do jjj=(j-int(slat))*(iperdeg)+1,(j-int(slat))*(iperdeg)+id2
        do iii=1,id2
        jloc=jjj-(j-int(slat))*(iperdeg)
		if (wrdlen .eq. 2) then
       			 topoin(iii,jjj)=topotmp(iii,id2-jloc+1)
		elseif (wrdlen .eq. 1) then
        		 topoin(iii,jjj)=topotmp_1(iii,id2-jloc+1)
		endif
        enddo
        enddo

        enddo ! enddo for j

        write(6,*) 'ilim, jlim: ', ilim, jlim


cc      come up with an estimate of how many input 30 s points should
cc      be included.  consider roughly enough to average over the size of
cc      the gridbox.

       limit=1

cc
cc      30" ~ 0.927 km
cc      10' ~ 18.5 km
cc

cc      limit represents an intentionally over-large search radius to be
cc      used below to specify which input points define the output topography

!        limit=int((testrad/2.)/0.927)+2

        resfac=(1./iperdeg)*(111.2)
        limit=int((testrad/2.)/resfac)+2

!        write(6,*) 'if possible, will search +/- ', limit

        gds=0

        gds(1)=0
        gds(2)=ilim
        gds(3)=jlim
        gds(4)=int(slat*1000)
        gds(5)=int(iiii*1000)
        gds(6)=128
        gds(7)=int((rnlat-1./(iperdeg))*1000)
        gds(8)=int((iiii+10-1./(iperdeg))*1000)
        gds(9)=int(1./(iperdeg)*1000)
        gds(10)=int(1./(iperdeg)*1000)

	write(6,*) 'gds= ', (gds(i),i=1,10)

	call gdswiz(gds,1,ilim*jlim,-9999.,xpts,ypts,
     &	topolon,topolat,nret,0,dum,dum)

	do j=1,jlim
	 do i=1,ilim
	  if (topolon(i,j) .gt. 180) then
	   topolon(i,j)=topolon(i,j)-360.
	  endif
         enddo
	enddo


!!!!
	if (categorical) then
	write(6,*) 'calling alt_categories ', pure_nn
		call alt_categories(im,jm,glatr,glonr,
     &		ilim,jlim,topoin,topolat,topolon,pure_nn,out3d,
     &		ncat,out2d,gds,limit,testrad)
	else
	write(6,*) 'calling alt_interp ', pure_nn

	if (path(len:len) .eq. 'u') then
		topogen=.true.
	else
		topogen=.false.
	endif
		call alt_interp(im,jm,glatr,glonr,
     &		ilim,jlim,topoin,topolat,topolon,pure_nn,out2d,slope_lt,
     &		slope_ln,stdev,gds,limit,testrad,ignoreval,topogen)
	endif

!	write(6,*) 'after this tile: '

!	do j=jm,1,-jm/20
!	write(6,633) (out2d(i,j),i=1,im,im/15)
!	enddo

  633	format(15(f5.0,1x))


!!!!!!!!
!!!!!!!!
	enddo ! enddo for iiii (looping over strips)
!!!!!!!!
!!!!!!!!

ccccccccccccc

	do j=1,jm
	do i=1,im
	if (out2d(i,j) .eq. -999.) then
	write(6,*) 'undefined at ', i,j
	write(6,*) 'glat,glon: ', glatr(i,j),glonr(i,j)

	rmax=-999.
	do jj=j-1,j+1
	do ii=i-1,i+1
		if (out2d(ii,jj) .gt. rmax) then
			rmax=out2d(ii,jj)
		endif
	enddo
	enddo

	write(6,*) 'redefined out2d to be ', rmax
	out2d(i,j)=rmax

	endif
	enddo
	enddo


        deallocate(topoin,topolat,topolon,xpts,ypts,dum)

        end subroutine alt_10by10_all

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine alt_categories(im,jm,glatr,glonr,
     &          ilim,jlim,datain,datalat,datalon,
     &		pure_nn,cat3d,
     &          ncat,cat2d,gds,limit,res)

	real glatr(im,jm),glonr(im,jm)
	real datalat(ilim,jlim),datalon(ilim,jlim)
	integer*2 datain(ilim,jlim)
	real cat2d(im,jm),cat3d(im,jm,ncat)

	integer gds(200), counter(ncat)
	logical pure_nn

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	write(6,*) 'search limit is : ', limit
	write(6,*) 'ncat= ', ncat
        do j=1,jm
        do i=1,im
                                                                                
c determine lat/lon of target (output) grid point
c
        glatd = glatr(i,j)
        glond = glonr(i,j)
                                                                                
        call ced_ij(glatd,glond,x,y,gds)
                                                                                
!! kludgy fix for dead space "between" strips (like -60.003 longitude)
         
        if (x .le. 1200.5) then
        ii=int(x+0.5)
        else
!       write(6,*) 'kludgy fix applied...'
        ii=int(x)
        endif
         
        ji=int(y+0.5)
         
!       ii,ji represents nearest neighbor from input array
         
        if (ii .ge. 1 .and. ii .le. ilim .and.
     +      ji .ge. 1 .and. ji .le. jlim ) then
         
	if (.not. pure_nn) then

        icount=0
        sum=0.
        counter=0
         
         
!       no strictly n.n. now.  test as many points within the limits
!       as possible.

	
        do jjj=ji-limit,ji+limit
           do iii=ii-limit,ii+limit
                                                                                
        if (jjj .lt. 1 .or. jjj .gt. jlim .or.
     &      iii .lt. 1 .or. iii .gt. ilim) goto 47
                                                                                
        call greatcir(glatd,glond,
     &          datalat(iii,jjj),datalon(iii,jjj),dist)
                                                                                

        if ( dist .le. res/2.) then
        counter(datain(iii,jjj))=counter(datain(iii,jjj))+1
        icount=icount+1
        endif
                                                                                
  47    continue
                                                                                
           enddo
        enddo
c
                                                                                
	if (icount .gt. 0) then

        imaxcount=-9
        do n=1,ncat
                                                                                
        cat3d(i,j,n)=float(counter(n))/icount
                                                                                
        if (counter(n) .gt. imaxcount) then
                imaxcount=counter(n)
                ibest=n
        endif
                                                                                
        enddo
                                                                                
        cat2d(i,j)=ibest

	if (mod(i,10) .eq. 0 .and. mod(j,10) .eq. 0) then
!	write(6,*) 'i,j,cat2d: ', i,j,cat2d(i,j)
!	write(6,*) 'icount,imaxcount ', icount,imaxcount
	endif
	
	endif

	else ! pure n.n. case

!!! dont believe this section will/should ever be used

	cat2d(i,j)=datain(ii,ji)
	
	do n=1,ncat
	if (cat2d(i,j) .eq. n) then
	  cat3d(i,j,n)=1.0
	else
	  cat3d(i,j,n)=0.0
	endif
	enddo
	
	endif
                                                                                
        else
! outside of this strip...ignore
                                                                                
        endif
                                                                                
                                                                                
        enddo ! enddos for im,jm loop
        enddo

	end subroutine alt_categories
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

         subroutine alt_hemi_all(im,jm,glatr,glonr,path,
     &     type,ncat,res,cat3d,cat2d,categorical,outland,useland)

	integer counter(ncat)

        parameter(r2d=57.29577951,dtr=0.017453293)

        character*12 fnameout
        integer gds(200)
        character*2 uname
	character*1 type
        character*8 filename
	character*7 geoloc(2)
        character*(*) path
        character*256 fullname
        character*256 title,unit_name

        real cat3d(im,jm,ncat)
        real cat2d(im,jm),x(im,jm),y(im,jm)
        real glatr(im,jm),glonr(im,jm)
        real dum2d(im,jm)
	real, allocatable:: datalat(:,:),datalon(:,:)
	real, allocatable:: xpts(:,:),ypts(:,:),dum(:,:)
	integer, allocatable:: idata(:,:,:)

!mp	respecify with proper "kind" type statment????
	integer*2, allocatable:: data2d(:,:)

	integer*2 ignoreval

!	
	real outland(im,jm) ! (land=1, water=0)

	real leftlon, rightlon

	logical have_left_east, have_left_west
	logical have_right_east, have_right_west
	logical read_west, read_east, read_both
	logical categorical, pure_nn, useland
	logical redo(im,jm), topogen

!
	cat3d=-99999.
	cat2d=-99999.
!



	have_left_east=.false.
	have_left_west=.false.
	have_right_east=.false.
	have_right_west=.false.

	do j=1,jm

	if (glonr(1,j) .lt. 0) then
	  have_left_west=.true.
	elseif (glonr(1,j) .gt. 0) then
	  have_left_east=.true.
	endif

	if (glonr(im,j) .lt. 0) then
	  have_right_west=.true.
	elseif (glonr(im,j) .gt. 0) then
	  have_right_east=.true.
	endif

	enddo


!	from this information can gather the straddling situation, and
!	devise a way to find west and east longitude


	read_west=.false.
	read_east=.false.
	read_both=.false.

	geoloc(1)='1234567'
	geoloc(2)='1234567'

	if (have_left_west .and. have_right_east) then
	   write(6,*) 'straddles pm'
	   read_both=.true.
	geoloc(1)='90s180w'
	geoloc(2)='90s000e'
	elseif (have_left_east .and. have_right_west) then
	   write(6,*) 'straddles dl'
	   read_both=.true.
	geoloc(1)='90s180w'
	geoloc(2)='90s000e'
	else
	   read_both=.false.
	   write(6,*) 'purely in one hemisphere'
	   if (have_left_west .and. have_right_west) then
	      read_west=.true.
	geoloc(1)='90s180w'
	   elseif (have_left_east .and. have_right_east) then
	      read_east=.true.
	geoloc(1)='90s000e'
	   endif
	endif

        call s_len(path,lenp)
        title=path(1:lenp)//'header'
              lentd=index(title,' ')-1
      call jclget(29,title(1:lentd),'formatted',1,istatus)
      if(istatus .ne. 1)then
         write(6,*)'warning: proc_geodat_tiles opening header: check'
     1            ,'geog paths and header file'
         return
      endif
      read(29,*)iblksizo,no,isbego,iwbego,rwoff,rsoff
      print *,'title=',title(1:lentd)
      print *,'rwoff,rsoff = ',rwoff,rsoff
      print *,'isbego,iwbego=',isbego,iwbego
      print *,'iblksizo,no=',iblksizo,no
      close(29)

	write(6,*) 'iblksizo= ', iblksizo

	ilim=no
	jlim=no

        write(6,*) 'ilim, jlim: ', ilim,jlim
	if (allocated(datalat)) deallocate(datalat)
	if (allocated(datalon)) deallocate(datalon)
	if (allocated(xpts)) deallocate(xpts)
	if (allocated(ypts)) deallocate(ypts)
	if (allocated(dum)) deallocate(dum)

        allocate(datalat(ilim,jlim))
        allocate(datalon(ilim,jlim))
        allocate(xpts(ilim,jlim))
        allocate(ypts(ilim,jlim))
        allocate(dum(ilim,jlim))

        do j=1,jlim
        do i=1,ilim
        xpts(i,j)=i
        ypts(i,j)=j
        enddo
        enddo

	if (read_both) then
	  ntiles=2
        else
	  ntiles=1
	endif	

!!!!
!!!!
	do nn=1,ntiles
!!!!
!!!!

        nn1=ilim
        nn2=ilim
	nn4=ncat

	if (allocated(idata)) deallocate(idata)
	if (allocated(data2d)) deallocate(data2d)

	write(6,*) 'allocated with dims: ', nn4,nn1,nn2
        allocate (idata(nn4,nn1,nn2))

	write(6,*) 'geoloc(nn): ', geoloc(nn)

	if (geoloc(nn)(4:7) .eq. '180w') then
		 wstart=-180
	elseif (geoloc(nn)(4:7) .eq. '000e') then
		 wstart=0
	else
		write(6,*) 'not right'
		stop
	endif

	write(6,*) 'wstart= ', wstart

	unit_name=path(1:lenp)//geoloc(nn)
	write(6,*) 'unit_name= ', unit_name

        call s_len(unit_name,len)

	i1=1
	i2=4

!	write(6,*) 'nn1*nn2*nn4: ', nn1*nn2*nn4

!endianissue???
        call read_binary_field(idata,i1,i2,nn1*nn2*nn4,
     &          unit_name,len)

!        n=1
!        write(6,*) 'ndata: ', n
!        do j=nn2,1,-nn2/35
!        write(6,633) (idata(n,i,j),i=1,nn1,nn1/20)
!        enddo


  633	format(40i3)

	
cc      put this tile into longitudinal strip


!        write(6,*) 'ilim, jlim: ', ilim, jlim


cc      come up with an estimate of how many input 30 s points should
cc      be included.  consider roughly enough to average over the size of
cc      the gridbox.

       limit=1

cc	limit represents an intentionally over-large search radius to be
cc	used below to specify which input points define the output point

!!	are these defs proper??

	deltlat=float(iblksizo)/no
	deltlon=float(iblksizo)/no

	write(6,*) 'no,iblksizo, deltlat: ', no,iblksizo, deltlat

        limit=int((res/2.)/(deltlat*111.2))+2

        gds=0
	
	slat=isbego+rsoff
	wlat=wstart+rwoff

	rnlat=90-rsoff
	elat=wlat+180

        gds(1)=0
        gds(2)=ilim
        gds(3)=jlim
        gds(4)=int(slat*1000)
        gds(5)=int(wlat*1000)
        gds(6)=128
        gds(7)=int((rnlat)*1000)
        gds(8)=int(elat*1000)
	gds(9)=int(deltlat*1000)
	gds(10)=int(deltlon*1000)

	write(6,*) 'deltalat, deltalon: ', deltlat, deltlon

!	write(6,*) 'gds= ', gds
!	write(6,*) 'ilim*jlim= ', ilim*jlim

	call gdswiz(gds,1,ilim*jlim,-9999.,xpts,ypts,
     &	datalon,datalat,nret,0,dum,dum)

	do j=1,jlim
	 do i=1,ilim
	  if (datalon(i,j) .gt. 180) then
	   datalon(i,j)=datalon(i,j)-360.
	  endif
         enddo
	enddo


	if (deltlat .lt. 1) then
	pure_nn=.false.
	else
	pure_nn=.true.
	endif


	allocate(data2d(ilim,jlim))

        if (categorical) then

        write(6,*) 'calling alt_categories!! ', pure_nn
	write(6,*) 'dont think this is proper!!!'

                call alt_categories(im,jm,glatr,glonr,
     &          ilim,jlim,idata,datalat,datalon,pure_nn,cat3d,
     &          ncat,cat2d,gds,limit,res)

        else

	do n=1,ncat

	do j=1,jlim
	do i=1,ilim
	data2d(i,j)=idata(n,i,j)
	enddo
	enddo

!	always want to ignore zeros for albedo, green frac, slope

	ignoreval=0
	
	topogen=.false.

                call alt_interp(im,jm,glatr,glonr,
     &          ilim,jlim,data2d,datalat,
     &		datalon,pure_nn,cat2d,dum2d,
     &          dum2d,dum2d,gds,limit,res,ignoreval,topogen)

	
	if (useland) then

! check if any output land points have a zero value for field.
! if so, assign an average of nearest non-zero value points within
! search radius

	isrch=3

	do j=1,jm
	do i=1,im

	redo(i,j)=.false.

	if (outland(i,j) .eq. 1 .and. cat2d(i,j) .eq. 0) then

!! trouble point

	ipts=0
	sum=0.

	do jj=j-isrch,j+isrch
	do ii=i-isrch,i+isrch

	if (ii .ge. 1 .and. ii .le. im .and. 
     &		jj .ge. 1 .and. jj .le. jm) then

	if (cat2d(ii,jj) .gt. 0) then
	sum=sum+cat2d(ii,jj)
	ipts=ipts+1
	endif

	endif

	enddo
	enddo ! end ii, jj loops

	if (ipts .gt. 0) then
!	  write(6,*) 'based on  ', ipts, 'valid points'
!	  write(6,*) 'derived new value of ', sum/ipts
	  cat2d(i,j)=sum/ipts
	else
	  write(6,*) 'found no valid points at ' , i,j
	  redo(i,j)=.true.
	endif

	endif

	enddo
	enddo


!new

	

	do j=1,jm
	 do i=1,im
	isrch=3
	  do while (redo(i,j) .and. isrch .lt. im/2)

! spiral outward as needed to find a valid value
	isrch=isrch+1

	do jj=j-isrch,j+isrch
	do ii=i-isrch,i+isrch

	if (ii .ge. 1 .and. ii .le. im .and. 
     &		jj .ge. 1 .and. jj .le. jm) then

	if (cat2d(ii,jj) .gt. 0) then
	  cat2d(i,j)=cat2d(ii,jj)
	  write(6,*) 'defined ', i,j, 'to be: ', cat2d(i,j)
	  write(6,*) 'found value at search : ', isrch
	  redo(i,j)=.false.
	  goto 99
	endif

	endif

	enddo
	enddo


  99	continue

	enddo

	if (redo(i,j)) then
	write(6,*) 'giving up on point...'
	cat2d(i,j)=0.25 ! a reasonable default for albedo/gfrac ???
	write(6,*) 'setting to default val ', cat2d(i,j)
	endif

	enddo
	enddo
	endif ! endif useland block

!endnew

	do j=1,jm
	do i=1,im
	cat3d(i,j,n)=cat2d(i,j)
	enddo
	enddo

	enddo ! ncat loop

  467	format(25(f3.0,1x))

!	each month will be interpolated, then dumped into the 3d
!	array.  will it work for slope/maxsnowalb?
	
        endif ! categorical/interp branch

	enddo ! for tiles

ccccccccccccc

        deallocate(datalat,datalon,xpts,ypts,dum)

        end subroutine alt_hemi_all

c +------------------------------------------------------------------+
	subroutine alt_interp(im,jm,glatr,glonr,
     &          ilim,jlim,datain,datalat,datalon,
     &		pure_nn,out2d,slope_lt,
     &          slope_ln,stdev,gds,limit,testrad,
     &		ignoreval,topogen)


	real glatr(im,jm),glonr(im,jm)
	real datalat(ilim,jlim),datalon(ilim,jlim)
	integer*2 datain(ilim,jlim),ignoreval
	real out2d(im,jm),slope_lt(im,jm),slope_ln(im,jm)
	real stdev(im,jm)

	integer gds(200)

	logical pure_nn, topogen

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!	write(6,*) 'testrad = ', testrad

!	write(6,*) 'datain for alt_interp, lims ', ilim, jlim
!	do j=jlim,1,-jlim/20
!	write(6,647) (datain(i,j),i=1,ilim,ilim/15)
!	enddo
  647	format(25i3)
	

        do j=1,jm
        do i=1,im
                                                                                
c determine lat/lon of target (output) grid point

        glatd = glatr(i,j)
        glond = glonr(i,j)
                                                                                
        call ced_ij(glatd,glond,x,y,gds)
                                                                                
!! kludgy fix for dead space "between" strips (like -60.003 longitude)

        if (x .le. id2+0.5) then
        ii=int(x+0.5)
        else
!       write(6,*) 'kludgy fix applied...'
        ii=int(x)
        endif
                                                                                
        ji=int(y+0.5)
                                                                                
!       ii,ji represents nearest neighbor from input data array
                                                                                
        if (ii .ge. 1 .and. ii .le. ilim .and.
     +      ji .ge. 1 .and. ji .le. jlim ) then
                                                                                
                if (.not. pure_nn) then
                                                                                
        icount=0
	ignorcnt=0
        sum=0.
        icnt_lt=0
        icnt_ln=0
        sum_lt=0.
        sum_ln=0.
        do jjj=ji-limit,ji+limit
           do iii=ii-limit,ii+limit
                                                                                
        if (jjj .lt. 1 .or. jjj .gt. jlim .or.
     &      iii .lt. 1 .or. iii .gt. ilim) goto 47
                                                                                
        call greatcir(glatd,glond,
     &          datalat(iii,jjj),datalon(iii,jjj),dist)

        if ( dist .le. testrad/2.) then
	  if (datain(iii,jjj) .ne. ignoreval) then
              sum=sum+real(datain(iii,jjj))
              icount=icount+1
	  else
              ignorcnt=ignorcnt+1
	  endif

!	
                                                                                
        if (iii+1 .le. ilim .and. jjj+1 .le. jlim ) then
        icnt_ln=icnt_ln+1
        icnt_lt=icnt_lt+1
        sum_ln = sum_ln + 0.5*(datain(iii+1,jjj)-datain(iii,jjj)+
     &                         datain(iii+1,jjj+1)-datain(iii,jjj+1))
        sum_lt = sum_lt + 0.5*(datain(iii,jjj+1)-datain(iii,jjj)+
     &                         datain(iii+1,jjj+1)-datain(iii+1,jjj))
        endif
                                                                                
        endif
                                                                                
  47    continue
           enddo
        enddo
c
                                                                                
        slope_ln(i,j)=(sum_ln/float(icnt_ln))/(1000.*testrad)
        slope_lt(i,j)=(sum_lt/float(icnt_lt))/(1000.*testrad)
        out2d(i,j)=sum/float(icount)

	if (icount .eq. 0) then ! force nearest neighbor
	out2d(i,j)=datain(ii,ji)

!!	in this case, can feel comfortable that is not topo
!!	so dont worry about stdev and slope computations

	endif


!!	if predominantly zero values, dominant land/soil type will be water
!!	so set topo to zero here.  looking for better topo/landmask 
!!	agreement

	if (topogen .and. ignorcnt .gt. icount  .and. 
     &		out2d(i,j) .gt. 0) then
	write(6,*) 'set topo to zero from value ', out2d(i,j)
	out2d(i,j)=0.
	endif

                                                                                
!!      now essentially redo to get standard deviation info

!! 1002 - not using below for anything...comment out?
                                                                                
        icount=0
        sum=0.
                                                                                
        do jjj=ji-limit,ji+limit
           do iii=ii-limit,ii+limit
                                                                                
        if (jjj .ge. 1 .and. jjj .le. jlim .and.
     &      iii .ge. 1 .and. iii .le. ilim) then
                                                                                
        call greatcir(glatd,glond,
     &          datalat(iii,jjj),datalon(iii,jjj),dist)
                                                                                
        if ( dist .le. testrad/2.) then
              sum=sum+(out2d(i,j)-real(datain(iii,jjj)))**2.
              icount=icount+1
        endif
                                                                                
        endif
                                                                                
                                                                                
        enddo
        enddo
                                                                                
        stdev(i,j)=(sum/float(icount))**(0.5)
                                                                                
!! end standard deviation
                else ! pure n.n. case
                                                                                
                stdev(i,j)=0.
                                                                                
                out2d(i,j)=datain(ii,ji)
                slope_ln(i,j)=0.
                slope_lt(i,j)=0.
                                                                                
                endif
                                                                                
        else
! outside of this strip...ignore
                                                                                
        endif
                                                                                
        enddo ! enddos for im,jm loop
        enddo

	end subroutine alt_interp


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine ced_ij(rlat,rlon,xpts,ypts,kgds)

        integer kgds(200)

        im=kgds(2)
        jm=kgds(3)
        rlat1=kgds(4)*1.e-3
        rlon1=kgds(5)*1.e-3
        rlat2=kgds(7)*1.e-3
        rlon2=kgds(8)*1.e-3
        iscan=mod(kgds(11)/128,2)
        jscan=mod(kgds(11)/64,2)
        nscan=mod(kgds(11)/32,2)
        hi=(-1.)**iscan
        hj=(-1.)**(1-jscan)
        dlon=hi*(mod(hi*(rlon2-rlon1)-1+3600,360.)+1)/(im-1)
        dlat=(rlat2-rlat1)/(jm-1)
        xmin=0
        xmax=im+1
        if(im.eq.nint(360/abs(dlon))) xmax=im+2
        ymin=0
        ymax=jm+1
        nret=0
        lrot=0
        fill=-999

            if(abs(rlon).le.360.and.abs(rlat).le.90) then
              xpts=1+hi*mod(hi*(rlon-rlon1)+3600,360.)/dlon
              ypts=1+(rlat-rlat1)/dlat

              if(xpts.ge.xmin.and.xpts.le.xmax.and.
     &           ypts.ge.ymin.and.ypts.le.ymax) then
                nret=nret+1
                if(lrot.eq.1) then
                  crot=1
                  srot=0
                endif
              else
                xpts=fill
                ypts=fill
              endif
            else
              xpts=fill
              ypts=fill
            endif


        return
        end

c --------------------

        subroutine gen_filename(lat,lon,name)

        integer lat,lon
        character*7 name
        character*1 latswitch,lonswitch,type
        character*2 latval
        character*3 lonval

        if (lat .ge. 0) latswitch='n'
        if (lat .lt. 0) latswitch='s'

        if (lon .ge. 0) lonswitch='e'
        if (lon .lt. 0) lonswitch='w'


         write(latval,'(i2.2)') abs(lat)
         write(lonval,'(i3.3)') abs(lon)

        name=latval//latswitch//lonval//lonswitch

        end subroutine gen_filename

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	subroutine greatcir(rlat1,rlon1,rlat2,rlon2,dist)

        parameter (pie=3.141592654,rconv=3.141592654/180.)
        parameter (erad=6.3712e+6)

c	formulation from steers (1965)

        costerm=cos(rlat2*rconv)*cos(rlat1*rconv)
        absdelong=abs(rlon2-rlon1)*rconv
        havdelong=(1.-cos(absdelong))/2.
        colatf=rconv*(90-rlat1)
        colatt=rconv*(90-rlat2)
        factor=amax1(colatf,colatt)-amin1(colatf,colatt)
        havfactor=(1.-cos(factor))/2.
        havft=havdelong*costerm+havfactor
        angdist=acos(1-2*havft)/rconv
        distm=angdist*(2*pie*erad)/360.
        dist=distm/1000.

	end subroutine greatcir
