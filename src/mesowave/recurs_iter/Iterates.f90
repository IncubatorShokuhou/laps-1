subroutine iterates(id,bkgd,ldf,nx,ny,ds,ncycles,nvlaps,nfic)

!*************************************************
!  this routine iteratively solves data analysis
!  problem.
!
!  history: feb. 2004 by yuanfu xie.
!           sep. 2004 by yuanfu xie penalizing div.
!*************************************************

  implicit none

  integer, intent(in) :: id,nx,ny,ncycles,nvlaps,nfic
  real,    intent(in) :: bkgd(nx,ny,ncycles,nvlaps)
  real,    intent(in) :: ldf(nx,ny),ds(3)

  integer :: iter,iobs,i,j,k,no_v,idp,nbqc,numo
  real    :: y0,b(2,3),rms,stdv(nvlaps),amx,amn,tmx,tmn

  real :: oo(4,mobs)
  oo = o

  ! unified analysis of velocity:
  idp = id
  if (id .eq. 201) idp = id+1

  s(1:n(1),1:n(2),1:n(3),id:idp) = 0.0

  do iter=1,nrf(id)

     a(1:n(1),1:n(2),1:n(3),id:idp) = 0.0

     ! qc: bound check and compute standard deviations:
     if (iter .eq. 1) then

	stdv(id) = 0.0
	no_v = 0

	nbqc = 0
	numo = 0
	amx = -1.0e6
	amn = 1.0e6

        do iobs=1,nobs

	   if ((id .eq. vid(iobs)) .and. &
	       (idx(1,iobs) .gt. nfic) .and. &
	       (idx(1,iobs) .lt. n(1)-nfic) .and. &
	       (idx(2,iobs) .gt. nfic) .and. &
               (idx(2,iobs) .lt. n(2)-nfic)) then

              y0 = 0.0

	      numo = numo+1
        
              b(1,1:3) = 1.0-coe(1:3,iobs)
              b(2,1:3) = coe(1:3,iobs)
              do k=1,2
		 if ((idx(3,iobs)+k-1 .ge. 1) .and. &
		     (idx(3,iobs)+k-1 .le. n(3))) then
                 do j=1,2
                    do i=1,2
                       y0 = y0 + bkgd(idx(1,iobs)+i-1-nfic, &
                                      idx(2,iobs)+j-1-nfic, &
                                      idx(3,iobs)+k-1,vid(iobs))* &
                                 b(i,1)*b(j,2)*b(k,3)
                    enddo   
                 enddo
		 endif
              enddo

	      stdv(id) = stdv(id)+(o(1,iobs)-y0)**2
	      no_v = no_v+1

	      if (abs(o(1,iobs)-y0) .gt. qc_cons(id)) then
		 print*,'threshold qc: ',o(1,iobs),y0,id,o(2:4,iobs)
		 o(1,iobs) = y0
		 w(iobs) = 0.0
	      endif

	      ! bad qc:
              if ((id .eq. 1) .or. (id .eq. 5)) then
                 if (abs(o(1,iobs)-y0) .gt. 10.0) then
	            !print*,'bad qc: ',o(1,iobs),y0, &
			! vid(iobs),iobs,o(2:4,iobs)
	            o(1,iobs) = y0
                    w(iobs) = 0.0
		    nbqc = nbqc+1
	         !else
		 !   print*,'god qc: ',o(1,iobs),y0, &
		 !	vid(iobs),iobs,o(2:4,iobs)
                 endif
	      endif

	   endif

        enddo

	print*,'total bad qc: ',nbqc,id,numo
	if (id .eq. 4) print*,'maximum diff: ',amx

	stdv(id) = sqrt(stdv(id)/no_v)

	!if (no_v .gt. 0) then
        !   print*,'standard deviation: ',stdv(id),id
	!else
	!   print*,'standard deviation: no observation'
	!endif

	! standard deviation check: with 4.0*stdv
	nbqc = 0
        do iobs=1,nobs

	   if ((id .eq. vid(iobs)) .and. &
	       (idx(1,iobs) .gt. nfic) .and. &
	       (idx(1,iobs) .lt. n(1)-nfic) .and. &
	       (idx(2,iobs) .gt. nfic) .and. &
               (idx(2,iobs) .lt. n(2)-nfic)) then

              y0 = 0.0
        
              b(1,1:3) = 1.0-coe(1:3,iobs)
              b(2,1:3) = coe(1:3,iobs)
              do k=1,2
		 if ((idx(3,iobs)+k-1 .ge. 1) .and. &
		     (idx(3,iobs)+k-1 .le. n(3))) then
                 do j=1,2
                    do i=1,2
                       y0 = y0 + bkgd(idx(1,iobs)+i-1-nfic, &
                                      idx(2,iobs)+j-1-nfic, &
                                      idx(3,iobs)+k-1,vid(iobs))* &
                                 b(i,1)*b(j,2)*b(k,3)
                    enddo   
                 enddo
		 endif
              enddo

              if (abs(o(1,iobs)-y0) .gt. 3.0*stdv(id)) then
	         !print*,'standard deviation qc: ', &
	    	 !   o(1,iobs),y0,vid(iobs),iobs,stdv(id)
	         nbqc = nbqc+1
	         o(1,iobs) = y0
                 w(iobs) = 0.0
	      !else
	      !  print*,'pass std: ',iobs,o(1,iobs),y0,3.0*stdv(id)
              endif

	   endif

        enddo

	print*,'number of std qced: ',nbqc,id,stdv(id),numo

     endif

     if (iter .gt. 1) then
        call minimize(id,ds)
     else
	a(nfic+1:n(1)-nfic,nfic+1:n(2)-nfic,1:n(3),id:idp) = &
                                 bkgd(1:nx,1:ny,1:n(3),id:idp)

	! fictitious points:
	do i=1,nfic
	   a(i,nfic+1:n(2)-nfic,1:n(3),id:idp) = &
                      a(nfic+1,nfic+1:n(2)-nfic,1:n(3),id:idp)
	   a(n(1)-nfic+i,nfic+1:n(2)-nfic,1:n(3),id:idp) = &
                      a(n(1)-nfic,nfic+1:n(2)-nfic,1:n(3),id:idp)
	enddo
	do i=1,nfic
	   a(1:n(1),i,1:n(3),id:idp) = a(1:n(1),nfic+1,1:n(3),id:idp)
	   a(1:n(1),n(2)-nfic+i,1:n(3),id:idp) = &
             a(1:n(1),n(2)-nfic,1:n(3),id:idp)
	enddo
     endif

     no_v = 0
     rms = 0.0

     do iobs=1,nobs

	if ((id .eq. vid(iobs)) .or. (idp .eq. vid(iobs))) then

           y0 = 0.0
        
           b(1,1:3) = 1.0-coe(1:3,iobs)
           b(2,1:3) = coe(1:3,iobs)
           do k=1,2
	      if ((idx(3,iobs)+k-1 .ge. 1) .and. &
		  (idx(3,iobs)+k-1 .le. n(3))) then
              do j=1,2
                 do i=1,2
                    y0 = y0 + a(idx(1,iobs)+i-1,idx(2,iobs)+j-1, &
                                idx(3,iobs)+k-1,vid(iobs))* &
                              b(i,1)*b(j,2)*b(k,3)
                 enddo
              enddo
	      endif
           enddo
        
           o(1,iobs) = o(1,iobs)-y0

	   !rms = rms + o(1,iobs)*o(1,iobs)
           !no_v = no_v + 1

	endif

     enddo

     !if (no_v .ne. 0) then
     !   write(11,*) 'rms: ', sqrt(rms/no_v),iter,al(1:3,id)
     !else
     !   write(11,*) 'rms: ', sqrt(rms),iter,al(1:3,id)
     !endif

     ! filter:
     if (iter .gt. 1) then
        if (id .eq. 6) then
	   al(1:3,id) = al(1:3,id)*0.9
        else
           al(1:3,id:idp) = al(1:3,id:idp)*0.9
	endif
     endif

     ! accumulate:
     s(1:n(1),1:n(2),1:n(3),id:idp) = &
	s(1:n(1),1:n(2),1:n(3),id:idp)+ &
        a(1:n(1),1:n(2),1:n(3),id:idp)

  enddo

  ! temporature problem:
  if (id .eq. 1) then

     ! compute the max/min differences from background:
     tmx = maxval(s(nfic+1:nfic+nx,nfic+1:nfic+ny,n(3),1)- &
	       bkgd(1:nx,1:ny,n(3),1))
     tmn = minval(s(nfic+1:nfic+nx,nfic+1:nfic+ny,n(3),1)- &
	       bkgd(1:nx,1:ny,n(3),1))

     print*,'max diff in t: ', tmx, tmn

     if ((tmx .gt. 15.0e0) .or. (tmn .lt. -15.0e0)) then

	open(unit=11,file='badata.dat')
	write(11,*) nx,ny,n,dm(1,3),dm(2,3)
	write(11,*) bkgd(1:n(1),1:n(2),1:n(3),1)
	write(11,*) ldf(1:nx,1:ny)
	write(11,*) nobs
	write(11,*) oo(1:4,1:nobs),vid(1:nobs),w(1:nobs)
	write(11,*) s(1:n(1),1:n(2),1:n(3),1)
	close(11)
	
        ! do j=1,ny
        !   do i=1,nx
	!      s(nfic+i,nfic+j,1:n(3),id) = &
        !        bkgd(i,j,1:n(3),id)
        !   enddo
        ! enddo

     endif
  endif

  ! land/water weight:
  if ((id .ne. 6) .and. (id .ne. 4)) then
     do j=1,ny
        do i=1,nx
  	   s(nfic+i,nfic+j,1:n(3),id:idp) = &
		ldf(i,j)*s(nfic+i,nfic+j,1:n(3),id:idp)+ &
                 (1.0-ldf(i,j))*bkgd(i,j,1:n(3),id:idp)
        enddo
     enddo
  else 
     ! do j=1,ny
     !    do i=1,nx
        ! s(nfic+i,nfic+j,1:n(3),id) = bkgd(i,j,1:n(3),id)
     !    enddo
     ! enddo
  endif

end subroutine iterates
