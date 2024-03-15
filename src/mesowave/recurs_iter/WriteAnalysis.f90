subroutine writeanalysis(a,n)

!==========================================================
!  this routine is to write the analysis out in laps format.
!
!  history: may. 2004 by yuanfu xie.
!==========================================================

  implicit none

  real, parameter :: temp0 = 273.16,capr = 287.0, cp=1004.0
  integer, intent(in) :: n(4)
  real,    intent(in) :: a(n(1),n(2),n(3),n(4))
  real :: te,ept	! theta_e
  real :: dxy(nx,ny,2),flu(n(1)),flv(n(2))
  real :: mjohn,r_missing_data

  character :: ext*3,varnames(nvlaps+10)*3,vunits(nvlaps+10)*3
  character :: lvl_coord(nvlaps+10)*3,comment(nvlaps+10)*125
  integer   :: lvl(nvlaps+10),istatus,it,i,j
  real      :: data(n(1)-2*nfic,n(2)-2*nfic,16)
  real	    :: pp(n(1),n(2))

  call get_directory('lsx', dir_s, len)
  ext = 'lsx'

  print*,'data directory: ',dir_s(1:len)

  !print*,'max u: ',maxval(a(1:n(1),1:n(2),n(3),2))
  !print*,'min u: ',minval(a(1:n(1),1:n(2),n(3),2))

  ! variable names and units:
  varnames(1) = 't  '	! temperature
  vunits(1) = 'k  '
  varnames(2) = 'u  '	! u
  vunits(2) = 'm/s'
  varnames(3) = 'v  '	! v
  vunits(3) = 'm/s'
  varnames(4) = 'vis' 	! visibility (m)
  ! comment(4) = ''
  vunits(4) = 'm'
  varnames(5) = 'td '	! dew point
  vunits(5) = 'k  '
  varnames(6) = 'p  '	! reduced pressure
  vunits(6) = 'pa  '
  comment(6)(1:21) = '0  m reduced pressure'
  varnames(7) = 'the'	! theta_e
  vunits(7) = 'k  '
  varnames(8) = 'div'	! divergence
  vunits(8) = '/s '
  varnames(9) = 'th '	! potential temperature
  vunits(9) = 'k  '
  varnames(10) = 'mrc '	! moisture convergence
  vunits(10) = '/s  '
  varnames(11) = 'pp  '	! reduced pressure change
  vunits(11) = 'pa  '
  comment(11)(1:28) = '0  m reduced pressure change'

  do it=1,nvlaps+9
     lvl(it) = 0
     lvl_coord(it) = 'agl'
  enddo

  print*,'n v n(4) = ',n(4)

  dxy = 5000.0

  do it=n(3)-1,n(3)

     ! gridpoints/meter:
     mjohn = 1.0/sqrt(grid_spacingx)

     ! divergence:
     data(1:nx,1:ny,n(4)+2) = &
	(a(nfic+2:n(1)-nfic+1,nfic+1:n(2)-nfic,it,2)- &
         a(nfic  :n(1)-nfic-1,nfic+1:n(2)-nfic,it,2))/grid_spacingx*0.5+ &
        (a(nfic+1:n(1)-nfic,nfic+2:n(2)-nfic+1,it,3)- &
         a(nfic+1:n(1)-nfic,nfic  :n(2)-nfic-1,it,3))/grid_spacingy*0.5

     ! convert to real*4:
     data(1:nx,1:ny,1:n(4)) = &
         a(1+nfic:n(1)-nfic,1+nfic:n(2)-nfic,it,1:n(4))

     ! theta_e: equivalent potential temperature
     do j=1,ny
	do i=1,nx
	   data(i,j,n(4)+1) = ept(data(i,j,1)-temp0, &
                                  data(i,j,5)-temp0, &
                                  data(i,j,4)/100.0)+temp0
	enddo
     enddo

     ! theta:
     data(1:nx,1:ny,n(4)+3) = data(1:nx,1:ny,1)* &
	                (100000.0/data(1:nx,1:ny,4))**(capr/cp)

     ! flux convergence:
     data(1:nx,1:ny,16) = data(1:nx,1:ny,6)/100.0
     call meso_anl(data(1,1,2),data(1,1,3),data(1,1,16),data(1,1,1), &
	           data(1,1,5),data(1,1,n(4)+3),dxy(1,1,1),dxy(1,1,2), &
	           data(1,1,n(4)+5),data(1,1,n(4)+4),data(1,1,n(4)+6), &
                   data(1,1,n(4)+7),data(1,1,n(4)+8),nx,ny)

     ! reduced pressure change:
     ! data(1:nx,1:ny,11) = 0.0
     call gridbarnes(a(1,1,it,6),n,n,pp)
     data(1:nx,1:ny,11) = a(1+nfic:n(1)-nfic,1+nfic:n(2)-nfic,it,6)-&
			  pp(1+nfic:n(1)-nfic,1+nfic:n(2)-nfic)

     call write_laps_data(istarttime+(it-1)*laps_cycle_time, &
                          dir_s,ext,nx,ny,n(4)+5,n(4)+5,varnames, &
                          lvl,lvl_coord,vunits,comment,data,istatus)
     write(6,*)' lsx file write completed, istatus = ',istatus,nx,ny

  enddo
  
end subroutine
