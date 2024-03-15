module prep_stmas4d

   use prmtrs_stmas

   private gtscaling
   public preprocss, grdmemalc, rdinitbgd, physcpstn

!**************************************************
!comment:
!   this module is used by stmas_core.f90 to do some preprocess.
!   subroutines:
!      preprocss: call subroutine of grdmemalc to allocate memorys to analysis field and call gtscaling to get scales of observations.
!      grdmemalc: memory allocate for analysis fields, background fields and the coordinates.
!      gtscaling: get scales of each observation.
!      rdinitbgd: get background filed on the current level grid from the initial background (with max grid number).
!      physcpstn: calculate penalty coefficents and make some scaling of control variables and coordinates.
!**************************************************

contains

   subroutine preprocss
!*************************************************
! data preprocess
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
      call grdmemalc
      call gtscaling
      return
   end subroutine preprocss

   subroutine grdmemalc
!*************************************************
! memory allocate for grd array
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, er
! --------------------
      allocate (www(numgrid(1), numgrid(2), numgrid(3), numgrid(4)), stat=er)
      if (er .ne. 0) stop 'www allocate wrong'
      allocate (cor(numgrid(1), numgrid(2)), stat=er)
      if (er .ne. 0) stop 'cor allocate wrong'
      allocate (xxx(numgrid(1), numgrid(2)), stat=er)
      if (er .ne. 0) stop 'xxx allocate wrong'
      allocate (yyy(numgrid(1), numgrid(2)), stat=er)
      if (er .ne. 0) stop 'yyy allocate wrong'
      allocate (deg(numgrid(1), numgrid(2)), stat=er)
      if (er .ne. 0) stop 'deg allocate wrong'
      allocate (den(numgrid(1), numgrid(2), numgrid(3), numgrid(4)), stat=er)
      if (er .ne. 0) stop 'den allocate wrong'
      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then
         allocate (zzz(numgrid(1), numgrid(2), numgrid(3), numgrid(4)), stat=er)
         if (er .ne. 0) stop 'zzz allocate wrong'
      elseif (ifpcdnt .eq. 1) then
         allocate (ppp(numgrid(3)), stat=er)
         if (er .ne. 0) stop 'ppp allocate wrong'
      end if

! yuanfu add numstat+1 grdbkgnd array for saving radar reflectivity genrerated low bound and upper bound:
      allocate (grdbkgnd(numgrid(1), numgrid(2), numgrid(3), numgrid(4), numstat + 2), stat=er)
      if (er .ne. 0) stop 'grdbkgnd allocate wrong'
      allocate (grdanals(numgrid(1), numgrid(2), numgrid(3), numgrid(4), numstat), stat=er)
      if (er .ne. 0) stop 'grdanals allocate wrong'
      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         do s = 1, numstat
            grdanals(i, j, k, t, s) = 0.0d0
            grdbkgnd(i, j, k, t, s) = 0.0d0
         end do
      end do
      end do
      end do
      end do
      numvars = numgrid(1)*numgrid(2)*numgrid(3)*numgrid(4)*numstat
      allocate (gradint(numgrid(1), numgrid(2), numgrid(3), numgrid(4), numstat), stat=er)
      if (er .ne. 0) stop 'gradint allocate wrong'
      return
   end subroutine grdmemalc

   subroutine gtscaling
!*************************************************
! allocate observation memory and read in data and scale
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      real, parameter :: sm = 1.0e-5
      integer  :: o, s, er, uu, vv, ww, no
!added by shuyuan liu 20101028 for adding the oi oi=obserror(o)**2  oi=1.0/oi
      real    ::  oi
! --------------------
      uu = u_cmpnnt
      vv = v_cmpnnt
      ww = w_cmpnnt

!jhui
!numstat+1
      allocate (scl(numstat + 1), stat=er)
      if (er .ne. 0) stop 'scl allocate wrong'
      do s = 1, numstat + 1
         scl(s) = 0.0d0
      end do

      if (nallobs .eq. 0) return
!!!!!!!!!!!!!!!!!!!!
      o = 0
!  do s=1,numstat
!    do no=1,nobstat(s)
      !     o=o+1
!      scl(s)=scl(s)+obsvalue(o)*obsvalue(o)
!    enddo
      ! enddo
      ! do s=numstat+1,numstat+2
!    do no=1,nobstat(s)
!      o=o+1
!      scl(uu)=scl(uu)+obsvalue(o)*obsvalue(o)
!      scl(vv)=scl(vv)+obsvalue(o)*obsvalue(o)
      !   enddo
!  enddo
!jhui
      ! do s=numstat+3,numstat+3
      !   do no=1,nobstat(s)
      !     o=o+1
      !     scl(numstat+1)=scl(numstat+1)+obsvalue(o)*obsvalue(o)
      !   enddo
      ! enddo
! statistic scales of each state variable based on observations.
!!!!!!!!!!adding the oi  modified by shuyuan liu 20101028
      o = 0
      do s = 1, numstat
         do no = 1, nobstat(s)
            o = o + 1
            oi = obserror(o)*obserror(o)
            oi = 1.0/oi
            scl(s) = scl(s) + obsvalue(o)*obsvalue(o)*oi
         end do
      end do
      do s = numstat + 1, numstat + 2
         do no = 1, nobstat(s)
            o = o + 1
            oi = obserror(o)*obserror(o)
            oi = 1.0/oi
            scl(uu) = scl(uu) + obsvalue(o)*obsvalue(o)*oi
            scl(vv) = scl(vv) + obsvalue(o)*obsvalue(o)*oi
         end do
      end do
      do s = numstat + 3, numstat + 3
         do no = 1, nobstat(s)
            o = o + 1
            oi = obserror(o)*obserror(o)
            oi = 1.0/oi
            scl(numstat + 1) = scl(numstat + 1) + obsvalue(o)*obsvalue(o)*oi
         end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do s = 1, numstat
         if (nobstat(s) .gt. 0 .and. scl(s) .ge. 1e-5) then
            if (s .eq. uu .or. s .eq. vv) then
               scl(s) = sqrt(scl(s)/(nobstat(s) + nobstat(numstat + 1) + nobstat(numstat + 2)))
            else
               scl(s) = sqrt(scl(s)/nobstat(s))
            end if
         elseif (nobstat(s) .eq. 0) then
            scl(s) = sl0(s)
         end if
      end do
      if (numdims .ge. 2) scl(uu) = max(scl(uu), scl(vv))
      if (numdims .ge. 2) scl(vv) = scl(uu)
      if (ww .ne. 0) scl(ww) = scl(uu)
!jhui
      if (nobstat(numstat + 3) .gt. 0) then  ! yuanfu: avoid division of zero
         do s = numstat + 3, numstat + 3
            scl(numstat + 1) = sqrt(scl(numstat + 1)/nobstat(s))
         end do
      end if

! scale the observations
      o = 0
      do s = 1, numstat
         do no = 1, nobstat(s)
            o = o + 1
            obsvalue(o) = obsvalue(o)/scl(s)
            obserror(o) = obserror(o)/scl(s)
         end do
      end do
      do s = numstat + 1, numstat + 2
         do no = 1, nobstat(s)
            o = o + 1
            obsvalue(o) = obsvalue(o)/scl(uu)
            obserror(o) = obserror(o)/scl(uu)
         end do
      end do
!jhui
      !do s=numstat+3,numstat+3  !changed by shuyuan 20100903
      s = numstat + 3
      do no = 1, nobstat(s)
         o = o + 1
         obsvalue(o) = obsvalue(o)/scl(numstat + 1)
         obserror(o) = obserror(o)/scl(numstat + 1)
      end do
      !enddo

      return
   end subroutine gtscaling

   subroutine rdinitbgd
!*************************************************
! read in initial background fields
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, nx, ny, nz, nt, i1, j1, k1, t1, tt
      real     :: r
! --------------------
      r = 287
      tt = temprtur

      if (numgrid(1) .ge. 2) then
         nx = (maxgrid(1) - 1)/(numgrid(1) - 1)
      else
         nx = 1
      end if
      if (numgrid(2) .ge. 2) then
         ny = (maxgrid(2) - 1)/(numgrid(2) - 1)
      else
         ny = 1
      end if
      if (numgrid(3) .ge. 2) then
         nz = (maxgrid(3) - 1)/(numgrid(3) - 1)
      else
         nz = 1
      end if
      if (numgrid(4) .ge. 2) then
         nt = (maxgrid(4) - 1)/(numgrid(4) - 1)
      else
         nt = 1
      end if
! interplate the background field onto the current analysis grid points.
      do t = 1, maxgrid(4), nt
      do k = 1, maxgrid(3), nz
      do j = 1, maxgrid(2), ny
      do i = 1, maxgrid(1), nx
         i1 = (i - 1)/nx + 1
         j1 = (j - 1)/ny + 1
         k1 = (k - 1)/nz + 1
         t1 = (t - 1)/nt + 1
!jhui
         do s = 1, numstat
            grdbkgnd(i1, j1, k1, t1, s) = grdbkgd0(i, j, k, t, s)/scl(s)
         end do
         ! make sure sh background positive:
         grdbkgnd(i1, j1, k1, t1, 5) = max(0.0, grdbkgnd(i1, j1, k1, t1, 5))
         ! sh lower and upper bound at the multigrid level:
         grdbkgnd(i1, j1, k1, t1, numstat + 1) = grdbkgd0(i, j, k, t, numstat + 1)/scl(humidity)
         grdbkgnd(i1, j1, k1, t1, numstat + 2) = grdbkgd0(i, j, k, t, numstat + 2)/scl(humidity)
         if (ifpcdnt .ne. 1) den(i1, j1, k1, t1) = dn0(i, j, k, t)     ! here the pressure is not aviable
        if(ifpcdnt.eq.1) den(i1,j1,k1,t1)=(ppp(k1)*scp(psl)+orivtcl)/r/((grdbkgnd(i1,j1,k1,t1,tt)+grdanals(i1,j1,k1,t1,tt))*scl(tt))
      end do
      end do
      end do
      end do

      return
   end subroutine rdinitbgd

   subroutine physcpstn
!*************************************************
! dealing with physical position and penalty coefficent and making some scaling
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, nx, ny, nz, nt, i1, j1, k1, t1, uu, vv, pp
      real     :: z1, z2
! --------------------
      uu = u_cmpnnt
      vv = v_cmpnnt
      pp = pressure
! physical position and coriolis frequency
      if (numgrid(1) .ge. 2) then
         nx = (maxgrid(1) - 1)/(numgrid(1) - 1)
      else
         nx = 1
      end if
      if (numgrid(2) .ge. 2) then
         ny = (maxgrid(2) - 1)/(numgrid(2) - 1)
      else
         ny = 1
      end if
      do j = 1, maxgrid(2), ny
      do i = 1, maxgrid(1), nx
         i1 = (i - 1)/nx + 1
         j1 = (j - 1)/ny + 1
         if (numgrid(1) .ge. 2) xxx(i1, j1) = xx0(i, j)
         if (numgrid(2) .ge. 2) yyy(i1, j1) = yy0(i, j)
         if (pnlt0pu .ge. 1.0e-10 .or. pnlt0pv .ge. 1.0e-10) cor(i1, j1) = cr0(i, j)
         if (numgrid(2) .ge. 2) deg(i1, j1) = dg0(i, j)
      end do
      end do

      if (numgrid(3) .ge. 2) then
         nz = (maxgrid(3) - 1)/(numgrid(3) - 1)
      else
         nz = 1
      end if
      if (numgrid(4) .ge. 2) then
         nt = (maxgrid(4) - 1)/(numgrid(4) - 1)
      else
         nt = 1
      end if

      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then          ! for sigma and heightcoordinate
         do t = 1, maxgrid(4), nt
         do k = 1, maxgrid(3), nz
         do j = 1, maxgrid(2), ny
         do i = 1, maxgrid(1), nx
            i1 = (i - 1)/nx + 1
            j1 = (j - 1)/ny + 1
            k1 = (k - 1)/nz + 1
            t1 = (t - 1)/nt + 1
            zzz(i1, j1, k1, t1) = zz0(i, j, k, t)
         end do
         end do
         end do
         end do
      elseif (ifpcdnt .eq. 1) then
         do k = 1, maxgrid(3), nz
            k1 = (k - 1)/nz + 1
            ppp(k1) = pp0(k)
         end do
      end if

! penalty coefficent
      do s = 1, numstat
         penal_x(s) = 0.0
         penal_y(s) = 0.0
         penal_z(s) = 0.0
         penal_t(s) = 0.0
         if (numgrid(1) .ge. 3) penal_x(s) = penal0x(s)*((maxgrid(1) - 2)*maxgrid(2)*maxgrid(3)*maxgrid(4)) &
                                             /((numgrid(1) - 2)*numgrid(2)*numgrid(3)*numgrid(4))
         if (numgrid(2) .ge. 3) penal_y(s) = penal0y(s)*((maxgrid(2) - 2)*maxgrid(1)*maxgrid(3)*maxgrid(4)) &
                                             /((numgrid(2) - 2)*numgrid(1)*numgrid(3)*numgrid(4))
         if (numgrid(3) .ge. 3) penal_z(s) = penal0z(s)*((maxgrid(3) - 2)*maxgrid(1)*maxgrid(2)*maxgrid(4)) &
                                             /((numgrid(3) - 2)*numgrid(1)*numgrid(2)*numgrid(4))
         if (numgrid(4) .ge. 3) penal_t(s) = penal0t(s)*((maxgrid(4) - 2)*maxgrid(1)*maxgrid(2)*maxgrid(3)) &
                                             /((numgrid(4) - 2)*numgrid(1)*numgrid(2)*numgrid(3))
      end do
      !  penal_x(numstat+1) = penal_x(1)
      !  penal_y(numstat+1) = penal_y(1)
      !  penal_z(numstat+1) = penal_z(1)
      !  penal_t(numstat+1) = penal_t(1)

      pnlt_pu = pnlt0pu*scl(uu)*scl(uu)/(2**(grdlevl - 1)) &
                /(numgrid(1)*numgrid(2)*numgrid(3)*numgrid(4))
      pnlt_pv = pnlt0pv*scl(vv)*scl(vv)/(2**(grdlevl - 1)) &
                /(numgrid(1)*numgrid(2)*numgrid(3)*numgrid(4))
      if (nallobs .ne. 0) then
         if (nobstat(uu) .eq. 0 .and. nobstat(pp) .eq. 0) pnlt_pu = 0.0d0
         if (nobstat(vv) .eq. 0 .and. nobstat(pp) .eq. 0) pnlt_pv = 0.0d0
      end if

      pnlt_hy = pnlt0hy*taul_hy**(grdlevl - 1) &
                /(numgrid(1)*numgrid(2)*numgrid(3)*numgrid(4))

! calculate scale for physical position and coriolis frequency
      scp(xsl) = abs((xxx(numgrid(1), numgrid(2)) - xxx(1, 1))/(numgrid(1) - 1))
      scp(ysl) = abs((yyy(numgrid(1), numgrid(2)) - yyy(1, 1))/(numgrid(2) - 1))
      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then         ! for sigma and heightcoordinate
         z2 = zzz(1, 1, numgrid(3), 1)
         z1 = zzz(1, 1, 1, 1)
         do i = 1, numgrid(1)
         do j = 1, numgrid(2)
         do t = 1, numgrid(4)
            z2 = max(z2, zzz(i, j, numgrid(3), t))
            z1 = min(z1, zzz(i, j, 1, t))
         end do
         end do
         end do
         scp(psl) = abs((z2 - z1)/(numgrid(3) - 1))
      elseif (ifpcdnt .eq. 1) then
         scp(psl) = abs((ppp(numgrid(3)) - ppp(1))/(numgrid(3) - 1))
      end if
      if (pnlt0pu .ge. 1.0e-10 .or. pnlt0pv .ge. 1.0e-10) scp(csl) = abs(cor(1, 1))
      do i = 1, numgrid(1)
      do j = 1, numgrid(2)
         if (pnlt0pu .ge. 1.0e-10 .or. pnlt0pv .ge. 1.0e-10) scp(csl) = max(scp(csl), abs(cor(i, j)))
      end do
      end do
      if (scp(csl) .lt. 1.0e-5) scp(csl) = 1.0 ! temporarily turn off scaling coroilis
      if (pnlt0pu .ge. 1.0e-10 .or. pnlt0pv .ge. 1.0e-10) scp(dsl) = abs(den(1, 1, 1, 1))
      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         if (pnlt0pu .ge. 1.0e-10 .or. pnlt0pv .ge. 1.0e-10) scp(dsl) = max(scp(dsl), abs(den(i, j, k, t)))
      end do
      end do
      end do
      end do
      if (scp(dsl) .lt. 1.0e-5) scp(dsl) = 1.0 ! temporarily turn off scaling geostrophic with density
! scale the physical position and coriolis frequency
      do i = 1, numgrid(1)
      do j = 1, numgrid(2)
         if (numgrid(1) .ge. 2) xxx(i, j) = (xxx(i, j) - oripstn(1))/scp(xsl)
         if (numgrid(2) .ge. 2) yyy(i, j) = (yyy(i, j) - oripstn(2))/scp(ysl)
         if (pnlt0pu .ge. 1.0e-10 .or. pnlt0pv .ge. 1.0e-10) cor(i, j) = cor(i, j)/scp(csl)
      end do
      end do
      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         if (pnlt0pu .ge. 1.0e-10 .or. pnlt0pv .ge. 1.0e-10) den(i, j, k, t) = den(i, j, k, t)/scp(dsl)
      end do
      end do
      end do
      end do
      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then         ! for sigma and height coordinate
         orivtcl = zzz(1, 1, 1, 1)
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            zzz(i, j, k, t) = (zzz(i, j, k, t) - orivtcl)/scp(psl)
         end do
         end do
         end do
         end do
      elseif (ifpcdnt .eq. 1) then                       ! for pressure coordinate
         orivtcl = ppp(1)
         do k = 1, numgrid(3)
            ppp(k) = (ppp(k) - orivtcl)/scp(psl)
         end do
      end if

      return
   end subroutine physcpstn

end module prep_stmas4d
