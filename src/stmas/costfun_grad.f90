module costfun_grad
!*************************************************
! calculate cost function and the relevent gradients
! history: departed from stmas4d_core module by zhongjie he, august 2007.
!*************************************************

   use prmtrs_stmas
   use wcompt_gradt
   use smcostf_grad
   use gsbcost_grad
   use hydcost_grad, only: hydrocost, hydrograd, hydrocost_xie, hydrograd_xie, hydrocost_shuyuan, hydrograd_shuyuan

   public costfunct, costgradt, costfunct1, costgradt1, costfunct2, costgradt2, cg_value, cg_grad

!***************************************************
!!comment:
!   this module is used by the module of stmas4d_core, the function is to calculate the value of costfunction and gradients.
!   subroutines:
!      costfunct: calculate the value of costfunction, while for the radial wind data, the vertical velocity is not considered.
!      costgradt: calculate the gradient of costfunctions correspond to costfunct.
!      costfunct1: calculate the value of costfunction, for the case that vertical velocitys are state variables.
!      costgradt1: calculate the gradient of costfunctions correspond to costfunct1.
!      costfunct2: just like the subroutine of costfunct, while the vertical velocity is considered, but not state variables.
!      costgradt2: calculate the gradient of costfunctions correspond to costfunct2.
!***************************************************
contains

   subroutine costfunct1
!*************************************************
! calculate cost function
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, o, no, m, n, uu, vv, ww
      integer  :: np(maxdims)
      real  :: ht, oi, hu, hv, hw, dg, dv
      real  :: ht0, hu0, hv0, hw0
      real  :: cs(numstat), cm, cb, ch
      real  :: d2r

      d2r = 3.14159/180.0

! --------------------
      uu = u_cmpnnt
      vv = v_cmpnnt
      ww = w_cmpnnt
! initialization
      costfun = 0.0
      do s = 1, numstat
         cs(s) = 0.0
      end do
      if (nallobs .eq. 0) return
      o = 0
      do s = 1, numstat
         do no = 1, nobstat(s)
            o = o + 1
            do n = 1, maxdims
               np(n) = obsidxpc(n, o)
            end do
            oi = obserror(o)**2
            oi = 1.0/oi
            ht = 0.0
            m = 0
            do t = np(4), min0(np(4) + 1, numgrid(4))
            do k = np(3), min0(np(3) + 1, numgrid(3))
            do j = np(2), min0(np(2) + 1, numgrid(2))
            do i = np(1), min0(np(1) + 1, numgrid(1))
               m = m + 1
               ht = ht + obscoeff(m, o)*grdanals(i, j, k, t, s)
            end do
            end do
            end do
            end do
            costfun = costfun + (ht - obsvalue(o))*oi*(ht - obsvalue(o))/(nobstat(s)*1.0)
            cs(s) = cs(s) + 0.5*(ht - obsvalue(o))*oi*(ht - obsvalue(o))/(nobstat(s)*1.0)
         end do
      end do
      s = numstat + 1                                  !   for radar data
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)**2
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*grdanals(i, j, k, t, ww)
         end do
         end do
         end do
         end do
         dg = obseinf1(no)
         dv = obseinf2(no)
         ht = hu*sin(d2r*dg)*cos(d2r*dv) + hv*cos(d2r*dg)*cos(d2r*dv) + hw*sin(d2r*dv)
         costfun = costfun + (ht - obsvalue(o))*oi*(ht - obsvalue(o))*obsradar
      end do
      s = numstat + 2                        ! for srmr data
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)**2
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         hu0 = 0.0
         hv0 = 0.0
         hw0 = 0.0
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*grdanals(i, j, k, t, ww)
            hu0 = hu0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, uu)
            hv0 = hv0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, vv)
            hw0 = hw0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, ww)
         end do
         end do
         end do
         end do
         ht = sqrt((hu + hu0)*(hu + hu0) + (hv + hv0)*(hv + hv0) + (hw + hw0)*(hw + hw0))
         costfun = costfun + (ht - obsvalue(o))*oi*(ht - obsvalue(o))*obs_sfmr
      end do

      costfun = 0.5d0*costfun

! smooth term
      cm = costfun
      call smothcost
      cm = costfun - cm

! geostrophic balance term
      if (grdlevl .le. endgslv) then
         cb = costfun
         if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then           ! for sigma and height coordinate
            call gsblncost_non_uniform_z
         elseif (ifpcdnt .eq. 1) then                         ! for pressure coordinate
!    call gsblncost_p_hs
            call cntnscost
         end if
         cb = costfun - cb
      end if

!  write(100,*)(cs(s),s=1,numstat),cm,cb

! hydrostatic condition term                 ! added by zhongjie he
      if (grdlevl .le. endhylv) then
         ch = costfun
         call hydrocost
         ch = costfun - ch
      end if

      return
   end subroutine costfunct1

   subroutine costgradt1
!*************************************************
! calculate gradient of cost function
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, o, no, m, n, uu, vv, ww
      integer  :: np(maxdims)
      real     :: ht, oi, hu, hv, hw, dg, dv
      real     :: ht0, hu0, hv0, hw0
      real     :: d2r

      d2r = 3.14159/180.0

! --------------------
      uu = u_cmpnnt
      vv = v_cmpnnt
      ww = w_cmpnnt
! initialization
      do s = 1, numstat
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            gradint(i, j, k, t, s) = 0.0
         end do
         end do
         end do
         end do
      end do
      if (nallobs .eq. 0) return
      o = 0
      do s = 1, numstat
         do no = 1, nobstat(s)
            o = o + 1
            do n = 1, maxdims
               np(n) = obsidxpc(n, o)
            end do
            oi = obserror(o)**2
            oi = 1.0/oi
            ht = 0.0
            m = 0
            do t = np(4), min0(np(4) + 1, numgrid(4))
            do k = np(3), min0(np(3) + 1, numgrid(3))
            do j = np(2), min0(np(2) + 1, numgrid(2))
            do i = np(1), min0(np(1) + 1, numgrid(1))
               m = m + 1
               ht = ht + obscoeff(m, o)*grdanals(i, j, k, t, s)
            end do
            end do
            end do
            end do
            m = 0
            do t = np(4), min0(np(4) + 1, numgrid(4))
            do k = np(3), min0(np(3) + 1, numgrid(3))
            do j = np(2), min0(np(2) + 1, numgrid(2))
            do i = np(1), min0(np(1) + 1, numgrid(1))
               m = m + 1
               gradint(i, j, k, t, s) = gradint(i, j, k, t, s) &
                                        + (ht - obsvalue(o))*oi*obscoeff(m, o)/(nobstat(s)*1.0)
            end do
            end do
            end do
            end do
         end do
      end do
      s = numstat + 1                         !  for radar data
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)**2
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*grdanals(i, j, k, t, ww)
         end do
         end do
         end do
         end do
         dg = obseinf1(no)
         dv = obseinf2(no)
         ht = hu*sin(d2r*dg)*cos(d2r*dv) + hv*cos(d2r*dg)*cos(d2r*dv) + hw*sin(d2r*dv)
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            gradint(i, j, k, t, uu) = gradint(i, j, k, t, uu) &
                                      + (ht - obsvalue(o))*oi*obscoeff(m, o)*sin(d2r*dg)*cos(d2r*dv)*obsradar
            gradint(i, j, k, t, vv) = gradint(i, j, k, t, vv) &
                                      + (ht - obsvalue(o))*oi*obscoeff(m, o)*cos(d2r*dg)*cos(d2r*dv)*obsradar
            gradint(i, j, k, t, ww) = gradint(i, j, k, t, ww) &
                                      + (ht - obsvalue(o))*oi*obscoeff(m, o)*sin(d2r*dv)*obsradar
         end do
         end do
         end do
         end do
      end do
      s = numstat + 2                !  for sfmr observation
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)*obserror(o)
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         hu0 = 0.0
         hv0 = 0.0
         hw0 = 0
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*grdanals(i, j, k, t, ww)
            hu0 = hu0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, uu)
            hv0 = hv0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, vv)
            hw0 = hw0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, ww)
         end do
         end do
         end do
         end do
         ht = sqrt((hu + hu0)*(hu + hu0) + (hv + hv0)*(hv + hv0) + (hw + hw0)*(hw + hw0))
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            gradint(i, j, k, t, uu) = gradint(i, j, k, t, uu) &
                                      + (ht - obsvalue(o))*oi*obscoeff(m, o)*(hu + hu0)/ht*obs_sfmr
            gradint(i, j, k, t, vv) = gradint(i, j, k, t, vv) &
                                      + (ht - obsvalue(o))*oi*obscoeff(m, o)*(hv + hv0)/ht*obs_sfmr
            gradint(i, j, k, t, ww) = gradint(i, j, k, t, ww) &
                                      + (ht - obsvalue(o))*oi*obscoeff(m, o)*(hw + hw0)/ht*obs_sfmr
         end do
         end do
         end do
         end do
      end do

! smooth term
      call smothgrad

! geostrophic balance term
      if (grdlevl .le. endgslv) then
         if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then       ! for sigma and height coordinate
            call gsblngrad_non_uniform_z
         elseif (ifpcdnt .eq. 1) then                     ! for pressure coordinate
!    call gsblngrad_p_hs
            call cntnsgrad
         end if
      end if

! hydrostatic condition term                 ! added by zhongjie he
      if (grdlevl .le. endhylv) then
         call hydrograd
      end if

      return
   end subroutine costgradt1

   subroutine costfunct2
!*************************************************
! calculate cost function
! history: august 2007, coded by wei li.
! history: modified by zhongjie he.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, o, no, m, n, uu, vv
      integer  :: np(maxdims)
      real  :: ht, oi, hu, hv, hw, dg, dv
      real  :: ht0, hu0, hv0
      real  :: cs(numstat + 3), cm, cb, ch, ref
      double precision :: rdr
! added by shuyuan 20100907  for calculate ref, qr,qs..
      real  ::rc     !   desity*rain water mixing ration
      real  ::sc     !   desity*snow water mixing ration
      real  :: segma     !variance  for ref_obs
      real  ::temp_ref, tempr, temps
      integer  :: rr, rs
      real  :: d2r

      d2r = 3.14159/180.0

      uu = u_cmpnnt
      vv = v_cmpnnt
!added by shuyuan 20100728
      rr = rour_cmpnnt
      rs = rous_cmpnnt

!  call wcompgernl        !  modified by zhongjie he.
! initialization
      costfun = 0.0
      do s = 1, numstat + 3
         cs(s) = 0.0
      end do
      if (nallobs .eq. 0) return
      o = 0
      do s = 1, numstat
         do no = 1, nobstat(s)
            o = o + 1
            do n = 1, maxdims
               np(n) = obsidxpc(n, o)
            end do
            oi = obserror(o)**2
            oi = 1.0/oi
            ht = 0.0
            m = 0
            do t = np(4), min0(np(4) + 1, numgrid(4))
            do k = np(3), min0(np(3) + 1, numgrid(3))
            do j = np(2), min0(np(2) + 1, numgrid(2))
            do i = np(1), min0(np(1) + 1, numgrid(1))
               m = m + 1
               ht = ht + obscoeff(m, o)*grdanals(i, j, k, t, s)
            end do
            end do
            end do
            end do
            costfun = costfun + (ht - obsvalue(o))*(ht - obsvalue(o)) !/(nobstat(s)*1.0)!!!modified by shuyuan   20101028
            cs(s) = cs(s) + 0.5*(ht - obsvalue(o))*(ht - obsvalue(o)) !/(nobstat(s)*1.0)!!!modified by shuyuan   20101028
         end do
      end do

      ! for radial wind velocity observations
      s = numstat + 1
      rdr = 0.0d0 !costfun
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)**2
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*www(i, j, k, t)
         end do
         end do
         end do
         end do
         dg = obseinf1(no)
         dv = obseinf2(no)
         ht = hu*sin(d2r*dg)*cos(d2r*dv) + hv*cos(d2r*dg)*cos(d2r*dv) + hw*sin(d2r*dv)
         ! costfun=costfun+(ht-obsvalue(o))*(ht-obsvalue(o)) !*obsradar!!!modified by shuyuan   20101028
         rdr = rdr + (ht - obsvalue(o))*(ht - obsvalue(o))*obsradar!!!modified by shuyuan   20101028
      end do
      print *, 'radar radial wind cost: ', rdr, ' obsradr= ', obsradar, ' : ', costfun, rdr
      cs(s) = rdr
      costfun = costfun + rdr

      ! for sfmr wind velocity observations
      s = numstat + 2
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)**2
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         hu0 = 0.0
         hv0 = 0.0
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*www(i, j, k, t)
            hu0 = hu0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, uu)
            hv0 = hv0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, vv)
         end do
         end do
         end do
         end do
         ht = sqrt((hu + hu0)*(hu + hu0) + (hv + hv0)*(hv + hv0) + hw*hw)
         costfun = costfun + (ht - obsvalue(o))*(ht - obsvalue(o))!*obs_sfmr!!!modified by shuyuan   20101028
      end do

      ! reflectivity:
      ref = costfun
! --------------------
! added by shuyuan 20100721
      ! check if rain analysis is requested:
      if (numstat .le. 5) goto 555
      segma = 1.
      s = numstat + 3
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)**2
         oi = 1.0/oi
         temp_ref = 0.0
         m = 0
         rc = 0.
         sc = 0.
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))

            m = m + 1
            rc = rc + obscoeff(m, o)*grdanals(i, j, k, t, rr)
         end do
         end do
         end do
         end do
         !changed 20100907  shuyuan

         ! temp_ref=(obsvalue(o)-43.1)/17.5
         ! temp_ref=(10.**temp_ref)              !rc  unit is g/m3
         temp_ref = obsvalue(o)
         temp_ref = (rc - temp_ref)*(rc - temp_ref)

         !costfun=costfun+temp_ref*oi/(segma*segma)/(nobstat(s)*1.0) !!!!20101025
         costfun = costfun + temp_ref!!!!20101025

         cs(s) = cs(s) + 0.5*temp_ref

      end do
      print *, 'radar reflectivity cost: ', costfun - ref, obsvalue(o - 1), nobstat(s), s, segma

      ! skip reflectivity:
555   continue

      costfun = 0.5d0*costfun

! smooth term
      cm = costfun
      call smothcost
      cm = costfun - cm

! geostrophic balance term
      cb = 0.0
      if (grdlevl .le. endgslv) then
         cb = costfun
         if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then       ! for sigma and height coordinate
            call gsblncost_non_uniform_z
         elseif (ifpcdnt .eq. 1) then                     ! for pressure coordinate
            call gsblncost_p_hs
         end if
         cb = costfun - cb
      end if

! hydrostatic condition term                 ! added by zhongjie he
      ch = 0.0
      if (grdlevl .le. endhylv) then
         ch = costfun
         ! hydrostatic is based on if rain/snow is analyzed:
         if (numstat .le. 5) call hydrocost_xie
         if (numstat .gt. 5) call hydrocost_shuyuan
         ch = costfun - ch
      end if

      write (*, 1) (cs(s), s=1, numstat + 3)
1     format('each state: ', 10e12.4)
      write (*, 2) cm, cb, ch
2     format('smoothing: ', e12.4, ' geostropic: ', e12.4, ' hydrostatic: ', e12.4)

      return
   end subroutine costfunct2

   subroutine costgradt2
!*************************************************
! calculate gradient of cost function
! history: august 2007, coded by wei li.
! history: january 2008, modified by zhongjie he
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, i1, i2, j1, j2, k1, s, o, no, m, n, uu, vv, zz, rr, rs
      integer  :: np(maxdims)
      real     :: ht, oi, hu, hv, hw, dg, dv, z1, z2
      real     :: ht0, hu0, hv0
      real     :: cc(ngptobs, nallobs)
!-------------------------------------------
! added by shuyuan 20100713  for calculate ref, qr,qs..
      real  ::rour, rc     !   desity*rain water mixing ration,rc=rour
      real  ::rous, sc    !   desity*snow water mixing ration,sc=rous
      real  :: c1, c2    ! coeff  for rour and rous
      real  :: segma     !variance  for ref_obs
      real  ::temp_ref, temp, temp1
      real  :: d2r

      d2r = 3.14159/180.0

!  declaration :
!                'cc' is the coefficent used to calculate gradient of w to control variable
! --------------------
      uu = u_cmpnnt
      vv = v_cmpnnt
      zz = pressure
!added by shuyuan 20100728
      rr = rour_cmpnnt
      rs = rous_cmpnnt
!  call wcompgernl        ! modified by zhongjie he
! initialization

      do s = 1, numstat
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            gradint(i, j, k, t, s) = 0.0
         end do
         end do
         end do
         end do
      end do
      if (nallobs .eq. 0) return
      o = 0
      do s = 1, numstat
         do no = 1, nobstat(s)
            o = o + 1
            do n = 1, maxdims
               np(n) = obsidxpc(n, o)
            end do
            oi = obserror(o)*obserror(o)
            oi = 1.0/oi
            ht = 0.0
            m = 0
            do t = np(4), min0(np(4) + 1, numgrid(4))
            do k = np(3), min0(np(3) + 1, numgrid(3))
            do j = np(2), min0(np(2) + 1, numgrid(2))
            do i = np(1), min0(np(1) + 1, numgrid(1))
               m = m + 1
               ht = ht + obscoeff(m, o)*grdanals(i, j, k, t, s)
            end do
            end do
            end do
            end do
            m = 0
            do t = np(4), min0(np(4) + 1, numgrid(4))
            do k = np(3), min0(np(3) + 1, numgrid(3))
            do j = np(2), min0(np(2) + 1, numgrid(2))
            do i = np(1), min0(np(1) + 1, numgrid(1))
               m = m + 1
               gradint(i, j, k, t, s) = gradint(i, j, k, t, s) &
                                        + (ht - obsvalue(o))*obscoeff(m, o)!!!modified by shuyuan   20101028
               ! +(ht-obsvalue(o))*obscoeff(m,o)/(nobstat(s)*1.0)
            end do
            end do
            end do
            end do
         end do
      end do
      s = numstat + 1                           ! for radial wind observations
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)*obserror(o)
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*www(i, j, k, t)
         end do
         end do
         end do
         end do

         dg = obseinf1(no)
         dv = obseinf2(no)
         ht = hu*sin(d2r*dg)*cos(d2r*dv) + hv*cos(d2r*dg)*cos(d2r*dv) + hw*sin(d2r*dv)
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            gradint(i, j, k, t, uu) = gradint(i, j, k, t, uu) &
                                      + (ht - obsvalue(o))*obscoeff(m, o)*sin(d2r*dg)*cos(d2r*dv)*obsradar
            ! +(ht-obsvalue(o))*oi*obscoeff(m,o)*sin(d2r*dg)*cos(d2r*dv)*obsradar!!!modified by shuyuan   20101028
            gradint(i, j, k, t, vv) = gradint(i, j, k, t, vv) &
                                      + (ht - obsvalue(o))*obscoeff(m, o)*cos(d2r*dg)*cos(d2r*dv)*obsradar
            !+(ht-obsvalue(o))*oi*obscoeff(m,o)*cos(d2r*dg)*cos(d2r*dv)*obsradar!!!modified by shuyuan   20101028
            cc(m, o) = (ht - obsvalue(o))*sin(d2r*dv)*obscoeff(m, o)*obsradar!!!modified by shuyuan   20101028
         end do
         end do
         end do
         end do
      end do

      call wwgradient(cc, s)  ! added by zhongjie he to calculate the gradient of www to u, v and z

      s = numstat + 2                !  for sfmr observation
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)*obserror(o)
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         hu0 = 0.0
         hv0 = 0.0
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*www(i, j, k, t)
            hu0 = hu0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, uu)
            hv0 = hv0 + obscoeff(m, o)*grdbkgnd(i, j, k, t, vv)
         end do
         end do
         end do
         end do
         ht = sqrt((hu + hu0)*(hu + hu0) + (hv + hv0)*(hv + hv0) + hw*hw)
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            gradint(i, j, k, t, uu) = gradint(i, j, k, t, uu) &
                                      + (ht - obsvalue(o))*obscoeff(m, o)*(hu + hu0)/ht
            ! +(ht-obsvalue(o))*oi*obscoeff(m,o)*(hu+hu0)/ht*obs_sfmr!!!modified by shuyuan   20101028
            gradint(i, j, k, t, vv) = gradint(i, j, k, t, vv) &
                                      + (ht - obsvalue(o))*obscoeff(m, o)*(hv + hv0)/ht
            !    +(ht-obsvalue(o))*oi*obscoeff(m,o)*(hv+hv0)/ht*obs_sfmr!!!modified by shuyuan   20101028
            cc(m, o) = (ht - obsvalue(o))*obscoeff(m, o)*hw/ht!*oi*obs_sfmr!!!modified by shuyuan   20101028
         end do
         end do
         end do
         end do
      end do

      call wwgradient(cc, s)
! for radar  reflectivity !  shuyuan 20100721
! --20100907------------------

      ! check if rain analysis is requested:
      if (numstat .le. 5) goto 555

      segma = 1.
      s = numstat + 3
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)*obserror(o)
         oi = 1.0/oi
         temp_ref = 0.0
         m = 0
         rc = 0.
         sc = 0.
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))

            m = m + 1
            rc = rc + obscoeff(m, o)*grdanals(i, j, k, t, rr)
         end do
         end do
         end do
         end do
         m = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            !changed 20100907 shuyuan  rc unit is  g/m3
            !temp_ref=(obsvalue(o)-43.1)/17.5
            !temp_ref=10.**temp_ref
            temp_ref = obsvalue(o)
            temp = (rc - temp_ref)/(segma*segma)
            temp1 = temp*obscoeff(m, o)!!*oi/(nobstat(s)*1.0)  !!!modified by shuyuan   20101028
            gradint(i, j, k, t, rr) = gradint(i, j, k, t, rr) + temp1
         end do
         end do
         end do
         end do
      end do

      ! skip reflectivity:
555   continue

! smooth term
      call smothgrad

! geostrophic balance term
      if (grdlevl .le. endgslv) then
         if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then     ! for sigma and height coordinate
            call gsblngrad_non_uniform_z
         elseif (ifpcdnt .eq. 1) then                   ! for pressure coordinate
            call gsblngrad_p_hs
         end if
      end if

! hydrostatic condition term                 ! added by zhongjie he
      if (grdlevl .le. endhylv) then
         ! hydrostatic is based on if rain/snow is analyzed:
         if (numstat .le. 5) call hydrograd_xie
         if (numstat .gt. 5) call hydrograd_shuyuan
      end if

      return
   end subroutine costgradt2

   subroutine test_fun2
!*************************************************
! show the difference between the analysis and the observations
! history: january 2008, coded by zhongjie he.

!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, o, no, m, n, uu, vv
      integer  :: np(maxdims)
      real  :: ht, oi, hu, hv, hw, dg, dv
      real  :: cs(numstat), cm, cb

      real  :: obs_x, obs_y, obs_z
      real  :: d2r

      d2r = 3.14159/180.0

! --------------------
      uu = u_cmpnnt
      vv = v_cmpnnt

      call wcompgernl        !  modified by zhongjie he.

! initialization
      costfun = 0.0
      do s = 1, numstat
         cs(s) = 0.0
      end do
      if (nallobs .eq. 0) return
      o = 0
      do s = 1, numstat
         do no = 1, nobstat(s)
            o = o + 1
            do n = 1, maxdims
               np(n) = obsidxpc(n, o)
            end do
            oi = obserror(o)**2
            oi = 1.0/oi
            ht = 0.0
            m = 0
            obs_x = 0
            obs_y = 0
            obs_z = 0
            do t = np(4), min0(np(4) + 1, numgrid(4))
            do k = np(3), min0(np(3) + 1, numgrid(3))
            do j = np(2), min0(np(2) + 1, numgrid(2))
            do i = np(1), min0(np(1) + 1, numgrid(1))
               m = m + 1
               ht = ht + obscoeff(m, o)*grdanals(i, j, k, t, s)
               obs_x = obs_x + obscoeff(m, o)*i
               obs_y = obs_y + obscoeff(m, o)*j
               obs_z = obs_z + obscoeff(m, o)*k
            end do
            end do
            end do
            end do
            costfun = costfun + (ht - obsvalue(o))*oi*(ht - obsvalue(o))
            cs(s) = cs(s) + 0.5*(ht - obsvalue(o))*oi*(ht - obsvalue(o))
         end do
      end do
      s = numstat + 1
      do no = 1, nobstat(s)
         o = o + 1
         do n = 1, maxdims
            np(n) = obsidxpc(n, o)
         end do
         oi = obserror(o)**2
         oi = 1.0/oi
         hu = 0.0
         hv = 0.0
         hw = 0.0
         m = 0
         obs_x = 0
         obs_y = 0
         obs_z = 0
         do t = np(4), min0(np(4) + 1, numgrid(4))
         do k = np(3), min0(np(3) + 1, numgrid(3))
         do j = np(2), min0(np(2) + 1, numgrid(2))
         do i = np(1), min0(np(1) + 1, numgrid(1))
            m = m + 1
            hu = hu + obscoeff(m, o)*grdanals(i, j, k, t, uu)
            hv = hv + obscoeff(m, o)*grdanals(i, j, k, t, vv)
            hw = hw + obscoeff(m, o)*www(i, j, k, t)
            obs_x = obs_x + obscoeff(m, o)*i
            obs_y = obs_y + obscoeff(m, o)*j
            obs_z = obs_z + obscoeff(m, o)*k
         end do
         end do
         end do
         end do
         dg = obseinf1(no)
         dv = obseinf2(no)
         ht = hu*sin(d2r*dg)*cos(d2r*dv) + hv*cos(d2r*dg)*cos(d2r*dv) + hw*sin(d2r*dv)
         costfun = costfun + (ht - obsvalue(o))*oi*(ht - obsvalue(o))*obsradar
      end do
      costfun = 0.5d0*costfun

! smooth term
      cm = costfun
      call smothcost
      cm = costfun - cm

! geostrophic balance term
      cb = costfun
      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then      ! for sigma and height coordinate
         call gsblncost_non_uniform_z
      elseif (ifpcdnt .eq. 1) then                    ! for pressure coordinate
         call gsblncost_p_hs
      end if
      cb = costfun - cb
!  write(100,*)(cs(s),s=1,numstat),cm,cb

      return
   end subroutine test_fun2

   subroutine cg_value(f, x, nm)
!*************************************************
! calculate the costfunction, used by minimizer_cg
! history: april 2008, coded by zhongjie he.
!*************************************************
      implicit none
! --------------------
      integer    :: i, j, k, t, s, n, nm
!  double precision :: x(nm),f
      real :: x(nm), f
      n = 0
      do s = 1, numstat
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            n = n + 1
            grdanals(i, j, k, t, s) = x(n)
         end do
         end do
         end do
         end do
      end do

      if (w_cmpnnt .ne. 0) then
         call costfunct1
      else
         call wcompgernl
         call costfunct2
      end if

      f = costfun

   end subroutine cg_value

   subroutine cg_grad(g, x, nm)
!*************************************************
! calculate the gradient, used by minimizer_cg
! history: april 2008, coded by zhongjie he.
!*************************************************
      implicit none
! --------------------
      integer    :: i, j, k, t, s, n, nm
!  double precision :: x(nm),g(nm)
      real :: x(nm), g(nm)

      n = 0
      do s = 1, numstat
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            n = n + 1
            grdanals(i, j, k, t, s) = x(n)
         end do
         end do
         end do
         end do
      end do
      if (w_cmpnnt .ne. 0) then
         call costgradt1
      else
         call wcompgernl
         call costgradt2
      end if
      n = 0
      do s = 1, numstat
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            n = n + 1
            g(n) = gradint(i, j, k, t, s)
         end do
         end do
         end do
         end do
      end do

   end subroutine cg_grad

end module costfun_grad
