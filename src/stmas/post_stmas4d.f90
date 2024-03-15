module post_stmas4d

   use prmtrs_stmas

   public pstprocss, grdmemrls
   private unscaling, retrieval, obsmemrls

!***************************************************
!!comment:
!   this module is used by the module of stmas4d_core to do some post processes.
!   subroutines:
!      pstprocss : the main subroutine of this module.
!      unscaling : unscaling every array.
!      retrieval : retrieval the analysis field for output.
!      obsmemrls : release observations memory
!      grdmemrls : release memory of background and analysis fields

contains

   subroutine pstprocss
!*************************************************
! post process
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
      print *, 'unscaling ........'
      call unscaling
      print *, 'retrieval ........'        ! for laps: add bkg in output by yuanfu
      call retrieval
      print *, 'obs memory releasing .......'
      call obsmemrls
      print *, 'grid memory releasing ........'
      call grdmemrls
      print *, 'end of post process --------'
      return
   end subroutine pstprocss

   subroutine unscaling
!*************************************************
! unscaling every array
! history: august 2007, coded by wei li.
!          december 2008, modified by yuanfu xie:
!                         change if statement:
!                        "if(nallobs.eq.o) return"
!                         to:
!                        "if(nallobs.eq.0) return"
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, o, no, uu
! --------------------
      uu = u_cmpnnt
! unscaling xxx, yyy, zzz or ppp, and cor

      do i = 1, numgrid(1)
      do j = 1, numgrid(2)
         if (pnlt0pu .ge. 1.0e-10 .or. pnlt0pv .ge. 1.0e-10) cor(i, j) = cor(i, j)*scp(csl)
         do k = 1, numgrid(3)
         do t = 1, numgrid(4)
            den(i, j, k, t) = den(i, j, k, t)*scp(dsl)
         end do
         end do
         if (numgrid(1) .ge. 2) xxx(i, j) = oripstn(1) + xxx(i, j)*scp(xsl)
         if (numgrid(2) .ge. 2) yyy(i, j) = oripstn(2) + yyy(i, j)*scp(ysl)
      end do
      end do
      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then      ! for sigma and height coordinate
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            zzz(i, j, k, t) = orivtcl + zzz(i, j, k, t)*scp(psl)
         end do
         end do
         end do
         end do
         close (1)
      elseif (ifpcdnt .eq. 1) then                    ! for pressure coordinate
         do k = 1, numgrid(3)
            ppp(k) = orivtcl + ppp(k)*scp(psl)
         end do
      end if

      if (nallobs .eq. 0) return        ! yuanfu: change o to 0
! unscaling the observations
      o = 0
      do s = 1, numstat
         do no = 1, nobstat(s)
            o = o + 1
            obsvalue(o) = obsvalue(o)*scl(s)
            obserror(o) = obserror(o)*scl(s)
         end do
      end do
      do s = numstat + 1, numstat + 2
         do no = 1, nobstat(s)
            o = o + 1
            obsvalue(o) = obsvalue(o)*scl(uu)
            obserror(o) = obserror(o)*scl(uu)
         end do
      end do
!jhui
      ! only if reflectivity is analyzed:
      if (numstat .le. 5) goto 55
      do s = numstat + 3, numstat + 3
         do no = 1, nobstat(s)
            o = o + 1
            obsvalue(o) = obsvalue(o)*scl(numstat + 1)
            obserror(o) = obserror(o)*scl(numstat + 1)
         end do
      end do
      ! skip rain and snow obs:
55    continue

! unscaling grid point variables
!jhui
      do s = 1, 5
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s)*scl(s)
            grdbkgnd(i, j, k, t, s) = grdbkgnd(i, j, k, t, s)*scl(s)
         end do
         end do
         end do
         end do
      end do
!-----------------------------
!unscale rour and rous
! added by shuyuan 20100818
      ! check if rain and snow are analyzed:
      if (numstat .le. 5) goto 555

      do s = 6, numstat
         do t = 1, numgrid(4)
         do k = 1, numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            grdanals(i, j, k, t, s) = grdanals(i, j, k, t, s)*scl(numstat + 1)
            grdbkgnd(i, j, k, t, s) = grdbkgnd(i, j, k, t, s)*scl(numstat + 1)
         end do
         end do
         end do
         end do
      end do

      ! skip scaling rain and snow:
555   continue

!------------------------------

      print *, 'specific: ', minval(grdanals(1:numgrid(1), 1:numgrid(2), 1:numgrid(3), 2, 5) + &
                                    grdbkgnd(1:numgrid(1), 1:numgrid(2), 1:numgrid(3), 2, 5))

      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         www(i, j, k, t) = www(i, j, k, t)*scl(uu)
      end do
      end do
      end do
      end do

      return
   end subroutine unscaling

   subroutine obsmemrls
!*************************************************
! release observations memory
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
      if (nallobs .ge. 1) then
         deallocate (obspostn)
         deallocate (obscoeff)
         deallocate (obsidxpc)
         deallocate (obsvalue)
         deallocate (obserror)
!    deallocate(obsstate)
         deallocate (obseinf1)
         deallocate (obseinf2)
!    deallocate(obseinf3)
         deallocate (nobstat)
      end if
      deallocate (penal_x)
      deallocate (penal_y)
      deallocate (penal_z)
      deallocate (penal_t)
      deallocate (penal0x)
      deallocate (penal0y)
      deallocate (penal0z)
      deallocate (penal0t)
      deallocate (scl)
      deallocate (sl0)
      return
   end subroutine obsmemrls

   subroutine grdmemrls
!*************************************************
! memory release for grd array
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
      deallocate (www)
      deallocate (cor)
      deallocate (xxx)
      deallocate (yyy)
      deallocate (deg)
      deallocate (den)
      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then
         deallocate (zzz)
      elseif (ifpcdnt .eq. 1) then
         deallocate (ppp)
      end if
      deallocate (grdbkgnd)
      deallocate (grdanals)
      deallocate (gradint)
      return
   end subroutine grdmemrls

   subroutine retrieval
!*************************************************
! retrieval the analysis field for output
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, ln
      character(len=200) :: dr
      integer  :: er
! --------------------
      if (temprtur .gt. 10000) then ! not needed for output in kelvin
         do s = temprtur, temprtur
            do t = 1, maxgrid(4)
            do k = 1, maxgrid(3)
            do j = 1, maxgrid(2)
            do i = 1, maxgrid(1)
               grdbkgd0(i, j, k, t, s) = grdbkgd0(i, j, k, t, s) - 273.15
            end do
            end do
            end do
            end do
         end do
      end if
!  call get_directory('bufr',dr,ln)
!  open(2,file=dr(1:ln)//'difference.dat')
      !open(2,file='./difference.dat')
      !write(2,*) numgrid(1:4),numstat
      !do s=1,numstat
      !  do t=1,numgrid(4)
      !  do k=1,numgrid(3)
      !  do j=1,numgrid(2)
      !  do i=1,numgrid(1)
      !    write(2,*) grdbkgd0(i,j,k,t,s) ! grdanals(i,j,k,t,s)
      !  enddo
      !  enddo
      !  enddo
      !  enddo
      !enddo
      !close(2)
!  open(2,file='background.dat')
!  do s=1,numstat
!    do t=1,numgrid(4)
!    do k=1,numgrid(3)
!    do j=1,numgrid(2)
!    do i=1,numgrid(1)
!      write(2,*)grdbkgd0(i,j,k,t,s)
!    enddo
!    enddo
!    enddo
!    enddo
!  enddo
!  close(2)

!jhui
      do s = 1, numstat
         do t = 1, maxgrid(4)
         do k = 1, maxgrid(3)
         do j = 1, maxgrid(2)
         do i = 1, maxgrid(1)
            grdbkgd0(i, j, k, t, s) = grdanals(i, j, k, t, s)         !+grdbkgd0(i,j,k,t,s)
         end do
         end do
         end do
         end do
      end do

!============= just for output by zhongjie he ==========
      !allocate(difftout(numgrid(1),numgrid(2),numgrid(3),numgrid(4),numstat),stat=er)
      !if(er.ne.0)stop 'difftout allocate wrong'
      !allocate(wwwout(numgrid(1),numgrid(2),numgrid(3),numgrid(4)),stat=er)
      !if(er.ne.0)stop 'wwwout allocate wrong'
      !do t=1,numgrid(4)
      !do k=1,numgrid(3)
      !do j=1,numgrid(2)
      !do i=1,numgrid(1)
      ! wwwout(i,j,k,t)=www(i,j,k,t)
      !  do s=1,numstat
      !    difftout(i,j,k,t,s)=grdanals(i,j,k,t,s)
      !  enddo
      !enddo
      !enddo
      !enddo
      !enddo

!  open(2,file='./wwwout.dat')
!    do t=1,numgrid(4)
!    do k=1,numgrid(3)
!    do j=1,numgrid(2)
!    do i=1,numgrid(1)
!      write(2,*)wwwout(i,j,k,t)
!    enddo
!    enddo
!    enddo
!    enddo
!  close(2)

!=======================================================

      return

   end subroutine retrieval

end module post_stmas4d
