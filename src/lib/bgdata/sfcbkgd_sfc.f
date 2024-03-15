      subroutine sfcbkgd_sfc(bgmodel,t,q,height,heightsfc
     &,td_sfc,tsfc,sh_sfc, ter,p,imx,jmx,kx,psfc,nx_pr,ny_pr)

c inputs are from laps analyzed 2- 3d fields
c     input       t        analyzed temperature     3d
c                 q        analyzed mixxing ratio   3d
c                 height   analyzed height          3d
c                 heightsfc analyzed height         2d
c                 ter      terrain                  2d
c                 imx      dimension e-w
c                 jmx      dimension n-s
c                 kx       number of vertical levels
c                 p        laps pressure levels     1d
c                 tsfc     !i/o input model t, output recomputed t
c                 sh_sfc   !i,   model sh
c                 td_sfc   !i/o, model td recomputed in routine
c                 qsfc     !i/o surface spec hum, input as q or computed internally
c
c recompute surface variables on hi-res terrain using 2d model surface data. 

c
c j. smart    09-22-98:	original version: this is used to compute
c                       sfc p for lgb when using nogaps1.0 deg since
c                       this field currently does not come with the
c                       model grids at afwa.
c    "        02-01-99: recompute sh_sfc with new psfc and tsfc for
c                       consistency.
c    "        11-18-99: put psfc,tsfc, and sh_sfc comps into subroutine
c kml:  changes made april 2004
c this subroutine now reads in heightsfc and td_sfc
c heightsfc and td_sfc are passed into subroutine compute_sfc_bgfields
c kml: end
c
      implicit none

      integer    i
      integer    imx
      integer    j
      integer    jmx
      integer    k
      integer    kx
      integer    it
      integer    bgmodel
      integer    nx_pr,ny_pr,ip,jp
      integer    istatus

      logical    lfndz

      real       height      ( imx , jmx, kx )
      real       heightsfc   ( imx, jmx )
      real       p           ( nx_pr,ny_pr,kx )
      real       psfc        ( imx , jmx )               ! i    (mb)
      real       tsfc        ( imx , jmx )               ! i/o  (k)
      real       q           ( imx , jmx , kx )
      real       rh3d        ( imx , jmx , kx )
      real       rh2d        ( imx , jmx )
      real       sh_sfc      ( imx,  jmx )               ! i    (kg/kg)
      real       td_sfc      ( imx,  jmx )               ! i/o  (k)
      real       t           ( imx , jmx , kx )
      real       esat
      real       ter         ( imx , jmx )
      real       make_td
      real       make_rh
      real       make_ssh
      real       t_ref,badflag
      real       p_mb
      real       r_missing_data
c
c if bgmodel = 6 or 8 then sh_sfc is indeed td sfc
c if bgmodel = 4      then sh_sfc is rh (wfo - ruc)
c if bgmodel = 9      then no surface fields input. compute all from 3d
c                     fields. q3d used. (nos - eta)
c otherwise sh_sfc is input with q
c 
      write(6,*)' subroutine sfcbkgd_sfc, bgmodel = ',bgmodel

      call get_r_missing_data(r_missing_data,istatus)

      t_ref=-132.0

      if(minval(td_sfc) .eq. r_missing_data .and.
     &   maxval(td_sfc) .eq. r_missing_data       )then
         write(6,*)' note: td_sfc has missing data values'
         write(6,*)' computing td_sfc from sh_sfc'

          do i = 1,imx
          do j = 1,jmx
              td_sfc(i,j)=make_td(psfc(i,j)/100.,tsfc(i,j)-273.15
     1                           ,sh_sfc(i,j)*1000.
     1                           ,t_ref)+273.15
          enddo ! j
          enddo ! i
      endif

      if(bgmodel.eq.3.or.bgmodel.eq.9)then
         write(6,*)' bgmodel is ',bgmodel,' convert 3d q to rh'
         do k=1,kx
            do j=1,jmx
            do i=1,imx
               ip = min(i,nx_pr)
               jp = min(j,ny_pr)
               p_mb = p(ip,jp,k) / 100.
               rh3d(i,j,k)=make_rh(p_mb,t(i,j,k)-273.15,q(i,j,k)*1000.
     +,t_ref)*100.
            enddo
            enddo
         enddo

         badflag=0.
         write(6,*)
     &   ' interp 3d t and rh to hi-res terrain, set sh_sfc array to rh'
         call interp_to_sfc(ter,rh3d,height,imx,jmx,kx,
     &                      badflag,sh_sfc) ! here sh_sfc temporarily is rh
         call interp_to_sfc(ter,t,height,imx,jmx,kx,badflag,
     &                      tsfc)
      endif

      write(6,*)' compute sfc fields using 2d model surface data'

      do j=1,jmx
      do i=1,imx

         ip = min(i,nx_pr)
         jp = min(j,ny_pr)

         k=0
         lfndz=.false.
         do while(.not.lfndz)
            k=k+1
            if(height(i,j,k).gt.ter(i,j))then
               lfndz=.true.
c
c
c recompute psfc with moisture consideration. snook's orig code.
c              it=sh_sfc(i,j)*100.         
c              it=min(45000,max(15000,it))
c              xe=esat(it)
c              qsfc=0.622*xe/(p_sfc-xe)
c              qsfc=qsfc/(1.+qsfc)
c --------------------------------------------------

               if(k.gt.0)then
                  call compute_sfc_bgfields_sfc(bgmodel,imx,jmx,kx,i,j,k
     &,ter(i,j),height,heightsfc(i,j),t,p,q,t_ref,psfc(i,j),tsfc(i,j),
     &sh_sfc(i,j),td_sfc(i,j),ip,jp,nx_pr,ny_pr)
c              else
c                 tsfc(i,j)=t(i,j,k)
c                 td_sfc(i,j)=make_td(p(k),t(i,j,k)-273.15,q(i,j,k)*1000.
c    &,t_ref)+273.15

               endif

            endif
            if(k.eq.kx)lfndz=.true.
         enddo
      enddo
      enddo 

      if(minval(sh_sfc) .eq. r_missing_data .and.
     &   maxval(sh_sfc) .eq. r_missing_data       )then
         write(6,*)' note: sh_sfc has missing data values'
      endif

      if(minval(td_sfc) .eq. r_missing_data .and.
     &   maxval(td_sfc) .eq. r_missing_data       )then
         write(6,*)' warning: td_sfc has missing data values'
      endif

      write(6,*)' returning from subroutine sfcbkgd_sfc'
      
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine compute_sfc_bgfields_sfc(bgm,nx,ny,nz,i,j,k,ter,height
     &,heightsfc,t,p,q,t_ref,psfc,tsfc,sh_sfc,td_sfc,ip,jp,nx_pr,ny_pr)
c
c j. smart 11/18/99 put existing code in subroutine for other process use in laps
c
c kml: changes made april 2004
c heightsfc (analyzed 2m height) and td_sfc are now read in
c dz is computed as the difference between heightsfc and terrain
c lapse rate adjustment is made to tsfc and td_sfc (2m level)...
c  ...rather than t and td2 (pressure levels)
c kml: end
      implicit none
                                
      integer nx,ny,nz
      integer i,j,k            !i, i,j,k coordinate of point for calculation
      integer bgm              !i, model type {if = 0, then sh_sfc input = qsfc}
      integer nx_pr,ny_pr,ip,jp
      integer istatus
                                                      
      real   p(nx_pr,ny_pr,nz) !i, pressure of levels (pa)
      real   ter               !i, terrain height at i,j
      real   t_ref             !i, reference temp for library moisture conv routines
                                                     
      real   height(nx,ny,nz)  !i, heights of pressure levels 3d
      real   heightsfc         !i, height of surface at i,j
      real   t(nx,ny,nz)       !i, temperatures 3d
      real   q(nx,ny,nz)       !i, specific humidity 3d
      real   psfc              !i/o, input bkgd model sfc p; output recalculated surface pressure, pa
      real   tsfc              !i/o input model t, output recomputed t
      real   sh_sfc            !i/o input model sh
      real   td_sfc            !i/o input model td, output recomputed td
      real   qsfc              !i/o surface spec hum, input as q or computed internally
                                                           
      real   tbar
      real   td1,td2
      real   p_mb,p_mb_p1,p_mb_m1                                 
      real   pdry   ! bls/iteris - needed for estimate of sfc pressure before moisture adjustment
      real   g,r
      real   ssh2,make_ssh,make_td
      real   dz,dzp,dtdz
      real   tvsfc,tvk,tbarv
      real   r_missing_data
                                                        
      parameter (g         = 9.8,
     &           r         = 287.04)
c
c if first guess values are missing data then return missing data
c
       call get_r_missing_data(r_missing_data,istatus)
       if(tsfc.lt.500.0.and.t(i,j,k).lt.500.0.and.
     .    td_sfc.lt.500.0)then
c
c first guess psfc without moisture consideration
c

          p_mb = p(ip,jp,k) / 100.
                                              
          tbar=(tsfc+t(i,j,k))*0.5
          dz=heightsfc-ter
c         psfc=p_mb*exp(g/(r*tbar)*dz)
          psfc=psfc/100.                                ! calcs done in mb

          ! bls/iteris - original code modified psfc directly, only to 
          ! do it again farther below, which was not correct and leading
          ! to wacky surface pressures in area of big dz
          ! orig:  psfc=psfc*exp(g/(r*tbar)*dz)
          pdry=psfc*exp(g/(r*tbar)*dz)
          
!         determine qsfc
          ! bls/iteris - where pdry is used below used to be psfc (incorrect)
          if(bgm.eq.0.or.bgm.eq.6.or.bgm.eq.8)then      ! sh_sfc is td
             qsfc=ssh2(pdry,tsfc-273.15
     &                  ,sh_sfc-273.15,t_ref)*.001      ! sh is kg/kg
          elseif(bgm.eq.3.or.bgm.eq.4.or.bgm.eq.9)then  ! sh_sfc is rh
             qsfc=make_ssh(pdry,tsfc-273.15
     &                      ,sh_sfc/100.,t_ref)*.001    ! sh is kg/kg
                                           
          else                                          ! sh_sfc is q (kg/kg)
             qsfc=sh_sfc
          endif
                                         
c pressure
          tvsfc=tsfc*(1.+0.608*qsfc)
          tvk=t(i,j,k)*(1.+0.608*q(i,j,k))
          tbarv=(tvsfc+tvk)*.5
c         psfc=(p_mb*exp(g/(r*tbarv)*dz))*100.  !return units = pa
          psfc=(psfc*exp(g/(r*tbarv)*dz))*100.  !return units = pa
          if(k.gt.1)then
             dzp=height(i,j,k)-height(i,j,k-1)
c temp
             dtdz=(t(i,j,k-1)-t(i,j,k))/dzp
             tsfc=tsfc+dtdz*dz
c dew point temp
             td2=make_td(p_mb,t(i,j,k)-273.15,q(i,j,k)*1000.,t_ref)
             p_mb_m1 = p(ip,jp,k-1) / 100.
             td1=make_td(p_mb_m1,t(i,j,k-1)-273.15,q(i,j,k-1)*1000.
     .,t_ref)
             td_sfc=td_sfc+((td1-td2)/dzp)*dz

             td_sfc = min(tsfc,td_sfc)
c
          else
             dzp=height(i,j,k+1)-height(i,j,k)
             p_mb_p1 = p(ip,jp,k+1) / 100.
             td2=make_td(p_mb_p1,t(i,j,k+1)-273.15
     &,q(i,j,k+1)*1000.,t_ref)+273.15
             td1=make_td(p_mb,t(i,j,k)-273.15
     &,q(i,j,k)*1000.,t_ref)+273.15
             dtdz=(t(i,j,k)-t(i,j,k+1))/dzp
             tsfc=tsfc+dtdz*dz
             tbar=(tsfc+t(i,j,k))*0.5
c            psfc=p_mb*exp(g/(r*tbar)*dz)
             psfc=psfc*exp(g/(r*tbar)*dz)
             td_sfc=td_sfc+((td1-td2)/dzp)*dz

             td_sfc = min(tsfc,td_sfc)
                                   
          endif
       else
          psfc =r_missing_data
          tsfc =r_missing_data
          td_sfc=r_missing_data
       endif
       return
       end
