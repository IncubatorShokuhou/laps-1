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
      subroutine dprep_eta_conusc(nx,ny,nz,
     +        ht, p, tp, uw, vw, rh, ht_sfc, p_sfc,
     +        rh_sfc, tp_sfc, uw_sfc, vw_sfc, mslp,istatus)
      implicit none
      real*4    cp,kappa
      parameter(cp=1004.,kappa=287./1004.)

      real*4 esat,es
      common /estab/esat(15000:45000),es(15000:45000)
      integer nx,ny,nz, i, j, k, it, istatus
      real xe, mrsat
      
      real mslp( nx,  ny),    ht( nx,  ny,  nz), 
     +     ht_sfc( nx,  ny),  p(nx,ny,nz), p_sfc( nx,  ny), 
     +     rh( nx,  ny,  nz), rh_sfc( nx,  ny), 
     +     tp( nx,  ny,  nz), tp_sfc( nx,  ny), 
     +     uw( nx,  ny,  nz), uw_sfc( nx,  ny), 
     +     vw( nx,  ny,  nz), vw_sfc( nx,  ny)

     
      real tmp(nx,ny,nz), factor,
     +     th( nx,  ny,  nz),  th_sfc( nx,  ny)

c *** convert surface and 3d temp to theta.
c *** compute exner function.
c *** convert sfc pressure and mean sea level pressure from pa to mb.
c *** convert surface and 3d rh (%) to mr.
c
      
c
c_______________________________________________________________________________
c  p is input as a single verticle column output as 3d array
c
      do k=1,nz
         p(1,1,k) = p(k,1,1)
      enddo
      
      do k=1,nz
         do j=1,ny
            do i=1,nx
               p(i,j,k)=p(1,1,k)
            enddo
         enddo
      enddo

      
      do k=1,nz
         do j=1,ny
            do i=1,nx


               it=tp(i,j,k)*100
               it=min(45000,max(15000,it))
               xe=esat(it)
               mrsat=0.00622*xe/(p(i,j,k)-xe) !assumes units of rh on next line is %
               rh(i,j,k)=rh(i,j,k)*mrsat


               factor=(1000./p(i,j,k))**kappa
               tp(i,j,k) = tp(i,j,k)*factor
               p(i,j,k) = cp/factor

            enddo
         enddo
      enddo
      print*,(p(2,2,k),k=1,nz)


      do j=1,ny
         do i=1,nx
            p_sfc(i,j)=p_sfc(i,j)*0.01
            mslp(i,j)=mslp(i,j)*0.01
            th_sfc(i,j)=tp_sfc(i,j)*(1000./p_sfc(i,j))**kappa
            it=tp_sfc(i,j)*100
            it=min(45000,max(15000,it))
            xe=esat(it)
            mrsat=0.00622*xe/(p_sfc(i,j)-xe) !assumes units of rh on next line is %
            rh_sfc(i,j)=rh_sfc(i,j)*mrsat
         enddo
      enddo

      print *,'sfc:',p_sfc(31,31),ht_sfc(31,31),tp_sfc(31,31)
     .   ,rh_sfc(31,31),uw_sfc(31,31),vw_sfc(31,31),mslp(31,31)
      print *,'ua :',ht(nx,ny,nz),tp(nx,ny,nz),rh(nx,ny,nz),
     .        uw(nx,ny,nz),vw(nx,ny,nz)
      print *,'ua :',ht(nx,ny,1),tp(nx,ny,1),rh(nx,ny,1),
     .        uw(nx,ny,1),vw(nx,ny,1)

      return
      end

