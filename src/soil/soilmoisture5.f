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
      subroutine soil_moisture(imax,jmax,
     &           laps_u,laps_v,laps_t,laps_td,
     &           laps_rain,laps_sc,laps_in,laps_wfz,
     &           laps_mwf,laps_mwf_pre,laps_wx,soiltype,    !inputs down to here
     &           griddry,laps_evap,laps_smc_3d,istatus)

c     laps soilmoisture subroutine. driver is program lsm (laps soil moisture).
c     original version implemented for fsl demonstration pc workstation by
c     chandran subramaniam
c     2/8/93
      
c     john smart 12/1/93: adapt the software to run in real time on the
c     			  unix platform.  set up laps standard i/o.
c     john smart 9/22/97: dynamic array mods

      integer*4 imax,jmax
      include 'soilm.inc'
c
c**** model12 is a soil moisture content model developed in june l986, and
c     based upon the remote sensing based watershed model developed at
c     the university of maryland from 1982 to 1985.  this is version 12.4,
c     as submitted to water resources bulletin in january 1989.
c     created by groves.


      dimension r(100)
      real 	ksat,lamda,in
      real*4    rmiss
      data 	day,sumr,in,oldwea/4*0./
      integer   istatus
      integer   icycle_time
      logical 	dry, griddry
      integer   soiltype(imax,jmax)
      real      laps_u(imax,jmax)     !u-component
      real      laps_v(imax,jmax)     !v-component
      real      laps_t(imax,jmax)     !temperature
      real      laps_td(imax,jmax)    !dewpoint temperature
      real      laps_sc(imax,jmax)    !snow cover - estimate - (0.0 to 1.0)
      real      laps_in(imax,jmax)    !infiltration
      real      laps_wx(imax,jmax)    !weather data
      real      laps_wfz(imax,jmax)   !wetting front depth, z
      real      laps_mwf(imax,jmax)   !wetting front moisture content
      real      laps_mwf_pre(imax,jmax)!previous wetting front moisture content
      real 	laps_rain(imax, jmax) !radar-estimated liq precip
      real      laps_evap(imax,jmax)  !amount of evaporation (calc within)
      real      laps_smc_3d(imax,jmax,3)!three layer soil moisture content

      call get_laps_cycle_time(icycle_time,istat)
      if(istat.ne.1)then
         write(6,*)'error getting laps cycle time'
         write(6,*)'warning: setting delt = 1.0'
         write(6,*)'assuming cycle time = 3600 sec'
         icycle_time = 3600
      endif
      delt = float(icycle_time)/3600.

      call get_r_missing_data(rmiss,istat)
      if(istat.ne.1)goto 1000

      write(6,*)'calculating pan evaporation'
      call calc_evap(imax,jmax,laps_u,laps_v,
     &               laps_t,laps_td,    !input to here.
     &               laps_evap,istatus)
      if (istatus.eq.1) then
 	 write(6,*)'completed calculating evaporation'
      else
         write(6,*)'error computing evaporation, returning'
         istatus = -1
         return
      endif

      write(6,*)' begin soil moisture calculation for grid'

      do j = 1, jmax
      do i = 1, imax

	isoil = soiltype(i,j)
        call soils(isoil,ksat,ths,thr,psif,psiae,lamda)
        call amc(thfc,isoil,ths)

        oldwea = laps_wx(i,j)	     
        in = laps_in(i,j)
        z = laps_wfz(i,j)
	thi = laps_mwf(i,j)
        xthi = laps_mwf_pre(i,j)

        if(in .ne.rmiss.and.  z .ne.rmiss.and.
     &     thi.ne.rmiss.and.xthi.ne.rmiss)then

        if((laps_sc(i,j).lt.snowthres).or.(laps_sc(i,j).gt.1.e30))then

          dry = .false.
          if((griddry).or.(laps_rain(i,j).lt.rainthres))dry=.true.
          if (dry) then
	    pan = laps_evap(i,j)
c	    if(oldwea.le.0.)xthi=thi     ! i'm not sure  about this statement
            call moist(xthi,in,xday,ksat,ths,thr,psiae,lamda,thi,
     1                 depth,pan,cumd,cumet,bal,z,rzst,oldwea)
	    call split(z,thi,xthi,thi1,thi2,thi3,hor1,hor2,hor3)
            oldwea = 1.0/24.
            in = 0.0
          else        ! rain 
            if(z.lt.dcm)xthi=(thi*z+xthi*(dcm-z))/dcm
            r(1) = laps_rain(i,j)
	    call entry(cumq,in,xthi,ksat,ths,psif,delt,r,depth,n)
	    z=in/(ths-xthi)
	    call split(z,ths,xthi,thi1,thi2,thi3,hor1,hor2,hor3)
            oldwea = -1./24
            thi = ths
          endif
          laps_mwf(i,j) = thi
          laps_mwf_pre(i,j) = xthi
          laps_wfz(i,j)= z
  	  laps_in(i,j) = in
	  laps_wx(i,j) = oldwea
          if(ths.eq.0.0)then
             laps_smc_3d(i,j,1) = 0.0
             laps_smc_3d(i,j,2) = 0.0
             laps_smc_3d(i,j,3) = 0.0
          else
	     thsper = ths / 100.0 
	     laps_smc_3d(i,j,1) = thi1 ! / thsper
	     laps_smc_3d(i,j,2) = thi2 ! / thsper
	     laps_smc_3d(i,j,3) = thi3 ! / thsper

          endif

        else

          laps_mwf(i,j)=ths
          laps_mwf_pre(i,j)=xthi
	  laps_wfz(i,j) = hor1*2.54
          laps_in(i,j)=0.1  !not sure about this
	  laps_wx(i,j) = 1./24. 

          if(ths.eq.0.0)then
             laps_smc_3d(i,j,1) = 0.0
             laps_smc_3d(i,j,2) = 0.0
             laps_smc_3d(i,j,3) = 0.0
          else
	     thsper = ths / 100.0 
	     laps_smc_3d(i,j,1) = ths ! / thsper
	     laps_smc_3d(i,j,2) = ths ! / thsper
	     laps_smc_3d(i,j,3) = ths ! / thsper
          endif

        endif

        endif

      enddo
      enddo	
      write(6,*)' completed sm calculation for grid'
      istatus = 1
c
100   write(6,*)'end of simulation'
      goto 1000
c
1000  return
      end

c ====================================================================

      subroutine soils(isoil,ksat,ths,thr,psif,psiae,lamda)
c**** subroutine provides default values for soil hydraulic parameters
      dimension satcon(6),respor(6),effpor(6),bubpr(6)
      real ksat,lam(6),lamda
      data satcon/10.,3.5,.65,1.3,.08,.03/
      data respor/.04,.05,.08,.10,.08,.11/
      data effpor/.39,.4,.39,.32,.4,.38/
      data bubpr/10.,15.,27.,15.,80.,130./
      data lam/.43,.38,.31,.23,.23,.20/
      the=effpor(isoil)
      ksat=satcon(isoil)/60.
      thr=respor(isoil)
      psibub=bubpr(isoil)
      lamda=lam(isoil)
      ths=thr+the
      psiae=psibub/2.
      psif=psiae*(2.+3.*lamda)/(1.+3.*lamda)
      return
      end
c
      subroutine amc(thi,isoil,ths)
c**** subroutine assigns initial soil moisture value
      dimension thfc(6),thwilt(6)
      data thwilt/.06,.075,.13,.16,.22,.27/
      data thfc/.125,.175,.25,.23,.37,.425/
      thi=thfc(isoil)
      if(thi.lt.thwilt(isoil))xthi=thwilt(isoil)
      if(thi.gt.ths)xthi=ths
      return
      end
c
      subroutine entry(ro,in,thi,ksat,ths,psif,delt,r,idepth,n)
c**** subroutine calculates infiltration volume for rainfall period
      dimension cumf(100),r(100),cumq(100),q(100)
      real ksat,in
      integer idepth
      if(thi.gt.ths)thi=ths
      smax=(ths-thi)*float(idepth)*2.54
      cumf(1)=in
c      cumf(1)=0.
      cumq(1)=0.
      fc=ksat*delt*60.
      do 10  it=2,n+1
	if(cumf(it-1).gt.smax) go to 2
	call infil(thi,delt,ksat,ths,psif,cumf(it-1),deli)
	go to 4
2	deli=fc
4	if(deli.gt.r(it-1)) deli=r(it-1)
	cumf(it)=cumf(it-1)+deli
	q(it)=r(it-1)-deli
	cumq(it)=cumq(it-1)+q(it)
	ro=cumq(it)
10	in=cumf(it)
      return
      end
c
      subroutine infil(thi,delt,ksat,ths,psif,in,deli)
c**** subroutine calculates incremental infiltration volume
      real ksat,in,itx
      dt=delt*60.
      bb=psif*(ths-thi)
      if(in.le.0.)in=0.01
      xe=exp(-in/bb)
      ti=(in-bb+bb*xe)/ksat
      tx=ti+dt
      itx=in
      do 10 m=1,5
	if(itx .gt. 5)goto 20
	xe=exp(-itx/bb)
10	itx=itx-(itx-bb+bb*xe-ksat*tx)/(1.-xe)
20    deli=itx-in
      if(in.eq..01)in=0.
      return
      end
c
      subroutine moist(xthi,xin,xday,ksat,ts,tr,ps,la,
     1theta,ird,pr,cumd,cumtr,bal,z,rzst,oldwea)
c**** subroutine simulates redistribution of soil moisture between storms
      common/com1/ptr,dtheta,thetaa,thetab
      common/com2/ks,ths,thr,psiae,lamda,rzd,t
      real ksat,ks,lamda,la
      data psia,psib/15500.,5166./
      ks=ksat/60.
      ths=ts
      thr=tr
      psiae=ps
      lamda=la
      rzd=ird*2.54
c      if(rzd.lt.50.)rzd=50.
      ptr=pr*2.54/86400.
c**** change time units to seconds
      tstart=0.
      t=tstart*86400.
      tout=xday*86400.
c**** initialize theta, z, storages
      if(oldwea.lt.0.)then
	 theta=ths
	 z=xin/(theta-xthi)
	 oldrzs=theta*rzd
	 if(z.lt.rzd)oldrzs=z*theta+(rzd-z)*xthi
      else
	 oldrzs=rzst
      endif
      cumd=0.
      cumtr=0.
c**** determine et thresholds
      psibub=psiae*2.
      the=ths-thr
      thetaa=the*(psibub/psia)**lamda+thr
c ******************************************
c  hard wire minimum theta to be lower here. this sets it
c    to the wilting point.
c  try 20% of ths first
      thetaa=.2*ths
c ******************************************
      thetab=the*(psibub/psib)**lamda+thr
c**** advance soil moisture from tstart to xday
      call intgrl(theta,xthi,z,cumd,cumtr,tout)
c**** calculate storages
      rzst=theta*rzd
      if(z.lt.rzd)rzst=z*theta+(rzd-z)*xthi
      bal=-cumtr-cumd-(rzst-oldrzs)
c     if(z.eq.0.)z=rzd
      return
      end
c
      subroutine intgrl(theta,xthi,z,cumd,cumtr,tout)
c**** subroutine integrates soil moisture variables
      common/com1/ptr,dtheta,thetaa,thetab
      common/com2/ks,ths,thr,psiae,lamda,rzd,t
      dimension x(5),dx(5),dx1(5),x1(5)
      real ks,lamda
      x(1)=theta
      x(2)=xthi
      x(3)=z
      x(4)=cumd
      x(5)=cumtr
      dtheta=1.
      dx1(1)=1.
      go to 10
1     continue
c**** check for small change in theta at wetting front
      if(z.eq.0.)go to 10
      if(dtheta.gt..03.and.z.lt.(.9*rzd))go to 10
      if(z.lt.rzd)x(1)=x(2)+x(3)*(x(1)-x(2))/rzd
      x(2)=x(1)
      x(3)=0.
c**** calculate dt
10    dt=amin1(-.004/dx1(1),tout-t)
      if(dt.lt.1.)dt=1.
      dt2=dt/2.
c**** integration routine
      call rate (x(1),x(2),x(3),dx1(1),dx1(2),dx1(3),dx1(4),dx1(5))
      do 20 i=1,5
20    x1(i)=x(i)+dt2*dx1(i)
      t=t+dt
      call rate (x(1),x(2),x(3),dx(1),dx(2),dx(3),dx(4),dx(5))
      do 30 i=1,5
30    x(i)=x(i)+dt2*(dx1(i)+dx(i))
      theta=x(1)
      xthi=x(2)
      z=x(3)
      cumd=x(4)
      cumtr=x(5)
      if(t.lt.tout)go to 1
      return
      end
c
      subroutine rate(theta,xthi,z,dthdt,dthidt,dzdt,drn,tr)
c**** subroutine calculates rates of change of soil moisture redistribution
      common/com1/ptr,dtheta,thetaa,thetab
      common/com2/ks,ths,thr,psiae,lamda,rzd,t
      real ks,lamda,lam23
      the=ths-thr
      c1=psiae/(1.+3.*lamda)
      lam23=3.+2./lamda
      tr1=0.
      if(z.eq.0.)go to 10
c**** determine rates above the wetting front
      dtheta=theta-xthi
      if (dtheta .eq. 0) then
        z = 0
	go to 10
      endif
      qz=ks*(c1*((theta-thr)/the)**(3.+1./lamda)/z+((theta-thr)
     1/the)**lam23)
      dzdt=qz/dtheta
      f1=1.
      if(theta.lt.thetab)f1=(theta-thetaa)/(thetab-thetaa)
      if(theta.lt.thetaa)f1=0.
      tr1=f1*ptr*z/rzd
      dthdt=-((qz)+tr1)/z
c**** determine et below the wetting front
10    f1=1.
      if(xthi.lt.thetab)f1=(xthi-thetaa)/(thetab-thetaa)
      if(xthi.lt.thetaa)f1=0.
      tr2=f1*ptr*(rzd-z)/rzd
c**** determine moisture loss rates
      tr=tr1+tr2
      drn=ks*((xthi-thr)/the)**lam23
      dthidt=-(drn+tr2)/(rzd-z)
      if(z.ne.0.)return
      dzdt=0.
      dthdt=dthidt
      return
      end
c
      subroutine split(z,thi,xthi,thi1,thi2,thi3,hor1,hor2,hor3)
c**** subroutine calculates moisture content for three soil horizons
      h1=hor1*2.54
      h2=hor2*2.54
      h3=hor3*2.54
      if(z.le.h1)then
	thi1=(z*thi+(h1-z)*xthi)/h1
	thi2=xthi
	thi3=xthi
	return
      else if (z.le.h2)then
	thi1=thi
	thi2=((z-h1)*thi+(h2-z)*xthi)/(h2-h1)
	thi3=xthi
	return
      else if(z.le.h3)then
	thi1=thi
	thi2=thi
	thi3=((z-h2)*thi+(h3-z)*xthi)/(h3-h2)
	return
      else
	thi1=thi
	thi2=thi
	thi3=thi
      endif
      return
      end	
