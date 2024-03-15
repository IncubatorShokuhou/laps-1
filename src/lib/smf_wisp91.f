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

      subroutine get_smf_1d(nk,cbase_m,ctop_m,itype_in
     1  ,height_laps_1d,pres_laps_1d,temp_laps_1d,
     1                          rlwc_laps,prob_laps,mode)

!     1991 steve albers (smith-feddes code adapted to laps)

!     the array pres_laps_1d begins at 900mb and decreases at 50mb steps.

!     real rlwc_laps(nk),prob_laps(nk)
!     real height_laps_1d(nk),pres_laps_1d(nk),temp_laps_1d(nk)

!     inputs
      integer nk              ! number of laps vertical levels
      real cbase_m            ! cloud base (meters msl)
      real ctop_m             ! cloud top (meters msl)
      integer itype_in        ! cloud type (now hardwired to 1) (stratus)
      real height_laps_1d(nk) ! vert array, heights (meters msl)
      real pres_laps_1d(nk)   ! vert array, laps pressure levels (mb)
      real temp_laps_1d(nk)   ! vert array, ambient temperature (deg k)
      integer mode            ! not used (can be removed)


!     outputs
      real rlwc_laps(nk)      ! vert array, lwc (g/m**3)
                                ! should be assigned for all laps levels within
                                ! the cloud layer.
      real prob_laps(nk)      ! vert array, not used right now.

      common /pant/ bot(4),top(4),prb(4),prt(4),btem(4),ttem(4)
      common /qant/ nlvl,heitl(260),pres(260),temp(260)
      common /iodev/k5,k6,k7,k8
c
c     dimension watr(4,2,260),xmvd(4),icld(8),ictp(4),icp(4)
      dimension watr(4,2,260),xmvd(4),ictp(4),icp(4)
      character*70 title
      character*2 cldtyp(4)
c     character*1 ans
      character*11 aindex(4,260),hindex(4,260)

      common /prof/ nobhgt,height(200),pressp(200),tempro(200)
c
c-------define the input/output logical devices
c
      k5 = 5
      k6 = 6
      k7 = 7
      k8 = 8

      bot(1) = cbase_m    ! bottom height m
      top(1) = ctop_m    ! top height m

      if(mode .eq. 1)then ! laps data is passed in through the call
          nobhgt = nk
          do i = 1,nobhgt
              height(i) = height_laps_1d(i)
              pressp(i) = pres_laps_1d(i)
              tempro(i) = temp_laps_1d(i)
          enddo ! i

      elseif(mode .eq. 2)then ! test mode
!         read in the test sounding
          call raob(nxflag,nvflag,maxobsr) ! this routine was modified for profs

      endif

      ndecks = 1

!     0        cldtyp(i) = '  '
!     1        cldtyp(i) = 'st'
!     2        cldtyp(i) = 'sc'
!     3        cldtyp(i) = 'cu'
!     4        cldtyp(i) = 'ns'
!     5        cldtyp(i) = 'ac'
!     6        cldtyp(i) = 'as'
!     7        cldtyp(i) = 'cs'
!     8        cldtyp(i) = 'ci'
!     9        cldtyp(i) = 'cc'
!    10        cldtyp(i) = 'cb'

      ictp(1) = itype_in     ! type of layer 1
      icp(1)  = 100   ! percent cloud cover of layer 1

      if(top(1) .le. bot(1))return

      call interp_smf(ndecks)
      call lwc(ndecks,icp,ictp,watr,xmvd) ! add this to link into .new routines

!     call a modified version of the output routine
c     call output(cldtyp,watr,xmvd,aindex,hindex,icp,ndecks,title)

      call out_laps
     1  (nk,cldtyp,watr,xmvd,aindex,hindex,icp,ndecks,title,
     1         rlwc_laps,prob_laps)

!     write(6,*)' height  lwc  %liquid'
      i = 1
      n1 = watr(i,1,1)
        iscript = 0
        do k = 3,iscript ! 260
            r_height = (k+n1-2)*100.
            if(watr(i,1,k) .gt. 0.)
     1  write(6,101)r_height,(watr(i,j,k),j=1,2)
101         format(1x,f8.1,'m ',2f7.3)
        enddo ! k

      return

      end
      subroutine raob(nxflag,nvflag,maxobsr)
c***********************************************************************
c                          subroutine raob
c***********************************************************************
c<begin>
c<identification>          name:  raob
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  usrint
c=======================================================================
c<description>
c    this subroutine reads in the vertical profile or upperair data
c    needed to run the smith-feddes model.
c=======================================================================
c<called routines>
c    hsort - (subroutine) eliminates redundant height levels and sorts
c            the raob data in ascending order of height
c    inkey - (subroutine) gets a character from the keyboard without
c            echo.  the routine is not available in source code form.
c=======================================================================
c<parameters>
c    formal declaration:
c       call raob(nxflag,nvflag,maxobsr)
c    input:
c       maxobsr  - (integer) the maximum number of upperair levels
c                  allowed
c    output:
c       nxflag   - (integer) a flag
c                  0 - the option to compute the icing index was not
c                      selected
c                  1 - the option to compute the icing index was
c                      selected
c       nvflag   - (integer) a flag
c                  0 - the vertical profile data is not yet available
c                      for use
c                  1 - the vertical profile data is available for use
c    output (in common):
c       nobhgt   - (integer) the number of levels in the temperature-
c                  pressure-height vertical profile
c       height   - (real) the height of each level in the profile
c       pressp   - (real) the pressure profile
c       tempro   - (real) the temperature profile
c    common:
c       k5       - (integer) the keyboard input unit
c       k6       - (integer) the screen output unit
c       k7       - (integer) the disk file unit
c=======================================================================
c<history>
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              original code.
c    09/14/89  asl    (505) 678-1570    elton p. avara
c              created subroutine hsort from part of this code
c    12/21/89  cut out screen i/o for profs version, all this does now
c              is essentially reads the sounding from a file.
c=======================================================================
c<end>
c***********************************************************************
c
c-------set the minimums and the maximums of the different input
c       parameters so there are no magic numbers floating around.
c
      real tempmin,tempmax,presmin,presmax,heitmin,heitmax
      parameter (tempmin  =   150.000)
      parameter (tempmax  =   349.999)
      parameter (presmin  =    50.000)
      parameter (presmax  =  1099.999)
      parameter (heitmin  =    50.000)
      parameter (heitmax  = 24999.999)
c
      common /prof/ nobhgt,height(200),pressp(200),tempro(200)
      common /iodev/k5,k6,k7,k8
c
      character*39 datafile
c     character*10 scrclr,answer
c     character*8 curpos
c     character*4 scrblk,hunits(2),tunits(3),punits(2)
      character*4 hunits(2),tunits(3),punits(2)
c     character*3 lerase
c     character*1 num(10),ans,ans2
      character*1 num(10)
c
      data num/'0','1','2','3','4','5','6','7','8','9'/
      data hunits/'(m) ','(ft)'/
      data tunits/'(k) ','(c) ','(f) '/
      data punits/'(mb)','(in)'/
      data nheit,ntemp,npres/3*1/
c
        datafile='sounding.dat'
        open(k7,file=datafile,err=19,status='old')
        read(k7,*) nobhgt
        do 21 j=1,nobhgt
   21     read(k7,*) height(j),tempro(j),pressp(j)
        close(k7)
c
c-------get the desired units for height, temperature, and pressure.
c
      nheit=1
      ntemp=1
      npres=1
c
c-------if vertical profile data are available from a file, convert the
c       data to the desired units.
c
   50 heitmn=heitmin
      heitmx=heitmax
c
      tempmn=tempmin
      tempmx=tempmax
c
      presmn=presmin
      presmx=presmax
c
c-------initialize the cursor position and the icing index computation
c       flag
c
      item=0
      nxflag=0
c
c-------set the vertical profile flag and return.
c
      nvflag=1
c
19    return
      end
      subroutine interp_smf(ndecks)
c***********************************************************************
c                          subroutine interp
c***********************************************************************
c<begin>
c<identification>          name:  interp
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  input
c=======================================================================
c<description>
c    this subroutine takes the temperature (tempro) and pressure
c    (pressp) profiles and interpolates between them to yield the
c    temperatures at the geometric layer bottoms (btem), 100 m height
c    levels (temp), and geometric layer tops (ttem); and the pressures
c    at the geometric layer bottoms (prb), 100 m height levels (pres),
c    and geometric layer tops (prt).
c=======================================================================
c<called routines>
c    none
c=======================================================================
c<parameters>
c    formal declaration:
c       call interp_smf(ndecks)
c    input:
c       ndecks   - (integer) the number of cloud decks (1-4)
c    input (in common):
c       nobhgt   - (integer) the number of levels in the temperature-
c                  pressure-height vertical profile
c       height   - (real) the height of each level in the profile
c       pressp   - (real) the pressure profile
c       tempro   - (real) the temperature profile
c       bot      - (real) the height of the (up to 4) cloud layer
c                  bottoms
c       top      - (real) the height of the (up to 4) cloud layer tops
c    output (in common):
c       prb      - (real) the pressures at the cloud layer bottoms
c       prt      - (real) the pressures at the cloud layer tops
c       btem     - (real) the temperatures at the cloud layer bottoms
c       ttem     - (real) the temperatures at the cloud layer tops
c       nlvl     - (integer) number of 100 m height levels
c       heitl    - (real) array of 100 m heights
c       pres     - (real) the pressure every 100 m
c       temp     - (real) the temperature every 100 m
c    common:
c       k6       - (integer) the screen output unit
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              cleaned up the code
c    05/03/89  asl    (505) 678-1570    elton p. avara
c              minor modifications to the code and added the variable
c              pres computations.
c    09/13/89  asl    (505) 678-1570    elton p. avara
c              modified the definitions of pres and temp to give values
c              every 100 m.  also added nlvl, prt, ttem, and heitl.
c              major modifications.
c=======================================================================
c<end>
c***********************************************************************
c
      common /pant/ bot(4),top(4),prb(4),prt(4),btem(4),ttem(4)
      common /qant/ nlvl,heitl(260),pres(260),temp(260)
      common /prof/ nobhgt,height(200),pressp(200),tempro(200)
      common /iodev/k5,k6,k7,k8
c
c-------initialize the height indices
c
      nlvl=0
      l=2
c
c-------find the greatest height level required
c
      xmaxht=0
      do 5 j=1,ndecks
    5 if(xmaxht.lt.top(j)) xmaxht=top(j)
c
c-------the first height level is ground level
c
      nlvl=nlvl+1
      heitl(nlvl)=height(1)
      pres(nlvl)=pressp(1)
      temp(nlvl)=tempro(1)
c
c-------find the next height above ground level which is a multiple of
c       100 m
c
      heitz=100*int(0.01*height(1))
      if(heitz.lt.height(1)) heitz=heitz+100.0
c
c-------increment the height by 100 meters and get the temperature
c       and pressure
c
   10 nlvl=nlvl+1
      heitz=heitz+100.0
      heitl(nlvl)=heitz
c
c-------check the height level against the maximum height required
c
      if(heitl(nlvl)-xmaxht.ge.100.0) then
        nlvl=nlvl-1
        goto 30
      endif
c
c-------loop over the temperature and pressure heights from the
c       radiosonde
c
      do 20 j=l,nobhgt
c
c-------now as soon as the observation height reaches the desired
c       height level, interpolate the temperature to that height.
c
      if (height(j).ge.heitl(nlvl)) then
        slope = (tempro(j)-tempro(j-1))/(height(j)-height(j-1))
        trcpt = -slope * height(j) + tempro(j)
        temp(nlvl) = slope * heitl(nlvl) + trcpt
c
c-------now interpolate to find the pressure
c
c       if(abs(slope).gt.1.0e-4) then
c
c--------------linear trend to temperature in vertical.  the number
c              0.034109527 is the acceleration due to gravity divided
c              by the gas constant
c
c         ratio = ((tempro(j-1) + slope * (heitl(nlvl)-height(j-1))) /
c    1    tempro(j-1))**(-0.034109527/slope)
c       else
c
c--------------use average layer temperature-avtemp-to calculate pressure
c
c
          avtemp = (tempro(j-1) + tempro(j))/2.
          ratio =exp(-0.034109527*(heitl(nlvl)-height(j-1))/avtemp)
c       endif
        pres(nlvl) = pressp(j-1) * ratio
c
c-------go get another height level
c
        l=j
        goto 10
c
c-------have we run out of upperair data?
c
      elseif (j.eq.nobhgt) then
        write(k6,*)'temperature and pressure profiles do not extend 100
     1m above highest cloud top.'
        write(k6,*)'the model needs this information to perform the calc
     1ulations.'
c
        stop
      endif
   20 continue
c
c-------loop over the geometric layers and get the layer bottom values
c       of temperature and pressure.
c
   30 l=2
      do 40 i=1,ndecks
c
c-------loop over the temperature and pressure heights from the
c       radiosonde
c
      do 50 j=l,nobhgt
c
c-------now as soon as the observation height reaches the desired
c       height level, interpolate the temperature to that height.
c
      if (height(j).ge.bot(i)) then
        slope = (tempro(j)-tempro(j-1))/(height(j)-height(j-1))
        trcpt = -slope * height(j) + tempro(j)
        btem(i) = slope * bot(i) + trcpt
c
c-------now interpolate to find the pressure
c
c       if(abs(slope).gt.1.0e-4) then
c
c--------------linear trend to temperature in vertical.  the number
c              0.034109527 is the acceleration due to gravity divided
c              by the gas constant
c
c         ratio = ((tempro(j-1) + slope * (bot(i)-height(j-1))) /
c    1    tempro(j-1))**(-0.034109527/slope)
c       else
c
c--------------use average layer temperature-avtemp-to calculate pressure
c
          avtemp = (tempro(j-1) + tempro(j))/2.
          ratio =exp(-0.034109527*(bot(i)-height(j-1))/avtemp)
c       endif
        prb(i) = pressp(j-1) * ratio
c
c-------go get another geometric layer
c
        l=j
        goto 40
c
c-------have we run out of upperair data?
c
      elseif (j.eq.nobhgt) then
        write(k6,*)'temperature and pressure profiles do not extend 100
     1m above highest cloud top.'
        write(k6,*)'the model needs this information to perform the calc
     1ulations.'
c
        stop
      endif
   50 continue
c
   40 continue
c
c-------loop over the geometric layers and get the layer top values
c       of temperature and pressure.
c
      l=2
      do 60 i=1,ndecks
c
c-------loop over the temperature and pressure heights from the
c       radiosonde
c
      do 70 j=l,nobhgt
c
c-------now as soon as the observation height reaches the desired
c       height level, interpolate the temperature to that height.
c
      if (height(j).ge.top(i)) then
        slope = (tempro(j)-tempro(j-1))/(height(j)-height(j-1))
        trcpt = -slope * height(j) + tempro(j)
        ttem(i) = slope * top(i) + trcpt
c
c-------now interpolate to find the pressure
c
c       if(abs(slope).gt.1.0e-4) then
c
c--------------linear trend to temperature in vertical.  the number
c              0.034109527 is the acceleration due to gravity divided
c              by the gas constant
c
c         ratio = ((tempro(j-1) + slope * (top(i)-height(j-1))) /
c    1    tempro(j-1))**(-0.034109527/slope)
c       else
c
c--------------use average temperature in layer-avtemp-to calculate pressure
c
          avtemp = (tempro(j-1) + tempro(j))/2.
          ratio =exp(-0.034109527*(top(i)-height(j-1))/avtemp)
c       endif
        prt(i) = pressp(j-1) * ratio
c
c-------go get another geometric layer
c
        l=j
        goto 60
c
c-------have we run out of upperair data?
c
      elseif (j.eq.nobhgt) then
        write(k6,*)'temperature and pressure profiles do not extend 100
     1m above highest cloud top.'
        write(k6,*)'the model needs this information to perform the calc
     1ulations.'
c
        stop
      endif
   70 continue
c
   60 continue
c
      return
      end


      subroutine al(lvlct,layr,ti,n1,n2,calw)
c***********************************************************************
c                          subroutine al
c***********************************************************************
c<begin>
c<identification>          name:  al
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  lw
c=======================================================================
c<description>
c    this subroutine computes moist adiabatic liquid water content of
c    non-cirriform cloud types
c=======================================================================
c<called routines>
c    tz     - (subroutine) calculates moist adiabatic lapse rate of
c             temperature and saturation vapor pressure
c=======================================================================
c<parameters>
c    formal declaration:
c       call al(lvlct,layr,ti,n1,n2,calw)
c    input:
c       lvlct    - (integer) cloud type of the cloud deck
c       layr     - (integer) cloud layer
c       ti       - (real) temperature at bottom of cloud layer (c)
c       n1       - (integer) the lowest 100 m height level within the
c                  cloud layer
c       n2       - (integer) the highest 100 m height level within the
c                  cloud layer
c    output:
c       calw     - (real) array of adiabatic liquid water content at
c                  each 100 m within the cloud layer
c    common:
c       bot      - (real) the height of the (up to 4) cloud layer
c                  bottoms
c       top      - (real) the height of the (up to 4) cloud layer tops
c       prb      - (real) the pressures at the cloud layer bottoms
c       prt      - (real) the pressures at the cloud layer tops
c       btem     - (real) the temperatures at the cloud layer bottoms
c       ttem     - (real) the temperatures at the cloud layer tops
c       nlvl     - (integer) number of 100 m height levels
c       heitl    - (real) array of 100 m heights
c       pres     - (real) the pressure every 100 m
c       temp     - (real) the temperature every 100 m
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              cleaned up the code.
c    05/03/89  asl    (505) 678-1570    elton p. avara
c              minor modifications to the code and added the variable
c              pres for calculating dpz.
c    09/12/89  asl    (505) 678-1570    elton p. avara
c              modified the definition of pres to give values every
c              100 m.  also added nlvl, top, btem, ttem, prt, temp,
c              and heitl.  major modifications.
c=======================================================================
c<end>
c***********************************************************************
c
c  author - c. william rogers
c
c  procedure:
c       compute moist adiabatic temperature lapse rate in s/r tz.  then
c  compute water condensed out in air parcel which rises 10m and cools
c  at lapse rate from s/r tz.  lwc is difference between saturated
c  conditions at bottom and top of 10m thick layer.  this procedure
c  is repeated in 10m steps or fraction thereof until lwc at each 100 m
c  in cloud is obtained.  all sublayer variables are for 10m increments
c  except at top of layer or bottom of next higher layer.  for
c  cumuliform cloud reduce lwc to account for entrainment.
c
c  reference:
c       rogers, c.w., j.t. hanley and e.j. mack, 1985: "updating the
c       smith-feddes model", final report contract no. n00228-84-c-3157
c       calspan report no. 7330-1.  calspan corp., p.o. box 400 buffalo,
c       new york 14225. (see section 3.5)
c
c  note:  all computations are carried out for 10m thick layers except
c         for shallower layers ending at bottom of next higher layer.
c
c  local glossary:
c       t      = bottom temperature of 10m layer (c) as input to s/r tz
c                and in computing tn for es2 computation.  it will be
c                (k) when returned by s/r tz and in computation of rho.
c       alw    = adiabatic liquid water content in cloud layer
c                bottom to bottom of next layer and at 10m intervals in
c                between
c       dp     = pressure change of 10m or less interval
c       dtz    = moist adiabatic lapse rate of temperature
c       dz     = 10m interval or less, the latter at bot(layr) and
c                top(layr)
c       es     = saturation vapor pressure at bottom of 10m layer
c       es2    = saturation vapor pressure at top of 10m or less layer
c       ht     = height of height level in cloud layer relative to cloud
c                deck base
c       igo    = control variable
c                if igo = 1 - have reached top of layer
c                       = 0 - have not reached top of layer
c       p      = pressure at bottom of 10m layer
c       rho    = absolute humidity at bottom of 10m layer
c       tn     = top temperature of 10m layer (c) in es2 computation
c                and (k) in rho2 computations
c       rho2   = absolute humidity at top of 10m or less layer
c       rs     = universal gas constant
c       wr     = gas constant for water vapor
c       y      = fraction of adiabatic lwc produced by entrainment
c       z      = height parameter within cloud layer used in s/r al
c
c-----------------------------------------------------------------------
c
      common /pant/ bot(4),top(4),prb(4),prt(4),btem(4),ttem(4)
      common /qant/ nlvl,heitl(260),pres(260),temp(260)
c
      dimension calw(260)
      logical nskip
c
      data a1,b1,c1,a2,b2,c2/8.4897, -13.2191, 4.7295, 10.357,
     1-28.2416, 8.8846 /
      data rs/8.313e07/
      data wr/2.16528e-7/
c
c-------initialize some variables
c
      alw = 0.0
      igo=0
      inc=0
      t=ti
c
      calw(1)=alw
      nskip=.false.
      if(n1.eq.0.or.n2.eq.0) nskip=.true.
c
c-------begin computing the lwc each 100 m within the cloud.
c
      deltaz = 100.
      nz=n1-1
      do 100 i=1,2600
c
c-------compute the lwc in 10 m (or less) height layers
c
      dz=deltaz
c
c-------compute the raob pressure and temperature gradients and height
c       increment for each 100 m height level.
c
      if(inc.eq.0) then
        nz=nz+1
        lvl=nz-n1+2
        if(nskip) then
          p=prb(layr)
          z=bot(layr)
cc        t=btem(layr)-273.16
          delz=top(layr)-bot(layr)
          dz=delz-deltaz*aint(delz/deltaz)
          dpz=(prt(layr)-prb(layr))/delz
          dtz=(ttem(layr)-btem(layr))/delz
          inc=1
        elseif(nz.eq.n1) then
          p=prb(layr)
          z=bot(layr)
cc        t=btem(layr)-273.16
          delz=heitl(n1)-bot(layr)
          dz=delz-deltaz*aint(delz/deltaz)
          dpz=(pres(n1)-prb(layr))/delz
          dtz=(temp(n1)-btem(layr))/delz
          inc=1
        elseif(nz.gt.n1.and.nz.le.n2) then
          p=pres(nz-1)
          z=heitl(nz-1)
cc        t=temp(nz-1)-273.16
          dpz=(pres(nz)-pres(nz-1))/(heitl(nz)-heitl(nz-1))
          dtz=(temp(nz)-temp(nz-1))/(heitl(nz)-heitl(nz-1))
          inc=1
        elseif(nz.gt.n2) then
          p=pres(n2-1)
          z=heitl(n2-1)
cc        t=temp(n2-1)-273.16

          dpz=(prt(layr)-pres(n2))/(top(layr)-heitl(n2))
          dtz=(ttem(layr)-temp(n2))/(top(layr)-heitl(n2))
          inc=1
        endif
      endif
c
c-------if at the top of the cloud layer, change the index to show it.
c
      if(top(layr)-z.lt.deltaz) then
        dz=top(layr)-z
        igo=1
      endif
c
c-------compute the moist adiabatic temperature lapse rate and
c       saturation vapor pressure
c
!     code for tz has been pulled in. note that dtz was changed to dmtz.
!     the value of t is also modified in that subroutine.
!     call tz(p,t,dmtz,es)
!     subroutine tz(p,t,dtz,es)

      g=980.
      cp=.239
      r=.06855
c
c  compute saturation vapor pressure at bottom of 10m sublayer
c
      es=6.112*exp(17.67*t/(t+243.5))
      xl=595.-0.5*t
c
c  compute differential of es with celsius temperature
c
      des=es*((17.67/(t+243.5))-(17.67*t/(t+243.5)**2))
      told = t
      t=t+273.16
c
c  compute moist adiabatic lapse rate of temperature
c
      dmtz=-g*((1.0+.621*es*xl/(p*r*t))/(cp+.621*xl*des/p))
      dmtz=dmtz*2.39e-6


c
c-------compute the mixing ratio and absolute humidity at the bottom of
c       the layer
c
      amr1=0.622*es/(p-es)
      rho=1000.*p/(2.8704*t)
c
c-------get the height, pressure, and moist adiabatic temperature (c) of
c       the top of the 10 m (or less) layer
c
      z=z+dz
      p=p+dpz*dz
      tn=(t-273.16)+dmtz*dz
cc    tn=(t-273.16)+dtz*dz
c
c-------compute the and saturation vapor pressure, mixing ratio, and
c       absolute humidity at the top of the layer
c
c-------try an approximate form for the saturation vapor pressure
c       at the higher level
      dtemp = tn - told
      es2 = es + des*dtemp
c     es2=6.112*exp(17.67*tn/(tn+243.5))
      tn=tn+273.16
      amr2=0.622*es2/(p-es2)
      rho2=1000.*p/(2.8704*tn)
c
c-------compute the increase in the adiabatic liquid water over the 10 m
c       ascent and add to the total liquid water content thus far.
c
      alw=alw+(amr1-amr2)*((rho+rho2)*0.5)
      t=tn-273.16
      calw(lvl)=alw
c
c-------if at the top of a 100 m height level, set the index to show it.
c
      if(igo.ne.1) then
        if(nskip) then
          goto 100
        else
          if(abs(z-heitl(nz)).gt.1.0) goto 100
        endif
      endif
      inc=0
c                                           !modifications nov/87
c-------time to reduce the calw due to entrainment.
c       there are two methods for reducing the calw.  the first is via
c       skatskii's curve which is used on the cumuliform clouds
c       (cb, sc, cu, and ac) and on the unknown cloud type.
c       the second is via warner's curve which is used on the stratiform
c       clouds (st, as, and ns).
c
c
c-------calculate the thickness of the cloud so far in kilometers.
c
      ht=(z-bot(layr))*.001
c
      if((lvlct.eq.10).or.(lvlct.eq.2).or.(lvlct.eq.3).or.(lvlct.eq.5)
     1.or.(lvlct.eq.25)) then
c
c-------reduce calw due to entrainment via skatskii's curve
c                                        !modifications oct/86
        if(ht.lt.0.3) then
          y = -1.667 * (ht - 0.6)
        elseif (ht.lt.1.0) then
          arg1 = b1 * b1 - 4.0 * a1 * (c1 - ht)
          y = (-b1 - sqrt(arg1)) / (2.0 * a1)
        elseif (ht.lt.2.9) then
          arg2 = b2 * b2 - 4.0 * a2 * (c2 - ht)
          y = (-b2 - sqrt(arg2)) / (2.0 * a2)
        else
          y = 0.26
        endif

c                                        !end of modifications oct/86
      elseif ((lvlct.eq.1).or.(lvlct.eq.4).or.(lvlct.eq.6)) then
c
c-------reduce calw due to entrainment via warner's curve.
c
        if(ht.le.0.032) then
          y = -11.0 * ht + 1.0
        elseif (ht.le.0.177) then
          y =  -1.4 * ht + 0.6915
        elseif (ht.le.0.726) then
          y =  -0.356 * ht + 0.505
        elseif (ht.le.1.5) then
          y =  -0.0608 * ht + 0.2912
        else
          y =   0.20
        endif
c
      else
c
c-------do nothing because cirriform clouds or clear and then we
c       shouldn't even be in this subroutine.
c
        y = 1.0
c
      endif
c                                        !end of modifications nov/87
c
c-------now calculate the reduced calw.
c
      calw(lvl) = calw(lvl) * y
c
c-------if at the top of the cloud layer, exit the loop.
c
      if(igo .eq. 1) goto 150
c
  100 continue
c
  150 return
      end
      subroutine tz(p,t,dtz,es)
c***********************************************************************
c                          subroutine tz
c***********************************************************************
c<begin>
c<identification>          name:  tz
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  al
c=======================================================================
c<description>
c    compute moist adiabatic lapse rate of temperature and saturation
c    vapor pressure at bottom of 10m layer.
c=======================================================================
c<called routines>
c    none
c=======================================================================
c<parameters>
c    formal declaration:
c       call tz(p,t,dtz,es)
c    input:
c       p        - (real) pressure at bottom of 10m layer
c       t        - (real) bottom temperature of 10m layer (c)
c    output:
c       dtz      - (real) moist adiabatic lapse rate based on variables
c                  at bottom of 10m layer
c       es       - (real) vapor pressure at bottom of 10m layer
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              cleaned up the code
c=======================================================================
c<end>
c***********************************************************************
c
c  author - c. william rogers
c
c  local glossary:
c       cp      = specific heat of air
c       g       = acceleration of gravity
c       xl      = latent heat of vaporization
c       r       = gas constant
c       des     = differential of es with temperature (c)
c
c  reference:
c       rogers, c.w., j.t. hanley and e.j. mack, 1985: "updating the
c       smith feddes model", final report contract no. n00228-84-c-3157
c       calspan report no. 7330-1.  calspan corp., p.o. box 400 buffalo,
c       new york 14225. (see section 3.5)
c
c-----------------------------------------------------------------------
c
      g=980.
      cp=.239
      r=.06855
c
c  compute saturation vapor pressure at bottom of 10m sublayer
c
      es=6.112*exp(17.67*t/(t+243.5))
      xl=595.-0.5*t
c
c  compute differential of es with celsius temperature
c
      des=es*((17.67/(t+243.5))-(17.67*t/(t+243.5)**2))
      t=t+273.16
c
c  compute moist adiabatic lapse rate of temperature
c
      dtz=-g*((1.0+.621*es*xl/(p*r*t))/(cp+.621*xl*des/p))
      dtz=dtz*2.39e-6
c
      return
      end
      subroutine mvd(ndecks,ictp,xmvd)
c***********************************************************************
c                          subroutine mvd
c***********************************************************************
c<begin>
c<identification>          name:  mvd
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  smfd
c=======================================================================
c<description>
c    determines the mean volume diameter of the droplets within the
c    cloud decks
c=======================================================================
c<called routines>
c    none
c=======================================================================
c<parameters>
c    formal declaration:
c       call mvd(ndecks,ictp,xmvd)
c    input:
c       ndecks   - (integer) the number of cloud decks
c       ictp     - (integer) code for the cloud type (up to 4 cloud
c                  decks)
c    output:
c       xmvd     - (real) mean volume diameter for each cloud
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code as part of main smfd routine
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              cleaned up the code
c    09/14/89  asl    (505) 678-1570    elton p. avara
c              created a separate subroutine for xmvd determination
c=======================================================================
c<end>
c***********************************************************************
c
      dimension ictp(4),xmvd(4)
c
c-------the mvds will be defined to a default representative value for
c       each of the different cloud types
c
      do 10 i = 1,ndecks
c                                                 ! no cloud
      if(ictp(i).eq.0) then
        xmvd(i) = 0.0
c                                                 ! st
      elseif (ictp(i).eq. 1) then
        xmvd(i) = 12.0
c                                                 ! sc
      elseif (ictp(i).eq. 2) then
        xmvd(i) = 10.0
c                                                 ! cu
      elseif (ictp(i).eq. 3) then
        xmvd(i) = 18.0
c                                                 ! ns
      elseif (ictp(i).eq. 4) then
        xmvd(i) = 12.0
c                                                 ! ac
      elseif (ictp(i).eq. 5) then
        xmvd(i) = 18.0
c                                                 ! as
      elseif (ictp(i).eq. 6) then
        xmvd(i) = 10.0
c                                                 ! cb
      elseif (ictp(i).eq.10) then
        xmvd(i) = 25.0
c                                                 ! cirriform
      else
        xmvd(i) = 10.0
c
      endif
   10 continue
c
      return
      end
      subroutine iindex(xmvd,watr,ndecks,aindex,hindex)
c***********************************************************************
c                          subroutine iindex
c***********************************************************************
c<begin>
c<identification>          name:  iindex
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  smfd
c=======================================================================
c<description>
c    this subroutine will determine the icing severity index for both
c    fixed-wing aircraft (aindex) and rotary-wing helicopters (hindex),
c    given the xmvd, watr, and height level temperatures temp, btem, and
c    ttem.
c=======================================================================
c<called routines>
c    none
c=======================================================================
c<parameters>
c    formal declaration:
c       call iindex(xmvd,watr,ndecks,aindex,hindex)
c    input:
c       xmvd     - (real) mean volume diameter for each cloud
c       watr     - (real) output array containing lwc every 100 m
c       ndecks   - (integer) the number of cloud decks
c    output:
c       aindex   - (character) icing severity index for fixed wing
c                  aircraft every 100 m within the cloud
c       hindex   - (character) icing severity index for helicopters
c                  every 100 m within the cloud
c    common:
c       bot      - (real) the height of the (up to 4) cloud layer
c                  bottoms
c       top      - (real) the height of the (up to 4) cloud layer tops
c       prb      - (real) the pressures at the cloud layer bottoms
c       prt      - (real) the pressures at the cloud layer tops
c       btem     - (real) the temperatures at the cloud layer bottoms
c       ttem     - (real) the temperatures at the cloud layer tops
c       nlvl     - (integer) number of 100 m height levels
c       heitl    - (real) array of 100 m heights
c       pres     - (real) the pressure every 100 m
c       temp     - (real) the temperature every 100 m
c       k6       - (integer) the screen output unit
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              cleaned up the code
c    05/03/89  asl    (505) 678-1570    elton p. avara
c              corrected a bug in determining severe lwc and minor
c              modifications to other code.
c    09/13/89  asl    (505) 678-1570    elton p. avara
c              modified the definition of temp, watr, aindex and hindex
c              to give values every 100 m.  also added nlvl.
c=======================================================================
c<end>
c***********************************************************************
c
      integer ndiams
      parameter   (ndiams = 6)

      real freeze,toocld
      parameter   (freeze = 273.15)
      parameter   (toocld = 243.15)
c
      common /pant/ bot(4),top(4),prb(4),prt(4),btem(4),ttem(4)
      common /qant/ nlvl,heitl(260),pres(260),temp(260)
      common /iodev/k5,k6,k7,k8
c
      character*11 aindex(4,260),hindex(4,260)
      dimension xmvd(4),watr(4,2,260),diam(ndiams),xlwcd(ndiams)
c
      data diam/10.0,15.0,20.0,30.0,40.0,50.0/
      data xlwcd/2.50,1.30,0.85,0.65,0.55,0.50/
c
c-------loop over the number of geometric layers which have
c       clouds
c
      do 300 i = 1,ndecks
c
c-------first thing we need to do is to determine the threshold of the
c       severe liquid water content (xlwcs)
c
      if(diam(1).le.xmvd(i).and.xmvd(i).le.diam(ndiams)) then
c
        do 10 k=2,ndiams
c
        if(xmvd(i).ge.diam(k-1).and.xmvd(i).le.diam(k)) then
          slope = (xlwcd(k) - xlwcd(k-1))/(diam(k) - diam(k-1))
          trcpt = -slope * diam(k) + xlwcd(k)
          xlwcs = slope * xmvd(i) + trcpt
        endif
c
   10   continue
c
      else
c
        write(k6,1111)'the calculated mvd = ',xmvd(i),' microns.'
! hongli jiang: change f3.1 to f4.1. 11/27/2013
 1111   format(1x,a,f4.1,a)
        write(k6,*)'this is not within the  possible values of 10-50 mic
     1rons.  thus, the'
        write(k6,*)'icing severity index is not available for this case.
     1'
c
        return
c
      endif
c
c-------now we can determine the thresholds for the moderate and light
c       liquid water contents, xlwcm and xlwcl.
c
      xlwcm = 0.5 * xlwcs
      xlwcl = 0.1 * xlwcs
c
c-------compute the icing severity index for each height level within
c       the cloud deck
c
      n1=watr(i,1,1)+0.001
      n2=watr(i,1,2)+0.001
c
      if(n1.lt.0.or.n2.lt.0) then
        aindex(i,1) = '    missing'
        hindex(i,1) = '    missing'
        goto 300
      elseif(n1.eq.0.or.n2.eq.0) then
        maxcnt=2
      else
        maxcnt=n2-n1+3
      endif
c
      do 200 l=1,maxcnt
c
c-------get the temperature for the height level
c
      if(l.eq.1) then
        t=btem(i)
      elseif(l.eq.maxcnt) then
        t=ttem(i)
      else
        t=temp(n1+l-2)
      endif
c
c-------now to determine the actual icing severity index for both the
c       fixed wing aircraft (aindex) and the rotary wing helicopters
c       (hindex), if the temperature is below freezing.
c
      if((toocld.lt.t).and.(t.le.freeze)) then
c
c-------first, the fixed wing (aindex)
c
        if(watr(i,1,l+2).lt.xlwcl) then
          aindex(i,l) = '  a1: trace'
        elseif (watr(i,1,l+2).lt.xlwcm) then
          aindex(i,l) = '  a2: light'
        elseif (watr(i,1,l+2).lt.xlwcs) then
c
          if(t.le.freeze - 5.0) then
            aindex(i,l) = 'a3:moderate'
          else
            aindex(i,l) = '  a2: light'
          endif
c
        else
c
          if(t.le.freeze - 5.0) then
            aindex(i,l) = ' a4: severe'
          elseif (t.le.freeze - 3.0) then
            aindex(i,l) = 'a3:moderate'
          else
            aindex(i,l) = '  a2: light'
          endif
c
        endif
c
c-------now for the rotary wing (hindex)
c
        if(watr(i,1,l+2).lt.xlwcl) then
          hindex(i,l) = '  h1: trace'
        elseif (watr(i,1,l+2).lt.xlwcm) then
c
          if(t.le.freeze - 5.0) then
            hindex(i,l) = '  h2: light'
          else
            hindex(i,l) = '  h1: trace'
          endif
c
        elseif (watr(i,1,l+2).lt.xlwcs) then
c
          if(t.le.freeze - 10.0) then
            hindex(i,l) = 'h3:moderate'
          elseif (t.le.freeze - 5.0) then
            hindex(i,l) = '  h2: light'
          else
            hindex(i,l) = '  h1: trace'
          endif
c
        else
c
          if(t.le.freeze - 10.0) then
            hindex(i,l) = ' h4: severe'
          elseif (t.le.freeze - 5.0) then
            hindex(i,l) = 'h3:moderate'
          else
            hindex(i,l) = '  h2: light'
          endif
c
        endif
c
      else
c
c-------the temperature is too warm or too cold to produce freezing
c
        aindex(i,l) = 'a0:no icing'
        hindex(i,l) = 'h0:no icing'
c
      endif
c
  200 continue
c
  300 continue
c
      return
      end
      subroutine cltype(ictp,ndecks,cldtyp)
c***********************************************************************
c                          subroutine cltype
c***********************************************************************
c<begin>
c<identification>          name:  cltype
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  smfd
c=======================================================================
c<description>
c    this subroutine will determine the two character abbreviation for
c    the cloud type from the 3dneph integer code for the cloud type.
c=======================================================================
c<called routines>
c    none
c=======================================================================
c<parameters>
c    formal declaration:
c       call cltype(ictp,ndecks,cldtyp)
c    input:
c       ictp     - (integer) code for the cloud type (up to 4 cloud
c                  decks)
c       ndecks   - (integer) the number of cloud decks
c    output:
c       cldtyp   - (character) cloud type in two ascii characters
c    common:
c       k6       - (integer) the screen output unit
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              cleaned up the code
c=======================================================================
c<end>
c***********************************************************************
c
      dimension ictp(4)
      character*2 cldtyp(4)
c
c-------loop over the number of cloud decks
c
      do 100 i = 1,ndecks
c
      if(ictp(i).eq.0) then
        cldtyp(i) = '  '
      elseif (ictp(i).eq. 1) then
        cldtyp(i) = 'st'
      elseif (ictp(i).eq. 2) then
        cldtyp(i) = 'sc'
      elseif (ictp(i).eq. 3) then
        cldtyp(i) = 'cu'
      elseif (ictp(i).eq. 4) then
        cldtyp(i) = 'ns'
      elseif (ictp(i).eq. 5) then
        cldtyp(i) = 'ac'
      elseif (ictp(i).eq. 6) then
        cldtyp(i) = 'as'
      elseif (ictp(i).eq. 7) then
        cldtyp(i) = 'cs'
      elseif (ictp(i).eq. 8) then
        cldtyp(i) = 'ci'
      elseif (ictp(i).eq. 9) then
        cldtyp(i) = 'cc'
      elseif (ictp(i).eq.10) then
        cldtyp(i) = 'cb'
      elseif (ictp(i).eq.25) then
        cldtyp(i) = '??'
      else
        write(6,*)'something is wrong with the cloud type array,'
        write(6,1001) i,ictp(i)
 1001   format(1x,'ictp array.  ictp(',i1,') = ',i2)
c
        stop
      endif
c
  100 continue
c
      return
      end
      subroutine output(cldtyp,watr,xmvd,aindex,hindex,icp,ndecks,title)
c***********************************************************************
c                          subroutine output
c***********************************************************************
c<begin>
c<identification>          name:  output
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  smfd
c=======================================================================
c<description>
c    this subroutine takes the output parameters and writes them out to
c    the screen or printer.
c=======================================================================
c<called routines>
c    inkey - (subroutine) gets a character from the keyboard without
c            echo.  the routine is not available in source code form.
c=======================================================================
c<parameters>
c    formal declaration:
c       call output(cldtyp,watr,xmvd,aindex,hindex,icp,ndecks,title)
c    input:
c       cldtyp   - (character) cloud type in two ascii characters
c       watr     - (real) output array containing lwc and percent liquid
c                  water every 100 m in the cloud
c       xmvd     - (real) mean volume diameter for each cloud
c       aindex   - (character) icing severity index for fixed wing
c                  aircraft every 100 m
c       hindex   - (character) icing severity index for helicopters
c                  every 100 m
c       icp      - (integer) percent (1-100%) cloud coverage for up to 4
c                  cloud decks
c       ndecks   - (integer) the number of cloud decks
c       title    - (character) the title to be displayed in the output
c    output:
c       none
c    common:
c       bot      - (real) the height of the (up to 4) cloud layer
c                  bottoms
c       top      - (real) the height of the (up to 4) cloud layer tops
c       prb      - (real) the pressures at the cloud layer bottoms
c       prt      - (real) the pressures at the cloud layer tops
c       btem     - (real) the temperatures at the cloud layer bottoms
c       ttem     - (real) the temperatures at the cloud layer tops
c       nlvl     - (integer) number of 100 m height levels
c       heitl    - (real) array of 100 m heights
c       pres     - (real) the pressure every 100 m
c       temp     - (real) the temperature every 100 m
c       k5       - (integer) the keyboard input unit
c       k6       - (integer) the screen output unit
c       k7       - (integer) the disk file unit
c       k8       - (integer) the printer output unit
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              major modifications to the code.  completely rewritten.
c    05/03/89  asl    (505) 678-1570    elton p. avara
c              minor modifications to the code.
c    09/14/89  asl    (505) 678-1570    elton p. avara
c              modified the definition of temp, watr, aindex and hindex
c              to give values every 100 m.  also added nlvl, heitl,
c              pres, prt, and ttem.  major modifications.
c    12/21/89  steve albers - modified for profs
c=======================================================================
c<end>
c***********************************************************************
c
      real       freeze
      parameter   (freeze = 273.15)
c
      common /pant/ bot(4),top(4),prb(4),prt(4),btem(4),ttem(4)
      common /qant/ nlvl,heitl(260),pres(260),temp(260)
      common /iodev/k5,k6,k7,k8
c
      dimension icp(4),watr(4,2,260),xmvd(4)
c
      character*70 title
c     character*14 filenm
      character*11 aindex(4,260),hindex(4,260)
c     character*10 scrclr
c     character*8 curpos
c     character*4 scrblk
c     character*3 lerase
      character*2 cldtyp(4)
c     character*1 num(10),ans
      character*1 num(10)
c
c     logical nprt,nfil
c
      data num/'0','1','2','3','4','5','6','7','8','9'/

    1 format(a,$)
c
c-------now to write the cloud data starting with the lowest cloud
c       to the highest cloud.
c
c-------first write the cloud deck data that is not height conditional
c       to the screen (and printer/disk file).
c
      write(k6,26) (int(top(j)),j=1,ndecks)
   26 format(' cloud top (m):',16x,4i12)
c
      write(k6,24) (int(bot(j)),j=1,ndecks)
   24 format(' cloud base (m):',15x,4i12)
c
      write(k6,30) (cldtyp(j),j=1,ndecks)
   30 format(' cloud type:',19x,4(10x,a2))
c
      write(k6,32) (xmvd(j),j=1,ndecks)
   32 format(' drop mvd (microns):',11x,4f12.1)
c
      write(k6,34) (icp(j),j=1,ndecks)
   34 format(' prob. of encountering cloud:',2x,4i12)
c
c-------now write the height conditional cloud data to the screen (and
c       printer/disk file).
c
      do 100 i = 1,ndecks
c
      n1=watr(i,1,1)+0.001
      n2=watr(i,1,2)+0.001
c
      if(n1.lt.0.or.n2.lt.0) then
        goto 100
      elseif(n1.eq.0.or.n2.eq.0) then
        maxcnt=2
      else
        maxcnt=n2-n1+3
      endif
c
c-------write the cloud layer heading
c
      write(k6,41) i
   41 format(/35x,'layer - ',i1//' height pressure  temp    lwc    prob
     1cloud    fix-wng icing    rot-wng icing'/'   (m)    (mb)     (c) (
     2gm/cu m) all liquid   severity index   severity index'/)
c
      do 200 l=1,maxcnt
c
c-------get the temperature for the height level and convert to celsius.
c
      if(l.eq.1) then
        z=bot(i)
        p=prb(i)
        t=btem(i)-freeze
      elseif(l.eq.maxcnt) then
        z=top(i)
        p=prt(i)
        t=ttem(i)-freeze
      else
        nz=n1+l-2
        z=heitl(nz)
        p=pres(nz)
        t=temp(nz)-freeze
      endif
c
      alw=watr(i,1,l+2)
c
c-------convert watr(2,i) from decimal to percent format
c
      prob=100.0*watr(i,2,l+2)
c
c-------write the data to the screen (and printer/disk file).
c
      write(k6,42) int(z),p,t,alw,prob,aindex(i,l),hindex(i,l)
 42   format(i7,f9.1,f7.1,f8.3,f10.1,8x,a11,6x,a11)

c
  200 continue
c
  100 continue
c
   99 return
      end
      subroutine out_laps
     1  (nk,cldtyp,watr,xmvd,aindex,hindex,icp,ndecks,title,
     1         rlwc_laps,prob_laps)
c=======================================================================
c<parameters>
c    formal declaration:
c       call output(cldtyp,watr,xmvd,aindex,hindex,icp,ndecks,title)
c    input:
c       cldtyp   - (character) cloud type in two ascii characters
c       watr     - (real) output array containing lwc and percent liquid
c                  water every 100 m in the cloud
c       xmvd     - (real) mean volume diameter for each cloud
c       aindex   - (character) icing severity index for fixed wing
c                  aircraft every 100 m
c       hindex   - (character) icing severity index for helicopters
c                  every 100 m
c       icp      - (integer) percent (1-100%) cloud coverage for up to 4
c                  cloud decks
c       ndecks   - (integer) the number of cloud decks
c       title    - (character) the title to be displayed in the output
c    output:
c       none
c    common:
c       bot      - (real) the height of the (up to 4) cloud layer
c                  bottoms
c       top      - (real) the height of the (up to 4) cloud layer tops
c       prb      - (real) the pressures at the cloud layer bottoms
c       prt      - (real) the pressures at the cloud layer tops
c       btem     - (real) the temperatures at the cloud layer bottoms
c       ttem     - (real) the temperatures at the cloud layer tops
c       nlvl     - (integer) number of 100 m height levels
c       heitl    - (real) array of 100 m heights
c       pres     - (real) the pressure every 100 m
c       temp     - (real) the temperature every 100 m
c       k5       - (integer) the keyboard input unit
c       k6       - (integer) the screen output unit
c       k7       - (integer) the disk file unit
c       k8       - (integer) the printer output unit
c=======================================================================
c    12/21/89  steve albers - modified for profs map lwc output to
c                             laps grid
c=======================================================================
c<end>
c***********************************************************************
c
      real       freeze
      parameter   (freeze = 273.15)
c
      common /pant/ bot(4),top(4),prb(4),prt(4),btem(4),ttem(4)
      common /qant/ nlvl,heitl(260),pres(260),temp(260)
      common /iodev/k5,k6,k7,k8
c
      dimension icp(4),watr(4,2,260),xmvd(4)
c
      real rlwc_laps(nk),prob_laps(nk)

      character*70 title
c     character*14 filenm
      character*11 aindex(4,260),hindex(4,260)
      character*2 cldtyp(4)
c     character*1 num(10),ans
      character*1 num(10)
c
c     logical nprt,nfil
c
      data num/'0','1','2','3','4','5','6','7','8','9'/

    1 format(a,$)

!     initialize laps arrays
!     do k = 1,nk
!         rlwc_laps(k) = 0.
!         prob_laps(k) = 0.
!     enddo

c
c-------now write the height conditional cloud data to the screen (and
c       printer/disk file).
c
      do 100 i = 1,ndecks
c
      n1=watr(i,1,1)+0.001
      n2=watr(i,1,2)+0.001
c
      if(n1.lt.0.or.n2.lt.0) then
        goto 100
      elseif(n1.eq.0.or.n2.eq.0) then
        maxcnt=2
      else
        maxcnt=n2-n1+3
      endif

      z_new = -999.
      do 200 l=1,maxcnt
c
c-------get the temperature for the height level and convert to celsius.
c
      if(l.eq.1) then
        z=bot(i)
        p=prb(i)
        t=btem(i)-freeze
      elseif(l.eq.maxcnt) then
        z=top(i)
        p=prt(i)
        t=ttem(i)-freeze
      else
        nz=n1+l-2
        z=heitl(nz)
        p=pres(nz)
        t=temp(nz)-freeze
      endif
c
      alw=watr(i,1,l+2)
c
c-------convert watr(2,i) from decimal to percent format
c
      prob=100.0*watr(i,2,l+2)
c
c-------write the data to the screen (and printer/disk file).
c

      alw_old = alw
      prob_old = prob
      z_old = z_new

      z_new = zcoord_of_pressure(p*100.)

      if(int(z_old) .ne. int(z_new) .and. z_old .ne. -999.)then
!         interpolate to laps grid point
          iz = int(z_new)
          frac = (z_new - float(iz)) / (z_new - z_old)
          rlwc_laps(iz) = alw_old * frac + alw * (1.0 - frac)
          prob_laps(iz) = prob_old * frac + prob * (1.0 - frac)
      endif


!      write(k6,42) int(z),p,t,alw,prob,aindex(i,l),hindex(i,l)
! 42   format(i7,f9.1,f7.1,f8.3,f10.1,8x,a11,6x,a11)

c
  200 continue
c
  100 continue
c
   99 return
      end

      subroutine lwc(ndecks,icp,ictp,watr,xmvd)
c***********************************************************************
c                          subroutine lwc
c***********************************************************************
c<begin>
c<identification>          name:  lwc
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  smfd
c=======================================================================
c<description>
c    this routine controls the computation of liquid water content (lwc)
c    and the probability the cloud is all liquid.  it calls the
c    microphysics portion of the code by rtneph cloud deck.
c=======================================================================
c<called routines>
c    lw     - (subroutine) calculates lwc and the probability the cloud
c             is all liquid for super-cooled clouds
c=======================================================================
c<parameters>
c    formal declaration:
c       call lwc(ndecks,icp,ictp,watr)
c    input:
c       ndecks   - (integer) number of cloud decks present
c       icp      - (integer) percent cloud amount in layers
c       ictp     - (integer) code for the cloud type (up to 4 cloud
c                  decks)
c    output:
c       watr     - (real) output array containing lwc and percent liquid
c                  water for every 100 m in the cloud layers
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              cleaned up the code
c    09/13/89  asl    (505) 678-1570    elton p. avara
c              modified the definition of watr to give values every
c              100 m.  minor modifications.
c=======================================================================
c<end>
c***********************************************************************
c
c  description:
c       loop 10 controls the stepping through the observed rtneph cloud
c  decks.  subroutine lw is called to calculate the liquid water content
c  (lwc) and thermodynamic phase (tdp).
c
c  local glossary:
c       irain     - no rain (1) flag
c       lay       - rtneph cloud deck index
c       ictp(4)   - cloud type by 3dneph cloud deck
c
c-----------------------------------------------------------------------
c
      dimension watr(4,2,260),ictp(4),icp(4),xmvd(4)
c
c ------------------- start of executable code -------------------------
c
c-------set the rain flag to no rain permenantly
c
      irain = 1
c
      do 10 lay=1,ndecks
c
c-------set the first two values of watr to indicate no clouds
c
      watr(lay,1,1)=-2.0
      watr(lay,1,2)=-2.0
c
c-------if the cloud type is 'clear', skip it.
c
      if(ictp(lay).ne.0) call lw(lay,irain,ictp(lay),watr,icp(lay)
     1                                          ,xmvd(lay))
c
  10  continue
c
      return
      end
      subroutine lw(lay,irain,lvlct,watr,neph,rmvd)
c***********************************************************************
c                          subroutine lw
c***********************************************************************
c<begin>
c<identification>          name:  lw
c                          type:  fortran-77 subroutine
c                      filename:  ice.for
c                        parent:  lwc
c=======================================================================
c<description>
c    this subroutine calculates the lwc and probability the cloud water
c    is all liquid for super-cooled clouds, and calls s/r al for each
c    height level within a cloud deck.
c=======================================================================
c<called routines>
c    al     - (subroutine) calculates moist adiabatic lwc of non-
c             cirriform cloud types
c=======================================================================
c<parameters>
c    formal declaration:
c       call lw(lay,irain,lvlct,watr,neph)
c    input:
c       lay      - (integer) cloud layer index (1-4)
c       irain    - (integer) rain (2) or no rain (1)
c       lvlct    - (integer) cloud type of the cloud deck
c       neph     - (integer) percent cloud amount in cloud deck
c    output:
c       watr     - (real) output array containing lwc and percent liquid
c                  water for every 100 m in the cloud layers
c    common:
c       bot      - (real) the height of the (up to 4) cloud layer
c                  bottoms
c       top      - (real) the height of the (up to 4) cloud layer tops
c       prb      - (real) the pressures at the cloud layer bottoms
c       prt      - (real) the pressures at the cloud layer tops
c       btem     - (real) the temperatures at the cloud layer bottoms
c       ttem     - (real) the temperatures at the cloud layer tops
c       nlvl     - (integer) number of 100 m height levels
c       heitl    - (real) array of 100 m heights
c       pres     - (real) the pressure every 100 m
c       temp     - (real) the temperature every 100 m
c       k6       - (integer) the screen output unit
c=======================================================================
c<history>
c    11/15/87  udri   (513) 229-3921    james k. luers
c              delivered basic source code
c    02/13/89  asl    (505) 678-1570    elton p. avara
c              cleaned up the code
c    05/03/89  asl    (505) 678-1570    elton p. avara
c              minor modifications to the code and added the variables
c              pres and top.
c    09/13/89  asl    (505) 678-1570    elton p. avara
c              modified the definitions of pres and watr to give
c              values every 100 m.  also added nlvl, ttem, temp, lay,
c              and heitl.  major modifications.
c=======================================================================
c<end>
c***********************************************************************
c
c  description:
c       authors - m. d. dykton and c. w. rogers.  loop 200 finds iloc
c                 and itmp and puts the calculated lwc in watr.
c
c  local glossary:
c       iloc       = (integer) location of point within cloud deck
c       itmp       = (integer) temperature converted to table index
c       lay        = (integer) cloud layer
c       permax(5,10,2)
c                  = (real) table of the percent of maximum lwc for a
c                    given cloud type (10) at intervals above the cloud
c                    base (5) for raining or non-raining clouds (2).
c                    see tn 74-4, figures 4, 5, & 6 on pp 6-8 and
c                    paragraph 4a, p 4.  iloc, lvlct, and irain are the
c                    indices for permax.
c       t          = (real) temperature in degrees kelvin
c       themax(10,10)
c                  = (real) table of maximum lwc in g/m**3 that can
c                    occur for a given cloud type (10) and for a given 5
c                    degree temperature interval (10).  see calspan
c                    final report, table 5, p 26.
c
c  references:
c       1) usafetac technical note 74-4, "a synoptic-scale model for
c          simulating condensed atmospheric moisture", april 1974, by
c          capt robert g. feddes.
c
c       2) rogers, c.w., hanley j.t. and e.j. mack, 1985: "updating the
c          smith-feddes model", final report contract no. n00228-84-c-
c          3157 calspan report no. 7330-1. calspan corp., p.o. box 400
c          buffalo, new york 14225.
c
c-----------------------------------------------------------------------
c
      common /pant/ bot(4),top(4),prb(4),prt(4),btem(4),ttem(4)
      common /qant/ nlvl,heitl(260),pres(260),temp(260)
      common /iodev/k5,k6,k7,k8
c
      dimension calw(260),watr(4,2,260),themax(10,10),permax(5,10,2)
c
      data themax/.2,.25,.30,.35,.4,.45,.50,.55,.6,.6,
     1            .35,.4,.45,.50,.55,.60,.65,3*.7,10*3.,.35,.4,
     2            .45,.5,.6,.6,.75,3*.90,.25,.30,.35,.40,.40,
     3            .45,.6,3*.7,.20,.25,.25,.3,.35,.40,.45,3*.5,
     4            3*.15,3*.2,4*.25,4*.1,3*.15,3*.2,
     5            4*.05,3*.1,3*.15,10*6.5/
c
      data permax/5*.4,.38,.62,.74,.58,.37,.40,.60,.80,
     1           .95,.74,.38,.62,.74,.58,.37,.38,.62,.74,.58,.37,
     2            20*.4,.37,.57,.76,.90,.82,.93,.77,.62,.47,.32,
     3           .96,.88,.74,.58,.37,1.9,1.6,1.,.5,.4,.96,.88,
     4           .74,.58,.37,.96,.88,.74,.58,.37,.93,.77,.62,
     5           .47,.32,.93,.77,.62,.47,.32,.93,.77,.62,.47,
     6           .32,.93,.77,.62,.47,.32,2.75,1.95,1.,.63,.47/
c
c ******************* start of executable code *************************
c
c  if cloud percentage is zero, skip microphysics.
c
      if(neph.eq.0) goto 100
c
c  set s/r al parameter
c
      ti=btem(lay)-273.16
c
c-------determine which 100 m height levels are within this cloud layer.
c
      n1=1
      n2=nlvl
      do 400 j=1,nlvl
      if(bot(lay).ge.heitl(j)) n1=j+1
      if(top(lay).gt.heitl(j)) n2=j
  400 continue
      nmax=n2-n1+3
c
      if(n1.gt.n2) then
        n1=0
        n2=0
        nmax=2
      endif
c
      dif=top(lay)-bot(lay)
c
c-------store the indices which indicate which 100 m height levels are
c       within the cloud in the first two elements of watr.
c
      watr(lay,1,1)=n1
      watr(lay,1,2)=n2
c
c-------if the cloud type is not cirriform, calculate the adiabatic
c       liquid water content within the cloud layer
c
      if(lvlct.lt.7.or.lvlct.gt.9) call al(lvlct,lay,ti,n1,n2,calw)
c
c-------get the liquid water content at all height levels within the
c       cloud layer
c
      do 200 l=1,nmax
c
c-------get the height and temperature for each height level within the
c       cloud layer
c
      if(l.eq.1) then
        z=bot(lay)
        t=btem(lay)
      elseif(l.lt.nmax) then
        z=heitl(n1+l-2)
        t=temp(n1+l-2)
      else
        z=top(lay)
        t=ttem(lay)
      endif
c
c  scale temperature into an index (itmp) for the array themax
c
      itmp = (t - 238.0) / 5.0
      if(itmp.gt.10) itmp = 10
      if(itmp.lt. 1) itmp = 1
c
c  calculate percent height of point above base within cloud
c  deck and scale it into an index for the table permax.
c
      hgtpct = (z - bot(lay)) / dif
      iloc = hgtpct * 5. + 1.
      if(iloc.lt.1) iloc = 1
      if(iloc.gt.5) iloc = 5
c
c  if not cirriform cloud, modify adiabatic cmc.
c
      if(lvlct .lt. 7 .or. lvlct .gt. 9) then
c
c  no reduction in cmc near cloud top if stratus type clouds
c  (st, as, and ns).  for cumuliform clouds reduce cmc near cloud top
c  if height greater than 80% of cloud depth
c
        if(lvlct.ne.1.and.lvlct.ne.4.and.lvlct.ne.6) then
          if(hgtpct.ge.0.8) then
            ll=lvlct
            if(ll.eq.25) ll=3
            deltop=5.0*(hgtpct-0.8)
            calw(l)=calw(l)*(1.0-deltop*(1.0-permax(iloc,ll,1)))
          endif
        endif
c
c  if precipitation is present, modify the liquid water content
c
        if(irain.eq.2) then
          perdif=permax(iloc,lvlct,2)-permax(iloc,lvlct,1)
          if(perdif .le. 0.0) then
c
c  decrease cmc in precipitation
c
            calw(l)=calw(l)*permax(iloc,lvlct,2)/permax(iloc,lvlct,1)
          else
c
c  increase cmc in precipitation
c
            pdif=perdif/permax(iloc,lvlct,2)
            calw(l)=calw(l)/(1.0-pdif)
          endif
        endif
c
        watr(lay,1,l+2)=calw(l)
c
c  calculate fractional probability cloud water is all liquid
c
        if( t .gt. 273. ) then
c
c  all cloud water is liquid if temperature > 273.
c
          watr(lay,2,l+2) = 1.
        elseif( t .lt. 233. ) then
c
c  all cloud water is ice if temperature < 233.
c
          watr(lay,2,l+2) = 0.
        else ! t is between 233 and 273
c
c  check cloud type if (233. <= temperature <= 273.)
c  there is super-cooled water
c
          if(lvlct.eq.3.or.lvlct.eq.10) then ! cu type cloud, 233 < t < 273
c
c  super-cooled unstable - cumulus or cumulonimbus clouds
c  -.03927 = (pi/2.) / -40.
c
            watr(lay,2,l+2)=cos(-.03927*(t-273.))
          else ! not cu type cloud (stratiform cloud), 233 < t < 273
c
c  super-cooled stable - stratiform clouds (0c to -30c)
c  cloud water is depleted by a temperature dependent hyperbolic
c  tangent function-no cloud water at -40c and all water at 0c
c  the "width" of the drop-off between all water and no water
c  in the cloud is controlled by the denominator of argument.
c  the middle of the dropoff is currently at -20c(253k) since
c  this value ensures mathematical consistentcy
c
c
c           arg=(t-253.)/7.5
c           fac=(tanh(arg)+1.0)/2.
c
c    11/6/90 change-use a linear factor for lwc depletion by glaciation
c    instead of the hyperbolic profile above (ramps from -10c to -30c)
c
            if(t.ge.263.) then
              fac=1.0
            elseif(t.le.243.) then
              fac=0.0
            else
              fac=(t-243.)/20.
            endif

            watr(lay,1,l+2)=watr(lay,1,l+2)*fac
            watr(lay,2,l+2) = exp(.0909091*(t-273.)) - .0263
          endif ! cu type cloud
        endif ! t

        if(t .le. 233.)watr(lay,1,l+2)=0. ! case for cu or stratiform cloud < 40c

      else ! cirroform cloud
c
c  if cirrifrom cloud do not use adiabatic lwc, use table.
c
        watr(lay,1,l+2) = themax(itmp,lvlct) * permax(iloc,lvlct,irain)
c
c  cloud is all ice
c
        watr(lay,2,l+2) = 0.0
      endif
c
c-------check to make sure that the lwc calculated by this program
c       is non-negative.  if it is negative then tell the user to
c       check their input values for errors.
c
      if(watr(lay,1,l+2).lt.0.0) then
        write(k6 ,*)'the model has calculated a negative lwc.'
        write(k6 ,*)'please check your input data for errors.'
        stop
      endif
c
  200 continue
      ndif=n2-n1
      nmid=(ndif)/2+2
      nmid1=nmid+1
      if((ndif/2)*2.eq.ndif) then
       tmid=temp(nmid)
      else
       tmid=(temp(nmid)+temp(nmid1))/2.
      endif
!     call dsd(lay,lvlct,irain,iloc,permax,watr,tmid,n1,n2,rmvd)
c
  100 return
      end

      subroutine dsd
     g              ( lay, lvlct, irain, iloc, permax,
     b                watr,tempk,n1,n2,mvd )
c
c    purpose
c         this subroutine calculates drop size distribution
c    parameter list
c         given
c              lvlct (i) = cloud type (1-10)
c              irain  (i) = no rain (1) or rain (2)
c              iloc   (i) = location of a calculation point
c                           within a cloud deck
c              lvl    (i) = cloud level
c              permax(5,10,2) (r) = percent of maximum lwc that can
c                                   a function of cloud type and
c                                   relative cloud position
c         yielded
c              rndrps(4,40) (r)    = drop size ditribution at
c                                    the mid-point of each cloud
c                                    deck
c                        mvd = mean volume diameter per layer
c
c    local glossary
c              drpcld (r) = radii at which cloud dsd is evaluated
c              cldrps (r) = number of drops/cm**3/micron radius interval
c                           centered at a particular radius
c              ilayr  (r) = layer number (same as lay)
c              lwcc   (r) = cloud condensed moisture
c              lwccl  (r) = cloud liquid content
c              lwcci  (r) = cloud ice content
c              lwcr   (r) = rain condensed moisture
c              lwcrl  (r) = rain liquid content
c              lwcri  (r) = rain ice content
c lwc=total condensed moisture, precipitation or no precipitation, liqui
c  cloud.  used in both no precipitation and precipitation cases.
c  tdp=total liquid condensed moisture, precipitation or no precipitatio
c  tdp=lwc for liquid cloud
c  tdp=0.0 for ice cloud
c  used only if all moisture is only in cloud
c              tempk  (r) = layer midpoint temperature (degrees kelvin)
c
c.......................................................................
c
c
      dimension  voli(40),cldrps(4,40),rndrps(10,4,2),
     1           permax(5,10,2), drpcld(40),watr(4,2,260)
      real       lwc, lwcc, lwcr, lwccl, lwcci, lwcrl, lwcri
      real mvd,mvdi(40)

      data drpcld/1.0,3.0,5.0,7.0,9.0,11.0,13.0,15.0,17.0,19.0,
     1 21.0,23.0,25.0,27.0,29.0,31.0,33.0,35.0,37.0,39.0,41.0,
     2 43.0,45.0,47.0,49.0,51.0,53.0,55.0,57.0,59.0,61.0,63.0,
     3 65.0,67.0,69.0,71.0,73.0,75.0,77.0,79.0/
c
c ******************** start of executable code ************************
c
c     presetting lwc variables
c
      lwc   = -99.0
      lwcc  = -99.0
      lwcr  = -99.0
      lwccl = -99.0
      lwcci = -99.0
      lwcrl = -99.0
      lwcri = -99.0
c
c   for the time being we will calculate a drop size distribution only
      ilayr=lay

c   at the mid-point of a cloud deck
      ndif=n2-n1
      nmid=ndif/2+1
      nmid1=nmid+1
      if((ndif/2)*2.eq.ndif) then
       lwc=watr(lay,1,nmid)
      else
       lwc=(watr(lay,1,nmid)+watr(lay,1,nmid1))/2.
      endif
      tdp=lwc
c  if no precipitation at this level all moisiture is cloud water
      if ( .not. ( irain .eq. 2 ) ) go to 20
c
c         it's precipitating; calculate the lwc percent difference
c         between precipitating & non-precipitating cloud types
 10   dif  = permax(iloc,lvlct,2) - permax(iloc,lvlct,1)
      pdif = dif / permax(iloc,lvlct,2)
      if (.not.( dif .le. 0. )) go to 60
c
c     all liquid is in cloud form - no precipitation
c
   20 lwccl = tdp
      lwcci = lwc - lwccl
c          cloud liquid content
c          cloud ice content
c
c     following computes and stores drop size distributions
c     corresponding to lwccl and lwcci
c
c
      if(lwccl .gt. 0.)then
          call ndrops
     g           ( lvlct, lwccl, drpcld, 40,lay,
     y             cldrps    )
c
c-----calculate the volumes of each of the radii categories
c
      totvol=0.0
      pi=acos(-1.0)
      do 25 k=1,40
         voli(k)=(4.0/3.0)*pi*(drpcld(k)*1e-4)**3
         iwatr=int(cldrps(lay,k))
         totvol=totvol+iwatr*voli(k)
  25  continue
c
c-----calculate the mean volume diameter for this layer
c
      do 35 k=1,40
         l=41-k
!        write(3,*) cldrps(lay,k)
         iwatr=int(cldrps(lay,l))
c         write(3,*)iwatr,voli(l),totvol
         mvdi(l)=iwatr*voli(l)/totvol
c         write(3,*)mvdi(l)
   35 continue
      totmvd=0.0

      mvd=0.0
      do 40 k=1,40
         l=41-k
         totmvd=totmvd+mvdi(l)
         if (totmvd.ge.0.5) then
             mvd=drpcld(l)
c            write(3,*)mvd
             go to 45
         endif
   40 continue
      endif ! if no lwc
c
   45 if(lwcci.eq.0.)go to 50
c
c     call ndrops
c    g           (8, lwcci, drpcld, 8,
c    y             watr(13,ilayr)     )
c
c
   50 continue
      return
c
c     following treats cases where some liquid is in rain form
c
 60   continue
c split lwc into rain and cloud water
      lwcr = pdif * lwc
      lwcc = lwc - lwcr
c
c     tempk is temperature at the layer midpoint
      if((lvlct .ge. 7) .and. (lvlct .le. 9)) go to 400
c          if temperature > 273. then all moisture is liquid
      if ( tempk .gt. 273. ) lwccl = 1.
c          if temperature < 233. then all moisture is ice
      if ( tempk .lt. 233. ) lwccl = 0.
c          skip if no super-cooled liquid is present
      if ( tempk .gt. 273. .or. tempk .lt. 233. ) go to 210
      if ( lvlct .eq. 3 .or. lvlct .eq. 10 ) go to 200
c
c     super-cooled stable -- stratiform clouds
      lwccl = exp(.0909091*(tempk-273.))-.0263
      go to 210
c
 200  continue
c     super-cooled unstable -- cumulus or cumulonimbus clouds
c     -.03927 = ( pi / 2. ) / -40.
      lwccl = cos (-.03927 * (tempk - 273.))
c
 210  continue
c  fractional probability cloud is all water
      water=lwccl
c set lwccl for liquid cloud processing
c  lwccl=1.0. insures lwccl converts to lwcc in 2nd statement hence.
      lwccl=1.0
      if(water .eq. 0.0) lwccl=0.0
c  liquid cloud variables set for s/r ndrops and s/r wrtout.
      lwccl = lwccl * lwcc
      lwcci = lwcc - lwccl
      watl = lwccl
      wati = lwcci
      go to 500
c  ice cloud variable set for s/r ndrops and s/r wrtout
c  if cirriform cloud then no cloud liquid content and probability of cl
c  condensed moisture being all liquid is zero
  400 watl=0.0
      wati=lwcc
      lwcci=lwcc
      lwccl=0.0
c
c     following computes and stores drop size distributions
c     corresponding to lwccl and lwcci
  500 continue
c
c
      if(lwccl .gt. 0.)then
          call ndrops
     g           ( lvlct, lwccl, drpcld, 40,lay,
     :             cldrps  )
c
c-----calculate the volumes of each of the radii categories
c
      totvol=0.0
      pi=acos(-1.0)
      do 525 k=1,40
         voli(k)=(4.0/3.0)*pi*(drpcld(k)*1e-4)**3
         iwatr=int(cldrps(lay,k))
         totvol=totvol+iwatr*voli(k)
c         twatr=watr(j,ilayr)
c         totvol=totvol+twatr*voli(k)
 525  continue
c
c-----calculate the mean volume diameter for this layer
c
      do 535 k=1,40
         l=41-k
c         twatr=watr(j,ilayr)
         iwatr=int(cldrps(lay,l))
c        write(3,*)iwatr,voli(l),totvol
c         mvdi(l)=twatr*voli(l)/totvol
         mvdi(l)=iwatr*voli(l)/totvol
c         write(3,*)mvdi(l)
  535 continue
      totmvd=0.0
      mvd=0.0
      do 540 k=1,40
         l=41-k
         totmvd=totmvd+mvdi(l)
         if (totmvd.ge.0.5) then
             mvd=drpcld(l)
c            write(3,*)mvd
             go to 45
         endif
  540 continue
      endif ! if no lwc
c
c
      if(lwcci.eq.0.)go to 90
c
c      call ndrops
c    g           ( 8, lwcci, drpcld,40,
c    y             watr(13,ilayr)  )
c
c
   90 continue
c     the factor qz is incorporated into the model to account for
c     the fact that "the differential fall velocities for liquid
c     and solid water will decrease the total condensed moisture
c     below the freezing level."  (tn 74-4, p 11, para. a)
c     in those layers whose heights are below the freezing level,
c     the condensed moisture is halved, for convective cloud types.
      qz = 1.
c
c     follow same computational procedure for rain as for clouds (above)
      if ( tempk .gt. 273. ) lwcrl = 1.
      if ( tempk .lt. 233. ) lwcrl = 0.
      if ( tempk .gt. 273. .and. lvlct .eq. 10 ) go to 230
      if ( tempk .gt. 273. .and. lvlct .eq. 3 ) go to 230
      if ( tempk .gt. 273. .or.  tempk .lt. 233. ) go to 240
      if(lvlct .eq. 3 .or. lvlct .eq. 10 ) go to 220
c
c     super-cooled stable -- non-convective so factor doesn't apply
      lwcrl = exp(.0909091 * (tempk-273.))-.0263
c  precipitation at temperatures below zero is frozen. stratiform clouds
      lwcrl=0.0
      go to 240
c
 220  continue
c     super-cooled unstable
      lwcrl =  cos(-.03927  *  (tempk - 273.))
c  preciptation at temperatures below zero is frozen.  cumuliform clouds
      lwcrl=0.0
c
 230  continue
c     branch here if convective clouds; test to see if temperature
c     above freezing (so that layer is below freezing point).
      if ( tempk .gt. 273.) qz = .5
c
 240  continue
c
      lwcrl = (lwcrl * qz) * lwcr
      lwcri = (lwcr * qz) - lwcrl
c
c     lwc = lwccl + lwcci + lwcrl + lwcri
c
c
c     compute drop size distribution (number/m**3/micron) at 300
c     micron intervals from 150 to 2850 microns, corresponding to ...
c
c     lwcrl (rain liquid):
      if ( .not. ( lwcrl .ne. 0 ) ) go to 105
      r = -150.
      do 100 ia = 1,10
           r = r+300.
  100      rndrps(ia,lay,1)= 20.*exp((-.004555*r)/lwcrl**.25)
c          see equation, p 8, tn 74-4 (which is from p 17, tn 74-1)
c     enddo
c
  105 continue
c
c     lwcri (rain ice):
      if (.not.( lwcri .ne. 0 )) go to 115
      r = -150.
      do 110 ia = 1,10
           r = r+300.
  110      rndrps(ia,lay,2)= 20. * exp((-.004555*r)/lwcri**.25)
c          see equation, p 8, tn 74-4
c     enddo
c
c
  115 return
      end
c
c

      subroutine ndrops(lvlct,lwc,dropsz,n,lay,
     1 cldrps)
cc
cc   the cloud drop size distribution calculated herein is based
cc   on berry and reinhardt(journal of atmospheric science, 1974)
cc   the total drop concentration must be specified.  for now,
cc   this is specified for all cloud types.
cc   later i intend to make this a function of cloud base vertical
cc   velocity and temperature.
cc
cc   parameter list
cc
cc    given
cc     lvlct(i) = cloud type
cc     lwc(r)    = liquid water content(gm/cm**3)
cc     dropsz(r) = number of drops per 2 micron diameter interval
cc     n(i)      = dimension of dropsz
cc     lay(i)    = cloud deck index
cc    yielded
cc     cldrps(r) = number of drops/cm**3/2 micron interval
cc
cc   -------------------------------------------------------------
      real pi,density
      parameter(pi=3.14159,density=1.0)
      real lwc,cldrps(4,40),dropsz(40),tnd(10)
      integer gmafct
      data gamma/2./
      data tnd/2.5e2,3.0e2,5.0e2,2.5e2,2.5e2,3.0e2,2.5e2,
     1 3.0e2,2.5e2,5.0e2/
cc
cc
      lwc=lwc*1e-6
      gmafct=1
      do 100 i=1,gamma
       gmafct=gmafct*i
 100  continue
cc
      g=(1+gamma)**(1+gamma)/float(gmafct)
cc
      ddiam=2.0e-4
      do 200 i=1,40
       dia=dropsz(i)/1e4
       dvol=pi*dia*dia*ddiam/2.
       vol=pi*dia**3/6.
       s=tnd(lvlct)*vol*density/lwc
       arg=(1+gamma)*s
       if(arg.gt.75) arg=75
       f=tnd(lvlct)**2*g*s**gamma*exp(-arg)/lwc
       cldrps(lay,i)=f*dvol*density
 200  continue
      return
      end




