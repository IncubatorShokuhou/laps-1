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
c
c
        subroutine read_obs_i(infile,maxstns,atime,num_meso,
     &   num_saos,num_sfc,stations,lat_s,lon_s,elev_s,wx_s,cover_s,
     &   hgt_ceil,hgt_low,t_s,td_s,dd_s,ff_s,ddg_s,ffg_s,pr_s,sr_s,
     &   istatus)
c
c*******************************************************************************
c
c       routine to read sao and mesonet surface obs written by the laps
c       'lapsdata' program.
c
c       changes:
c               p.a. stamus     11-29-88        original version.
c                               02-01-90        version for interactive mdat.
c                               02-14-91        add solar radiation.
c
c       input/output:
c
c          variable      var type    i/o    description
c         ----------    ----------  -----  -------------
c          infile          a*70       i     directory where input data is.
c          maxstns          i         i     max number of stations allowed
c          atime           a*24       o     data time: dd-mmm-yyyy hh:mm
c          num_meso         i         o     number of mesonet stations in file
c          num_saos         i         o     number of sao stations in file
c          num_sfc          i         o     total number of surface obs.
c          stations        a*3 a      o     array of the station names
c          lat_s            ra        o     latitude of the stations
c          lon_s            ra        o     longitude of the stations
c          elev_s           ra        o     elevation of the stations (m)
c          wx_s            a*8 a      o     array of observed weather
c          cover_s          ra        o     cloud cover (tenths)
c          hgt_ceil         ra        o     ceiling height (m)
c          hgt_low          ra        o     height lowest cloud (m)
c          t_s              ra        o     temperature (f)
c          td_s             ra        o     dewpoint (f)
c          dd_s             ra        o     wind direction (deg)
c          ff_s             ra        o     wind speed (kt)
c          ddg_s            ra        o     gust wind direction (deg)
c          ffg_s            ra        o     gust wind speed (kt)
c          pr_s             ra        o     pressure variable - see note #2.
c          sr_s             ra        o     solar radiation.
c          istatus          i         o     status flag: 1 = normal
c                                                       -1 = file not found
c
c       user notes:
c
c       1.  arrays should be dimensioned 'maxstns' in the calling program,
c           with maxstns >= 60 for this routine.
c
c       2.  the pressure variable for the mesonet data is station pressure
c           in millibars.  the pressure variable for saos is altimeter setting
c           in millibars.  the mesonet stations are always the first 'num_meso'
c           stations in the pr_s array, and saos are the rest.  no corrections,
c           changes, reductions, etc., are made to the data in this routine.
c
c       3.  infile is now passed in, complete and ready for the open statement.
c
c*******************************************************************************
c
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cover_s(maxstns), hgt_ceil(maxstns), hgt_low(maxstns)
        real t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns), ffg_s(maxst
     1ns)
c
        character stations(maxstns)*3, wx_s(maxstns)*8
        character atime*24, infile*256
c
        istatus = 0
c
c.....  open the obs data file.  check for a
c.....  'file not found', and notify the user if necessary.
c
        open(1,iostat=ios,file=infile,status='old')
        if(ios .eq. 29) then            !file not found
          print *,
     &    ' +++++ obs file (sao & mesonet data) not available. +++++'
          istatus = -1
          return
        endif
c
c.....  now read the time and number of stations in the file, then read
c.....  the data.
c
        read(1,901) atime,num_meso,num_sfc
901     format(1x,a24,2i6)
c
        do k=1,num_sfc
          read(1,902)stations(k),lat_s(k),lon_s(k),elev_s(k),wx_s(k),
     &           cover_s(k),hgt_ceil(k),hgt_low(k),t_s(k),td_s(k),
     &           dd_s(k),ff_s(k),ddg_s(k),ffg_s(k),pr_s(k),sr_s(k)
        enddo !k
902     format(1x,a3,2f7.2,1x,f5.0,1x,a8,1x,f5.1,2(1x,f7.1),1x,2f6.1,
     &         1x,4(1x,f5.0),1x,f6.1,1x,f6.1)
c
        num_saos = num_sfc - num_meso
        istatus = 1                     ! normal return
c
        return
        end
