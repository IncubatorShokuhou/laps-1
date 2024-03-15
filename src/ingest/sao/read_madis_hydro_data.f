c
c  subroutine to read the file "madis - hydrological surface" 
c
      subroutine read_madis_hydro_netcdf(nf_fid, icchecknum, 
     +     qcchecknum, maxstaticids, ninventorybins, recnum, 
     +     elevation, latitude, longitude, precip12hr, precip12hrqcd, 
     +     precip1hr, precip1hrqcd, precip24hr, precip24hrqcd, 
     +     precip3hr, precip3hrqcd, precip5min, precip5minqcd, 
     +     precip6hr, precip6hrqcd, precipaccum, precipaccumqcd, 
     +     riverflow, riverstage, ict, qct, dataprovider, 
     +     handbook5id, homewfo, precip12hrdd, precip1hrdd, 
     +     precip24hrdd, precip3hrdd, precip5mindd, precip6hrdd, 
     +     precipaccumdd, providerid, rawmessage, staticids, 
     +     stationid, stationname, stationtype, observationtime, 
     +     receivedtime, reporttime, riverreportchangetime, 
     +     filtersetnum, firstinbin, firstoverflow, globalinventory, 
     +     invtime, inventory, isoverflow, lastinbin, lastrecord, 
     +     nstaticids, numericwmoid, precip12hrica, precip12hricr, 
     +     precip12hrqca, precip12hrqcr, precip1hrica, precip1hricr, 
     +     precip1hrqca, precip1hrqcr, precip24hrica, precip24hricr, 
     +     precip24hrqca, precip24hrqcr, precip3hrica, precip3hricr, 
     +     precip3hrqca, precip3hrqcr, precip5minica, precip5minicr, 
     +     precip5minqca, precip5minqcr, precip6hrica, precip6hricr, 
     +     precip6hrqca, precip6hrqcr, precipaccumica, 
     +     precipaccumicr, precipaccumqca, precipaccumqcr, 
     +     prevrecord, secondsstage1_2, secondsstage3,badflag)
c
      include 'netcdf.inc'
      integer icchecknum, qcchecknum, maxstaticids, ninventorybins, 
     +     recnum,nf_fid, nf_vid, nf_status
      integer filtersetnum, firstinbin(ninventorybins), firstoverflow,
     +     globalinventory, invtime(recnum), inventory(maxstaticids),
     +     isoverflow(recnum), lastinbin(ninventorybins),
     +     lastrecord(maxstaticids), nstaticids,
     +     numericwmoid(recnum), precip12hrica(recnum),
     +     precip12hricr(recnum), precip12hrqca(recnum),
     +     precip12hrqcr(recnum), precip1hrica(recnum),
     +     precip1hricr(recnum), precip1hrqca(recnum),
     +     precip1hrqcr(recnum), precip24hrica(recnum),
     +     precip24hricr(recnum), precip24hrqca(recnum),
     +     precip24hrqcr(recnum), precip3hrica(recnum),
     +     precip3hricr(recnum), precip3hrqca(recnum),
     +     precip3hrqcr(recnum), precip5minica(recnum),
     +     precip5minicr(recnum), precip5minqca(recnum),
     +     precip5minqcr(recnum), precip6hrica(recnum),
     +     precip6hricr(recnum), precip6hrqca(recnum),
     +     precip6hrqcr(recnum), precipaccumica(recnum),
     +     precipaccumicr(recnum), precipaccumqca(recnum),
     +     precipaccumqcr(recnum), prevrecord(recnum),
     +     secondsstage1_2(recnum), secondsstage3(recnum)
      real elevation(recnum), latitude(recnum), longitude(recnum),
     +     precip12hr(recnum), precip12hrqcd( qcchecknum, recnum),
     +     precip1hr(recnum), precip1hrqcd( qcchecknum, recnum),
     +     precip24hr(recnum), precip24hrqcd( qcchecknum, recnum),
     +     precip3hr(recnum), precip3hrqcd( qcchecknum, recnum),
     +     precip5min(recnum), precip5minqcd( qcchecknum, recnum),
     +     precip6hr(recnum), precip6hrqcd( qcchecknum, recnum),
     +     precipaccum(recnum), precipaccumqcd( qcchecknum, recnum),
     +     riverflow(recnum), riverstage(recnum)
      double precision observationtime(recnum), receivedtime(recnum),
     +     reporttime(recnum), riverreportchangetime(recnum)
      character precip5mindd(recnum)
      character precip12hrdd(recnum)
      character*72 ict(icchecknum)
      character precip1hrdd(recnum)
      character*11 dataprovider(recnum)
      character*60 qct(qcchecknum)
      character*11 handbook5id(recnum)
      character precip24hrdd(recnum)
      character precipaccumdd(recnum)
      character*51 stationname(recnum)
      character precip3hrdd(recnum)
      character*11 stationtype(recnum)
      character precip6hrdd(recnum)
      character*256 rawmessage(recnum)
      character*24 staticids(maxstaticids)
      character*12 providerid(recnum)
      character*11 stationid(recnum)
      character*4 homewfo(recnum)

      write(6,*)' subroutine read_madis_hydro_netcdf: recnum = ',recnum

c   variables of type real
c
c     variable        netcdf long name
c     elevation     "elevation"
c
      nf_status=nf_inq_varid(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for elevation'
       print *,'set elevation to badflag'
       elevation = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,elevation)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for elevation'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - elevation'
       else
        call ck_array_real(elevation,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - elevation'
       else
        call ck_array_real(elevation,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     latitude      "latitude"
c
      nf_status=nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitude'
       print *,'set latitude to badflag'
       latitude = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,latitude)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitude'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - latitude'
       else
        call ck_array_real(latitude,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - latitude'
       else
        call ck_array_real(latitude,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     longitude     "longitude"
c
      nf_status=nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for longitude'
       print *,'set longitude to badflag'
       longitude = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,longitude)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitude'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - longitude'
       else
        call ck_array_real(longitude,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - longitude'
       else
        call ck_array_real(longitude,recnum,valmis
     1                    ,badflag)
       endif
      endif

c
c     variable        netcdf long name
c     precip12hr    "12 hour precip accumulation "
c
      nf_status=nf_inq_varid(nf_fid,'precip12hr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip12hr'
       print *,'set precip12hr to badflag'
       precip12hr = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip12hr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip12hr'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip12hr'
       else
        call ck_array_real(precip12hr,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip12hr'
       else
        call ck_array_real(precip12hr,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip12hrqcd "12-hr precip amount qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'precip12hrqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip12hrqcd'
       print *,'set precip12hrqcd to badflag'
       precip12hrqcd = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip12hrqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip12hrqcd'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip12hrqcd'
       else
        call ck_array_real(precip12hrqcd,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip12hrqcd'
       else
        call ck_array_real(precip12hrqcd,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip1hr     "1 hour precip accumulation "
c
      nf_status=nf_inq_varid(nf_fid,'precip1hr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip1hr'
       print *,'set precip1hr to badflag'
       precip1hr = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip1hr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip1hr'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip1hr'
       else
        call ck_array_real(precip1hr,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip1hr'
       else
        call ck_array_real(precip1hr,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip1hrqcd  "1-hr precip amount qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'precip1hrqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip1hrqcd'
       print *,'set precip1hrqcd to badflag'
       precip1hrqcd = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip1hrqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip1hrqcd'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip1hrqcd'
       else
        call ck_array_real(precip1hrqcd,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip1hrqcd'
       else
        call ck_array_real(precip1hrqcd,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip24hr    "24 hour precip accumulation "
c
      nf_status=nf_inq_varid(nf_fid,'precip24hr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip24hr'
       print *,'set precip24hr to badflag'
       precip24hr = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip24hr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip24hr'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip24hr'
       else
        call ck_array_real(precip24hr,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip24hr'
       else
        call ck_array_real(precip24hr,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip24hrqcd "24-hr precip amount qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'precip24hrqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip24hrqcd'
       print *,'set precip24hrqcd to badflag'
       precip24hrqcd = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip24hrqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip24hrqcd'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip24hrqcd'
       else
        call ck_array_real(precip24hrqcd,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip24hrqcd'
       else
        call ck_array_real(precip24hrqcd,recnum,valmis
     1                    ,badflag)
       endif
      endif

c
c     variable        netcdf long name
c     precip3hr     "3 hour precip accumulation "
c
      nf_status=nf_inq_varid(nf_fid,'precip3hr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip3hr'
       print *,'set precip3hr to badflag'
       precip3hr = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip3hr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip3hr'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip3hr'
       else
        call ck_array_real(precip3hr,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip3hr'
       else
        call ck_array_real(precip3hr,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip3hrqcd  "3-hr precip amount qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'precip3hrqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip3hrqcd'
       print *,'set precip3hrqcd to badflag'
       precip3hrqcd = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip3hrqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip3hrqcd'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip3hrqcd'
       else
        call ck_array_real(precip3hrqcd,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip3hrqcd'
       else
        call ck_array_real(precip3hrqcd,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip5min    "5 minute precip accumulation "
c
      nf_status=nf_inq_varid(nf_fid,'precip5min',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip5min'
       print *,'set precip5min to badflag'
       precip5min = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip5min)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip5min'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip5min'
       else
        call ck_array_real(precip5min,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip5min'
       else
        call ck_array_real(precip5min,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip5minqcd "5-min precip amount qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'precip5minqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip5minqcd'
       print *,'set precip5minqcd to badflag'
       precip5minqcd = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip5minqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip5minqcd'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip5minqcd'
       else
        call ck_array_real(precip5minqcd,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip5minqcd'
       else
        call ck_array_real(precip5minqcd,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip6hr     "6 hour precip accumulation "
c
      nf_status=nf_inq_varid(nf_fid,'precip6hr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip6hr'
       print *,'set precip6hr to badflag'
       precip6hr = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip6hr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip6hr'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip6hr'
       else
        call ck_array_real(precip6hr,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip6hr'
       else
        call ck_array_real(precip6hr,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precip6hrqcd  "6-hr precip amount qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'precip6hrqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip6hrqcd'
       print *,'set precip6hrqcd to badflag'
       precip6hrqcd = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precip6hrqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip6hrqcd'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precip6hrqcd'
       else
        call ck_array_real(precip6hrqcd,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precip6hrqcd'
       else
        call ck_array_real(precip6hrqcd,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precipaccum   "precip accumulation with an unknown time period"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccum',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccum'
       print *,'set precipaccum to badflag'
       precipaccum = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precipaccum)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccum'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precipaccum'
       else
        call ck_array_real(precipaccum,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precipaccum'
       else
        call ck_array_real(precipaccum,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     precipaccumqcd"precip amount qc departures"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccumqcd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccumqcd'
       print *,'set precipaccumqcd to badflag'
       precipaccumqcd = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,precipaccumqcd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccumqcd'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - precipaccumqcd'
       else
        call ck_array_real(precipaccumqcd,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - precipaccumqcd'
       else
        call ck_array_real(precipaccumqcd,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     riverflow     "river flow"
c
      nf_status=nf_inq_varid(nf_fid,'riverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for riverflow'
       print *,'set riverflow to badflag'
       riverflow = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,riverflow)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for riverflow'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - riverflow'
       else
        call ck_array_real(riverflow,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - riverflow'
       else
        call ck_array_real(riverflow,recnum,valmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     riverstage    "river stage"
c
      nf_status=nf_inq_varid(nf_fid,'riverstage',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for riverstage'
       print *,'set riverstage to badflag'
       riverstage = badflag
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,riverstage)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for riverstage'
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'_fillvalue',valfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - riverstage'
       else
        call ck_array_real(riverstage,recnum,valfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_real(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - riverstage'
       else
        call ck_array_real(riverstage,recnum,valmis
     1                    ,badflag)
       endif
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c     filtersetnum  
c
      nf_status=nf_inq_varid(nf_fid,'filtersetnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for filtersetnum'
       print *,'set filtersetnum to -99'
       filtersetnum = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,filtersetnum)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for filtersetnum'
       endif
      endif
c
c     variable        netcdf long name
c     firstinbin    
c
      nf_status=nf_inq_varid(nf_fid,'firstinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for firstinbin'
       print *,'set firstinbin to -99'
       firstinbin = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,firstinbin)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for firstinbin'
       endif
      endif
c
c     variable        netcdf long name
c     firstoverflow 
c
      nf_status=nf_inq_varid(nf_fid,'firstoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for firstoverflow'
       print *,'set firstoverflow to -99'
       firstoverflow = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,firstoverflow)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for firstoverflow'
       endif
      endif
c
c     variable        netcdf long name
c     globalinventory
c
      nf_status=nf_inq_varid(nf_fid,'globalinventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for globalinventory'
       print *,'set globalinventory to -99'
       globalinventory = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,globalinventory)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for globalinventory'
       endif
      endif
c
c     variable        netcdf long name
c     invtime       
c
      nf_status=nf_inq_varid(nf_fid,'invtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for invtime'
       print *,'set invtime to -99'
       invtime = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,invtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for invtime'
       endif
      endif
c
c     variable        netcdf long name
c     inventory     
c
      nf_status=nf_inq_varid(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for inventory'
       print *,'set inventory to -99'
       inventory = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,inventory)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for inventory'
       endif
      endif
c
c     variable        netcdf long name
c     isoverflow    
c
      nf_status=nf_inq_varid(nf_fid,'isoverflow',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for isoverflow'
       print *,'set isoverflow to -99'
       isoverflow = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,isoverflow)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for isoverflow'
       endif
      endif
c
c     variable        netcdf long name
c     lastinbin     
c
      nf_status=nf_inq_varid(nf_fid,'lastinbin',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for lastinbin'
       print *,'set lastinbin to -99'
       lastinbin = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,lastinbin)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for lastinbin'
       endif
      endif
c
c     variable        netcdf long name
c     lastrecord    
c
      nf_status=nf_inq_varid(nf_fid,'lastrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for lastrecord'
       print *,'set lastrecord to -99'
       lastrecord = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,lastrecord)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for lastrecord'
       endif
      endif
c
c     variable        netcdf long name
c     nstaticids    
c
      nf_status=nf_inq_varid(nf_fid,'nstaticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for nstaticids'
       print *,'set nstaticids to -99'
       nstaticids = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,nstaticids)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for nstaticids'
       endif
      endif

c
c     variable        netcdf long name
c     numericwmoid  "numeric wmo identification"
c
      nf_status=nf_inq_varid(nf_fid,'numericwmoid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for numericwmoid'
       print *,'set numericwmoid to -99'
       numericwmoid = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,numericwmoid)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for numericwmoid'
       endif
      endif
c
c     variable        netcdf long name
c     precip12hrica "12-hr precip amount ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip12hrica',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip12hrica'
       print *,'set precip12hrica to -99'
       precip12hrica = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip12hrica)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip12hrica'
       endif
      endif
c
c     variable        netcdf long name
c     precip12hricr "12-hr precip amount ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip12hricr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip12hricr'
       print *,'set precip12hricr to -99'
       precip12hricr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip12hricr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip12hricr'
       endif
      endif
c
c     variable        netcdf long name
c     precip12hrqca "12-hr precip amount qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip12hrqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip12hrqca'
       print *,'set precip12hrqca to -99'
       precip12hrqca = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip12hrqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip12hrqca'
       endif
      endif
c
c     variable        netcdf long name
c     precip12hrqcr "12-hr precip amount qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip12hrqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip12hrqcr'
       print *,'set precip12hrqcr to -99'
       precip12hrqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip12hrqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip12hrqcr'
       endif
      endif
c
c     variable        netcdf long name
c     precip1hrica  "1-hr precip amount ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip1hrica',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip1hrica'
       print *,'set precip1hrica to -99'
       precip1hrica = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip1hrica)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip1hrica'
       endif
      endif
c
c     variable        netcdf long name
c     precip1hricr  "1-hr precip amount ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip1hricr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip1hricr'
       print *,'set precip1hricr to -99'
       precip1hricr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip1hricr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip1hricr'
       endif
      endif
c
c     variable        netcdf long name
c     precip1hrqca  "1-hr precip amount qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip1hrqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip1hrqca'
       print *,'set precip1hrqca to -99'
       precip1hrqca = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip1hrqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip1hrqca'
       endif
      endif
c
c     variable        netcdf long name
c     precip1hrqcr  "1-hr precip amount qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip1hrqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip1hrqcr'
       print *,'set precip1hrqcr to -99'
       precip1hrqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip1hrqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip1hrqcr'
       endif
      endif
c
c     variable        netcdf long name
c     precip24hrica "24-hr precip amount ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip24hrica',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip24hrica'
       print *,'set precip24hrica to -99'
       precip24hrica = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip24hrica)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip24hrica'
       endif
      endif
c
c     variable        netcdf long name
c     precip24hricr "24-hr precip amount ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip24hricr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip24hricr'
       print *,'set precip24hricr to -99'
       precip24hricr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip24hricr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip24hricr'
       endif
      endif
c
c     variable        netcdf long name
c     precip24hrqca "24-hr precip amount qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip24hrqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip24hrqca'
       print *,'set precip24hrqca to -99'
       precip24hrqca = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip24hrqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip24hrqca'
       endif
      endif
c
c     variable        netcdf long name
c     precip24hrqcr "24-hr precip amount qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip24hrqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip24hrqcr'
       print *,'set precip24hrqcr to -99'
       precip24hrqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip24hrqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip24hrqcr'
       endif
      endif
c
c     variable        netcdf long name
c     precip3hrica  "3-hr precip amount ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip3hrica',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip3hrica'
       print *,'set precip3hrica to -99'
       precip3hrica = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip3hrica)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip3hrica'
       endif
      endif
c
c     variable        netcdf long name
c     precip3hricr  "3-hr precip amount ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip3hricr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip3hricr'
       print *,'set precip3hricr to -99'
       precip3hricr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip3hricr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip3hricr'
       endif
      endif
c
c     variable        netcdf long name
c     precip3hrqca  "3-hr precip amount qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip3hrqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip3hrqca'
       print *,'set precip3hrqca to -99'
       precip3hrqca = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip3hrqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip3hrqca'
       endif
      endif
c
c     variable        netcdf long name
c     precip3hrqcr  "3-hr precip amount qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip3hrqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip3hrqcr'
       print *,'set precip3hrqcr to -99'
       precip3hrqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip3hrqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip3hrqcr'
       endif
      endif
c
c     variable        netcdf long name
c     precip5minica "5-min precip amount ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip5minica',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip5minica'
       print *,'set precip5minica to -99'
       precip5minica = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip5minica)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip5minica'
       endif
      endif
c
c     variable        netcdf long name
c     precip5minicr "5-min precip amount ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip5minicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip5minicr'
       print *,'set precip5minicr to -99'
       precip5minicr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip5minicr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip5minicr'
       endif
      endif
c
c     variable        netcdf long name
c     precip5minqca "5-min precip amount qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip5minqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip5minqca'
       print *,'set precip5minqca to -99'
       precip5minqca = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip5minqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip5minqca'
       endif
      endif
c
c     variable        netcdf long name
c     precip5minqcr "5-min precip amount qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip5minqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip5minqcr'
       print *,'set precip5minqcr to -99'
       precip5minqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip5minqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip5minqcr'
       endif
      endif
c
c     variable        netcdf long name
c     precip6hrica  "6-hr precip amount ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip6hrica',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip6hrica'
       print *,'set precip6hrica to -99'
       precip6hrica = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip6hrica)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip6hrica'
       endif
      endif
c
c     variable        netcdf long name
c     precip6hricr  "6-hr precip amount ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip6hricr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip6hricr'
       print *,'set precip6hricr to -99'
       precip6hricr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip6hricr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip6hricr'
       endif
      endif
c
c     variable        netcdf long name
c     precip6hrqca  "6-hr precip amount qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precip6hrqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip6hrqca'
       print *,'set precip6hrqca to -99'
       precip6hrqca = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip6hrqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip6hrqca'
       endif
      endif
c
c     variable        netcdf long name
c     precip6hrqcr  "6-hr precip amount qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'precip6hrqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip6hrqcr'
       print *,'set precip6hrqcr to -99'
       precip6hrqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precip6hrqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip6hrqcr'
       endif
      endif
c
c     variable        netcdf long name
c     precipaccumica"precip amount ic applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccumica',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccumica'
       print *,'set precipaccumica to -99'
       precipaccumica = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precipaccumica)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccumica'
       endif
      endif
c
c     variable        netcdf long name
c     precipaccumicr"precip amount ic results word"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccumicr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccumicr'
       print *,'set precipaccumicr to -99'
       precipaccumicr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precipaccumicr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccumicr'
       endif
      endif
c
c     variable        netcdf long name
c     precipaccumqca"precip amount qc applied word"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccumqca',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccumqca'
       print *,'set precipaccumqca to -99'
       precipaccumqca = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precipaccumqca)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccumqca'
       endif
      endif
c
c     variable        netcdf long name
c     precipaccumqcr"precip amount qc results word"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccumqcr',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccumqcr'
       print *,'set precipaccumqcr to -99'
       precipaccumqcr = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,precipaccumqcr)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccumqcr'
       endif
      endif

c
c     variable        netcdf long name
c     prevrecord    
c
      nf_status=nf_inq_varid(nf_fid,'prevrecord',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for prevrecord'
       print *,'set prevrecord to -99'
       prevrecord = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,prevrecord)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for prevrecord'
       endif
      endif
c
c     variable        netcdf long name
c     secondsstage1_2
c
      nf_status=nf_inq_varid(nf_fid,'secondsstage1_2',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for secondsstage1_2'
       print *,'set secondsstage1_2 to -99'
       secondsstage1_2 = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,secondsstage1_2)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for secondsstage1_2'
       endif
      endif
c
c     variable        netcdf long name
c     secondsstage3 
c
      nf_status=nf_inq_varid(nf_fid,'secondsstage3',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for secondsstage3'
       print *,'set secondsstage3 to -99'
       secondsstage3 = -99
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,secondsstage3)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for secondsstage3'
       endif
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c     observationtime"time of observation"
c
      nf_status=nf_inq_varid(nf_fid,'observationtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for observationtime'
       print *,'set observationtime to badflag'
       observationtime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,observationtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for observationtime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - observationtime'
       else
        call ck_array_dble(observationtime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - observationtime'
       else
        call ck_array_dble(observationtime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     receivedtime  "time data was received"
c
      nf_status=nf_inq_varid(nf_fid,'receivedtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for receivedtime'
       print *,'set receivedtime to badflag'
       receivedtime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,receivedtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for receivedtime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - receivedtime'
       else
        call ck_array_dble(receivedtime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - receivedtime'
       else
        call ck_array_dble(receivedtime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     reporttime    "time data was processed by the provider"
c
      nf_status=nf_inq_varid(nf_fid,'reporttime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for reporttime'
       print *,'set reporttime to badflag'
       reporttime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,reporttime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for reporttime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - reporttime'
       else
        call ck_array_dble(reporttime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - reporttime'
       else
        call ck_array_dble(reporttime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif
c
c     variable        netcdf long name
c     riverreportchangetime"time of last new river stage/flow rpt"
c
      nf_status=nf_inq_varid(nf_fid,'riverreportchangetime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for riverreportchangetime'
       print *,'set riverreportchangetime to badflag'
       riverreportchangetime = badflag
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,riverreportchangetime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for riverreportchangetime'
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'_fillvalue',dvalfil)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for _fillvalue - riverreportchangetime'
       else
        call ck_array_dble(riverreportchangetime,recnum,dvalfil
     1                    ,badflag)
       endif
       nf_status=nf_get_att_double(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
     1         ,' for missing_value - riverreportchangetime'
       else
        call ck_array_dble(riverreportchangetime,recnum,dvalmis
     1                    ,badflag)
       endif
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c     ict           "list of possible ic checks"
c
      nf_status=nf_inq_varid(nf_fid,'ict',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for ict'
       print *,'set ict to " "'
       ict = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,ict)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for ict'
       endif
      endif
c
c     variable        netcdf long name
c     qct           "list of possible qc checks"
c
      nf_status=nf_inq_varid(nf_fid,'qct',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for qct'
       print *,'set qct to " "'
       qct = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,qct)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for qct'
       endif
      endif
c
c     variable        netcdf long name
c     dataprovider  "ldad data provider"
c
      nf_status=nf_inq_varid(nf_fid,'dataprovider',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dataprovider'
       print *,'set dataprovider to " "'
       dataprovider = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,dataprovider)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dataprovider'
       endif
      endif
c
c     variable        netcdf long name
c     handbook5id   "handbook5 id (afos or shef id)"
c
      nf_status=nf_inq_varid(nf_fid,'handbook5id',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for handbook5id'
       print *,'set handbook5id to " "'
       handbook5id = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,handbook5id)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for handbook5id'
       endif
      endif
c
c     variable        netcdf long name
c     homewfo       "home wfo id"
c
      nf_status=nf_inq_varid(nf_fid,'homewfo',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for homewfo'
       print *,'set homewfo to " "'
       homewfo = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,homewfo)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for homewfo'
       endif
      endif
c
c     variable        netcdf long name
c     precip12hrdd  "12-hr precip amount qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precip12hrdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip12hrdd'
       print *,'set precip12hrdd to " "'
       precip12hrdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precip12hrdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip12hrdd'
       endif
      endif
c
c     variable        netcdf long name
c     precip1hrdd   "1-hr precip amount qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precip1hrdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip1hrdd'
       print *,'set precip1hrdd to " "'
       precip1hrdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precip1hrdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip1hrdd'
       endif
      endif
c
c     variable        netcdf long name
c     precip24hrdd  "24-hr precip amount qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precip24hrdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip24hrdd'
       print *,'set precip24hrdd to " "'
       precip24hrdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precip24hrdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip24hrdd'
       endif
      endif
c
c     variable        netcdf long name
c     precip3hrdd   "3-hr precip amount qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precip3hrdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip3hrdd'
       print *,'set precip3hrdd to " "'
       precip3hrdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precip3hrdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip3hrdd'
       endif
      endif
c
c     variable        netcdf long name
c     precip5mindd  "5-min precip amount qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precip5mindd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip5mindd'
       print *,'set precip5mindd to " "'
       precip5mindd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precip5mindd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip5mindd'
       endif
      endif
c
c     variable        netcdf long name
c     precip6hrdd   "6-hr precip amount qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precip6hrdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precip6hrdd'
       print *,'set precip6hrdd to " "'
       precip6hrdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precip6hrdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precip6hrdd'
       endif
      endif
c
c     variable        netcdf long name
c     precipaccumdd "precip amount qc summary value"
c
      nf_status=nf_inq_varid(nf_fid,'precipaccumdd',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for precipaccumdd'
       print *,'set precipaccumdd to " "'
       precipaccumdd = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,precipaccumdd)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for precipaccumdd'
       endif
      endif
c
c     variable        netcdf long name
c     providerid    "data provider station id"
c
      nf_status=nf_inq_varid(nf_fid,'providerid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for providerid'
       print *,'set providerid to " "'
       providerid = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,providerid)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for providerid'
       endif
      endif
c
c     variable        netcdf long name
c     rawmessage    "raw text ldad hydro report"
c
      nf_status=nf_inq_varid(nf_fid,'rawmessage',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for rawmessage'
       print *,'set rawmessage to " "'
       rawmessage = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,rawmessage)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for rawmessage'
       endif
      endif
c
c     variable        netcdf long name
c     staticids     
c
      nf_status=nf_inq_varid(nf_fid,'staticids',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for staticids'
       print *,'set staticids to " "'
       staticids = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,staticids)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for staticids'
       endif
      endif
c
c     variable        netcdf long name
c     stationid     "alphanumeric station id"
c
      nf_status=nf_inq_varid(nf_fid,'stationid',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationid'
       print *,'set stationid to " "'
       stationid = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,stationid)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationid'
       endif
      endif
c
c     variable        netcdf long name
c     stationname   "alphanumeric station name"
c
      nf_status=nf_inq_varid(nf_fid,'stationname',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationname'
       print *,'set stationname to " "'
       stationname = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,stationname)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationname'
       endif
      endif
c
c     variable        netcdf long name
c     stationtype   "ldad station type"
c
      nf_status=nf_inq_varid(nf_fid,'stationtype',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for stationtype'
       print *,'set stationtype to " "'
       stationtype = ' '
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,stationtype)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for stationtype'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
