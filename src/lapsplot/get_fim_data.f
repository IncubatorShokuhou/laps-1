      subroutine get_fim_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,pres_p
     +                   ,clwc_p
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer latitude, longitude, time,nf_fid, nf_vid, nf_status
      real pres_p(nx_l,ny_l,21)
      real clwc_p(nx_l,ny_l,21)
c
c  open netcdf file for reading
c
      nf_status=nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),filename
        istatus=0
        return
      endif
c
c  fill all dimension values
c
c
c get size of latitude
c
      nf_status=nf_inq_dimid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim latitude'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,latitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim latitude'
      endif
c
c get size of longitude
c
      nf_status=nf_inq_dimid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim longitude'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,longitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim longitude'
      endif
c
c get size of time
c
      nf_status=nf_inq_dimid(nf_fid,'time',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim time'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,time)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim time'
      endif
      call read_fim_data(nf_fid, latitude, longitude, time,
     +     i4time_sys, ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, pres_p, clwc_p, istatus)

      write(6,*)' dimensions: ',nx_l,ny_l,longitude,latitude

      write(6,*)' range of clwc_p (kg/m^3)',minval(clwc_p)
     +                                     ,maxval(clwc_p)

      return
      end
c
c
      subroutine read_fim_data(nf_fid, latitude, longitude, time,
     +     i4time_sys, ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, pres_p, clwc_p, istatus)


      use constants_laps, only: r

      include 'netcdf.inc'
      integer latitude, longitude, time,nf_fid, nf_vid, nf_status

      real pres_p(longitude,latitude,21)
      real clwc_p(longitude,latitude,21)
      real rho_p(longitude,latitude,21)
      real buff_p(longitude,latitude,21)

      real pres(longitude,latitude,65)
      real clwc(longitude,latitude,65)

      real clwmr_10hybridlevel( longitude,  latitude, time),
     +     clwmr_11hybridlevel( longitude,  latitude, time),
     +     clwmr_12hybridlevel( longitude,  latitude, time),
     +     clwmr_13hybridlevel( longitude,  latitude, time),
     +     clwmr_14hybridlevel( longitude,  latitude, time),
     +     clwmr_15hybridlevel( longitude,  latitude, time),
     +     clwmr_16hybridlevel( longitude,  latitude, time),
     +     clwmr_17hybridlevel( longitude,  latitude, time),
     +     clwmr_18hybridlevel( longitude,  latitude, time),
     +     clwmr_19hybridlevel( longitude,  latitude, time),
     +     clwmr_1hybridlevel( longitude,  latitude, time),
     +     clwmr_20hybridlevel( longitude,  latitude, time),
     +     clwmr_21hybridlevel( longitude,  latitude, time),
     +     clwmr_22hybridlevel( longitude,  latitude, time),
     +     clwmr_23hybridlevel( longitude,  latitude, time),
     +     clwmr_24hybridlevel( longitude,  latitude, time),
     +     clwmr_25hybridlevel( longitude,  latitude, time),
     +     clwmr_26hybridlevel( longitude,  latitude, time),
     +     clwmr_27hybridlevel( longitude,  latitude, time),
     +     clwmr_28hybridlevel( longitude,  latitude, time),
     +     clwmr_29hybridlevel( longitude,  latitude, time),
     +     clwmr_2hybridlevel( longitude,  latitude, time),
     +     clwmr_30hybridlevel( longitude,  latitude, time),
     +     clwmr_31hybridlevel( longitude,  latitude, time),
     +     clwmr_32hybridlevel( longitude,  latitude, time),
     +     clwmr_33hybridlevel( longitude,  latitude, time),
     +     clwmr_34hybridlevel( longitude,  latitude, time),
     +     clwmr_35hybridlevel( longitude,  latitude, time),
     +     clwmr_36hybridlevel( longitude,  latitude, time),
     +     clwmr_37hybridlevel( longitude,  latitude, time),
     +     clwmr_38hybridlevel( longitude,  latitude, time),
     +     clwmr_39hybridlevel( longitude,  latitude, time),
     +     clwmr_3hybridlevel( longitude,  latitude, time),
     +     clwmr_40hybridlevel( longitude,  latitude, time),
     +     clwmr_41hybridlevel( longitude,  latitude, time),
     +     clwmr_42hybridlevel( longitude,  latitude, time),
     +     clwmr_43hybridlevel( longitude,  latitude, time),
     +     clwmr_44hybridlevel( longitude,  latitude, time),
     +     clwmr_45hybridlevel( longitude,  latitude, time),
     +     clwmr_46hybridlevel( longitude,  latitude, time),
     +     clwmr_47hybridlevel( longitude,  latitude, time),
     +     clwmr_48hybridlevel( longitude,  latitude, time),
     +     clwmr_49hybridlevel( longitude,  latitude, time),
     +     clwmr_4hybridlevel( longitude,  latitude, time),
     +     clwmr_50hybridlevel( longitude,  latitude, time),
     +     clwmr_51hybridlevel( longitude,  latitude, time),
     +     clwmr_52hybridlevel( longitude,  latitude, time),
     +     clwmr_53hybridlevel( longitude,  latitude, time),
     +     clwmr_54hybridlevel( longitude,  latitude, time),
     +     clwmr_55hybridlevel( longitude,  latitude, time),
     +     clwmr_56hybridlevel( longitude,  latitude, time),
     +     clwmr_57hybridlevel( longitude,  latitude, time),
     +     clwmr_58hybridlevel( longitude,  latitude, time),
     +     clwmr_59hybridlevel( longitude,  latitude, time),
     +     clwmr_5hybridlevel( longitude,  latitude, time),
     +     clwmr_60hybridlevel( longitude,  latitude, time),
     +     clwmr_61hybridlevel( longitude,  latitude, time),
     +     clwmr_62hybridlevel( longitude,  latitude, time),
     +     clwmr_63hybridlevel( longitude,  latitude, time),
     +     clwmr_64hybridlevel( longitude,  latitude, time),
     +     clwmr_6hybridlevel( longitude,  latitude, time),
     +     clwmr_7hybridlevel( longitude,  latitude, time),
     +     clwmr_8hybridlevel( longitude,  latitude, time),
     +     clwmr_9hybridlevel( longitude,  latitude, time),
     +     pres_10hybridlevel( longitude,  latitude, time),
     +     pres_11hybridlevel( longitude,  latitude, time),
     +     pres_12hybridlevel( longitude,  latitude, time),
     +     pres_13hybridlevel( longitude,  latitude, time),
     +     pres_14hybridlevel( longitude,  latitude, time),
     +     pres_15hybridlevel( longitude,  latitude, time),
     +     pres_16hybridlevel( longitude,  latitude, time),
     +     pres_17hybridlevel( longitude,  latitude, time),
     +     pres_18hybridlevel( longitude,  latitude, time),
     +     pres_19hybridlevel( longitude,  latitude, time),
     +     pres_1hybridlevel( longitude,  latitude, time),
     +     pres_20hybridlevel( longitude,  latitude, time),
     +     pres_21hybridlevel( longitude,  latitude, time),
     +     pres_22hybridlevel( longitude,  latitude, time),
     +     pres_23hybridlevel( longitude,  latitude, time),
     +     pres_24hybridlevel( longitude,  latitude, time),
     +     pres_25hybridlevel( longitude,  latitude, time),
     +     pres_26hybridlevel( longitude,  latitude, time),
     +     pres_27hybridlevel( longitude,  latitude, time),
     +     pres_28hybridlevel( longitude,  latitude, time),
     +     pres_29hybridlevel( longitude,  latitude, time),
     +     pres_2hybridlevel( longitude,  latitude, time),
     +     pres_30hybridlevel( longitude,  latitude, time),
     +     pres_31hybridlevel( longitude,  latitude, time),
     +     pres_32hybridlevel( longitude,  latitude, time),
     +     pres_33hybridlevel( longitude,  latitude, time),
     +     pres_34hybridlevel( longitude,  latitude, time),
     +     pres_35hybridlevel( longitude,  latitude, time),
     +     pres_36hybridlevel( longitude,  latitude, time),
     +     pres_37hybridlevel( longitude,  latitude, time),
     +     pres_38hybridlevel( longitude,  latitude, time),
     +     pres_39hybridlevel( longitude,  latitude, time),
     +     pres_3hybridlevel( longitude,  latitude, time),
     +     pres_40hybridlevel( longitude,  latitude, time),
     +     pres_41hybridlevel( longitude,  latitude, time),
     +     pres_42hybridlevel( longitude,  latitude, time),
     +     pres_43hybridlevel( longitude,  latitude, time),
     +     pres_44hybridlevel( longitude,  latitude, time),
     +     pres_45hybridlevel( longitude,  latitude, time),
     +     pres_46hybridlevel( longitude,  latitude, time),
     +     pres_47hybridlevel( longitude,  latitude, time),
     +     pres_48hybridlevel( longitude,  latitude, time),
     +     pres_49hybridlevel( longitude,  latitude, time),
     +     pres_4hybridlevel( longitude,  latitude, time),
     +     pres_50hybridlevel( longitude,  latitude, time),
     +     pres_51hybridlevel( longitude,  latitude, time),
     +     pres_52hybridlevel( longitude,  latitude, time),
     +     pres_53hybridlevel( longitude,  latitude, time),
     +     pres_54hybridlevel( longitude,  latitude, time),
     +     pres_55hybridlevel( longitude,  latitude, time),
     +     pres_56hybridlevel( longitude,  latitude, time),
     +     pres_57hybridlevel( longitude,  latitude, time),
     +     pres_58hybridlevel( longitude,  latitude, time),
     +     pres_59hybridlevel( longitude,  latitude, time),
     +     pres_5hybridlevel( longitude,  latitude, time),
     +     pres_60hybridlevel( longitude,  latitude, time),
     +     pres_61hybridlevel( longitude,  latitude, time),
     +     pres_62hybridlevel( longitude,  latitude, time),
     +     pres_63hybridlevel( longitude,  latitude, time),
     +     pres_64hybridlevel( longitude,  latitude, time),
     +     pres_65hybridlevel( longitude,  latitude, time),
     +     pres_6hybridlevel( longitude,  latitude, time),
     +     pres_7hybridlevel( longitude,  latitude, time),
     +     pres_8hybridlevel( longitude,  latitude, time),
     +     pres_9hybridlevel( longitude,  latitude, time)
!     double precision latitude(latitude), longitude(longitude),
!    +     time(time)

!     declarations for 'write_zzz' call
!     integer iwmostanum(recnum)
!     character a9time_ob_r(recnum)*9
      logical l_closest_time, l_closest_time_i, l_in_domain
      real*4 lat_a(nx_l,ny_l)
      real*4 lon_a(nx_l,ny_l)
      real*4 topo_a(nx_l,ny_l)

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting r_missing_data'
          return
      endif
      call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_domain_perimeter'
          return
      endif

      write(6,*)' call read_fim_netcdf ',latitude,longitude,time

      call read_fim_netcdf(nf_fid, ! latitude, longitude, time, 
     +     clwmr_10hybridlevel, clwmr_11hybridlevel, 
     +     clwmr_12hybridlevel, clwmr_13hybridlevel, 
     +     clwmr_14hybridlevel, clwmr_15hybridlevel, 
     +     clwmr_16hybridlevel, clwmr_17hybridlevel, 
     +     clwmr_18hybridlevel, clwmr_19hybridlevel, 
     +     clwmr_1hybridlevel, clwmr_20hybridlevel, 
     +     clwmr_21hybridlevel, clwmr_22hybridlevel, 
     +     clwmr_23hybridlevel, clwmr_24hybridlevel, 
     +     clwmr_25hybridlevel, clwmr_26hybridlevel, 
     +     clwmr_27hybridlevel, clwmr_28hybridlevel, 
     +     clwmr_29hybridlevel, clwmr_2hybridlevel, 
     +     clwmr_30hybridlevel, clwmr_31hybridlevel, 
     +     clwmr_32hybridlevel, clwmr_33hybridlevel, 
     +     clwmr_34hybridlevel, clwmr_35hybridlevel, 
     +     clwmr_36hybridlevel, clwmr_37hybridlevel, 
     +     clwmr_38hybridlevel, clwmr_39hybridlevel, 
     +     clwmr_3hybridlevel, clwmr_40hybridlevel, 
     +     clwmr_41hybridlevel, clwmr_42hybridlevel, 
     +     clwmr_43hybridlevel, clwmr_44hybridlevel, 
     +     clwmr_45hybridlevel, clwmr_46hybridlevel, 
     +     clwmr_47hybridlevel, clwmr_48hybridlevel, 
     +     clwmr_49hybridlevel, clwmr_4hybridlevel, 
     +     clwmr_50hybridlevel, clwmr_51hybridlevel, 
     +     clwmr_52hybridlevel, clwmr_53hybridlevel, 
     +     clwmr_54hybridlevel, clwmr_55hybridlevel, 
     +     clwmr_56hybridlevel, clwmr_57hybridlevel, 
     +     clwmr_58hybridlevel, clwmr_59hybridlevel, 
     +     clwmr_5hybridlevel, clwmr_60hybridlevel, 
     +     clwmr_61hybridlevel, clwmr_62hybridlevel, 
     +     clwmr_63hybridlevel, clwmr_64hybridlevel, 
     +     clwmr_6hybridlevel, clwmr_7hybridlevel, 
     +     clwmr_8hybridlevel, clwmr_9hybridlevel, 
     +     pres_10hybridlevel, pres_11hybridlevel, 
     +     pres_12hybridlevel, pres_13hybridlevel, 
     +     pres_14hybridlevel, pres_15hybridlevel, 
     +     pres_16hybridlevel, pres_17hybridlevel, 
     +     pres_18hybridlevel, pres_19hybridlevel, pres_1hybridlevel, 
     +     pres_20hybridlevel, pres_21hybridlevel, 
     +     pres_22hybridlevel, pres_23hybridlevel, 
     +     pres_24hybridlevel, pres_25hybridlevel, 
     +     pres_26hybridlevel, pres_27hybridlevel, 
     +     pres_28hybridlevel, pres_29hybridlevel, pres_2hybridlevel, 
     +     pres_30hybridlevel, pres_31hybridlevel, 
     +     pres_32hybridlevel, pres_33hybridlevel, 
     +     pres_34hybridlevel, pres_35hybridlevel, 
     +     pres_36hybridlevel, pres_37hybridlevel, 
     +     pres_38hybridlevel, pres_39hybridlevel, pres_3hybridlevel, 
     +     pres_40hybridlevel, pres_41hybridlevel, 
     +     pres_42hybridlevel, pres_43hybridlevel, 
     +     pres_44hybridlevel, pres_45hybridlevel, 
     +     pres_46hybridlevel, pres_47hybridlevel, 
     +     pres_48hybridlevel, pres_49hybridlevel, pres_4hybridlevel, 
     +     pres_50hybridlevel, pres_51hybridlevel, 
     +     pres_52hybridlevel, pres_53hybridlevel, 
     +     pres_54hybridlevel, pres_55hybridlevel, 
     +     pres_56hybridlevel, pres_57hybridlevel, 
     +     pres_58hybridlevel, pres_59hybridlevel, pres_5hybridlevel, 
     +     pres_60hybridlevel, pres_61hybridlevel, 
     +     pres_62hybridlevel, pres_63hybridlevel, 
     +     pres_64hybridlevel, pres_65hybridlevel, pres_6hybridlevel, 
     +     pres_7hybridlevel, pres_8hybridlevel, pres_9hybridlevel, 
     +     latitude, longitude, time)
c
c the netcdf variables are filled - your zzz write call may go here
c

      pres(:,:,1) = pres_1hybridlevel(:,:,1)
      pres(:,:,2) = pres_2hybridlevel(:,:,1)
      pres(:,:,3) = pres_3hybridlevel(:,:,1)
      pres(:,:,4) = pres_4hybridlevel(:,:,1)
      pres(:,:,5) = pres_5hybridlevel(:,:,1)
      pres(:,:,6) = pres_6hybridlevel(:,:,1)
      pres(:,:,7) = pres_7hybridlevel(:,:,1)
      pres(:,:,8) = pres_8hybridlevel(:,:,1)
      pres(:,:,9) = pres_9hybridlevel(:,:,1)
      pres(:,:,10) = pres_10hybridlevel(:,:,1)
      pres(:,:,11) = pres_11hybridlevel(:,:,1)
      pres(:,:,12) = pres_12hybridlevel(:,:,1)
      pres(:,:,13) = pres_13hybridlevel(:,:,1)
      pres(:,:,14) = pres_14hybridlevel(:,:,1)
      pres(:,:,15) = pres_15hybridlevel(:,:,1)
      pres(:,:,16) = pres_16hybridlevel(:,:,1)
      pres(:,:,17) = pres_17hybridlevel(:,:,1)
      pres(:,:,18) = pres_18hybridlevel(:,:,1)
      pres(:,:,19) = pres_19hybridlevel(:,:,1)
      pres(:,:,20) = pres_20hybridlevel(:,:,1)
      pres(:,:,21) = pres_21hybridlevel(:,:,1)
      pres(:,:,22) = pres_22hybridlevel(:,:,1)
      pres(:,:,23) = pres_23hybridlevel(:,:,1)
      pres(:,:,24) = pres_24hybridlevel(:,:,1)
      pres(:,:,25) = pres_25hybridlevel(:,:,1)
      pres(:,:,26) = pres_26hybridlevel(:,:,1)
      pres(:,:,27) = pres_27hybridlevel(:,:,1)
      pres(:,:,28) = pres_28hybridlevel(:,:,1)
      pres(:,:,29) = pres_29hybridlevel(:,:,1)
      pres(:,:,30) = pres_30hybridlevel(:,:,1)
      pres(:,:,31) = pres_31hybridlevel(:,:,1)
      pres(:,:,32) = pres_32hybridlevel(:,:,1)
      pres(:,:,33) = pres_33hybridlevel(:,:,1)
      pres(:,:,34) = pres_34hybridlevel(:,:,1)
      pres(:,:,35) = pres_35hybridlevel(:,:,1)
      pres(:,:,36) = pres_36hybridlevel(:,:,1)
      pres(:,:,37) = pres_37hybridlevel(:,:,1)
      pres(:,:,38) = pres_38hybridlevel(:,:,1)
      pres(:,:,39) = pres_39hybridlevel(:,:,1)
      pres(:,:,40) = pres_40hybridlevel(:,:,1)
      pres(:,:,41) = pres_41hybridlevel(:,:,1)
      pres(:,:,42) = pres_42hybridlevel(:,:,1)
      pres(:,:,43) = pres_43hybridlevel(:,:,1)
      pres(:,:,44) = pres_44hybridlevel(:,:,1)
      pres(:,:,45) = pres_45hybridlevel(:,:,1)
      pres(:,:,46) = pres_46hybridlevel(:,:,1)
      pres(:,:,47) = pres_47hybridlevel(:,:,1)
      pres(:,:,48) = pres_48hybridlevel(:,:,1)
      pres(:,:,49) = pres_49hybridlevel(:,:,1)
      pres(:,:,50) = pres_50hybridlevel(:,:,1)
      pres(:,:,51) = pres_51hybridlevel(:,:,1)
      pres(:,:,52) = pres_52hybridlevel(:,:,1)
      pres(:,:,53) = pres_53hybridlevel(:,:,1)
      pres(:,:,54) = pres_54hybridlevel(:,:,1)
      pres(:,:,55) = pres_55hybridlevel(:,:,1)
      pres(:,:,56) = pres_56hybridlevel(:,:,1)
      pres(:,:,57) = pres_57hybridlevel(:,:,1)
      pres(:,:,58) = pres_58hybridlevel(:,:,1)
      pres(:,:,59) = pres_59hybridlevel(:,:,1)
      pres(:,:,60) = pres_60hybridlevel(:,:,1)
      pres(:,:,61) = pres_61hybridlevel(:,:,1)
      pres(:,:,62) = pres_62hybridlevel(:,:,1)
      pres(:,:,63) = pres_63hybridlevel(:,:,1)
      pres(:,:,64) = pres_64hybridlevel(:,:,1)
      pres(:,:,65) = pres_65hybridlevel(:,:,1)

      clwc(:,:,1) = clwmr_1hybridlevel(:,:,1)
      clwc(:,:,2) = clwmr_2hybridlevel(:,:,1)
      clwc(:,:,3) = clwmr_3hybridlevel(:,:,1)
      clwc(:,:,4) = clwmr_4hybridlevel(:,:,1)
      clwc(:,:,5) = clwmr_5hybridlevel(:,:,1)
      clwc(:,:,6) = clwmr_6hybridlevel(:,:,1)
      clwc(:,:,7) = clwmr_7hybridlevel(:,:,1)
      clwc(:,:,8) = clwmr_8hybridlevel(:,:,1)
      clwc(:,:,9) = clwmr_9hybridlevel(:,:,1)
      clwc(:,:,10) = clwmr_10hybridlevel(:,:,1)
      clwc(:,:,11) = clwmr_11hybridlevel(:,:,1)
      clwc(:,:,12) = clwmr_12hybridlevel(:,:,1)
      clwc(:,:,13) = clwmr_13hybridlevel(:,:,1)
      clwc(:,:,14) = clwmr_14hybridlevel(:,:,1)
      clwc(:,:,15) = clwmr_15hybridlevel(:,:,1)
      clwc(:,:,16) = clwmr_16hybridlevel(:,:,1)
      clwc(:,:,17) = clwmr_17hybridlevel(:,:,1)
      clwc(:,:,18) = clwmr_18hybridlevel(:,:,1)
      clwc(:,:,19) = clwmr_19hybridlevel(:,:,1)
      clwc(:,:,20) = clwmr_20hybridlevel(:,:,1)
      clwc(:,:,21) = clwmr_21hybridlevel(:,:,1)
      clwc(:,:,22) = clwmr_22hybridlevel(:,:,1)
      clwc(:,:,23) = clwmr_23hybridlevel(:,:,1)
      clwc(:,:,24) = clwmr_24hybridlevel(:,:,1)
      clwc(:,:,25) = clwmr_25hybridlevel(:,:,1)
      clwc(:,:,26) = clwmr_26hybridlevel(:,:,1)
      clwc(:,:,27) = clwmr_27hybridlevel(:,:,1)
      clwc(:,:,28) = clwmr_28hybridlevel(:,:,1)
      clwc(:,:,29) = clwmr_29hybridlevel(:,:,1)
      clwc(:,:,30) = clwmr_30hybridlevel(:,:,1)
      clwc(:,:,31) = clwmr_31hybridlevel(:,:,1)
      clwc(:,:,32) = clwmr_32hybridlevel(:,:,1)
      clwc(:,:,33) = clwmr_33hybridlevel(:,:,1)
      clwc(:,:,34) = clwmr_34hybridlevel(:,:,1)
      clwc(:,:,35) = clwmr_35hybridlevel(:,:,1)
      clwc(:,:,36) = clwmr_36hybridlevel(:,:,1)
      clwc(:,:,37) = clwmr_37hybridlevel(:,:,1)
      clwc(:,:,38) = clwmr_38hybridlevel(:,:,1)
      clwc(:,:,39) = clwmr_39hybridlevel(:,:,1)
      clwc(:,:,40) = clwmr_40hybridlevel(:,:,1)
      clwc(:,:,41) = clwmr_41hybridlevel(:,:,1)
      clwc(:,:,42) = clwmr_42hybridlevel(:,:,1)
      clwc(:,:,43) = clwmr_43hybridlevel(:,:,1)
      clwc(:,:,44) = clwmr_44hybridlevel(:,:,1)
      clwc(:,:,45) = clwmr_45hybridlevel(:,:,1)
      clwc(:,:,46) = clwmr_46hybridlevel(:,:,1)
      clwc(:,:,47) = clwmr_47hybridlevel(:,:,1)
      clwc(:,:,48) = clwmr_48hybridlevel(:,:,1)
      clwc(:,:,49) = clwmr_49hybridlevel(:,:,1)
      clwc(:,:,50) = clwmr_50hybridlevel(:,:,1)
      clwc(:,:,51) = clwmr_51hybridlevel(:,:,1)
      clwc(:,:,52) = clwmr_52hybridlevel(:,:,1)
      clwc(:,:,53) = clwmr_53hybridlevel(:,:,1)
      clwc(:,:,54) = clwmr_54hybridlevel(:,:,1)
      clwc(:,:,55) = clwmr_55hybridlevel(:,:,1)
      clwc(:,:,56) = clwmr_56hybridlevel(:,:,1)
      clwc(:,:,57) = clwmr_57hybridlevel(:,:,1)
      clwc(:,:,58) = clwmr_58hybridlevel(:,:,1)
      clwc(:,:,59) = clwmr_59hybridlevel(:,:,1)
      clwc(:,:,60) = clwmr_60hybridlevel(:,:,1)
      clwc(:,:,61) = clwmr_61hybridlevel(:,:,1)
      clwc(:,:,62) = clwmr_62hybridlevel(:,:,1)
      clwc(:,:,63) = clwmr_63hybridlevel(:,:,1)
      clwc(:,:,64) = clwmr_64hybridlevel(:,:,1)
      clwc(:,:,65) = 0. ! clwmr_65hybridlevel(:,:,1)

!     vertical interpolation
      call vinterp_sub(r_missing_data,nx_l,ny_l,nx_l,ny_l
     .                ,1,nx_l,1,ny_l,21,65
     .                ,pres_p,pres,clwc,clwc_p)

      write(6,*)' range of input clwc_p',minval(clwc_p)
     +                                  ,maxval(clwc_p)

!     fix units to kg/kg as advertised
      clwc_p = clwc_p * .001

      write(6,*)' range of clwc_p (kg/kg)',minval(clwc_p)
     +                                    ,maxval(clwc_p)

!     convert mixing ratio to density (content)
      rho_p(:,:,:) = pres_p(:,:,:) / (r * 273.15)
      clwc_p(:,:,:) = clwc_p(:,:,:) * rho_p(:,:,:)

      write(6,*)' applying 180 degree offset to clwc_p and pres_p'

!     apply 180 degree roll to clwc_p
      do i = 1,nx_l
         iroll = i+nx_l/2
         write(6,*)' i/iroll',i,iroll
         if(iroll .gt. nx_l)iroll = iroll - nx_l   
         buff_p(i,:,:) = clwc_p(iroll,:,:)
      enddo ! i
      clwc_p = buff_p

!     apply 180 degree roll to pres_p
      do i = 1,nx_l
         iroll = i+nx_l/2
         if(iroll .gt. nx_l)iroll = iroll - nx_l   
         buff_p(i,:,:) = pres_p(iroll,:,:)
      enddo ! i
      pres_p = buff_p
      
      write(6,*)' range of rho_p (kg/m^3)',minval(rho_p)
     +                                    ,maxval(rho_p)

      return
      end
c
c  subroutine to read the file 
c
      subroutine read_fim_netcdf(nf_fid, ! latitude, longitude, time, 
     +     clwmr_10hybridlevel, clwmr_11hybridlevel, 
     +     clwmr_12hybridlevel, clwmr_13hybridlevel, 
     +     clwmr_14hybridlevel, clwmr_15hybridlevel, 
     +     clwmr_16hybridlevel, clwmr_17hybridlevel, 
     +     clwmr_18hybridlevel, clwmr_19hybridlevel, 
     +     clwmr_1hybridlevel, clwmr_20hybridlevel, 
     +     clwmr_21hybridlevel, clwmr_22hybridlevel, 
     +     clwmr_23hybridlevel, clwmr_24hybridlevel, 
     +     clwmr_25hybridlevel, clwmr_26hybridlevel, 
     +     clwmr_27hybridlevel, clwmr_28hybridlevel, 
     +     clwmr_29hybridlevel, clwmr_2hybridlevel, 
     +     clwmr_30hybridlevel, clwmr_31hybridlevel, 
     +     clwmr_32hybridlevel, clwmr_33hybridlevel, 
     +     clwmr_34hybridlevel, clwmr_35hybridlevel, 
     +     clwmr_36hybridlevel, clwmr_37hybridlevel, 
     +     clwmr_38hybridlevel, clwmr_39hybridlevel, 
     +     clwmr_3hybridlevel, clwmr_40hybridlevel, 
     +     clwmr_41hybridlevel, clwmr_42hybridlevel, 
     +     clwmr_43hybridlevel, clwmr_44hybridlevel, 
     +     clwmr_45hybridlevel, clwmr_46hybridlevel, 
     +     clwmr_47hybridlevel, clwmr_48hybridlevel, 
     +     clwmr_49hybridlevel, clwmr_4hybridlevel, 
     +     clwmr_50hybridlevel, clwmr_51hybridlevel, 
     +     clwmr_52hybridlevel, clwmr_53hybridlevel, 
     +     clwmr_54hybridlevel, clwmr_55hybridlevel, 
     +     clwmr_56hybridlevel, clwmr_57hybridlevel, 
     +     clwmr_58hybridlevel, clwmr_59hybridlevel, 
     +     clwmr_5hybridlevel, clwmr_60hybridlevel, 
     +     clwmr_61hybridlevel, clwmr_62hybridlevel, 
     +     clwmr_63hybridlevel, clwmr_64hybridlevel, 
     +     clwmr_6hybridlevel, clwmr_7hybridlevel, 
     +     clwmr_8hybridlevel, clwmr_9hybridlevel, 
     +     pres_10hybridlevel, pres_11hybridlevel, 
     +     pres_12hybridlevel, pres_13hybridlevel, 
     +     pres_14hybridlevel, pres_15hybridlevel, 
     +     pres_16hybridlevel, pres_17hybridlevel, 
     +     pres_18hybridlevel, pres_19hybridlevel, pres_1hybridlevel, 
     +     pres_20hybridlevel, pres_21hybridlevel, 
     +     pres_22hybridlevel, pres_23hybridlevel, 
     +     pres_24hybridlevel, pres_25hybridlevel, 
     +     pres_26hybridlevel, pres_27hybridlevel, 
     +     pres_28hybridlevel, pres_29hybridlevel, pres_2hybridlevel, 
     +     pres_30hybridlevel, pres_31hybridlevel, 
     +     pres_32hybridlevel, pres_33hybridlevel, 
     +     pres_34hybridlevel, pres_35hybridlevel, 
     +     pres_36hybridlevel, pres_37hybridlevel, 
     +     pres_38hybridlevel, pres_39hybridlevel, pres_3hybridlevel, 
     +     pres_40hybridlevel, pres_41hybridlevel, 
     +     pres_42hybridlevel, pres_43hybridlevel, 
     +     pres_44hybridlevel, pres_45hybridlevel, 
     +     pres_46hybridlevel, pres_47hybridlevel, 
     +     pres_48hybridlevel, pres_49hybridlevel, pres_4hybridlevel, 
     +     pres_50hybridlevel, pres_51hybridlevel, 
     +     pres_52hybridlevel, pres_53hybridlevel, 
     +     pres_54hybridlevel, pres_55hybridlevel, 
     +     pres_56hybridlevel, pres_57hybridlevel, 
     +     pres_58hybridlevel, pres_59hybridlevel, pres_5hybridlevel, 
     +     pres_60hybridlevel, pres_61hybridlevel, 
     +     pres_62hybridlevel, pres_63hybridlevel, 
     +     pres_64hybridlevel, pres_65hybridlevel, pres_6hybridlevel, 
     +     pres_7hybridlevel, pres_8hybridlevel, pres_9hybridlevel, 
     +     latitude, longitude, time)
c
      include 'netcdf.inc'
      integer latitude, longitude, time,nf_fid, nf_vid, nf_status

      real clwmr_10hybridlevel( longitude,  latitude, time),
     +     clwmr_11hybridlevel( longitude,  latitude, time),
     +     clwmr_12hybridlevel( longitude,  latitude, time),
     +     clwmr_13hybridlevel( longitude,  latitude, time),
     +     clwmr_14hybridlevel( longitude,  latitude, time),
     +     clwmr_15hybridlevel( longitude,  latitude, time),
     +     clwmr_16hybridlevel( longitude,  latitude, time),
     +     clwmr_17hybridlevel( longitude,  latitude, time),
     +     clwmr_18hybridlevel( longitude,  latitude, time),
     +     clwmr_19hybridlevel( longitude,  latitude, time),
     +     clwmr_1hybridlevel( longitude,  latitude, time),
     +     clwmr_20hybridlevel( longitude,  latitude, time),
     +     clwmr_21hybridlevel( longitude,  latitude, time),
     +     clwmr_22hybridlevel( longitude,  latitude, time),
     +     clwmr_23hybridlevel( longitude,  latitude, time),
     +     clwmr_24hybridlevel( longitude,  latitude, time),
     +     clwmr_25hybridlevel( longitude,  latitude, time),
     +     clwmr_26hybridlevel( longitude,  latitude, time),
     +     clwmr_27hybridlevel( longitude,  latitude, time),
     +     clwmr_28hybridlevel( longitude,  latitude, time),
     +     clwmr_29hybridlevel( longitude,  latitude, time),
     +     clwmr_2hybridlevel( longitude,  latitude, time),
     +     clwmr_30hybridlevel( longitude,  latitude, time),
     +     clwmr_31hybridlevel( longitude,  latitude, time),
     +     clwmr_32hybridlevel( longitude,  latitude, time),
     +     clwmr_33hybridlevel( longitude,  latitude, time),
     +     clwmr_34hybridlevel( longitude,  latitude, time),
     +     clwmr_35hybridlevel( longitude,  latitude, time),
     +     clwmr_36hybridlevel( longitude,  latitude, time),
     +     clwmr_37hybridlevel( longitude,  latitude, time),
     +     clwmr_38hybridlevel( longitude,  latitude, time),
     +     clwmr_39hybridlevel( longitude,  latitude, time),
     +     clwmr_3hybridlevel( longitude,  latitude, time),
     +     clwmr_40hybridlevel( longitude,  latitude, time),
     +     clwmr_41hybridlevel( longitude,  latitude, time),
     +     clwmr_42hybridlevel( longitude,  latitude, time),
     +     clwmr_43hybridlevel( longitude,  latitude, time),
     +     clwmr_44hybridlevel( longitude,  latitude, time),
     +     clwmr_45hybridlevel( longitude,  latitude, time),
     +     clwmr_46hybridlevel( longitude,  latitude, time),
     +     clwmr_47hybridlevel( longitude,  latitude, time),
     +     clwmr_48hybridlevel( longitude,  latitude, time),
     +     clwmr_49hybridlevel( longitude,  latitude, time),
     +     clwmr_4hybridlevel( longitude,  latitude, time),
     +     clwmr_50hybridlevel( longitude,  latitude, time),
     +     clwmr_51hybridlevel( longitude,  latitude, time),
     +     clwmr_52hybridlevel( longitude,  latitude, time),
     +     clwmr_53hybridlevel( longitude,  latitude, time),
     +     clwmr_54hybridlevel( longitude,  latitude, time),
     +     clwmr_55hybridlevel( longitude,  latitude, time),
     +     clwmr_56hybridlevel( longitude,  latitude, time),
     +     clwmr_57hybridlevel( longitude,  latitude, time),
     +     clwmr_58hybridlevel( longitude,  latitude, time),
     +     clwmr_59hybridlevel( longitude,  latitude, time),
     +     clwmr_5hybridlevel( longitude,  latitude, time),
     +     clwmr_60hybridlevel( longitude,  latitude, time),
     +     clwmr_61hybridlevel( longitude,  latitude, time),
     +     clwmr_62hybridlevel( longitude,  latitude, time),
     +     clwmr_63hybridlevel( longitude,  latitude, time),
     +     clwmr_64hybridlevel( longitude,  latitude, time),
     +     clwmr_6hybridlevel( longitude,  latitude, time),
     +     clwmr_7hybridlevel( longitude,  latitude, time),
     +     clwmr_8hybridlevel( longitude,  latitude, time),
     +     clwmr_9hybridlevel( longitude,  latitude, time),
     +     pres_10hybridlevel( longitude,  latitude, time),
     +     pres_11hybridlevel( longitude,  latitude, time),
     +     pres_12hybridlevel( longitude,  latitude, time),
     +     pres_13hybridlevel( longitude,  latitude, time),
     +     pres_14hybridlevel( longitude,  latitude, time),
     +     pres_15hybridlevel( longitude,  latitude, time),
     +     pres_16hybridlevel( longitude,  latitude, time),
     +     pres_17hybridlevel( longitude,  latitude, time),
     +     pres_18hybridlevel( longitude,  latitude, time),
     +     pres_19hybridlevel( longitude,  latitude, time),
     +     pres_1hybridlevel( longitude,  latitude, time),
     +     pres_20hybridlevel( longitude,  latitude, time),
     +     pres_21hybridlevel( longitude,  latitude, time),
     +     pres_22hybridlevel( longitude,  latitude, time),
     +     pres_23hybridlevel( longitude,  latitude, time),
     +     pres_24hybridlevel( longitude,  latitude, time),
     +     pres_25hybridlevel( longitude,  latitude, time),
     +     pres_26hybridlevel( longitude,  latitude, time),
     +     pres_27hybridlevel( longitude,  latitude, time),
     +     pres_28hybridlevel( longitude,  latitude, time),
     +     pres_29hybridlevel( longitude,  latitude, time),
     +     pres_2hybridlevel( longitude,  latitude, time),
     +     pres_30hybridlevel( longitude,  latitude, time),
     +     pres_31hybridlevel( longitude,  latitude, time),
     +     pres_32hybridlevel( longitude,  latitude, time),
     +     pres_33hybridlevel( longitude,  latitude, time),
     +     pres_34hybridlevel( longitude,  latitude, time),
     +     pres_35hybridlevel( longitude,  latitude, time),
     +     pres_36hybridlevel( longitude,  latitude, time),
     +     pres_37hybridlevel( longitude,  latitude, time),
     +     pres_38hybridlevel( longitude,  latitude, time),
     +     pres_39hybridlevel( longitude,  latitude, time),
     +     pres_3hybridlevel( longitude,  latitude, time),
     +     pres_40hybridlevel( longitude,  latitude, time),
     +     pres_41hybridlevel( longitude,  latitude, time),
     +     pres_42hybridlevel( longitude,  latitude, time),
     +     pres_43hybridlevel( longitude,  latitude, time),
     +     pres_44hybridlevel( longitude,  latitude, time),
     +     pres_45hybridlevel( longitude,  latitude, time),
     +     pres_46hybridlevel( longitude,  latitude, time),
     +     pres_47hybridlevel( longitude,  latitude, time),
     +     pres_48hybridlevel( longitude,  latitude, time),
     +     pres_49hybridlevel( longitude,  latitude, time),
     +     pres_4hybridlevel( longitude,  latitude, time),
     +     pres_50hybridlevel( longitude,  latitude, time),
     +     pres_51hybridlevel( longitude,  latitude, time),
     +     pres_52hybridlevel( longitude,  latitude, time),
     +     pres_53hybridlevel( longitude,  latitude, time),
     +     pres_54hybridlevel( longitude,  latitude, time),
     +     pres_55hybridlevel( longitude,  latitude, time),
     +     pres_56hybridlevel( longitude,  latitude, time),
     +     pres_57hybridlevel( longitude,  latitude, time),
     +     pres_58hybridlevel( longitude,  latitude, time),
     +     pres_59hybridlevel( longitude,  latitude, time),
     +     pres_5hybridlevel( longitude,  latitude, time),
     +     pres_60hybridlevel( longitude,  latitude, time),
     +     pres_61hybridlevel( longitude,  latitude, time),
     +     pres_62hybridlevel( longitude,  latitude, time),
     +     pres_63hybridlevel( longitude,  latitude, time),
     +     pres_64hybridlevel( longitude,  latitude, time),
     +     pres_65hybridlevel( longitude,  latitude, time),
     +     pres_6hybridlevel( longitude,  latitude, time),
     +     pres_7hybridlevel( longitude,  latitude, time),
     +     pres_8hybridlevel( longitude,  latitude, time),
     +     pres_9hybridlevel( longitude,  latitude, time)
!     double precision latitude(latitude), longitude(longitude),
!    +     time(time)

      write(6,*)' read_fim_netcdf ',latitude,longitude,time

c   variables of type real
c
c     variable        netcdf long name
c     clwmr_10hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_10hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_10hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_10hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_10hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_11hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_11hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_11hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_11hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_11hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_12hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_12hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_12hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_12hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_12hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_13hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_13hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_13hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_13hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_13hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_14hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_14hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_14hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_14hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_14hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_15hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_15hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_15hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_15hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_15hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_16hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_16hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_16hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_16hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_16hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_17hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_17hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_17hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_17hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_17hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_18hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_18hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_18hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_18hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_18hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_19hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_19hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_19hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_19hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_19hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_1hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_1hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_1hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_1hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_1hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_20hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_20hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_20hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_20hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_20hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_21hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_21hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_21hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_21hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_21hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_22hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_22hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_22hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_22hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_22hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_23hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_23hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_23hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_23hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_23hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_24hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_24hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_24hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_24hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_24hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_25hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_25hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_25hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_25hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_25hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_26hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_26hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_26hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_26hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_26hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_27hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_27hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_27hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_27hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_27hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_28hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_28hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_28hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_28hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_28hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_29hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_29hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_29hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_29hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_29hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_2hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_2hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_2hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_2hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_2hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_30hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_30hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_30hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_30hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_30hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_31hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_31hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_31hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_31hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_31hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_32hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_32hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_32hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_32hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_32hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_33hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_33hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_33hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_33hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_33hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_34hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_34hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_34hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_34hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_34hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_35hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_35hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_35hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_35hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_35hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_36hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_36hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_36hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_36hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_36hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_37hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_37hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_37hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_37hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_37hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_38hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_38hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_38hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_38hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_38hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_39hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_39hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_39hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_39hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_39hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_3hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_3hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_3hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_3hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_3hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_40hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_40hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_40hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_40hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_40hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_41hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_41hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_41hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_41hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_41hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_42hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_42hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_42hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_42hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_42hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_43hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_43hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_43hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_43hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_43hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_44hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_44hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_44hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_44hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_44hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_45hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_45hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_45hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_45hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_45hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_46hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_46hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_46hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_46hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_46hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_47hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_47hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_47hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_47hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_47hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_48hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_48hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_48hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_48hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_48hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_49hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_49hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_49hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_49hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_49hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_4hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_4hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_4hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_4hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_4hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_50hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_50hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_50hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_50hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_50hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_51hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_51hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_51hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_51hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_51hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_52hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_52hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_52hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_52hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_52hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_53hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_53hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_53hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_53hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_53hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_54hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_54hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_54hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_54hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_54hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_55hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_55hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_55hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_55hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_55hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_56hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_56hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_56hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_56hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_56hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_57hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_57hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_57hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_57hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_57hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_58hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_58hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_58hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_58hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_58hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_59hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_59hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_59hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_59hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_59hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_5hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_5hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_5hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_5hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_5hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_60hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_60hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_60hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_60hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_60hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_61hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_61hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_61hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_61hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_61hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_62hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_62hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_62hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_62hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_62hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_63hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_63hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_63hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_63hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_63hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_64hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_64hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_64hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_64hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_64hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_6hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_6hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_6hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_6hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_6hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_7hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_7hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_7hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_7hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_7hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_8hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_8hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_8hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_8hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_8hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     clwmr_9hybridlevel"cloud mixing ratio"
c
      nf_status=nf_inq_varid(nf_fid,'clwmr_9hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clwmr_9hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clwmr_9hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clwmr_9hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_10hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_10hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_10hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_10hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_10hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_11hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_11hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_11hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_11hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_11hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_12hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_12hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_12hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_12hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_12hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_13hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_13hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_13hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_13hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_13hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_14hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_14hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_14hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_14hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_14hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_15hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_15hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_15hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_15hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_15hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_16hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_16hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_16hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_16hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_16hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_17hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_17hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_17hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_17hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_17hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_18hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_18hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_18hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_18hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_18hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_19hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_19hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_19hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_19hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_19hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_1hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_1hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_1hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_1hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_1hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_20hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_20hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_20hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_20hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_20hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_21hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_21hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_21hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_21hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_21hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_22hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_22hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_22hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_22hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_22hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_23hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_23hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_23hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_23hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_23hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_24hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_24hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_24hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_24hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_24hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_25hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_25hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_25hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_25hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_25hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_26hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_26hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_26hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_26hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_26hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_27hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_27hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_27hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_27hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_27hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_28hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_28hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_28hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_28hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_28hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_29hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_29hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_29hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_29hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_29hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_2hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_2hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_2hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_2hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_2hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_30hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_30hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_30hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_30hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_30hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_31hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_31hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_31hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_31hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_31hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_32hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_32hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_32hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_32hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_32hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_33hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_33hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_33hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_33hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_33hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_34hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_34hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_34hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_34hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_34hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_35hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_35hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_35hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_35hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_35hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_36hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_36hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_36hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_36hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_36hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_37hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_37hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_37hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_37hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_37hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_38hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_38hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_38hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_38hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_38hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_39hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_39hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_39hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_39hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_39hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_3hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_3hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_3hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_3hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_3hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_40hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_40hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_40hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_40hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_40hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_41hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_41hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_41hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_41hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_41hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_42hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_42hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_42hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_42hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_42hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_43hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_43hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_43hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_43hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_43hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_44hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_44hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_44hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_44hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_44hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_45hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_45hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_45hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_45hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_45hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_46hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_46hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_46hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_46hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_46hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_47hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_47hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_47hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_47hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_47hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_48hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_48hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_48hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_48hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_48hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_49hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_49hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_49hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_49hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_49hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_4hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_4hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_4hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_4hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_4hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_50hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_50hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_50hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_50hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_50hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_51hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_51hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_51hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_51hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_51hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_52hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_52hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_52hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_52hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_52hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_53hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_53hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_53hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_53hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_53hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_54hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_54hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_54hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_54hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_54hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_55hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_55hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_55hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_55hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_55hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_56hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_56hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_56hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_56hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_56hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_57hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_57hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_57hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_57hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_57hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_58hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_58hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_58hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_58hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_58hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_59hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_59hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_59hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_59hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_59hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_5hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_5hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_5hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_5hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_5hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_60hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_60hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_60hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_60hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_60hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_61hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_61hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_61hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_61hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_61hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_62hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_62hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_62hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_62hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_62hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_63hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_63hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_63hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_63hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_63hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_64hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_64hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_64hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_64hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_64hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_65hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_65hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_65hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_65hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_65hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_6hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_6hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_6hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_6hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_6hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_7hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_7hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_7hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_7hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_7hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_8hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_8hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_8hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_8hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_8hybridlevel'
       endif
      endif
c
c     variable        netcdf long name
c     pres_9hybridlevel"pressure"
c
      nf_status=nf_inq_varid(nf_fid,'pres_9hybridlevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pres_9hybridlevel'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pres_9hybridlevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pres_9hybridlevel'
       endif
      endif

c   variables of type int
c

c   variables of type double
c
c
c     variable        netcdf long name
c     latitude      "latitude"
c
      if(.false.)then

      nf_status=nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitude'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,latitude)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitude'
       endif
      endif
c
c     variable        netcdf long name
c     longitude     "longitude"
c
      nf_status=nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for longitude'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,longitude)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitude'
       endif
      endif
c
c     variable        netcdf long name
c     time          "verification time generated by wgrib2 function verftime()"
c
      nf_status=nf_inq_varid(nf_fid,'time',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for time'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,time)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for time'
       endif
      endif
 
      endif

c   variables of type char
c

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
