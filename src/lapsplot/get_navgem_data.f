      subroutine get_navgem_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l,nz_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,pres_p
     +                   ,clwc_p
     +                   ,cice_p
     +                   ,rain_p
     +                   ,snow_p
     +                   ,seaice,snow_depth
     +                   ,lun_out
     +                   ,istatus)

      use mem_allsky, only: aod_3d

      include 'netcdf.inc'

      character*(*) filename

      integer latitude, level, longitude,nf_fid, nf_vid, nf_status

      parameter (lvlp=21)

      real ht_p(nx_l,ny_l,lvlp)
      real pres_p(nx_l,ny_l,lvlp)
      real clwc_p(nx_l,ny_l,lvlp)
      real cice_p(nx_l,ny_l,lvlp)
      real rain_p(nx_l,ny_l,lvlp)
      real snow_p(nx_l,ny_l,lvlp)
      real coarse_ext_p(nx_l,ny_l,lvlp)
      real fine_ext_p(nx_l,ny_l,lvlp)
      real total_ext_p(nx_l,ny_l,lvlp)
      real seaice(nx_l,ny_l)
      real snow_depth(nx_l,ny_l)

      write(6,*)' entered get_navgem_data'

      if(nz_l .ne. lvlp)then
        write(6,*)' incorrect vertical levels',nz_l,lvlp
        istatus=0
        return
      endif

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
c get size of level
c
      nf_status=nf_inq_dimid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim level'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,level)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim level'
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
      call read_navgem_data(nf_fid, latitude, level, longitude,
     +     i4time_sys, ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, clwc_p, cice_p, pres_p, 
     +     coarse_ext_p, fine_ext_p, total_ext_p, lvlp,
     +     seaice, snow_depth, istatus)

      if(.true.)then
         write(6,*)' use navgem aerosols'
         aod_3d(:,:,:) = total_ext_p(:,:,:)
      endif

      return
      end
c
c
      subroutine read_navgem_data(nf_fid, latitude, level, longitude,
     +     i4time_sys, ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, clw_p, ciw_p, pres_p,
     +     coarse_ext_p, fine_ext_p, total_ext_p, lvlp,
     +     seaice, snow_depth, istatus)


      include 'netcdf.inc'
      integer latitude, level, longitude,nf_fid, nf_vid, nf_status


      real abf_aot( longitude, latitude), abf_ext( longitude, 
     +     latitude, level), abf_mass_concentration( longitude, 
     +     latitude, level), masl( longitude,  latitude, level),
     +     ss_aot( longitude, latitude), ss_ext( longitude, 
     +     latitude, level), ss_mass_concentration( longitude, 
     +     latitude, level), cf( longitude,  latitude, level), ciw(
     +     longitude,  latitude, level), clw( longitude,  latitude,
     +     level), coarse_aot( longitude, latitude), coarse_ext(
     +     longitude,  latitude, level),
     +     coarse_mode_mass_concentration( longitude,  latitude,
     +     level), dust_aot( longitude, latitude), dust_ext(
     +     longitude,  latitude, level), dust_mass_concentration(
     +     longitude,  latitude, level), dz( longitude,  latitude,
     +     level), fine_aot( longitude, latitude), fine_ext(
     +     longitude,  latitude, level),
     +     fine_mode_mass_concentration( longitude,  latitude,
     +     level), height_agl( longitude,  latitude, level),
     +     latitude_a(latitude), level_a(level), longitude_a(longitude),
     +     pressure( longitude,  latitude, level), ps( longitude,
     +     latitude), q( longitude,  latitude, level), rh( longitude,
     +     latitude, level), sigma_a(level), sigma_b(level),
     +     smoke_aot( longitude, latitude), smoke_ext( longitude, 
     +     latitude, level), smoke_mass_concentration( longitude, 
     +     latitude, level), temperature( longitude,  latitude,
     +     level), total_aot( longitude, latitude), total_ext(
     +     longitude,  latitude, level), tprec( longitude, latitude),
     +     u( longitude,  latitude, level), v( longitude,  latitude,
     +     level), w( longitude,  latitude, level), wind_speed(
     +     longitude,  latitude, level),
     +     seaice(longitude, latitude), snow_depth(longitude, latitude)
c
      real pres_p(longitude,latitude,lvlp)
      real ciw_p(longitude,latitude,lvlp)
      real clw_p(longitude,latitude,lvlp)
      real coarse_ext_p(longitude,latitude,lvlp)
      real fine_ext_p(longitude,latitude,lvlp)
      real total_ext_p(longitude,latitude,lvlp)

!     declarations for 'write_nvg' call
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

      call read_navgem_netcdf(nf_fid, latitude, level, longitude, 
     +     abf_aot, abf_ext, abf_mass_concentration, masl, ss_aot, 
     +     ss_ext, ss_mass_concentration, cf, ciw, clw, coarse_aot, 
     +     coarse_ext, coarse_mode_mass_concentration, dust_aot, 
     +     dust_ext, dust_mass_concentration, dz, fine_aot, fine_ext, 
     +     fine_mode_mass_concentration, height_agl, 
     +     latitude_a, level_a, longitude_a,
     +     pressure, ps, q, rh, sigma_a, sigma_b, 
     +     smoke_aot, smoke_ext, smoke_mass_concentration, 
     +     temperature, total_aot, total_ext, tprec, u, v, w, 
     +     wind_speed,seaice,snow_depth)
c
c the netcdf variables are filled - your nvg write call may go here
c

      write(6,*)' range of input ciw',minval(ciw)
     +                               ,maxval(ciw)

      write(6,*)' range of input clw',minval(clw)
     +                               ,maxval(clw)

      write(6,*)' range of input pressure',minval(pressure)
     +                                    ,maxval(pressure)

      write(6,*)' range of simgrid pres_p',minval(pres_p)
     +                                    ,maxval(pres_p)

!     vertical interpolation
      call vinterp_sub(r_missing_data,nx_l,ny_l,nx_l,ny_l
     .                ,1,nx_l,1,ny_l,lvlp,level
     .                ,pres_p,pressure,ciw,ciw_p)

!     vertical interpolation
      call vinterp_sub(r_missing_data,nx_l,ny_l,nx_l,ny_l
     .                ,1,nx_l,1,ny_l,lvlp,level
     .                ,pres_p,pressure,clw,clw_p)

!     vertical interpolation
      call vinterp_sub(r_missing_data,nx_l,ny_l,nx_l,ny_l
     .                ,1,nx_l,1,ny_l,lvlp,level
     .                ,pres_p,pressure,coarse_ext,coarse_ext_p)

!     vertical interpolation
      call vinterp_sub(r_missing_data,nx_l,ny_l,nx_l,ny_l
     .                ,1,nx_l,1,ny_l,lvlp,level
     .                ,pres_p,pressure,fine_ext,fine_ext_p)

!     vertical interpolation
      call vinterp_sub(r_missing_data,nx_l,ny_l,nx_l,ny_l
     .                ,1,nx_l,1,ny_l,lvlp,level
     .                ,pres_p,pressure,total_ext,total_ext_p)

      write(6,*)' range of interp ciw_p',minval(ciw_p)
     +                                  ,maxval(ciw_p)

      write(6,*)' range of interp clw_p',minval(clw_p)
     +                                  ,maxval(clw_p)

      write(6,*)' range of interp coarse_ext_p',minval(coarse_ext_p)
     +                                         ,maxval(coarse_ext_p)

      write(6,*)' range of interp fine_ext_p',minval(fine_ext_p)
     +                                       ,maxval(fine_ext_p)

      write(6,*)' range of interp total_ext_p',minval(total_ext_p)
     +                                        ,maxval(total_ext_p)

      return
      end
c
c  subroutine to read the file 
c
      subroutine read_navgem_netcdf(nf_fid, latitude, level, 
     +     longitude, abf_aot, abf_ext, abf_mass_concentration, masl, 
     +     ss_aot, ss_ext, ss_mass_concentration, cf, ciw, clw, 
     +     coarse_aot, coarse_ext, coarse_mode_mass_concentration, 
     +     dust_aot, dust_ext, dust_mass_concentration, dz, fine_aot, 
     +     fine_ext, fine_mode_mass_concentration, height_agl, 
     +     latitude_a, level_a, longitude_a,
     +     pressure, ps, q, rh, sigma_a, 
     +     sigma_b, smoke_aot, smoke_ext, smoke_mass_concentration, 
     +     temperature, total_aot, total_ext, tprec, u, v, w, 
     +     wind_speed,seaice,snow_depth)
c
      include 'netcdf.inc'
      integer latitude, level, longitude,nf_fid, nf_vid, nf_status

      real abf_aot( longitude, latitude), abf_ext( longitude, 
     +     latitude, level), abf_mass_concentration( longitude, 
     +     latitude, level), masl( longitude,  latitude, level),
     +     ss_aot( longitude, latitude), ss_ext( longitude, 
     +     latitude, level), ss_mass_concentration( longitude, 
     +     latitude, level), cf( longitude,  latitude, level), ciw(
     +     longitude,  latitude, level), clw( longitude,  latitude,
     +     level), coarse_aot( longitude, latitude), coarse_ext(
     +     longitude,  latitude, level),
     +     coarse_mode_mass_concentration( longitude,  latitude,
     +     level), dust_aot( longitude, latitude), dust_ext(
     +     longitude,  latitude, level), dust_mass_concentration(
     +     longitude,  latitude, level), dz( longitude,  latitude,
     +     level), fine_aot( longitude, latitude), fine_ext(
     +     longitude,  latitude, level),
     +     fine_mode_mass_concentration( longitude,  latitude,
     +     level), height_agl( longitude,  latitude, level),
     +     latitude_a(latitude), level_a(level), longitude_a(longitude),
     +     pressure( longitude,  latitude, level), ps( longitude,
     +     latitude), q( longitude,  latitude, level), rh( longitude,
     +      latitude, level), sigma_a(level), sigma_b(level),
     +     smoke_aot( longitude, latitude), smoke_ext( longitude, 
     +     latitude, level), smoke_mass_concentration( longitude, 
     +     latitude, level), temperature( longitude,  latitude,
     +     level), total_aot( longitude, latitude), total_ext(
     +     longitude,  latitude, level), tprec( longitude, latitude),
     +     u( longitude,  latitude, level), v( longitude,  latitude,
     +     level), w( longitude,  latitude, level), wind_speed(
     +     longitude,  latitude, level),
     +     seaice(longitude, latitude), snow_depth(longitude, latitude)

c
c   variables of type real
c
c     variable        netcdf long name
c     abf_aot       "anthropogenic and biogenic fine aerosol optical thickness at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'abf_aot',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for abf_aot'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,abf_aot)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for abf_aot'
       endif
      endif
c
c     variable        netcdf long name
c     abf_ext       "anthropogenic and biogenic fine aerosol extinction at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'abf_ext',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for abf_ext'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,abf_ext)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for abf_ext'
       endif
      endif
c
c     variable        netcdf long name
c     abf_mass_concentration"anthropogenic and biogenic fine_mass_concentration"
c
      nf_status=nf_inq_varid(nf_fid,'abf_mass_concentration',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for abf_mass_concentration'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,abf_mass_concentration)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for abf_mass_concentration'
       endif
      endif
c
c     variable        netcdf long name
c     masl          "meters above sea level"
c
      nf_status=nf_inq_varid(nf_fid,'masl',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for masl'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,masl)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for masl'
       endif
      endif
c
c     variable        netcdf long name
c     ss_aot        "sea salt aot"
c
      nf_status=nf_inq_varid(nf_fid,'ss_aot',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for ss_aot'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,ss_aot)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for ss_aot'
       endif
      endif
c
c     variable        netcdf long name
c     ss_ext        "sea salt  aerosol extinction at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'ss_ext',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for ss_ext'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,ss_ext)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for ss_ext'
       endif
      endif
c
c     variable        netcdf long name
c     ss_mass_concentration"sea_salt_mass_concentration"
c
      nf_status=nf_inq_varid(nf_fid,'ss_mass_concentration',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for ss_mass_concentration'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,ss_mass_concentration)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for ss_mass_concentration'
       endif
      endif
c
c     variable        netcdf long name
c     cf            "cloud fraction"
c
      nf_status=nf_inq_varid(nf_fid,'cf',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for cf'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,cf)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for cf'
       endif
      endif
c
c     variable        netcdf long name
c     ciw           "cloud ice water"
c
      nf_status=nf_inq_varid(nf_fid,'ciw',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for ciw'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,ciw)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for ciw'
       endif
      endif
c
c     variable        netcdf long name
c     clw           "cloud liquid water"
c
      nf_status=nf_inq_varid(nf_fid,'clw',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for clw'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,clw)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for clw'
       endif
      endif
c
c     variable        netcdf long name
c     coarse_aot    "coarse mode  aerosol optical thickness at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'coarse_aot',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for coarse_aot'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,coarse_aot)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for coarse_aot'
       endif
      endif
c
c     variable        netcdf long name
c     coarse_ext    "coarse mode  aerosol extinction at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'coarse_ext',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for coarse_ext'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,coarse_ext)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for coarse_ext'
       endif
      endif
c
c     variable        netcdf long name
c     dust_aot      "dust aerosol optical thickness at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'dust_aot',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dust_aot'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,dust_aot)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dust_aot'
       endif
      endif
c
c     variable        netcdf long name
c     dust_ext      "dust  aerosol extinction at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'dust_ext',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dust_ext'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,dust_ext)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dust_ext'
       endif
      endif
c
c     variable        netcdf long name
c     dust_mass_concentration"dust_mass_concentration"
c
      nf_status=nf_inq_varid(nf_fid,'dust_mass_concentration',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dust_mass_concentration'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,dust_mass_concentration)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dust_mass_concentration'
       endif
      endif
c
c     variable        netcdf long name
c     dz            "layer thickness"
c
      nf_status=nf_inq_varid(nf_fid,'dz',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dz'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,dz)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dz'
       endif
      endif
c
c     variable        netcdf long name
c     fine_aot      "fine mode  aerosol optical thickness at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'fine_aot',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for fine_aot'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,fine_aot)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for fine_aot'
       endif
      endif
c
c     variable        netcdf long name
c     fine_ext      "fine mode  aerosol extinction at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'fine_ext',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for fine_ext'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,fine_ext)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for fine_ext'
       endif
      endif
c
c     variable        netcdf long name
c     height_agl    "height above ground level"
c
      nf_status=nf_inq_varid(nf_fid,'height_agl',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for height_agl'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,height_agl)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for height_agl'
       endif
      endif
c
c     variable        netcdf long name
c     latitude      "latitude"
c
      nf_status=nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for latitude'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,latitude_a)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for latitude'
       endif
      endif
c
c     variable        netcdf long name
c     level         "model levels"
c
      nf_status=nf_inq_varid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for level'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,level_a)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for level'
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
       nf_status=nf_get_var_real(nf_fid,nf_vid,longitude_a)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for longitude'
       endif
      endif
c
c     variable        netcdf long name
c     pressure      "pressue"
c
      nf_status=nf_inq_varid(nf_fid,'pressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pressure'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pressure)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pressure'
       endif
      endif
c
c     variable        netcdf long name
c     ps            "surface pressure"
c
      nf_status=nf_inq_varid(nf_fid,'ps',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for ps'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,ps)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for ps'
       endif
      endif
c
c     variable        netcdf long name
c     q             "specific humidity"
c
      nf_status=nf_inq_varid(nf_fid,'q',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for q'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,q)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for q'
       endif
      endif
c
c     variable        netcdf long name
c     rh            "relative humidity"
c
      nf_status=nf_inq_varid(nf_fid,'rh',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for rh'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,rh)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for rh'
       endif
      endif
c
c     variable        netcdf long name
c     sigma_a       "sigma_a"
c
      nf_status=nf_inq_varid(nf_fid,'sigma_a',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for sigma_a'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,sigma_a)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for sigma_a'
       endif
      endif
c
c     variable        netcdf long name
c     sigma_b       "sigma_b"
c
      nf_status=nf_inq_varid(nf_fid,'sigma_b',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for sigma_b'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,sigma_b)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for sigma_b'
       endif
      endif
c
c     variable        netcdf long name
c     smoke_aot     "smoke aot"
c
      nf_status=nf_inq_varid(nf_fid,'smoke_aot',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for smoke_aot'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,smoke_aot)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for smoke_aot'
       endif
      endif
c
c     variable        netcdf long name
c     smoke_ext     "smoke  aerosol extinction at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'smoke_ext',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for smoke_ext'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,smoke_ext)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for smoke_ext'
       endif
      endif
c
c     variable        netcdf long name
c     smoke_mass_concentration"smoke_mass_concentration"
c
      nf_status=nf_inq_varid(nf_fid,'smoke_mass_concentration',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for smoke_mass_concentration'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,smoke_mass_concentration)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for smoke_mass_concentration'
       endif
      endif
c
c     variable        netcdf long name
c     temperature   "temperature"
c
      nf_status=nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for temperature'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,temperature)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for temperature'
       endif
      endif
c
c     variable        netcdf long name
c     total_aot     "total aot"
c
      nf_status=nf_inq_varid(nf_fid,'total_aot',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for total_aot'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,total_aot)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for total_aot'
       endif
      endif
c
c     variable        netcdf long name
c     total_ext     "total  aerosol extinction at 550nm"
c
      nf_status=nf_inq_varid(nf_fid,'total_ext',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for total_ext'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,total_ext)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for total_ext'
       endif
      endif
c
c     variable        netcdf long name
c     tprec         "total precipitation"
c
      nf_status=nf_inq_varid(nf_fid,'tprec',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for tprec'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,tprec)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for tprec'
       endif
      endif
c
c     variable        netcdf long name
c     u             "eastward wind"
c
      nf_status=nf_inq_varid(nf_fid,'u',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for u'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,u)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for u'
       endif
      endif
c
c     variable        netcdf long name
c     v             "northward wind"
c
      nf_status=nf_inq_varid(nf_fid,'v',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for v'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,v)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for v'
       endif
      endif
c
c     variable        netcdf long name
c     w             "omega"
c
      nf_status=nf_inq_varid(nf_fid,'w',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for w'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,w)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for w'
       endif
      endif
c
c     variable        netcdf long name
c     wind_speed    "wind speed (cup)"
c
      nf_status=nf_inq_varid(nf_fid,'wind_speed',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for wind_speed'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,wind_speed)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for wind_speed'
       endif
      endif

c
c     variable        netcdf long name
c     seaice
c
      nf_status=nf_inq_varid(nf_fid,'seaice',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for seaice'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,seaice)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for seaice'
       endif
      endif

c
c     variable        netcdf long name
c     snow_depth    
c
      nf_status=nf_inq_varid(nf_fid,'snow_depth',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for snow_depth'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,snow_depth)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for snow_depth'
       endif
      endif


c   variables of type int
c

c   variables of type double
c


c   variables of type char
c

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
