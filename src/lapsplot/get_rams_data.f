      subroutine get_rams_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l,nz_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lat,lon
     +                   ,pres_3d
     +                   ,ht_p
     +                   ,clwc_p
     +                   ,cice_p
     +                   ,rain_p
     +                   ,snow_p
     +                   ,aext_p
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer x, y, z,nf_fid, nf_vid, nf_status
      real pres_3d(nx_l,ny_l,nz_l)
      real ht_p(nx_l,ny_l,nz_l)
      real clwc_p(nx_l,ny_l,nz_l)
      real cice_p(nx_l,ny_l,nz_l)
      real rain_p(nx_l,ny_l,nz_l)
      real snow_p(nx_l,ny_l,nz_l)
      real aext_p(nx_l,ny_l,nz_l)
      real glat(nx_l,ny_l)
      real glon(nx_l,ny_l)
      real lat(nx_l,ny_l)
      real lon(nx_l,ny_l)
      real gri(nx_l,ny_l)
      real grj(nx_l,ny_l)
      real buf_3d(nx_l,ny_l,nz_l)
      real buf_2d(nx_l,ny_l)
      logical l_hinterp /.false./
      logical wrapped /.false./
c
c  open netcdf file for reading
c
      write(6,*)' opening ',trim(filename)

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
c get size of x
c
      nf_status=nf_inq_dimid(nf_fid,'x',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,x)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
      endif
c
c get size of y
c
      nf_status=nf_inq_dimid(nf_fid,'y',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim y'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,y)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim y'
      endif
c
c get size of z
c
      nf_status=nf_inq_dimid(nf_fid,'z',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,z)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z'
      endif

      if(x .eq. nx_l .and. y .eq. ny_l .and. z .eq. nz_l)then
        write(6,*)' model grid dimensions match laps grid'
      else
        write(6,*)' model grid deviates from laps grid'
        write(6,*)x,y,z 
        write(6,*)nx_l,ny_l,nz_l
        stop
      endif

      call read_rams_data(nf_fid, x, y, z, i4time_sys,
     +     ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, glat, glon, ht_p, clwc_p,
     +     cice_p, rain_p, snow_p, aext_p, istatus)

      if(istatus .eq. 1)then
        write(6,*)' model corners'
        write(6,*)' ll ',glat(1,1),glon(1,1)
        write(6,*)' ul ',glat(1,ny_l),glon(1,ny_l)
        write(6,*)' lr ',glat(nx_l,1),glon(nx_l,1)
        write(6,*)' ur ',glat(nx_l,ny_l),glon(nx_l,ny_l)

        write(6,*)' laps corners'
        write(6,*)' ll ',lat(1,1),lon(1,1)
        write(6,*)' ul ',lat(1,ny_l),lon(1,ny_l)
        write(6,*)' lr ',lat(nx_l,1),lon(nx_l,1)
        write(6,*)' ur ',lat(nx_l,ny_l),lon(nx_l,ny_l)

        ricen = (float(x)+1.)/2.
        rjcen = (float(y)+1.)/2.
        call bilinear_laps(ricen,rjcen,nx_l,ny_l,glat,glat_cen)
        call bilinear_laps(ricen,rjcen,nx_l,ny_l,glon,glon_cen)
        call bilinear_laps(ricen,rjcen,nx_l,ny_l,lat,lat_cen)
        call bilinear_laps(ricen,rjcen,nx_l,ny_l,lon,lon_cen)
        write(6,*)' rams center is ',glat_cen,glon_cen
        write(6,*)' laps center is ',lat_cen,lon_cen
      else
        write(6,*)' istatus from read_rams_data is ',istatus
        return
      endif

      if(l_hinterp)then
        call latlon_to_grij_new(lat,lon,nx_l,ny_l,
     1                          glat,glon,gri,grj,x,y,
     1                          istatus_grij)

        write(6,*)' rams gri,grj at laps corners'
        write(6,*)' ll ',gri(1,1),grj(1,1)
        write(6,*)' ul ',gri(1,ny_l),grj(1,ny_l)
        write(6,*)' lr ',gri(nx_l,1),grj(nx_l,1)
        write(6,*)' ur ',gri(nx_l,ny_l),grj(nx_l,ny_l)

        call stats_missing(pres_3d,x,y,nz_l,r_missing_data,nmiss)
        write(6,*)' nmiss for pres_3d is ',nmiss
        buf_3d = pres_3d
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,pres_3d
     1                       ,wrapped)
        call stats_missing(pres_3d,x,y,nz_l,r_missing_data,nmiss)
        write(6,*)' nmiss for pres_3d is ',nmiss

        stop

        buf_3d = ht_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,ht_p
     1                       ,wrapped)

        buf_3d = clwc_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,clwc_p
     1                       ,wrapped)

        buf_3d = cice_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,cice_p
     1                       ,wrapped)

        buf_3d = rain_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,rain_p
     1                       ,wrapped)

        buf_3d = rain_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,rain_p
     1                       ,wrapped)

        buf_3d = snow_p
        call hinterp_field_3d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_3d,snow_p
     1                       ,wrapped)

        buf_2d = glat
        call hinterp_field_2d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_2d,glat
     1                       ,wrapped)

        buf_2d = glon
        call hinterp_field_2d(x,y,nx_l,ny_l,nz_l,gri,grj,buf_2d,glon
     1                       ,wrapped)

        write(6,*)' remapped glat/glon arrays at laps corners'
        write(6,*)' ll ',glat(1,1),glon(1,1)
        write(6,*)' ul ',glat(1,ny_l),glon(1,ny_l)
        write(6,*)' lr ',glat(nx_l,1),glon(nx_l,1)
        write(6,*)' ur ',glat(nx_l,ny_l),glon(nx_l,ny_l)

      endif

      return
      end
c
c
      subroutine read_rams_data(nf_fid, x, y, z, i4time_sys,
     +     ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, glat, glon, ht, lwc,
     +     ice, rai, sno, aext, istatus)


      include 'netcdf.inc'
      integer x, y, z,nf_fid, nf_vid, nf_status
      integer ht_fcinv(z), ice_fcinv(z), lwc_fcinv(z), rai_fcinv(z),
     +     sno_fcinv(z)
      real glat( x, y), glon( x, y), ht( x,  y, z), ice( x,  y, z),
     +     level(z), lwc( x,  y, z), rai( x,  y, z), rtgt_ascii( x,
     +     y), sno( x,  y, z), topt( x, y), topt_ascii( x, y)
      real aext(x,y,z)

!     declarations for 'write_zzz' call
      logical l_closest_time, l_closest_time_i, l_in_domain
      real*4 lat_a(nx_l,ny_l)
      real*4 lon_a(nx_l,ny_l)
      real*4 topo_a(nx_l,ny_l)

      write(6,*)' subroutine read_rams_data'

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting r_missing_data'
          return
      endif

!     call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat_a,lon_a,
!    1            topo_a,1.0,rnorth,south,east,west,istatus)
!     if(istatus .ne. 1)then
!         write(6,*)' error in get_domain_perimeter'
!         return
!     endif

      write(6,*)' calling read_rams_netcdf'

      call read_rams_netcdf(nf_fid, x, y, z, glat, glon, ht, ice, 
     +     level, lwc, rai, aext, rtgt_ascii, sno, topt, topt_ascii, 
     +     ht_fcinv, ice_fcinv, lwc_fcinv, rai_fcinv, sno_fcinv)
c
c the netcdf variables are filled - your zzz write call may go here
c
!     initial loop through obs to get times and stanums
      where (lwc .eq. r_missing_data)lwc = 0.
      where (ice .eq. r_missing_data)ice = 0.
      where (rai .eq. r_missing_data)rai = 0.
      where (sno .eq. r_missing_data)sno = 0.

!     fill in missing heights
      do iz = 1,z-1
          where (ht(:,:,iz) .eq. r_missing_data)ht(:,:,iz) = iz * 100.
      enddo 

      where (ht(:,:,z) .eq. r_missing_data)ht(:,:,z) = 20000.

      return
      end
c
c  subroutine to read the file "rams microphysics interpolated to constant pressure surfaces" 
c
      subroutine read_rams_netcdf(nf_fid, x, y, z, glat, glon, ht, 
     +     ice, level, lwc, rai, aext, rtgt_ascii, sno, topt, topt_ascii, 
     +     ht_fcinv, ice_fcinv, lwc_fcinv, rai_fcinv, sno_fcinv)
c
      include 'netcdf.inc'
      integer x, y, z,nf_fid, nf_vid, nf_status
      integer ht_fcinv(z), ice_fcinv(z), lwc_fcinv(z), rai_fcinv(z),
     +     sno_fcinv(z)
      real glat( x, y), glon( x, y), ht( x,  y, z), ice( x,  y, z),
     +     level(z), lwc( x,  y, z), rai( x,  y, z), rtgt_ascii( x,
     +     y), sno( x,  y, z), topt( x, y), topt_ascii( x, y)



c   variables of type real
c
c     variable        netcdf long name
c     glat          
c
      nf_status=nf_inq_varid(nf_fid,'glat',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for glat'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,glat)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for glat'
       endif
      endif
c
c     variable        netcdf long name
c     glon          
c
      nf_status=nf_inq_varid(nf_fid,'glon',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for glon'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,glon)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for glon'
       endif
      endif
c
c     variable        netcdf long name
c     ht            
c
      nf_status=nf_inq_varid(nf_fid,'ht',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for ht'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,ht)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for ht'
       endif
      endif
c
c     variable        netcdf long name
c     ice           
c
      nf_status=nf_inq_varid(nf_fid,'ice',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for ice'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,ice)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for ice'
       endif
      endif
c
c     variable        netcdf long name
c     level         
c
      nf_status=nf_inq_varid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for level'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,level)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for level'
       endif
      endif
c
c     variable        netcdf long name
c     lwc           
c
      nf_status=nf_inq_varid(nf_fid,'lwc',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for lwc'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,lwc)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for lwc'
       endif
      endif
c
c     variable        netcdf long name
c     rai           
c
      nf_status=nf_inq_varid(nf_fid,'rai',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for rai'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,rai)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for rai'
       endif
      endif
c
c     variable        netcdf long name
c     aext
c
      nf_status=nf_inq_varid(nf_fid,'aext',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for aext'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,aext)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for aext'
       endif
      endif
c
c     variable        netcdf long name
c     rtgt_ascii    
c
      nf_status=nf_inq_varid(nf_fid,'rtgt_ascii',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for rtgt_ascii'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,rtgt_ascii)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for rtgt_ascii'
       endif
      endif
c
c     variable        netcdf long name
c     sno           
c
      nf_status=nf_inq_varid(nf_fid,'sno',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for sno'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,sno)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for sno'
       endif
      endif
c
c     variable        netcdf long name
c     topt          
c
      nf_status=nf_inq_varid(nf_fid,'topt',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for topt'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,topt)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for topt'
       endif
      endif
c
c     variable        netcdf long name
c     topt_ascii    
c
      nf_status=nf_inq_varid(nf_fid,'topt_ascii',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for topt_ascii'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,topt_ascii)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for topt_ascii'
       endif
      endif


c   variables of type double
c


c   variables of type char
c

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif
c
      return
      end

      subroutine stats_missing(array,ni,nj,nk,r_missing_data,nmiss)

      real array(ni,nj,nk)

      nmiss = 0

      do k = 1,nk
      do i = 1,ni
      do j = 1,nj
          if(array(i,j,k) .eq. r_missing_data)then
              nmiss = nmiss + 1
          endif
      enddo ! j
      enddo ! i
      enddo ! k

      return
      end
