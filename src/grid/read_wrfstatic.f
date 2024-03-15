
      subroutine read_wrfstatic(ni,nj,lat,lon,filename,topo,istatus)
      include 'netcdf.inc'
      integer time, south_north, west_east, nf_fid, nf_vid, nf_status
      character*5 path
      character*255 filename
      data path/'data/'/

      real lat(ni,nj),lon(ni,nj),topo(ni,nj)
c
c  open netcdf file for reading
c
      istatus = 1

      call s_len(filename,lenp)
      write(6,*)' read_wrfstatic - file is ',filename(1:lenp)

      nf_status = nf_open(filename(1:lenp),nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open geo_em.d01.nc'
        istatus = 0
        return
      endif
c
c  fill all dimension values
c
c
c get size of time
c
      nf_status = nf_inq_dimid(nf_fid,'time',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim time'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,time)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim time'
      endif
c
c get size of south_north
c
      nf_status = nf_inq_dimid(nf_fid,'south_north',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim south_north'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,south_north)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim south_north'
      endif
c
c get size of west_east
c
      nf_status = nf_inq_dimid(nf_fid,'west_east',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim west_east'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,west_east)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim west_east'
      endif

      write(6,*)' wps dims are ',west_east,south_north
      if(west_east .ne. ni .or. south_north .ne. nj)then
          write(6,*)' error: wps dims differ from laps'
          write(6,*)
     1    ' only matching grids are currently supported: should be '
     1    ,ni,nj
          istatus = 0
          return
      endif

      call read_wrfstatic_sub(nf_fid , time, south_north, west_east
     1                       ,topo,lat,lon,istatus)                     

      return
      end
c
c
      subroutine read_wrfstatic_sub(nf_fid, time, south_north
     1                    ,west_east,hgt_m, xlat_m, xlong_m,istatus)
      include 'netcdf.inc'
      integer time, south_north, west_east, nf_fid, nf_vid, nf_status
      real hgt_m( west_east,  south_north, time), 
     +   xlat_m( west_east,  south_north, time), 
     +   xlong_m( west_east,  south_north, time)
      call read_netcdf(nf_fid , time, south_north, west_east,
     +    hgt_m, xlat_m, xlong_m)
c
c the netcdf variables are filled - your code goes here
c
      return
      end

      subroutine read_netcdf(nf_fid , time, south_north, west_east,
     +    hgt_m, xlat_m, xlong_m)
      include 'netcdf.inc'
      integer time, south_north, west_east, nf_fid, nf_vid, nf_status

      real hgt_m( west_east,  south_north, time), 
     +   xlat_m( west_east,  south_north, time), 
     +   xlong_m( west_east,  south_north, time)
c
c     variable        netcdf long name
c      hgt_m        
c
      nf_status = nf_inq_varid(nf_fid,'hgt_m',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var hgt_m'
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,hgt_m)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ hgt_m '
      endif
c
c     variable        netcdf long name
c      xlat_m       
c
      nf_status = nf_inq_varid(nf_fid,'xlat_m',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var xlat_m'
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,xlat_m)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ xlat_m '
      endif
c
c     variable        netcdf long name
c      xlong_m      
c
      nf_status = nf_inq_varid(nf_fid,'xlong_m',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var xlong_m'
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,xlong_m)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ xlong_m '
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
