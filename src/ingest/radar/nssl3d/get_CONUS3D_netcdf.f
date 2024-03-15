      subroutine get_conus3d_netcdf(filename
     &          ,nx_dim,ny_dim,nz_dim
     &          ,xxmin,xxmax,yymin,yymax, mrefl, istat) 
c *
c * a subroutine to read 1 km resolution nssl 3d in netcdf format
c * 4/19/2004 dongsoo kim, original ingest routine
c *
c * input
      character filename*132
      integer*4 nx_dim, ny_dim, nz_dim
c * output
      real*4   mrefl(nx_dim, ny_dim, nz_dim)
      real*4   xxmin, xxmax, yymin, yymax
      integer  ncid
c
      include '/usr/local/apps/netcdf-3.4/include/netcdf.inc'
c
      status=nf_open(filename,nf_nowrite,ncid)
      if (status .ne. nf_noerr) then
           istat = -1
           print *, filename,' is not avail'
           return
      endif
c
c  statements to fill image                          
c
      nf_status = nf_inq_varid(ncid,'mrefl_mosaic',nf_vid)
      nf_status = nf_get_var_real(ncid,nf_vid,mrefl)
c
c  read global attributes
c
      nf_status = nf_get_att_real(ncid,nf_global,'xmin', xxmin)
      nf_status = nf_get_att_real(ncid,nf_global,'xmax', xxmax)
      nf_status = nf_get_att_real(ncid,nf_global,'ymin', yymin)
      nf_status = nf_get_att_real(ncid,nf_global,'ymax', yymax)
c
      nf_status = nf_close(ncid)
c
      return
      end
