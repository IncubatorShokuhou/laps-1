      subroutine get_conus3d_hdr(filename
     &          ,nx_dim, ny_dim, nz_dim
     &          ,istat)
c *
c * to read header info of conus3d netcdf format digital data
c * 4/20/04 dongsoo kim,  original
c
c * input
      character filename*132
c * output
      integer nx_dim, ny_dim, nz_dim, istat
c*************************************
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
c * statements to fill x, y and levels
c
      nf_status = nf_inq_dimid(ncid,'x',dimid)
      nf_status = nf_inq_dimlen(ncid,dimid,nx_dim)

      nf_status = nf_inq_dimid(ncid,'y',dimid)
      nf_status = nf_inq_dimlen(ncid,dimid,ny_dim)

      nf_status = nf_inq_dimid(ncid,'levels',dimid)
      nf_status = nf_inq_dimlen(ncid,dimid,nz_dim)
c *
      nf_status = nf_close(ncid)
c *
      return
      end
