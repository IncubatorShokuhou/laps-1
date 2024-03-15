
c *
c * to read header info of wideband image data
c *
      subroutine get_wideband_hdr(filename
     &          ,rad_dim, z_dim, v_dim
     &          ,sitelat,sitelon,sitealt,vcp
     &          ,istat)
c
      parameter (rad_max = 999, ang_max = 999)
c * output
      integer rad_dim, z_dim, v_dim, vcp
      real    sitelat,sitelon,sitealt
c*************************************
      character filename*72
c
      include 'netcdf.inc'
c
      status=nf_open(filename,nf_nowrite,ncid)
      if (status .ne. nf_noerr) then
           istat = -1
c           print *, filename,' is not avail'
           return
      endif
c
c * statements to fill nx and ny
c
      nf_status = nf_inq_dimid(ncid,'radial',nf_vid)
      nf_status = nf_inq_dimlen(ncid,nf_vid,rad_dim)
      nf_status = nf_inq_dimid(ncid,'z_bin',nf_vid)
      nf_status = nf_inq_dimlen(ncid,nf_vid,z_dim)
      nf_status = nf_inq_dimid(ncid,'v_bin',nf_vid)
      nf_status = nf_inq_dimlen(ncid,nf_vid,v_dim)

      nf_status = nf_inq_varid(ncid,'vcp',nf_vid)
      nf_status = nf_get_var_int(ncid,nf_vid,vcp)
      nf_status = nf_inq_varid(ncid,'sitelat',nf_vid)
      nf_status = nf_get_var_real(ncid,nf_vid,sitelat)
      nf_status = nf_inq_varid(ncid,'sitelon',nf_vid)
      nf_status = nf_get_var_real(ncid,nf_vid,sitelon)
      nf_status = nf_inq_varid(ncid,'sitealt',nf_vid)
      nf_status = nf_get_var_real(ncid,nf_vid,sitealt)
c
      nf_status = nf_close(ncid)
      istat = 0
c *
      return
      end
