
      subroutine get_wideband_netcdf(filename
     &          ,rad_dim,z_dim,v_dim
     &          ,image_z,image_v,image_w
     &          ,azim_ang, elev_ang, resolv
     &          ,istat)
c
c * a subroutine to read wideband level ii data
c   00 nov 2001  dongsoo kim    original
c 
      integer rad_max, ang_max
      parameter (rad_max = 999, ang_max = 999)
      integer rad_dim, z_dim, v_dim
      byte  z(z_dim,rad_dim),w(v_dim,rad_dim)
     &     ,v(v_dim,rad_dim)
c
c * output
c
      real  azim_ang(rad_dim), elev_ang(rad_dim), resolv
      integer*2  image_z(460,rad_max)
     &     ,image_v(920,rad_max)
     &     ,image_w(920,rad_max)
c
c *************************************
c
      character filename*72
c
      include 'netcdf.inc'
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
      nf_status = nf_inq_varid(ncid,'z',nf_vid)
      nf_status = nf_get_var_int1(ncid,nf_vid,z)
      nf_status = nf_inq_varid(ncid,'v',nf_vid)
      nf_status = nf_get_var_int1(ncid,nf_vid,v)
      nf_status = nf_inq_varid(ncid,'w',nf_vid)
      nf_status = nf_get_var_int1(ncid,nf_vid,w)
c
      nf_status = nf_inq_varid(ncid,'radialazim',nf_vid)
      nf_status = nf_get_var_real(ncid,nf_vid,azim_ang)
      nf_status = nf_inq_varid(ncid,'radialelev',nf_vid)
      nf_status = nf_get_var_real(ncid,nf_vid,elev_ang)
      nf_status = nf_inq_varid(ncid,'resolutionv',nf_vid)
      nf_status = nf_get_var_real(ncid,nf_vid,resolv)

      nf_status = nf_close(ncid)
c
      do j=1,rad_dim
      do i=1,z_dim
         image_z(i,j) = iand(z(i,j),255)
      enddo
      enddo
c
      do j=1,rad_dim
      do i=1,v_dim
         image_v(i,j) = iand(v(i,j),255)
         image_w(i,j) = iand(w(i,j),255)
      enddo
      enddo
c
      return
      end
