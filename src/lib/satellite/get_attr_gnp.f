      subroutine get_attribute_gnp(nf_fid,centrallat,centrallon,
     &stdlat,stdlon,latnxny,lonnxny,latdxdy,londxdy,dx,dy,nx,ny,
     &istatus)
c
c  open netcdf file for reading 
c
      character dummy*31
      integer   istatus
      integer   dim_id
      integer   nx,ny
      integer   rcode
      integer   nf_fid
      integer   lenf
      integer   nf_status
      integer   nf_attid
      integer   nf_attnum
      real      centrallat
      real      centrallon
      real      stdlat
      real      stdlon
      real      dx,dy
      real      latnxny
      real      lonnxny
      real      latdxdy
      real      londxdy

      include 'netcdf.inc'

c     lenf=index(filename,' ')-1
c     print*,'opening ',filename(1:lenf)
c     nf_status = nf_open(filename,nf_nowrite,nf_fid)

c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'nf_open ',filename(1:lenf)
c       istatus=-1
c       goto 100
c     endif
c
c get variable id of projection parameters
      nf_status = nf_inq_varid(nf_fid,'lambert_projection',nf_projid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lambert_projection'
      endif

c get attributes, both global and from 'lambert_projection'      

c     
c get centrallat
c
      nf_attid=nf_global
      nf_status = nf_inq_attid(nf_fid,nf_attid,'tile_center_latitude'
     +                        ,nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'tile_center_latitude attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=nf_get_att_real(nf_fid,nf_attid,'tile_center_latitude'
     +                         ,centrallat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'tile_center_latitude'
         istatus=-1
         goto 100
      endif
c
c get centrallon
c
      nf_attid=nf_global
      nf_status = nf_inq_attid(nf_fid,nf_attid,'tile_center_longitude'
     +                        ,nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'tile_center_longitude attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=nf_get_att_real(nf_fid,nf_attid,'tile_center_longitude'
     +                         ,centrallon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'tile_center_longitude'
        istatus=-1
        goto 100
      endif
c
c get standard_parallel
c
      nf_attid=nf_projid
      nf_status = nf_inq_attid(nf_fid,nf_attid,'standard_parallel'
     +                        ,nf_attnum)
      if(nf_status.ne.nf_noerr)	then
         print*, nf_strerror(nf_status)
         print*, 'standard_parallel attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'standard_parallel'
     +                           ,stdlat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'standard_parallel'
         istatus=-1
         goto 100
      endif
c
c get stdlon
c
      nf_attid=nf_projid
      nf_status = nf_inq_attid(nf_fid,nf_attid
     +                       ,'longitude_of_central_meridian',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'longitude_of_central_meridian attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid
     +                          ,'longitude_of_central_meridian',stdlon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'longitude_of_central_meridian'
        istatus=-1
        goto 100
      endif
c
c get latnxny
c
      nf_attid=nf_global
      nf_status = nf_inq_attid(nf_fid,nf_attid,'latnxny',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'latnxny attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'latnxny',latnxny)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'latnxny'
         istatus=-1
         goto 100
      endif
c
c get lonnxny
c
      nf_attid=nf_global
      nf_status = nf_inq_attid(nf_fid,nf_attid,'lonnxny',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'lonnxny attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'lonnxny',lonnxny)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'lonnxny'
        istatus=-1
        goto 100
      endif
c
c get latdxdy
c
      nf_attid=nf_global
      nf_status = nf_inq_attid(nf_fid,nf_attid,'latdxdy',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'latdxdy attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'latdxdy',latdxdy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'latdxdy'
         istatus=-1
         goto 100
      endif
c
c get londxdy
c
      nf_attid=nf_global
      nf_status = nf_inq_attid(nf_fid,nf_attid,'londxdy',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'londxdy attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'londxdy',londxdy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'londxdy'
        istatus=-1
        goto 100
      endif
c
c get resolution-x
c
      nf_attid=nf_global
      nf_status = nf_inq_attid(nf_fid,nf_attid,'pixel_x_size',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'pixel_x_size attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'pixel_x_size',dx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'pixel_x_size'
        istatus=-1
        goto 100
      endif
c
c get resolution-y
c
      nf_attid=nf_global
      nf_status = nf_inq_attid(nf_fid,nf_attid,'pixel_y_size',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'pixel_y_size attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'pixel_y_size',dy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'pixel_y_size'
        istatus=-1
        goto 100
      endif
c
c get x dimension
c
      nf_attid=nf_global
      dim_id = ncdid(nf_fid, 'x', rcode)
      if(rcode.ne.0)then
         write(6,*)'error getting x id code - returning'
         istatus=-1
         goto 100
      endif

      call ncdinq(nf_fid,dim_id,dummy,nx,rcode)
      if(rcode.ne.0)then
         write(6,*)'error getting x dimension - nx'
         istatus=-1
         goto 100
      endif
c
c get y dimension
c
      nf_attid=nf_global
      dim_id = ncdid(nf_fid, 'y', rcode)
      if(rcode.ne.0)then
         write(6,*)'error getting y id code - returning'
         istatus=-1
         goto 100
      endif

      call ncdinq(nf_fid,dim_id,dummy,ny,rcode)
      if(rcode.ne.0)then
         write(6,*)'error getting y dimension - ny'
         istatus=-1
         goto 100
      endif
c
      istatus=0
100   return
      end
