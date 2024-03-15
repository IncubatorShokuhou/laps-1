      subroutine get_attribute_vol(nf_fid,stationlatitude
     &  ,stationlongitude,stationelevationinmeters,station 
!    &  ,latnxny,lonnxny,latdxdy,londxdy,dx,dy,nx,ny
     &  ,istatus)
c
c  obtain attributes, still should add radar name
c
      character dummy*31
      character station*31
      integer   istatus
      integer   dim_id
      integer   nx,ny
      integer   rcode
      integer   nf_fid
      integer   lenf
      integer   nf_status
      integer   nf_attid
      integer   nf_attnum
      real      stationlatitude
      real      stationlongitude
      real      rstationelevationinmeters
      real      rlon00
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
c get stationlatitude
c
      nf_attid=0
      nf_attnum=0
      nf_status = nf_inq_attid(nf_fid,nf_attid,'stationlatitude'
     .                        ,nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'stationlatitude attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=nf_get_att_real(nf_fid,nf_attid,'stationlatitude'
     .                         ,stationlatitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'stationlatitude'
         istatus=-1
         goto 100
      endif
c
c get stationlongitude
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'stationlongitude'
     .                        ,nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'stationlongitude attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=nf_get_att_real(nf_fid,nf_attid,'stationlongitude'
     .                         ,stationlongitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'stationlongitude'
        istatus=-1
        goto 100
      endif
c
c get stationelevationinmeters
c
      nf_status = nf_inq_attid(nf_fid,nf_attid
     .                        ,'stationelevationinmeters',nf_attnum)
      if(nf_status.ne.nf_noerr)	then
         print*, nf_strerror(nf_status)
         print*, 'stationelevationinmeters attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid
     .                           ,'stationelevationinmeters'
     .                           ,stationelevationinmeters)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'stationelevationinmeters'
         istatus=-1
         goto 100
      endif

c
c get station     
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'station',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'station attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_text(nf_fid,nf_attid,'station',station)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'station'
        istatus=-1
        goto 100
      endif

      return
c
c get latnxny
c
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
      nf_status = nf_inq_attid(nf_fid,nf_attid,'dxkm',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'dxkm attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'dxkm',dx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dxkm'
        istatus=-1
        goto 100
      endif
c
c get resolution-y
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'dykm',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'dykm attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'dykm',dy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dykm'
        istatus=-1
        goto 100
      endif
c
c get x dimension
c
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
c get x dimension
c
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
