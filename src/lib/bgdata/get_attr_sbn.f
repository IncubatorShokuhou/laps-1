      subroutine get_attribute_sbn(cdfname,centrallat,centrallon,
     &rlat00,rlon00,latnxny,lonnxny,latdxdy,londxdy,dx,dy,nx,ny,
     &rotation,projname,istatus)
c
c  open netcdf file for reading
c
      character cdfname*200
      character dummy*31
      character projname*30
      integer   istatus
      integer   dim_id
      integer   nx,ny
      integer   rcode
      integer   nf_fid
      integer   lenf,lenp
      integer   nf_status
      real      centrallat
      real      centrallon
      real      rlat00
      real      rlon00
      real      dx,dy
      real      latnxny
      real      lonnxny
      real      latdxdy
      real      londxdy
      real      rotation

      include 'netcdf.inc'

      print*,'in get_attribute_sbn'
      istatus = -1

      call s_len(cdfname,lenp)
      nf_status = nf_open(cdfname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', cdfname(1:lenp)
        return
      endif
c
c get centrallat
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'centrallat',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'centrallat attribute id'
         goto 100
      endif

      nf_status=nf_get_att_real(nf_fid,nf_attid,'centrallat',centrallat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'centrallat'
         goto 100
      endif
c
c get centrallon
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'centrallon',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'centrallon attribute id'
         goto 100
      endif

      nf_status=nf_get_att_real(nf_fid,nf_attid,'centrallon',centrallon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'centrallon'
        goto 100
      endif
c
c get lat00
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'lat00',nf_attnum)
      if(nf_status.ne.nf_noerr)	then
         print*, nf_strerror(nf_status)
         print*, 'lat00 attribute id'
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'lat00',rlat00)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'lat00'
         goto 100
      endif
c
c get lon00
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'lon00',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'lon00 attribute id'
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'lon00',rlon00)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'lon00'
        goto 100
      endif
c
c get latnxny
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'latnxny',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'latnxny attribute id'
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'latnxny',latnxny)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'latnxny'
         goto 100
      endif
c
c get lonnxny
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'lonnxny',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'lonnxny attribute id'
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'lonnxny',lonnxny)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'lonnxny'
        goto 100
      endif
c
c get latdxdy
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'latdxdy',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'latdxdy attribute id'
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'latdxdy',latdxdy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'latdxdy'
         goto 100
      endif
c
c get londxdy
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'londxdy',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'londxdy attribute id'
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'londxdy',londxdy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'londxdy'
        goto 100
      endif
c
c get resolution-x
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'dxkm',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'dxkm attribute id'
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'dxkm',dx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dxkm'
        goto 100
      endif
c
c get resolution-y
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'dykm',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'dykm attribute id'
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'dykm',dy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dykm'
        goto 100
      endif
c
c get x dimension
c
      dim_id = ncdid(nf_fid, 'x', rcode)
      if(rcode.ne.0)then
         write(6,*)'error getting x id code - returning'
         goto 100
      endif

      call ncdinq(nf_fid,dim_id,dummy,nx,rcode)
      if(rcode.ne.0)then
         write(6,*)'error getting x dimension - nx'
         goto 100
      endif
c
c get x dimension
c
      dim_id = ncdid(nf_fid, 'y', rcode)
      if(rcode.ne.0)then
         write(6,*)'error getting y id code - returning'
         goto 100
      endif

      call ncdinq(nf_fid,dim_id,dummy,ny,rcode)
      if(rcode.ne.0)then
         write(6,*)'error getting y dimension - ny'
         goto 100
      endif

c
c get resolution-y
c
	nf_status = nf_inq_attid(nf_fid,nf_attid,'rotation',nf_attnum)
	if(nf_status.ne.nf_noerr) then
	 print*, nf_strerror(nf_status)
	 print*, 'rotation: attribute id'
	 goto 100
	endif

	nf_status = nf_get_att_real(nf_fid,nf_attid,'rotation'
     .,rotation)
	if(nf_status.ne.nf_noerr) then
	print *, nf_strerror(nf_status)
	print *,'rotation'
	goto 100
	endif
c
c get projection name
c
        nf_status = nf_inq_attid(nf_fid,nf_attid,'projname',nf_attnum)
        if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'projname: attribute id'
         goto 100
        endif

        nf_status = nf_get_att_text(nf_fid,nf_attid,'projname'
     .,projname)
        if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'projname'
        goto 100
        endif

c
      istatus=1
100   return
      end
