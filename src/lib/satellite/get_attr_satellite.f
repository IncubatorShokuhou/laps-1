      subroutine get_attribute_wfo(nf_fid,centrallat,centrallon,
     &rlat00,rlon00,latnxny,lonnxny,latdxdy,londxdy,dx,dy,nx,ny,
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
      real      rlat00
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
c get centrallat
c
      nf_attid=0
      nf_attnum=0
      nf_status = nf_inq_attid(nf_fid,nf_attid,'centrallat',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'centrallat attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=nf_get_att_real(nf_fid,nf_attid,'centrallat',centrallat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'centrallat'
         istatus=-1
         goto 100
      endif
c
c get centrallon
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'centrallon',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'centrallon attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=nf_get_att_real(nf_fid,nf_attid,'centrallon',centrallon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'centrallon'
        istatus=-1
        goto 100
      endif
c
c get lat00
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'lat00',nf_attnum)
      if(nf_status.ne.nf_noerr)	then
         print*, nf_strerror(nf_status)
         print*, 'lat00 attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'lat00',rlat00)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'lat00'
         istatus=-1
         goto 100
      endif
c
c get lon00
c
      nf_status = nf_inq_attid(nf_fid,nf_attid,'lon00',nf_attnum)
      if(nf_status.ne.nf_noerr) then
         print*, nf_strerror(nf_status)
         print*, 'lon00 attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = nf_get_att_real(nf_fid,nf_attid,'lon00',rlon00)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'lon00'
        istatus=-1
        goto 100
      endif
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
c get y dimension
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
c
c============================================================
c
      function c_afwa_fname(csatid,chtype)

      character csatid*(*)
      character chtype*(*)
      character cs1*1
      character cs2*2
      character cs3*2
      character c_afwa_fname*(*)

      cs1='i'
      if(chtype.eq.'vis')cs1='v'

      if(csatid.eq.'meteos')then
         c_afwa_fname=csatid//cs1//'1_'//chtype
      else
         cs2=csatid(1:2)
         cs3=csatid(5:6)
         c_afwa_fname='u'//cs2//cs3//cs1//'1_'//chtype
      endif

      return 
      end

c
c ===========================================================
c
      function nw_vis_line_gwc(chtype,decimation,bescnfc,fsci)

      integer   nw_vis_line_gwc
      integer   fsci
      integer   bescnfc
      integer   decimation
      character chtype*3

      if(chtype.eq.'vis')then
         nw_vis_line_gwc=(bescnfc*decimation+fsci)
      else
         nw_vis_line_gwc=(bescnfc*decimation+fsci)*4 
      endif

      return
      end
c
c ===========================================================
c
      function nw_vis_pix_gwc(chtype,decimation,bepixfc,goalpha)

      integer    nw_vis_pix_gwc
      integer    bepixfc
      integer    decimation
      real       goalpha
      character  chtype*3

      if(chtype.eq.'vis')then
         nw_vis_pix_gwc=(bepixfc*decimation/4+int(goalpha))*8
      else
         nw_vis_pix_gwc=(bepixfc*decimation+int(goalpha))*8
      endif 

      return
      end
