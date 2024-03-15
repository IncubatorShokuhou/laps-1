      subroutine get_unidata_model_id(filename,cmodel,ivaltimes,ntbg
     &,istatus)

      implicit none
      include 'netcdf.inc'
      character*132 cmodel
      character*200 filename
      character*132 model
      integer ntbg,istatus
      integer ivaltimes(ntbg)
      integer nf_fid,nf_vid,nf_status
c
c  open netcdf file for reading
c
      istatus = 1
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', filename
        return
      endif

      nf_status = nf_inq_varid(nf_fid,'model',nf_vid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in var model'
         return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,model)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif
      nf_status=nf_inq_varid(nf_fid,'valtimeminusreftime',nf_vid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif
      nf_status=nf_get_vara_int(nf_fid,nf_vid,1,ntbg,ivaltimes)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      if(model(1:3).ne.cmodel(1:3))then
         print*,'mismatch between model and cmodel'
         print*,model,cmodel
      endif

      istatus = 0

      return
      end

      subroutine get_unidata_dims(cdfname,cmodel
     +,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv
     +,nzbg_ww,n_valtimes,istatus)

      implicit none
      include 'netcdf.inc'
      integer slen, nf_status,nf_fid, i, istat
      integer nf_vid
      character*132 cmodel
      character*200 cdfname

      integer nxbg,nybg
      integer nzbg_ht
      integer nzbg_tp
      integer nzbg_sh
      integer nzbg_uv
      integer nzbg_ww
      integer n_valtimes 
      integer record
      integer nvars
      integer ivaltimes(100)
      integer istatus

      integer ncid,itype,ndims
      integer j,k,kk,lc,nclen,lenc
      integer dimlen
      character*16 cvars(10)

      integer dimids(10)
      integer idims(10,10)
      integer nattr
      integer nf_attid,nf_attnum
      character*13 fname9_to_wfo_fname13, fname13
c
      istatus = 0
      call s_len(cdfname,slen)
c
c get size of n_valtimes
c

      call get_nvaltimes_unidata(cdfname,n_valtimes,ivaltimes,istatus)
      if(istatus.ne.1) then
         print *,'error: get_nvaltimes '
         return
      endif
c
c
c  open netcdf file for reading
c
      nf_status = nf_open(cdfname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open unidata'
      endif
c
c get everything for each variable
c
      call s_len(cmodel,nclen)
 
      nvars = 6 
      if ( (cmodel(1:nclen) .eq. 'ruc_iso') .or. 
     +     (cmodel(1:nclen) .eq. 'gfs_iso') ) then
        cvars(1)='z               '
        cvars(2)='rh              '
        cvars(3)='t               '
        cvars(4)='u               ' 
        cvars(5)='v               '
        cvars(6)='omega           '
      elseif(cmodel(1:nclen) .eq. 'ruc_hyb')  then
        cvars(1)='z_hybr          '
        cvars(2)='hum_mix_hybr    '
        cvars(3)='vptmp_hybr      '
        cvars(4)='u_hybr          '
        cvars(5)='v_hybr          '
        cvars(6)='omega_hybr      '
      else
        print *, "unsupported unidata model: ", cmodel
        istatus = 0
        return
      endif

      do i=1,nvars

         nf_status = nf_inq_varid(nf_fid, cvars(i),nf_vid)
         nf_status = nf_inq_var(nf_fid,nf_vid,cvars(i)
     +,itype,ndims,dimids,nattr)
         print *, cvars(i),nf_vid,ndims
         do j=1,ndims
            nf_status = nf_inq_dimlen(nf_fid,dimids(j),dimlen)
            idims(j,i)= dimlen
         enddo

         if(i.eq.1)then
            nxbg = idims(1,i)
            nybg = idims(2,i)
            nzbg_ht=idims(3,i)
         elseif(i.eq.3)then
            nzbg_tp=idims(3,i)
         elseif(i.eq.2)then
            nzbg_sh=idims(3,i)
         elseif(i.eq.4)then
            nzbg_uv=idims(3,i)
         elseif(i.eq.6)then
            nzbg_ww=idims(3,i)
         endif

      enddo

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 1
      return 
      end

      subroutine get_nvaltimes_unidata(cdfname,nvaltimes,ivaltimes,
     +    istatus)

      implicit none

      include 'netcdf.inc'

      integer       nf_fid,nf_status,nf_vid
      integer, intent(out)::       nvaltimes
      integer, intent(out)::    ivaltimes(100)
      integer       ireftimes(100)
      integer,allocatable ::       ilocaltimes(:)
      integer,intent(out)::   istatus
      integer t
      character*200,intent(in) :: cdfname

c     logical       l2

      istatus=0

c     l2=.false.
c     inquire(file=cdfname,opened=l2)
c     if(.not.l2)then

      nf_status = nf_open(cdfname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'nf_open: get_nvaltimes'
      endif

c     endif
         
c switched from n_valtimes to record to handle incomplete files
c lw 7-9-03
c     nf_status = nf_inq_dimid(nf_fid,'n_valtimes',nf_vid)
      nf_status = nf_inq_dimid(nf_fid,'record',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim n_valtimes'
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nvaltimes)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim n_valtimes'
        return
      endif
      nf_status=nf_inq_varid(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif
      allocate(ilocaltimes(nvaltimes))
      nf_status=nf_get_vara_int(nf_fid,nf_vid,1,nvaltimes,ilocaltimes)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif

      ! get the reftimes
      nf_status=nf_inq_varid(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif
      nf_status=nf_get_vara_int(nf_fid,nf_vid,1,nvaltimes,ireftimes)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif

      do t = 1, nvaltimes

        ivaltimes(t) = (ilocaltimes(t) - ireftimes(t)) * 3600
      enddo
      deallocate(ilocaltimes)
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'nf_close: get_nvaltimes'
         return
      endif

      istatus=1
      return
      end
c
c ------------------------------------------------------------
      subroutine read_unidata_iso(cdfname,af,cmodel,
     .nxbg,nybg,nzbght,nzbgtp,nzbgsh,nzbguv,nzbgww,
     .prbght,prbgsh,prbguv,prbgww,
     .ht,tp,sh,uw,vw,ww,
     .ht_sfc,pr_sfc,uw_sfc,vw_sfc,sh_sfc,tp_sfc,mslp,
     .ctype,istatus)
c
      implicit none
c
      include 'netcdf.inc'
      include 'bgdata.inc'

c     integer ncid, lenstr, ntp, nvdim, nvs, ndsize
      integer model_out
      integer ncid

c     model_out=1  => lga
c     model_out=2  => dprep

      integer ndims ,dimids(nf_max_var_dims)
      integer itype,nattr

      integer nxbg,nybg
      integer nzbght
      integer nzbgtp
      integer nzbgsh
      integer nzbguv
      integer nzbgww
      integer nzunidata
      integer ntbg
      integer rcode
      integer ivaltimes(100)
      integer ind2m, ind10m
      logical lcmpsfcq
c
      real, intent(out)  ::   mslp(nxbg,nybg)

c *** 3d output arrays.
c
      real, intent(out)  :: prbght(nxbg,nybg,nzbght)
      real, intent(out)  :: prbgsh(nxbg,nybg,nzbgsh)
      real, intent(out)  :: prbguv(nxbg,nybg,nzbguv)
      real, intent(out)  :: prbgww(nxbg,nybg,nzbgww)
      real, intent(out)  ::     ht(nxbg,nybg,nzbght)
      real, intent(out)  ::     tp(nxbg,nybg,nzbgtp)
      real, intent(out)  ::     sh(nxbg,nybg,nzbgsh)
      real, intent(out)  ::     uw(nxbg,nybg,nzbguv)
      real, intent(out)  ::     vw(nxbg,nybg,nzbguv)
      real, intent(out)  ::     ww(nxbg,nybg,nzbgww)
      real, intent(out)  ::     ht_sfc(nxbg,nybg)
      real, intent(out)  ::     tp_sfc(nxbg,nybg)
      real, intent(out)  ::     sh_sfc(nxbg,nybg)
      real, intent(out)  ::     uw_sfc(nxbg,nybg)
      real, intent(out)  ::     vw_sfc(nxbg,nybg)
      real, intent(out)  ::     pr_sfc(nxbg,nybg)
 
c
      real              ::  prbg(nzbght)

      integer start(10),count(10)
 
      integer i,j,k,n,ip,jp,ii,jj,it,kk
      integer istatus,slen,lent
      integer ibdht,ibdtp,ibduv,ibdsh,ibdww
c
      character*9   fname,oldfname,model
      character*5   ctype
      character*4   af
      character*16  cvar
      character*2   gproj
      character*200 cdfname
      character*132 cmodel
c
      real   xe,mrsat
      real   make_ssh

      integer nf_vid,nn,nf_status
      real cp,rcp, factor
      parameter (cp=1004.,rcp=287./cp)
c
c_______________________________________________________________________________
c
      interface

        subroutine read_netcdf_real(nf_fid,fname,n1,f
     +,start,count,istatus)
          integer n1
          integer nf_fid
          integer istatus
          integer start(10),count(10)
          real    f(n1)
          character*(*) fname
        end subroutine
      end interface
c
c -------------------------------------------------------

      print*,'here: read_unidata_iso'

      istatus = 1

      call s_len(cdfname,slen)

      print*,'cdfname: ',cdfname(1:slen)

      call get_nvaltimes_unidata(cdfname,ntbg,ivaltimes,istatus)
c
      print*,'opening cdf file: ',cdfname(1:slen)

      rcode = nf_open(cdfname,nf_nowrite,ncid)
      if(rcode.ne.nf_noerr) then
         print *, nf_strerror(rcode)
         print *,'nf_open ',cdfname(1:slen)
         return
      endif

      read(af,'(i4)') nn

      n=1
      do while(n.lt.ntbg.and.ivaltimes(n)/3600.ne. nn)
         n=n+1
      enddo
      if(ivaltimes(n)/3600.ne.nn) then

         print*,'error: no record valid at requested time '
         print*,'ntbg/nn/af/n/ivaltimes(n) ',ntbg,' ',nn,' ',af,
     &' ',n,' ',ivaltimes(n)

         rcode= nf_close(ncid)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var: ',cmodel
            return
         endif

         goto 999

      else

         print*,'found valid record at ivaltime'
         print*,'ntbg/nn/af/n/ivaltimes(n) ',ntbg,' ',nn,' ',af,
     &' ',n,' ',ivaltimes(n) 
         print*
      endif

      ! get the pressure levels for this data
      nf_status = nf_inq_varid(ncid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in level '
        istatus = 0
        return
      endif
      nf_status = nf_get_var_real(ncid,nf_vid,prbg)
     
      ! get index for 2m and 10m winds
      ind2m = 1
      ind10m = 2
  
      ! set some indices
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbght
      start(4)=n
      count(4)=1

      print*,'read ht'
      cvar='z'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),ht,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (z): ',cmodel
         else
            print *,'missing ht data detected: return'
         endif
         print*
         return
      endif

c
c ****** statements to fill tp.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgtp
      start(4)=n
      count(4)=1

      print*,'read tp'
      cvar='t'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),tp,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (t): ',cmodel
         else
            print *,'missing t data detected: return'
         endif
         print*
         return
      endif

c
c ****** statements to fill rh.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgsh
      start(4)=n
      count(4)=1

      print*,'read rh'
      cvar='rh'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),sh,start
     +     ,count,rcode)
     
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (rh): ',cmodel
         else
c          rh is missing above 100 mb in gfs
           if (cmodel .eq. 'gfs_iso') then
            print *, "filling rh at top levels!"
             do j=1,nybg
             do i=1,nxbg
             do k = 1,nzbgsh
               if ((sh(i,j,k) .lt. 0.).or. 
     +             (sh(i,j,k) .gt. 200.)) then
                 sh(i,j,k) = 1.0
               endif
             enddo
             enddo
             enddo
           else
             print *,'missing rh data detected: return'
             return
           endif
         endif
      endif
c
c ****** statements to fill uw. 
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv
      start(4)=n
      count(4)=1

      print*,'read uw'
      cvar='u'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),uw,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (uw): ',cmodel
         else
            print *,'missing u data detected: return'
         endif
         print*
         return
      endif

c
c ****** statements to fill vw.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv
      start(4)=n
      count(4)=1

      print*,'read vw'
      cvar='v'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),vw,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (vw): ',cmodel
         else
            print *,'missing v data detected: return'
         endif
         print*
         return
      endif
c
c ****** statements to fill ww.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgww
      start(4)=n
      count(4)=1

      print*,'read ww'
      cvar='omega'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),ww  
     +  ,start,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (ww): ',cmodel
         else
c           ww is missing above 100 mb in gfs
           if (cmodel .eq. 'gfs_iso') then
             do j=1,nybg
             do i=1,nxbg
             do k = 1,nzbgww
               if ((ww(i,j,k) .lt. -1000.) .or.
     +             (ww(i,j,k) .gt. 1000)) then
                 ww(i,j,k) = 0.
               endif
             enddo
             enddo
             enddo
           else
             print *,'missing ww data detected: return'
             print *,'missing ww data detected: continue without'
             print *,'filling ww with 0.0'
             ww(:,:,:) = 0.0
           endif
         endif
      endif

c   get 2m t and rh, 10m u and v
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind2m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read tp_sfc'
      lcmpsfcq=.true.
      cvar='t_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,tp_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (tp_sfc): ',cmodel
         endif
      endif
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read ht_sfc'
      lcmpsfcq=.true.
      cvar='z_sfc'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,ht_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (ht_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind2m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read sh_sfc'
      lcmpsfcq=.true.
      cvar='rh_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,sh_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (sh_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind10m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read uw_sfc'
      lcmpsfcq=.true.
      cvar='u_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,uw_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (uw_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind10m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read vw_sfc'
      lcmpsfcq=.true.
      cvar='v_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,vw_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (vw_sfc): ',cmodel
         endif
      endif

c
c get sfc pressure field
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=1
      start(4)=n
      count(4)=1

      print*,'read p'

      lcmpsfcq=.true.
      
      cvar='p_sfc'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,pr_sfc,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (p): ',cmodel
         endif
         print *,'missing sfc p data detected'
         print*,' -> continue without; compute in sfcbkgd'
         lcmpsfcq=.false.
         print*
c        return
      endif
c
c get mslp (this field name differs from one model to the other)
c
      if(cmodel(1:3).eq.'ruc')then
         cvar='pm_msl'
      else
         cvar='p_msl'
      endif
      call read_netcdf_real(ncid,cvar,nxbg*nybg,mslp
     +           ,start,count,rcode)
      if(rcode.ne.nf_noerr) then
        if(rcode.gt.-61)then
          print *, nf_strerror(rcode)
          print *,'in nf_get_var (emsp): ',cmodel
        else
          print*,'error status returned from read_netcdf_real'
        endif
        print *,'missing emsp data detected'
        print*
        return
      endif

c
c *** close netcdf file.
c
      rcode= nf_close(ncid)
      if(rcode.ne.nf_noerr) then
         print *, nf_strerror(rcode)
         print *,'in nf_get_var: ',cmodel
         return
      endif
c
ccc      endif

c
c *** fill ouput arrays.
c *** convert rh to sh.
c
      print*,'load prbg arrays'
      do j=1,nybg
      do i=1,nxbg
         prbgsh(i,j,:)=prbg(:)
         prbght(i,j,:)=prbg(:)
         prbguv(i,j,:)=prbg(:)
         prbgww(i,j,:)=prbg(:)
      enddo
      enddo

c for laps-lgb

      call s_len(ctype,lent)

      print*,'ctype ',ctype(1:lent)

c since we may not have pr_sfc at this point lets not
c do this comp here and instead let sfcbkgd do it later.
      if(ctype(1:lent).eq.'lapsb' .and. lcmpsfcq)then

         print*,' computing sfc q '
         ibdtp=0
         ibduv=0
         do j=1,nybg
         do i=1,nxbg
            if(tp_sfc(i,j).lt.missingflag.and.
     .         tp_sfc(i,j).gt.150.)       then
c
c make sfc q from rh
c
              it=int(tp_sfc(i,j)*100)
              it=min(45000,max(15000,it))
              xe=esat(it)
              mrsat=0.00622*xe/(pr_sfc(i,j)*0.01-xe)
              sh_sfc(i,j)=sh_sfc(i,j)*mrsat
              sh_sfc(i,j)=sh_sfc(i,j)/(1.+sh_sfc(i,j))
            else
              ibdtp=ibdtp+1
            endif
            if(uw_sfc(i,j).gt.500..or.vw_sfc(i,j).gt.500.)then
              ibduv=ibduv+1
            endif
         enddo
         enddo

         print*,'done computing sfc q'

      endif


      if(ibdtp.gt.0.or.ibduv.gt.0)then
         print*,'found bad sfc data (tp/uv) ',ibdtp,ibduv
         return
      endif


      ibdht=0
      ibdtp=0
      ibduv=0
      do i=1,nxbg
      do j=1,nybg

         do k=1,nzbght              
            if (ht(i,j,k) .ge. 99999.)then
                ibdht=ibdht+1
            endif
         enddo
         do k=1,nzbgtp
            if(tp(i,j,k).lt.100.or.tp(i,j,k).ge.missingflag)then
               ibdtp=ibdtp+1
            endif
         enddo
         do k=1,nzbguv
            if(abs(uw(i,j,k)).ge. 500. .or.
     .         abs(vw(i,j,k)).ge. 500.)then
                ibduv=ibduv+1
            endif
         enddo
      enddo
      enddo

      if(ibdht.gt.0.or.ibdtp.gt.0.or.ibduv.gt.0)then
         print*,'found bad 3d data (ht/t/uv) ',ibdht,ibdtp,ibduv
         print*,'return to read_bgdata'
         return
      endif

      ibdsh=0
      if(ctype(1:lent).eq.'lapsb')then
         print*,'computing 3d sh'
         print*,'nzbgsh/nzbgtp = ',nzbgsh,nzbgtp
         do j=1,nybg
         do i=1,nxbg
         do k=1,nzbgsh
           if (sh(i,j,k).lt.200)then
             sh(i,j,k)= make_ssh(prbgsh(i,j,k)
     .,tp(i,j,k)-273.15,sh(i,j,k)/100., -132.)*0.001  !kg/kg 

           else
            ibdsh=ibdsh+1
           endif
         enddo
         enddo
         enddo

      endif

      print*,'done computing 3d sh'

      if(ibdsh.gt.0)then
         print*,'found bad rh 3d data',ibdsh
         print*,'return to read_bgdata'
         return
      endif


c this for dprep ... ingnore for now!
c -----------------------------------
      if(.false. .and. model_out.eq.2) then
c
c compute exner and convert temp to theta
c
         do k=1,nzbght
            do j=1,nybg
               do i=1,nxbg
                  if(tp(i,j,k).ne.missingflag) then
                     factor=(1000./prbght(i,j,k))**rcp
                     tp(i,j,k) = tp(i,j,k)*factor
                     prbght(i,j,k) = cp/factor
                  endif
               enddo
            enddo
         enddo

      endif


c     if(istatus_211 .eq. 0) then
c       print*, 'no valid data found for',fname, af
c       return
c     endif

c      istatus = 1

      if(0.eq.1) then
 900     print*,'error: bad dimension specified in netcdf file'
         print*, (count(i),i=1,4)
         istatus=-1
      endif
 999  return
      end

c -----------------------------------------------------------

      subroutine get_unidata_grid(filename,cmodel,nx,ny
     &,stdlat1,stdlat2,lov,la1,lo1,la2,lo2,dx,dy,gproj, istatus)
c
      implicit none
      include 'netcdf.inc'

      real    stdlat1,stdlat2
      real    lov
      real    la1,lo1
      real    la2,lo2
      real    dx, dy
      character*2 gproj
      character*100 grid_type
      integer  nx, ny, nf_fid, nf_vid, nf_status
      character filename*200
      character cmodel*(*)
      character*2, dimname
      integer istatus

      istatus = 0
c
c  open netcdf file for reading
c
      print *, "opening:  ", trim(filename)
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', filename
        return
      endif
c
c     variable        netcdf long name
c     grid_type
c
      nf_status = nf_inq_varid(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_type'
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,grid_type)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_type'
      endif
                                                                            
      if (grid_type(1:7) .eq. "lambert") then
        gproj = "lc"
      elseif(grid_type(1:18) .eq. "latitude/longitude") then
        gproj = "ll"
      else
        print *, "found unsupported grid_type:",trim(grid_type)
        print *, "in get_unidata_grid"
        istatus = 0
        return
      endif

c     nx
      if (gproj .ne. "ll") then
        dimname = "nx"
      else
        dimname = "ni"
      endif
      nf_status = nf_inq_varid(nf_fid,dimname,nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ',dimname
        istatus = 0
        return
      else  
        nf_status = nf_get_var_int(nf_fid,nf_vid,nx)
      endif

c     ny
      if (gproj .ne. "ll") then
        dimname = "ny"
      else
        dimname = "nj"
      endif

      nf_status = nf_inq_varid(nf_fid,dimname,nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ',dimname
        istatus = 0
        return
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,ny)
      endif


c
c     variable        netcdf long name
c      la1          "first latitude"
c
      nf_status = nf_inq_varid(nf_fid,'la1',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la1'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,la1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la1'
      endif
c
c     variable        netcdf long name
c      lo1          "first longitude"
c
        nf_status = nf_inq_varid(nf_fid,'lo1',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo1'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lo1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo1'
      endif
      if (lo1 .gt. 180.) lo1 = lo1 - 360.
c
c     variable        netcdf long name
c      lov          "orientation of grid"
c
      if (gproj .ne. "ll") then
        nf_status = nf_inq_varid(nf_fid,'lov',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lov'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,lov)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lov - lov'
        endif
        if (lov .gt. 180.) lov = lov - 360.
c
ci      variable        netcdf long name
c       latin1
c
        nf_status = nf_inq_varid(nf_fid,'latin1',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var latin1'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,stdlat1)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var latin1 - stdlat1'
        endif
c
c       variable        netcdf long name
c         latin2
c
        nf_status = nf_inq_varid(nf_fid,'latin2',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var latin2'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,stdlat2)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var latin2 - stdlat2'
        endif
      endif
c
c     variable        netcdf long name
c     dx 
c
      if (gproj .ne. "ll") then
        dimname = "dx"
      else 
        dimname = "di"
      endif 
      nf_status = nf_inq_varid(nf_fid,dimname,nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dx'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,dx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dx'
      endif
c
c     variable        netcdf long name
c     dy
c
      if (gproj .ne. "ll") then
        dimname = "dy"
      else
        dimname = "dj"
      endif
      nf_status = nf_inq_varid(nf_fid,dimname,nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dy'
      endif
       nf_status = nf_get_var_real(nf_fid,nf_vid,dy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dy'
      endif

      ! set or compute corner points

      if (cmodel(1:3) .eq. "ruc") then
        la2=55.481
        lo2=-57.381
      elseif(gproj .eq. "ll") then
c
c       variable        netcdf long name
c       la2          "last latitude"
c
        nf_status = nf_inq_varid(nf_fid,'la2',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var la2'
        endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,la2)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var la2'
        endif
c
c       variable        netcdf long name
c       lo2          "last longitude"
c
        nf_status = nf_inq_varid(nf_fid,'lo2',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lo1'
        endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lo2)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lo1'
        endif
        if (lo2 .gt. 180.) lo2 = lo2 - 360.

      else
        print *, "need to add computation of la2,lo2"
        print *, " in get_unidata_grid"
        istatus = 0
        return
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        return
      endif

      istatus=1 
      return
      end
c
c ------------------------------------------------------------
      subroutine read_unidata_ruc_hyb(cdfname,af,cmodel,
     .nxbg,nybg,nzbght,nzbgtp,nzbgsh,nzbguv,nzbgww,
     .prbght,prbgsh,prbguv,prbgww,
     .ht,tp,sh,uw,vw,ww,
     .ht_sfc,pr_sfc,uw_sfc,vw_sfc,sh_sfc,tp_sfc,mslp,
     .ctype,istatus)
c
      implicit none
c
      include 'netcdf.inc'
      include 'bgdata.inc'

c     integer ncid, lenstr, ntp, nvdim, nvs, ndsize
      integer model_out
      integer ncid

c     model_out=1  => lga
c     model_out=2  => dprep

      integer ndims ,dimids(nf_max_var_dims)
      integer itype,nattr

      integer nxbg,nybg
      integer nzbght
      integer nzbgtp
      integer nzbgsh
      integer nzbguv
      integer nzbgww
      integer nzunidata
      integer ntbg
      integer rcode
      integer ivaltimes(100)
      integer ind2m, ind10m
      logical lcmpsfcq
c
      real, intent(out)  ::   mslp(nxbg,nybg)

c *** 3d output arrays.
c
      real, intent(out)  :: prbght(nxbg,nybg,nzbght)
      real, intent(out)  :: prbgsh(nxbg,nybg,nzbgsh)
      real, intent(out)  :: prbguv(nxbg,nybg,nzbguv)
      real, intent(out)  :: prbgww(nxbg,nybg,nzbgww)
      real, intent(out)  ::     ht(nxbg,nybg,nzbght)
      real, intent(out)  ::     tp(nxbg,nybg,nzbgtp)
      real, intent(out)  ::     sh(nxbg,nybg,nzbgsh)
      real, intent(out)  ::     uw(nxbg,nybg,nzbguv)
      real, intent(out)  ::     vw(nxbg,nybg,nzbguv)
      real, intent(out)  ::     ww(nxbg,nybg,nzbgww)
      real, intent(out)  ::     ht_sfc(nxbg,nybg)
      real, intent(out)  ::     tp_sfc(nxbg,nybg)
      real, intent(out)  ::     sh_sfc(nxbg,nybg)
      real, intent(out)  ::     uw_sfc(nxbg,nybg)
      real, intent(out)  ::     vw_sfc(nxbg,nybg)
      real, intent(out)  ::     pr_sfc(nxbg,nybg)
 
c
      integer start(10),count(10)
 
      integer i,j,k,n,ip,jp,ii,jj,it,kk
      integer istatus,slen,lent
      integer ibdht,ibdtp,ibduv,ibdsh,ibdww
c
      character*9   fname,oldfname,model
      character*5   ctype
      character*4   af
      character*16  cvar
      character*2   gproj
      character*200 cdfname
      character*132 cmodel
c
      real tv
      integer nf_vid,nn,nf_status,ic,jc
      real cp,rcp, factor
      parameter (cp=1004.,rcp=287./cp)
c
c_______________________________________________________________________________
c
      interface

        subroutine read_netcdf_real(nf_fid,fname,n1,f
     +,start,count,istatus)
          integer n1
          integer nf_fid
          integer istatus
          integer start(10),count(10)
          real    f(n1)
          character*(*) fname
        end subroutine
      end interface
c
c -------------------------------------------------------

      print*,'here: read_unidata_ruc_hyb'

      istatus = 1

      call s_len(cdfname,slen)

      print*,'cdfname: ',cdfname(1:slen)

      call get_nvaltimes_unidata(cdfname,ntbg,ivaltimes,istatus)
c
      print*,'opening cdf file: ',cdfname(1:slen)

      rcode = nf_open(cdfname,nf_nowrite,ncid)
      if(rcode.ne.nf_noerr) then
         print *, nf_strerror(rcode)
         print *,'nf_open ',cdfname(1:slen)
         return
      endif

      read(af,'(i4)') nn

      n=1
      do while(n.lt.ntbg.and.ivaltimes(n)/3600.ne. nn)
         n=n+1
      enddo
      if(ivaltimes(n)/3600.ne.nn) then

         print*,'error: no record valid at requested time '
         print*,'ntbg/nn/af/n/ivaltimes(n) ',ntbg,' ',nn,' ',af,
     &' ',n,' ',ivaltimes(n)

         rcode= nf_close(ncid)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var: ',cmodel
            return
         endif

         goto 999

      else

         print*,'found valid record at ivaltime'
         print*,'ntbg/nn/af/n/ivaltimes(n) ',ntbg,' ',nn,' ',af,
     &' ',n,' ',ivaltimes(n) 
         print*
      endif

      ! get the pressure levels for this data
      nf_status = nf_inq_varid(ncid,'p_hybr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in level '
        istatus = 0
        return
      endif
      nf_status = nf_get_var_real(ncid,nf_vid,prbght)
 
      ! get index for 2m and 10m winds
      ind2m = 1
      ind10m = 2
  
      ! set some indices
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbght
      start(4)=n
      count(4)=1

      print*,'read ht'
      cvar='z_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),ht,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (z): ',cmodel
         else
            print *,'missing ht data detected: return'
         endif
         print*
         return
      endif

c
c ****** statements to fill tp.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgtp
      start(4)=n
      count(4)=1

      print*,'read vptmp'
      cvar='vptmp_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),tp,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (t): ',cmodel
         else
            print *,'missing t data detected: return'
         endif
         print*
         return
      endif

c
c ****** statements to fill sh.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgsh
      start(4)=n
      count(4)=1

      print*,'read qv'
      cvar='hum_mix_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),sh,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (qv): ',cmodel
         else
            print *,'missing qv data detected: return'
         endif
         print*
         return
      endif
c
c ****** statements to fill uw. 
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv
      start(4)=n
      count(4)=1

      print*,'read uw'
      cvar='u_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),uw,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (uw): ',cmodel
         else
            print *,'missing u data detected: return'
         endif
         print*
         return
      endif

c
c ****** statements to fill vw.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv
      start(4)=n
      count(4)=1

      print*,'read vw'
      cvar='v_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),vw,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (vw): ',cmodel
         else
            print *,'missing v data detected: return'
         endif
         print*
         return
      endif
c
c ****** statements to fill ww.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgww
      start(4)=n
      count(4)=1

      print*,'read omega'
      cvar='omega_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),ww  
     +  ,start,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (ww): ',cmodel
         else
            print *,'missing ww data detected: continue without'
            print *,'filling ww with 0.0'
            ww(:,:,:) = 0.0
         endif
         print*
      endif

c   get 2m t and rh, 10m u and v
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind2m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read t_2m'
      lcmpsfcq=.true.
      cvar='t_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,tp_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (tp_sfc): ',cmodel
         endif
      endif

c  no surface height...so use lowest hybrid level heigh
      ht_sfc = ht(:,:,1)

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind2m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read sh_sfc'
      lcmpsfcq=.true.
      cvar='hum_mix_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,sh_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (sh_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind10m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read uw_sfc'
      lcmpsfcq=.true.
      cvar='u_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,uw_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (uw_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind10m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read vw_sfc'
      lcmpsfcq=.true.
      cvar='v_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,vw_sfc,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (vw_sfc): ',cmodel
         endif
      endif

c
c get sfc pressure field (use lowest hybrid level)
c
c
      pr_sfc = prbght(:,:,1)

c get mslp (this field name differs from one model to the other)
c
      if(cmodel(1:3).eq.'ruc')then

         cvar='pm_msl'
         call read_netcdf_real(ncid,cvar,nxbg*nybg,mslp
     +           ,start,count,rcode)
         if(rcode.ne.nf_noerr) then
            if(rcode.gt.-61)then
               print *, nf_strerror(rcode)
               print *,'in nf_get_var (emsp): ',cmodel
            else
               print*,'error status returned from read_netcdf_real'
            endif
            print *,'missing emsp data detected'
            print*
            return
         endif

      endif

c
c *** close netcdf file.
c
      rcode= nf_close(ncid)
      if(rcode.ne.nf_noerr) then
         print *, nf_strerror(rcode)
         print *,'in nf_get_var: ',cmodel
         return
      endif
c
ccc      endif

c
c *** fill ouput arrays.
c ***   convert theta-v to temperature
c ***   convert mixing ratio to specific humidity
c ***   convert 3d pressure level arrays to from pa to mb
c ***   load remaining pressure arrays
c
c for laps-lgb

      call s_len(ctype,lent)

      print*,'ctype ',ctype(1:lent)
     
      do j=1,nybg
        do i=1,nxbg
           do k = 1, nzbght

c             conver pressure to mb and fill additional arrays
              prbght(i,j,k) = prbght(i,j,k) * 0.01
              prbgsh(i,j,k) = prbght(i,j,k)
              prbguv(i,j,k) = prbght(i,j,k)
              prbgww(i,j,k) = prbght(i,j,k)

c             convert theta-v into virtual temp
              tv = tp(i,j,k) * (prbgsh(i,j,k)*0.001)**rcp
c             convert virtual temp into temp
              tp(i,j,k) = tv/(1.0+0.61*sh(i,j,k))

c             convert mixrat into spechum
              sh(i,j,k) = sh(i,j,k)/(1.+sh(i,j,k))

           enddo
c          also convert surface mixrat into spechum
           sh_sfc(i,j)=sh_sfc(i,j)/(1.+sh_sfc(i,j))
         enddo
       enddo

      ic = nxbg/2
      jc = nybg/2
      do k = 1, nzbght
        print 850, k,prbght(ic,jc,k),ht(ic,jc,k),tp(ic,jc,k), 
     +              sh(ic,jc,k),uw(ic,jc,k),vw(ic,jc,k),ww(ic,jc,k)
      enddo
! honglig jiang: change from f7.5 to f8.5. 11/27/2013
 850  format(i2,1x,f6.1,1x,f7.1,1x,f5.1,1x,f8.5,1x,f5.1,1x,f5.1,1x,
     +   f8.5)
    
      istatus =  1
      if(0.eq.1) then
 900     print*,'error: bad dimension specified in netcdf file'
         print*, (count(i),i=1,4)
         istatus=-1
      endif
 999  return
      end
