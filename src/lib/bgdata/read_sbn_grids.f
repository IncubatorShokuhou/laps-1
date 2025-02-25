      subroutine get_sbn_model_id(filename,cmodel,ivaltimes,ntbg
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

      subroutine get_sbn_dims(cdfname,cmodel
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

c     integer ntp, nvdim, nvs, lenstr, ndsize
c     integer ntp, nvdim, nvs
c     character*31 dummy
      
c     integer id_fields(5), vdims(10)
c     data id_fields/1,4,7,10,13/

      integer ncid,itype,ndims
      integer j,k,kk,lc,nclen,lenc
      integer dimlen
      character*10 cvars(10)

      integer dimids(10)
      integer idims(10,10)
      integer nattr
      integer nf_attid,nf_attnum
      character*13 fname9_to_wfo_fname13, fname13
c     linda wharton 10/27/98 removed several commented out lines:
c        print *,'ndsize = ', ndsize
c        value ndsize not set anywhere in this subroutine
c
      istatus = 0
      call s_len(cdfname,slen)
c
c get size of n_valtimes
c

      call get_nvaltimes(cdfname,n_valtimes,ivaltimes,istatus)
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
        print *,'nf_open rucsbn'
      endif
c
c get size of record
c
c     nf_status = nf_inq_dimid(nf_fid,'record',nf_vid)
c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'dim record'
c       return
c     endif
c     nf_status = nf_inq_dimlen(nf_fid,nf_vid,record)
c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'dim record'
c       return
c     endif
c
c get size of x
c
c     nf_status = nf_inq_dimid(nf_fid,'x',nf_vid)
c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'dim x'
c       return
c     endif
c     nf_status = nf_inq_dimlen(nf_fid,nf_vid,nxbg)
c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'dim x'
c       return
c     endif
c
c get size of y
c
c     nf_status = nf_inq_dimid(nf_fid,'y',nf_vid)
c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'dim y'
c       return
c     endif
c     nf_status = nf_inq_dimlen(nf_fid,nf_vid,nybg)
c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'dim y'
c       return
c     endif
c
c get everything for each variable
c
      call s_len(cmodel,nclen)

      nvars=8
      cvars(1)='gh'
      cvars(2)='rh'
      cvars(3)='t'
      cvars(4)='uw'
      cvars(5)='vw'
      cvars(6)='pvv'
      cvars(7)='p'  !sfc pressure
      if(cmodel(1:nclen).eq.'ruc40_native')then
         cvars(8)='mmsp'
      elseif(cmodel(1:nclen).eq.'eta48_conus'.or.
     .       cmodel(1:7).eq.'mesoeta')then
         cvars(8)='emsp'
      elseif(cmodel(1:nclen).eq.'avn_sbn_cyleq')then
         cvars(8)='pmsl'
      endif

      do i=1,nvars

         nf_status = nf_inq_varid(nf_fid, cvars(i),nf_vid)
         nf_status = nf_inq_var(nf_fid,nf_vid,cvars(i)
     +,itype,ndims,dimids,nattr)

         do j=1,ndims
            nf_status = nf_inq_dimlen(nf_fid,dimids(j),dimlen)
            idims(j,i)= dimlen
         enddo

c        nf_status = nf_inq_attid(nf_fid,nf_vid,'_n3d',nf_attnum)
c        if(nf_status.ne.nf_noerr) then
c           print*, nf_strerror(nf_status)
c           print*, 'attribute id ', cvars(i)
c           return
c        endif

         nf_status=nf_get_att_int(nf_fid,nf_vid,'_n3d',idims(3,i))
         if(nf_status.ne.nf_noerr) then
            print *, nf_strerror(nf_status)
            print *,'get attribute ',cvars(i)
            return
         endif


         if(cvars(i).eq.'gh')then
            nxbg = idims(1,i)
            nybg = idims(2,i)
            nzbg_ht=idims(3,i)
         elseif(cvars(i).eq.'t')then
            nzbg_tp=idims(3,i)
         elseif(cvars(i).eq.'rh')then
            nzbg_sh=idims(3,i)
         elseif(cvars(i).eq.'uw')then
            nzbg_uv=idims(3,i)
         elseif(cvars(i).eq.'pvv')then
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

      subroutine get_nvaltimes(cdfname,nvaltimes,ivaltimes,istatus)

      implicit none

      include 'netcdf.inc'

      integer       nf_fid,nf_status,nf_vid
      integer       nvaltimes
      integer       ivaltimes(100)
      integer       istatus
      character*200 cdfname

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
      nf_status=nf_inq_varid(nf_fid,'valtimeminusreftime',nf_vid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif
      nf_status=nf_get_vara_int(nf_fid,nf_vid,1,nvaltimes,ivaltimes)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ model '
         return
      endif

c     if(.not.l2)then

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
      subroutine read_sbn_grids(cdfname,af,cmodel,
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
      integer nzsbn
      integer ntbg
      integer rcode
      integer ivaltimes(100)

      logical lcmpsfcq
c
c *** sfc output arrays.
c
      real, intent(out)  :: pr_sfc(nxbg,nybg)
      real, intent(out)  :: uw_sfc(nxbg,nybg)
      real, intent(out)  :: vw_sfc(nxbg,nybg)
      real, intent(out)  :: sh_sfc(nxbg,nybg)
      real, intent(out)  :: tp_sfc(nxbg,nybg)
      real, intent(out)  :: ht_sfc(nxbg,nybg)
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
c
c     real ::  prbg_ht(100)
c     real ::  prbg_sh(100)
c     real ::  prbg_uv(100)
c     real ::  prbg_ww(100)

      real, allocatable ::  prbg_ht(:)
      real, allocatable ::  prbg_sh(:)
      real, allocatable ::  prbg_uv(:)
      real, allocatable ::  prbg_ww(:)


c in the near future this array will be used to read
c the sbn grids so that we avoid some inconsistencies
c with array allocations in the vertical interpolation.

      real, allocatable ::  data(:,:,:)

      integer start(10),count(10)
 
      integer i,j,k,n,ip,jp,ii,jj,it,kk
      integer istatus,slen,lent
      integer kskp
      integer ibdht,ibdtp,ibduv,ibdsh,ibdww
c
      character*9   fname,oldfname,model
      character*5   ctype
      character*4   af
      character*4   cvar
      character*2   gproj
      character*200 cdfname
      character*132 cmodel
c
      real   xe,mrsat
      real   make_ssh
c
c *** common block variables for lambert-conformal grid.
c
c     integer nx_lc,ny_lc,nz_lc  !no. of lc domain grid points
c     real   lat1,lat2,lon0,       !lambert-conformal std lat1, lat, lon
c    .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
c     common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
c     real   lon0_lc
c     real   lat1_lc,lat2_lc

      integer nf_vid,nn
      real cp,rcp, factor
      parameter (cp=1004.,rcp=287./cp)
c
ccc      save htn,tpn,rhn,uwn,vwn,prn,oldfname
c_______________________________________________________________________________
c
      interface
        subroutine get_prbg(nf_fid,nlvls,cvar,cmodel
     +,pr_levels_bg)
          integer        nf_fid
          integer        nlvls
          character*4    cvar
          character*132  cmodel
          real  ::       pr_levels_bg(nlvls)
        end subroutine

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

      print*,'here: read_sbn_grids'

      istatus = 1

      call s_len(cdfname,slen)

      print*,'cdfname: ',cdfname(1:slen)

      call get_nvaltimes(cdfname,ntbg,ivaltimes,istatus)
c
c *** open the netcdf file.
c
      print*,'opening cdf file: ',cdfname(1:slen)

      rcode = nf_open(cdfname,nf_nowrite,ncid)
      if(rcode.ne.nf_noerr) then
         print *, nf_strerror(rcode)
         print *,'nf_open ',cdfname(1:slen)
         return
      endif

      read(af,'(i4)') nn

      rcode=nf_inq_varid(ncid,'valtimeminusreftime',nf_vid)
      if(rcode.ne.nf_noerr) then
         print *, nf_strerror(rcode)
         print *,'in nf_get_var: ',cmodel
         return
      endif
      rcode=nf_get_vara_int(ncid,nf_vid,1,ntbg,ivaltimes)
      if(rcode.ne.nf_noerr) then
         print *, nf_strerror(rcode)
         print *,'in nf_get_var: ',cmodel
         return
      endif

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

c     if(.not.allocated(prbg_ht))allocate(prbg_ht(mxlvls))
c     if(.not.allocated(prbg_sh))allocate(prbg_sh(mxlvls))
c     if(.not.allocated(prbg_uv))allocate(prbg_uv(mxlvls))
c     if(.not.allocated(prbg_ww))allocate(prbg_ww(mxlvls))

      print*,'allocating data array for i/o'

      if(.not.allocated(data))allocate (data(nxbg,nybg,60)) 
c
      kskp=1
      if(cmodel.eq.'avn_sbn_cyleq')kskp=2
      if(cmodel(1:7).eq.'mesoeta')kskp=0


      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbght
      start(4)=n
      count(4)=1

      print*,'read ht ',ncid,nxbg*nybg*count(3),start,count,rcode
      print*,'nxbg/nybg/nzbght ',nxbg,nybg,nzbght
      cvar='gh'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),data,start
     +     ,count,rcode)
      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (gh): ',cmodel
         else
            print *,'missing ht data detected in read_sbn_grids: return'
         endif
         print*
         return
      endif

      do k=1,nzbght
         ht(:,:,k)=data(:,:,k)
      enddo

      ht_sfc = ht(:,:,1)
c
c get the pressures for this variable
c
      cvar='gh'
      allocate(prbg_ht(nzbght))
      call get_prbg(ncid,nzbght,cvar,cmodel,prbg_ht)
c
c ****** statements to fill tp.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgtp+kskp
      start(4)=n
      count(4)=1

      print*,'read tp'
      cvar='t'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),data,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (t): ',cmodel
         else
            print *,'missing t data detected in read_sbn_grids: return'
         endif
         print*
         return
      endif

      tp_sfc = data(:,:,1)
      k=0
      do kk=kskp+1,nzbgtp+kskp
         k=k+1
         tp(:,:,k)=data(:,:,kk) 
      enddo

c
c ****** statements to fill rh.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgsh+kskp
      start(4)=n
      count(4)=1

      print*,'read rh'
      cvar='rh'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),data,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (rh): ',cmodel
         else
            print *,'missing rh data detected in read_sbn_grids: return'
         endif
         print*
         return
      endif
      sh_sfc = data(:,:,1)
      k=0
      do kk=kskp+1,nzbgsh+kskp
         k=k+1
         sh(:,:,k)=data(:,:,kk)
      enddo
c
c get the pressures for this variable
c
      allocate (prbg_sh(nzbgsh))
      call get_prbg(ncid,nzbgsh,cvar,cmodel,prbg_sh)
c
c ****** statements to fill uw. 
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv+kskp
      start(4)=n
      count(4)=1

      print*,'read uw'
      cvar='uw'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),data,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (uw): ',cmodel
         else
            print *,'missing u data detected in read_sbn_grids: return'
         endif
         print*
         return
      endif

      uw_sfc = data(:,:,1)
      k=0
      do kk=kskp+1,nzbguv+kskp
         k=k+1
         uw(:,:,k)=data(:,:,kk)
      enddo

c
c get the pressures for this variable
c
      allocate (prbg_uv(nzbguv))
      call get_prbg(ncid,nzbguv,cvar,cmodel,prbg_uv)

c
c ****** statements to fill vw.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv+kskp
      start(4)=n
      count(4)=1

      print*,'read vw'
      cvar='vw'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),data,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (vw): ',cmodel
         else
            print *,'missing v data detected in read_sbn_grids: return'
         endif
         print*
         return
      endif
      vw_sfc = data(:,:,1)
      k=0
      do kk=kskp+1,nzbguv+kskp
         k=k+1
         vw(:,:,k)=data(:,:,kk)
      enddo
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
      cvar='pvv'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),data
     +  ,start,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (ww): ',cmodel
         else
            print *,'missing ww data detected: continue without'
            print *,'filling ww with 0.0'
            data(:,:,1:nzbgww) = 0.0
         endif
         print*
      endif
      do k=1,nzbgww
         ww(:,:,k)=data(:,:,k)
      enddo
c
c get the pressures for this variable
c
      allocate (prbg_ww(nzbgww))
      call get_prbg(ncid,nzbgww,cvar,cmodel,prbg_ww)

      deallocate (data)

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
      
      cvar='p'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,pr_sfc,start
     +     ,count,rcode)

      if(rcode.ne.nf_noerr) then
         if(rcode.gt.-61)then
            print *, nf_strerror(rcode)
            print *,'in nf_get_var (p): ',cmodel
         endif
         print *,'missing sfc p data detected in read_sbn_grids'
         print*,' -> continue without; compute in sfcbkgd'
         lcmpsfcq=.false.
         print*
c        return
      endif
c
c get mslp (this field name differs from one model to the other)
c
      if(cmodel.eq.'eta48_conus'.or.
     .   cmodel(1:7).eq.'mesoeta')then

         cvar='emsp'
         call read_netcdf_real(ncid,cvar,nxbg*nybg,mslp
     +           ,start,count,rcode)
         if(rcode.ne.nf_noerr) then
            if(rcode.gt.-61)then
               print *, nf_strerror(rcode)
               print *,'in nf_get_var (emsp): ',cmodel
            else
               print*,'error status returned from read_netcdf_real'
            endif
            print *,'missing emsp data detected in read_sbn_grids'
            print*
            return
         endif

      elseif(cmodel.eq.'ruc40_native')then

         cvar='mmsp'
         call read_netcdf_real(ncid,cvar,nxbg*nybg,mslp
     +           ,start,count,rcode)
         if(rcode.ne.nf_noerr) then
            if(rcode.gt.-61)then
               print *, nf_strerror(rcode)
               print *,'in nf_get_var (mmsp): ',cmodel
            else
               print*,'error status returned from read_netcdf_real'
            endif
            print *,'missing mmsp data detected in read_sbn_grids'
            print*
            return
         endif

      elseif(cmodel.eq.'avn_sbn_cyleq')then

         cvar='pmsl'
         call read_netcdf_real(ncid,cvar,nxbg*nybg,mslp
     +           ,start,count,rcode)
         if(rcode.ne.nf_noerr) then
            if(rcode.gt.-61)then
               print*, nf_strerror(rcode)
               print*,'in nf_get_var (pmsl): ',cmodel
            else
               print*,'error status returned from read_netcdf_real'
            endif
            print*,'missing pmsl data detected in read_sbn_grids'
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
c *** convert rh to sh.
c
      print*,'load prbg arrays'
      do j=1,nybg
      do i=1,nxbg
         prbgsh(i,j,:)=prbg_sh(:)
         prbght(i,j,:)=prbg_ht(:)
         prbguv(i,j,:)=prbg_uv(:)
         prbgww(i,j,:)=prbg_ww(:)
      enddo
      enddo

      deallocate (prbg_ht,prbg_sh,prbg_uv,prbg_ww)

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
         do k=1,nzbgsh
            do kk=1,nzbgtp
               if(prbght(1,1,kk).eq.prbgsh(1,1,k))then
                  goto500
               endif
            enddo
            print*,'------------------------------------'
            print*,'error: no match prbght and prbgsh'
            print*,'       in read_sbn_grids.'
            print*,'return with no data.'
            print*,'------------------------------------'
            goto 999 
500         do i=1,nxbg
            do j=1,nybg
               if (sh(i,j,k).lt.200)then
                sh(i,j,k)= make_ssh(prbgsh(i,j,k)
     .,tp(i,j,kk)-273.15,sh(i,j,k)/100., -132.)*0.001  !kg/kg 

c                  it=tp(i,j,k)*100
c                  it=min(45000,max(15000,it))
c                  xe=esat(it)
c                  mrsat=0.00622*xe/(prbgsh(i,j,k)-xe)
c                  sh(i,j,k)=sh(i,j,k)*mrsat
c                  sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))

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

      istatus = 0

      if(0.eq.1) then
 900     print*,'error: bad dimension specified in netcdf file'
         print*, (count(i),i=1,4)
         istatus=-1
      endif
 999  return
      end

c -----------------------------------------------------------

      subroutine get_prbg(nf_fid,nlvls,cvar,cmodel
     +,prlevels_bg)

      implicit none

      integer nf_fid
      integer nf_vid
      integer nf_status
      integer nlvls
      integer level
      integer lenc,lc,nclen
      integer k,kk
     
      character      cvar*4
      character      cmodel*132
      character      clvln10(100)*10
      character      clvln11(100)*11
      character      cnewvar*10
      character      ctmp*10
 
      real prlevels_bg(nlvls)
      real pr_levels_bg(100)

      include 'netcdf.inc'

      call s_len(cvar,lenc)
      cnewvar=cvar(1:lenc)//'levels'
      nf_status = nf_inq_varid(nf_fid,cnewvar,nf_vid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in var id: ',cnewvar
         return
      endif

      clvln11 = '           '
      clvln10 = '          '

      call s_len(cmodel,nclen)
      if(cmodel(1:nclen).eq.'avn_sbn_cyleq')then
         nf_status = nf_get_var_text(nf_fid,nf_vid,clvln11)
      else
         nf_status = nf_get_var_text(nf_fid,nf_vid,clvln10)
         if(nf_status.eq.nf_noerr)then
            do k=1,100
               clvln11(k)=clvln10(k)
            enddo
         endif
      endif
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in var: ',clvln11
         return
      endif
      pr_levels_bg = 0.0
      kk=0
      do k=1,100
         call s_len(clvln11(k),lc)
         if(clvln11(k)(1:2).eq.'mb')then
            kk=kk+1
            ctmp=trim(clvln11(k)(4:10))
            read(ctmp,'(i4.4)')level
            pr_levels_bg(kk)=float(level)
         endif
      enddo
      do k=1,kk
         prlevels_bg(k)=pr_levels_bg(k)
      enddo

      return
      end
