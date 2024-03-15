c-----------------------------------------------------------------------
      subroutine readnet_nam(ncid,imx,imy,imz,iblev,
     +                    gh_sfc,gh,rh_2mfh,rh,rh_lbis,t_2mfh,
     +                    t,t_lbis,uw_10mfh,uw,uw_lbis,vw,
     +                    vw_lbis,vw_10mfh,av,pvv,p_sfc,heli,
     +                    cape_sfc,cape_lbis,cin_sfc,cin,
     +                    bli_lbis,pli_lbis,pw,emspmsl,prmsl,
     +                    cp_sfc,tp_sfc,isolevel,boundrylevel,
     +                    valtime,reftime,origin,model,grid_type,
     +                    j_dim,i_dim,ni,nj,la1,la2,lo1,lo2,
     +                    di,dj,intlat1,intlat2,lon0,start,
     +                    count,vdims,strbuf,num_of_ens)
c-----------------------------------------------------------------------

      include 'netcdf.inc'

c     define variables.
c     variable ids run sequentially from 1 to nvars =   50

      parameter      (nvars=50)          ! number of variables
      parameter      (nrec=1)           ! change this to generalize
      parameter      (ndims=7)          ! number of dimensions
      parameter      (record=1)
      parameter      (namelen=132)
      parameter      (nav=1)  

      integer      rcode                     ! error code
      integer      recdim                    ! record dimension
      integer        imx,imy,imz,iblev
      integer        num_of_ens   

      real         gh_sfc(imx,imy,nrec)
      real         gh(imx,imy,imz,nrec)
      real         rh_2mfh(imx,imy,nrec)
      real         rh(imx,imy,imz,nrec)
      real         rh_lbis(imx,imy,iblev,nrec)
      real         t_2mfh(imx,imy,nrec)
      real         t(imx,imy,imz,nrec)
      real         t_lbis(imx,imy,iblev,nrec)
      real         uw_10mfh(imx,imy,nrec)
      real         uw(imx,imy,imz,nrec)
      real         uw_lbis(imx,imy,iblev,nrec)
      real         vw(imx,imy,imz,nrec)
      real         vw_lbis(imx,imy,iblev,nrec)
      real         vw_10mfh(imx,imy,nrec)
      real         av(imx,imy,imz,nrec)
      real         pvv(imx,imy,imz,nrec)
      real         p_sfc(imx,imy,nrec)
      real         heli(imx,imy,nrec)
      real         cape_sfc(imx,imy,nrec)
      real         cape_lbis(imx,imy,nrec)
      real         cin_sfc(imx,imy,nrec)
      real         cin(imx,imy,nrec)
      real         bli_lbis(imx,imy,nrec)
      real         pli_lbis(imx,imy,nrec)
      real         pw(imx,imy,nrec)
      real         emspmsl(imx,imy,nrec)
      real         prmsl(imx,imy,nrec)
      real         cp_sfc(imx,imy,nrec)
      real         tp_sfc(imx,imy,nrec)

      integer      isolevel(imz)
      integer      boundrylevel(iblev)

      real*8         valtime(nrec)
      real*8         reftime(nrec)
      character*1    origin(namelen) 
      character*1    model(namelen) 

      character*1    grid_type(namelen,nrec)
      character*1    j_dim(namelen,nrec)
      character*1    i_dim(namelen,nrec)

      integer*2      version
      integer*2      ni(nav) 
      integer*2      nj(nav) 

      real         la1(nav) 
      real         la2(nav) 
      real         lo1(nav) 
      real         lo2(nav)
      real         di(nav) 
      real         dj(nav) 
      real         intlat1(nav) 
      real         intlat2(nav) 
      real         lon0(nav) 

      integer      start(ndims)            ! hyperslab starting index
      integer      count(ndims)            ! hyperslab count from start
      integer        vdims(ndims)            ! max # of var dims
      character*1024 strbuf                  ! string buffer for var
                                             !  and attr names



c     get info on the record dimension for this file.
      call ncinq(ncid,ndims,nvars,ngatts,recdim,rcode)
      call ncdinq(ncid,recdim,strbuf,nrecs,rcode)
c     nrecs now contains the # of records for this file

c     retrieve data for gh_sfc variable.
      call ncvinq(ncid,   1,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   1,start,count,gh_sfc,rcode)

c     retrieve data for gh variable.
      call ncvinq(ncid,   2,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   2,start,count,gh,rcode)

c     retrieve data for rh_2mfh variable.
      call ncvinq(ncid,   3,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   3,start,count,rh_2mfh,rcode)

c     retrieve data for rh variable.
      call ncvinq(ncid,   4,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   4,start,count,rh,rcode)

c     retrieve data for rh_lbis variable.
      call ncvinq(ncid,   5,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   5,start,count,rh_lbis,rcode)

c     retrieve data for t_2mfh variable.
      call ncvinq(ncid,   6,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   6,start,count,t_2mfh,rcode)

c     retrieve data for t variable.
      call ncvinq(ncid,   7,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   7,start,count,t,rcode)

c     retrieve data for t_lbis variable.
      call ncvinq(ncid,   8,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   8,start,count,t_lbis,rcode)

c     retrieve data for uw_10mfh variable.
      call ncvinq(ncid,   9,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,   9,start,count,uw_10mfh,rcode)

c     retrieve data for uw variable.
      call ncvinq(ncid,  10,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  10,start,count,uw,rcode)

c     retrieve data for uw_lbis variable.
      call ncvinq(ncid,  11,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  11,start,count,uw_lbis,rcode)

c     retrieve data for vw variable.
      call ncvinq(ncid,  12,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  12,start,count,vw,rcode)

c     retrieve data for vw_lbis variable.
      call ncvinq(ncid,  13,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  13,start,count,vw_lbis,rcode)

c     retrieve data for vw_10mfh variable.
      call ncvinq(ncid,  14,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  14,start,count,vw_10mfh,rcode)

c     retrieve data for av variable.
      call ncvinq(ncid,  15,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  15,start,count,av,rcode)

c     retrieve data for pvv variable.
      call ncvinq(ncid,  16,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  16,start,count,pvv,rcode)

c     retrieve data for p_sfc variable.
      call ncvinq(ncid,  17,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  17,start,count,p_sfc,rcode)

c     retrieve data for heli variable.
      call ncvinq(ncid,  18,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  18,start,count,heli,rcode)

c     retrieve data for cape_sfc variable.
      call ncvinq(ncid,  19,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  19,start,count,cape_sfc,rcode)

c     retrieve data for cape_lbis variable.
      call ncvinq(ncid,  20,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  20,start,count,cape_lbis,rcode)

c     retrieve data for cin_sfc variable.
      call ncvinq(ncid,  21,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  21,start,count,cin_sfc,rcode)

c     retrieve data for cin variable.
      call ncvinq(ncid,  22,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  22,start,count,cin,rcode)

c     retrieve data for bli_lbis variable.
      call ncvinq(ncid,  23,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  23,start,count,bli_lbis,rcode)

c     retrieve data for pli_lbis variable.
      call ncvinq(ncid,  24,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  24,start,count,pli_lbis,rcode)

c     retrieve data for pw variable.
      call ncvinq(ncid,  25,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  25,start,count,pw,rcode)

c     retrieve data for emspmsl variable.
      call ncvinq(ncid,  26,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  26,start,count,emspmsl,rcode)

c     retrieve data for prmsl variable.
      call ncvinq(ncid,  27,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  27,start,count,prmsl,rcode)

c     retrieve data for cp_sfc variable.
      call ncvinq(ncid,  28,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  28,start,count,cp_sfc,rcode)

c     retrieve data for tp_sfc variable.
      call ncvinq(ncid,  29,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  29,start,count,tp_sfc,rcode)

c     retrieve data for isolevel variable.
      call ncvinq(ncid,  30,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  30,start,count,isolevel,rcode)

c     retrieve data for boundrylevel variable.
      call ncvinq(ncid,  31,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  31,start,count,boundrylevel,rcode)

c     retrieve data for valtime variable.
      call ncvinq(ncid,  32,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  32,start,count,valtime,rcode)

c     retrieve data for reftime variable.
      call ncvinq(ncid,  33,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  33,start,count,reftime,rcode)

c     retrieve data for origin variable.
      call ncvinq(ncid,  34,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  34,start,count,origin,lenstr,rcode)

c     retrieve data for model variable.
      call ncvinq(ncid,  35,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  35,start,count,model,lenstr,rcode)

c     retrieve data for version variable.
      call ncvinq(ncid,  36,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  36,start,count,version,rcode)

c     retrieve data for grid_type variable.
      call ncvinq(ncid,  37,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  37,start,count,grid_type,lenstr,rcode)

c     retrieve data for j_dim variable.
      call ncvinq(ncid,  38,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  38,start,count,j_dim,lenstr,rcode)

c     retrieve data for i_dim variable.
      call ncvinq(ncid,  39,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgtc(ncid,  39,start,count,i_dim,lenstr,rcode)

c     retrieve data for ni variable.
      call ncvinq(ncid,  40,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  40,start,count,ni,rcode)

c     retrieve data for nj variable.
      call ncvinq(ncid,  41,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  41,start,count,nj,rcode)

c     retrieve data for la1 variable.
      call ncvinq(ncid,  42,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  42,start,count,la1,rcode)

c     retrieve data for la2 variable.
      call ncvinq(ncid,  43,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  43,start,count,la2,rcode)

c     retrieve data for lo1 variable.
      call ncvinq(ncid,  44,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  44,start,count,lo1,rcode)

c     retrieve data for lo2 variable.
      call ncvinq(ncid,  45,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  45,start,count,lo2,rcode)

c     retrieve data for di variable.
      call ncvinq(ncid,  46,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  46,start,count,di,rcode)

c     retrieve data for dj variable.
      call ncvinq(ncid,  47,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  47,start,count,dj,rcode)

c     retrieve data for intlat1 variable.
      call ncvinq(ncid,  48,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  48,start,count,intlat1,rcode)

c     retrieve data for intlat2 variable.
      call ncvinq(ncid,  49,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  49,start,count,intlat2,rcode)

c     retrieve data for lon0 variable.
      call ncvinq(ncid,  50,strbuf,nctype,nvdim,vdims,nvatts,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),strbuf,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      end do
      call ncvgt(ncid,  50,start,count,lon0,rcode)

c
c     begin writing statements to use the data.
c

      return
      end
