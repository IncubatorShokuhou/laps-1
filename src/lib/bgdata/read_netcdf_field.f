      subroutine read_netcdf_real(nf_fid,fname,n1,f,start,count,istatus)

      implicit none

      include 'netcdf.inc'
      include 'bgdata.inc'
      integer n1,i, nf_fid, nf_vid,istatus,nf_status
      integer start(10),count(10)
      real f(n1) , nfmissing
      character*(*) fname

      print*,'here: in read_netcdf_real a, n1=',n1

      istatus=0
      nf_status = nf_inq_varid(nf_fid,fname,nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ', fname
        istatus = 1
        return
      endif

      if(start(1).eq.0.and.count(1).eq.0) then
         print*,'here: in read_netcdf_real b1'
         nf_status = nf_get_var_real(nf_fid,nf_vid,f)
      else
         print*,'here: in read_netcdf_real b2',start,count
         nf_status = nf_get_vara_real(nf_fid,nf_vid,start,count,f)
      endif
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ ', fname
        istatus = 1
        return
      endif

      print*,'here: in read_netcdf_real c'

      if(fname.ne.'isolevel')then
         nf_status = nf_get_att_real(nf_fid,nf_vid,'_fillvalue'
     .,nfmissing)
         if(nf_status.ne.nf_noerr) then
            print *, nf_strerror(nf_status)
         endif
      endif
      do i=1,n1
         if(f(i).eq.nfmissing) then
            f(i)=missingflag
            istatus=istatus-1
         endif
      enddo
      
      return
      end
