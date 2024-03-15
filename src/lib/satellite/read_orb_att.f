c
c get size of num_att
c
c     nf_status = nf_inq_dimid(nf_fid,'num_att',nf_vid)
c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'dim num_att'
c     endif
c     nf_status = nf_inq_dimlen(nf_fid,nf_vid,num_att)
c     if(nf_status.ne.nf_noerr) then
c       print *, nf_strerror(nf_status)
c       print *,'dim num_att'
c     endif
c     call main_sub(nf_fid , num_att)

c     end
c
c
      subroutine read_orb_att(c_filespec,csatid, num_att, orb_att,
     &istatus)
      include 'netcdf.inc'
c     integer num_att, nf_fid, nf_vid, nf_status
      integer num_att, nf_fid, nf_status
      character*6 csatid
c     character*7 sat_name
      character*(*) c_filespec
      character*255 cfname
      double precision orb_att(num_att)

      istatus = 1
c
c  open netcdf file for reading
c
      cfname=c_filespec//'/'//csatid//'_orbatt.dat'
      n=index(cfname,' ')-1
      nf_status = nf_open(cfname,nf_nowrite,nf_fid)

      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',cfname(1:n)
        return
      endif

      call read_netcdf(nf_fid , num_att, orb_att,istatus)
c
c the netcdf variables are filled - your code goes here
c
      return
      end

      subroutine read_netcdf(nf_fid , num_att, orb_att,istatus)
      include 'netcdf.inc'
      integer num_att, nf_fid, nf_vid, nf_status

      character*7 sat_name
      double precision orb_att(num_att)

      istatus = 1
c
c     variable        netcdf long name
c      sat_name
      nf_status = nf_inq_varid(nf_fid,'sat_name',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sat_name'
        return
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,sat_name)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ sat_name '
        return
      endif
      write(6,*)'sat name in netcdf file ',sat_name
c
c     variable        netcdf long name
c      orb_att      "orbit attitudes" 
c
        nf_status = nf_inq_varid(nf_fid,'orb_att',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var orb_att'
        return
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,orb_att)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ orb_att '
        return
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 0
      return
      end
