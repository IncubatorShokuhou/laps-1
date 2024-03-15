      subroutine read_nowrad_cdf(ctype,cfname_in,lines,elems,
     + dlat,dlon,la1,lo1,la2,lo2,centerlon,toplat,validtime,
     + dx,dy,lov, latin, image, istatus)

      include 'netcdf.inc'
      integer elems, lines,nf_fid, nf_vid, nf_status
      integer nx, ny, image( elems, lines), validtime
      real dx, dy, la1, la2, lo1, lo2, centerlon, difflon,
     +     dlat, dlon, toplat, lov, latin

      character*200 x_dim
      character*200 y_dim
      character*200 grid_name
      character*200 grid_type
      character*(*)  cfname_in
      character*(*)  ctype
c
      call read_nowrad_head(ctype,cfname_in,nx,ny,validtime,
     + dx, dy, la1, la2, lo1, lo2, centerlon, difflon,dlat,dlon,
     + toplat, grid_name, grid_type,x_dim,y_dim,lov,latin)

      call read_wsi_image(cfname_in,elems,lines,image)

      return
      end
c
c  ===============================================================
c  subroutine to read the file "wsi nowrad data header" 
c  ===============================================================
c
      subroutine read_nowrad_head(ctype,cfname_in,nx,ny, 
     +     validtime, dx, dy, la1, la2, lo1, lo2, centerlon, difflon, 
     +     radsperelem, radsperline, toplat, grid_name, grid_type, 
     +     x_dim, y_dim, lov, latin)
c
      include 'netcdf.inc'
      integer elems, lines,nf_fid, nf_vid, nf_status
      integer nx, ny, validtime
      real dx, dy, la1, la2, lo1, lo2, centerlon, difflon,
     +     radsperelem, radsperline, toplat, lov, latin

      real*8 valtime
      integer attlen
      character c_atvalue*80
      character*(*) x_dim
      character*(*) y_dim
      character*(*) grid_name
      character*(*) grid_type
      character*(*) cfname_in
      character*(*) ctype
c
c  open netcdf file for reading
c
      call s_len(cfname_in,nc)
      nf_status = nf_open(cfname_in,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', cfname_in(1:nc)
        return
      endif

c   variables of type real
c
c     variable        netcdf long name
c      dx           
c
        nf_status = nf_inq_varid(nf_fid,'dx',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dx'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,dx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dx'
      endif

      call ncainq(nf_fid,nf_vid,'units',itype,attlen,istatus)
      call ncagtc(nf_fid,nf_vid,'units',c_atvalue,attlen,istatus)
      if(istatus.ne.0)then
         write(6,*)'error getting attribute - dx'
      endif

c
c     variable        netcdf long name
c      dy           
c
        nf_status = nf_inq_varid(nf_fid,'dy',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dy'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,dy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dy'
      endif
c
c wsi-type is degrees. wfo-type may be kilometers. dx and dy must = m.
c
      if(c_atvalue(1:4).eq.'kilo')then
         dx=dx*1000.
         dy=dy*1000.
      elseif(c_atvalue(1:4).eq.'degr')then
         dx=dx*111.1*1000.
         dy=dy*111.1*1000.
      endif
c
c     variable        netcdf long name
c      la1          
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
c      la2          
c
        nf_status = nf_inq_varid(nf_fid,'la2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la2'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,la2)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var la2'
        endif
      endif
c
c     variable        netcdf long name
c      lo1          
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
c
c     variable        netcdf long name
c      lo2          
c
        nf_status = nf_inq_varid(nf_fid,'lo2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo2'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,lo2)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lo2'
        endif
      endif

c ******** wsi type variables ***********

      if(ctype.eq.'wsi')then

c
c     variable        netcdf long name
c      centerlon    "center longitude"
c
          nf_status = nf_inq_varid(nf_fid,'centerlon',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var centerlon'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,centerlon)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var centerlon'
        endif
c
c     variable        netcdf long name
c      difflon      "difference longitude"
c
          nf_status = nf_inq_varid(nf_fid,'difflon',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var difflon'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,difflon)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var difflon'
        endif
c
c     variable        netcdf long name
c      radsperelem  "radians per element"
c
          nf_status = nf_inq_varid(nf_fid,'radsperelem',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var radsperelem'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,radsperelem)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var radsperelem'
        endif
c
c     variable        netcdf long name
c      radsperline  "radians per line"
c
          nf_status = nf_inq_varid(nf_fid,'radsperline',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var radsperline'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,radsperline)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var radsperline'
        endif
c
c     variable        netcdf long name
c      toplat       "top latitude"
c
          nf_status = nf_inq_varid(nf_fid,'toplat',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var toplat'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,toplat)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var toplat'
        endif

c ------------------------------

      elseif(ctype.eq.'wfo')then

c
c     variable        netcdf long name
c      latin       "latitude of tangent"
c
          nf_status = nf_inq_varid(nf_fid,'latin',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var latin'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,latin)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var latin'
        endif
c
c     variable        netcdf long name
c      lov       "longitude of vertical"
c
          nf_status = nf_inq_varid(nf_fid,'lov',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lov'
        endif
          nf_status = nf_get_var_real(nf_fid,nf_vid,lov)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var lov'
        endif

      else

        print*,'invalid radar type in ctype ',ctype

      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      nx           "number of x points"
c
        nf_status = nf_inq_varid(nf_fid,'nx',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nx'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,nx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nx'
      endif
c
c     variable        netcdf long name
c      ny           "number of y points"
c
        nf_status = nf_inq_varid(nf_fid,'ny',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ny'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,ny)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ny'
      endif

      if(ctype.eq.'wsi')then
c
c     variable        netcdf long name
c      validtime    "valid time"
c
          nf_status = nf_inq_varid(nf_fid,'validtime',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var validtime'
        endif
          nf_status = nf_get_var_int(nf_fid,nf_vid,validtime)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var validtime'
        endif

      elseif(ctype.eq.'wfo')then

c
c     variable        netcdf long name
c      validtime    "valid time"
c
        nf_status=nf_inq_varid(nf_fid,'valtime',nf_vid)
        if(nf_status.ne.0)then
           nf_status=nf_inq_varid(nf_fid,'validtime',nf_vid)
           if(nf_status.ne.0)then
              write(6,*)'error getting variable - valtime'
              return
           endif
        endif
        nf_status=nf_get_var1_double(nf_fid,nf_vid,1,valtime)
        validtime=int(valtime)

      endif

c   variables of type char
c
c
c     variable        netcdf long name
c      grid_name    "grid name"
c
        nf_status = nf_inq_varid(nf_fid,'grid_name',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_name'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,grid_name)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_name'
      endif
c
c     variable        netcdf long name
c      grid_type    "grid type"
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
c
c     variable        netcdf long name
c      x_dim        "longitude dimension"
c
        nf_status = nf_inq_varid(nf_fid,'x_dim',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var x_dim'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,x_dim)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var x_dim'
      endif
c
c     variable        netcdf long name
c      y_dim        "latitude dimension"
c
        nf_status = nf_inq_varid(nf_fid,'y_dim',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var y_dim'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,y_dim)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var y_dim'
      endif
c
c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
c
c ================================================================
c
      subroutine read_wsi_image(cfname_in,elems,lines,image)

      include 'netcdf.inc'

      integer elems, lines,nf_fid, nf_vid, nf_status
      integer image(elems,lines)
      character*(*) cfname_in
c
c  open netcdf file for reading
c
      call s_len(cfname_in,nc)
      nf_status = nf_open(cfname_in,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', cfname_in(1:nc)
        return
      endif
c
c     variable        netcdf long name
c      image        "image pixel values"
c
        nf_status = nf_inq_varid(nf_fid,'image',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var image'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,image)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var image'
      endif

c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
c
c ============================================================
c routine to convert netcdf data base (bytes) to nowrad scaled
c integers.
c ============================================================

      subroutine cvt_wsi_nowrad(ctype,elems,lines,image
     +,nsites_present,present_site_loc_i,present_site_loc_j
     +,nsites_absent,absent_site_loc_i,absent_site_loc_j
     +,max_sites,istatus)
c
c
      implicit none

      integer elems,lines
      integer image(elems,lines)

      integer max_sites

      integer nsites_present
      real    present_site_loc_i(max_sites)
      real    present_site_loc_j(max_sites)

      integer nsites_absent
      real    absent_site_loc_i(max_sites)
      real    absent_site_loc_j(max_sites)

      integer istatus
      integer imax_image_value,imin_image_value,i,j
      integer icount_out, icount_bad
      integer bad_data_flag
      parameter(bad_data_flag=255) 
 
      character ctype*(*)

      if(ctype.eq.'wfo')then
         do j=1,lines
         do i=1,elems
            if(image(i,j).lt.0)image(i,j)=256+image(i,j)
            image(i,j)=image(i,j)/16
         enddo
         enddo
      endif
 
      imax_image_value = 0
      imin_image_value = 255
      icount_bad=0
      icount_out=0
      do j=1,lines
         do i=1,elems
            if(image(i,j) .gt. imax_image_value)
     +           imax_image_value = image(i,j)
            if(image(i,j) .lt. imin_image_value)
     +           imin_image_value = image(i,j)
            if(image(i,j) .ne. bad_data_flag)then
               if((image(i,j).gt.64).or.(image(i,j).lt.0))then
                  icount_out = icount_out +1
                  image(i,j) = bad_data_flag
c     write(6,*) i, j, i_value
               endif

            else
               icount_bad=icount_bad+1
               image(i,j) = bad_data_flag
            endif
         enddo
      enddo
c
      write(6,*)'number of bad data points (> ',bad_data_flag,' )'
      write(6,*)'prior to calling c_scan_adjust: ',icount_bad
      write(6,*)'max value found in image array: ',imax_image_value
      write(6,*)'min value found in image array: ',imin_image_value
      write(6,*)
      write(6,*)'data found out-of-bounds (icount_out) ',icount_out
      if(icount_out.gt.0)then
         istatus = -1
         return
      endif


c       new c version of scan_adjust implemented 05-apr-96
ccc      call c_scan_adjust(image,lines,elems,bad_data_flag)

      nsites_present=0
      nsites_absent =0
      if(ctype.eq.'wsi')then
         do j=1,lines
         do i=1,elems
            if(image(i,j).ne.bad_data_flag)then
               if(image(i,j).gt.15.and.image(i,j).lt.32)then
                  nsites_present=nsites_present+1
                  present_site_loc_i(nsites_present)=float(i)
                  present_site_loc_j(nsites_present)=float(j)
               elseif(image(i,j).ge.32 .and.image(i,j).lt.48)then
                  nsites_absent=nsites_absent+1
                  absent_site_loc_i(nsites_absent)=float(i)
                  absent_site_loc_j(nsites_absent)=float(j)
               endif
               image(i,j)=mod(image(i,j),16)
            endif
         enddo
         enddo
      endif

 9999 return
      end
c
c ----------------------------------------------------------
c
      subroutine read_nowrad_dims(ctype,cfname,elems,lines)

      include 'netcdf.inc'

      integer elems,lines
      character*(*) cfname
      character*(*) ctype
      character*10  cvarname(2)

c
c  open netcdf file for reading
c
      call s_len(cfname,nc)
      nf_status = nf_open(cfname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', cfname(1:nc)
        return
      endif
      cvarname(1)='x'
      cvarname(2)='y'
      if(ctype.eq.'wsi')then
         cvarname(1)='elems'
         cvarname(2)='lines'
      endif
      call s_len(cvarname(1),il)
c
c get size of elems
c
      nf_status = nf_inq_dimid(nf_fid,cvarname(1)(1:il),nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim elems'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,elems)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim elems'
      endif
c
c get size of lines
c
      nf_status = nf_inq_dimid(nf_fid,cvarname(2)(1:il),nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim lines'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,lines)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim lines'
      endif
c
c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end

