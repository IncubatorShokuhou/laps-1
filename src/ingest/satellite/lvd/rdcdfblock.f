      subroutine rdblock_line_elem(csat_id,csat_type,chtype,
     &ncid,varid,n_elems,n_lines,data,istatus)
c
c routine designed to read 3-d satellite sounding and image netcdf data with variable
c line and element dimensions.
c
      implicit none

      include 'netcdf.inc'
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      integer      rcode
c
      integer n_elems,n_lines,nch
      real    data(n_elems, n_lines)
      integer*2 data_int2(n_elems, n_lines)

      integer n,nn
      integer varid,ncid
      integer istatus
      integer ndsize_ch
      integer ndsize_x
      integer ndsize_y
      integer dim_id_x
      integer dim_id_y
      integer dim_id_k
      integer istart
      integer iend
      integer jstart
      integer jend

      integer start(10)
      integer count(10)
      character*31 dummy
      character*40 c_data_type
      character*6  csat_id
      character*3  csat_type
      character*3  chtype
      integer i,j,ii,jj 
c
c **************************************************************************
c
      istatus = -1
c
c get starting and ending line/pixel values for the satellite data in question
c
      nn=index(chtype,' ')-1
      if(nn.le.0)nn=3

      call getsat_attributes(csat_id,csat_type,chtype,
     &istart,iend,jstart,jend,ndsize_x,ndsize_y,istatus)
      if(istatus.ne.1)goto 1000
      istatus = -1
c
c qc step
c -------
      if(jend.gt.ndsize_y)then
         write(6,*)'jend > ndsize_y - rdcdfblock.f ',jend,ndsize_y
         write(6,*)'returing: istatus = ',istatus
         return
      endif

      if(iend.gt.ndsize_x)then
         write(6,*)'iend > ndsize_x ',iend,ndsize_x
         write(6,*)'returning to main, istatus = ',istatus
         return
      endif

      start(1)=istart
      count(1)=iend-istart+1
      start(2)=jstart
      count(2)=jend-jstart+1
      if(csat_type.eq.'cdf'.or.
     +   csat_type.eq.'wfo')then
         start(3)=1
         count(3)=1
      elseif(csat_type.eq.'ncp')then
         start(3)=1
         count(3)=1
         start(4)=1
         count(4)=1
      endif
c
c read line. switch here discriminates 1-byte versus 2-byte data.
c
      if(csat_type.eq.'rll' .or. csat_type.eq.'gnp')then
          write(6,*)'ncid/varid = ',ncid,varid
          write(6,*)'start=',start(1:2),' count=',count(1:2)
          rcode=nf_get_vara_int2(ncid,varid,start,count,data_int2)
          if (rcode .ne. nf_noerr) then
              print *, nf_strerror(rcode)
          endif
          write(6,*)'center pixel i2: ',data_int2(n_elems/2,n_lines/2)
          write(6,*)'rdblock_line_elem i2 data range: '
     1              ,minval(data_int2),maxval(data_int2)
          data(:,:) = data_int2(:,:)
      else
          rcode=nf_get_vara_real(ncid,varid,start,count,data)
          if (rcode .ne. nf_noerr) then
              print *, nf_strerror(rcode)
          endif
      endif

      write(6,*)'rdblock_line_elem data range: '
     1          ,minval(data),maxval(data)

      if(csat_type.ne.'rll')then
          do j=1,n_lines
          do i=1,n_elems
              if(data(i,j).lt.0)data(i,j)=data(i,j)+256.
          enddo
          enddo
      endif

      write(6,*)'center pixel r4: ',data(n_elems/2,n_lines/2)
c
      istatus = 1
1000  return
      end
