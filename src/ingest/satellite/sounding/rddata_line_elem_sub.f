      subroutine rddata_line_elem(ncid,imax,jmax,nch
     &,sounding,istatus)
c
c routine designed to read 3-d satellite sounding netcdf data with variable
c line and element dimensions. should work for 2d data with variable dimensions
c if nch = 1
c
      implicit none
      include 'netcdf.inc'
      integer   rcode
c
      integer   imax,jmax,nch
      integer   sounding(imax, jmax, nch)
      integer   sounding_rec(imax)

      integer i,j,k
      integer varid,ncid
      integer istatus
      integer dim_id_x
      integer dim_id_y
      integer dim_id_k
      integer i4dum

      integer start(10)
      integer count(10)
      character*31 dummy
c
c **************************************************************************
c
      istatus = 1
c
c now ready to read the data
c
      start(1)=1
      count(2)=1
      count(3)=1
c
c get sounding variable id
c
      rcode=nf_inq_varid(ncid,'sounding',varid)
      if(rcode.ne.0)then
         write(6,*)'error getting sounding varid'
         istatus = -1
         return
      endif

      do k = 1,nch  !# of channels in sounding database (= 19).

         write(6,*)'reading channel ',k

         start(3)=k

         do j = 1,jmax

            count(1)=imax
            start(2)=j

c read line
      rcode=nf_get_vara_int(ncid,varid,start,count,sounding_rec)
            if(rcode.ne.0)then 
               print*,'error reading sounding database ',rcode
               istatus = -1
               print*,'returning to main, istatus = ',istatus
               return
            endif
c
c load line into the 3d output array
c
            sounding(:,j,k) = sounding_rec(:)

         enddo
      enddo
c -------------------------------------
c
      write(6,*)'sucessful reading sounding - rddata_line_elem '

      return
      end
