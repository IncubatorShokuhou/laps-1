      subroutine get_cloud_drift_madis 
     ~           (i4time_sys, i4_window, filename, istatus)     

      implicit none
      include 'netcdf.inc'
      integer, intent(in)  :: i4time_sys, i4_window
      character*(*),intent(in) ::  filename
      integer, intent(out) :: istatus
      integer, parameter  :: lun_cdw = 11

      integer :: i,nobs, ncid,nf_vid,nf_status
      real  r_missing_data
      real, allocatable :: oblat(:),oblon(:),
     ~                     pressure(:),winddir(:), windspd(:)

      integer,parameter::double=selected_real_kind(p=13,r=200) 
      real(kind=double), allocatable :: validtime(:)
      character(len=1),allocatable ::winddirdd(:),windspddd(:)
      integer, parameter :: unix2i4 = 315619200
      character(len=9)  :: a9timeobs
      integer :: obtime, i4dif 
      integer :: nkept,nreject
      call get_r_missing_data(r_missing_data,istatus)
      if ( istatus .ne. 1 )  then
         write (6,*) 'error getting r_missing_data'
         return
      endif
      print *, "opening ", trim(filename)
      ! open the netcdf file
      nf_status = nf_open(filename,nf_nowrite,ncid)
      if (nf_status .ne. nf_noerr) then
        print *, "problem opening madis netcdf file."
        print *, "filename:",trim(filename)
        istatus = 0
        return
      endif
      ! get the number of records
      nf_status = nf_inq_dimid(ncid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'problem getting recnum'
        print *, ' in get_cloud_drift_madis'
        istatus = 0
        return
      endif
      nf_status = nf_inq_dimlen(ncid,nf_vid,nobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim record'
        istatus = 0
        return
      endif
c
      print *, " madis cloud drift: nobs = ", nobs

      print *, " allocating arrays..."
      allocate(oblat(nobs))
      allocate(oblon(nobs))
      allocate(validtime(nobs))
      allocate(pressure(nobs))
      allocate(winddir(nobs))
      allocate(winddirdd(nobs))
      allocate(windspd(nobs))
      allocate(windspddd(nobs))
 
      print *, "reading data" 
      ! read in the data
      nf_status = nf_inq_varid(ncid,"oblat",nf_vid)
      if (nf_status .ne. nf_noerr) then
        print *," problem getting oblat"
        istatus = 0
        goto 900
      endif
      nf_status = nf_get_var_real(ncid,nf_vid,oblat)

      nf_status = nf_inq_varid(ncid,"oblon",nf_vid)
      if (nf_status .ne. nf_noerr) then
        print *," problem getting oblon"
        istatus = 0
        goto 900
      endif
      nf_status = nf_get_var_real(ncid,nf_vid,oblon)

      nf_status = nf_inq_varid(ncid,"validtime",nf_vid)
      if (nf_status .ne. nf_noerr) then
        print *," problem getting validtime"
        istatus = 0
        goto 900
      endif
      nf_status = nf_get_var_double(ncid,nf_vid,validtime)

      nf_status = nf_inq_varid(ncid,"pressure",nf_vid)
      if (nf_status .ne. nf_noerr) then
        print *," problem getting pressure"
        istatus = 0
        goto 900
      endif
      nf_status = nf_get_var_real(ncid,nf_vid,pressure)

      nf_status = nf_inq_varid(ncid,"winddir",nf_vid)
      if (nf_status .ne. nf_noerr) then
        print *," problem getting winddir"
        istatus = 0
        goto 900
      endif
      nf_status = nf_get_var_real(ncid,nf_vid,winddir)

      nf_status = nf_inq_varid(ncid,"windspd",nf_vid)
      if (nf_status .ne. nf_noerr) then
        print *," problem getting windspd"
        istatus = 0
        goto 900
      endif
      nf_status = nf_get_var_real(ncid,nf_vid,windspd)

      nf_status = nf_inq_varid(ncid,"winddirdd",nf_vid)
      if (nf_status .ne. nf_noerr) then
        print *," problem getting winddirdd"
        istatus = 0
        goto 900
      endif
      nf_status = nf_get_var_text(ncid,nf_vid,winddirdd)

      nf_status = nf_inq_varid(ncid,"windspddd",nf_vid)
      if (nf_status .ne. nf_noerr) then
        print *," problem getting windspddd"
        istatus = 0
        goto 900
      endif
      nf_status = nf_get_var_text(ncid,nf_vid,windspddd)

 
      call open_ext(lun_cdw,i4time_sys,'cdw',istatus)
      nreject =0
      nkept  = 0
      do i= 1,nobs
c       if this ob passes the quality flag
         obtime = nint(validtime(i))+unix2i4
         i4dif = i4time_sys - obtime
         if (abs(i4dif) .le. i4_window) then
           if (windspddd(i) .ne. "c" .or. 
     ~         windspd(i) .lt. 0. .or.
     ~         windspd(i) .gt. 125.) then
             windspd(i) = r_missing_data
           endif
           if (winddirdd(i) .ne. "c" .or.
     ~         winddir(i) .lt. 0. .or.
     ~         winddir(i) .gt. 360.) then
             winddir(i) = r_missing_data
           endif
           call c_time2fname(nint(validtime(i)),a9timeobs)

           write (lun_cdw,'(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)') 
     ~                 oblat(i), oblon(i), pressure(i), 
     ~                 winddir(i), windspd(i), a9timeobs

           nkept = nkept+1
         else
           nreject = nreject +1
         endif
      enddo
      print *, "total obs/# kept/#rejected:",nobs,nkept,nreject
      istatus = 1
900   deallocate(oblat)
      deallocate(oblon)
      deallocate(validtime)
      deallocate(pressure)
      deallocate(winddir)
      deallocate(winddirdd)
      deallocate(windspd)
      deallocate(windspddd)
      return
      end
