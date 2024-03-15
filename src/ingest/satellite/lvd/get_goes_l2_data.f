      subroutine get_goes_l2_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_offset
     +                   ,maxchannel,max_files,nchannels
     +                   ,csat_id,csat_type,chtype
     +                   ,path
     +                   ,image_13,n_ir_elem,n_ir_lines
     +                   ,image_07
     +                   ,image_02,n_vis_elem,n_vis_lines
     +                   ,i4time_goes,istatus)                      ! o

      include 'netcdf.inc'

      character*(*) path
      character*255 filename

      integer nelem, nlines,nf_fid, nf_vid, nf_status
      real image_13(n_ir_elem,n_ir_lines)
      real image_07(n_ir_elem,n_ir_lines)
      real image_02(n_vis_elem,n_vis_lines)

      character*6  csat_id
      character*3  csat_type
      character*3  chtype(maxchannel),chtype_local
      character*9  a9time
      character*13 cfname13,cvt_i4time_wfo_fname13

      write(6,*)' subroutine get_goes_l2_data ',n_ir_elem,n_ir_lines

      i4time_goes = ((i4time_sys + 0)
     1            / ilaps_cycle_time) * ilaps_cycle_time + i4time_offset

!     read ir channel 13 data
      call make_fnam_lp(i4time_goes,a9time,istatus)
      filename = trim(path)//'/'//a9time//'_10p4.nc'
      write(6,*)'gr2/13 fname is ',trim(filename)
c
c  open netcdf file for reading
c
      nf_status=nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),filename
        istatus=0
        return
      endif
c
c  fill all dimension values
c
c
c get size of nelem
c
      nf_status=nf_inq_dimid(nf_fid,'x',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nelem'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,nelem)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nelem'
      endif
c
c get size of nlines
c
      nf_status=nf_inq_dimid(nf_fid,'y',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nlines'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,nlines)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nlines'
      endif

      chtype_local = '11u'
      call read_goes_np_data(nf_fid,n_ir_elem,n_ir_lines,
     +     csat_id,csat_type,chtype_local,
     +     ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, image_13, istatus)

      nf_status=nf_close(nf_fid)
      
      if(.true.)then

!     read ir channel 07 data
      call make_fnam_lp(i4time_goes,a9time,istatus)
      filename = trim(path)//'/'//a9time//'_4u.nc'
      write(6,*)'gr2/07 fname is ',trim(filename)
c
c  open netcdf file for reading
c
      nf_status=nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),filename
        istatus=0
!       return
      endif
c
c  fill all dimension values
c
c
c get size of nelem
c
      nf_status=nf_inq_dimid(nf_fid,'x',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nelem'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,nelem)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nelem'
      endif
c
c get size of nlines
c
      nf_status=nf_inq_dimid(nf_fid,'y',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nlines'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,nlines)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nlines'
      endif

      chtype_local = '39u'
      call read_goes_np_data(nf_fid,n_ir_elem,n_ir_lines,
     +     csat_id,csat_type,chtype_local,
     +     ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, image_07, istatus)

      nf_status=nf_close(nf_fid)

!     read vis channel 02 data
      call make_fnam_lp(i4time_goes,a9time,istatus)
      filename = trim(path)//'/'//a9time//'_vis.nc'
      write(6,*)'gr2/02 fname is ',trim(filename)
c
c  open netcdf file for reading
c
      nf_status=nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),filename
        istatus=0
        return
      endif
c
c  fill all dimension values
c
c
c get size of nelem
c
      nf_status=nf_inq_dimid(nf_fid,'x',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nelem'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,nelem)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nelem'
      endif
c
c get size of nlines
c
      nf_status=nf_inq_dimid(nf_fid,'y',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nlines'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,nlines)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim nlines'
      endif

      chtype_local = 'vis'
      call read_goes_np_data(nf_fid,n_vis_elem,n_vis_lines,
     +     csat_id,csat_type,chtype_local,
     +     ilaps_cycle_time, nx_l, ny_l, i4time_earliest,
     +     i4time_latest, lun_out, image_02, istatus)

      endif

      nf_status=nf_close(nf_fid)

      write(6,*)' call qc_goes_vis'

      call qc_goes_vis(image_02,n_vis_elem,n_vis_lines)

      istatus = 1

      return
      end
