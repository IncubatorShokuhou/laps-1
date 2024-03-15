module wrfutil

   implicit none

   save

   integer :: ncid

end module

!===============================================================================

subroutine get_wrf_dims(fname, nx, ny, nz, istatus)

   use wrfutil

   implicit none

   include 'netcdf.inc'

   integer :: nx, ny, nz, nid, icode, nrec, istatus
   character(len=*) :: fname
   logical :: there

! istatus: 1=good return, 0=error/bad return
   istatus = 1  ! assume good return

! open wrf file, and leave open for future use.

   inquire (file=trim(fname), exist=there)
   if (there) then
      icode = nf_open(trim(fname), nf_nowrite, ncid)
      if (ncid <= 0) then
         print *, 'could not open wrf file: ', trim(fname)
         istatus = 0
         return
      else
         print *, 'opened wrf file: ', trim(fname)
      end if
   else
      print *, 'could not find wrf file: ', trim(fname)
      istatus = 0
      return
   end if

   icode = nf_inq_dimid(ncid, 'time', nid)
   icode = nf_inq_dimlen(ncid, nid, nrec)
   print *, 'time (number of records): ', nrec

! verify that file has 1 or more records (dimension "time")
   if (nrec .eq. 0) then
      print *, 'error: wrf file contains no data: ', trim(fname)
      print *, 'error: cannot create output files...stopping!'
      istatus = 0
      return
   end if

! read wrf grid dimensions.

   icode = nf_inq_dimid(ncid, 'west_east', nid)
   icode = nf_inq_dimlen(ncid, nid, nx)
   print *, 'west_east: ', nx
   icode = nf_inq_dimid(ncid, 'south_north', nid)
   icode = nf_inq_dimlen(ncid, nid, ny)
   print *, 'south_north: ', ny
   icode = nf_inq_dimid(ncid, 'bottom_top', nid)
   icode = nf_inq_dimlen(ncid, nid, nz)
   print *, 'bottom_top: ', nz

   return
end

!===============================================================================

subroutine fill_wrf_grid

   use lfmgrid
   use wrfutil
   use constants

   implicit none

   include 'netcdf.inc'

   integer :: nid, icode, mapproj, i, j, k
   real, allocatable, dimension(:, :) :: ncon_pcp_tot
   real, allocatable, dimension(:, :, :) :: fld3d

! fill native map projection settings.

   icode = nf_get_att_int(ncid, nf_global, 'map_proj', mapproj)
   icode = nf_get_att_real(ncid, nf_global, 'dx', ngrid_spacingx)
   icode = nf_get_att_real(ncid, nf_global, 'dy', ngrid_spacingy)
   icode = nf_get_att_real(ncid, nf_global, 'truelat1', ntruelat1)
   icode = nf_get_att_real(ncid, nf_global, 'truelat2', ntruelat2)
   icode = nf_get_att_real(ncid, nf_global, 'stand_lon', nstdlon)

   select case (mapproj)
   case (1)
      nprojection = 'lambert conformal'
   case (2)
      nprojection = 'polar stereographic'
   case (3)
      nprojection = 'mercator'
   end select

! allocate local variables.

   allocate (ncon_pcp_tot(nx, ny))

! read model data.

   if (l_process_uv) then ! u,v
      print *, ' reading wrf u/v ', large_ngrid, large_pgrid, l_process_uv
      icode = nf_inq_varid(ncid, 'u', nid)
      if (icode .eq. 0) then
         allocate (fld3d(nx + 1, ny, nz))
         icode = nf_get_var_real(ncid, nid, fld3d)
         print *, 'u: ', icode
         do i = 1, nx - 1
            nusig(i, :, :) = (fld3d(i, :, :) + fld3d(i + 1, :, :))*0.5
         end do
         do i = nx, nx
            nusig(i, :, :) = (fld3d(i, :, :) + fld3d(i, :, :))*0.5
         end do
         deallocate (fld3d)
      end if

      icode = nf_inq_varid(ncid, 'v', nid)
      if (icode .eq. 0) then
         allocate (fld3d(nx, ny + 1, nz))
         icode = nf_get_var_real(ncid, nid, fld3d)
         print *, 'v: ', icode
         do j = 1, ny - 1
            nvsig(:, j, :) = (fld3d(:, j, :) + fld3d(:, j + 1, :))*0.5
         end do
         do j = ny, ny
            nvsig(:, j, :) = (fld3d(:, j, :) + fld3d(:, j, :))*0.5
         end do
         deallocate (fld3d)
      end if
   else
      print *, ' skipping read of wrf u/v ', large_ngrid, large_pgrid, l_process_uv
   end if

   icode = nf_inq_varid(ncid, 'p', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, npsig)
      icode = nf_inq_varid(ncid, 'pb', nid)
      if (icode .eq. 0) then
         allocate (fld3d(nx, ny, nz))
         icode = nf_get_var_real(ncid, nid, fld3d)
         npsig = npsig + fld3d
         deallocate (fld3d)
      else
         npsig = rmsg
      end if
   end if

   icode = nf_inq_varid(ncid, 't', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, ntsig)
      ntsig = (ntsig + 300.)*(npsig/p0)**kappa
   end if

   icode = nf_inq_varid(ncid, 'qvapor', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nmrsig)

   icode = nf_inq_varid(ncid, 'qcloud', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, ncldliqmr_sig)

   icode = nf_inq_varid(ncid, 'qice', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, ncldicemr_sig)

   icode = nf_inq_varid(ncid, 'qrain', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nrainmr_sig)

   icode = nf_inq_varid(ncid, 'qsnow', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nsnowmr_sig)

   icode = nf_inq_varid(ncid, 'qgraup', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, ngraupelmr_sig)

   if (c_m2z .eq. 'wrf') then
      print *, ' reading wrf refl_10cm'
      icode = nf_inq_varid(ncid, 'refl_10cm', nid)
      if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nrefl_sig)
   end if

   if (l_process_w) then ! w
      print *, ' reading wrf w ', large_ngrid, large_pgrid, l_process_w
      allocate (fld3d(nx, ny, nz + 1))
      icode = nf_inq_varid(ncid, 'w', nid)
      if (icode .eq. 0) then
         icode = nf_get_var_real(ncid, nid, fld3d)
         do k = 1, nz
            nwsig(:, :, k) = (fld3d(:, :, k) + fld3d(:, :, k + 1))*0.5
         end do
      end if
      deallocate (fld3d)
   else
      print *, ' skipping read of wrf w ', large_ngrid, large_pgrid, l_process_w
   end if ! l_process_w

   if (.true.) then ! z
      print *, ' reading wrf z'
      allocate (fld3d(nx, ny, nz + 1))
      icode = nf_inq_varid(ncid, 'ph', nid)
      if (icode .eq. 0) then
         icode = nf_get_var_real(ncid, nid, fld3d)
         do k = 1, nz
            nzsig(:, :, k) = (fld3d(:, :, k) + fld3d(:, :, k + 1))*0.5
         end do
         icode = nf_inq_varid(ncid, 'phb', nid)
         if (icode .eq. 0) then
            icode = nf_get_var_real(ncid, nid, fld3d)
            do k = 1, nz
               nzsig(:, :, k) = (nzsig(:, :, k) + (fld3d(:, :, k) + fld3d(:, :, k + 1))*0.5)/grav
            end do
         else
            nzsig = rmsg
         end if
      end if
      deallocate (fld3d)
   else
      print *, ' skipping read of wrf z'
   end if ! large_ngrid

   icode = nf_inq_varid(ncid, 'tsk', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nground_t)

   icode = nf_inq_varid(ncid, 'psfc', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, npsfc)

   icode = nf_inq_varid(ncid, 't2', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, ntsfc)

   icode = nf_inq_varid(ncid, 'q2', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nmrsfc)

   icode = nf_inq_varid(ncid, 'u10', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nusfc)

   icode = nf_inq_varid(ncid, 'v10', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nvsfc)

   icode = nf_inq_varid(ncid, 'rainc', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, ncon_pcp_tot)
   else
      ncon_pcp_tot = 0.
   end if

   icode = nf_inq_varid(ncid, 'rainnc', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, npcp_tot)
   else
      npcp_tot = 0.
   end if

   icode = nf_inq_varid(ncid, 'hgt', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nzsfc)

   icode = nf_inq_varid(ncid, 'xlat', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nlat)
   print *, 'xlat: ', icode, maxval(nlat)

   icode = nf_inq_varid(ncid, 'xlong', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nlon)
   print *, 'xlong: ', icode, maxval(nlon)

   icode = nf_inq_varid(ncid, 'swdown', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nswdown)

   icode = nf_inq_varid(ncid, 'glw', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nlwdown)

   icode = nf_inq_varid(ncid, 'olr', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nlwout)

   icode = nf_inq_varid(ncid, 'lh', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nlhflux)

   icode = nf_inq_varid(ncid, 'grdflx', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nshflux)

   icode = nf_inq_varid(ncid, 'pblh', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, npblhgt)
      if (maxval(npblhgt) <= 1.) npblhgt = rmsg
   end if

   icode = nf_close(ncid)

! fill total precip and convert from mm to m.

   npcp_tot = (npcp_tot + ncon_pcp_tot)*0.001

   deallocate (ncon_pcp_tot)

   if (fcsttime == 0.) then
      nlwdown = rmsg
      nswdown = rmsg
      nshflux = rmsg
      nlhflux = rmsg
   end if

   return
end
