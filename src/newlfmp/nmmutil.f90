module nmmutil

   implicit none

   save

   integer :: ncid

end module

!===============================================================================

subroutine get_nmm_dims(fname, nx, ny, nz)

   use nmmutil

   implicit none

   include 'netcdf.inc'

   integer :: nx, ny, nz, nid, icode
   character(len=*) :: fname
   logical :: there

! open nmm file, and leave open for future use.

   inquire (file=trim(fname), exist=there)
   if (there) then
      icode = nf_open(trim(fname), nf_nowrite, ncid)
      if (ncid <= 0) then
         print *, 'could not open nmm file: ', trim(fname)
         stop
      end if
   else
      print *, 'could not find nmm file: ', trim(fname)
      stop
   end if

! read nmm grid dimensions.

   icode = nf_inq_dimid(ncid, 'west_east', nid)
   icode = nf_inq_dimlen(ncid, nid, nx)
   icode = nf_inq_dimid(ncid, 'south_north', nid)
   icode = nf_inq_dimlen(ncid, nid, ny)
   icode = nf_inq_dimid(ncid, 'bottom_top', nid)
   icode = nf_inq_dimlen(ncid, nid, nz)

   return
end

!===============================================================================

subroutine fill_nmm_grid

   use lfmgrid
   use nmmutil
   use constants

   implicit none

   include 'netcdf.inc'

   integer :: nid, icode, k
   real, allocatable, dimension(:, :) :: fld2d
   real, allocatable, dimension(:, :, :) :: fld3d

   nprojection = 'rotated lat-lon'
   icode = nf_get_att_real(ncid, nf_global, 'cen_lat', ntruelat1)
   icode = nf_get_att_real(ncid, nf_global, 'cen_lon', nstdlon)
   icode = nf_get_att_real(ncid, nf_global, 'dx', ngrid_spacingx)
   icode = nf_get_att_real(ncid, nf_global, 'dy', ngrid_spacingy)

   icode = nf_inq_varid(ncid, 'fis', nid)
   icode = nf_get_var_real(ncid, nid, nzsfc)
   nzsfc(:, :) = nzsfc(:, :)/grav
   icode = nf_inq_varid(ncid, 'glat', nid)
   icode = nf_get_var_real(ncid, nid, nlat)
   icode = nf_inq_varid(ncid, 'glon', nid)
   icode = nf_get_var_real(ncid, nid, nlon)

   allocate (fld3d(nx, ny, nz + 1))
   icode = nf_inq_varid(ncid, 'pint', nid)
   icode = nf_get_var_real(ncid, nid, fld3d)
   icode = nf_inq_varid(ncid, 't', nid)
   icode = nf_get_var_real(ncid, nid, ntsig)
   icode = nf_inq_varid(ncid, 'q', nid)
   icode = nf_get_var_real(ncid, nid, nmrsig)

! compute height and pressure at half-levels.
   call nmm_height(nx, ny, nz, fld3d, ntsig, nmrsig, nzsfc, npsig, nzsig)

   icode = nf_inq_varid(ncid, 'u', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nusig)
   icode = nf_inq_varid(ncid, 'v', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nvsig)
   icode = nf_inq_varid(ncid, 'w', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, fld3d)
      do k = 1, nz
         nwsig(:, :, k) = (fld3d(:, :, k) + fld3d(:, :, k + 1))*0.5
      end do
      nwsfc(:, :) = fld3d(:, :, 1)
   end if

   deallocate (fld3d)

! nmm cloud water fields:
!   cwm is total condensate.
!   f_ice is fraction of cwm that is ice.
!   liquid water content is cwm - f_ice*cwm.
!   f_rain is fraction of liquid water content that is rain.
!   f_rimef is ratio of total ice to unrimed ice (>=1).
!   snow cannot be differentiated from ice in nmm output.

   icode = nf_inq_varid(ncid, 'cwm', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, ncldicemr_sig) ! cwm
      icode = nf_inq_varid(ncid, 'f_ice', nid)
      if (icode .eq. 0) then
         icode = nf_get_var_real(ncid, nid, nsnowmr_sig) ! f_ice
         icode = nf_inq_varid(ncid, 'f_rain', nid)
         if (icode .eq. 0) then
            icode = nf_get_var_real(ncid, nid, nrainmr_sig) ! f_rain
            icode = nf_inq_varid(ncid, 'f_rimef', nid)
            if (icode .eq. 0) then
               icode = nf_get_var_real(ncid, nid, ngraupelmr_sig) ! f_rimef
!           ncldliqmr_sig=ncldicemr_sig-ncldicemr_sig*nsnowmr_sig
               ncldliqmr_sig = ncldicemr_sig*(1.0 - nsnowmr_sig)

               nrainmr_sig = ncldliqmr_sig*nrainmr_sig

               ncldicemr_sig = ncldicemr_sig*nsnowmr_sig

               ngraupelmr_sig = ncldicemr_sig/ngraupelmr_sig

               nsnowmr_sig = rmsg
            else
               ncldliqmr_sig = rmsg
               nrainmr_sig = rmsg
               ncldicemr_sig = rmsg
               ngraupelmr_sig = rmsg
               nsnowmr_sig = rmsg
            end if
         end if
      end if
   end if

! accumulate precip.

!icode=nf_inq_varid(ncid,'cuprec',nid)
!if (icode .eq. 0) then
!   icode=nf_get_var_real(ncid,nid,npcp_tot)
!else
!   npcp_tot=0.
!endif

! in the nmm, acprec is the total precipitation.
! that is, acprec is the sum of cumulus and explicit precipitation.
   npcp_tot = 0.
   allocate (fld2d(nx, ny))
   icode = nf_inq_varid(ncid, 'acprec', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, fld2d)
      npcp_tot = npcp_tot + fld2d
   end if
   deallocate (fld2d)

   icode = nf_inq_varid(ncid, 'pshltr', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, npsfc)
   icode = nf_inq_varid(ncid, 'tshltr', nid)
   if ((icode .eq. 0) .and. (minval(npsfc) < rmsg)) then
      icode = nf_get_var_real(ncid, nid, ntsfc)
      ntsfc = ntsfc*(npsfc/p0)**kappa
   end if
   icode = nf_inq_varid(ncid, 'qshltr', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nmrsfc)
   icode = nf_inq_varid(ncid, 'u10', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nusfc)
   icode = nf_inq_varid(ncid, 'v10', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nvsfc)
   icode = nf_inq_varid(ncid, 'tground', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nground_t)
   icode = nf_inq_varid(ncid, 'sfcshx', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nshflux)
   icode = nf_inq_varid(ncid, 'sfclhx', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nlhflux)
   icode = nf_inq_varid(ncid, 'pblh', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, npblhgt)
   icode = nf_inq_varid(ncid, 'rlwtoa', nid)
   if (icode .eq. 0) icode = nf_get_var_real(ncid, nid, nlwout)

   if (fcsttime == 0) then
      nshflux = rmsg
      nlhflux = rmsg
      npblhgt = rmsg
   end if
   if (maxval(npblhgt) < 1.) npblhgt = rmsg

   return
end

!===============================================================================

subroutine nmm_height(nx, ny, nz, npri, ntp, nmr, sht, npr, nht)

   use constants

   implicit none

!real, parameter :: r=287.,g=9.8,rog=r/g

   integer :: nx, ny, nz, i, j, k

   real :: alnp(nx, ny)
   real, dimension(nx, ny) :: sht
   real, dimension(nx, ny, nz) :: ntp, nmr, npr, nht
   real, dimension(nx, ny, nz + 1) :: npri, nhti

! compute pressure at half levels in pa.

   do k = 1, nz
      npr(:, :, k) = (npri(:, :, k) + npri(:, :, k + 1))*0.5
   end do

! calculate height at full levels using hypsometric eqn.

   nhti(:, :, 1) = sht(:, :)
   do k = 1, nz
      nhti(:, :, k + 1) = nhti(:, :, k) + rog*(ntp(:, :, k) + 0.61*nmr(:, :, k)) &
                          *alog(npri(:, :, k)/npri(:, :, k + 1))
   end do

! logarithmic interpolation to half-levels.

   do k = 1, nz
      alnp(:, :) = alog((npr(:, :, k)/npri(:, :, k))) &
                   /alog(npri(:, :, k + 1)/npri(:, :, k))
      nht(:, :, k) = nhti(:, :, k) + (nhti(:, :, k + 1) - nhti(:, :, k))*alnp(:, :)
   end do

   return
end
