!dis
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis
!dis

module wfo_models

   ! module that contains routines and variable pertaining to the setup
   ! of the background model grid being processed by sbnprep.

   use map_utils
   implicit none
   include 'netcdf.inc'
   integer                :: nfstatus
   integer                :: maxlevels = 100
   integer                :: maxtimes = 100

   public open_wfofile, close_wfofile, get_wfomodel_var_levels, &
      get_wfomodel_fcsttimes, get_wfomodel_proj, get_wfomodel_var_inv, &
      read_wfomodel_data, get_wfomodel_topo
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine open_wfofile(nfname, nfid, istatus)

      ! given a netcdf file name, this routine opens it and returns the
      ! integer netcdf file handle.

      implicit none
      character(len=256), intent(in)  :: nfname  ! input file name
      integer, intent(out)             :: nfid    ! output unit number
      integer, intent(out)             :: istatus ! status flag (1=succes)

      istatus = 1
      nfstatus = nf_open(nfname, nf_nowrite, nfid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem with netcdf file:', trim(nfname)
         print *, 'netcdf error = ', nfstatus
         istatus = 0
      end if
      return
   end subroutine open_wfofile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine close_wfofile(nfid, istatus)

      ! closes an open netcdf file given the integer netcdf file handle.

      implicit none
      integer, intent(in)            :: nfid
      integer, intent(out)           :: istatus

      istatus = 1
      nfstatus = nf_close(nfid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem closing netcdf file.'
         print *, 'netcdf error = ', nfstatus
         istatus = 0
      end if
      return
   end subroutine close_wfofile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_wfomodel_proj(nfid, mname, proj, istatus)

      ! populates a projection information data structure defined
      ! by module_map_utils.f for an already open wfo netcdf model
      ! file.  use open_wfofile to open the file before calling this
      ! routine.  getting projection info is a bit tricky, because the awips
      ! model files have often been clipped and do not give the grid spacing
      ! at the actual true latitude of the projection.  hence, the model name,
      ! which corresponds to the last two elements of the directory the model
      ! file is located in (e.g., "conus212/mesoeta"), is used to make some
      ! educated guesses.

      implicit none
      integer, intent(in)           :: nfid
      character(len=132), intent(in) :: mname
      type(proj_info), intent(out)  :: proj  ! declared via "use map_utils"
      integer, intent(out)          :: istatus

      integer                       :: nx, ny
      integer                       :: vid, attnum, attid
      real                          :: truelat1, truelat2
      real                          :: stdlon
      real                          :: latsw, lonsw
      real                          :: latne, lonne
      real                          :: latne_c, lonne_c
      real                          :: difflat, difflon
      real                          :: dx, dy
      real                          :: truedx
      real                          :: latdxdy, londxdy
      character(len=132)            :: gproj

      istatus = 1
      gproj(1:132) = ' '
      ! get the horizontal dimensions
      nfstatus = nf_inq_dimid(nfid, 'x', vid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem getting x dimension id'
         istatus = 0
         return
      end if

      nfstatus = nf_inq_dimlen(nfid, vid, nx)
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem getting x dimension'
         istatus = 0
         return
      end if

      nfstatus = nf_inq_dimid(nfid, 'y', vid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem getting y dimension id'
         istatus = 0
         return
      end if

      nfstatus = nf_inq_dimlen(nfid, vid, ny)
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem getting y dimension'
         istatus = 0
         return
      end if

      ! get projection info
      nfstatus = nf_inq_attid(nfid, attid, 'projname', attnum)
      nfstatus = nf_get_att_text(nfid, attid, 'projname', gproj)
      nfstatus = nf_inq_attid(nfid, attid, 'centrallat', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'centrallat', truelat1)
      nfstatus = nf_inq_attid(nfid, attid, 'centrallon', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'centrallon', stdlon)
      nfstatus = nf_inq_attid(nfid, attid, 'rotation', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'rotation', truelat2)
      nfstatus = nf_inq_attid(nfid, attid, 'lat00', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'lat00', latsw)
      nfstatus = nf_inq_attid(nfid, attid, 'lon00', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'lon00', lonsw)
      nfstatus = nf_inq_attid(nfid, attid, 'latnxny', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'latnxny', latne)
      nfstatus = nf_inq_attid(nfid, attid, 'lonnxny', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'lonnxny', lonne)
      nfstatus = nf_inq_attid(nfid, attid, 'dxkm', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'dxkm', dx)
      nfstatus = nf_inq_attid(nfid, attid, 'dykm', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'dykm', dy)
      nfstatus = nf_inq_attid(nfid, attid, 'latdxdy', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'latdxdy', latdxdy)
      nfstatus = nf_inq_attid(nfid, attid, 'londxdy', attnum)
      nfstatus = nf_get_att_real(nfid, attid, 'londxdy', londxdy)

      ! determine grid spacing at true latitude.  we have to do some
      ! level of hard coding for efficiency here, because the sbn data feed
      ! does not provide grid spacing at the true latitude.  rather, it provides
      ! a "self-computed" grid spacing at the approximate center of the domain.  so,
      ! we will make a guess at the true dx, then use the mapping routines to
      ! verify that we get the correct coordinate conversion when using the "guessed"
      ! value. down the road, we could put in some iterative method to do this using
      ! the provided dx/dy values as a starting point.

      if (mname(1:8) .eq. 'conus211') then
         !  this is the ncep 211 grid, which is not really 48km
         truedx = 81270.50  ! meters at 25.0 n
         call map_set(proj_lc, latsw, lonsw, truedx, stdlon, truelat1, truelat2, &
                      nx, ny, proj)
      else if (mname(1:8) .eq. 'conus212') then
         truedx = 40635.25
         call map_set(proj_lc, latsw, lonsw, truedx, stdlon, truelat1, truelat2, &
                      nx, ny, proj)
      else if (mname(1:8) .eq. 'conus215') then
         truedx = 20317.625
         call map_set(proj_lc, latsw, lonsw, truedx, stdlon, truelat1, truelat2, &
                      nx, ny, proj)
      else if (mname(1:6) .eq. 'latlon') then
         ! true dx for this grid is in lat/lon increment
         truedx = 1.25  ! degrees
         ! for latlon projections, the latitude increment is stored in truelat1 and
         ! the longitudinal increment is put into stdlon.  for the avn data on sbn,
         ! the origin in this data set is listed as the nw corner, but the
         ! array is actually set up with the origin at the southwest.
         call map_set(proj_latlon, latsw, lonsw, 0., truedx, truedx, 0., nx, ny, proj)
      else
         print *, 'model type not yet supported: ', trim(mname)
         istatus = 0
         return
      end if

      ! verify that the computed upper-right corner matches the specified upper right
      ! corner
      call ij_to_latlon(proj, float(nx), float(ny), latne_c, lonne_c)
      print *, 'comp latne/lonne = ', latne_c, lonne_c
      print *, 'spec latne/lonne = ', latne, lonne
      difflat = abs(latne_c - latne)
      difflon = abs(lonne_c - lonne)
      if ((difflat .gt. 0.001) .or. (difflon .gt. 0.001)) then
         print *, 'problem with projection information:'
         print *, 'specified lat/lon at nx/ny = ', latne, lonne
         print *, 'computed lat/lon at nx/ny  = ', latne_c, lonne_c
         istatus = 0
         return
      end if

      ! if we made it this far, print out some diagnostic information.
      print *, 'wfo model map projection info for ', trim(mname)
      print *, '------------------------------------------------------------'
      print *, 'projection type:              ', trim(gproj)
      print *, 'southwest lat/lon:            ', proj%lat1, proj%lon1
      print *, 'dlat/dlon (latlon proj only): ', proj%dlat, proj%dlon
      print *, 'standard lon:                 ', proj%stdlon
      print *, 'truelat1/truelat2:            ', proj%truelat1, proj%truelat2
      print *, 'hemi (1 for nh, -1 for sh):   ', proj%hemi
      print *, 'cone factor:                  ', proj%cone
      print *, 'pole point (i/j):             ', proj%polei, proj%polej
      return
   end subroutine get_wfomodel_proj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_wfomodel_var_levels(nfid, varname, varid, n_levels, levels_c, &
                                      n_plevels, plevels_r, pbot_ind, ptop_ind, &
                                      havesfc, sfc_ind, &
                                      istatus)

      ! given a file handle for a previously opened netcdf wfo
      ! awips-format model file and a character variable name,
      ! this routine will return various pieces of information
      ! about the variable, including number of total levels,
      ! an array of level ids, number of pressure levels, array
      ! of pressure levels in pa, and a flag as to whether or
      ! not there is a surface value.  it also returns the netcdf
      ! integer variable id.

      implicit none
      integer, intent(in)            :: nfid
      character(len=10), intent(in)   :: varname
      integer, intent(out)           :: varid
      integer, intent(out)           :: n_levels
      character(len=10), intent(out)  :: levels_c(maxlevels)
      integer, intent(out)           :: n_plevels
      real, intent(out)              :: plevels_r(maxlevels)
      integer, intent(out)           :: ptop_ind
      integer, intent(out)           :: pbot_ind
      integer, intent(out)           :: sfc_ind
      logical, intent(out)           :: havesfc
      integer, intent(out)           :: istatus

      character(len=32)              :: levname
      integer                        :: levid
      integer                        :: dimid(4)
      character(len=10)              :: sfc_level
      character(len=10)              :: level_txt
      integer                        :: k, kp
      real                           :: press_mb
      character(len=2)               :: dummy2
      istatus = 1
      n_levels = 0
      n_plevels = 0
      havesfc = .false.
      varid = -1
      pbot_ind = -1
      ptop_ind = -1
      sfc_ind = -1

      if ((trim(varname) .eq. 't') .or. &
          (trim(varname) .eq. 'rh')) then
         sfc_level = 'fhag 2    '
      else if ((trim(varname) .eq. 'uw') .or. &
               (trim(varname) .eq. 'vw')) then
         sfc_level = 'fhag 10   '
      else if ((trim(varname) .eq. 'emsp') .or. &
               (trim(varname) .eq. 'pmsl') .or. &
               (trim(varname) .eq. 'mmsp')) then
         sfc_level = 'msl       '
      else
         sfc_level = 'unknown   '
      end if

      ! get the variable id for the variable requested
      nfstatus = nf_inq_varid(nfid, varname, varid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'variable not found: ', varname
         istatus = 0
         return
      end if

      ! build the name of the level variable
      levname = trim(varname)//'levels'
      nfstatus = nf_inq_varid(nfid, levname, levid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'no levels variable found for ', trim(levname)
         istatus = 0
         return
      end if
      nfstatus = nf_inq_vardimid(nfid, levid, dimid)
      nfstatus = nf_inq_dimlen(nfid, dimid(2), n_levels)
      if (n_levels .gt. 0) then
         nfstatus = nf_get_var_text(nfid, levid, levels_c)
         ! scan for number of levels on pressure surfaces
         do k = 1, n_levels
            level_txt = levels_c(k)
            if (level_txt(1:2) .eq. 'mb') then
               n_plevels = n_plevels + 1

               ! assume the pressure level data is contiguous from
               ! the lower atmosphere upward in the array
               if (pbot_ind .le. 0) pbot_ind = k
            else if (level_txt .eq. sfc_level) then
               havesfc = .true.
               if (sfc_ind .le. 0) sfc_ind = k
            end if
         end do
         if (n_plevels .gt. 0) then
            ptop_ind = pbot_ind + n_plevels - 1
            kp = 1
            do k = 1, n_levels
               level_txt = levels_c(k)
               if (level_txt(1:2) .eq. 'mb') then
                  read (level_txt, '(a2,1x,f7.0)') dummy2, press_mb
                  plevels_r(kp) = press_mb*100.
                  kp = kp + 1
               end if
            end do
         end if
      end if
      print *, 'level info for variable ', trim(varname)
      print *, '-------------------------------------------------'
      print *, 'num total levels:        ', n_levels
      print *, 'level ids: ', levels_c(1:n_levels)
      print *, 'num pressure levels:     ', n_plevels
      print *, 'pressure levels (pa):', plevels_r(1:n_plevels)
      print *, 'found surface value:     ', havesfc
      print *, 'kbotp/ktopp/ksfc:   ', pbot_ind, ptop_ind, sfc_ind
      return
   end subroutine get_wfomodel_var_levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_wfomodel_fcsttimes(nfid, ntimes, fcstsec, istatus)

      ! given the integer file handle of a previously opened
      ! awips netcdf model file, this routine returns the  number
      ! of output forecast times and each valtime-reftime in seconds.

      implicit none

      integer, intent(in)                :: nfid
      integer, intent(out)               :: ntimes
      integer, intent(out)               :: fcstsec(maxtimes)
      integer, intent(out)               :: istatus

      integer                            :: timeid
      istatus = 1

      nfstatus = nf_inq_dimid(nfid, 'n_valtimes', timeid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'could not obtain n_valtimes variable id'
         istatus = 0
         return
      end if
      nfstatus = nf_inq_dimlen(nfid, timeid, ntimes)
      nfstatus = nf_inq_varid(nfid, 'valtimeminusreftime', timeid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'could not obtain valtimeminusreftime'
         istatus = 0
         return
      end if
      nfstatus = nf_get_vara_int(nfid, timeid, 1, ntimes, fcstsec(1:ntimes))
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem obtaining the valid times'
         istatus = 0
      end if
      print *, 'forecast hours found: ', fcstsec(1:ntimes)/3600
      return
   end subroutine get_wfomodel_fcsttimes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_wfomodel_data(nfid, vid, proj, time_ind, &
                                 start_lev, stop_lev, data, istatus)

      ! subroutine to read a specific variable (given by integer
      ! variable id) from a specific file (given by integer file id)
      ! from the start/stop values for each dimension of the array.   the
      ! routine presumes you have already called open_wfofile and
      ! get_wfomodel_var_levels for the variable you want to obtain.

      implicit none
      integer, intent(in)                   :: nfid
      integer, intent(in)                   :: vid
      type(proj_info), intent(in)            :: proj
      integer, intent(in)                    :: time_ind
      integer, intent(in)                    :: start_lev
      integer, intent(in)                    :: stop_lev
      real, intent(out)                      :: data(proj%nx, proj%ny, &
                                                     stop_lev - start_lev + 1)
      integer, intent(out)                  :: istatus

      integer                               :: start_ind(4)
      integer                               :: count_ind(4)

      istatus = 1
      start_ind(4) = time_ind
      count_ind(4) = 1
      start_ind(3) = start_lev
      count_ind(3) = stop_lev - start_lev + 1
      start_ind(2) = 1
      count_ind(2) = proj%ny
      start_ind(1) = 1
      count_ind(1) = proj%nx

      nfstatus = nf_get_vara_real(nfid, vid, start_ind, count_ind, data)
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem obtaining variable for vid = ', vid
         istatus = 0
      end if
      return
   end subroutine read_wfomodel_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_wfomodel_topo(nfid, proj, topo, istatus)

      ! this routine will obtain the topography data from
      ! an awips model file.  it assumes you have already opened
      ! the file to get nfid and called the get_wfomodel_proj routine
      ! to get the proj structure and that topo is already allocated

      implicit none
      integer, intent(in)           :: nfid
      type(proj_info), intent(in)    :: proj
      real, intent(out)             :: topo(proj%nx, proj%ny)
      integer, intent(out)          :: istatus
      integer                       :: topoid
      istatus = 1

      nfstatus = nf_inq_varid(nfid, 'statictopo', topoid)
      if (nfstatus .ne. nf_noerr) then
         print *, nf_strerror(nfstatus)
         print *, 'nf_inq_varid statictopo'
         istatus = 0
         return
      end if

      nfstatus = nf_get_var_real(nfid, topoid, topo)
      if (nfstatus .ne. nf_noerr) then
         print *, nf_strerror(nfstatus)
         print *, 'nf_get_var_real statictopo'
         istatus = 0
         return
      end if

      return
   end subroutine get_wfomodel_topo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_wfomodel_var_inv(nfid, varname, nlevs, ntimes, inv, istatus)

      ! subroutine to return the inventory of a specific variable
      ! from an awips netcdf model bigfile.  the inventory is returned
      ! as a 2d array of logical flags dimensioned (n_levels,ntimes)

      implicit none
      integer, intent(in)              :: nfid
      character(len=10), intent(in)     :: varname
      integer, intent(in)               :: nlevs
      integer, intent(in)               :: ntimes
      logical, intent(out)              :: inv(nlevs, ntimes)
      integer, intent(out)              :: istatus

      character(len=1), allocatable     :: invchar(:, :)
      character(len=16)                :: invname
      integer                          :: invid

      istatus = 1
      invname = trim(varname)//'inventory'
      nfstatus = nf_inq_varid(nfid, invname, invid)
      if (nfstatus .ne. nf_noerr) then
         print *, 'no variable found: ', trim(invname)
         istatus = 0
         return
      end if

      allocate (invchar(nlevs, ntimes))
      nfstatus = nf_get_var_text(nfid, invid, invchar)
      if (nfstatus .ne. nf_noerr) then
         print *, 'problem getting inventory.'
         istatus = 0
         return
      end if

      inv(:, :) = .false.
      where (invchar .eq. '1') inv = .true.
      deallocate (invchar)
      return
   end subroutine get_wfomodel_var_inv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module wfo_models

