module lapsdata
!!
!! configures the list of laps variables and their associated grib-2
!! table entries

   implicit none

   ! 3d isobaric variable metadata
   type var_iso3d
      integer              :: qbal ! 0=use standard, 1=use balance output
      character(len=3)     :: ext  ! subdirectory and extension
      character(len=3)     :: var  ! netcdf variable name
      real                 :: pmin ! minimum pressure level to output
      real                 :: pmax ! maximum pressure level to output
      real                 :: conv_fac ! used to convert units to grib standard
      integer              :: scale_fac ! scale factor when gribing (power of 10)
      integer              :: discipline
      integer              :: cat_id
      integer              :: parm_id
   end type var_iso3d

   type var_2d
      character(len=3)     :: ext
      character(len=3)     :: var
      real                 :: conv_fac
      integer              :: scale_fac
      integer              :: lev1_type
      integer              :: lev1_scale
      integer              :: lev1_value
      integer              :: lev2_type
      integer              :: lev2_scale
      integer              :: lev2_value
      integer              :: discipline
      integer              :: cat_id
      integer              :: parm_id
   end type var_2d

   type var_a2d
      character(len=3)     :: ext
      character(len=3)     :: var
      real                 :: conv_fac
      integer              :: scale_fac
      integer              :: lev1_type
      integer              :: lev1_scale
      integer              :: lev1_value
      integer              :: lev2_type
      integer              :: lev2_scale
      integer              :: lev2_value
      integer              :: discipline
      integer              :: cat_id
      integer              :: parm_id

      integer              :: eyear
      integer              :: emon
      integer              :: eday
      integer              :: ehour
      integer              :: emin
      integer              :: esec
      integer              :: ntimes
      integer              :: ntimes_miss
      integer              :: stattype
      integer              :: periodtype
      integer              :: etime_unit
      integer              :: etime_value
   end type var_a2d

   integer, parameter     ::  max_iso = 40
   integer, parameter     ::  max_2d = 50
   integer, parameter     ::  max_a2d = 10

   type(var_iso3d)        :: meta_iso3d(max_iso)
   type(var_2d)           :: meta_2d(max_2d)
   type(var_a2d)          :: meta_accum2d(max_a2d)
   integer                :: n_iso, n_2d, n_accum2d

contains

   subroutine get_data_config(laps_data_root, vname)

      implicit none
      character(len=*), intent(in)    :: laps_data_root
      character(len=*), intent(in)    :: vname
      character(len=256)             :: vtab
      logical                        :: file_exists
      integer, parameter              :: lun = 10
      character(len=80)              :: fline
      integer                        :: fstat, line_len
      integer  :: qbal, scale_fac, discipline, cat_id, parm_id
      integer  :: lev1_type, lev1_scale, lev1_value
      integer  :: lev2_type, lev2_scale, lev2_value
      real     :: conv_fac, pmin, pmax
      character(len=3) :: var, ext
      vtab = trim(laps_data_root)//'/static/'//trim(vname)
      inquire (file=vtab, exist=file_exists)

      if (.not. file_exists) then
         print *, "-- config file not found:  ", trim(vtab)
         stop 'no config file'
      end if

      n_iso = 0
      n_2d = 0
      n_accum2d = 0
      open (file=vtab, unit=lun, access='sequential', form='formatted', status='old')
      readfile: do
         read (lun, '(a)', iostat=fstat) fline
         if (fstat .eq. 0) then
            line_len = len_trim(fline)
            if (fline(1:2) .eq. '3d') then
               read (fline(3:line_len), *, iostat=fstat) qbal, ext, var, pmin, pmax, conv_fac, scale_fac, &
                  discipline, cat_id, parm_id
               if (fstat .eq. 0) then
                  n_iso = n_iso + 1
                  if (n_iso .gt. max_iso) then
                     print *, "number of 3d variables requested exceeds max_iso parameter!"
                     stop 'module_lapsdata.f90'
                  end if
                  meta_iso3d(n_iso)%qbal = qbal
                  meta_iso3d(n_iso)%ext = ext
                  meta_iso3d(n_iso)%var = var
                  meta_iso3d(n_iso)%pmin = pmin
                  meta_iso3d(n_iso)%pmax = pmax
                  meta_iso3d(n_iso)%conv_fac = conv_fac
                  meta_iso3d(n_iso)%scale_fac = scale_fac
                  meta_iso3d(n_iso)%discipline = discipline
                  meta_iso3d(n_iso)%cat_id = cat_id
                  meta_iso3d(n_iso)%parm_id = parm_id
               else
                  print *, "-- format problem with 3d variable line: ", trim(fline)
               end if

            elseif (fline(1:2) .eq. '2d') then
               read (fline(3:line_len), *, iostat=fstat) ext, var, conv_fac, scale_fac, &
                  lev1_type, lev1_scale, lev1_value, lev2_type, lev2_scale, lev2_value, &
                  discipline, cat_id, parm_id
               if (fstat .eq. 0) then
                  n_2d = n_2d + 1
                  if (n_2d .gt. max_2d) then
                     print *, "number of 2d variables requested exceeds max_2d parameter!"
                     stop 'module_lapsdata.f90'
                  end if
                  meta_2d(n_2d)%ext = ext
                  meta_2d(n_2d)%var = var
                  meta_2d(n_2d)%conv_fac = conv_fac
                  meta_2d(n_2d)%scale_fac = scale_fac
                  meta_2d(n_2d)%lev1_type = lev1_type
                  meta_2d(n_2d)%lev1_scale = lev1_scale
                  meta_2d(n_2d)%lev1_value = lev1_value
                  meta_2d(n_2d)%lev2_type = lev2_type
                  meta_2d(n_2d)%lev2_scale = lev2_scale
                  meta_2d(n_2d)%lev2_value = lev2_value
                  meta_2d(n_2d)%discipline = discipline
                  meta_2d(n_2d)%cat_id = cat_id
                  meta_2d(n_2d)%parm_id = parm_id

               else
                  print *, "-- format problem with 2d variable line: ", trim(fline)
               end if

            elseif (fline(1:3) .eq. 'a2d') then
               read (fline(4:line_len), *, iostat=fstat) ext, var, conv_fac, scale_fac, &
                  lev1_type, lev1_scale, lev1_value, lev2_type, lev2_scale, lev2_value, &
                  discipline, cat_id, parm_id
               if (fstat .eq. 0) then
                  n_accum2d = n_accum2d + 1

                  if (n_accum2d .gt. max_a2d) then
                     print *, "number of accumulated 2d variables requested exceeds max_a2d parameter!"
                     stop 'module_lapsdata.f90'
                  end if

                  meta_accum2d(n_accum2d)%ext = ext
                  meta_accum2d(n_accum2d)%var = var
                  meta_accum2d(n_accum2d)%conv_fac = conv_fac
                  meta_accum2d(n_accum2d)%scale_fac = scale_fac
                  meta_accum2d(n_accum2d)%lev1_type = lev1_type
                  meta_accum2d(n_accum2d)%lev1_scale = lev1_scale
                  meta_accum2d(n_accum2d)%lev1_value = lev1_value
                  meta_accum2d(n_accum2d)%lev2_type = lev2_type
                  meta_accum2d(n_accum2d)%lev2_scale = lev2_scale
                  meta_accum2d(n_accum2d)%lev2_value = lev2_value
                  meta_accum2d(n_accum2d)%discipline = discipline
                  meta_accum2d(n_accum2d)%cat_id = cat_id
                  meta_accum2d(n_accum2d)%parm_id = parm_id
                  ! octet 35-40, reset to current runtime in laps2grib program
                  meta_accum2d(n_accum2d)%eyear = 0
                  meta_accum2d(n_accum2d)%emon = 0
                  meta_accum2d(n_accum2d)%eday = 0
                  meta_accum2d(n_accum2d)%ehour = 0
                  meta_accum2d(n_accum2d)%emin = 0
                  ! octet 41-43, 47-49
                  meta_accum2d(n_accum2d)%esec = 0
                  meta_accum2d(n_accum2d)%ntimes = 1      !42 indicator of time unit (see octet 49)
                  meta_accum2d(n_accum2d)%ntimes_miss = 0 !43 missing precip causes accumulation to reset
                  meta_accum2d(n_accum2d)%stattype = 1 !47 1=accumulation, 0=avg
                  meta_accum2d(n_accum2d)%periodtype = 2 !48 2=succesive runs have same start time of accum
                  meta_accum2d(n_accum2d)%etime_unit = 0 !49 0=min, 1=hour, 2=day
                  ! octet 50, reset to current runtime in laps2grib program
                  meta_accum2d(n_accum2d)%etime_value = 60 !50 length of time range (see octet 49)

               else
                  print *, "-- format problem with accumulated 2d variable line: ", trim(fline)
               end if

            end if

         elseif (fstat .lt. 0) then
            exit readfile
         end if

      end do readfile

      if ((n_iso .le. 0) .and. (n_2d .le. 0)) then
         print *, "-- no variables configured in config file. using default values."
         call config_data_static
      end if
      return
   end subroutine get_data_config

   subroutine config_data_static

      ! this sets up the list of variables to
      ! process.  eventually, we should make this some sort
      ! of user-editable table that gets read in.

      implicit none

      ! initialize the counter for the 3d variables
      n_iso = 0

      ! fill in the metadata for each variable, then increment the counter

      ! temperature
      n_iso = n_iso + 1
      meta_iso3d(n_iso)%qbal = 0
      meta_iso3d(n_iso)%ext = 'lt1'
      meta_iso3d(n_iso)%var = 't3 '
      meta_iso3d(n_iso)%pmin = 1000.
      meta_iso3d(n_iso)%pmax = 105000.
      meta_iso3d(n_iso)%conv_fac = 1.
      meta_iso3d(n_iso)%scale_fac = 1
      meta_iso3d(n_iso)%discipline = 0
      meta_iso3d(n_iso)%cat_id = 0
      meta_iso3d(n_iso)%parm_id = 0

      ! height
      n_iso = n_iso + 1
      meta_iso3d(n_iso)%qbal = 0
      meta_iso3d(n_iso)%ext = 'lt1'
      meta_iso3d(n_iso)%var = 'ht '
      meta_iso3d(n_iso)%pmin = 1000.
      meta_iso3d(n_iso)%pmax = 105000.
      meta_iso3d(n_iso)%conv_fac = 1.
      meta_iso3d(n_iso)%scale_fac = 0
      meta_iso3d(n_iso)%discipline = 0
      meta_iso3d(n_iso)%cat_id = 3
      meta_iso3d(n_iso)%parm_id = 5

      ! rh wrt liquid
      n_iso = n_iso + 1
      meta_iso3d(n_iso)%qbal = 0
      meta_iso3d(n_iso)%ext = 'lh3'
      meta_iso3d(n_iso)%var = 'rhl'
      meta_iso3d(n_iso)%pmin = 1000.
      meta_iso3d(n_iso)%pmax = 105000.
      meta_iso3d(n_iso)%conv_fac = 1.
      meta_iso3d(n_iso)%scale_fac = 1
      meta_iso3d(n_iso)%discipline = 0
      meta_iso3d(n_iso)%cat_id = 1
      meta_iso3d(n_iso)%parm_id = 1

      ! u wind
      n_iso = n_iso + 1
      meta_iso3d(n_iso)%qbal = 0
      meta_iso3d(n_iso)%ext = 'lw3'
      meta_iso3d(n_iso)%var = 'u3 '
      meta_iso3d(n_iso)%pmin = 1000.
      meta_iso3d(n_iso)%pmax = 105000.
      meta_iso3d(n_iso)%conv_fac = 1.
      meta_iso3d(n_iso)%scale_fac = 1
      meta_iso3d(n_iso)%discipline = 0
      meta_iso3d(n_iso)%cat_id = 2
      meta_iso3d(n_iso)%parm_id = 2

      ! v wind
      n_iso = n_iso + 1
      meta_iso3d(n_iso)%qbal = 0
      meta_iso3d(n_iso)%ext = 'lw3'
      meta_iso3d(n_iso)%var = 'v3 '
      meta_iso3d(n_iso)%pmin = 1000.
      meta_iso3d(n_iso)%pmax = 105000.
      meta_iso3d(n_iso)%conv_fac = 1.
      meta_iso3d(n_iso)%scale_fac = 1
      meta_iso3d(n_iso)%discipline = 0
      meta_iso3d(n_iso)%cat_id = 2
      meta_iso3d(n_iso)%parm_id = 3

      ! 2 d variables
      n_2d = 0

      ! mslp
      n_2d = n_2d + 1
      meta_2d(n_2d)%ext = 'lsx'
      meta_2d(n_2d)%var = 'msl'
      meta_2d(n_2d)%conv_fac = 1.
      meta_2d(n_2d)%scale_fac = 0
      meta_2d(n_2d)%lev1_type = 1
      meta_2d(n_2d)%lev1_scale = 0
      meta_2d(n_2d)%lev1_value = 0
      meta_2d(n_2d)%lev2_type = 255
      meta_2d(n_2d)%lev2_scale = 255
      meta_2d(n_2d)%lev2_value = 255
      meta_2d(n_2d)%discipline = 0
      meta_2d(n_2d)%cat_id = 3
      meta_2d(n_2d)%parm_id = 1

      ! psfc
      n_2d = n_2d + 1
      meta_2d(n_2d)%ext = 'lsx'
      meta_2d(n_2d)%var = 'ps '
      meta_2d(n_2d)%conv_fac = 1.
      meta_2d(n_2d)%scale_fac = 0
      meta_2d(n_2d)%lev1_type = 1
      meta_2d(n_2d)%lev1_scale = 0
      meta_2d(n_2d)%lev1_value = 0
      meta_2d(n_2d)%lev2_type = 255
      meta_2d(n_2d)%lev2_scale = 255
      meta_2d(n_2d)%lev2_value = 255
      meta_2d(n_2d)%discipline = 0
      meta_2d(n_2d)%cat_id = 3
      meta_2d(n_2d)%parm_id = 0

      ! 2m temp
      n_2d = n_2d + 1
      meta_2d(n_2d)%ext = 'lsx'
      meta_2d(n_2d)%var = 't  '
      meta_2d(n_2d)%conv_fac = 1.
      meta_2d(n_2d)%scale_fac = 1
      meta_2d(n_2d)%lev1_type = 103
      meta_2d(n_2d)%lev1_scale = 0
      meta_2d(n_2d)%lev1_value = 2
      meta_2d(n_2d)%lev2_type = 255
      meta_2d(n_2d)%lev2_scale = 255
      meta_2d(n_2d)%lev2_value = 255
      meta_2d(n_2d)%discipline = 0
      meta_2d(n_2d)%cat_id = 0
      meta_2d(n_2d)%parm_id = 0

      ! 2m dewpoint
      n_2d = n_2d + 1
      meta_2d(n_2d)%ext = 'lsx'
      meta_2d(n_2d)%var = 'td '
      meta_2d(n_2d)%conv_fac = 1.
      meta_2d(n_2d)%scale_fac = 1
      meta_2d(n_2d)%lev1_type = 103
      meta_2d(n_2d)%lev1_scale = 0
      meta_2d(n_2d)%lev1_value = 2
      meta_2d(n_2d)%lev2_type = 255
      meta_2d(n_2d)%lev2_scale = 255
      meta_2d(n_2d)%lev2_value = 255
      meta_2d(n_2d)%discipline = 0
      meta_2d(n_2d)%cat_id = 0
      meta_2d(n_2d)%parm_id = 6

      ! 2m rh
      n_2d = n_2d + 1
      meta_2d(n_2d)%ext = 'lsx'
      meta_2d(n_2d)%var = 'rh '
      meta_2d(n_2d)%conv_fac = 1.
      meta_2d(n_2d)%scale_fac = 0
      meta_2d(n_2d)%lev1_type = 103
      meta_2d(n_2d)%lev1_scale = 0
      meta_2d(n_2d)%lev1_value = 2
      meta_2d(n_2d)%lev2_type = 255
      meta_2d(n_2d)%lev2_scale = 255
      meta_2d(n_2d)%lev2_value = 255
      meta_2d(n_2d)%discipline = 0
      meta_2d(n_2d)%cat_id = 1
      meta_2d(n_2d)%parm_id = 1

      ! 10m u
      n_2d = n_2d + 1
      meta_2d(n_2d)%ext = 'lsx'
      meta_2d(n_2d)%var = 'u  '
      meta_2d(n_2d)%conv_fac = 1.
      meta_2d(n_2d)%scale_fac = 0
      meta_2d(n_2d)%lev1_type = 103
      meta_2d(n_2d)%lev1_scale = 0
      meta_2d(n_2d)%lev1_value = 10
      meta_2d(n_2d)%lev2_type = 255
      meta_2d(n_2d)%lev2_scale = 255
      meta_2d(n_2d)%lev2_value = 255
      meta_2d(n_2d)%discipline = 0
      meta_2d(n_2d)%cat_id = 2
      meta_2d(n_2d)%parm_id = 2

      ! 10m v
      n_2d = n_2d + 1
      meta_2d(n_2d)%ext = 'lsx'
      meta_2d(n_2d)%var = 'v  '
      meta_2d(n_2d)%conv_fac = 1.
      meta_2d(n_2d)%scale_fac = 0
      meta_2d(n_2d)%lev1_type = 103
      meta_2d(n_2d)%lev1_scale = 0
      meta_2d(n_2d)%lev1_value = 10
      meta_2d(n_2d)%lev2_type = 255
      meta_2d(n_2d)%lev2_scale = 255
      meta_2d(n_2d)%lev2_value = 255
      meta_2d(n_2d)%discipline = 0
      meta_2d(n_2d)%cat_id = 2
      meta_2d(n_2d)%parm_id = 3

      return
   end subroutine config_data_static

end module lapsdata
