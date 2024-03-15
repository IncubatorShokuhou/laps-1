
                subroutine start_up_mix(ni, nj, niveau, intliq700)

          use lfmgrid, only: hzsig,htsig,hpsig,hmrsig,husig,hvsig,hwsig,hrefl_sig,hgraupelmr_sig,hcldliqmr_sig,intliqwater,llat,llon

                   w_to_omegaf(w, pres, scaleht) = -(w*pressure_pa)/scale_height

                !!-------------------------------------------------
                !! defined structure needed by read_config_file()
                !!-------------------------------------------------
                   type mystruct
                      integer :: ip2;                                  !! ip2
                      integer :: npas;                                                         !! npas
                      integer :: yyymmddhh;                                        !! yyyymmddhh
                      integer :: hh_init;                                                !! hh_init
                      integer :: hh_fcst;                                                !! hh_prevu
                      integer :: nbr_level_2keep;                        !! nbr_niveau_voulu
                      integer :: nbr_level_avail;                        !! nbr_niveau_disponible
                      real :: resolution;                                                !! resolution
                      real :: thres_meso;                                                !! seuil_meso
                      real :: thres_tornade;                                        !! seuil_tornade
                      real :: thres_severe_updraft(0:1)         !! seuil_severe_updraft
                      real :: thres_couplet_tourbillon(0:1, 0:1) !! seuil_couplet_tourbillon
                      character*10  :: region                                        !! region
                      character*200 :: ext_file_out;                !! ext_fichier_out
                      character*200 :: ext_file_in;                        !! ext_fichier_in
                      character*200 :: path_file_in;                !! path_fichier_in
                      character*200 :: path_file_out;                !! path_fichier_out
                      character*100 :: varname3d;                        !! nom_var3d
                      character*50  :: varname2d;                        !! nom_var2d
                   end type mystruct

                   real dx(ni, nj), dy(ni, nj), vort_2d(ni, nj), div_2d(ni, nj)
                   real intliq700(ni, nj)

                !!----------------------------------------------------
                !! defined const needed throughout the program
                !!----------------------------------------------------
                   integer, parameter :: debug = 1                  !! debug == 0 to not print any of the steps / debug ==1 to print every steps
                   integer, parameter :: valeur_bidon = -99999     ! flag value for skipping variable

                !!-------------------------------------------------------
                !! define 3d variable position in the 3d array called
                !! 3d_fields. the 3d_fields array has 3 dimensions
                !! (3d_fields[field_position][vertical_level][horizontal_grid_point])
                !!-------------------------------------------------------
                   integer, parameter :: gz = 0  !! geopotential height (dam)
                   integer, parameter :: tt = 1  !! air temperatutre (c)
                   integer, parameter :: hu = 2  !! specific humidity (kg/kg)
                   integer, parameter :: uv = 3  !! wind speed (knot)
                   integer, parameter :: ww = 4  !! vertical velocity (pa/s)
                   integer, parameter :: wd = 5  !! wind direction (meteorological degree)
                   integer, parameter :: qr = 6  !! relative vorticity (s^-1)
                   integer, parameter :: zet = 7  !! equivalent total reflectivity (dbz)
                   integer, parameter :: qjt1 = 8  !! graupel mixing ratio (kg/kg)
                   integer, parameter :: qit1 = 9  !! ice mixing ratio (kg/kg)
                   integer, parameter :: qht1 = 10 !! hail mixing ratio (kg/kg)
                   integer, parameter :: slw = 11 !! supercooled liquid water (kg/m^3)
                   integer, parameter :: dmh = 12 !! mean hail diameter (m)

                   integer, parameter :: nbr_var_3d = 13

                !!-------------------------------------------------------
                !! define 2d variable position in the 2d array called
                !! 2d_fields. the 3d_fields array has 2 dimensions
                !! (2d_fields[field_position][horizontal_grid_point])
                !!-------------------------------------------------------
                   integer, parameter :: sweat = 0 !! sweat index
                   integer, parameter :: lwc = 1 !! integrated liquid water content from the surface upward
                   integer, parameter :: lwc_700 = 2 !! integrated liquid water content from 700mb upward

                   integer, parameter :: nbr_var_2d = 3

                !!------------------------------------------------------------
                !! define the position of the severe weather variablesin the
                !! array called tv_var. the tv_var array has 2 dimensions
                !! (tv_var[sw_variable_position][horizontal_grid_point])
                !!-----------------------------------------------------------
                   integer, parameter :: hauteur_maxr_abovet0 = 0  !! height of the maximun recflectivity above the freezing level
                   integer, parameter :: maxr_abovet0 = 1  !! maximun recflectivity above the freezing level
                   integer, parameter :: couplet_tourbillon = 2  !! vorticity couplets
                   integer, parameter :: cisaillement_3km = 3  !! 0-3km wind shear
                   integer, parameter :: cisaillement_6km = 4  !! 0-6km wind shear
                   integer, parameter :: couplet_positif = 5  !! positive vorticity couplet
                   integer, parameter :: couplet_negatif = 6  !! negative vorticity couplet
                   integer, parameter :: zonesurplomb = 7  !! weak echo region
                   integer, parameter :: severe_updraft = 8  !! severe updraft
                   integer, parameter :: vertical_moisture_flux = 9  !! vertical moisture flux
                   integer, parameter :: mesocyclone = 10 !! mesocyclone
                   integer, parameter :: meso_base = 11 !! height of the mesocyclone's base
                   integer, parameter :: vef = 12 !! bounded weak echo region
                   integer, parameter :: tornade = 13 !! tornado
                   integer, parameter :: grosseur_grele_max = 14 !! maximum hail size in a column
                   integer, parameter :: grele_sfc = 15 !! hail size reaching the surface
                   integer, parameter :: supercooled_liquid_water = 16 !! integrated supercooled liquid water
                   integer, parameter :: vertical_ice_flux = 17 !! vertical ice flux
                   integer, parameter :: vertical_moisture_flux_new = 18 !! vertical moisture flux (new)

                   integer, parameter :: nbr_variabletv = 19

                !!-------------------------------------------
                !! defined variables needed in the program
                !!-------------------------------------------
!               integer, parameter :: nbr_niveau  = 10
!               integer, parameter :: ni          = 25
!                integer, parameter :: nj          = 25

                   integer :: i, j, k, m, ierr

                   real fields3d(0:(ni*nj) - 1, 0:nbr_niveau - 1, 0:nbr_var_3d - 1)
                   real fields3d_c(0:nbr_var_3d - 1, 0:nbr_niveau - 1, 0:(ni*nj) - 1)
                   real fields2d(0:(ni*nj) - 1, 0:nbr_var_2d - 1)
                   real fields2d_c(0:nbr_var_2d - 1, 0:(ni*nj) - 1)
                   real tv_var(0:nbr_variabletv - 1, 0:(ni*nj) - 1)
                   real probtv(0:(ni*nj) - 1)

                   character*6 file_name

                   type(mystruct):: paramconfig; 
                   file_name = "tv.cfg"

                !!--------------------------------
                !! initialize/populate the array
                !!--------------------------------
                   scaleht = 8000.
                   call get_grid_spacing_array(lat, lon, ni, nj, dx, dy)

!               3d fields
                   do k = 1, niveau

                      call get_grid_spacing_array(llat, llon, ni, nj, dx, dy)
                      call vortdiv(husig(:, :, k), hvsig(:, :, k), vort_2d, div_2d, ni, nj, dx, dy) ! qr

                      do j = 1, nj
                      do i = 1, ni
                         ij = (j - 1)*ni + (i - 1)
                         fields3d(ij, k - 1, gz) = hzsig(i, j, k)
                         fields3d(ij, k - 1, tt) = htsig(i, j, k)
                         fields3d(ij, k - 1, hu) = hmrsig(i, j, k) ! mixing ratio is used for now

                         speed = sqrt(husig(i, j, k)**2 + hvsig(i, j, k)**2)
                         fields3d(ij, k - 1, uv) = speed
                         fields3d(ij, k - 1, ww) = w_to_omegaf(hwsig(i, j, k), hwsig(i, j, k), scaleht)
                         if (speed .gt. 0.) then
                            fields3d(ij, k - 1, wd) = atan3d(-husig(i, j, k), -hvsig(i, j, k))
                         else
                            fields3d(ij, k - 1, wd) = 0.
                         end if

                         fields3d(ij, k - 1, qr) = vort_2d(i, j)
                         fields3d(ij, k - 1, zet) = hrefl_sig(i, j, k)
                         fields3d(ij, k - 1, qjt1) = hgraupelmr_sig(i, j, k)
                         fields3d(ij, k - 1, qit1) = 0.
                         fields3d(ij, k - 1, qht1) = 0.

                         if (htsig(i, j, k) .gt. 273.15) then
                            fields3d(ij, k - 1, slw) = 0.
                         else
                            fields3d(ij, k - 1, slw) = hcldliqmr_sig(i, j, k)
                         end if

                         fields3d(ij, k - 1, dmh) = valeur_bidon

                      end do ! i
                      end do ! j

                   end do ! k

                   fields2d = 0
                   do j = 1, nj
                   do i = 1, ni
                      ij = (j - 1)*ni + (i - 1)
                      fields2d(ij, sweat) = valeur_bidon
                      fields2d(ij, lwc) = intliqwater(i, j)
                      fields2d(ij, lwc_700) = intliq700(i, j)
                   end do ! i
                   end do ! j

                !!-------------------------------------------
                !! read the configuration file
                !!-------------------------------------------
                   ierr = read_config_file(%val(debug), file_name, paramconfig)

                !!-----------------------------------------
                !! re-arrange the multi-dimentional fortran
                !! to have the c array format
                !!-----------------------------------------
                   fields3d_c = reshape(fields3d, shape(fields3d_c), order=(/3, 2, 1/))

                   fields2d_c = transpose(fields2d)

             !!-------------------------------------------
                !! find the severe weather structural
                !! elements and the environmental ingredients
                !!-------------------------------------------
                   ierr = potentielle_temps_violent(%val(debug),%val(nbr_niveau), &
                                                    %val(ni),%val(nj), paramconfig,%val(nbr_var_3d), &
                                                    fields3d_c,%val(nbr_variabletv), tv_var)

                !!-------------------------------------------
                !! compute the severe thunderstorm intensity
                !! for each detected cell
                !!-------------------------------------------
                   ierr = calculprobtempsviolent(%val(debug), &
                                                 %val(ni),%val(nj), paramconfig,%val(nbr_var_2d), &
                                                 %val(nbr_variabletv), fields2d_c, tv_var, probtv)

                   return
                end

