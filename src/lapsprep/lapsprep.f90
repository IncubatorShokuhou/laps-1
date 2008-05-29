!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis
 
  PROGRAM lapsprep   
    !
    ! PURPOSE
    ! =======
    ! Prepares LAPS analysis data for ingest by various NWP model pre-processors.
    ! Currently supports MM5V3 (outputs PREGRID v3 format), WRF (outputs 
    ! gribprep format), and RAMS 4.3 (outputs RALPH2 format).   
    !
    ! ARGUMENTS
    ! =========
    !  LAPS_valid_time   - Optional command line argument of form YYJJJHHMM
    !                      specifying time for which to build output
    !                      If not present, the program will use the latest
    !                      available analysis based on LAPS systime.dat
    !
    ! REMARKS
    ! =======
    !  1. You must set the LAPS_DATA_ROOT environment variable before running.
    !  2. Other program controls in lapsprep.nl
    !
    ! HISTORY
    ! =======
    ! 28 Nov 2000 -- Original -- Brent Shaw 
    !    (based on lapsreader program originally developed by Dave Gill of
    !     NCAR to support MM5)

    ! Module declarations

    USE constants
    USE setup
    USE laps_static
    USE lapsprep_mm5
    USE lapsprep_wrf
    USE lapsprep_rams
    USE lapsprep_netcdf

    ! Variable Declarations

    IMPLICIT NONE

    ! Declarations for use of NetCDF library

    INCLUDE "netcdf.inc" 
    INTEGER :: cdfid , rcode
    INTEGER :: zid
    INTEGER :: z 
    INTEGER , DIMENSION(4) :: start , count
    INTEGER , DIMENSION(2) :: startc, countc
    INTEGER :: vid
    CHARACTER (LEN=132) :: dum 

    ! Arrays for data
    REAL , ALLOCATABLE , DIMENSION (:,:,:) :: u , v , t , rh , ht, &   
                                             lwc,rai,sno,pic,ice, sh, mr, w, & 
                                             virtual_t, rho,lcp
    REAL , ALLOCATABLE , DIMENSION (:,:)   :: slp , psfc, snocov, d2d,tskin
    REAL , ALLOCATABLE , DIMENSION (:)     :: p
    REAL , PARAMETER                       :: tiny = 1.0e-30
    
    ! Miscellaneous local variables
                                        
    INTEGER :: out_loop, loop , var_loop , i, j, k, kbot,istatus
    LOGICAL :: file_present
    REAL    :: rhmod, lwcmod, shmod, icemod
    REAL    :: rhadj
    REAL    :: lwc_limit
    REAL    :: hydrometeor_scale_pcp, hydrometeor_scale_cld

    ! Some stuff for JAX to handle lga problem
    ! with constant mr above 300 mb
    LOGICAL :: jaxsbn
    REAL  :: weight_top, weight_bot, newsh
    REAL, EXTERNAL ::make_rh 
    INTEGER :: k300
    jaxsbn = .false.
     
    ! Beginning of code

    CALL get_systime(i4time,laps_file_time,istatus)

    if (istatus.ne.1)then
     	write(*,*)'bad istatus from get_systime'
	stop
    endif

    PRINT *, 'LAPS_FILE_TIME = ', laps_file_time

    READ(laps_file_time, '(I2.2,I3.3,I2.2,I2.2)') valid_yyyy, valid_jjj, &
                                                   valid_hh, valid_min
    IF (valid_yyyy.LT.80) THEN
      valid_yyyy = 2000 + valid_yyyy
    ELSE
      valid_yyyy = 1900 + valid_yyyy
    ENDIF

    PRINT '(2A)', 'Running LAPSPREP using A9_time of: ', laps_file_time
  
    ! Get the LAPS_DATA_ROOT from the environment.  

    CALL GETENV('LAPS_DATA_ROOT', laps_data_root)
    PRINT '(2A)', 'LAPS_DATA_ROOT=',laps_data_root

    !  Get the namelist items (from the setup module).

    CALL read_namelist

    ! Get the static information (projection, dimensions,etc.)
 
    PRINT '(A)', 'Getting horizontal grid specs from static file.'
    CALL get_horiz_grid_spec(laps_data_root)

    ! Now that we have LAPS grid info, set up the hydrometeor scaling
    ! factor, which scales the concentrations of hydormeteors for this
    ! grid spacing.  We assume the values from LAPS are approprate on
    ! a grid with radar scaling (approx. 2km)

!beka    hydrometeor_scale = 2./dx  ! dx is in km
       hydrometeor_scale_pcp = hydrometeor_scale_factor_pcp/dx  ! dx is in km
       hydrometeor_scale_cld = hydrometeor_scale_factor_cld/dx  ! dx is in km

	print *,'beka',hydrometeor_scale_pcp, hydrometeor_scale_factor_pcp
	print *,'cj',hydrometeor_scale_cld, hydrometeor_scale_factor_cld

    !  Loop through each of the requested extensions for this date.  Each of the
    !  extensions has a couple of the variables that we want.

    PRINT '(A)', 'Starting Loop for each LAPS file'
    file_loop : DO loop = 1 , num_ext
      
      PRINT *, 'Looking for ',ext(loop)
      !  If this is a microphysical species but not doing 
      !  a hotstart, then cycle over this file.

      IF (((TRIM(ext(loop)).EQ.'lwc').OR.(TRIM(ext(loop)).EQ.'lcp')) .AND. &
          (.NOT.hotstart) ) THEN
        CYCLE file_loop
      ENDIF

      !  Build the input file name.   the input file.

      IF ((TRIM(ext(loop)) .NE. 'lw3' ).AND. &
          (TRIM(ext(loop)) .NE. 'lt1' ).AND. &
          (TRIM(ext(loop)) .NE. 'lq3' ).AND. &
          (TRIM(ext(loop)) .NE. 'lh3' )) THEN
        input_laps_file = TRIM(laps_data_root) //'/lapsprd/' // &
            TRIM(ext(loop)) // '/' // laps_file_time // '.' // &
            TRIM(ext(loop))
      ELSE
        IF (balance) THEN
          input_laps_file = TRIM(laps_data_root) //'/lapsprd/balance/' // &
          TRIM(ext(loop)) // '/' // laps_file_time // '.' // &
          TRIM(ext(loop))
        ELSE
          input_laps_file = TRIM(laps_data_root) //'/lapsprd/' // &
              TRIM(ext(loop)) // '/' // laps_file_time // '.' // &
              TRIM(ext(loop)) 
        ENDIF
      ENDIF
      PRINT *, 'Opening: ', input_laps_file

      ! Determine if the file exists

      INQUIRE (FILE=TRIM(input_laps_file), EXIST=file_present)
      IF (.NOT.file_present) THEN
        IF( (ext(loop).EQ.'lt1').OR. &
            (ext(loop).EQ.'lw3').OR. &
            (ext(loop).EQ.'lh3').OR. &
            (ext(loop).EQ.'lsx') ) THEN 
          PRINT '(A)', 'Mandatory file not available:' ,input_laps_file
          STOP 'not_enough_data'
        ELSE IF ( (ext(loop).EQ.'lq3') .OR. &
               (ext(loop).EQ.'lwc') ) THEN
          PRINT '(A)', 'File not available, cannot do hotstart.'
          hotstart = .false.
          CYCLE file_loop
        ELSE IF (ext(loop).EQ.'l1s') THEN
          PRINT '(A)', 'File not available, cannot do snowcover.'
          CYCLE file_loop
        ELSE
          PRINT '(A)', 'File not available, but not mandatory.'
          CYCLE file_loop
        ENDIF
      ENDIF

      ! Open the netcdf file and get the vertical dimension

      cdfid = NCOPN ( TRIM(input_laps_file) , NCNOWRIT , rcode )

      zid = NCDID ( cdfid , 'z' , rcode )
      CALL NCDINQ ( cdfid , zid , dum , z , rcode )

      IF ( ( ext(loop) .EQ. 'lsx' ) .OR. &
           ( ext(loop) .EQ. 'lm2') ) THEN
         z2 = z
      ELSE
         z3 = z
      END IF
      
      IF ( loop .EQ. 1 ) THEN

      ! ALLOCATE space for all of the variables in this data set.  Note
      ! that some ofthe 3d variables are allocated by z+1 instead of just z, as
      ! we are going to put the surface values into these arrays as well.

        ALLOCATE ( u   ( x , y , z3 + 1 ) )
        ALLOCATE ( v   ( x , y , z3 + 1 ) )
        ALLOCATE ( w   ( x , y , z3 + 1 ) )
        ALLOCATE ( t   ( x , y , z3 + 1 ) )
        ALLOCATE ( rh  ( x , y , z3 + 1 ) )
        ALLOCATE ( ht  ( x , y , z3 + 1 ) )
        ALLOCATE ( slp ( x , y         ) )
        ALLOCATE ( psfc (x , y         ) )
        ALLOCATE ( d2d ( x , y         ) )
        ALLOCATE ( tskin ( x , y       ) )
        ALLOCATE ( p   (         z3 + 1 ) ) 
        ! The following variables are not "mandatory"
        ALLOCATE ( lwc ( x , y , z3 ) ) 
        ALLOCATE ( rai ( x , y , z3 ) ) 
        ALLOCATE ( sno ( x , y , z3 ) ) 
        ALLOCATE ( pic ( x , y , z3 ) ) 
        ALLOCATE ( ice ( x , y , z3 ) )
        ALLOCATE ( snocov ( x , y ) ) 
        ALLOCATE ( lcp ( x , y , z3 ) ) 
        ! The following variables are only used
        ! for converting non-mandatory cloud variables
        ! to mixing ratio values 
        ALLOCATE ( rho ( x , y , z3 ) )
        ALLOCATE ( virtual_t ( x , y , z3 ) )
        ALLOCATE ( sh ( x , y , z3 ) )
        ALLOCATE ( mr ( x , y , z3+1 ) )

        ! Initialize the non-mandatory values

        lwc(:,:,:) = -999.
        rai(:,:,:) = -999.
        sno(:,:,:) = -999.
        pic(:,:,:) = -999.
        snocov(:,:) = -999.
        lcp(:,:,:) = 0.
      END IF

      IF       ( ext(loop) .EQ. 'lh3' ) THEN

        ! Loop over the number of variables for this data file.

        var_lh3 : DO var_loop = 1 , num_cdf_var(loop)

          ! Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          CALL NCVGT ( cdfid , vid , start , count , rh , rcode )

          !  Do this just once for pressure.

          vid = NCVID ( cdfid , 'level' , rcode )

          CALL NCVGT ( cdfid , vid , 1 , z3 , p , rcode )
 
          ! Set the pressure level of the lowest level of our
          ! pressure array as 2001 mb to flag the surface
          p(z3+1) = 2001

        END DO var_lh3 

      ELSE IF ( ext(loop) .EQ. 'lq3' )THEN

        !  Loop over the number of variables for this data file.
        var_lq3 : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          CALL NCVGT ( cdfid , vid , start , count , sh , rcode )

        END DO var_lq3

      ELSE IF ( ext(loop) .EQ. 'lsx' ) THEN

        !  Loop over the number of variables for this data file.

        var_lsx : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , 1 , 1 /)

          IF      ( cdf_var_name(var_loop,loop) .EQ. 'u  ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , u  (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'v  ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , v  (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'vv ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , w  (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 't  ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , t  (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'rh ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , rh (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'mr ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , mr (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'msl' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , slp           , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'ps ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , psfc          , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'tgd') THEN
            CALL NCVGT ( cdfid , vid , start , count , tskin         , rcode )
          END IF

        END DO var_lsx

        ! Convert sfc mixing ratio from g/kg to kg/kg.

        mr(:,:,z3+1)=mr(:,:,z3+1)*0.001

      ELSE IF ( ext(loop) .EQ. 'lm2' ) THEN

        var_l1s : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , 1 , 1 /)

          IF      ( cdf_var_name(var_loop,loop) .EQ. 'sc ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , snocov        , rcode )
          END IF

        END DO var_l1s                             
        
      ELSE IF ( ext(loop) .EQ. 'lt1' ) THEN

        !  Loop over the number of variables for this data file.

        var_lt1 : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          IF      ( cdf_var_name(var_loop,loop) .EQ. 't3 ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , t  , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'ht ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , ht , rcode )
          END IF

        END DO var_lt1

      ELSE IF ( ext(loop) .EQ. 'lw3' ) THEN

        !  Loop over the number of variables for this data file.

        var_lw3 : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          IF      ( cdf_var_name(var_loop,loop) .EQ. 'u3 ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , u , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'v3 ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , v , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'om ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , w , rcode )
          END IF

        END DO var_lw3

      ELSE IF (( ext(loop) .EQ. 'lwc' ).AND.(hotstart)) THEN

        !  Loop over the number of variables for this data file.

        var_lwc1 : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          IF      ( cdf_var_name(var_loop,loop) .EQ. 'lwc' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , lwc, rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'rai' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , rai , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'sno' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , sno , rcode ) 
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'ice' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , ice , rcode ) 
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'pic' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , pic , rcode ) 
          END IF

        END DO var_lwc1    

      ELSE IF (( ext(loop) .EQ. 'lcp' ).AND.(hotstart)) THEN

        !  Loop over the number of variables for this data file.

        var_lvc : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          IF      ( cdf_var_name(var_loop,loop) .EQ. 'lcp' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , lcp, rcode )
            print *, 'Got cloud cover...min/max = ',minval(lcp),maxval(lcp)
          END IF

        END DO var_lvc  

      END IF

    END DO file_loop

    ! Compute mixing ratio from spec hum.
    ! Fill missing values with sfc value.

    k300 = 0
    do k= 1,z3
      if (p(k) .eq. 300.) k300 = k
    enddo
    if (k300 .eq. 0) THEN
      print *, "Could not find k300!"
      stop
    endif
    do k=1,z3
    do j=1,y
    do i=1,x

      if ((jaxsbn).and.(p(k).LT.300.)) then
           weight_bot = (p(k) - 50) / (250)
           weight_top = 1.0 - weight_bot
           newsh = weight_bot * sh(i,j,k300) + &
                       weight_top * tiny 
               
           newsh = MIN(sh(i,j,k),newsh)
            
           ! Make sure sh does not exceed 
           ! ice saturation value
           CALL saturate_ice_points(t(i,j,k), &
                                    p(k),1.0, &
                                    shmod,rhmod)
           sh(i,j,k) = MIN(shmod, newsh)
           rh(i,j,k) = make_rh(p(k),t(i,j,k)-273.15, &
             sh(i,j,k)*1000., -132.) * 100.
      endif
     
      if (sh(i,j,k) .ge. 0. .and. sh(i,j,k) .lt. 1.) then
        mr(i,j,k)=sh(i,j,k)/(1.-sh(i,j,k))
      else
        mr(i,j,k)=mr(i,j,z3+1)
      endif
!      if (u(i,j,k) .eq. 1.e-30 .or. abs(u(i,j,k)) .gt. 200.) u(i,j,k)=u(i,j,z3+1)
!     if (v(i,j,k) .eq. 1.e-30 .or. abs(v(i,j,k)) .gt. 200.) v(i,j,k)=v(i,j,z3+1)
    enddo
    enddo
    enddo

    !  Set the lowest level of the geopotential height to topographic height

    ht(:,:,z3+1) = topo 

    IF (hotstart) THEN

      ! If this is a hot start, then we need to convert the microphysical
      ! species from mass per volume to mass per mass (mixing ratio).  This
      ! requires that we compute the air density from virtual temperature
      ! and divide each species by the air density.
       
      ! Compute virtual temperature from mixing ratio and temperature
      virtual_t(:,:,:)=( 1. + 0.61*mr(:,:,1:z3))*t(:,:,1:z3)
 
      ! Compute density from virtual temperature and gas constant for dry air
      DO k = 1, z3
        rho(:,:,k) = p(k)*100. / (rdry * virtual_t(:,:,k))
      ENDDO

      ! For each of the species, ensure they are not "missing".  If missing
      ! then set their values to 0.000.  Otherwise, divide by the density to 
      ! convert from concentration to mixing ratio.

      IF (MAXVAL(lwc) .LT. 99999.) THEN

        ! Scale lwc for grid spacing
        lwc = lwc*hydrometeor_scale_cld
        ! Cap lwc to autoconversion rate for liquid to rain
        WHERE(lwc .GT. autoconv_lwc2rai) lwc = autoconv_lwc2rai
        ! Convert lwc concentration to mixing ratio
        lwc(:,:,:) = lwc(:,:,:)/rho(:,:,:)   ! Cloud liquid mixing ratio

        ! Convert lwc mixing ratio to vapor mixing ratio
        IF (lwc2vapor_thresh .GT. 0.) THEN
          DO k=1,z3
            DO j=1,y
              DO i=1,x  
                IF( (lcp(i,j,k).GE.lcp_min).AND.&
                    (t(i,j,k).GE.263.0).AND.&
                    (lwc(i,j,k).GT.lwc_min))THEN  
                !IF (lwc(i,j,k).GT.0.00010) THEN
                  !CALL lwc2vapor(lwc(i,j,k),sh(i,j,k),t(i,j,k), &
                  !               p(k),lwc2vapor_thresh, &
                  !               lwcmod,shmod,rhmod)

                  CALL saturate_lwc_points(sh(i,j,k),t(i,j,k), &
                                           p(k),lwc2vapor_thresh, &
                                           shmod,rhmod)
                  ! Update moisture arrays

                  sh(i,j,k) = shmod
                  rh(i,j,k) = rhmod
                  mr(i,j,k) = shmod/(1.-shmod)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSE
        PRINT *,'Missing cloud liquid, setting values to 0.0'
        lwc(:,:,:) = 0.0
      ENDIF

      IF (MAXVAL(rai) .LT. 99999.) THEN
        rai(:,:,:) = rai(:,:,:) * hydrometeor_scale_pcp
        rai(:,:,:) = rai(:,:,:)/rho(:,:,:)   ! Rain mixing ratio
      ELSE
        PRINT *, 'Missing rain, setting values to 0.0' 
        rai(:,:,:) = 0.0
      ENDIF

      IF (MAXVAL(sno) .LT. 99999.) THEN
        sno(:,:,:) = sno(:,:,:) * hydrometeor_scale_pcp
        sno(:,:,:) = sno(:,:,:)/rho(:,:,:)   ! Snow mixing ratio
      ELSE
        PRINT *, 'Missing snow, setting values to 0.0'    
        sno(:,:,:) = 0.0
      ENDIF

      IF (MAXVAL(ice) .LT. 99999.) THEN 
        ! Limit ice to autoconversion threshold

        ice(:,:,:) = ice(:,:,:) * hydrometeor_scale_cld
        WHERE(ice .GT. autoconv_ice2sno) ice = autoconv_ice2sno
        ice(:,:,:) = ice(:,:,:)/rho(:,:,:)   ! Ice mixing ratio
      
        ! Convert ice mixing ratio to vapor mixing ratio
        ! (for now ice conversion tied to lwc2vapor_thresh)

        IF (lwc2vapor_thresh .GT. 0.) THEN
          DO k=1,z3
            DO j=1,y
              DO i=1,x
                IF ((lcp(i,j,k).GE.lcp_min).AND. &
                    (t(i,j,k).LT.263.).AND. &
                    (ice(i,j,k).GT.ice_min)) THEN  

                      ! Update moisture arrays
                       CALL saturate_ice_points(t(i,j,k), &
                                           p(k),lwc2vapor_thresh, &
                                           shmod,rhmod)
                  sh(i,j,k) = shmod
                  rh(i,j,k) = rhmod
                  mr(i,j,k) = sh(i,j,k)/(1.-sh(i,j,k))
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSE
        PRINT *, 'Missing ice, setting values to 0.0' 
        ice(:,:,:) = 0.0
      ENDIF

      IF (MAXVAL(pic) .LT. 99999.) THEN
         pic(:,:,:) = pic(:,:,:)*hydrometeor_scale_pcp
        pic(:,:,:) = pic(:,:,:)/rho(:,:,:)   ! Graupel (precipitating ice) mixing rat.
      ELSE
        PRINT *, 'Missing pice, setting values to 0.0' 
        pic(:,:,:) = 0.0
      ENDIF

      ! Convert 3d omega from Pa/s to m/s, or fill with sfc value if missing.

      do k=1,z3
      do j=1,y
      do i=1,x
        if (w(i,j,k) .eq. 1.e-30 .or. abs(w(i,j,k)) .gt. 100.) then
          w(i,j,k)=w(i,j,z3+1)
        else
          w(i,j,k)=-w(i,j,k)/(rho(i,j,k)*g)
        endif
      enddo
      enddo
      enddo


    ENDIF

    ! If make_sfc_uv set, then replace surface winds with 
    ! winds interpolated from the 3D field.
    IF (make_sfc_uv) THEN
      PRINT *, 'Creating surface u/v from 3D field...'
      DO j = 1,y
        DO i = 1,x
          kbot = 0
          get_lowest: DO k = z3,1,-1
            IF (ht(i,j,k) .GT. topo(i,j)) THEN
              kbot = k
              EXIT get_lowest
            ENDIF
          ENDDO get_lowest
          IF (kbot .NE. 0) THEN
            u(i,j,z3+1) = u(i,j,kbot)
            v(i,j,z3+1) = v(i,j,kbot)
          ELSE
            print *, 'Problem finding kbot.'
            STOP
          ENDIF
        ENDDO
      ENDDO
  
    ENDIF

    ! Loop over the each desired output format

    DO out_loop = 1, num_output
      ! Now it is time to output these arrays.  The arrays are ordered
      !  as (x,y,z).  The origin is the southwest corner at the top of the 
      ! atmosphere for the 3d arrays, where the last layer (z3+1) contains    
      ! the surface information.  This is where you would insert a call
      ! to a custom output routine.

      select_output: SELECT CASE (output_format(out_loop))
        CASE ('mm5 ')
          CALL output_pregrid_format(p, t, ht, u, v, rh, slp, &
                            lwc, rai, sno, ice, pic,snocov, tskin)

        CASE ('wrf ')
          CALL output_gribprep_format(p, t, ht, u, v, rh, slp, psfc,&
                             lwc, rai, sno, ice, pic,snocov, tskin)
        CASE ('wps ')
          CALL output_metgrid_format(p, t, ht, u, v, rh, slp, psfc,&
                             lwc, rai, sno, ice, pic,snocov, tskin)
     
        CASE ('rams') 
          CALL output_ralph2_format(p,u,v,t,ht,rh,slp,psfc,snocov,tskin)
        CASE ('sfm ')
          PRINT '(A)', 'Support for SFM (RAMS 3b) coming soon...check back later!'

        CASE ('cdf ')
          CALL output_netcdf_format(p,ht,t,mr,u,v,w,slp,psfc,lwc,ice,rai,sno,pic)

        CASE DEFAULT
          PRINT '(2A)', 'Unrecognized output format: ', output_format
          PRINT '(A)', 'Recognized formats include mm5, rams, wrf, sfm, and cdf'

      END SELECT select_output
    ENDDO 
    PRINT '(A)', 'LAPSPREP Complete.'

  END program lapsprep
  


