module prmtrs_stmas

   ! laps parameters
   character*9          :: lapsast            ! laps ascii time
   integer              :: lapsi4t            ! laps i4time
   integer              :: icycle             ! laps cycle time
   integer              :: uniform            ! 1 uniform vertical grid

   real                 :: rmissing           ! laps default missing value (yuanfu)
   real                 :: badsfcdt           ! laps default bad surface (yuanfu)
   real, allocatable :: latitude(:, :)      ! laps latitude of the grid (yuanfu)
   real, allocatable :: longitud(:, :)      ! laps longitude of the grid (yuanfu)
   real, allocatable :: topogrph(:, :)      ! laps longitude of the grid (yuanfu)
   ! real    ,allocatable :: lapsradw(:)        ! laps radial wind (yuanfu)
   real, allocatable :: z_fcstgd(:)        ! vertical level of model or laps

   ! end of laps parameters.

   ! gps processing parameters:
   ! specific gas constant for water vapor:
   real, parameter :: rv_gas = 461
   ! coefficients of refractivity at microwave frequencies:
   real, parameter :: k2_rfc = 22.1
   real, parameter :: k3_rfc = 373900

   ! bufr parameters:
   real, parameter :: bufrmiss = 10.0e10 ! missing value in real for bufr data

   integer, parameter :: maxdims = 4          ! number of total dimension
   integer, parameter :: tmgobs_channel = 32  ! channel for tmg file
   integer, parameter :: pigobs_channel = 33  ! channel for tmg file
   integer              :: bgthobs            ! not use background until nobsgrid less than bgthobs
   integer              :: grdlevl            ! grid level
   integer              :: numgrid(maxdims)   ! number of grid for every dimension
   integer              :: inigrid(maxdims)   ! initial number of grid for every dimension
   integer              :: maxgrid(maxdims)   ! max number of grid for every dimension
   integer              :: ntmpgrd(maxdims)   ! number of old grid for every dimension
   integer              :: fcstgrd(maxdims)   ! number of old grid for every dimension
   integer              :: numdims            ! number of dimension valid
   integer              :: ngptobs            ! number of grid point per observation
   integer              :: numstat            ! number of state
   integer              :: nallobs            ! number of total observation valid
   integer, allocatable :: nobstat(:)         ! number of variety of observation valid
   integer              :: numvars            ! number of total control variables valid
   integer              :: ifbkgnd            ! whether use background or not, 1 is yes, 0 is not
   integer              :: ifbound            ! whether use bound or not, 1 is yes, 0 is not
   integer              :: ifpcdnt            ! whether use pure pressure coordinate, 1 is yes, 0 is not
   integer              :: if_test            ! whether run for test case, 1 is for the test case, added by zhongjie he
   integer              :: fnstgrd            ! the finest grid level
   integer              :: ifrepet            ! whethere run the mutigrid frame in monotonous or repeatedly, 0 for monotonous, 1 for repeatedly
   integer              :: itrepet            ! the times to repeat

   real, allocatable :: penal0x(:)         ! initial x direction penalty coefficent
   real, allocatable :: penal0y(:)         ! initial y direction penalty coefficent
   real, allocatable :: penal0z(:)         ! initial z direction penalty coefficent
   real, allocatable :: penal0t(:)         ! initial t direction penalty coefficent
   real, allocatable :: penal_x(:)         ! x direction penalty coefficent
   real, allocatable :: penal_y(:)         ! y direction penalty coefficent
   real, allocatable :: penal_z(:)         ! z direction penalty coefficent
   real, allocatable :: penal_t(:)         ! t direction penalty coefficent
   real                 :: pnlt0pu            ! initial geostrophic balance penalty coefficent, for p and u
   real                 :: pnlt0pv            ! initial geostrophic balance penalty coefficent, for p and v
   real                 :: pnlt_pu            ! geostrophic balance penalty coefficent, for p and u
   real                 :: pnlt_pv            ! geostrophic balance penalty coefficent, for p and v
   real                 :: grdspac(maxdims)   ! grid spacing
   real                 :: oripstn(maxdims)   ! original position of the study domain
   real                 :: costfun            ! cost function
   real, allocatable :: gradint(:, :, :, :, :) ! gradient of cost function
   character(len=100)   :: obsfile            ! observation file name

   integer              :: cosstep            ! the iterate steps before middle grid level
   integer              :: midgrid            ! middle grid level where the iterate step changed
   integer              :: finstep            ! the iterate steps after middle grid level

   integer              :: u_cmpnnt           ! u component index
   integer              :: v_cmpnnt           ! v component index
   integer              :: w_cmpnnt           ! w component index
   integer              :: pressure           ! pressure index
   integer              :: temprtur           ! temperature index
   integer              :: humidity           ! specific humidity index
!added by juxiang for cloud optical depth
   integer              :: cloudice           ! cloud ice index
   integer              :: cloudwat           ! cloud water index
   integer              :: raincont           ! rain content index
   integer              :: snowcont           ! snow content index
   integer              :: grapcont           ! graupel content index

!added by shuyuan 20100721 for radar reflectivity
   integer              :: rour_cmpnnt        ! density*rain water mixing ratio
   integer              :: rous_cmpnnt        ! density*snow water mixing ratio
!
   integer              :: xsl                ! x coordinate index
   integer              :: ysl                ! y coordinate index
   integer              :: psl                ! pressure index
   integer              :: csl                ! coriolis force index
   integer              :: dsl                ! density index

   ! integer ,allocatable :: idxradwn(:,:)      ! grid indices of radial wind

   real                 :: xytrans            ! coefficient used to translate the x and y coordinate to meters
   real                 :: z_trans            ! coefficient used to translate the z coordinate to meters

   real, allocatable :: www(:, :, :, :)       ! w component
   real, allocatable :: xxx(:, :)           ! distance of horizontal dimension 1
   real, allocatable :: yyy(:, :)           ! distance of horizontal dimension 2
   real, allocatable :: zzz(:, :, :, :)       ! height of general vertical coordinate
   real, allocatable :: ppp(:)             ! pressure vertical coordinate
   real, allocatable :: cor(:, :)           ! coriolis frequency
   real, allocatable :: deg(:, :)           ! direction angle
   real, allocatable :: den(:, :, :, :)       ! density
   real, allocatable :: scl(:)             ! scaling the physical variable
   real, allocatable :: sl0(:)             ! default scaling the physical variable
   real                 :: scp(5)             ! scaling xxx,yyy,zzz or ppp,cor,den
   real                 :: orivtcl
   real, allocatable :: xx0(:, :)           ! distance of horizontal dimension 1
   real, allocatable :: yy0(:, :)           ! distance of horizontal dimension 2
   real, allocatable :: zz0(:, :, :, :)       ! height of general vertical coordinate
   real, allocatable :: zzb(:)       ! sigma level for vertical coordinate
   real, allocatable :: pp0(:)             ! pressure vertical coordinate
   real, allocatable :: ppm(:)             ! multigrid pressure vertical coordinate
   real, allocatable :: cr0(:, :)           ! coriolis frequency
   real, allocatable :: dg0(:, :)           ! direction angle
   real, allocatable :: dn0(:, :, :, :)       ! density
! for observation
   real                 :: obsradar           ! the observation ratio of convention to radar data
   real                 :: obs_sfmr           ! the observation ratio of convention to sfmr data
   real                 :: obsref           ! the observation ratio of ref
   real, allocatable :: obspostn(:, :)      ! coordinate of observation point
   real, allocatable :: obscoeff(:, :)      ! interpolation coefficent
   integer, allocatable :: obsidxpc(:, :)      ! index of every node per cell
   integer, allocatable :: obsstate(:)        ! number of observation for every state
   real, allocatable :: obsvalue(:)        ! observation value
   real, allocatable :: obserror(:)        ! observation error
   real, allocatable :: obseinf1(:)        ! observation extra information 1
   real, allocatable :: obseinf2(:)        ! observation extra information 2
   real, allocatable :: obseinf3(:)        ! observation extra information 3

   ! gpsmet wetdelay data:
   integer, parameter :: max_gps = 5000                ! temporarily hard coded  !!changefrom 2000 to 5000 20100525 liu
   integer :: num_gps, i4t_gps                        ! i4t_gps: i4time of jan 1, 1970
   real    :: gps_xyt(3, max_gps), gps_elv(max_gps), gps_tim(max_gps), &
              gps_err(max_gps), gps_tpw(max_gps), gps_wet(max_gps)

   ! surface pressure:
   real, allocatable :: p_sfc_f(:, :)

! for grid point variable
   real, allocatable :: grdbkgd0(:, :, :, :, :)! grid background value
   real, allocatable :: grdbkgnd(:, :, :, :, :)! grid background value
   real, allocatable :: grdanals(:, :, :, :, :)! grid analysis value
   real, allocatable :: tmpanals(:, :, :, :, :)! grid analysis value

   real, allocatable :: wwwout(:, :, :, :)    ! w component for output
   real, allocatable :: difftout(:, :, :, :, :)! grid analysis value for outpur
   real, allocatable :: z_maxgid(:)        ! vertical level of finest analysis grid

   real                 :: pnlt0hy            ! initial hydrostatic condition penalty coefficent. addied by zhongjie he
   real                 :: pnlt_hy            ! hydrostatic condition penalty coefficent. addied by zhongjie he
   real                 :: taul_hy            ! reduction coefficent of hydrostatic condition penalty term, addied by zhongjie he
   real                 :: endhylv            ! the grid level after which the hydrostatic condition penalty term is ommitted. added by zhongjie he
   real                 :: endgslv            ! the grid level after which the geostrophic balance penalty term is ommitted. added by zhongjie he
   real                 :: limit_3            ! the limitation of height to decide in which range the observation is aviable. added by zhongjie he
   real                 :: limit_4            ! the limitation of time to decide in which range the observation is aviable. added by zhongjie he

   ! analysis domain by yuanfu xie:
   integer              :: itime2(2)          ! analysis time window

end module prmtrs_stmas
