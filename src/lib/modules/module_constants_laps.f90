module constants_laps

  ! constants used by laps     

      real,       parameter        :: cm2inch  =   1.0 / 2.54
      real,       parameter        :: m2inch   =   100. / 2.54
      real,       parameter        :: inch2m   =   .0254
      real,       parameter        :: cp       =   1005.7
      real,       parameter        :: grav     =      9.81
      integer,    parameter        :: imiss    = -99999
      integer,    parameter        :: ismth    =      0
      integer,    parameter        :: jsmth    =      0
      real,       parameter        :: lapse    =      6.5e-03
      real,       parameter        :: lv       =      2.5e+06
      real,       parameter        :: m2feet   =      3.281
      real,       parameter        :: mps2knts =      1.944
      real,       parameter        :: mps2mph  =      2.237
      real,       parameter        :: pi       =      3.1415927
      real,       parameter        :: p0       =   100000.0
      real,       parameter        :: r_gas    =    287.04
      real,       parameter        :: r        =    r_gas
      real,       parameter        :: rv       =    461.5
      real,       parameter        :: t0       =    273.15
      real,       parameter        :: xmiss    = -99999.9
      real,       parameter        :: cpog     = cp / grav
      real,       parameter        :: cpor     = cp / r
      real,       parameter        :: deg2rad  = pi / 180.0
      real,       parameter        :: e        = r / rv
      real,       parameter        :: gor      = grav / r
      real,       parameter        :: kappa    = r / cp
      real,       parameter        :: rad2deg  = 180.0 / pi
      real,       parameter        :: rog      = r / grav
      real,       parameter        :: rvolv    = rv / lv
      real,       parameter        :: rmissing = 1.e9
      integer,    parameter        :: imissing = nint(rmissing)
  
end module constants_laps


