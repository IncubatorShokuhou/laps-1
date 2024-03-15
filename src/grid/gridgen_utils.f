      subroutine compute_latlon(nnxp,nnyp,n_staggers,mdlat,mdlon
     +,deltax,xtn,ytn,lats,lons,istatus)


      implicit none

      integer nnxp,nnyp
      integer n_staggers
      integer i,j,k,nc

      real    deltax,deltay
      real    mdlat,mdlon
      real    stagger_ew,stagger_ns
      real    erad

c these returned to gridgen
      real    xtn(nnxp,n_staggers)
      real    ytn(nnyp,n_staggers)

c these are the anchor points (sw corner) of the grid
      real    xmn1,ymn1

c     real    lat(nnxp,nnyp)
c     real    lon(nnxp,nnyp)

c a-c staggers contained within these arrays.
      real    lats(nnxp,nnyp,n_staggers)
      real    lons(nnxp,nnyp,n_staggers)

      character c_dataroot*255
      character c10_grid_fname*10

      integer istatus,itstatus,ishow_timer

      print*,'calculate lat/lon at stagger grid points.'

      deltay=deltax

c     call get_grid_center(mdlat,mdlon,istatus)
c     if(istatus .ne. 1)then
c        write(6,*)' error returned: get_grid_center'
c        return
c     endif

      call get_earth_radius(erad,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error calling get_earth_radius'
         return
      endif

      do k=1,n_staggers

         if(k.eq.1)then

c get x/y for lower left corner
            call polar_gp(mdlat,mdlon,xmn1,ymn1,deltax,deltay,
     +  nnxp,nnyp,1.)
            stagger_ns=0
            stagger_ew=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn(1,k),ytn(1,k))
 
         elseif(k.eq.2)then    !this is a-stagger (.5 e-w stagger) 

            stagger_ew=0.5*deltax
            stagger_ns=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn(1,k),ytn(1,k))

         elseif(k.eq.3)then !this is the b-stagger (.5 n-s stagger)

            call get_xytn(deltax,deltay,nnxp,nnyp,0,0,xmn1,ymn1
     +,xtn(1,k),ytn(1,k))
            stagger_ew=0
            stagger_ns=0.5*deltay 
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn(1,k),ytn(1,k)) 

         elseif(k.eq.4)then !this is the c-stagger (.5 both n-s and e-w)

            call get_xytn(deltax,deltay,nnxp,nnyp,0,0,xmn1,ymn1
     +,xtn(1,k),ytn(1,k))
            stagger_ew=0.5*deltax
            stagger_ns=0.5*deltay
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn(1,k),ytn(1,k))

         endif
c
c*****************************************************************
c*  convert it to lat/lon using the library routines.            *
         itstatus=ishow_timer()

         do j = 1,nnyp
           do i = 1,nnxp

             call xy_to_latlon(xtn(i,k),ytn(j,k),erad ! ,90.,std_lon     
     1                                     ,lats(i,j,k),lons(i,j,k))

c            print *,'i,j,xtn,ytn,pla,lplo=',i,j,xtn,ytn,pla,plo
           enddo
           if(j .eq. (j/50)*50)then
             print*,'completed row ',j
           endif
         enddo

      enddo

      print*,'completed compute_latlon routine'
 
      itstatus=ishow_timer()

      istatus = 1

      return
      end

c ----------------------------------------------

      subroutine get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)

      implicit none

      integer  nnxp,nnyp
      integer  i,j

      real    deltax,deltay
      real    stagger_ew,stagger_ns
      real    xmn1,ymn1
      real    xmn(nnxp)
      real    ymn(nnyp)
      real    xtn(nnxp)
      real    ytn(nnyp)
      real    ridiff,rjdiff

      xmn(1)=xmn1+stagger_ew
      ymn(1)=ymn1+stagger_ns

      do 600 i=2,nnxp
         ridiff = float(i-1)
         xmn(i)=xmn(1)+ridiff*deltax
 600  continue
      xmn(nnxp)=2*xmn(nnxp-1)-xmn(nnxp-2)
 
      do 610 j=2,nnyp
         rjdiff = float(j-1)
         ymn(j)=ymn(1)+rjdiff*deltay
 610  continue
      ymn(nnyp)=2*ymn(nnyp-1)-ymn(nnyp-2)
 
      do 650 i=2,nnxp
         xtn(i)=.5*(xmn(i)+xmn(i-1))
 650  continue
      xtn(1)=1.5*xmn(1)-.5*xmn(2)

      do 660 j=2,nnyp
         ytn(j)=.5*(ymn(j)+ymn(j-1))
 660  continue
      ytn(1)=1.5*ymn(1)-.5*ymn(2)

      write(6,*)' results at end of get_xytn:'
      write(6,*)' deltax/deltay = ',deltax,deltay
      write(6,*)' xmn bounds = ',xmn(1),xmn(nnxp)
      write(6,*)' ymn bounds = ',ymn(1),ymn(nnyp)
      write(6,*)' xtn bounds = ',xtn(1),xtn(nnxp)
      write(6,*)' ytn bounds = ',ytn(1),ytn(nnyp)

      return
      end
c
c -------------------------------------------------------------
c
      subroutine get_map_factor_grid(nx,ny,n_staggers
     +,rlats,rlons ,rmap_factors,istatus)

      implicit none

      integer nx,ny
      integer n_staggers
      integer i,j,k
      integer istatus

      real    rlats(nx,ny,n_staggers)
      real    rlons(nx,ny,n_staggers)
      real    rmap_factors(nx,ny,n_staggers)
      real    sigma

      do k=1,n_staggers
      do j=1,ny
      do i=1,nx
         call get_sigma(rlats(i,j,k),rlons(i,j,k)
     +,sigma,istatus)
         if(istatus.ne.1)goto 999
         rmap_factors(i,j,k)=sigma
      enddo
      enddo
      enddo

      return

999   print*,'error returned: get_sigma'

      return
      end
c
c ------------------------------------------------------------
c
      subroutine get_coriolis_components(nx,ny,lat,coriolis_parms)
c
      include 'trigd.inc'
      implicit none

      integer nx,ny
      integer i,j

      real    lat(nx,ny)
      real    coriolis_parms(nx,ny,2)

      real    omega_ear
      data    omega_ear/7.292e-5/

      do j=1,ny
      do i=1,nx
         coriolis_parms(i,j,1)=2*omega_ear*sind(lat(i,j))
         coriolis_parms(i,j,2)=2*omega_ear*cosd(lat(i,j))
      enddo
      enddo

      return
      end
c
c ----------------------------------------------------------
c
      subroutine get_projrot_grid(nx,ny,lat,lon,projrot_grid
     +,istatus)

      include 'trigd.inc'
      implicit none

      integer nx,ny
      integer istatus
      integer i,j
      real    r_missing_data
      real    projrot_deg
      real    projrot_latlon
      real    projrot_grid(nx,ny,2)
      real    lat(nx,ny),lon(nx,ny)

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'error returned: get_r_missing_data'
         return
      endif

      do j=1,ny
      do i=1,nx
         projrot_grid(i,j,1)=r_missing_data
         projrot_grid(i,j,2)=r_missing_data
         projrot_deg=projrot_latlon(lat(i,j),lon(i,j),istatus)
         if(istatus.eq.1)then
            projrot_grid(i,j,1)=sind(projrot_deg)
            projrot_grid(i,j,2)=cosd(projrot_deg)
         endif
      enddo
      enddo

      return
      end
c
c--------------------------------------------------------
c
      subroutine get_static_albedo(nx,ny,lat,lon,landfrac
     +,static_albedo,istatus)

      implicit none

      integer nx,ny
      real    lat(nx,ny),lon(nx,ny)
      real    landfrac(nx,ny)
      real    static_albedo(nx,ny)
      real    water_albedo_cmn
      real    r_missing_data
      integer istatus
      integer i,j
      integer nwater

      data water_albedo_cmn/0.08/

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'error returned: get_r_missing_data'
         return
      endif

      nwater = 0
      do j=1,ny
      do i=1,nx
         if(landfrac(i,j).le.0.01)then
            nwater = nwater+1
            static_albedo(i,j)=water_albedo_cmn
         else
            static_albedo(i,j)=r_missing_data
         endif
      enddo
      enddo
      if(nwater .gt. 0)then
         print*,'static_albedo ',water_albedo_cmn,' used ',
     +'at ',nwater,' grid points '
      else
         print*,'no water grid points for water_albedo',
     +' in this domain'
      endif

      return
      end
c
c --------------------------------------------------------
c
      subroutine bilinear_interp(i,j,imax,jmax,array_2d,result)
c
c this used only for getting topography on the c-staggered grid

      implicit none
      integer i,j,imax,jmax
      real result
      real array_2d(imax,jmax)
      real z1,z2,z3,z4
      real fraci,fracj

      fraci = 0.5
      fracj = 0.5

      z1=array_2d(i  , j  )
      z2=array_2d(i-1, j  )
      z3=array_2d(i-1, j-1)
      z4=array_2d(i  , j-1)

      result= z1+(z2-z1)*fraci+(z4-z1)*fracj
     1     - (z2+z4-z3-z1)*fraci*fracj

      return
      end

c
c--------------------------------------------------------
c
      subroutine get_gridgen_var(nf,ngrids,var,comment)

      implicit none

      integer        nf,ngrids,i,j
      character*(*)  var(nf)
      character*(*)  comment(nf)
      character*2    cat

      write(6,*)' get_gridgen_var: ngrids = ',ngrids

      if(ngrids.eq.38)then  !this is the laps analysis section

         var(1)    = 'lat'
         var(2)    = 'lon'
         var(3)    = 'avg'
         var(4)    = 'ldf'
         var(5)    = 'use'
         var(6)    = 'alb'  !now used for max snow alb 2-20-03 js.
         var(7)    = 'std'
         var(8)    = 'sln'
         var(9)    = 'slt'
         var(10)   = 'stl'
         var(11)   = 'sbl'
         var(12)   = 'lnd'  ! land-water mask based on usgs landuse
         i=12
         do j=1,12
            write(cat,'(i2.2)')j
            if(cat(1:1).eq.' ')cat(1:1)='0'
            if(cat(2:2).eq.' ')cat(2:2)='0'
            var(i+j)= 'g'//cat   ! vegetation greenness fraction
         enddo

         var(25)='tmp'
         i=25
         do j=1,12
            write(cat,'(i2.2)')j
            if(cat(1:1).eq.' ')cat(1:1)='0'
            if(cat(2:2).eq.' ')cat(2:2)='0'
            var(i+j)= 'a'//cat   ! monthly albedo
         enddo

         var(ngrids)   = 'zin'
 
         comment(1) = 'lat: from model by j. smart/ s. albers 2-03\0'
         comment(2) = 'lon: from model by j. smart/ s. albers 2-03\0'
         comment(3) = 'average terrain elevation (m) \0'
         comment(4) = 'land fraction: derived from usgs land use \0'
         comment(5) = 'land use dominant category (usgs 24 category) \0'
         comment(6) = 'maximum snow albedo; defined over land only \0'
         comment(7) = 'standard deviation of elevation data (m)\0'
         comment(8) = 'mean longitudinal terrain slope (m/m)\0'
         comment(9) = 'mean latitudinal terrain slope (m/m)\0'
         comment(10)= 'top layer (0-30cm) dominant category soiltype\0'
         comment(11)= 'bot layer (30-90cm) dominant category soiltype\0'
         comment(12)= 'land-water mask (0=water; 1 otherwise) \0'

         i=12
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= 'vegetation greenness fraction: mon = '//cat
         enddo

         comment(25)='mean annual soil temp (deg k)'
         i=25
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= 'climatological albedo: mon = '//cat
         enddo


         comment(ngrids)='\0'

      elseif(ngrids.eq.112)then   !this is the wrfsi section

         var(1)    = 'lat'  ! non-staggered (analysis-grid) lats
         var(2)    = 'lon'  ! non-staggered (analysis-grid) lons
         var(3)    = 'laa'  ! a-stagger (.5*deltax (e-w)) lats
         var(4)    = 'loa'  ! a-stagger (.5*deltax (e-w)) lons
         var(5)    = 'lab'  ! b-stagger (.5*deltay (n-s)) lats
         var(6)    = 'lob'  ! b-stagger (.5*deltay (n-s)) lons
         var(7)    = 'lac'  ! c-stagger (.5*deltax and .5*deltay) lats
         var(8)    = 'loc'  ! c-stagger (.5*deltay and .5*deltay) lons

         var(9)    = 'avg'  ! topo (m) on analysis-grid

         var(10)   = 'ldf'  ! land fraction
         var(11)   = 'use'  ! landuse dominant category
         var(12)   = 'lnd'  ! land-water mask based on usgs landuse
         var(13)   = 'stl'  ! soiltype top layer dominant category
         var(14)   = 'sbl'  ! soiltype bot layer dominant category

         i=14
         do j=1,24
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'u'//cat
         enddo

         i=39
         var(i)     = 'spr'  ! sin(projection rotation) from true
         var(i+1)   = 'cpr'  ! cos(projection rotation) from true
         var(i+2)   = 'mfl'  ! map factor analysis grid
         var(i+3)   = 'mfa'  ! map factor a-stagger grid
         var(i+4)   = 'mfb'  ! map factor b-stagger grid
         var(i+5)   = 'mfc'  ! map factor c-stagger grid
         var(i+6)   = 'cph'  ! horizontal component of coriolis parameter
         var(i+7)   = 'cpv'  ! vertical component of coriolis parameter

         var(i+8)   = 'alb'  ! maximum snow albedo defined over land only
         var(i+9)   = 'std'  ! standard deviation of elevation data (m)
         var(i+10)  = 'sln'  ! terrain slope; longitudinal component (m/m)
         var(i+11)  = 'slt'  ! terrain slope; latitudinal component (m/m)
         var(i+12)  = 'avc'  ! topo (m) on c-stagger grid

         i=51
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 't'//cat   !top layer (0-30cm) soiltype (% dist)
         enddo
         i=67
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'b'//cat   ! bot layer (30-90cm) soiltype (% dist)
         enddo
         i=83
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'g'//cat   ! vegetation greenness fraction
         enddo

         var(96)='tmp'

         i=96
         do j=1,12
            write(cat,'(i2.2)')j
            if(cat(1:1).eq.' ')cat(1:1)='0'
            if(cat(2:2).eq.' ')cat(2:2)='0'
            var(i+j)= 'a'//cat ! monthly albedo
         enddo

         var(109) = 'slp'      !terrain slope index, dominant category
         var(110) = 'gnx'      !max greenness fraction
         var(111) = 'gnn'      !min greenness fraction

         var(ngrids)   = 'zin'

         comment(1) = 'made from model by j. smart/ s. albers 2-03\0'
         comment(2) = 'made from model by j. smart/ s. albers 2-03\0'
         comment(3) = 'a-stagger grid latitudes \0'
         comment(4) = 'a-stagger grid longitudes for wrf_si \0'
         comment(5) = 'b-stagger grid latitudes for wrf_si \0'
         comment(6) = 'b-stagger grid longitudes for wrf_si \0'
         comment(7) = 'c-stagger grid latitudes for wrf_si \0'
         comment(8) = 'c-stagger grid longitudes for wrf_si \0'
         comment(9) = 'average terrain elevation (m) \0'
         comment(10)= 'land fraction computed from usgs landuse\0'
         comment(11)= 'land use dominant category \0'
         comment(12)= 'land-water mask (0=water; 1 otherwise) \0'
         comment(13)= 'soiltype top layer dominant category \0'
         comment(14)= 'soiltype bot layer dominant category \0'

         i=14
         do j=1,24
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'u'//cat
          comment(i+j)= '% dist land use category '//cat//' \0'
         enddo

         i=39
         comment(i)= 'sin of projection rotation (rad) \0'
         comment(i+1)= 'cosine of projection rotation (rad) \0'
         comment(i+2)= 'map factor analysis grid \0'
         comment(i+3)= 'map factor a-stagger grid \0'
         comment(i+4)='map factor b-stagger grid \0'
         comment(i+5)= 'map factor c-stagger grid \0'
         comment(i+6)= 'horizontal component coriolis parameter \0'
         comment(i+7)= 'vertical component coriolis parameter \0'

         comment(i+8)= 'maximum snow albedo (%) over land only \0'
         comment(i+9)= 'standard deviation of elevation data (m)\0'
         comment(i+10)= 'mean longitudinal terrain slope (m/m)\0'
         comment(i+11)= 'mean latitudinal terrain slope (m/m)\0'
         comment(i+12)= 'average terrain elevation (c-stagger) (m)\0'

         i=51
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)='% dist top layer soiltype category '//cat
         enddo

         i=67
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= '% dist bot layer soiltype category '//cat
         enddo

         i=83
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= 'vegetation greenness fraction: mon = '//cat
         enddo

         comment(96)='1 degree mean annual soiltemp (deg k)'
         i=96
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= 'climatological albedo: mon = '//cat
         enddo

         comment(109)=  'terrain slope index'
         comment(110)=  'annual max greenness fraction'
         comment(111)=  'annual min greenness fraction'

         comment(ngrids)= '\0'

	elseif (ngrids .eq. 103) then ! rotlat case

         var(1)   = 'lah'  ! mass(h)-point lats
         var(2)   = 'loh'  ! mass(h)-point lons
         var(3)   = 'lav'  ! wind(v)-point lats
         var(4)   = 'lov'  ! wind(v)-point lons
         var(5)   = 'ldf'  ! land fraction
         var(6)   = 'use'  ! landuse dominant category
         var(7)   = 'lnd'  ! land-water mask based on usgs landuse
         var(8)   = 'stl'  ! soiltype top layer dominant category
         var(9)   = 'sbl'  ! soiltype bot layer dominant category
     
         i=10
                                                                        
         var(i)     = 'spr'  ! sin(projection rotation) from true
         var(i+1)   = 'cpr'  ! cos(projection rotation) from true
         var(i+2)   = 'cph'  ! horizontal component of coriolis parameter
         var(i+3)   = 'cpv'  ! vertical component of coriolis parameter
         var(i+4)   = 'alb'  ! maximum snow albedo
         var(i+5)   = 'std'  ! standard deviation of elevation data (m)
         var(i+6)  = 'sln'  ! terrain slope; longitudinal component (m/m)
         var(i+7)  = 'slt'  ! terrain slope; latitudinal component (m/m)
         var(i+8)  = 'avc'  ! topo (m) on mass points 

        i=18
         do j=1,24
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'u'//cat
         enddo
                                                                         
         i=42
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 't'//cat   !top layer (0-30cm) soiltype (% dist)
         enddo
         i=58
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'b'//cat   ! bot layer (30-90cm) soiltype (% dist)
         enddo
         i=74
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'g'//cat   ! vegetation greenness fraction
         enddo
                                                                          
         var(87)='tmp'
                                                                              
         i=87
         do j=1,12
            write(cat,'(i2.2)')j
            if(cat(1:1).eq.' ')cat(1:1)='0'
            if(cat(2:2).eq.' ')cat(2:2)='0'
            var(i+j)= 'a'//cat   ! monthly albedo
         enddo
                                                                           
         var(100) = 'slp'      !terrain slope index, dominant category
         var(101) = 'gnx'      !max greenness fraction
         var(102) = 'gnn'      !min greenness fraction
         var(ngrids)   = 'zin'

         comment(1) = 'h-point (mass) latitudes \0'
         comment(2) = 'h-point (mass) longitudes \0'
         comment(3) = 'v-point (mass) latitudes \0'
         comment(4) = 'v-point (mass) longitudes \0'
         comment(5)= 'land fraction  \0'
         comment(6)= 'land use dominant category \0'
         comment(7)= 'land-water mask (0=water; 1 otherwise) \0'
         comment(8)= 'soiltype top layer dominant category \0'
         comment(9)= 'soiltype bot layer dominant category \0'
                                                                                         
         i=10
         comment(i)= 'sin of projection rotation (rad) \0'
         comment(i+1)= 'cosine of projection rotation (rad) \0'
         comment(i+2)= 'horizontal component coriolis parameter \0'
         comment(i+3)= 'vertical component coriolis parameter \0'
         comment(i+4)= 'maximum snow albedo (%) over land only \0'
         comment(i+5)= 'standard deviation of elevation data (m)\0'
         comment(i+6)= 'mean longitudinal terrain slope (m/m)\0'
         comment(i+7)= 'mean latitudinal terrain slope (m/m)\0'
         comment(i+8)= 'average terrain elevation (c-stagger) (m)\0'
                                                                               
         i=18
         do j=1,24
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          var(i+j)= 'u'//cat
          comment(i+j)= '% dist land use category '//cat//' \0'
         enddo

         i=42
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)='% dist top layer soiltype category '//cat
         enddo

         i=58
         do j=1,16
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= '% dist bot layer soiltype category '//cat
         enddo

         i=74
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= 'vegetation greenness fraction: mon = '//cat
         enddo

         comment(87)='1 degree mean annual soiltemp (deg k)'
         i=87
         do j=1,12
          write(cat,'(i2.2)')j
          if(cat(1:1).eq.' ')cat(1:1)='0'
          if(cat(2:2).eq.' ')cat(2:2)='0'
          comment(i+j)= 'climatological albedo: mon = '//cat
         enddo

         comment(100)=  'terrain slope index'
         comment(101)=  'annual max greenness fraction'
         comment(102)=  'annual min greenness fraction'
         comment(ngrids)= '\0'

      endif
      return
      end

c ********************************************************************

	subroutine read_dem(unit_no,unit_name,nn1,nn2,i1,i2
     &,data,istat)

	implicit none

	integer countx,county,unit_no,nn1,nn2

        character  cdata(nn1,nn2)*2
        integer,   allocatable :: idata(:,:)
	real       data(nn1,nn2)
	integer len, lend, i1, i2, ia, istat
        real multiplier
c       logical l1,l2
	character*(*) unit_name
        character*1   ctiletype

c	open(unit_no,file=unit_name,status='old',access='direct',
c	. recl=nn2*nn1*2)
c	inquire(unit_no,exist=l1,opened=l2)
c	read(unit_no,rec=1) idata

	call s_len(unit_name,len)
        call get_directory_length(unit_name,lend)
        ctiletype=unit_name(lend+1:lend+1)

        if(.not.allocated(idata))then
           allocate (idata(nn1,nn2),stat=istat)
           if(istat.ne.0)then
              print*,'unable to allocate idata array: read_dem'
              print*,'nn1/nn2/istat: ',nn1,nn2,istat
              return
           endif
        endif


        multiplier=1.0
	if(ctiletype.eq.'t'.or.ctiletype.eq.'u')then
           if(nn1.eq.1201)then
              open(unit_no,file=unit_name,status='old',
     .form='unformatted')
              read(unit_no)data
              close(unit_no)
           else
              call read_binary_field(cdata,i1,i2,nn1*nn2,unit_name,len)

              do county=1,nn2
              do countx=1,nn1
               idata(countx,county) = ia (cdata(countx,county),2,0)
              enddo
              enddo
              if(ctiletype.eq.'t')multiplier=.01  !(for t data these are temps * 100)
           endif
        else
           call read_binary_field(idata,i1,i2,nn1*nn2,unit_name,len)
        endif

        if(nn1.le.1200)then
	do county=1,nn2
	do countx=1,nn1
	 if(idata(countx,county).eq.-9999) idata(countx,county)=0
	  data(countx,county)=float(idata(countx,nn2-county+1))
     &*multiplier
c sg97 initial data (dem format) starts in the lower-left corner;
c sg97 this format is wrapped around to have upper-left corner as its start.
c
c js00 some machines do not account for signed integers
	   if(data(countx,county).ge.15535.0)
     &data(countx,county)=data(countx,county)-65535

	enddo
	enddo
        endif
 
        if(allocated (idata))deallocate(idata)

ccc	 close(unit_no)
	return
	end
c ********************************************************************

        subroutine read_dem_g(unit_no,unit_name,nn1,nn2,nn3,nn4
     &,nofr,i1,i2,data,istat)

        implicit none
        integer  countx,county,countz
        integer  unit_no,nn1,nn2,nn3,nn4,nofr
        integer  len, lend, i1, i2, i
        integer  istat

        real     data(nn1,nn2,nn3,nn4)
        integer, allocatable ::  idata(:,:,:)

c       logical  l1,l2

        character*(*) unit_name
        character*1   ctype

c       open(unit_no,file=unit_name,status='old',access='direct',
c       . recl=nn2*nn1*2)
c       inquire(unit_no,exist=l1,opened=l2)
c       read(unit_no,rec=1) idata

        if(.not.allocated(idata))then
           print*,'allocate idata in read_dem_g'
           allocate (idata(nn4,nn1,nn2),stat=istat)
           if(istat.ne.0)then
              print*,'unable to allocate idata array: read_dem_g'
              print*,'nn1/nn2/nn4/istat: ',nn1,nn2,nn4,istat
              return
           endif
        endif

        call s_len(unit_name,len)
        call get_directory_length(unit_name,lend)
        ctype=unit_name(lend+1:lend+1)

        print*,'read_dem_g: tile type = ',ctype

        call read_binary_field(idata,i1,i2,nn1*nn2*nn4,unit_name,len)

c       if(nn1.ne.1250 .and. nn2.ne.1250)then
        if(ctype.ne.'a' .and. 
     &     ctype.ne.'g' .and.
     &     ctype.ne.'i' .and.
     &     ctype.ne.'m')then

           do county=1,nn2
           do countx=1,nn1
           do countz=1,nn4

              if(idata(countz,countx,county).eq.-9999)
     &idata(countz,countx,county)=0

              data(countx,county,nofr,countz)=
     &float(idata(countz,countx,nn2-county+1))

c sg97 initial data (dem format) starts in the lower-left corner;
c sg97 this format is wrapped around to have upper-left corner as its start.


           enddo
           enddo
           enddo

        else   !new greenfrac data starts at 90s

           do county=1,nn2
           do countx=1,nn1
           do countz=1,nn4

              data(countx,county,nofr,countz)=
     &float(idata(countz,countx,county))

           enddo
           enddo
           enddo

        endif

c we'll resurrect the actual categories later
c but for now we want categories 1-9 for
c terrain slope index.
        if(ctype == 'i')then
           do county=1,nn2
           do countx=1,nn1
              if(data(countx,county,1,1).eq.13)then
                 data(countx,county,1,1)=8
              elseif(data(countx,county,1,1).eq.0)then
                 data(countx,county,1,1)=9
              endif
           enddo
           enddo
c          where(data .eq. 13)data = 8
c          where(data .eq. 0)data = 9
        endif

        if(allocated(idata))deallocate(idata)

ccc      close(unit_no)
        return
        end

c +------------------------------------------------------------------+
	subroutine jcl
	character*(*) filenm,formt

c	-------------------------------------------------------
	entry jclget(iunit,filenm,formt,iprnt,istatus)
c
c	  this routine access an existing file with the file name of
c	    filenm and assigns it unit number iunit.
c
		if(iprnt.eq.1) then
	print*,' opening input unit ',iunit,' file name ',filenm
	print*,'		format  ',formt
	endif

	open(iunit,status='old',file=filenm,form=formt,err=1)

	istatus=1
	return

 1    istatus = 0
	return

	end
!
!----------------------------------------------------------------------
        function ia(chr,n,ispval)                         
!                                                              
!  purpose: to convert a n-bytes character (chr) to integer ia. 
!        ** the integer data file is saved as a n-byte character
!           data file. this function is used to recover the    
!           character data to the integer data.               
!                                                            
!  n      --- the number of bytes in chr                    
!  ispval --- default value for the negative integer.      
!                                                       
        character*2 :: chr                                
        integer  n, ii1, ii2, jj, isn, m, nbit, mshft, ia2, ispval
        integer  bit_1, bit_2                            
!                                                    
        bit_1 = o'200'     ! binary '10000000'        
        bit_2 = o'377'     ! binary '11111111'       
        ia    = 0                                   
!                                                
        ii1 = ichar(chr(1:1))                     
        if(ii1 < 0) ii1=ii1+256

! .. get the sign -- isn=0 positive, isn=1 negative:
        jj  = iand(ii1,bit_1)                        
        isn = ishft(jj,-7)                          
!                                                
! .. for negative number:
!    because the negative integers are represented by the supplementary
!    binary code inside machine.
!                              
        if (isn.eq.1) then    
          do m = n+1,4   
           nbit = (m-1)*8   
           jj = ishft(bit_2,nbit)
           ia = ieor(jj,ia)     
          end do                
        endif                   
!                              
!   .. get the byte from chr: 
        do m = 1,n          
         ii2 = ichar(chr(m:m)) 
         if(ii2 < 0) ii2=ii2+256
         mshft = (n-m)*8      
         ia2   = ishft(ii2,mshft)
!   .. the abs(integer):          
         ia = ieor(ia,ia2)     
        end do                 
!                              
        if (ia.lt.0) ia = ispval
!                            
        return                
        end
c
c--------------------------------------------------------------------
c
       subroutine polar_gp(lat,lon,x,y,dx,dy,nx,ny,dir)
c
      include 'trigd.inc'
       real lat,lon,x,y,dx,dy,
     1        erad,tlat,tlon                                      ! ,plat,plon,
     1        xdif,ydif
c
       integer   nx,ny
       integer   idir  !positive (1.) going from center lat/lon to sw x/y;
c                       negative (-1.) going from sw lat/lon to center x/y.
c
       rad=3.141592654/180.

       call get_earth_radius(erad,istatus)
       if(istatus .ne. 1)then
           write(6,*)' error calling get_earth_radius'
           stop
       endif

!      calculate xy coordinates at domain center
       call latlon_to_xy(lat,lon,erad,xdif,ydif)

       x=xdif+(1.-dir*float(nx)/2.)*dx
       y=ydif+(1.-dir*float(ny)/2.)*dy
 
       return
 
       end
c
c------------------------------------------------------------
c
      subroutine blend_topo(nnxp,nnyp,lats,lons
     1,topt_10,topt_10_s,topt_10_ln,topt_10_lt
     1,topt_30,topt_30_s,topt_30_ln,topt_30_lt
     1,topt_out,topt_out_s,topt_out_ln,topt_out_lt)

      implicit none
      integer  nnxp,nnyp
      integer  i,j
      integer  icount_10
      integer  icount_30
      integer  icount_ramp

      real     alat1n
      real     alat2n
      real     alat1s
      real     alat2s
      real     nboundary
      real     sboundary
      real     frac10
      real     width

      real     lats(nnxp,nnyp)
      real     lons(nnxp,nnyp)
      real     topt_out(nnxp,nnyp)
      real     topt_out_s(nnxp,nnyp)
      real     topt_out_ln(nnxp,nnyp)
      real     topt_out_lt(nnxp,nnyp)
      real     topt_10(nnxp,nnyp)
      real     topt_10_s(nnxp,nnyp)
      real     topt_10_ln(nnxp,nnyp)
      real     topt_10_lt(nnxp,nnyp)
      real     topt_30(nnxp,nnyp)
      real     topt_30_s(nnxp,nnyp)
      real     topt_30_ln(nnxp,nnyp)
      real     topt_30_lt(nnxp,nnyp)


      do i = 1,nnxp
      do j = 1,nnyp

! select 30s or 10m topo data for this grid point (or a blend)

!              check whether 30s data is missing or zero
         if(topt_30(i,j) .eq. 1e30 .or. topt_30(i,j) .eq. 0.
!    1                              .or.
!                  are we in the pittsburgh data hole?
!    1            (lats(i,j) .gt. 39.7 .and. lats(i,j) .lt. 41.3 .and.
!    1             lons(i,j) .gt.-79.3 .and. lons(i,j) .lt.-77.7)
!
     1                                                      )then 

!                  use 10 min data
            topt_out(i,j) = topt_10(i,j)
            topt_out_s(i,j)=topt_10_s(i,j)
            topt_out_ln(i,j)=topt_10_ln(i,j)
            topt_out_lt(i,j)=topt_10_lt(i,j)
            icount_10 = icount_10 + 1

         else ! use 30s data, except ramp to 10m if near data boundary

! determine the northern boundary of the 30s data at this lon
            if(lons(i,j).ge.-129..and.
     +         lons(i,j).le.-121.)then       
               nboundary = 51.
            elseif(lons(i,j).ge.-121..and.
     +             lons(i,j).le.-120.)then     
                   nboundary = 51. - lons(i,j) - (-121.)
            elseif(lons(i,j).ge.-120..and.
     +             lons(i,j).le.-118.)then     
                   nboundary = 50.
            elseif(lons(i,j).ge.-118..and.
     +             lons(i,j).le.-117.)then     
                   nboundary = 50. + lons(i,j) - (-118.)
            elseif(lons(i,j).ge.-117..and.
     +             lons(i,j).le. -89.)then     
                   nboundary = 51.
            elseif(lons(i,j).ge. -89..and.
     +             lons(i,j).le. -85.)then     
                   nboundary = 50.
            elseif(lons(i,j).ge. -85..and.
     +             lons(i,j).le. -83.)then     
                   nboundary = 49.
            elseif(lons(i,j).ge. -83..and.
     +             lons(i,j).le. -81.)then     
                   nboundary = 48.
            elseif(lons(i,j).ge. -81..and.
     +             lons(i,j).le. -73.)then     
                   nboundary = 46.
            elseif(lons(i,j).ge. -73..and.
     +             lons(i,j).le. -67.)then     
                   nboundary = 47.
            elseif(lons(i,j).ge. -67.)then     
                   nboundary = 46.
            else
                   nboundary = 51.
            endif

            alat1n = nboundary - 0.3
            alat2n = nboundary - 0.1

! determine the southern boundary of the 30s data at this lon
            if    (lons(i,j) .le. -127.)then         
                   sboundary = 49. 
            elseif(lons(i,j) .le. -126.)then         
                   sboundary = 48. 
            elseif(lons(i,j) .le. -125.)then         
                   sboundary = 40. 
            elseif(lons(i,j) .le. -124.)then         
                   sboundary = 37. 
            elseif(lons(i,j) .le. -123.)then         
                   sboundary = 36. 
            elseif(lons(i,j) .le. -122.)then         
                   sboundary = 35. 
            elseif(lons(i,j) .le. -120.)then         
                   sboundary = 33. 
            elseif(lons(i,j) .le. -118.)then     
                   sboundary = 32. 
            elseif(lons(i,j) .le. -107.)then     
                   sboundary = 30. 
            elseif(lons(i,j) .le. -103.)then     
                   sboundary = 28. 
            elseif(lons(i,j).ge.-103. .and.
     +             lons(i,j).le.-102.)then       
                   sboundary = 25. +  (-102. - lons(i,j)) * 3.
            elseif(lons(i,j).ge.-102. .and.
     +             lons(i,j).le. -99.)then       
                   sboundary = 25.
            elseif(lons(i,j).ge.-99.  .and.
     +             lons(i,j).le. -98.)then       
                   sboundary = 24. +  ( -98. - lons(i,j))
            elseif(lons(i,j).ge.-98. )then       
                   sboundary = 24.
            endif

            alat1s = sboundary + 0.3
            alat2s = sboundary + 0.1

! decide whether to use 30s or 10m data (or a blend)

            if  (  lats(i,j) .ge. alat2n)then    ! use 10m data
                   topt_out(i,j) = topt_10(i,j)
                   topt_out_s(i,j)=topt_10_s(i,j)
                   topt_out_ln(i,j)=topt_10_ln(i,j)
                   topt_out_lt(i,j)=topt_10_lt(i,j)
                   icount_10 = icount_10 + 1

            elseif(lats(i,j) .ge. alat1n .and. 
     1             lats(i,j) .le. alat2n)then

! between alat1n and alat2n,        use weighted average

                   width = alat2n - alat1n
                   frac10 = (lats(i,j) - alat1n) / width
                   topt_out(i,j) = topt_10(i,j) * frac10 
     1                           + topt_30(i,j) * (1. - frac10)
                   topt_out_s(i,j) = topt_10_s(i,j) * frac10
     1                             + topt_30_s(i,j) * (1. - frac10)
                   topt_out_ln(i,j) = topt_10_ln(i,j) * frac10
     1                              + topt_30_ln(i,j) * (1. - frac10)
                   topt_out_lt(i,j) = topt_10_lt(i,j) * frac10
     1                              + topt_30_lt(i,j) * (1. - frac10)
                   icount_ramp = icount_ramp + 1

                   if(icount_ramp .eq. (icount_ramp/5) * 5 )then       
                      write(6,*)
                      write(6,*)'in blending zone, nboundary = '
     1                                       ,nboundary,alat1n,alat2n       
                      write(6,*)'lat/lon/frac =',lats(i,j)
     1                               ,lons(i,j) ,frac10
                      write(6,*)'topt_30      =',topt_30(i,j)
                      write(6,*)'topt_10      =',topt_10(i,j)
                      write(6,*)'topt_out     =',topt_out(i,j)
                   endif

            elseif(lats(i,j) .ge. alat1s .and. 
     1             lats(i,j) .le. alat1n)then
                   topt_out(i,j) = topt_30(i,j)
                   topt_out_s(i,j)=topt_30_s(i,j)
                   topt_out_ln(i,j)=topt_30_ln(i,j)
                   topt_out_lt(i,j)=topt_30_lt(i,j)
                   icount_30 = icount_30 + 1       ! use 30s data

            elseif(lats(i,j) .ge. alat2s .and. 
     1             lats(i,j) .le. alat1s)then

! between alat1s and alat2s,        use weighted average

                   width = alat1s - alat2s
                   frac10 = (alat1s - lats(i,j)) / width
                   topt_out(i,j) = topt_10(i,j) * frac10
     1                           + topt_30(i,j) * (1. - frac10)
                   topt_out_s(i,j) = topt_10_s(i,j) * frac10 
     1                             + topt_30_s(i,j) * (1. - frac10)
                   topt_out_ln(i,j) = topt_10_ln(i,j) * frac10
     1                              + topt_30_ln(i,j) * (1. - frac10)
                   topt_out_lt(i,j) = topt_10_lt(i,j) * frac10
     1                              + topt_30_lt(i,j) * (1. - frac10)
                   icount_ramp = icount_ramp + 1

                   if(icount_ramp .eq. (icount_ramp/5) * 5 )then       
                      write(6,*)
                      write(6,*)'in blending zone, sboundary = '
     1                                   ,sboundary,alat1s,alat2s       
                      write(6,*)'lat/lon/frac =',lats(i,j)
     1                           ,lons(i,j), frac10
                      write(6,*)'topt_30      =',topt_30(i,j)
                      write(6,*)'topt_10      =',topt_10(i,j)
                      write(6,*)'topt_out     =',topt_out(i,j)
                   endif

            elseif(lats(i,j) .le. alat2s)then    
                   topt_out(i,j) = topt_10(i,j)    ! use 10m data
                   topt_out_s(i,j)=topt_10_s(i,j)
                   topt_out_ln(i,j)=topt_10_ln(i,j)
                   topt_out_lt(i,j)=topt_10_lt(i,j)
                   icount_10 = icount_10 + 1

            else
                   write(6,*)' software error in gridgen_model.f'
                   write(6,*)' lat/lon = ',lats(i,j),lons(i,j)
                   stop

            endif ! test to see if we blend the data

         endif ! 30s data check

      enddo ! j
      enddo ! i 

      return
      end
c
c ---------------------------------------------------------------
c
      subroutine get_meanlattemp(path_to_tiles,temp,istat)
c
      implicit      none

      integer       istat
      integer       n,i
      integer       ldir
      character*(*) path_to_tiles

      real          temp(180)

      istat=0
      call s_len(path_to_tiles,ldir)
      open(22,file=path_to_tiles(1:ldir)//'/latmeantemp.dat'
     &,form='formatted',status='old',iostat=istat)
      if(istat.ne.0) then
	write(6,*) 'insert bogus temp of 280'
	temp=280.0
	istat=1
	return
!	goto 3
      endif
      do i=1,180
         read(22,222,err=4)temp(i)
      enddo
 
      close(22)

      istat=1

      return

c     print*,'rmeantemp 1/2/3/4/5; ',rmeantemp(1),rmeantemp(45),rmeantemp(90)&
c          ,rmeantemp(135),rmeantemp(180)

222   format(1x,f6.2)

  3   print*,'error: opening latmeantemp file '
      print*,'path_to_tiles: ',path_to_tiles(1:ldir+3),ldir
      return
  4   print*,'error: reading latmeantemp file '

      return
      end
c
c ---------------------------------------------------------
c
      subroutine eval_localization(cstaticdir,nest,localize
     .,cgrid_fname,la1_dom,lo1_dom,istatus)

c routine performs the following:
c 1. tests if the static file has been generated for this domain
c 2. tests if the static variables consistent with the current namelist specs
c
c if 1 or 2 is false then we either localize or re-localize
c the domain; "localize" is returned indicating such (true or false).
c
      implicit  none

      integer   nest,ifl
      integer   nx,ny
      integer   nx_dom,ny_dom
      integer   lf,ldir,len_cfl
      integer   istatus

      character cstaticdir*200
      character cstaticfile*200
      character cgrid_fname*10
      character c6_maproj*6
      character c8_maproj*8
      character cfl*3
      character cnest*2

      logical   localize
      logical   static_exists

      real      dx,dy,la1,lo1,lov,latin1,latin2
      real      grid_spacing_dom_m,grid_spacing_m
      real      la1_dom,lo1_dom

      call getenv('force_localization',cfl)
      call s_len(cfl,len_cfl)

      if(len_cfl.eq.1)then
         print*,'environment variable force_localization= ',cfl
         read(cfl,'(i1.1)')ifl
         if(ifl.eq.nest)then
            localize = .true.
            return
         endif
      elseif(len_cfl.gt.1)then
         call downcase(cfl,cfl)
         if(cfl.eq.'all')then
            localize = .true.
            return
         else
            print*,'unknown force_localization setting ',cfl
            stop
         endif
      endif

      if(cgrid_fname.eq."wrfsi")then
         localize=.false.
      else
         localize=.true.
         return
      endif

      call get_c6_maproj(c6_maproj,istatus)

      call s_len(cstaticdir,ldir)

      call rd_static_attr(cstaticdir,nest,cgrid_fname
     .,nx, ny, dx, dy, la1, lo1, latin1, latin2, lov
     .,c8_maproj,istatus)

      if(c8_maproj.eq.'lambert'.and. c6_maproj.ne.'lambrt')then
         print*,'static file map-proj differs from namelist'
         print*,'**** relocalize this domain ****'
         localize=.true.
         return
      elseif(c8_maproj.eq.'polar'.and. c6_maproj.ne.'plrstr')then
         print*,'static file map-proj differs from namelist'
         print*,'**** relocalize this domain ****'
         localize=.true.
         return
      elseif(c8_maproj.eq.'mercator'.and. c6_maproj.ne.'merctr')then
         print*,'static file map-proj differs from namelist'
         print*,'**** relocalize this domain ****'
         localize=.true.
         return
      endif

      write(cnest,'(i2.2)')nest
      cstaticfile=trim(cstaticdir)//'static.'//trim(cgrid_fname)
      cstaticfile=trim(cstaticfile)//'.d'//cnest

      if (istatus .ne. 1) then
        print *,' eval localization: did not read wrf static file'
        print *,' status = ',istatus
        localize = .true.
        return
      end if 

      if(la1.ne.la1_dom)localize=.true.
      if(lo1.ne.lo1_dom)localize=.true.

      call get_grid_dim_xy(nx_dom,ny_dom,istatus) 
      if(nx.ne.nx_dom)localize=.true.
      if(ny.ne.ny_dom)localize=.true.
c
      call get_grid_spacing(grid_spacing_dom_m,istatus)
      if(dx.ne.grid_spacing_dom_m)localize=.true.

      return
      end
c
c --------------------------------------------------------
c
      subroutine rd_static_attr(cstaticdir,nest
     .,cgrid_fname,nx,ny,dx,dy,la1,lo1,latin1
     .,latin2,lov,c8_maproj,istatus)

c routine performs the following:
c 1. tests if the static file has been generated for this domain or nest?
c 2. tests if the static variables consistent with the current namelist specs
c
c if 1 or 2 is false then we need to either localize or re-localize
c the domain and "localize" is returned indicating such (true or false).
c
      implicit  none

      integer   nest
      integer   nx,ny
      integer   istatus

      character staticfile*200
      character cstaticdir*200
      character cgrid_fname*10
      character c8_maproj*8
      character cnest*2

      logical   localize
      logical   static_exists

      real      dx,dy,la1,lo1,lov,latin1,latin2
      real      grid_spacing_dom_m,grid_spacing_m

      istatus = 1

      staticfile=trim(cstaticdir)//'static.'//cgrid_fname
      if(trim(cgrid_fname).eq.'wrfsi')then
         write(cnest,'(i2.2)')nest
         staticfile=trim(staticfile)//'.d'//cnest
      endif

      inquire(file=staticfile, exist=static_exists)

      if (static_exists) then

        print*,'static file exists: rd_static_attr'
        print*,'static filename: ',trim(staticfile)

        call rd_static_attr_sub(staticfile, nx, ny
     .,la1, latin1, latin2, lo1, lov, dx, dy
     .,c8_maproj,istatus)

        if (istatus .ne. 1) then
	write(6,*) '2nd time'
           print '(a,i5)', ' error reading wrf static file: ',istatus
           return
        end if

        if (lov .lt. -180.) lov = lov + 360.
        if (lov .gt. 180.) lov = lov - 360.
        if (lo1 .lt. -180.) lo1 = lo1 + 360.
        if (lo1 .gt. 180.) lo1 = lo1 - 360.

      else

        print '(a)', 'static file not found: ', trim(staticfile)
        istatus = 0

      endif

      return
      end
c
c --------------------------------------------------------------------------
c
      subroutine gridcompare(nx,ny,ii,datain,datalm,istatus)

      implicit none

      integer nx,ny
      integer i,ii,j
      integer cntww
      integer cntwnw
      integer cntnww
      integer cntnwnw
      integer istatus

      real datain(nx,ny)  !input data to compare to land mask
      real datalm(nx,ny)  !land mask
      real rmsng,thresh

      call get_r_missing_data(rmsng,istatus)

      cntwnw =0
      cntww  =0
      cntnww =0
      cntnwnw=0

      if(ii == 1)then
         thresh=0.0                 !terrain
         print*,'array comparison:   terrain'
      elseif(ii == 2 .or. ii==3)then
         thresh=14                  !soil texture
         print*,'array comparison:   soil texture'
      elseif(ii == 4)then 
         thresh=0.0                 !max greenness
         print*,'array comparison:   greeness: mo 6'  ! max greenness'
      elseif(ii == 5)then
         thresh=0.0                 !min greenness
         print*,'array comparison:   min greenness'
      elseif(ii == 6)then
         thresh=rmsng               !deep soil temp
         print*,'array comparison:   deep soil temp'
      elseif(ii == 7)then
         thresh=0.0                 !terrain slope index
         print*,'array comparison:   terrain slope index'
      elseif(ii == 8)then
         thresh=0.08                !albedo: month 6
         print*,'array comparison:   albedo: month 6'
      elseif(ii == 9)then
         thresh=0.08                !max snow albedo
         print*,'array comparison:   max snow albedo'
      elseif(ii == 10)then
         thresh=16                  !dominant cat landuse
         print*,'array comparison:   dominant cat landuse'
      endif


      do j=1,ny
      do i=1,nx

         if(datalm(i,j)    .eq.0 .and. datain(i,j).ne.thresh)then
            cntwnw=cntwnw+1
         elseif(datalm(i,j).eq.0 .and. datain(i,j).eq.thresh)then
            cntww=cntww+1
         elseif(datalm(i,j).eq.1 .and. datain(i,j).eq.thresh)then
            cntnww=cntnww+1
         elseif(datalm(i,j).eq.1 .and. datain(i,j).ne.thresh)then
            cntnwnw=cntnwnw+1
         endif

      enddo
      enddo

      print*,'============================================'
      print*,'============================================'
      print*,'land mask=water, data array = water        : ',cntww
      print*,'land mask=not water, data array = not water: ',cntnwnw
      print*,'land mask=water, data array = not water    : ',cntwnw
      print*,'land mask=not water, data array = water    : ',cntnww
      print*

      return
      end
