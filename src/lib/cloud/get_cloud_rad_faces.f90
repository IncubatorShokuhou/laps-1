     
     subroutine get_cloud_rad_faces(             &
                obj_alt,obj_azi,                 & ! i
                solalt,solazi,                   & ! i 
                clwc_3d,cice_3d,rain_3d,snow_3d, & ! i
                topo_a,                          & ! i
                ni,nj,nk,idb,jdb,                & ! i
                heights_3d,                      & ! i 
                transm_2t,transm_3d,transm_4d)     ! o

     include 'trigd.inc'

!    calculate 3d radiation field looping through each of the 6 faces in
!    the domain.

     use mem_namelist, only: r_missing_data,earth_radius,grid_spacing_m &
                            ,aod,aero_scaleht,angstrom_exp_a,redp_lvl
     use mem_allsky, only: uprad_4d ! (upward spectral irradiance)
     use mem_allsky, only: ext_g, nc
     use mem_allsky, only: aod_3d   ! (extinction coefficient)            ! i
     use mem_allsky, only: mode_aero_cld
     include 'rad.inc' ! e.g. for ext_o, o3_du

     trans(od) = exp(-min(od,80.))
     scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! x/scurve range is 0-1

!    for each face, set start/end points to loop in the i,j,k dimension
!    index of 1 points to the 1st element, 2 points to the last element
!                  t  b  w  e  n  s
     real i1(6)  / 1, 1, 1, 2, 1, 1/            
     real i2(6)  / 2, 2, 1, 2, 2, 2/           
     real j1(6)  / 1, 1, 1, 1, 2, 1/            
     real j2(6)  / 2, 2, 2, 2, 2, 1/           
     real k1(6)  / 2, 1, 1, 1, 1, 1/
     real k2(6)  / 2, 1, 2, 2, 2, 2/
     character*6 cfaces(6) /'top','bottom','west','east','north','south'/

     real heights_3d(ni,nj,nk)
     real clwc_3d(ni,nj,nk)  ! kg/m**3
     real cice_3d(ni,nj,nk)  ! kg/m**3
     real rain_3d(ni,nj,nk)  ! kg/m**3
     real snow_3d(ni,nj,nk)  ! kg/m**3
     real b_alpha_3d(ni,nj,nk) ! m**-1          
     real transm_3d(ni,nj,nk) ! direct transmission plus forward scattered
     real transm_3t(ni,nj,nk) ! terrain 3d transmission
     real transm_2t(ni,nj)    ! terrain 2d transmission
     real transm_4d(ni,nj,nk,nc) ! color information added
     integer kterr(ni,nj)

     real obj_alt(ni,nj),obj_azi(ni,nj)
     real topo_a(ni,nj),terr_max_path_a(ni,nj)
     real projrot(ni,nj)

     real heights_1d(nk)

     real sprad_to_nl(nc)
     real trans_c(nc)
     real bi_coeff(2,2),tri_coeff(2,2,2),b_alpha_cube(2,2,2)
     equivalence (s,scurr) ! needed only during transition from s to scurr

     integer htlutlo,htluthi
     parameter (htlutlo = -1000)
     parameter (htluthi = 30000)
     real htlut(htlutlo:htluthi)

     integer*8 nsteps,nfacesteps

!    these parameters can be obtained/updated from the 'cloud_rad' module
!    backscattering efficiencies
     real, parameter :: bksct_eff_clwc    = .063
     real, parameter :: bksct_eff_cice    = .14
     real, parameter :: bksct_eff_rain    = .063
     real, parameter :: bksct_eff_snow    = .14
     real, parameter :: bksct_eff_graupel = .30
     real, parameter :: bksct_eff_aero    = .125 

!    scattering efficiencies
     real, parameter :: q_clwc    = 2.0
     real, parameter :: q_cice    = 2.0
     real, parameter :: q_rain    = 1.0
     real, parameter :: q_snow    = 1.0
     real, parameter :: q_graupel = 1.0

!    densities
     real, parameter :: rholiq     =   1e3 ! kilograms per cubic meter
     real, parameter :: rhosnow    = .07e3 ! kilograms per cubic meter
     real, parameter :: rhograupel = .50e3 ! kilograms per cubic meter

!    effective radii
     real, parameter :: reff_clwc    = .000007 ! m
     real, parameter :: reff_cice    = .000034 ! m
     real, parameter :: reff_rain    = .000750 ! m
     real, parameter :: reff_snow    = .004000 ! m
     real, parameter :: reff_graupel = .010000 ! m

     logical l_same_point

     clwc2alpha = 1.5 / (rholiq  * reff_clwc)
     cice2alpha = 1.5 / (rholiq  * reff_cice)
     rain2alpha = 1.5 / (rholiq  * reff_rain)
     snow2alpha = 1.5 / (rhosnow * reff_snow)
     pice2alpha = 1.5 / (rhograupel * reff_graupel)

     twi_alt = -4.5

!    initialize transm_3d
     write(6,*)' initialize transm_3d and kterr'
     do k = 1,nk
         transm_3d(:,:,k) = r_missing_data
         where(topo_a(:,:) .ge. heights_3d(:,:,k)) ! below terrain
             transm_3d(:,:,k) = 0.
             kterr(:,:) = k
         endwhere
     enddo ! k

     write(6,*)' call get_terrain_shadows...'
     call get_terrain_shadows(heights_3d,topo_a,obj_alt,obj_azi,kterr,ni,nj,nk,idb,jdb,transm_3t,transm_2t)

     transm_4d = 0.                   

     do ic = 1,nc
       call nl_to_sprad(1.,1,wa(ic),sprad)
       sprad_to_nl(ic) = 1. / sprad
     enddo ! ic

     if(mode_aero_cld .lt. 3)then
       b_alpha_3d = clwc_3d * clwc2alpha * bksct_eff_clwc &
                  + cice_3d * cice2alpha * bksct_eff_cice &
                  + rain_3d * rain2alpha * bksct_eff_rain &
                  + snow_3d * snow2alpha * bksct_eff_snow 
     else
       b_alpha_3d = clwc_3d * clwc2alpha * bksct_eff_clwc &
                  + cice_3d * cice2alpha * bksct_eff_cice &
                  + rain_3d * rain2alpha * bksct_eff_rain &
                  + snow_3d * snow2alpha * bksct_eff_snow &
                  + aod_3d               * bksct_eff_aero
     endif

     write(6,*)' subroutine get_cloud_rad_faces: solar alt/az ',solalt,solazi

     ntot = 0
     nsteps = 0
     refraction = 0.5 ! initialize

     heights_1d(:) = heights_3d(ni/2,nj/2,:)

     do k = 1,nk-1    
       htlow  = heights_1d(k)
       hthigh = heights_1d(k+1)
       do klut = nint(htlow),nint(hthigh)
         rklut = klut
         if(klut .ge. htlutlo .and. klut .le. htluthi)then
           htlut(klut) = float(k) + (rklut-htlow)/(hthigh-htlow)
         endif
       enddo ! klut
     enddo ! k

     dhmin = heights_1d(2)  - heights_1d(1)
     dhmax = heights_1d(nk) - heights_1d(1)
     write(6,*)' dhmin/dhmax = ',dhmin,dhmax

     faceperim = 0.4
     facestepij = 0.80 ! 0.85 causes missing points at 15.5 deg (top face)
     raysteps = 0.5 * grid_spacing_m

     refr_mn = 0.0

!    different criteria might be used at high altitude depending on how far
!    away the limb is, in turn related to horz_dep for the observer
     if(solalt .ge. twi_alt)then ! daylight or early twilight

      do if = 1,6
       idebug = 1 ! initialize for first ray trace
       write(6,*)
       is = 1 + (ni-1)*(i1(if)-1)
       ie = 1 + (ni-1)*(i2(if)-1)
       js = 1 + (nj-1)*(j1(if)-1)
       je = 1 + (nj-1)*(j2(if)-1)
       ks = 1 + (nk-1)*(k1(if)-1)
       ke = 1 + (nk-1)*(k2(if)-1)

       write(6,*)cfaces(if),' face'
       write(6,*)'rangei = ',is,ie  
       write(6,*)'rangej = ',js,je  

       nnew = 0
       nfacesteps = 0

       ris = is-faceperim; rie = ie+faceperim
       rjs = js-faceperim; rje = je+faceperim
       rks = ks          ; rke = ke

!      exception needed if obs_alt = 90.
       objalt = obj_alt(nint(ris),nint(rjs))
       facesteph = min(max(grid_spacing_m*tand(objalt),dhmin),dhmax) ! * 0.3

       write(6,*)'rangek = ',ks,ke,facesteph

       do rit = ris,rie,facestepij
       do rjt = rjs,rje,facestepij

!      do rkt = rks,rke,facestepk
       do htt = heights_1d(ks),heights_1d(ke),facesteph
         rkt = htlut(nint(htt))

         it = rit; jt = rjt; kt = rkt
         id = nint(rit); jd = nint(rjt); kd = nint(rkt)

!        kl = max(min(int(rkt),nk-1),1) ; kh = kl + 1
!        fk = rkt - kl
!        htt = heights_1d(kl) * (1.-fk) + heights_1d(kh) * fk

         ihit_terrain = 0
 
         slast = 0
         b_alpha_new = 0.
         btau = 0.

!        we presently do ray marching a constant distance intervals.
!        it may be more efficient to have successive steps march 
!        to the next 3d grid box boundary. we can then interpolate
!        from these "end points" to form the needed integrated values
!        at these grid box boundaries. we'd also want to know the
!        integrated values at the location where the ray comes closest
!        to the center of the grid box that is traverses.

!        when doing bilinear or trilinear interpolation, we're looking at
!        squares/cubes that have grid points lying on the vertices and the
!        'int' operation is used. when assigning the radiance values to the
!        grid, we are considering the nearest grid point using the 'nint'
!        operation. here the cube is centered on a grid point.

!        intersections with terrain are considered by finding maxima in the
!        bilinearly interpolated terrain field, relative to the ray height.
!        these maxima are thought to be located along lines connecting two
!        adjacent or diagonally adjacent grid points.

!        start ray trace at this point
         if(idebug .eq. 1)write(6,1)if,it,jt,kt
1        format('if,it,jt,kt',4i3)

         do ls = 0,10000 ! max number of ray segments

           if(ls .eq. 0)then ! values at start of trace                
             objalt = obj_alt(id,jd) + refr_mn
             objazi = obj_azi(id,jd)
             dids = -sind(objazi)*cosd(objalt)/grid_spacing_m
             djds = -cosd(objazi)*cosd(objalt)/grid_spacing_m
             dhtds = -sind(objalt)                                 
             dxyds =  cosd(objalt)
           else            ! ls > 0, thus after first looping
             slast = scurr ! initialized down below
             htlast = ht
           endif

!          if(idebug .eq. 1)write(6,*)'dids/djds/dhtds = ',dids,djds,dhtds

!          calculate 'scurr' as the total path traversed so far by the ray.
!          the incremental path length ('ds') will be calculated later on.
           if(objalt .gt. 15. .and. if .le. 2)then ! march by height levels
             nksteps = 2
             rkmarch = rkt - float(ls)/float(nksteps)
             if(rkmarch .le. 0.)goto 10 ! going outside the domain
             kmarch = nint(rkmarch) ! kl,kh
             htmarch = heights_1d(max(kmarch,1))
             scurr = (htt - htmarch) / (-dhtds)

             if(scurr .lt. slast)then
               write(6,*)' error s<slast 1:',s,slast,ls,htt,htmarch,rkmarch
               stop
             endif

           elseif(objalt .lt. 4.0 .and. if .eq. 1)then
             scurr = float(ls) * grid_spacing_m

           else
             scurr = float(ls) * raysteps

             if(scurr .lt. slast)then
               write(6,*)' error s<slast 2:',s,slast,ls,raysteps
               stop
             endif

           endif

           ri = rit + dids*s
           rj = rjt + djds*s

!          update the height based on the slope of the ray with a correction
!          for earth curvature. the approximation is made that the curve of
!          the earth can be given by a quadratic (parabolic) expression.
!          ht = htt + dhtds*s + (dxyds*s)**2 / (2.0*earth_radius)

!          this can be used as part of a refraction strategy
           refk = 0.179
           ht = htt + dhtds*s + (dxyds*s)**2 / ((2.0/(1.-refk))*earth_radius)

!          if(ht .le. float(htluthi) .and. ht .ge. float(htlutlo))then
           if(ht .le.  1.0*(htluthi) .and. ht .ge.  1.0*(htlutlo))then
             rk = htlut(nint(ht))
!          elseif(ht .ge. htlutlo-500.)then
!            rk = 1.0
           else
             if(idebug .eq. 1)write(6,*)' ht outside lut',ht,heights_1d(1)
             goto 10
           endif

           idlast = id; jdlast = jd; kdlast = kd
           id = nint(ri); jd = nint(rj); kd = nint(rk)

           if(id .lt. 1 .or. id .gt. ni &
         .or. jd .lt. 1 .or. jd .gt. nj &
         .or. kd .lt. 1 .or. kd .gt. nk)then ! outside domain

             if(idebug .eq. 1)write(6,2)s,ri,rj,rk,ht             
2            format('s/ri/rj/rk/ht = ',5f9.2,' outside box')
             goto 10 
           else                              ! inside domain
             if(transm_3d(id,jd,kd) .ne. r_missing_data)then
               if(id .eq. idlast .and. jd .eq. jdlast .and. kd .eq. kdlast)then
                 if(idebug .eq. 1)write(6,3)s,ri,rj,rk,ht,id,jd,kd             
3                format('s/ri/rj/rk/ht = ',5f9.2,' same march point',3i5)
                 goto 9 ! experimental speedup
               else
                 if(idebug .eq. 1)write(6,4)s,ri,rj,rk,ht,id,jd,kd             
4                format('s/ri/rj/rk/ht = ',5f9.2,' already assigned',3i5)
                 goto 9
               endif
             else ! transm_3d is missing
               if(idebug .eq. 1)write(6,5)s,ri,rj,rk,ht,id,jd,kd             
5              format('s/ri/rj/rk/ht = ',5f9.2,' new',13x,3i5)
               nnew = nnew + 1
             endif

!            if(if .eq. 3 .and. jt .eq. ni/2 .and. kt .eq. nk/2)then
!              alt_theo2 = -atand((ht-htlast) / (s-slast))
!              write(6,7)id,jd,s,obj_alt(id,jd),alt_theo2
!7             format('key ray',2i5,f9.1,2f9.4)
!            endif

!            valid trace (even if already assigned)
             illast = il; jllast = jl; kllast = kl

             il = max(min(int(ri),ni-1),1); fi = ri - il; ih=il+1
             jl = max(min(int(rj),nj-1),1); fj = rj - jl; jh=jl+1
             kl = max(min(int(rk),nk-1),1); fk = rk - kl; kh=kl+1

             if(il .eq. illast .and. jl .eq. jllast .and. kl .eq. kllast)then
               l_same_point = .true.
             else
               l_same_point = .false.
             endif

             if(ht - topo_a(id,jd) .le. 1000. .and. ihit_terrain .ne. 1 .and. .false.)then

!              interpolate to get topography at fractional grid point

!              if(il .le. 0)then
!                write(6,*)' il bounds check',il,id,ri,rj
!                stop
!              endif               

!              if(jl .le. 0)then
!                write(6,*)' jl bounds check',jl,jd,ri,rj
!                stop
!              endif               

               bi_coeff(1,1) = (1.-fi) * (1.-fj)
               bi_coeff(2,1) = fi      * (1.-fj)
               bi_coeff(1,2) = (1.-fi) *     fj 
               bi_coeff(2,2) = fi      *     fj 

               topo_bilin = sum(bi_coeff(:,:) * topo_a(il:ih,jl:jh))

               if(ht .lt. topo_bilin)then
                 ihit_terrain = 1
               endif
             endif

             if(transm_3t(id,jd,kd) .eq. 0.)then
!              transm_3d(id,jd,kd) = 0.0                      
!              transm_4d(id,jd,kd,:) = 0.0                      
             else ! trace free of terrain
               ds = s - slast
               b_alpha_last = b_alpha_new

               tri_coeff(1,1,1) = (1.-fi) * (1.-fj) * (1.- fk)
               tri_coeff(2,1,1) = fi      * (1.-fj) * (1.- fk)
               tri_coeff(1,2,1) = (1.-fi) *     fj  * (1.- fk)
               tri_coeff(1,1,2) = (1.-fi) * (1.-fj) *      fk
               tri_coeff(1,2,2) = (1.-fi) *     fj  *      fk
               tri_coeff(2,1,2) = fi      * (1.-fj) *      fk
               tri_coeff(2,2,1) = fi      *     fj  * (1.- fk)
               tri_coeff(2,2,2) = fi      *     fj  *      fk

               if(l_same_point .eqv. .false.)then
                 b_alpha_cube(:,:,:) = b_alpha_3d(il:ih,jl:jh,kl:kh)
               endif

               b_alpha_new = sum(tri_coeff(:,:,:)*b_alpha_cube(:,:,:))
               b_alpha_new = max(b_alpha_new,0.) ! clean up far edge extrapolation
               if(ls .gt. 0)then
!                b_alpha_m = 0.5 * (b_alpha_3d(id,jd,kd) + b_alpha_3d(idlast,jdlast,kdlast))
                 b_alpha_m = 0.5 * (b_alpha_new + b_alpha_last)
                 dbtau = ds * b_alpha_m
                 btau = btau + dbtau
                 albedo = btau / (1. + btau)       
               else
                 b_alpha_last = b_alpha_new
                 albedo = 0.
               endif

               if(b_alpha_new .lt. 0.)then
                 write(6,*)' error in b_alpha_new',b_alpha_new,ri,rj,rk,fi,fj,fk
                 write(6,*)'b_alpha_3d',b_alpha_3d(il:ih,jl:jh,kl:kh)
                 stop
               endif
               if(albedo .gt. 1.0 .or. albedo .lt. 0.0)then
                 write(6,*)' error in albedo ',albedo,btau,dbtau,ds,b_alpha_m,b_alpha_new,b_alpha_last,s,slast
                 stop
               endif

               transm_3d(id,jd,kd) = 1. - albedo              
               if(idebug .eq. 1)write(6,8)transm_3d(id,jd,kd),id,jd,kd
8              format('transm_3d     = ',5x,f9.4,48x,3i5)
             endif

             nfacesteps = nfacesteps + 1

           endif

9          continue

         enddo ! ls

10       continue

         idebug = 0

       enddo ! k
       enddo ! j
       enddo ! i

       write(6,*)' number new/steps assigned on this face =',nnew,nfacesteps

       ntot = ntot + nnew
       nsteps = nsteps + nfacesteps

       i4_elapsed = ishow_timer()

      enddo ! if

     else  ! solalt < twi_alt
      write(6,*)' solalt < twi_alt, raytrace not needed: ',solalt,twi_alt
      i4_elapsed = ishow_timer()

     endif ! solalt

     transm_3d(:,:,:) = transm_3d(:,:,:) * transm_3t(:,:,:)

     npts = ni*nj*nk
     fractot = float(ntot)/float(npts)
     write(6,*)' total is ',ntot, 'potential pts is',npts,' frac',fractot
     write(6,*)' nsteps is ',nsteps,'frac is',float(nsteps)/float(npts)

     arg1 = minval(transm_3d)
     arg2 = maxval(transm_3d)
     if(arg1 .lt. 0. .or. arg2 .gt. 1.0)then
       write(6,*)' warning: range of transm_3d = ',arg1,arg2
!      stop
     else
       write(6,*)' range of transm_3d = ',arg1,arg2
     endif

     i4_elapsed = ishow_timer()

     nshadow = 0
     if(solalt .ge. twi_alt)then ! daylight or early twilight
       imiss = 0
       do k = 1,nk
!        patm_k = exp(-heights_1d(k)/8000.)
         patm_k = ztopsa(heights_1d(k)) / 1013.
         topo = redp_lvl ! generic topo value
         ht_agl = heights_1d(k) - topo

!        see http://mintaka.sdsu.edu/gf/explain/atmos_refr/dip.html
         if(ht_agl .gt. 0.)then                               
           horz_dep_d = sqrt(2.0 * ht_agl / earth_radius) * 180./3.14
         else
           horz_dep_d = 0.
         endif

         obj_alt_last = r_missing_data

         do i = 1,ni; do j = 1,nj
           if(i .eq. idb .and. j .eq. jdb .and. k .eq. k)then
             iverbose = 1
           else
             iverbose = 0
           endif

           if(iverbose .eq. 1)then
             write(6,*)' here iverbose 1 ',i,j,k,heights_1d(k),ht_agl,transm_3d(i,j,k)
           endif

           if(transm_3d(i,j,k) .eq. r_missing_data)then
             imiss = imiss + 1
             if(imiss .le. 10 .or. iverbose .eq. 1)then
               write(6,*)' missing at ',i,j,k
             endif
             if(i .eq. 10 * (i/10) .and. j .eq. 10 * (j/10))then
               write(6,*)' missing at ',i,j,k
             endif
           elseif(transm_3d(i,j,k) .eq. 0.)then
             nshadow = nshadow + 1
             if(iverbose .eq. 1)then
               write(6,*)' shadow at ',i,j,k,heights_3d(i,j,k)
             endif
!          elseif(transm_3d(i,j,k) .lt. 0.)then
!            write(6,*)' error: transm_3d < 0',i,j,k,transm_3d(i,j,k)
!            stop
!          elseif(transm_3d(i,j,k) .gt. 1.)then
!            write(6,*)' error: transm_3d > 1',i,j,k,transm_3d(i,j,k)
!            stop
           else ! calculate transm_4d

!            direct illumination of the cloud is calculated here
!            indirect illumination is factored in via 'scat_frac'
             obj_alt_thr = .00 ! abs(obj_alt(i,j)) * .00
             if(abs(obj_alt(i,j) - obj_alt_last) .gt. obj_alt_thr .or. iverbose .eq. 1)then
!              ag = airmassf(cosd(90. - max(obj_alt(i,j),-3.0)),patm_k)
               ag = airmassf(90.-obj_alt(i,j), patm_k)

               if(.true.)then
                 aero_refht = redp_lvl
                 obj_alt_app = obj_alt(i,j) + refraction
                 call get_airmass(obj_alt_app,heights_3d(i,j,k),patm_k & ! i 
                                   ,aero_refht,aero_scaleht &   ! i
                                   ,earth_radius,iverbose &     ! i
                                   ,agdum,aodum,aa,refr_deg)    ! o
               else
                 aa = 0.
               endif

               obj_alt_last = obj_alt(i,j)
               refraction = refr_deg 
             endif                                                     

             obj_alt_cld = obj_alt(i,j) + horz_dep_d + refraction

             if(iverbose .eq. 1)then
               write(6,20)heights_3d(i,j,k),obj_alt(i,j),refraction,horz_dep_d,obj_alt_cld
20             format(' ht/alt/refr/hzdp/obj_alt_cld',5f9.3)
             endif

!            estimate solar extinction/reddening by rayleigh scattering
!            at this cloud altitude
             if(obj_alt_cld .lt. -0.25)then ! (early) twilight cloud lighting
!              twi_int = .1 * 10.**(+obj_alt_cld * 0.4) ! magnitudes per deg
               twi_int = 0.
               rint = twi_int
               gint = twi_int
               bint = twi_int

             else ! low daylight sun
               scat_frac = 1.00
               do ic = 1,nc
                 od_g = ag * ext_g(ic) * scat_frac

                 ext_a(ic) = (wa(ic)/.55)**(-angstrom_exp_a)
                 od_a = aa * ext_a(ic) * aod

                 od_o = (o3_du/300.) * ext_o(ic) * aodum

                 trans_c(ic) = trans(od_g + od_o + od_a)

                 if(iverbose .eq. 1)then
                   write(6,21)k,ic,obj_alt(i,j)
21                 format(' k/ic/objalt ',i4,i3,f9.2)
                   write(6,22)ag,agdum,aa,aodum,od_g,od_a,trans_c(ic)
22                 format(' ag/agdm/aa/aodm/od_g/od_a/trn',7f9.4)
                 endif
               enddo

!              fraction of solar disk (approximate)
               if(obj_alt_cld .gt. 0.25)then
                 sol_occ = 1.0
               else
                 occfrac = (obj_alt_cld - (-0.25)) / 0.5             
                 sol_occ = scurve(occfrac)             
               endif

               rint = trans_c(1) * sol_occ
               gint = trans_c(2) * sol_occ       
               bint = trans_c(3) * sol_occ              
             endif  

             transm_4d(i,j,k,1) = transm_3d(i,j,k) * rint
             transm_4d(i,j,k,2) = transm_3d(i,j,k) * gint
             transm_4d(i,j,k,3) = transm_3d(i,j,k) * bint

             if(iverbose .eq. 1)then
               write(6,23)sol_occ,rint,gint,bint,transm_4d(i,j,k,:)
23             format(' solocc/rint/gint/bint/transm_4d',f8.3,2x,3f10.5,2x,3f10.5)
             endif

           endif
         enddo ; enddo ! ij
       enddo ! k
     else ! nighttime: use red channel for sfc lighting
       do k = 1,nk       
         transm_4d(:,:,k,1) = (uprad_4d(:,:,k,2) / (2.*pi)) * sprad_to_nl(2)
       enddo ! k
       transm_3d(:,:,:) = 0.
     endif

     write(6,*)' nshadow = ',nshadow

     if(fractot .lt. 1.0)then
         write(6,*)' warning: missing points in get_cloud_rad_faces',fractot
!        stop
     endif

     where(transm_3d .eq. r_missing_data)transm_3d = 0.

     i4_elapsed = ishow_timer()

     return
     end

!    notes:
!       above threshold of 15 degrees has banding
!       at 2.2,4 degrees solalt works well - runs slow in top face
!       at 1330utc -0.45 deg looks ok
!       at 1315utc -3.19 deg good illumination - some banding
!                  turning off ag for transm_4d still has banding

     subroutine get_terrain_shadows(heights_3d,topo_a,obj_alt,obj_azi,kterr,ni,nj,nk,idb,jdb & ! i
                                   ,transm_3t,transm_2t)                                       ! o

     include 'trigd.inc'

     use mem_namelist, only: r_missing_data,earth_radius,grid_spacing_m

     real heights_3d(ni,nj,nk)
     real transm_3t(ni,nj,nk) ! terrain 3d transmission
     real transm_2t(ni,nj)    ! terrain 2d transmission
     real topo_a(ni,nj),terr_max_path_a(ni,nj)
     real obj_alt(ni,nj),obj_azi(ni,nj)
     integer kterr(ni,nj)

     write(6,*)' subroutine get_terrain_shadows'
     transm_3t = 1.0
     transm_2t = r_missing_data

     terr_max = maxval(topo_a)
     terr_min = minval(topo_a)
     terr_max_path_a = 100000. ! initialize

!    initial terrain effects
     nsteps = 0
     nloop_tot = 0
     do k = 1,nk
       i4_elapsed = ishow_timer()
       nloop = 0

       kunder = max(k-1,1)
       
       do i = 1,ni
       do j = 1,nj

         if(k .eq. kterr(i,j))then
           htstart = topo_a(i,j)
         else
           htstart = heights_3d(i,j,k)
         endif

         objalt = obj_alt(i,j)
         objazi = obj_azi(i,j)

         if(topo_a(i,j) .gt. heights_3d(i,j,k))then             ! below terrain
!          transm_3d(i,j,k) = 0.
           transm_3t(i,j,k) = 0.
         elseif(k .ge. 2 .and. transm_3t(i,j,k-1) .eq. 1.0)then ! underneath level already illuminated
           continue 
!        elseif(htstart .gt. terr_max_path(i,j) .and. objalt .gt. 0.)then
!          continue	   
         else                                                   ! above terrain and underneath level is dark
           dids = sind(objazi)*cosd(objalt)/grid_spacing_m
           djds = cosd(objazi)*cosd(objalt)/grid_spacing_m
           dhtds =  sind(objalt)                                 
           dxyds =  cosd(objalt)
           tan_suntop = tand(objalt+0.25)
           tan_sunbot = tand(objalt-0.25)
	   terr_tanalt_max = -1e10
           nloop = nloop + 1
           terr_max_path = topo_a(i,j)
           do is = 1,10000
              nsteps = nsteps + 1
              s = float(is) * grid_spacing_m
              dxy = dxyds * s
              rinew = float(i) + dids * s
              rjnew = float(j) + djds * s
              inew = int(rinew)
              jnew = int(rjnew)
              if(inew .lt. 1 .or. inew .ge. ni .or. jnew .lt. 1 .or. jnew .ge. nj)then ! outside domain
                 goto 101
	      endif

              curvat_ht = dxy**2 / (2.0*earth_radius)
              htnew = htstart + dhtds * s + curvat_ht

              terrht = topo_a(int(rinew),int(rjnew)) ! replace with bilinear interp
!             terr_max_path = max(terr_max_path,terrht)
              terr_tanalt     = ((terrht   - curvat_ht) - htstart) / dxy
              terr_tanalt_pot = ((terr_max - curvat_ht) - htstart) / dxy
              terr_tanalt_max = max(terr_tanalt_max,terr_tanalt)

              if(i .eq. idb .and. j .eq. jdb)then
                 write(6,103)is,s,inew,jnew,htstart,htnew,terrht    &
		            ,terr_tanalt,terr_tanalt_max,tan_suntop &
		            ,terr_tanalt_pot,tan_sunbot
103              format(i6,f8.0,2i5,3x,3f9.2,3x,3f9.4,3x,2f9.4)
              endif

              if(htnew .gt. terr_max .and. objalt .gt. 0.)then ! ray now above max terrain with source above horizon
                 goto 101
              endif

              if(htnew .lt. terr_min .and. objalt .lt. 0.)then ! ray now below min terrain with source below horizon
                 goto 101
              endif

              if(terr_tanalt_max .gt. tan_suntop)then          ! terrain has now completely covered light source
                 goto 101
              endif

              if(terr_tanalt_pot .lt. tan_sunbot)then          ! highest potential terrain beyond location doesn't cover light source
                 goto 101
              endif
           enddo ! is

101        continue

!          special handling at kterr
           if(k .ne. kterr(i,j))then
              if(terr_tanalt_max .gt. tan_suntop)then
                  transm_3t(i,j,k) = 0.          
              endif
           else ! k = kterr
              transm_3t(i,j,k) = 0.          
              if(terr_tanalt_max .gt. tan_suntop)then
                  transm_2t(i,j) = 0.          
              else
                  transm_2t(i,j) = 1.          
              endif
           endif

           if(i .eq. idb .and. j .eq. jdb)then
              write(6,102)k,kterr(i,j),is,inew,jnew,transm_3t(i,j,k),terr_tanalt_max,tan_suntop ! ,terr_max_path_a(i,j)
102           format(' initial terrain effects for k =',2i5,3i6,f9.3,3f9.3)
              if(k .eq. kterr(i,j))then
                 write(6,*)' transm_2t = ',transm_2t(i,j)
              endif
           endif
 
         endif ! non-trivial point

!        level is immediately above the terrain. we can consider integrating this by substituting
!        the terrain height for one of the 3d heights.
!        if(heights_3d(i,j,kunder) .lt. topo_a(i,j) .and. heights_3d(i,j,k) .ge. topo_a(i,j))then
!             transm_2t(i,j) = transm_3t(i,j,k)
!        endif

       enddo ! j
       enddo ! i
       nloop_tot = nloop_tot + nloop
       write(6,*)' nloop/ht for level ',k,nloop,heights_3d(idb,jdb,k)
     enddo ! k

     write(6,*)' nloop_tot/nsteps = ',nloop_tot,nsteps

     return
     end

     
