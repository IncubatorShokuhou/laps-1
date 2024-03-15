cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
        subroutine normalize_brightness(i4time,lat,lon,image,ni,nj
     1  ,sublat_d,sublon_d,range_m,l_national,iskip_bilin,r_missing_data
     1  ,lun,i_dir,phase_angle_d,specular_ref_angle_d,emission_angle_d
     1  ,istatus)


c***normalize a vis satellite image for solar angle.

c       s. albers          feb 94       original version
c       j. smart           jan 96       output phase angle array
c       j. smart           mar 97       added r_missing_data to argument list.
c       s. albers          may 97       no longer using 0 as missing data
c                                       value. added error trapping when
c                                       emission angle < 0.
c       s. albers          jun 97       error handling for 'iskip_bilin'
c       j. smart           mar 99       put emission angle code in subroutine.

c argument      i/o       type                    description
c --------      ---       ----    ------------------------------------------------
c i4time         i        i*4     i4time of image (1960 reference)
c lat            i        r*4  a  latitude array of image pixels (degrees)
c lon            i        r*4  a  longitude array of image pixels (degrees)
c image         i/o       r*4  a  remapped image array.
c ni             i        i*4     i dimension of image
c nj             i        i*4     j dimension of image
c sublat_d       i        r*4     latitude of satellite subpoint (degrees)
c sublon_d       i        r*4     longitude...                   (degrees)
c range_m        i        r*4     distance of spacecraft from center of earth (m)
c l_national     i        l       .true. for conus / .false. for laps
c iskip_bilin    i        i*4     >5 for conus / ~1 for laps
c                                 controls subsampling of the grid for the more
c                                 time-consuming calculations. larger values
c                                 give faster runtimes and less accuracy. 
c                                 smaller values give slower runtimes and
c                                 greater accuracy, and also risks an error
c                                 message about dimensions. in general, the 
c                                 product of the satellite data resolution and
c                                 'iskip_bilin' should lie between about 10km 
c                                 and 30km.
c lun            i        i*4     logical unit # for logging output
c i_dir          i        i*4     direction of normalization (-1,0,+1)
c phase_angle_d  o        r*4  a  phase angle (sparse array if iskip_bilin > 1)
c specular_ref_angle_d o  r*4  a  distance from specular reflection pt to sun
      include 'trigd.inc'
c istatus        o        i*4     standard status return.
        integer         ni,nj
c***parameter list variables
        real          lat(ni,nj),lon(ni,nj)
        real          sublat_d,sublon_d,range_m
        integer       i4time,istatus
        real          image(ni,nj)
        real          image_in(ni,nj)
        logical         l_national

c***local variables
        integer maxlut
        real rpd
        parameter   (maxlut=2500,         ! largest array size allowed
     1          rpd=3.1415926536/180.)

        real*8  tx,ty,tz,rx,ry,rz,
     1          satx,saty,satz,r8phase_angle_r,
     1          r8_ref_angle_r,
     1          stx_n,sty_n,stz_n,arg_dum,
     1          tx_n,ty_n,tz_n,
     1          refx,refy,refz

! note that solar_factor and phase_factor originally were declared to maxlut
        real  normfac,imgtmp,                                                                               
     1          solar_factor(ni,nj),solar_alt_d(ni,nj),
     1          xfrac,yfrac,
     1          sf_ul,sf_ur,sf_lr,sf_ll,sf_u,sf_l,s_f,weight,
     1                pf_ul,pf_ur,pf_lr,pf_ll,pf_u,pf_l,p_f,rbril,rbrih,
     1          phase_factor(ni,nj),phase_angle_d(ni,nj),
     1          sat_radius,emission_angle_d(ni,nj),
     1          specular_ref_angle_d(ni,nj)       

        real rbril_a(ni,nj),rbrih_a(ni,nj)

        integer nilut,njlut,i,j,img_i(maxlut),img_j(maxlut),
     1            ilut,jlut,ispace,jspace, ni2, nj2

        character*9 a9time

        call make_fnam_lp (i4time,a9time,istatus)
        if(istatus .ne. 1)return

        write(lun,*)' begin normalization routine for ',a9time       
        write(lun,*)' sub lat/lon/rng = ',sublat_d,sublon_d,range_m
        write(lun,*)' iskip_bilin = ',iskip_bilin

        nilut = (ni-2) / iskip_bilin + 2
        njlut = (nj-2) / iskip_bilin + 2

        write(lun,*)' ni/nj/nilut/njlut = ',ni,nj,nilut,njlut
        write(lun,*)' corner sw ',lat(1,1),lon(1,1)
        write(lun,*)' corner se ',lat(ni,1),lon(ni,1)
        write(lun,*)' corner nw ',lat(1,nj),lon(1,nj)
        write(lun,*)' corner ne ',lat(ni,nj),lon(ni,nj)

        call zero(phase_angle_d,ni,nj)
        call zero(emission_angle_d,ni,nj)
        call zero(specular_ref_angle_d,ni,nj)

        iwrite = 0

c***where's the sun?
!       use cartesian coordinates with the x-axis through the prime
!       meridian and distances in au.
        rlat = 0.
        rlon = 0.
        call solar_position(rlat,rlon,i4time,solar_alt_deg       
     1                     ,solar_dec_d,hr_angle_d)
        solar_range = 1.
        solar_sublon_d = -hr_angle_d
        rx = cosd(solar_sublon_d) * cosd(solar_dec_d) * solar_range
        ry = sind(solar_sublon_d) * cosd(solar_dec_d) * solar_range
        rz = sind(solar_dec_d)                        * solar_range

!   satellite location (in au)
        au_m = 149000000.
        sat_radius = range_m / au_m

        satx = cosd(sublon_d) * cosd(sublat_d) * sat_radius
        saty = sind(sublon_d) * cosd(sublat_d) * sat_radius
        satz = sind(sublat_d)                  * sat_radius

!       write(lun,*)
!    1 '   i    j    alt    emis   pha    pf   vis  spec      lat    lon'       

c***fill the solar brightness and phase angle arrays
        normfac=sind(58.)       ! normalized sun angle

        if(nilut .gt. maxlut .or. njlut .gt. maxlut)then
            write(lun,*)'warning: insufficient dimension for maxlut'
            write(lun,*)'maxlut = ',maxlut
            write(lun,*)'nilut = ',nilut
            write(lun,*)'njlut = ',njlut

            rskip_new = max(float(nilut)/float(maxlut)
     1                     ,float(njlut)/float(maxlut)) 
     1                    * float(iskip_bilin)
            iskip_new = int(rskip_new) + 1

            write(lun,*)'input parameter iskip_bilin currently equals'   
     1                 ,iskip_bilin
            write(lun,*)'try increasing iskip_bilin to approx '
     1                 ,iskip_new      
            write(lun,*)'secondary alternative is to increase maxlut'
     1                 ,' declared as a parameter in normalize.f'
            write(lun,*)
     1      'setting istatus to zero, returning without normalizing'

            istatus = 0
            return
        else
            write(lun,*)' dimension for maxlut',maxlut,nilut,njlut
        endif

        nj2 = nj
        ni2 = ni
        do jlut=1,njlut
         j = ((jlut-1) * iskip_bilin) + 1
         
         j = min(j,nj2)
         img_j(jlut) = j

         do ilut=1,nilut

          i = ((ilut-1) * iskip_bilin) + 1
          i = min(i,ni2)
          img_i(ilut) = i

          call solar_position(lat(i,j),lon(i,j),i4time
     1                                  ,solar_alt_d(i,j)
     1                                  ,solar_dec_d,hr_angle_d)

c   reduce and limit correction at terminator
          if(solar_alt_d(i,j).lt.8.)then
           solar_alt_d(i,j)=max(8.-(8.-solar_alt_d(i,j))*.5,3.8)
          endif

          solar_factor(ilut,jlut)=log(sind(solar_alt_d(i,j))/normfac)

c   compute emission angle (emission_angle_d = satellite angular altitude)
 
          call sat_angular_alt(sat_radius,lat(i,j),lon(i,j)
     .,satx,saty,satz,tx,ty,tz,emission_angle_d(i,j),istatus)

          if(emission_angle_d(i,j) .lt. 0.)then
              istatus = 0

              if(image(i,j) .ne. r_missing_data)then ! valid image data

                if(iwrite .le. 1)then
                  write(6,*)
     1            ' warning, emission_angle_d = ', emission_angle_d(i,j)
                  write(6,*)
     1            ' you are normalizing beyond the earths limb'     
                  write(6,*)
     1            ' check your satellite subpoint and lat/lons'
                  write(6,*)'i,j,lat(i,j),lon(i,j),sublat_d,sublon_d'        
                  write(6,*)i,j,lat(i,j),lon(i,j),sublat_d,sublon_d
                  write(6,*)'image counts is ',image(i,j)
                  iwrite = iwrite + 1
                elseif(iwrite .le. 10)then
                  write(6,*)
     1'warning, emission_angle_d < 0 (i/j/lat/lon/e): ',i,j,lat(i,j)
     1,lon(i,j),emission_angle_d(i,j)
                  iwrite=iwrite+1
                endif

              endif ! valid image pixel intensity
              emission_angle_d(i,j) = 0.0
  
          endif ! emission angle < 0 (looking beyond the limb)

c   compute phase angle
          call anglevectors(tx-satx,ty-saty,tz-satz,rx,ry,rz
     1                     ,r8phase_angle_r)
          phase_angle_d(i,j) = 180. - r8phase_angle_r / rpd

c   compute specular reflection angle
          stx_n = tx-satx
          sty_n = ty-saty
          stz_n = tz-satz
          call normalize(stx_n,sty_n,stz_n,arg_dum)

          tx_n = -tx
          ty_n = -ty
          tz_n = -tz
          call normalize(tx_n,ty_n,tz_n,arg_dum)

          call deviate_ray(-1d0,stx_n,sty_n,stz_n,tx_n,ty_n,tz_n
     1                         ,refx,refy,refz)
          call anglevectors(refx,refy,refz,rx,ry,rz,r8_ref_angle_r)
          specular_ref_angle_d(i,j) = 180. - r8_ref_angle_r / rpd


          if(i .lt. 40)then ! for testing

          phase_factor(ilut,jlut) = cosd(phase_angle_d(i,j))**6
!    1                    * sind(emission_angle_d(i,j))**1.0
!    1                    * sind(solar_alt_d(i,j))**1.0
     1                    * sind(emission_angle_d(i,j))**0.5
     1                    * sind(solar_alt_d(i,j))**0.25

          else
          phase_factor(ilut,jlut) = cosd(phase_angle_d(i,j))**6
!    1                    * sind(emission_angle_d(i,j))**1.0
!    1                    * sind(solar_alt_d(i,j))**1.0
     1                    * sind(emission_angle_d(i,j))**0.5
     1                    * sind(solar_alt_d(i,j))**0.125


          endif


!         phase_factor(ilut,jlut) = 0.    ! disable phase factor correction

!         perhaps the phase factor should be modified if the solar altitude
!         is low (< 23 deg). this would probably be designed to darken
!         the land areas at low phase angle and low solar altitude.


          if(ilut .eq. ilut/20*20 .and. jlut .eq. jlut/20*20)then

              imgtmp=image(i,j)
              if(imgtmp.eq.r_missing_data)imgtmp=-99.00
              if(iskip_bilin .ne. 1)then
                  write(lun,50)ilut,jlut,solar_alt_d(i,j)
     1                  ,emission_angle_d(i,j),phase_angle_d(i,j)
     1                  ,phase_factor(ilut,jlut),imgtmp
     1                  ,specular_ref_angle_d(i,j)
     1                  ,lat(i,j),lon(i,j)
 50               format(1x,2i5,f7.2,f7.2,f7.2,f6.2,2f6.1,2f7.2,' __')  
              endif
          endif

         enddo
        enddo

        if(iwrite.gt.0)then
           write(6,*)'warning: found ',iwrite,' pts with emission',
     +' angle < 0.0'
        endif
        write(lun,*)' normalization lookup tables complete'

c***apply the solar brightness normalization to the image
        if(iskip_bilin .eq. 1)write(lun,*)
     1   '   i    j    alt    emis   pha    '
     1  ,'pf   visi  viso  spec   lat  lon'     
        jlut=1
        jspace=img_j(jlut+1)-img_j(jlut)

        do j=1,nj

         if(j .gt. img_j(jlut+1))then
          jlut=jlut+1
          jspace=img_j(jlut+1)-img_j(jlut)
         endif

         yfrac=float(j-img_j(jlut))/jspace

         ilut=1
         ispace=img_i(ilut+1)-img_i(ilut)

         do i=1,ni

          if(i .gt. img_i(ilut+1))then
           ilut=ilut+1
           ispace=img_i(ilut+1)-img_i(ilut)
          endif

          xfrac=float(i-img_i(ilut))/ispace

c   bilinearly interpolate the normalization factor from the surrounding points.
          sf_ul=solar_factor(ilut  ,jlut  )
          sf_ur=solar_factor(ilut+1,jlut  )
          sf_lr=solar_factor(ilut+1,jlut+1)
          sf_ll=solar_factor(ilut  ,jlut+1)

          sf_u=sf_ul+(sf_ur-sf_ul)*xfrac
          sf_l=sf_ll+(sf_lr-sf_ll)*xfrac
          s_f=sf_u+(sf_l-sf_u)*yfrac

c   bilinearly interpolate the phase factor from the surrounding points.
          pf_ul=phase_factor(ilut  ,jlut  )
          pf_ur=phase_factor(ilut+1,jlut  )
          pf_lr=phase_factor(ilut+1,jlut+1)
          pf_ll=phase_factor(ilut  ,jlut+1)

          pf_u=pf_ul+(pf_ur-pf_ul)*xfrac
          pf_l=pf_ll+(pf_lr-pf_ll)*xfrac
          p_f=pf_u+(pf_l-pf_u)*yfrac

c   greater (lesser) abs. values of rbrih and rbril will brighten (darken) the
c   high and low ends, respectively.  rbril is modified to make darker land more
c   uniform in brightness.

!         original result
!         phase_const1   = 20.                                
!         ph_const2_l = 0. 
!         ph_const2_h = 0. 

!         new result
          phase_const1   = 0.                                
          ph_const2_l = 10. 
          ph_const2_h = 10. 

          rbrih=-60.*s_f   -   phase_const1*p_f
          rbrih_a(i,j) = rbrih

          if(l_national)then
            weight=18.
          else ! local type scales
            if(s_f.gt.-1.0)then
             weight=30.
            else
             weight=max(30.-12.*(-1.0-s_f),20.)
            endif
          endif

          rbril=-weight*s_f  -   phase_const1*p_f
          rbril_a(i,j) = rbril

          if(image(i,j) .ne. r_missing_data)then
            imgtmp = image(i,j)
            image_in(i,j) = image(i,j)
            if(i_dir .eq. +1)then ! regular normalization
              call stretch(68.-rbril,220.-rbrih,68.,220.,image(i,j))
              call stretch(40.+ph_const2_l*p_f, 114.+ph_const2_h*p_f,
     1                     40.,                 114.,  image(i,j))
            elseif(i_dir .eq. -1)then
              call stretch(68.,220.,68.-rbril,220.-rbrih,image(i,j))
            endif

            if(iskip_bilin .eq. 1)then ! print output visible counts
              if(i .eq. i/20*20 .and. j .eq. j/20*20)then
                write(lun,61)i,j,solar_alt_d(i,j)
     1                  ,emission_angle_d(i,j),phase_angle_d(i,j)
     1                  ,phase_factor(ilut,jlut),imgtmp,image(i,j)
     1                  ,specular_ref_angle_d(i,j)
     1                  ,lat(i,j),lon(i,j)
!               write(lun,60)image(i,j)
!60             format(1x,36x,f7.1)
 61             format(1x,2i5,f7.2,f7.2,f7.2,f6.2,3f6.1,2f7.2,' __')  
              endif
            endif

          endif

         enddo ! j
        enddo ! i

!       examine arrays
!       phamin = 9.9
!       do i = 1,ni
!       do j = 1,nj
!           phamin = min(phamin,phase_factor(i,j))
!           if(phase_factor(i,j) .le. 0.005)then
!               write(6,*)' phase factor zero ',i,j,solar_alt_d(i,j)
!    1                                             ,phase_angle_d(i,j)
!    1                                             ,phase_factor(i,j)
!    1                                             ,' ___'
!           else
!               write(6,*)' phase factor zero ',i,j,solar_alt_d(i,j)
!    1                                             ,phase_angle_d(i,j)
!    1                                             ,phase_factor(i,j)
!           endif
!       enddo ! j
!       enddo ! i

!       where(phase_factor(:,:) .lt. 0.01)
!           write(6,*)' phase_factor is zero ___'
!       end where

!       write(6,*)' phamin = ',phamin,minval(phase_factor)

        write(lun,70)minval(solar_alt_d),maxval(solar_alt_d)  
     1              ,minval(solar_factor),maxval(solar_factor)
70      format(1x,'solar alt / factor range:   ',2f8.2,3x,2f8.2)

        if(phase_const1 .gt. 0.)then
            write(lun,71)minval(phase_angle_d),maxval(phase_angle_d)
     1                  ,minval(phase_factor) * phase_const1
     1                  ,maxval(phase_factor) * phase_const1
            write(lun,72) 68.-maxval(rbril_a), 68.-minval(rbril_a)  
     1                  ,220.-maxval(rbrih_a),220.-minval(rbrih_a)
        else ! assume ph_const2 > 0.
            write(lun,71)minval(phase_angle_d),maxval(phase_angle_d)
     1                  ,minval(phase_factor) * ph_const2_l
     1                  ,maxval(phase_factor) * ph_const2_l
            write(lun,72) 40. + minval(phase_factor) * ph_const2_l, 
     1                    40. + maxval(phase_factor) * ph_const2_l,
     1                   114. + minval(phase_factor) * ph_const2_h,
     1                   114. + maxval(phase_factor) * ph_const2_h 
        endif

71      format(1x,'phase angle / rbri range:   ',2f8.2,3x,2f8.2)

72      format(1x,'stretch l h / range:        ',2f8.2,3x,2f8.2)

        write(lun,75)minval(image_in),maxval(image_in)         
     1              ,minval(image),maxval(image)          
75      format(1x,'image in/out range:         ',2f8.2,3x,2f8.2)

        istatus = 1
        write(lun,*)' normalization complete'

        return

        end
c-------------------------------------------------------------------------------
        subroutine stretch(il,ih,jl,jh,rarg)

        implicit        none

        real          a,b,il,ih,jl,jh,rarg

        a = (jh - jl) / (ih - il)
        b =  jl - il * a

        rarg = a * rarg + b
        rarg = max(min(255.,rarg),0.)

        return
        end


        subroutine deviate_ray(l1,x5,y5,z5,x3,y3,z3,x8,y8,z8)

	implicit real*8(a-z)

        c9 = x3*x5 + y3*y5 + z3*z5
        c8 = sqrt(1.d0 - l1**2 * (1.d0 - c9**2))
        m1 = c8 - l1*c9

        x8 = l1*x5 + m1*x3                    ! x dir cosine of outgoing ray
        y8 = l1*y5 + m1*y3                    ! y dir cosine of outgoing ray
        z8 = l1*z5 + m1*z3                    ! z dir cosine of outgoing ray

	return
	end

c-----------------------------------------------------------------------
c=======================================================================

      subroutine sat_angular_alt(sat_radius,lat,lon
     .,satx,saty,satz,tx,ty,tz,emission_angle_d,istatus)
c
c computes satellite altitude (degrees above horizon) for
c a given earth lat/lon. code taken from src/lib/normalize.f
c by albers, s.
c
c code put in subroutine for use in lvd satellite ingest for
c determining the polar extent of useable satellite data.
c smart, j. 3/10/99
c
      include 'trigd.inc'

      implicit none

      integer  i,j,istatus
 
      real     rpd,radius_earth_m
      parameter (rpd = 3.1415926536/180.,
     1           radius_earth_m = 6378137.)

      real  lat,lon
      real  sat_radius,au_m

!     coordinates are equatorial and relative to prime meridian
      real*8  tx,ty,tz,       ! o (earth surface)                   
     1        satx,saty,satz, ! i (satellite relative to earth center - au)
     1        stx,sty,stz,
     1        r8emission_angle_r     

      real  emission_angle_d

!     real cosd, sind

c================================================================
        au_m = 149000000.
c   compute equatorial coordinates of point on sfc of earth
        tx = cosd(lon) * cosd(lat) * radius_earth_m / au_m
        ty = sind(lon) * cosd(lat) * radius_earth_m / au_m
        tz = sind(lat)             * radius_earth_m / au_m

c       satellite relative to surface (au)
        stx = satx-tx
        sty = saty-ty
        stz = satz-tz

c   compute emission angle (emission_angle_d = satellite angular altitude)
        call anglevectors(stx,sty,stz,tx,ty,tz
     1                            ,r8emission_angle_r)
        emission_angle_d = 90. - r8emission_angle_r / rpd

        return
        end

       subroutine refl_to_albedo(reflectance,solalt,land_albedo    ! i
     1                          ,cloud_albedo)                     ! o

       include 'trigd.inc'

       real land_albedo

       solalt_eff = max(solalt,6.0)

       reflectance_land = land_albedo * sind(solalt_eff)
       reflectance_air  = 0.11

!      cloud albedo = 1 if reflectance = sind(solalt_eff)       
!      cloud albedo = 0 if reflectance = reflectance_land+reflectance_air

       zeropoint = reflectance_land + reflectance_air
       onepoint = sind(solalt_eff)

       cloud_albedo = (reflectance-zeropoint) / (onepoint-zeropoint)
       cloud_albedo = min(max(cloud_albedo,0.),1.)
       
       return
       end
