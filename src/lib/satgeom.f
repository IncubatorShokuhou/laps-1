        subroutine satgeom(i4time,lat,lon,ni,nj
     1  ,sublat_d_a,sublon_d_a,range_m,r_missing_data,phase_angle_d
     1  ,specular_ref_angle_d,emission_angle_d,azimuth_d,istatus)


        include 'trigd.inc'

        angdif(x,y)=mod(x-y+540.,360.)-180.

c istatus        o        i*4     standard status return.
        integer         ni,nj
c***parameter list variables
        real          lat(ni,nj),lon(ni,nj)
        real          sublat_d_a(ni,nj),sublon_d_a(ni,nj)
        real          sublat_d,sublon_d,range_m ! sat from earth center
        integer       i4time,istatus

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
     1          solar_alt_d(ni,nj),
     1          xfrac,yfrac,
     1          sf_ul,sf_ur,sf_lr,sf_ll,sf_u,sf_l,s_f,weight,
     1                pf_ul,pf_ur,pf_lr,pf_ll,pf_u,pf_l,p_f,rbril,rbrih,
     1          phase_factor(ni,nj),phase_angle_d(ni,nj),
     1          sat_radius,emission_angle_d(ni,nj),
     1          azimuth_d(ni,nj),
     1          specular_ref_angle_d(ni,nj)

        real rbril_a(ni,nj),rbrih_a(ni,nj)

        integer nilut,njlut,i,j,img_i(maxlut),img_j(maxlut),
     1            ilut,jlut,ispace,jspace, ni2, nj2

        character*9 a9time

        call make_fnam_lp (i4time,a9time,istatus)
        if(istatus .ne. 1)return

        lun = 6
        write(lun,*)' begin satgeom at ',a9time

        call zero(phase_angle_d,ni,nj)
        call zero(emission_angle_d,ni,nj)
        call zero(specular_ref_angle_d,ni,nj)

        iwrite = 0

c***where's the sun?
!       use cartesian coordinates with the x-axis through the prime
!       meridian and distances in au (geocentric equatorial).
        rlat = 0.
        rlon = 0.
        call solar_position(rlat,rlon,i4time,solar_alt_deg
     1                     ,solar_dec_d,hr_angle_d)
        solar_range = 1.
        solar_sublon_d = -hr_angle_d
        rx = cosd(solar_sublon_d) * cosd(solar_dec_d) * solar_range
        ry = sind(solar_sublon_d) * cosd(solar_dec_d) * solar_range
        rz = sind(solar_dec_d)                        * solar_range

!   satellite location (in au - geocentric equatorial)
        au_m = 149000000.
        sat_radius = range_m / au_m

        write(lun,*)'    i    j    alt    emis   pha    pf   vis  spec'

c***fill the solar brightness and phase angle arrays
        normfac=sind(58.)       ! normalized sun angle

        do j = 1,nj
         do i = 1,ni

          sublat_d = sublat_d_a(i,j)
          sublon_d = sublon_d_a(i,j)

          satx = cosd(sublon_d) * cosd(sublat_d) * sat_radius
          saty = sind(sublon_d) * cosd(sublat_d) * sat_radius
          satz = sind(sublat_d)                  * sat_radius

c   compute emission angle (emission_angle_d = satellite angular altitude)

          call sat_angular_alt(sat_radius,lat(i,j),lon(i,j)
     .,satx,saty,satz,tx,ty,tz,emission_angle_d(i,j),istatus)

          if(emission_angle_d(i,j) .lt. 0.)then
              istatus = 0
              emission_angle_d(i,j) = 0.0
          endif ! emission angle < 0 (looking beyond the limb)

!         topocentric equatorial vector of satellite relative to gridpoint
          dx=satx-tx
          dy=saty-ty
          dz=satz-tz

          intvl = max(ni/21,1)
          if(j .eq. nj/2 .and.
     1         (i .eq. (i/intvl)*intvl .or.
     1         (emission_angle_d(i,j) .gt. 0.0
     1          .and. emission_angle_d(i,j) .le. 10.0)) )then
              idebug = 1
          else
              idebug = 0
          endif

          if(idebug .eq. 1)then
              write(6,*)
              write(6,11)lon(i,j),sublon_d,dx,dy,dz
     1                  ,emission_angle_d(i,j)
11            format('/lon/sub/dx/dy/dz/emis',6f11.4)            
          endif

!         rotate this vector around z axis to get local cartesian coordinates
          call rotate_z(dx,dy,dz,-lon(i,j))

!         convert cartesian coordinates to dec and ha
          call xyz_to_polar_d(dx,dy,dz,dec,angle,r)

          if(angle .gt. 180.)then
              angle = angle - 360.
          endif

          ha = -angle

          if(idebug .eq. 1)then
              write(6,*)'rotated dx/dy/dz,dec,ha',dx,dy,dz,dec,ha
          endif

!         convert dec and ha to alt/az
          call equ_to_altaz_d(dec,ha,lat(i,j),alt,azimuth_d(i,j))

          if(idebug .eq. 1)then
              write(6,*)'alt/az',alt,azimuth_d(i,j)
          endif

          goto500

          call solar_position(lat(i,j),lon(i,j),i4time
     1                                  ,solar_alt_d(i,j)
     1                                  ,solar_dec_d,hr_angle_d)

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

500       continue

         enddo
        enddo

c-----------------------------------------------------------------------
c=======================================================================

        return
        end

