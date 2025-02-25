
       
        subroutine get_starglow(i4time,alt_a,azi_a,elong_a,minalt,maxalt,minazi,maxazi,rlat,rlon,alt_scale,azi_scale,horz_dep,l_zod,glwmid,glow_stars)

        include 'trigd.inc'

        use mem_allsky, only: nc
        include 'rad.inc'

!       http://arxiv.org/pdf/astro-ph/9706111.pdf
!       http://adsabs.harvard.edu/abs/1986pasp...98..364g
!       1 nl = 1.31e6  photons (550nm) sec**-1 cm**-2 sterad**-1
!       1 nl = 1.31e10 photons (550nm) sec**-1  m**-2 sterad**-1
!       1 watt = 2.76877e18 photons (550nm) / sec
!       1 nl = 4.7313e-9 watts m**-2 sterad**-1 (550nm)
!       wm2sr_to_nl(x) = 2.113e8 * x     ! (550nm)
!       wm2sr_to_nl(x) = 2145.707e9 * x
        wm2srum_to_s10(x) = x / 1.261e-8 ! 550nm
        s10_to_wm2srum(x) = x * 1.261e-8 ! 550nm
        wm2_to_wm2nm(x)     = 0.00136 * x ! solar 550nm approx
        wm2sr_to_wm2srnm(x) = 0.00136 * x ! solar 550nm approx
        wm2sr_to_wm2srum(x) = x * 1.36 ! solar 550nm approx
        wm2srum_to_wm2sr(x) = x / 1.36 ! solar 550nm approx
!       wm2srum_to_nl(x) = wm2sr_to_nl(wm2srum_to_wm2sr(x)) ! sol 550nm
!       s10_to_nl(x) = wm2srum_to_nl(s10_to_wm2srum(x)) ! solar 550nm

        real alt_a(minalt:maxalt,minazi:maxazi)
        real azi_a(minalt:maxalt,minazi:maxazi)
        real elong_a(minalt:maxalt,minazi:maxazi)
        real glow_stars(3,minalt:maxalt,minazi:maxazi) ! log nl
        real rad_stars(3,minalt:maxalt,minazi:maxazi) ! nl

!       local arrays
        real dec_a(minalt:maxalt,minazi:maxazi)
        real ha_a(minalt:maxalt,minazi:maxazi)
        real ra_a(minalt:maxalt,minazi:maxazi)
        real gallat_a(minalt:maxalt,minazi:maxazi)
        real gallon_a(minalt:maxalt,minazi:maxazi)
        real helioeclipticlat_a(minalt:maxalt,minazi:maxazi)
        real helioeclipticlon_a(minalt:maxalt,minazi:maxazi)

        character*20 c_conv

        parameter (nstars = 120000)
        real dec_d(nstars),ra_d(nstars),mag_stars(nstars),bmv(nstars)
        real alt_stars(nstars),azi_stars(nstars),ext_mag(nstars),lst_deg
        real angdif
        real*8 dangdif,jed,r8lon,lst,has(nstars),phi,als,azs,ras,decr,xx,yy
        real*8 dalt,dazi,dha,ddec,dra,sol_meanlon,sol_meananom

        integer*8 n_star_loops

        character*20 starnames(nstars)
        logical l_zod, l_psf /.true./

        dangdif(xx,yy)=dmod(xx-yy+9.4247779607694d0,6.2831853071796d0)-3.1415926535897932d0
        angdifd(x,y)=mod(x-y+540.,360.)-180.
!       addmags(a,b)=log10(10.**(-a*0.4) + 10.**(-b*0.4)) * (-2.5)
        addlogs(x,y) = log10(10.**x + 10.**y)

!       according to allen's "astrophysical quantities" page 134, sky brightness is
!       given in terms of the number of 10th magnitude stars per square degree in the visual band:

!       air glow at zenith....................40                 (6th mag per square degree)
!       zodiacal light........................100                (5th mag per square degree)
!       faint stars...........................320 on milky way
!                  ........................... 95 off milky way
!       diffuse galactic light................. 20
!       total sky brightness at zenith........290                (130 is 4.5 mag per square deg)

        write(6,*)' subroutine get_starglow...'
        write(6,*)' lat/lon = ',rlat,rlon             
        write(6,*)' l_zod = ',l_zod

        phi = rlat * rpd
        r8lon = rlon          

        call i4time_to_jd(i4time,jed,istatus)
        call sidereal_time(jed,r8lon,lst)

        lst_deg = lst / rpd
        write(6,*)' sidereal time (deg) = ',lst_deg               

!       obtain stars data
        call read_hyg(nstars,ns,dec_d,ra_d,mag_stars,bmv,starnames,istat_cat)
        if(istat_cat .eq. 0)then
          call read_bsc(nstars,ns,dec_d,ra_d,mag_stars,bmv,starnames)
        endif

        do is = 1,ns        
          ras = ra_d(is) * rpd
          decr=dec_d(is) * rpd

          has(is)=dangdif(lst,ras)
          call equ_to_altaz_r(decr,has(is),phi,als,azs)

          alt_stars(is) = als / rpd
          azi_stars(is) = azs / rpd
          if(is .le. 10)then
            write(6,5)is,starnames(is),dec_d(is),ra_d(is),has(is)/rpd,alt_stars(is),azi_stars(is),mag_stars(is),bmv(is)
5           format('dec/ra/ha/al/az/mg/bmv',i5,1x,a20,1x,f9.1,4f9.3,3f7.1)
          endif
        enddo ! is

        i4_elapsed = ishow_timer()

!       obtain planets data (placed at end of star arrays)
        do iobj = 2,6
            ns = ns + 1
            write(6,*)' call sun_planet for ',iobj,ns 
            call sun_planet(i4time,iobj,rlat,rlon,dec_d(ns),ra_d(ns),alt_stars(ns),azi_stars(ns),elgms_r4,mag_stars(ns))
            has(ns) = 0. ! this isn't being used since we already have alt/azi
            starnames(ns) = 'planet'
            bmv(ns) = 0.
        enddo ! iobj

        i4_elapsed = ishow_timer()

        do is = 1,ns        
          patm = 1.0     

          if(alt_stars(is) .gt. 0.0)then
            call calc_extinction(90.          ,patm,airmass,zenext)
            call calc_extinction(alt_stars(is),patm,airmass,totexto)
            ext_mag(is) = 0.0 ! totexto - zenext 
          else
            ext_mag(is) = 0.0
          endif

          if(is .le. 100 .or. (is .gt. 100 .and. is .ge. ns-4))then
            write(6,7)is,starnames(is),dec_d(is),ra_d(is),has(is)/rpd,alt_stars(is),azi_stars(is),mag_stars(is),bmv(is),ext_mag(is)
7           format('dec/ra/ha/al/az/mg/bmv/ext',i5,1x,a20,1x,f9.1,4f9.3,3f7.1)
          endif
        enddo ! is

!       zodiacal light
!       http://www.ing.iac.es/astronomy/observing/conditions/skybr/skybr.html
!       http://aas.aanda.org/articles/aas/pdf/1998/01/ds1449.pdf
!       http://download.springer.com/static/pdf/80/art%253a10.1186%252fbf03352136.pdf?auth66=1427760212_d8fc6d84b219c957974dd2b221bc226a&ext=.pdf
!       http://adsabs.harvard.edu/full/1975a%26a....38..405d

        glow_stars = 1.0 ! initialize to cosmic background level

        sqarcsec_per_sqdeg = 3600.**2

!       coordinate conversions
        nalt = (maxalt-minalt) + 1
        nazi = (maxazi-minazi) + 1
        c_conv = 'altaz2decra'
        call angcoord_2d(c_conv,nalt,nazi,lst_deg,arg2,alt_a,azi_a,rlat,dec_a,ra_a)

        c_conv = 'decra2galactic'
        call angcoord_2d(c_conv,nalt,nazi,arg1,arg2,dec_a,ra_a,argphi,gallat_a,gallon_a)

!       http://en.wikipedia.org/wiki/position_of_the_sun
        sol_meanlon  = 280.460 + 0.9856474 * (jed-2451545.0d0)
        sol_meananom = 357.528 + 0.9856003 * (jed-2451545.0d0)
        sollon = sol_meanlon + 1.915 * sind(sol_meananom) + .02 * sind(2. * sol_meananom)
        sollon = mod(sollon,360.)
        write(6,*)' sollon (computed) ',sollon
        c_conv = 'decra2helioecliptic'
        call angcoord_2d(c_conv,nalt,nazi,sollon,arg2,dec_a,ra_a,argphi,helioeclipticlat_a,helioeclipticlon_a)

        if(l_zod .eqv. .true.)then
          write(6,*)' computing milky way / zodiacal light / corona',horz_dep
        else
          write(6,*)' skipping milky way / zodiacal light / corona'
        endif

        glwmid_mag = (glwmid - 2.925) * 0.4
        pix_mag = log10((.25/alt_scale)**2) * 2.5
        ref_mag = 1.3 + pix_mag - glwmid_mag ! higher value has larger fainter stars

!       write(6,*)'star area/radius (deg/vpix)', size_glow_sqdg_st,strad_dg,strad_dg/alt_scale
        write(6,*)'glwmid_mag / pix_mag / ref_mag',glwmid_mag,pix_mag,ref_mag

        n_anti_calls = 0
        n_star_loops = 0

        do ialt = minalt,maxalt
        do jazi = minazi,maxazi

           altg = alt_a(ialt,jazi)
           azig = azi_a(ialt,jazi)

           if(ialt .eq. (ialt/10)*10 .and. jazi .eq. (ialt/10)*10)then
             write(6,*)' altg,azig/elong',ialt,jazi,altg,azig,elong_a(ialt,jazi) ! test
           endif
          
           if(altg .ge. -horz_dep)then

             if(l_zod .eqv. .true.)then

!             milky way (mainly from stars from 6-16 magnitude)
              explatterm = (abs(gallat_a(ialt,jazi)/20.))**1.5
              explonterm = (angdifd(gallon_a(ialt,jazi),0.)/140.)**2
              s10gal = 25. + 250. * exp(-explatterm) * exp(-explonterm) 
              if(s10gal .gt. 0.)then
                vgal = s10_to_magsecsq(s10gal)
                glow_gal = log10(v_to_b(vgal))
              else
                glow_gal = 0.
              endif

!             zodiacal light s10 units 140 - 90 sin(beta) ecliptic latitude
              elat = helioeclipticlat_a(ialt,jazi)
              elon = helioeclipticlon_a(ialt,jazi)
              if(elon .gt. 180.)elon = 360. - elon
              s10zod_pole = 60. ! high ecliptic latitude
              s10geg_ecliptic =   40. * cosd((elon-180.)/2.)**150.0
!             s10geg_ecliptic =  400. * cosd((elon-180.)/2.)**150.0 ! test
              geg_factor = 1. + (40./140.) * cosd((elong_a(ialt,jazi)-180.)/2.)**150.0
              s10zod_ecliptic =  900. * cosd((elon)/2.0001)**7.7 + 140. ! + s10geg_ecliptic
!             s10zod_ecliptic = 4000. * cosd((elon)/2.)**16.0 + 150. +
!             s10geg_ecliptic
!             note the gegenschein of s10 = 40 can be added near
!             helioecliptic 180. to yield a total of 180 s10 units.
              eclp = 1.9 - 0.9*abs(elon/180.)
              eclp = 1.0  ! control ecliptic sharpness
              ecllatterm = abs(sind(elat))**eclp
              s10zod = 10.**(log10(s10zod_pole) * ecllatterm + &
                             log10(s10zod_ecliptic) * (1.-ecllatterm))
              s10zod = s10zod * geg_factor
!             s10zod =            (s10zod_pole) * ecllatterm + &
!                                 (s10zod_ecliptic) * (1.-ecllatterm) 
!             s10zod = s10zod*10. ! debug for visibility

!             f/k corona section (inner)
!             http://arxiv.org/pdf/0909.1722.pdf (far corona in fig 1.)
!             https://ase.tufts.edu/cosmos/print_images.asp?id=28
!             https://www.terrapub.co.jp/journals/eps/pdf/5006_07/50060493.pdf
!             "interplanetary dust" edited by b. grun
              sqarcsec_per_sqdeg = 3600.**2
              size_glow_sqdg_sm = 0.2    ! sun/moon area           
              delta_mag = log10(size_glow_sqdg_sm*sqarcsec_per_sqdeg)*2.5
              distr = sqrt(elat**2 + elon**2)
              distr2 = max(distr,0.40) ! account for pixel size
              srad = distr2/0.25
              spowerk =  5.7 +       srad  * 2.0
              spowerf =  7.5 + log10(srad) * 2.0
!             spowerf =  8.0 + log10(srad) * 2.0
!             spowerf2 = 8.3 +       srad  * 0.16
              smag_eff = -26.7 + 2.5 * min(spowerk,spowerf)
              rmag_per_sqarcsec = smag_eff + delta_mag
              elong = min(sqrt(elon**2 + elat**2),180.)
              coeff = 3.5 - 2.5*cosd(elong/2)**180.
              if(elong .gt. 0.)then
!               pcoeff = 1.0 - elon**2/(elon**2 + elat**2)
!               pcoeff = abs(elat) / sqrt(elon**2 + elat**2)
!               the last exponent controls ecliptic sharpness
                pcoeff = (abs(elat) / sqrt(elon**2 + elat**2))**1.3 
              else
                pcoeff = 0.
              endif
              fexp = -2.5*(1.-pcoeff) + (-2.8*pcoeff)
              f_wm2srum = coeff * srad**fexp
              f_s10 = wm2srum_to_s10(f_wm2srum)
!             f_nl  = wm2srum_to_nl(f_wm2srum)

!             merge original zodiacal light with corona
              s10zod2 = max(s10zod,f_s10)

              sbu = s10zod2 * 1.25
              vzod = s10_to_magsecsq(s10zod2)
              glow_zod = log10(v_to_b(vzod))

              do ic = 1,3
                glow_stars(ic,ialt,jazi) = log10(10.**glow_stars(ic,ialt,jazi) + 10.**glow_zod + 10.**glow_gal)             
!               glow_stars(ic,ialt,jazi) = log10(10.**glow_stars(ic,ialt,jazi))             
              enddo ! ic

              iwrite = 0

              if((azig .eq. 0.0 .or. azig .eq. 90. .or. azig .eq. 180. .or. azig .eq. 270.) .and. &
                  (altg .eq. 0.0 .or. altg .eq. 10. .or. altg .eq. 20. .or. altg .eq. 30. .or. altg .eq. 40. .or. altg .eq. 60. .or. altg .eq. 90.) &
                                                                                                                            )then
                iwrite = 1
              endif

              if(elong_a(ialt,jazi) .gt. (180. - alt_scale))then
                iwrite = 2
              endif

              if(iwrite .ge. 1)then
                i4_elapsed = ishow_timer()
                write(6,*)
                if(iwrite .eq. 2)write(6,*)' gegenschein point'

                write(6,71)altg,azig,dec_a(ialt,jazi),ra_a(ialt,jazi) &
                          ,gallat_a(ialt,jazi),gallon_a(ialt,jazi) &
                          ,helioeclipticlat_a(ialt,jazi),helioeclipticlon_a(ialt,jazi)
71              format(' alt-azi/dec-ra/gallat-lon/helioecllat-lon',4(2x,2f8.2))
                write(6,*)'angdif/explonterm/s10gal',angdifd(gallon_a(ialt,jazi),0.),explonterm,s10gal
                write(6,72)elat,elon,srad,spowerk,spowerf,f_wm2srum,s10zod,f_s10,s10zod2,f_wm2srum*1e8,coeff,fexp
72              format(' elat/elon',5f9.2,f10.6,4f10.0,2f9.2)

                write(6,*)'s10geg_ecl/s10zod_ecl/s10zod',s10geg_ecliptic,s10zod_ecliptic,s10zod
                write(6,*)'s10zod/s10gal/glowstars',s10zod2,s10gal,glow_stars(2,ialt,jazi)
              endif

             endif ! l_zod
           endif ! alt > -20.
        enddo ! jazi
        enddo ! ialt

        do is = 1,ns        
          ialt_star = nint(alt_stars(is)/alt_scale)
          jazi_star = nint(azi_stars(is)/azi_scale)

          minalt_star = max(ialt_star - 16,minalt)
          maxalt_star = min(ialt_star + 16,maxalt)

          minazi_star = max(jazi_star - 16,minazi)
          maxazi_star = min(jazi_star + 16,maxazi)

          do ialt = minalt_star,maxalt_star
          do jazi = minazi_star,maxazi_star

            altg = alt_a(ialt,jazi)
            azig = azi_a(ialt,jazi)

            if(altg .ge. -horz_dep)then

!             threshold box to consider a star
              alt_cos = min(altg,89.)

                n_star_loops = n_star_loops + 1

!               consider varying the size with the magnitude
!               size_glow_sqdg_st = 0.1  ! final polar kernel size?
!               size_glow_sqdg_st = alt_scale**2
!               strad_dg = sqrt(size_glow_sqdg_st) / 3.14159
                
                strad_dg = alt_scale * 0.4

                if(mag_stars(is) .le. ref_mag)then
!                  ratio_mid = 10.**(glwmid - 2.925)
                   ratio_bri = 10.**((ref_mag-mag_stars(is))*0.4)
                   strad_dg = strad_dg * sqrt(ratio_bri)
                endif

!               if the psf is gaussian then we can note the area under the
!               curve is 2 * pi * amplitude * sigmax * sigmay                
                if(.not. l_psf)then
                   strad_dg = min(strad_dg,0.22)
                   thr_dist = 1.5
                else
                   strad_dg = min(strad_dg,0.33)
                   thr_dist = 3.5
                endif
                size_glow_sqdg_st = 3.14159 * strad_dg**2

!               radii of ellipse edge
!               alt_dist = alt_scale / 2.0          * 1.3
!               azi_dist = alt_dist / cosd(alt_cos) * 1.3   
                alt_dist = strad_dg
                azi_dist = strad_dg / cosd(alt_cos)

!               place star in center of grid box
                alts = nint(alt_stars(is)/alt_scale) * alt_scale
                azis = nint(azi_stars(is)/azi_scale) * azi_scale
                if(abs(alts-altg) .le. thr_dist*alt_dist .and. abs(azis-azig) .le. thr_dist*azi_dist)then
                    if(is .le. 20 .or. is .ge. ns-4 .or. is .eq. 272 .or. mag_stars(is) .le. -1.0)then
                      diste = sqrt( ((alts-altg)/alt_dist)**2 + ((azis-azig)/azi_dist)**2) 
                      write(6,81)is,mag_stars(is),strad_dg,size_glow_sqdg_st,ialt,jazi,alt_dist,azi_dist,alts-altg,azis-azig,diste
81                    format(/' sanity 1',i5,f9.1,2f9.3,2i5,5f9.1)
                      iverb = 1
                    else
                      iverb = 0
                    endif

                    radius_pix = strad_dg/alt_scale ! 0.25 / alt_scale ! 0.25 is radius_deg of moon
                    ricen = (azis-azig)/azi_scale
                    rjcen = (alts-altg)/alt_scale
                    aspect_ratio = 1. / cosd(min(abs(alts),89.))
                    n_anti_calls = n_anti_calls + 1
                    if(l_psf)then
                      dist_dg = sqrt((ricen/aspect_ratio)**2 + rjcen**2) * alt_scale
                      dist_dg = max(dist_dg,alt_scale*0.3) ! rough gaussian average of central pixel
                      frac_lit1 = 0.5 * exp(-((dist_dg/strad_dg)**2))
                    else
                      call antialias_ellipse(radius_pix,ricen,rjcen,aspect_ratio,frac_lit1,0,0,iverb)
                    endif

                    if(frac_lit1 .gt. 0.)then
                        delta_mag = log10(size_glow_sqdg_st*sqarcsec_per_sqdeg/frac_lit1)*2.5
                        rmag_per_sqarcsec = mag_stars(is) + delta_mag                  

!                       convert to log nanolamberts
                        glow_nl = log10(v_to_b(rmag_per_sqarcsec))
                    else
                        delta_mag = 99.
                        rmag_per_sqarcsec = 99.
                        glow_nl = 0.
                    endif

                    do ic = 1,3
                      if(ic .eq. 1)then
                        glow_delt = (bmv(is)-0.6) * 0.7 * 0.4
                      elseif(ic .eq. 2)then
                        glow_delt = 0.
                      elseif(ic .eq. 3)then
                        glow_delt = -(bmv(is)-0.6) * 0.4
                      endif

                      glow_nl = glow_nl + glow_delt

                      glow_stars(ic,ialt,jazi) = log10(10.**glow_stars(ic,ialt,jazi) + 10.**glow_nl)             

                      if(iverb .eq. 1)then
                        if(ic .eq. 2)then
!                         write(6,*)' sanity 2',is,alt_stars(is),alts,azis
                          write(6,91)is,mag_stars(is),bmv(is),rmag_per_sqarcsec,delta_mag,glow_nl,glow_stars(ic,ialt,jazi),frac_lit1
91                        format( &
                          ' mag/bmv/rmag_per_sqsec/dmag/glow_nl/glow_stars/flt1 = ',i4,2f8.2,4f10.3,f8.3)             
                        endif
                      endif

                    enddo ! ic
                endif ! within star kernel
            endif ! alt > -horz_dep
          enddo ! jazi
          enddo ! ialt
        enddo ! is

        write(6,*)' number of antialias calls is ',n_anti_calls

        rad_stars(:,:,:) = 10.**glow_stars(:,:,:)
 
        return
        end 

        subroutine read_stars(nstars,ns,dec_d,ra_d,mag_stars,starnames)

        real dec_d(nstars),ra_d(nstars),mag_stars(nstars)
        character*20 starnames(nstars)

        character*150 static_dir,filename
        character*120 cline

        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/stars.dat'          

        is = 0
        open(51,file=filename,status='old')
1       read(51,2,err=3,end=9)cline
2       format(a)
3       is = is + 1
!       write(6,*)cline

        if(is .le. 0)then
            do ic = 1,80
                write(6,*)' char ',ic,cline(ic:ic)
            enddo ! ic
        endif

        read(cline,4,err=5)ih,im,dec_d(is),mag_stars(is)
4       format(49x,i2,1x,i2,1x,f5.0,28x,f5.0)
5       continue
!       read(cline(66:70),*,err=6)mag_stars(is)
6       continue

!       convert coordinates
        ra_d(is) = float(ih)*15. + float(im)/4.

        starnames(is) = cline(30:49) 
!       write(6,*)' name is: ',starnames(is)
!       write(6,*)' mag string is: ',cline(91:95)
!       write(6,*)'ih/im',ih,im,dec_d(is),ra_d(is),mag_stars(is)

        goto 1

9       close(51)

        ns = is

        write(6,*)' read_stars completing with # stars of ',ns

        return
        end

        subroutine read_bsc(nstars,ns,dec_d,ra_d,mag_stars,bmv,starnames)

!       http://tdc-www.harvard.edu/catalogs/bsc5.html

        real dec_d(nstars),ra_d(nstars),mag_stars(nstars),bmv(nstars)
        character*20 starnames(nstars)

        character*150 static_dir,filename
        character*120 cline
        character*1 c1_dec

        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/bsc5.dat'          

        is = 0
        open(51,file=filename,status='old')
1       read(51,2,err=3,end=9)cline
2       format(a)
3       is = is + 1
!       write(6,*)cline

        if(is .le. 0)then
            do ic = 1,80
                write(6,*)' char ',ic,cline(ic:ic)
            enddo ! ic
        endif

        read(cline,4,err=5)ih,im,c1_dec,id,idm,mag_stars(is),bmv(is)
4       format(75x,i2,i2,4x,a1,i2,i2,14x,f5.0,2x,f5.2)
5       continue
!       read(cline(66:70),*,err=6)mag_stars(is)
6       continue

!       convert coordinates
        if(c1_dec .eq. '+')then
            rsign = 1.0
        else
            rsign = -1.0
        endif
        dec_d(is) = rsign * (float(id) + float(idm)/60.)
        ra_d(is) = float(ih)*15. + float(im)/4.

        starnames(is) = cline(5:14) 

        if(dec_d(is) .eq. 0. .and. ra_d(is) .eq. 0)then ! qc
            is = is - 1
            iqc = 0
        else
            iqc = 1
        endif

        if(mag_stars(is) .lt. 2.0 .and. iqc .eq. 1)then ! magnitude limit
!           write(6,*)' name is: ',
!           write(6,*)' mag string is: ',cline(91:95)
            write(6,13)is,starnames(is),ih,im,dec_d(is),ra_d(is),mag_stars(is),bmv(is)
13          format(i4,1x,a10,' ih/im',2i4,2f7.2,' mag/bmv',2f7.1)
        endif

        if(mag_stars(is) .gt. 99.0)then ! magnitude limit
            is = is - 1
        endif

        goto 1

9       close(51)

        ns = is

        write(6,*)' read_bsc catalog completing with # stars of ',ns

        return
        end

        subroutine read_hyg(nstars,ns,dec_d,ra_d,mag_stars,bmv,starnames,istatus)

!       http://www.projectrho.com/public_html/starmaps/catalogues.php
!       http://astronexus.com/node/34

        real dec_d(nstars),ra_d(nstars),mag_stars(nstars),bmv(nstars)
        character*20 starnames(nstars)

        character*150 static_dir,filename
        character*250 cline
        character*1 c1_dec
        character*20 cdum5,cdum6,cdum7,cdum16

        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/hygdata_v3.csv'          

        is = 0
        open(51,file=filename,status='old',err=990)

        write(6,*)' reading hyg catalog'
        
        read(51,*) ! header line
        read(51,*) ! sol

1       continue
        cdum7 = ''
!       read(51,*,err=3,end=9)cline
        read(51,*,err=3,end=9)dum1,dum2,dum3,dum4,cdum5,cdum6,cdum7,dum8,dum9,dum10,dum11,dum12,dum13,dum14,dum15,cdum16,dum17,dum18
2       format(a)
3       is = is + 1

        dec_d(is) = dum9
        ra_d(is) = dum8 * 15.

        starnames(is) = cdum7
        mag_stars(is) = dum14
        bmv(is) = dum17

        if(dec_d(is) .eq. 0. .or. ra_d(is) .eq. 0.)then ! qc
            is = is - 1
            iqc = 0
        elseif(mag_stars(is) .eq. 0. .or. mag_stars(is) .lt. -2.0)then ! qc
            is = is - 1
            iqc = 0
        elseif(mag_stars(is) .gt. 1e37)then ! magnitude limit
            is = is - 1
            iqc = 0
        else
            iqc = 1
        endif

!       check failure to read all the variables for betelgeuse at line 27919
        if(is .le. 100)then
            write(6,*)' dum8/9/14 ',nint(dum1),is,dum8,dum9,dum14,iqc
        endif

        if(mag_stars(is) .lt. 2.0 .and. iqc .eq. 1)then ! magnitude limit
!           write(6,*)' name is: ',
!           write(6,*)' mag string is: ',cline(91:95)
            write(6,*)' dum8/9/14 ',nint(dum1),is,dum8,dum9,dum14,iqc
            write(6,13)is,starnames(is),dec_d(is),ra_d(is),mag_stars(is),bmv(is)
13          format(i6,1x,a10,2f7.2,' mag/bmv',2f7.1)
        endif

        goto 1

9       close(51)

        ns = is

        write(6,*)' read_hyg catalog completing with # stars of ',ns
        istatus = 1
        return

990     write(6,*)' read_hyg catalog unavailable'
        istatus = 0
        
        end
       
        subroutine get_glow_obj(i4time,alt_a,azi_a,minalt,maxalt,minazi,maxazi &
                               ,alt_scale,azi_scale &
                               ,htmsl,patm & 
                               ,alt_obj_in,azi_obj_in,mag_obj,l_obsc &
                               ,l_phase,rill,va &
                               ,alt_obj2,azi_obj2,emag &
                               ,diam_deg,horz_dep,glow_obj)

        use mem_namelist, only: r_missing_data,earth_radius,aero_scaleht,redp_lvl
        use cloud_rad, only: ghi_zen_toa
        use mem_allsky, only: day_int0

        include 'trigd.inc'

        parameter (pi = 3.14159265)

        real alt_obj_in ! i (true altitude uncorrected for refraction)
        real alt_a(minalt:maxalt,minazi:maxazi)
        real azi_a(minalt:maxalt,minazi:maxazi)
        real glow_obj(minalt:maxalt,minazi:maxazi) ! log nl

        real ext_mag,lst_deg,mag_obj
        real*8 angdif,jed,lst,has,als,azs,ras,x,y,decr

        character*20 starnames

        logical l_obsc, l_phase
        real maginterp,magtobri0,bri0tomag
        magtobri0(a) = 10.**(-a*0.4)
        bri0tomag(a) = -log10(a) * 2.5
        maginterp(a,b,f) = bri0tomag(magtobri0(a) * (1.-f) + magtobri0(b) * f)
        addlogs(a,b) = log10(10.**a + 10.**b)

        angdif(x,y)=dmod(x-y+9.4247779607694d0,6.2831853071796d0)-3.1415926535897932d0

        write(6,*)' subroutine get_glow_obj...'
        write(6,21)alt_obj_in,azi_obj_in,mag_obj,diam_deg,l_obsc,l_phase,emag,horz_dep
21      format('   alt/az/mag/diam/lobsc/lphase/emag/hrzdp = ',4f9.3,2l2,f7.4,f9.4)

!       we may want to efficiently apply refraction here to the object 
!       altitude to convert from true altitude to apparent. it may be more 
!       accurate and slower to calculate refraction down below instead.

        iverbose = 0
        call       get_airmass(alt_obj_in,htmsl,patm &   ! i 
                              ,redp_lvl,aero_scaleht &   ! i
                              ,earth_radius,iverbose &   ! i
                              ,ag,ao,aa,refr_deg)        ! o

        write(6,*)' refr_deg (informational) = ',refr_deg

!       place object in center of grid box
        if(diam_deg .eq. 0.)then
            alt_obj = nint(alt_obj_in/alt_scale) * alt_scale
            azi_obj = nint(azi_obj_in/azi_scale) * azi_scale
        else
            alt_obj = alt_obj_in
            azi_obj = azi_obj_in
        endif
 
        rpd = 3.14159265 / 180.

        i4_elapsed = ishow_timer()

!       calculate extinction
!       patm = 1.0     

        if(alt_obj .gt. 0.0)then
            call calc_extinction(90.    ,patm,airmass,zenext)
            call calc_extinction(alt_obj,patm,airmass,totexto)
            ext_mag = totexto - zenext 
        else
            ext_mag = 0.0
        endif
        ext_mag = 0. ! test of attenuation in 'get_sky_rgb'

        write(6,7)alt_obj,azi_obj,mag_obj,ext_mag
7       format('  al/az/mg/ext',2f9.3,2f9.1)

!       glow_obj = 1.0 ! cosmic background level

        sqarcsec_per_sqdeg = 3600.**2

!       if(mag_obj .lt. -5.)then ! sun or moon (0.5 degree diameter)
        if(diam_deg .ge. 0.25)then ! sun or moon (0.5 degree diameter)
            alt_dist = diam_deg / 2.0
        else                     ! star or planet (point object)
            alt_dist = alt_scale / 2.0
        endif

        radius_deg = diam_deg / 2.0

        iwrite = 0

        sum_fraclit = 0.
        sum_fraclit_sr = 0.
    
        do ialt = minalt,maxalt

!         this is apparent altitude, and we can use refraction to convert it
!         to true altitude.
          altg_app = alt_a(ialt,minazi)
          iverbose = 0
          call       get_airmass(altg_app,htmsl,patm &     ! i 
                                ,redp_lvl,aero_scaleht &   ! i
                                ,earth_radius,iverbose &   ! i
                                ,ag,ao,aa,refr_deg)        ! o
          altg1 = altg_app - refr_deg ! true altitude of grid
          
          if(refr_deg .gt. .100 .or. (refr_deg .gt. 0. .and. htmsl .gt. 100e3 .and. (l_phase .eqv. .false.) ) )then
            write(6,31)ialt,altg_app,altg1,refr_deg,ag
31          format('altg_app/altg1/refr_deg/ag',i8,3f10.4,f9.2)
          endif

          do jazi = minazi,maxazi

            size_glow_sqdg = 1.0  ! alt/az grid
            size_glow_sqdg = 0.1  ! final polar kernel size
            size_glow_sqdg = 0.3  ! empirical middle ground 

            if(altg1 .ge. -90.)then 
              altg = altg1
              azig = azi_a(ialt,jazi)
            else ! refracted below the nadir
              altg = -180. - altg1
              azig = mod(azi_a(ialt,jazi) + 180.,360.)
            endif

            if(altg .ge. -2. .or. altg_app .ge. -horz_dep)then

                alt_cos = min(altg,89.)
                azi_dist = alt_dist / cosd(alt_cos)         

!               calculate distance in object/grid radii in cyl projection
!               calculate distance in degrees in cyl projection
                if(.false.)then
                distd = sqrt(((alt_obj-altg))**2 + ((azi_obj-azig)*cosd(alt_cos))**2)
                    distr = sqrt(((alt_obj-altg)/alt_dist)**2 + ((azi_obj-azig)/azi_dist)**2)
                else
                    call great_circle(alt_obj,azi_obj,altg,azig,distd,gcbearing)
                    distr = distd / radius_deg
                endif

                if(.not. l_phase)then
                  if(htmsl .ge. 100e3 .and. (jazi .eq. (jazi/360)*360) )then
                    write(6,41)ialt,jazi,alt_obj,azi_obj,altg,azig,distd,distr
41                  format(i8,i6,2f8.3,2x,2f8.3,4x,2f8.3)
                  endif
                endif

!               initialize
                frac_lit = 1.0
                frac_lit2 = 0.0
                grid_frac_obj = 1.0
                iblock=0
                rmag_per_sqarcsec = r_missing_data

!               if(diam_deg .ge. 0.75 .and. distr .le. radius_deg)then  ! solar corona
                if(diam_deg .ge. 0.75)then  ! solar corona
                  if(distd .le. 8.0)then
                    iblock = 1
                    size_glow_sqdg = 0.2    ! sun/moon area           
                    delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
!                   distr = sqrt(((alt_obj-alt))**2 + ((azi_obj-azi)*cosd(alt_cos))**2)

!                   test a call to antialias_ellipse
                    if(.false.)then
                      if(distd .ge. 0.25)then ! corona
                        frac_lit = 0.0
                      else                    ! dark side of moon
                        frac_lit = 1.0
                      endif
                    else
                      radius_pix = 0.25 / alt_scale ! 0.25 is radius_deg of moon
                      ricen = (azi_obj-azig)/azi_scale
                      rjcen = (alt_obj-altg)/alt_scale
                      aspect_ratio = 1. / cosd(min(abs(alt_obj),89.))
                      if(abs(alt_obj-altg) .le. 0.15)then                  
!                       write(6,*)' calling antialias_ellipse with iverbose = 1'
                        call antialias_ellipse(radius_pix,ricen,rjcen,aspect_ratio,frac_lit,0,0,0)
                      else
!                       write(6,*)' calling antialias_ellipse with iverbose = 0'
                        call antialias_ellipse(radius_pix,ricen,rjcen,aspect_ratio,frac_lit,0,0,0)
                      endif
                    endif

                    if(frac_lit .lt. 1.0)then    ! possible corona location
                        distd2 = max(distd,0.25+0.5*alt_scale) ! account for pixel size 
                        srad = distd2 / 0.25
!                       https://ase.tufts.edu/cosmos/print_images.asp?id=28
                        spowerk  =  -5.7 -      (srad-1.) * 2.4
!                       spowerk  =  -5.7 -      (srad-1.) * 2.0
!                       spowerf  =  -7.5 - log10(srad   ) * 2.0
!                       spowerf  =  -8.5 -       srad     * 0.1
!                       spowerf  =  -8.3 -       srad     * 0.16
!                       http://arxiv.org/pdf/0909.1722.pdf
                        spowerf1 =  -6.7 -      (srad-1.) * 1.6 
!                       spowerf2 =  -8.3 -       srad     * 0.16
!                       https://ase.tufts.edu/cosmos/print_images.asp?id=28
                        spowerf2 =  -7.7 - log10(srad   ) * 2.3
                        spowerf =     max(spowerf1,spowerf2)                 
                        spowerc = addlogs(spowerk,spowerf)
                        smag_effc = -26.7 - 2.5 * spowerc
                    else ! corona is blocked completely by moon
                        smag_effc = 0.
                    endif
                    smag_effm = -2.5
                    smag_eff = maginterp(smag_effc,smag_effm,frac_lit)

                    rmag_per_sqarcsec = smag_eff + delta_mag 
                    if(abs(alt_obj-altg) .le. 0.15)then                  
!                   if(.true.)then                                      
                        write(6,81)ialt,jazi,altg,alt_obj,azig,azi_obj,distd,spowerk,spowerf,spowerc,smag_eff,rmag_per_sqarcsec,frac_lit
81                      format(' cor: altg,alt_obj,azig,azi_obj,distd,spk,spf,spc,smag,rmag,flit =',2i5,5f9.3,6f7.2)
                    endif
                  else
                    iblock = 2
                    rmag_per_sqarcsec = r_missing_data
                  endif
                 
                elseif(diam_deg .ge. 0.05 .and. diam_deg .lt. 0.75 .and. distd .le. 1.0)then ! regular sun or moon
                    iblock = 3
                    size_glow_sqdg = pi * radius_deg**2 ! sun/moon area           

!                   calculate fraction of grid box illuminated by object
                    if(htmsl .gt. 100000e3)then ! very simple anti-aliasing 
                      distr_thrl = 1.0 - 0. ! consider size of pixel relative to object radius
                      distr_thrh = 1.0 + 0.
                      if(distr .le. distr_thrl)then ! object radii
                        frac_lit1 = 1.0
                      elseif(distr .gt. distr_thrh)then
                        frac_lit1 = 0.0
                      else
                        frac_lit1 = 1.5 - distr
                      endif
                    elseif(l_phase)then
                      iblock = 30    
                      if(distr .le. 1.5)then ! potentially lit
                        iverbose = 2  
                        radius_pix = radius_deg / alt_scale
                        ricen = (azi_obj-azig)/azi_scale
                        rjcen = (alt_obj-altg)/alt_scale
                        aspect_ratio = 1. / cosd(min(abs(alt_obj),89.))
                        write(6,*)' calling antialias_phase with iverbose = ',iverbose,distd,distr,va
                        call antialias_phase(radius_pix,ricen,rjcen,aspect_ratio,alt_scale,azi_scale,va,rill,frac_lit1m,0,0,iverbose)
                        call antialias_phase(radius_pix,ricen,rjcen,aspect_ratio,alt_scale,azi_scale,va,1.0 ,frac_lit1e,0,0,0)
                        earthshine = .00001
                        frac_lit1 = frac_lit1m + earthshine * frac_lit1e
                        write(6,86)ialt,jazi,frac_lit1,frac_lit1m,frac_lit1e
86                      format(' frac_lit1 returned ',2i6,3f9.4)
                        size_glow_sqdg = size_glow_sqdg * rill           
                      else
                        frac_lit1 = 0.
                      endif    
                    else           ! more accurate anti-aliasing scheme
                      radius_pix = radius_deg / alt_scale
                      ricen = (azi_obj-azig)/azi_scale
                      rjcen = (alt_obj-altg)/alt_scale
                      aspect_ratio = 1. / cosd(min(abs(alt_obj),89.))

!                     lots of writes can happen with solar eclipses from dscovr
!                     if(iwrite .le. 100)then
                      if(iwrite .le. 100 .and. distd .lt. (radius_deg + alt_scale * 1.42))then
!                     if(jazi .eq. minazi)then
!                     if(.true.)then
                        iwrite = iwrite + 1
                        write(6,87)jazi,iwrite,distd,altg,azig,radius_deg,alt_scale,radius_pix
87                      format(' call antialias_ellipse: dist/alt/az/raddeg/alt_scl/radpix = ',i5,i5,3f8.3,3f10.5)
                        iverb = 2
                      else
                        iverb = 0
                      endif

!                     obtain fraction of pixel that is lit (frac_lit1)
                      call antialias_ellipse(radius_pix,ricen,rjcen,aspect_ratio,frac_lit1,0,0,iverb)
                      area_pix = alt_scale * azi_scale
                    endif

!                   consider second "object" for obscuration or phase
                    if(l_obsc)then 
                      rill = 0.5
                      radius_pix2 = radius_pix
                      ricen = (azi_obj2-azig)/azi_scale
                      rjcen = (alt_obj2-altg)/alt_scale
                      aspect_ratio = 1. / cosd(min(abs(alt_obj),89.))
                      call antialias_ellipse(radius_pix2,ricen,rjcen &
                                            ,aspect_ratio,frac_lit2,0,0,0)

!                     thresholded and random strategies are being tested. 
!                     it might work better for the 'antialias_ellipse' 
!                     routine to consider both objects simultaneously.
                      if(.false.)then
!                       frac_lit = min(frac_lit1,(1.-frac_lit2))
                        frac_lit = frac_lit1 * (1.-frac_lit2) 
                        if(frac_lit2 .gt. frac_lit1)then
                          frac_lit = 0.
                        else
                          frac_lit = frac_lit1 - frac_lit2
                        endif
                      else ! consider approx overlap as f(emag)
                        aov1 = max((frac_lit1 + frac_lit2 - 1.0),0.)
                        aov2 = frac_lit1 * frac_lit2
                        aov3 = frac_lit2
                        fr1 = max(1.-abs(emag-0.0)/0.5,0.)
                        fr2 = max(1.-abs(emag-0.5)/0.5,0.)
                        fr3 = max(1.-abs(emag-1.0)/0.5,0.)
                        approx_overlap = aov1*fr1 + aov2*fr2 + aov3*fr3
                        frac_lit = max(frac_lit1-approx_overlap,0.)
                        write(6,88)aov1,aov2,aov3,fr1,fr2,fr3,approx_overlap,frac_lit
88                      format(' aov',3f7.2,'  fr',3f7.2,'  overlap/frac_lit',2f7.4)
                      endif
                    else ! l_obsc = f
                      frac_lit = frac_lit1 
                    endif

                    sum_fraclit    = sum_fraclit + frac_lit

                    pix_area_sqdg = area_pix * cosd(altg)
                    pix_area_sr = pix_area_sqdg * (pi/180.)**2
                    sum_fraclit_sr = sum_fraclit_sr + frac_lit * pix_area_sr

!                   write(6,89)ialt,jazi,altg_app,altg,alt_obj,azig,azi_obj,distd,frac_lit1,frac_lit2,frac_lit
89                  format(' altga-t/alt_obj/azig/azi_obj/distd/frac_lit =',i6,i5,6f9.3,2f7.2,f7.4)

!                   determine surface brightness of object averaged over the pixel
                    if(frac_lit .gt. 0.)then
!                       average pixel surface brightness accounting for fraction of pixel that is lit
                        delta_mag = log10((size_glow_sqdg*sqarcsec_per_sqdeg)/frac_lit)*2.5
                        rmag_per_sqarcsec = mag_obj + delta_mag                  
                    else
                        rmag_per_sqarcsec = r_missing_data
                    endif

                elseif(distd .le. radius_deg)then ! star
                    iblock = 4
                    delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
                    rmag_per_sqarcsec = mag_obj + delta_mag                  
                else
                    iblock = 5
                    rmag_per_sqarcsec = r_missing_data
                endif

                if(rmag_per_sqarcsec .ne. r_missing_data)then
!                 convert to nanolamberts. this has an empirical correction for now based on
!                 sunlight spread over the spherical solid angle being 3e9nl. max nl should be 5.538e14
                  glow_nl = v_to_b(rmag_per_sqarcsec) / 1.22062
                  glow_obj(ialt,jazi) = addlogs(glow_obj(ialt,jazi),log10(glow_nl))

!                 if(.true.)then                          
!                 if(abs(alt_obj-altg) .le. 0.15 .and. iwrite .le. 100)then                  
                  if(iblock.eq.30 .or. iblock .eq. 3)then                  
                      write(6,90)ialt,jazi,altg_app,altg,alt_obj,azig,azi_obj,distd,frac_lit1,frac_lit2,frac_lit
90                    format(' altga-t/alt_obj/azig/azi_obj/distd/area/frclit =',2i5,6f9.3,2f7.2,f7.4)

                      if(glow_nl .lt. 1e5)then
                          write(6,91)ialt,jazi,diam_deg,rmag_per_sqarcsec,delta_mag,glow_nl,glow_obj(ialt,jazi)
91                        format(' rmag_per_sqarcsec/dmag/glow_nl/glow_obj =     ',2i5,f5.2,2f10.3,e12.4,f11.2)
                      else
                          write(6,92)ialt,jazi,diam_deg,rmag_per_sqarcsec,delta_mag,glow_nl,glow_obj(ialt,jazi)
92                        format(' rmag_per_sqarcsec/dmag/glow_nl/glow_obj =     ',2i5,f5.2,2f10.3,e12.4,f11.2,' ***glow***')
                      endif
                  endif
                else ! missing rmag_per_sqarcsec
                  if(distd .le. 0.15)then
                    write(6,*)' warning: no glow assigned when close to object ',iblock
                  endif
                endif ! within star kernel
            endif ! alt > -2.
          enddo ! jazi
        enddo ! ialt

!       calculate solar surface brightness for reference
        sun_radius_deg = 0.5334 / 2.
        sun_area_sqdg = pi * sun_radius_deg**2
        sun_area_sr = sun_area_sqdg * (pi/180.)**2
        sun_area_sphere = sun_area_sr / (4. * pi)
        sun_nl = day_int0 / sun_area_sphere
        sun_radiance_calc = ghi_zen_toa / sun_area_sr
        sun_radiance      = 2.009e7 ! w/m^2/sr (via wikipedia)

!       solar radiance (wikipedia) is 2.009e7 w/m^2/sr

        write(6,93)diam_deg,mag_obj,delta_mag,sun_area_sr,sun_nl,log10(sun_nl)
93      format(' diam/mag/dmag/sun_sr/sun_nl/glow_sun_nl = ',3f9.3,2e12.4,f9.3,' ***glow***')

        write(6,*)' solar radiance is ',sun_radiance_calc,sun_radiance

        write(6,*)' sum_fraclit (pixels / sr) = ',sum_fraclit,sum_fraclit_sr

        write(6,*)' returning from get_glow_obj...'
        write(6,*)

        return
        end 

        subroutine get_twi_glow_ave(glow,alt_a,azi_a,ni,nj &
                                   ,sol_alt,sol_az,twi_glow_ave)

        angdif(x,y)=mod(x-y+540.,360.)-180.

        real glow(ni,nj)
        real alt_a(ni,nj)
        real azi_a(ni,nj)

        cnt = 0.
        sum = 0.

!       compare with -22. - sol_alt (as an integrated magnitude)
!       variable would be 'twi_mag'
        do i = 1,ni
        do j = 1,nj
            diffaz = angdif(azi_a(i,j),sol_az)
            if( abs(diffaz) .le. 90. .and. &
               (alt_a(i,j)  .le. 20. .and. alt_a(i,j) .ge. 0.) )then            
                cnt = cnt + 1.0
                sum = sum + 10.**glow(i,j)
            endif
        enddo ! j
        enddo ! i

        if(cnt .gt. 0.)then
            twi_glow_ave = log10(sum/cnt)
        else ! tested at -6.2 degrees
            write(6,*)' warning in get_twi_glow_ave - using empirical function'
            twi_glow_ave = 8.8 + (sol_alt*0.3)
            write(6,*)' sol_alt/sol_az/twi_glow_ave = ' &
                       ,sol_alt,sol_az,twi_glow_ave
        endif

        return
        end

        subroutine get_sky_rad_ave(rad,alt_a,azi_a,ni,nj &
                                  ,sky_rad_ave)

        include 'trigd.inc'

        real rad(ni,nj)
        real alt_a(ni,nj)
        real azi_a(ni,nj)

        write(6,*)' rad range is ',minval(rad),maxval(rad)

        alt_top = alt_a(ni,1)

!       average over window
        cnt = 0.
        sum = 0.
        do i = 1,ni
          cosi = cosd(alt_a(i,1))
          if(alt_a(i,1) .ge. 0.)then ! above horizontal
            do j = 1,nj
                cnt = cnt + (1.0      * cosi)
                sum = sum + (rad(i,j) * cosi)     
            enddo ! j
          endif
        enddo ! i

        if(cnt .gt. 0.)then
            sky_rad_ave_wdw = sum/cnt
        else
            write(6,*)' error in get_sky_rad_ave'
            sky_rad_ave_wdw = 10. 
        endif

!       average on top row
        cnt = 0.
        sum = 0.
        do i = ni,ni
          cosi = cosd(alt_a(i,1))
          if(alt_a(i,1) .ge. 0.)then ! above horizontal
            do j = 1,nj
                cnt = cnt + (1.0      * cosi)
                sum = sum + (rad(i,j) * cosi)     
            enddo ! j
          endif
        enddo ! i

        if(cnt .gt. 0.)then
            sky_rad_ave_top = sum/cnt
        else
            write(6,*)' warning in get_sky_rad_ave'
            write(6,*)' alt_top = ',alt_top
            sky_rad_ave_top = rad(ni,1)
        endif

!       if needed, the average on the top row is extrapolated for the
!       portion of the sky above the camera field of view (when alt_top 
!       is less than 90 degrees).
        if(alt_top .gt. 0.)then
          write(6,*)' get_sky_rad_ave: top above horizon ',alt_top
          area_frac = sind(alt_top)
          write(6,*)' area_frac/wdw/top = ',area_frac,sky_rad_ave_wdw,sky_rad_ave_top
          sky_rad_ave = sky_rad_ave_wdw * area_frac + sky_rad_ave_top * (1.-area_frac)
        else
          write(6,*)' get_sky_rad_ave: top below horizon ',alt_top
          sky_rad_ave = sky_rad_ave_wdw
        endif

        return
        end

        subroutine rad_2_irrad_down(rad,alt_a,azi_a,elong_a,ni,nj,alt_scale,azi_scale,irrad,beam,direct,diffuse)

!       subroutine name relates to primary output of ghi, with other outputs such as dhi and dni.          

!       input units can then be solar relative (to solar irradiance spread out over the full
!       sphere - 4 pi steradians). the (wideband or spectral) irradiance integrates to 1.
!       if the input (wideband or spectral) radiance has a value of 2. over the upper
!       hemisphere - 2 pi steradians).

        include 'trigd.inc'
        use mem_allsky, only: day_int0

        parameter (pi=3.14159265)
        parameter (rpd = pi/180.)

        real rad(ni,nj),irrad
        real alt_a(ni,nj)
        real azi_a(ni,nj)
        real elong_a(ni,nj)

        real*8 cnt,sum_ghi,sum_dni,sum_diffuse

        sincosint(x) = -0.5 * cosd(x)**2

        alt_top = alt_a(ni,1)

!       average over window
        cnt_ghi = 0.
        cnt_dni = 0.
        cnt_diffuse = 0.

        sum_ghi = 0.
        sum_dni = 0.
        sum_bhi = 0.
        sum_diffuse = 0.

        do i = 1,ni
          cosi = cosd(alt_a(i,1))
          sini = sind(alt_a(i,1))
          solidangle_pix = (alt_scale*rpd) * (azi_scale*rpd) * cosi 

          if(alt_a(i,1) .ge. 0.)then ! above horizontal
            do j = 1,nj
                sum_ghi         = sum_ghi     + (rad(i,j) * solidangle_pix * sini)     

                if(elong_a(i,j) .le. 2.5)then ! cone of acceptance for solar instrumentation
                    sum_dni     = sum_dni     + (rad(i,j) * solidangle_pix)     
                    sum_bhi     = sum_bhi     + (rad(i,j) * solidangle_pix * sini)     
                else
                    sum_diffuse = sum_diffuse + (rad(i,j) * solidangle_pix * sini)     
                endif

                cnt_ghi = cnt_ghi + (solidangle_pix * sini)
                cnt_dni = cnt_dni + (solidangle_pix       )
                cnt_diffuse = cnt_ghi

            enddo ! j
          endif ! alt > 0
        enddo ! i

        if(cnt_ghi .gt. 0.)then
            sky_rad_sum_wdw = sum_ghi
        else
            write(6,*)' error in rad_2_irrad_down'
            sky_rad_sum_wdw = 10. 
        endif

!       average on top row
        if(alt_top  .lt. 90.)then
          cnt = 0.
          sum = 0.
          do i = ni,ni
            cosi = cosd(alt_a(i,1))
            sini = sind(alt_a(i,1))
            if(alt_a(i,1) .ge. 0.)then ! above horizontal
              do j = 1,nj
                sum = sum + (rad(i,j) * cosi * sini)     
                cnt = cnt + (1.0      * cosi)
              enddo ! j
            endif
          enddo ! i

          if(cnt .gt. 0.)then
            sky_rad_ave_top = sum/cnt
          else
            write(6,*)' warning in get_sky_rad_ave'
            write(6,*)' alt_top = ',alt_top
            sky_rad_ave_top = rad(ni,1)
          endif
        endif

!       if needed, the average on the top row is extrapolated for the
!       portion of the sky above the camera field of view (when alt_top 
!       is less than 90 degrees).

!       'cnt_ghi' will integrate (via the sin effect) to pi, half of the
!       2 pi steradians subtended by the upper hemisphere.

        if(alt_top .eq. 90.)then
          irrad   = sky_rad_sum_wdw / (4. * pi) ! (2. * cnt_ghi)
          beam    = sum_bhi         / (4. * pi) ! (2. * cnt_ghi)
          direct  = sum_dni         / (4. * pi) ! (2. * cnt_dni)
          diffuse = sum_diffuse     / (4. * pi) ! (2. * cnt_ghi)
          write(6,*)' rad_2_irrad_down: rad | solrel | reflectance sample = ',rad(ni,1),rad(ni,1)*day_int0,rad(ni,1)/2.
          rad_max = maxval(rad)
          write(6,*)' rad_2_irrad_down: rad | solrel | reflectance max    = ',rad_max,rad_max*day_int0,rad_max/2.
          write(6,1)alt_top,sky_rad_sum_wdw,irrad,beam,direct,diffuse,cnt_ghi,cnt_dni
1         format('  rad_2_irrad_down: top at zenith ',f9.1,5f9.5,' cnt ',2f9.4)
        elseif(alt_top .gt. 0.)then
          write(6,*)' rad_2_irrad_down: top above horizon ',alt_top
          area_frac = (sincosint(alt_top) - sincosint(0.)) / (sincosint(90.) - sincosint(0.))
          sky_rad_sum_top = sky_rad_ave_top * (1.-area_frac) * 2. * pi
          write(6,*)' area_frac/wdw/top = ',area_frac,sky_rad_sum_wdw,sky_rad_sum_top
          irrad = sky_rad_sum_wdw + sky_rad_sum_top
        else
          write(6,*)' get_irrad: top below horizon ',alt_top
          irrad = sky_rad_ave_wdw
        endif

!       sp_irrad = sp_irrad * 2. * pi ! solid angle above horizontal

        return
        end

        subroutine apply_rel_extinction(rmaglim,alt,od_zen)

        include 'trigd.inc'

        z = 90. - alt
        airmass = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))
        od_slant = od_zen * (airmass - 1.0)
        delta_mag = od_slant * 1.086
        rmaglim = rmaglim - delta_mag

!       write(6,1)airmass,od_zen,od_slant,delta_mag,rmaglim 
1       format(6x,'airmass,od_zen,od_slant,delta_mag,rmaglim ',5f9.2)

        return
        end

        subroutine nl_to_sprad(nl,nc,wa,sprad)

        real nl(nc)       ! i (spectral radiance in scaled solar relative units)
        real wa(nc)       ! i (wavelength in microns)
        real fa(nc)       ! l (solar spectral irradiance w/(m**2 nm))
        real sprad(nc)    ! o (spectral radiance in w/(m**2 nm sr)) units)

        parameter (pi=3.14159265)

!       constant 3e9      ! l (scaling constant, approx solar illuminance 
                          !    spread over a spherical solid angle [nl])

        iverbose = 0

!       obtain solar spectral radiance
        call get_fluxsun(wa,nc,iverbose,fa)

        do ic = 1,nc
            sprad(ic) = (nl(ic)/3e9) * fa(ic) / (4. * pi)
        enddo ! ic

        continue

        return
        end

        subroutine get_auroraglow(i4time,alt_a,azi_a,minalt,maxalt,minazi,maxazi,rlat,rlon,alt_scale,azi_scale,glow_aurora)

        real glow_aurora(3,minalt:maxalt,minazi:maxazi) ! log nl
        real alt_a(minalt:maxalt,minazi:maxazi)
        real azi_a(minalt:maxalt,minazi:maxazi)

        return
        end
