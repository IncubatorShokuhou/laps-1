
        subroutine nl_to_rgb(rad,glwref,contrast,cntref,wa,iverbose,xyz,rc,gc,bc)

        parameter (nc = 3)

        real wa(nc)                    ! wavelength (um)
!       data wa    /.615,.546,.450/

!       luminance of sunlight from 1au distance uniformly scattered 
!       throughout a 360 degree sphere. this will be defined as the
!       spectral radiance matching that of the solar spectrum.
        parameter (day_int0 = 3e9)

        include 'wac.inc'

        real luminance ! candela / m**2
        real luma_of_counts
        real xyz(3)
        integer init /0/
        save init,fasun

        parameter (day_int = day_int0)

!       exponential clipping of gamma function
!       expgamma(x) = 1. - exp(-x**(1./gamma))
!       expgaminv(x) = (-(log(1.-x)))**gamma ! inverse of expgamma
!       radhi(glwref,cntref) = (10.**glwref) / (expgaminv(cntref/255.))
!       rad_to_counts2(rad1) = expgamma(rad1/radhi(glwref,cntref)) * 255.

!       straight gamma function
        fgamma(x) = x**(1./gamma)
        fgaminv(x) = x**gamma ! inverse of expgamma
        radhi(glwref,cntref) = (10.**glwref) / (fgaminv(cntref/255.))
        rad_to_counts2(rad1) = fgamma(rad1/radhi(glwref,cntref)) * 255.

!       convert nl to rintensity

!       radiance units are watts/(m**2 sr nm)

!       day_int0 is luminance of sunlight from 1au distance uniformly 
!       scattered throughout a 360 degree sphere. this will be defined as
!       the spectral radiance matching that of the solar spectrum.
!       note solar constant is 1361.5 watts/m**2
        real rad(nc) 
        real rel_solar(nc), cdm2(nc)
        real fasun(nct), rel_solar_extrap(nct), farad(nct)
        real nl_int

        if(init .eq. 0)then
            write(6,*)' initializing a call to get_fluxsun ',init
            call get_fluxsun(wa_tri,nct,1,fasun)
            init = 1
        endif

        if(maxval(rad(:)) .eq. 0.)then
            rc = 0.; gc = 0.; bc = 0.
            return
        endif

!       nl_int = rgb_to_y(rad(1),rad(2),rad(3))

        rel_solar(:) = rad(:) / day_int
        rel_solar_luminance = rgb_to_y(rel_solar(1),rel_solar(2) &
                                      ,rel_solar(3))
        if(iverbose .eq. 1)then
          write(6,*)
          write(6,*)'rel_solar = ',rel_solar
          write(6,*)'rel_solar_luminance = ',rel_solar_luminance
        endif

!       1 lambert = 3183.0988 candela/m**2        (luminance)
!       1 lux = 1 lumen/m**2                      (illuminance)
!       1 lux = 1 lumen/m**2                      (luminous exitance)
!       1 candela = 1/683.002 watts/sr (at 550nm) (luminous intensity)
!       1 watt = 683 lumens (555nm source)        (radiant flux)
!       1 watt =  93 lumens (solar spectrum)      (radiant flux)
!       1 watt/m**2                               (irradiance)
!       1 lumen = 1 cd*sr                         (luminous flux)
!       1 candela = 1 lumen per steradian         (luminous intensity)
        cdm2(:) = (rad(:)/1e9) * 3183.0988 ! nl to candela/m**2

!       extrapolate color            
        call extrap_color(rel_solar,wa,rel_solar_extrap) ! dimensionless
        farad(:) = fasun(:) * rel_solar_extrap(:)        ! spect irradiance
        if(iverbose .eq. 1)then
           write(6,*)' ict    wa    rel_solar    fa     farad = '
           do ict = 1,nct
              write(6,11)ict,wa_tri(ict),rel_solar_extrap(ict),fasun(ict),farad(ict)
11            format(i4,f9.4,f9.4,f9.4,f9.4)
           enddo ! ict
        endif

!       convert sprad to xyz
        call get_tricolor(farad,iverbose,wa,xx,yy,zz,x,y,z,luminance)
        if(iverbose .eq. 1)write(6,*)'xyzfarad = ',x,y,z
        if(iverbose .eq. 1)write(6,*)'luminance farad = ',luminance
        xyz(1) = x; xyz(2) = y; xyz(3) = z

!       convert xyz to (linear) rgb
        ct = 5780.
        if(iverbose .eq. 1)write(6,*)'ct = ',ct
        call xyztosrgb(x,y,z,ct,r,g,b)             

!       convert rintensity to rgb (gamma correction)
        call linearrgb_to_counts(r,g,b,rc,gc,bc)

        if(iverbose .eq. 1)write(6,*)'glwref = ',glwref

!       if(iverbose .eq. 1)write(6,*)'desired rad = ',desired_rad

        luma_of_counts = rgb2luma(rc,gc,bc)
        if(iverbose .eq. 1)write(6,*)'luma_of_counts = ',luma_of_counts

!       solar_counts = 240.
        solar_counts = rad_to_counts2(day_int)
        if(iverbose .eq. 1)then 
            write(6,*)'day_int = ',day_int
            write(6,*)'glwref = ',glwref
            write(6,*)'cntref = ',cntref
            write(6,*)'solar_counts = ',solar_counts
        endif

!       desired_luma = 240. * rel_solar_luminance**(1./gamma)
!       desired_luma = solar_counts + log10(rel_solar_luminance) * contrast
        desired_luma = rad_to_counts2(day_int*rel_solar_luminance)
        if(iverbose .eq. 1)write(6,*)'desired luma = ',desired_luma
        scale_luma = desired_luma / luma_of_counts
        if(iverbose .eq. 1)write(6,*)'scale_luma = ',scale_luma
        rc = rc * scale_luma       
        gc = gc * scale_luma       
        bc = bc * scale_luma       
        
        return
        end

        subroutine extrap_color(a_nc,wa,a_nct)

        parameter (nc = 3)

        real wa(nc)                    ! wavelength (um)

!       luminance of sunlight from 1au distance uniformly scattered 
!       throughout a 360 degree sphere. this will be defined as the
!       spectral radiance matching that of the solar spectrum.
        parameter (day_int0 = 3e9)

        include 'wac.inc'

        real a_nc(nc), a_nct(nct)

!       assuming nc = 3, find polynomial fit
        x1 = log(wa(1))  ; x2 = log(wa(2))  ; x3 = log(wa(3))
        y1 = log(a_nc(1)); y2 = log(a_nc(2)); y3 = log(a_nc(3))

        a = ((y2-y1)*(x1-x3) + (y3-y1)*(x2-x1)) / &
            ((x1-x3)*(x2**2-x1**2) + (x2-x1)*(x3**2-x1**2))
        b = ((y2 - y1) - a*(x2**2 - x1**2)) / (x2-x1)
        c = y1 - a*x1**2 - b*x1

        a_nct(:) = exp(a*log(wa_tri(:))**2 + b*log(wa_tri(:)) + c)
 
        return
        end

        subroutine rad_to_xyz(rad_in,x,y,z)

        real rad_in(3),rad(3)

        rad(1) = rad_in(1) + 0.2 * rad_in(3)
        rad(2) = rad_in(2)
        rad(3) = rad_in(3)

        x = rad(1) / sum(rad)
        y = rad(2) / sum(rad)
        z = rad(3) / sum(rad)

        return
        end

        subroutine linearrgb_to_counts(r,g,b,rc,gc,bc)

        parameter (gamma = 2.2)

        rr = r ** (1./gamma)
        gg = g ** (1./gamma)
        bb = b ** (1./gamma)

        rc = rr*255. 
        gc = gg*255. 
        bc = bb*255. 
      
        return
        end

        subroutine xyztorgb(x,y,z,r,g,b)

        xx = x; yy = y; zz = z

        r=  (0.41847   *xx)-(0.15866  *yy)-(0.082835*zz)
        g= -(0.091169  *xx)+(0.25243  *yy)+(0.015708*zz)
        b=  (0.00092090*xx)-(0.0025498*yy)+(0.17860 *zz)

        return
        end

        subroutine rgbtoxyz(r,g,b,x,y,z) 

        x = (.49    * r + .31    * g + .20    * b) / 0.17697
        y = (.17697 * r + .81240 * g + .01063 * b) / 0.17697
        z = (             .01    * g + .99    * b) / 0.17697

        return
        end

        subroutine xyztosrgb(x,y,z,ct,r,g,b)

        colfrac = (ct - 5000.) / (6500. - 5000.)
 
        r65 =  3.2406*x -1.5372*y -0.4986*z
        g65 = -0.9689*x +1.8758*y +0.0415*z
        b65 =  0.0557*x -0.2040*y +1.0570*z

        r50 =  3.1339*x -1.6169*y -0.4906*z
        g50 = -0.9785*x +1.9160*y +0.0333*z
        b50 =  0.0720*x -0.2290*y +1.4057*z

        r = r65 * colfrac + r50 * (1.-colfrac)
        g = g65 * colfrac + g50 * (1.-colfrac)
        b = b65 * colfrac + b50 * (1.-colfrac)

!       clip at 0. when outside the gamut
        r = max(r,0.)
        g = max(g,0.)
        b = max(b,0.)

        return
        end

        subroutine srgbtoxyz(r,g,b,ct,x,y,z)
 
        colfrac = (ct - 5000.) / (6500. - 5000.)

        x65 =  0.4124*r + 0.3576*g + 0.1805*b
        y65 =  0.2126*r + 0.7152*g + 0.0722*b
        z65 =  0.0193*r + 0.1192*g + 0.9505*b
 
        x50 =  0.4361*r + 0.3851*g + 0.1430*b
        y50 =  0.2225*r + 0.7169*g + 0.0606*b
        z50 =  0.0139*r + 0.0971*g + 0.7141*b

        x = x65 * colfrac + x50 * (1.-colfrac)
        y = y65 * colfrac + y50 * (1.-colfrac)
        z = z65 * colfrac + z50 * (1.-colfrac)
 
        return
        end

        subroutine get_fluxsun(wa,nc,iverbose,fa)

!       returns solar spectral irradiance for selected wavelengths

        real fa(nc) ! watts/(m**2 nm)
        real wa(nc) ! wavelengths in microns

        if(iverbose .eq. 1)write(6,*)'  subroutine get_fluxsun',nc,wa

!       note that 5900k is the correlated color temperature of the sun
!       giving the best fit to a blackbody curve shape. the total energy
!       emitted though is related to the effective temperature of 5780k.
!       we can thus consider the emissivity to be related to this difference.
        bbts = 5780.
        ccts = 5900. 
        emissivity = (bbts / ccts)**4.
        ds = 1.0 * 1.496e13 ! au to cm
        rs = 1.0 * 6.96e10  ! solar radii to cm
        xc = 0.; yc = 0.; zc = 0.

        do ic = 1,nc   
          w = wa(ic) * 1e-4 ! microns to cm  
          w_ang = wa(ic) * 10000.
          bb = (.0000374/(w**5.)) / (exp(1.43/(w*ccts))-1.)
          if(iverbose .eq. 1)write(6,*)'ic/w/bb',ic,w,bb
          fa(ic)=((rs/ds)**2)*bb*emissivity*1e-8 ! erg/cm2/s/a
          fa(ic)=fa(ic) * .01                    ! convert to w/m2/nm
        enddo ! ic

        if(iverbose .eq. 1)write(6,*)'fa sun (w/m2/nm) is ',fa

        return
        end

        subroutine get_tricolor(fa,iverbose,wa,xc,yc,zc,x,y,z,luminance)

        parameter (nc = 3)

        real wa(nc)                    ! wavelength (um)

!       luminance of sunlight from 1au distance uniformly scattered 
!       throughout a 360 degree sphere. this will be defined as the
!       spectral radiance matching that of the solar spectrum.
        parameter (day_int0 = 3e9)

        include 'wac.inc'

!       integrate spectral irradiance with the color matching functions

        real fa(nct)   ! spectral irradiance array vs wavelength (input)
        real luminance ! candela / m**2                    (output = yc)
        character*255 static_dir

!       3500 to 8000 angstroms in 10 steps (color matching functions)
        real x1(nct)!/0, .014, .336, .005, .433, 1.062, .283, .011, 0, 0/
        real y1(nct)!/0, .0004, .038, .323, .995, .631, .107, .004, .0001, 0/
        real z1(nct)!/0, .068, 1.773, .272, .009, 0,0,0,0,0/
        integer init/0/
        save x1,y1,z1,init

        if(init .eq. 0)then
          call get_directory('static',static_dir,len_dir)
          open(11,file=trim(static_dir)//'/cie2.txt',status='old',form='formatted')
!         open(11,file='cie2.txt',status='old')
          do ict = 1,nct
            read(11,*)inm,x1(ict),y1(ict),z1(ict)
          enddo ! ict
          init = 1
        endif

        xc = 0.; yc = 0.; zc = 0.

!       integrate with trapezoidal rule
        do ic = 2,nct   
          icm = ic-1
          w = wa_tri(ic) * 1e-4 ! microns to cm  
          w_ang = wa_tri(ic) * 10000.

          xs = 0.5 * (x1(icm)*fa(icm) + x1(ic)*fa(ic))
          ys = 0.5 * (y1(icm)*fa(icm) + y1(ic)*fa(ic))
          zs = 0.5 * (z1(icm)*fa(icm) + z1(ic)*fa(ic))
  
          xc = xc + xs
          yc = yc + ys
          zc = zc + zs

          if(iverbose .eq. 1)then
            write(6,1)ic,w_ang,fa(ic),xs,ys,zs,xc,yc,zc
1           format(' ic,wa,fa',i4,f8.0,f10.5,' xyzs',3f10.5,' xyzc',3f10.5)
          endif
 
        enddo ! ic

        x = xc / (xc+yc+zc)
        y = yc / (xc+yc+zc)
        z = zc / (xc+yc+zc)

        luminance = yc

        if(iverbose .eq. 1)write(6,11)x,y,z,luminance
11      format(' xyzl tricolor is ',3f10.6,f12.0)

        return
        end

      subroutine hsl_to_rgb(hue,sat,rintens,red,grn,blu)

      real luma

!     hue is 0:r, 1:b, 2:g, 3:r

      red1 = max(1.0 - abs(hue - 0.0),0.0)
      red2 = max(1.0 - abs(hue - 3.0),0.0)
      red = max(red1,red2)
      grn = max(1.0 - abs(hue  - 2.0),0.0)
      blu = max(1.0 - abs(hue  - 1.0),0.0)

!     normalize to the max intensity
      colmax = max(red,grn,blu)
      if(colmax .gt. 0.)then
          red = red/colmax
          grn = grn/colmax
          blu = blu/colmax
      endif

      red = (red*sat) + 1.0*(1.0-sat)
      grn = (grn*sat) + 1.0*(1.0-sat)
      blu = (blu*sat) + 1.0*(1.0-sat)

      red = red * rintens
      grn = grn * rintens
      blu = blu * rintens

!     scale the rgb intensity according to the luma calculation

      luma = .30 * red + .59 * grn + .11 * blu

      if(luma .gt. 0.)then
          ratio_corr = rintens / luma
      else
          ratio_corr = 1.0
      endif

      red = red * ratio_corr
      grn = grn * ratio_corr
      blu = blu * ratio_corr

      return
      end
