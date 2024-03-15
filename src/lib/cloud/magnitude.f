        subroutine magnitude(iplan,       ! planet number
     1                       isat,        ! satellite number (0 for the planet)
     1                       ex,ey,ez,    ! sun - observer radius vector (au)
     1                       px,py,pz,    ! sun - satellite radius vector (au)
     1                       mag,diam_sec)! magnitude, diameter in arcsec

        implicit real*8(a-z)

        include '../../include/astparms.for'

        integer iplan,isat,len_ast_root,len_mag_file

        dimension v10(0:13,0:70)
        dimension r_km(0:13,0:70)
        character*10 nm(0:13,0:70)

        character*150 c_ast_root,c_mag_file,c_line

!     saturns rings (zx,zy,zz are already present in utilparms.for)
      data xx,xy,xz/-.6216568507d0,0.7832897037d0,0.d0/
      data yx,yy,yz/-.7779524722d0,-.6174209639d0,0.1165388524d0/
!     data zx,zy,zz/-.0912836831d0,-.0724471759d0,-.9931861335d0/

        other_corr = 0.

        dx=px-ex
        dy=py-ey
        dz=pz-ez
        rlc=sqrt(px*px+py*py+pz*pz)
        delta=sqrt(dx*dx+dy*dy+dz*dz)

        if(iplan .eq. 0)then
            mag = v10(iplan,isat) + 5.*dlog10(delta)
            goto400
        endif

        if(delta .eq. 0.)then
            write(6,*)' delta = 0.'
            return
        endif

        call phase(ex,ey,ez,px,py,pz,phase_angle_rad,r_ill)

        phase_angle_deg = phase_angle_rad * dpr

        call get_physical_data(v10,r_km,nm)       

        call phase_func(iplan,isat,phase_angle_deg
     1                 ,phase_corr)       

!       saturn
        if(iplan .eq. 6)then
            if(isat .eq. 0)then
                sinb=-(zx*dx+zy*dy+zz*dz)/delta
                other_corr = ! +2.5210143d0*du
     1                         -abs(2.6d0*sinb)+1.25d0*sinb*sinb
            else if(isat .eq. 8)then ! iapetus
                theta = 0.! angle from western elongation
                h = 0.571 - .429 * sin(theta)
                other_corr = -5. * dlog10(h/.571)
            endif
        endif

        mag = v10(iplan,isat) + 5.*dlog10(rlc*delta)
     1                      + phase_corr + other_corr

        radii = (r_km(iplan,isat)/km_per_au) / delta ! radius / distance

400     diam_sec = asin(radii) * 2d0 * dpr * 3600d0

        call correct_near(radii,phase_angle_deg,mag,mag_new)
        mag = mag_new

!       write(6,*)' mag,mag_new ',mag,mag_new

        if(iplan .gt. 0 .and. .false.)then
            if(isat .eq. 0)then
                write(6,*)' sat   name       v10 phase_angle ',
     1                  'phase_corr other    mag      diam   geo alb'
                write(13,*)' sat   name       v10 phase_angle ',
     1                  'phase_corr other    mag      diam   geo alb'
            endif

            planet_brightness = 10.**(-v10(iplan,isat)*0.4)
            solar_brightness  = 10.**(-v10(0,0)*0.4)
            radius10_rad = r_km(iplan,isat) / km_per_au

!d               geoalb = planet_brightness
!d      1              / r_km(iplan,isat)**2. * 1.00e6

                geoalb = planet_brightness
     1                 / (solar_brightness * radius10_rad **2)

                write( 6,111)isat,nm(iplan,isat),v10(iplan,isat)
     1                  ,phase_angle_deg,phase_corr,other_corr
     1                  ,mag,diam_sec,geoalb
                write(13,111)isat,nm(iplan,isat),v10(iplan,isat)
     1                  ,phase_angle_deg,phase_corr,other_corr
     1                  ,mag,diam_sec,geoalb
        endif
111     format(i4,3x,a10,f6.2,5f9.2,f9.4)

        return
        end


        subroutine get_ast_root(c_ast_root_out,len_ast_root)

        character*(150) c_ast_root
        character*(*)   c_ast_root_out

        data c_ast_root /'..'/

        c_ast_root_out = c_ast_root
        len_ast_root = 2

        return
        end

        subroutine phase_ramp(p1,p2,a1,a2,phase_angle)

        if(phase_angle .le. p1)then
            continue
        elseif(phase_angle .ge. p2)then
            p1 = p2
        else
            frac = (phase_angle - a1) / (a2 - a1)
            p1 = p1 * (1.0 - frac) + p2 * frac
        endif

        return
        end

        subroutine get_physical_data(v10,r_km,nm)       

        implicit real*8(a-z)

        integer i,j,len_ast_root,len_mag_file,len_dir

        dimension v10(0:13,0:70)  ! opposition magnitude
        dimension r_km(0:13,0:70) ! radius of object
        character*10 nm(0:13,0:70)! name

        character*150 c_ast_root,c_mag_file,c_line

        data init/0/

        character*150 static_dir,filename
        character*120 cline

        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/magnitude.dat'          

        if(init .eq. 1)go to 102

!       nm = '          '

!       read i, j, v10(i,j), r_km(i,j), nm(i,j)
        c_mag_file = filename
        call s_len(filename,len_mag_file)

        open(15,file=c_mag_file(1:len_mag_file),status='old')
 50     read(15,51,err=101,end=101)c_line
 51     format(a)

        if(c_line(1:1) .eq. '!')then
            goto50
        endif

        read(c_line,61)i, j
        read(c_line,62)varg, opp_corr, r_km(i,j), nm(i,j)

        v10(i,j) = varg - opp_corr

 61     format(i10,i4)
 62     format(18x,f6.0,4x,f5.0,4x,f10.0,4x,a)

        go to 50

 101    init = 1

 102    continue

        return
        end


        subroutine phase_func(iplan,isat,phase_angle_deg
     1                       ,phase_corr)       

        implicit real*8(a-z)

        include '../../include/astparms.for'

        integer iplan,isat,len_ast_root,len_mag_file

!       dimension v10(0:13,0:70)
!       dimension r_km(0:13,0:70)
!       character*10 nm(0:13,0:70)

        para(x,x1) = (x-x1)**2

        phase_corr = 0.

!       default phase corrections
        arg = phase_angle_deg
        icy_corr = 0.0323  * arg - .00066 * arg**2                ! ganymede

        arg = phase_angle_deg / 100.
        rocky_corr  = 4.98 * arg - 4.88 * arg**2 + 3.02 * arg**3  ! mercury
        cloudy_corr = 0.09 * arg + 2.39 * arg**2 - 0.65 * arg**3  ! venus

!       earth
        if(iplan .eq. 1)then 
            arg = phase_angle_deg / 100.
            if(isat .eq. 0)
     1        phase_corr = 1.30 * arg - 0.19 * arg**2 + 0.48 * arg**3
            if(isat .eq. 1)then
              if(.true.)then
                phase_corr = 3.05 * arg - 1.02 * arg**2 + 1.05 * arg**3
!               phase2 = 10.1 - ((180. - phase_angle_deg) * 0.11)
                phase2 = 10.1 - ((180. - phase_angle_deg) * 0.26)
                phase_corr = max(phase_corr,phase2)
              else ! shevchenko (1980)
                phase_corr = 3.38 * arg - 1.07 * arg**2 + 0.99 * arg**3
              endif
            endif
        endif

!       mercury
        if(iplan .eq. 2)then 
            arg = phase_angle_deg / 100.
!           phase_corr = 3.80 * arg - 2.73 * arg**2 + 2.00 * arg**3 ! orig
            phase_corr = 4.98 * arg - 4.88 * arg**2 + 3.02 * arg**3 ! hilton 2005
        endif

!       venus
        if(iplan .eq. 3)then 
            arg = phase_angle_deg / 100.

!           original
!           phase_corr = 0.09 * arg + 2.39 * arg**2 - 0.65 * arg**3

!           hilton ~2005
!           if(phase_angle_deg .le. 163.6)then
!               phase_corr = 1.03 * arg + 0.57 * arg**2 + 0.13 * arg**3
!           else
!               phase_corr = 5.45 - 1.02 * arg
!           endif

!           mallama 2006 (smoother inflection by s. albers 2017)
            arg = phase_angle_deg
            aint = 163.72d0
            ahw = 4d0
            acoeff = .03d0 / ahw
            if(phase_angle_deg .le. aint)then
                phase_corr =       -1.04396e-3*arg +3.68682e-4*arg**2
     1                                             -2.81374e-6*arg**3
     1                                             +8.93796e-9*arg**4
                if(phase_angle_deg .ge. aint-ahw)then
                   phase_corr =
     1             phase_corr - para(phase_angle_deg,aint-ahw)*acoeff
                endif
            else
                phase_corr = 236.058  -2.81914*arg +8.39034e-3*arg**2
                phase_corr = phase_corr + 4.38
                if(phase_angle_deg .le. aint+ahw)then
                   phase_corr =
     1             phase_corr - para(phase_angle_deg,aint+ahw)*acoeff
                endif
            endif

!           write(6,*)'phase/phase_corr',phase_angle_deg,phase_corr
        endif

!       mars
        if(iplan .eq. 4)then
            if(isat .eq. 0)phase_corr = .016 * phase_angle_deg
!           if(isat .ge. 1)phase_corr = rocky_corr
        endif

!       jupiter
        if(iplan .eq. 5)then
            arg = phase_angle_deg
            if(isat .eq. 0)phase_corr = .005 * phase_angle_deg
            if(isat .eq. 1)phase_corr = 0.046  * arg - .0010  * arg**2 ! io
            if(isat .eq. 2)phase_corr = 0.0312 * arg - .00125 * arg**2 ! europa
            if(isat .eq. 3)phase_corr = 0.0323 * arg - .00066 * arg**2 ! ganymede
            if(isat .eq. 4)phase_corr = 0.078  * arg - .00274 * arg**2 ! callisto
            if(isat .ge. 5)phase_corr = rocky_corr
        endif

!       saturn
        if(iplan .eq. 6)then
            if(isat .eq. 0)then
                phase_corr = .044 * phase_angle_deg

            else if(isat .eq. 6)then ! titan
                phase_corr = cloudy_corr

            else if(isat .eq. 8)then ! iapetus
                phase_corr = 0.5 * (rocky_corr + icy_corr)

            else if(isat .le. 18 .and. isat .ne. 9)then
                phase_corr = icy_corr
                
            else if(isat .gt. 18 .or. isat .eq. 9)then
                phase_corr = rocky_corr
                
            endif
        endif

!       uranus
        if(iplan .eq. 7)then
            if(isat .eq. 0)then
                phase_corr = cloudy_corr
            else if(isat .ge. 1 .and. isat .le. 5)then
                phase_corr = icy_corr
            else if(isat .ge. 6)then
                phase_corr = rocky_corr
            endif
        endif

!       neptune
        if(iplan .eq. 8)then
            if(isat .eq. 0)then
                phase_corr = cloudy_corr
            else if(isat .eq. 1)then
                phase_corr = icy_corr
            else if(isat .ge. 2)then
                phase_corr = rocky_corr
            endif
        endif

!       pluto
        if(iplan .eq. 9)then
            if(isat .le. 1)then
                phase_corr = icy_corr
            else if(isat .ge. 2)then
                phase_corr = rocky_corr                
            endif
        endif

        return
        end


        subroutine correct_near(radii,phase_angle_deg,mag_orig,mag_new)       

        include 'trigd.inc'

        implicit real*8(a-z)

        include '../../include/astparms.for'

        dradii = 1.0 / radii   ! distance / radius

        if(dradii .lt. 1.0)then 
            mag_new = 1e37
            return
        endif

        if(dradii .le. 1d3)then
            mag_new = mag_orig
            return
        endif

        sin_phase_angle = sind(phase_angle_deg)

        if(sin_phase_angle .eq. 0.)then
            mag_new = mag_orig
            return
        endif

        frac_ill = 0.5 + (cosd(phase_angle_deg)/2.)

!       determine (distance/radius) within which terminator is invisible
        r_term_vis = 1.0 / sin_phase_angle

!       assume illuminated frac varies approximately as inverse distance
        if(dradii .le. r_term_vis)then ! terminator is invisible
            if(frac_ill .lt. .50)then
                mag_new = 1e37
                return
            else
                frac_ill_new = 1.0
            endif

        else                           ! terminator is visible
            arg_ill = r_term_vis / dradii ! zero means terminator is visible
                                          ! one  means terminator is invisible        

            if(frac_ill .lt. .50)then
                frac_ill_new = (1.0 - arg_ill) * frac_ill + arg_ill * 0.         
            else
                frac_ill_new = (1.0 - arg_ill) * frac_ill + arg_ill * 1.
            endif

        endif        

        ratio = frac_ill_new / frac_ill

        mag_corr = -(log10(ratio)) * 2.5

        mag_new = mag_orig + mag_corr

        return
        end
