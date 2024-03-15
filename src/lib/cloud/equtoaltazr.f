
        subroutine equ_to_altaz_r(dec,ha,phi,alt,az)

        implicit real*8(a-z)

        sindec = dsin(dec)
        cosdec = dcos(dec)
        sinphi = dsin(phi)
        cosphi = dcos(phi)
        cosha  = dcos(ha)

        alt=asin (sinphi*sindec+cosphi*cosdec*cosha)
        az =acos((cosphi*sindec-sinphi*cosdec*cosha)/dcos(alt))

        if(ha .gt. 0 .and. ha .lt. 3.1415926535897932d0)then !    0 to +180
            az = 6.2831853071796d0 - az
        endif

        if(ha .lt. -3.1415926535897932d0)then                ! -360 to -180
            az = 6.2831853071796d0 - az
        endif

        return
        end

        subroutine angcoord_2d(c_conv,ni,nj,arg1,arg2,lat1,lon1,phi
     1                                               ,lat2,lon2)

!       general purpose angular coordinate transformations in 2d
!       this version works with single precision arguments in degrees

        include 'trigd.inc'

        parameter (pi = 3.1415926535897932e0)
        parameter (rpd = pi/180.)

        angdif(x,y)=mod(x-y+9.4247779607694e0,6.2831853071796e0)
     1                     -3.1415926535897932e0

        atan3df(x,y) = mod((atan2(x,y)/rpd+360e0),360e0)

        character*20 c_conv

        real lat1(ni,nj),lon1(ni,nj),lat2(ni,nj),lon2(ni,nj) ! degrees
        real lat2s,lon2arg,lst,lstmra,ll0,l0

        write(6,*)' subroutine angcoord_2d ',c_conv

        if(trim(c_conv) .eq. 'altaz2decra')then
          lst = arg1

          do i = 1,ni
          do j = 1,nj
            z = 90. - lat1(i,j)
            azi = lon1(i,j) - 180. ! counted from meridian

            sindec = cosd(z)*sind(phi)-sind(z)*cosd(phi)*cosd(azi)
            sindec = max(min(sindec,1.0),-1.0)

            cosha = (cosd(z)*cosd(phi)+sind(z)*sind(phi)*cosd(azi))
     1                 /cosd(dec)
            cosha = min(max(cosha,-1.),+1.)

            cosdec = sqrt(1.-sindec**2)
            if(cosdec .ne. 0.)then
              sinha = (sind(azi)*sind(z)) / cosdec
            else
              sinha = 0.
            endif

            dec = asind(sindec)
            ha = atan3df(sinha,cosha)
            ra = lst - ha
           
            lat2(i,j) = dec
            lon2(i,j) = ra

          enddo ! j
          enddo ! i

        elseif(trim(c_conv) .eq. 'decra2galactic')then
          ra0 = 282.86 ! ra of ascending node of galactic equator 
          l0 = 32.93 ! galactic longitude of ascending node
          decngp = 27.13

          do i = 1,ni
          do j = 1,nj
            dec = lat1(i,j)
            ra = lon1(i,j)

            sinb = sind(dec)*sind(decngp) 
     1           - cosd(dec)*cosd(decngp)*sind(ra-ra0)
            sinb = max(min(sinb,1.0),-1.0)

            cosb = sqrt(1.-sinb**2)
            cosll0 = cosd(ra-ra0)*cosd(dec)/cosb
            sinll0 = (sind(dec)*cosd(decngp)
     1             + cosd(dec)*sind(decngp)*sind(ra-ra0)) / cosb
            ll0 = atan3df(sinll0,cosll0)

            lat2(i,j) = asind(sinb)
            lon2(i,j) = ll0 + l0
          enddo ! j
          enddo ! i
            
        elseif(trim(c_conv) .eq. 'decra2helioecliptic')then
          obl = 23.446
          sollam = arg1
          do i = 1,ni
          do j = 1,nj
            dec = lat1(i,j)
            ra = lon1(i,j)
            sinb = sind(dec)*cosd(obl)-cosd(dec)*sind(obl)*sind(ra)
            sinb = max(min(sinb,1.0),-1.0)
            cosb = sqrt(1.-sinb**2)
            coslam = cosd(ra)*cosd(dec)/cosb
            sinlam = (sind(dec)*sind(obl)+cosd(dec)*cosd(obl)*sind(ra))
     1                         /cosb
            lat2(i,j) = asind(sinb)
            lon2(i,j) = atan3df(sinlam,coslam) - sollam
          enddo ! j
          enddo ! i
        endif

        where(lon2 .lt. 0.)
            lon2 = lon2 + 360.
        end where

        where(lon2 .gt. 360.)
            lon2 = lon2 - 360.
        end where

        return
        end
