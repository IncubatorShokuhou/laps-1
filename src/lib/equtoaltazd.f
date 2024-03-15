
        subroutine equ_to_altaz_d(dec,ha,phi,alt,az)

        include 'trigd.inc'

        implicit real(a-z)

        sindec = sind(dec)
        cosdec = cosd(dec)
        sinphi = sind(phi)
        cosphi = cosd(phi)
        cosha  = cosd(ha)

        alt=asind (sinphi*sindec+cosphi*cosdec*cosha)
        cosarg = (cosphi*sindec-sinphi*cosdec*cosha)/cosd(alt)
        cosarg = min(max(cosarg,-1.),+1.)
        az =acosd(cosarg)

        if(ha .gt. 0. .and. ha .lt. 180.)az = 360.0 - az

        return
        end

