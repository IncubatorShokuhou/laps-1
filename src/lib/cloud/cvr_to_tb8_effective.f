
        subroutine cvr_to_tb8_effective(kcld,temp_3d,klaps,i,j,ni,nj,a
     1                                 ,f,ilyr,cvr,cld_hts,t_gnd_k
     1                                 ,heights_3d,t_effective,nlyr
     1                                 ,idebug,istatus)

!       forward model for multiple cloud layers

        real thresh_cvr
        parameter (thresh_cvr = 0.1)

        integer max_cld_lyrs                   
        parameter (max_cld_lyrs = 100)

        real a(max_cld_lyrs)     ! cloud fractions of layers
        real f(max_cld_lyrs)     ! apparent "cross-section" of cloud layers seen from above
        integer ik(max_cld_lyrs) ! height level representative of cloud layers
        integer ilyr(kcld)       ! layer index for each cloud lvl (needs kcloud)
        real cvr(kcld)           ! cloud cover from analysis
        real temp_3d(ni,nj,klaps),temp_lyr(max_cld_lyrs)
        real heights_3d(ni,nj,klaps)
        real cld_hts(kcld)

        if(idebug .eq. 1)then
            write(6,*)' subroutine cvr_to_tb8_effective at:',i,j
            write(6,*)'n,temp_lyr(n),f(n),rsum...'
        endif

!       convert from cloud cover to discreet cloud layer indices (cvr to a)
        nlyr = 0
        ilyr = 0
        do k = kcld-1,1,-1
            if(cvr(k)   .ge. thresh_cvr .and.
     1         cvr(k+1) .lt. thresh_cvr)then      ! top of layer
                nlyr = nlyr + 1
                a(nlyr) = cvr(k)
                ik(nlyr) = k
                sum_lyr = cvr(k)
                nlevels_lyr = 1

             else 
                if(nlyr .ge. 1)then
                    if(cvr(k) .gt. a(nlyr))then   ! max within layer
!                       a(nlyr) = cvr(k)
                        ik(nlyr) = k
                        sum_lyr = sum_lyr + cvr(k)
                        nlevels_lyr = nlevels_lyr + 1
                        a(nlyr) = sum_lyr / float(nlevels_lyr) ! update mean
                    endif
                endif
             endif

             if(cvr(k) .ge. thresh_cvr)then       ! still within layer
                 ilyr(k) = nlyr
             else                                 ! below layer
                 ilyr(k) = 0
             endif

        enddo ! k

!       get temperatures of the layers
        do n = 1,nlyr
            k = ik(n)
            z_temp = height_to_zcoord2(cld_hts(k),heights_3d,ni,nj,klaps
     1                                          ,i,j,istatus)

            if(istatus .ne. 1)then ! layer is above height grid?
                write(6,*)' note: istatus=0 in cvr_to_tb8_effective'       
                write(6,*)' i,j,k,cld_hts(k)',i,j,k,cld_hts(k)
                return
            endif

            z_temp = max(1.,min(z_temp,float(klaps)-.001))
            iz_temp = int(z_temp)
            frac = z_temp - iz_temp
            temp_lyr(n) = temp_3d(i,j,iz_temp)    * (1. - frac)
     1                 +  temp_3d(i,j,iz_temp+1)  * frac
        enddo ! n

!       add a layer for the ground
        nlyr = nlyr + 1
        a(nlyr) = 1.0
        temp_lyr(nlyr) = t_gnd_k

!       convert cloud layer fractions to "cross-section" seen from satellite
!       this solves for the f array given the a array
        a(1) = min(a(1),1.0)
        f(1) = a(1)
        sumf = f(1)
        if(nlyr .ge. 2)then
            do n = 2,nlyr
                a(n) = min(a(n),1.0)
                f(n) = a(n) * (1.0 - sumf)
                sumf = sumf + f(n)
            enddo ! n
        endif ! nlyr

!       calculate total radiance from all cloud layers + ground
        rsum = 0
        do n = 1,nlyr
!           if(idebug .eq. 1)then
!               write(6,*)n,temp_lyr(n),f(n),' call temp_to_rad'
!           endif
            rsum = rsum + temp_to_rad(temp_lyr(n)) * f(n)
            if(idebug .eq. 1)then
                write(6,*)n,temp_lyr(n),f(n),rsum
            endif
        enddo ! n

!       convert to effective temperature and compare to observed brightness temp
        t_effective = rad_to_temp(rsum)

        istatus = 1
        return
        end

        function temp_to_rad(temp)

        data init /0/
        save init

        data iwrite /0/
        save iwrite

        if(init .eq. 0)then
            init = 1
            nsat = 3
            write(6,*)' call plnkiv'
            call plnkiv(nsat)
        endif

        temp_to_rad = temp
        temp_to_rad = vplanc(temp,8)

        if(iwrite .eq. 0)then
            iwrite = 1
            write(6,*)' temp,temp_to_rad ',temp,temp_to_rad
        endif

        return
        end

        function rad_to_temp(rad)

        data init /0/
        save init

        data iwrite /0/
        save iwrite

        if(init .eq. 0)then
            init = 1
            nsat = 3
            call plnkiv(nsat)
        endif

        rad_to_temp = rad
        rad_to_temp = vbrite(rad,8)

        if(iwrite .eq. 0)then
            iwrite = 1
            write(6,*)' rad,rad_to_temp ',rad,rad_to_temp
        endif

        return
        end
