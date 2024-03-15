cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis
        subroutine comp_laps_maps(ulaps,vlaps,umaps,vmaps,ni,nj,nk
     1                         ,r_missing_data,rms)

        real ulaps(ni,nj,nk),vlaps(ni,nj,nk)
        real umaps(ni,nj,nk),vmaps(ni,nj,nk)

        nobs_tot = 0
        residualu_tot = 0.
        residualv_tot = 0.

        write(6,2)
2       format(/'      comparing laps to maps'/
     1  1x,'   i   j   k      laps            maps    '
     1     ,'       diff')

        do k = 1,nk

            nobs = 0
            residualu = 0.
            residualv = 0.
            rmaxerrsq = 0.
            sumu = 0.
            sumv = 0.

            do j = 1,nj
            do i = 1,ni

                if(ulaps(i,j,k) .ne. r_missing_data
     1          .and. umaps(i,j,k) .ne. r_missing_data )then

                  nobs = nobs + 1

                  diffu = umaps(i,j,k) - ulaps(i,j,k)
                  diffv = vmaps(i,j,k) - vlaps(i,j,k)

                  residualu = residualu + diffu ** 2
                  residualv = residualv + diffv ** 2

                  sumu = sumu + diffu
                  sumv = sumv + diffv

                  if(diffu**2+diffv**2 .gt. rmaxerrsq)then
                      rmaxerrsq = diffu**2+diffv**2
                      imaxerr = i
                      jmaxerr = j
                  endif

!                  write(6,101)i,j,k
!       1               ,grid_laps_u(i,j,k),grid_laps_v(i,j,k),
!       1               umaps(i,j,k),vmaps(i,j,k),diffu,diffv,rmaxerrsq
!101                format(1x,3i4,3(2x,2f7.1),f7.1)

                endif

            enddo ! j
            enddo ! i


!           do the rms just for this level
            if(nobs .gt. 0)then
                rmsu = sqrt(residualu/float(nobs))
                rmsv = sqrt(residualv/float(nobs))
                rmeanu = sumu / float(nobs)
                rmeanv = sumv / float(nobs)
            else
                rmsu = 0.
                rmsv = 0.
                rmeanu = 0.
                rmeanv = 0.                
            endif

            rms  = sqrt(rmsu**2 + rmsv**2)

            if(rmsu**2 - rmeanu**2 .ge. 0.)then
                rmsamu = sqrt(rmsu**2 - rmeanu**2)
            else
                rmsamu = 0.
            endif

            if(rmsv**2 - rmeanv**2 .ge. 0.)then
                rmsamv = sqrt(rmsv**2 - rmeanv**2)
            else
                rmsamv = 0.
            endif

            write(6 ,102)k,nobs,rmsu,rmsv,rms,sqrt(rmaxerrsq),imaxerr
     1                    ,jmaxerr,rmeanu,rmeanv,rmsamu,rmsamv
!           write(15,102)k,nobs,rmsu,rmsv,rms,sqrt(rmaxerrsq),imaxerr
!    1                    ,jmaxerr,rmeanu,rmeanv,rmsamu,rmsamv
102         format(
     1      ' rms: lp/mp (k,n,rmsu,rmsv,rms,mxerr,ubr,vbr,rmamu,rmamv)='       
     1                                      ,i3,i7,3f5.1,f5.1,2i4,4f5.1)

            residualu_tot = residualu_tot + residualu
            residualv_tot = residualv_tot + residualv
            nobs_tot = nobs_tot + nobs


        enddo ! k


!       do the total rms for the 3d array
        if(nobs_tot .gt. 0)then
            rmsu = sqrt(residualu_tot/nobs_tot)
            rmsv = sqrt(residualv_tot/nobs_tot)
        else
            rmsu = 0.
            rmsv = 0.
        endif

        rms  = sqrt(rmsu**2 + rmsv**2)

        k = 0
        write(6,*)' rms for entire array'
        write(6 ,102)k,nobs,rmsu,rmsv,rms
        write(15,102)k,nobs,rmsu,rmsv,rms

        return

        end

