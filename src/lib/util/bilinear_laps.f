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

        subroutine bilinear_laps(ri,rj,imax,jmax,array_2d,result)

cdoc    interpolate 2-d array to find the field value at a fractional grid
cdoc    point.

        real array_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in bilinear_laps: stop'
            stop
        endif

        i = int(ri)
        if(i .eq. imax)i=i-1

        j = int(rj)
        if(j .eq. jmax)j=j-1

        if(i .ge. 1 .and. i .le. imax .and.
     1     j .ge. 1 .and. j .le. jmax) then

            fraci = ri - i
            fracj = rj - j

            z1=array_2d(i  , j  )
            z2=array_2d(i+1, j  )
            z3=array_2d(i+1, j+1)
            z4=array_2d(i  , j+1)

            if(  z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                result =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                - (z2+z4-z3-z1)*fraci*fracj

            else
                result = r_missing_data

            endif

        else
            result = r_missing_data

        endif

        return
        end

        subroutine bilinear_interp_extrap(ri,rj,imax,jmax
     1                                   ,array_2d,result,istatus)

cdoc    interpolate 2-d array to find the field value at a fractional grid
cdoc    point. this one allows you to extrapolate very slightly outside the 
cdoc    grid.

        real array_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in bilinear_interp_extrap'
            return
        endif

        j = nint(rj)
        i = nint(ri)
        if(j .lt. jmax .and. j .gt. 1 .and.
     &     i .lt. imax .and. i .gt. 1)then

           fraci = ri - i
           fracj = rj - j
           if(j.gt.rj)fracj=j-rj
           if(i.gt.ri)fraci=i-ri

c standard bilinear interpolation

              z1=array_2d(i  , j  )
              z2=array_2d(i+1, j  )
              z3=array_2d(i+1, j+1)
              z4=array_2d(i  , j+1)

              if(    z1 .ne. r_missing_data
     1         .and. z2 .ne. r_missing_data
     1         .and. z3 .ne. r_missing_data
     1         .and. z4 .ne. r_missing_data)then

                  result =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                   - (z2+z4-z3-z1)*fraci*fracj
              else
                  result = r_missing_data
              endif

        elseif(j .gt. jmax .or. j .lt. 1 .or.
     &         i .gt. imax .or. i .lt. 1)then
                  result = r_missing_data

        elseif(j .eq. jmax .or. j .eq. 1)then
              if(i .eq. imax .or. i .eq. 1)then
                 result = array_2d(i,j)
              else
                 frac1 = 1-(ri-int(ri))
                 frac2 = 1-frac1
                 z1=array_2d(i  , j  )
                 z2=array_2d(i+1, j  )
                 result = z1*frac1+z2*frac2
              endif

        elseif(i .eq. imax .or. i .eq. 1)then
              if(j .eq. jmax .or. j .eq. 1)then
                 result = array_2d(i,j)
              else
                 frac1 = 1-(rj-int(rj))
                 frac2 = 1-frac1
                 z1=array_2d(i  , j  )
                 z2=array_2d(i, j+1  )
                 result = z1*frac1+z2*frac2
              endif

        endif

        istatus = 1

        return
        end

        subroutine bilinear_laps_3do(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                             ,nz_laps,array_3d,result)

cdoc    interpolate 3-d array to find the field values at fractional grid
cdoc    points.

        real array_3d(imax,jmax,nz_laps)                 ! i
        real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! i
        real result(nx_laps,ny_laps,nz_laps)             ! o

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in bilinear_laps: stop'
            stop
        endif

        do k  = 1,nz_laps
        do il = 1,nx_laps
        do jl = 1,ny_laps

          ri = ri_a(il,jl)

          i = int(ri)
          if(i .eq. imax)i=i-1

          rj = rj_a(il,jl)

          j = int(rj)
          if(j .eq. jmax)j=j-1

          if(i .ge. 1 .and. i .le. imax .and.
     1       j .ge. 1 .and. j .le. jmax) then

            fraci = ri - i
            fracj = rj - j

            z1=array_3d(i  , j  ,k)
            z2=array_3d(i+1, j  ,k)
            z3=array_3d(i+1, j+1,k)
            z4=array_3d(i  , j+1,k)

            if(  z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                result(il,jl,k) =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                          - (z2+z4-z3-z1)*fraci*fracj

            else
                result(il,jl,k) = r_missing_data

            endif

          else
            result(il,jl,k) = r_missing_data

          endif

        enddo ! j
        enddo ! i  
        enddo ! k

        return
        end
 
        subroutine bilinear_laps_3df(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                             ,nz_laps,array_3d,result)

cdoc    interpolate 3-d array to find the field values at fractional grid
cdoc    points.

        real array_3d(imax,jmax,nz_laps)                 ! i
        real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! i
        real result(nx_laps,ny_laps,nz_laps)             ! o

        real z1(nz_laps),z2(nz_laps),z3(nz_laps),z4(nz_laps)

        write(6,*)' subroutine bilinear_laps_3df...'

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in bilinear_laps: stop'
            stop
        endif

        do jl = 1,ny_laps
        do il = 1,nx_laps

          ri = ri_a(il,jl)

          i = int(ri)
          if(i .eq. imax)i=i-1

          rj = rj_a(il,jl)

          j = int(rj)
          if(j .eq. jmax)j=j-1

          if(i .ge. 1 .and. i .le. imax .and.
     1       j .ge. 1 .and. j .le. jmax) then

            fraci = ri - i
            fracj = rj - j
            fracij = fraci*fracj

            z1(:)=array_3d(i  , j  ,:)
            z2(:)=array_3d(i+1, j  ,:)
            z3(:)=array_3d(i+1, j+1,:)
            z4(:)=array_3d(i  , j+1,:)

!           result(il,jl,:) 
!    1                  =  z1(:)+(z2(:)-z1(:))*fraci+(z4(:)-z1(:))*fracj
!    1                  - (z2(:)+z4(:)-z3(:)-z1(:))*fracij

            result(il,jl,:) 
     1                  =  z1(:) * (1.0 - fraci - fracj + fracij)             
     1                  +  z2(:) * (      fraci         - fracij) 
     1                  +  z3(:) * (                      fracij) 
     1                  +  z4(:) * (              fracj - fracij) 

            where(abs(result(il,jl,:)) .gt. 1e10)
     1          result(il,jl,:) = r_missing_data
  
          else
            result(il,jl,:) = r_missing_data

          endif

        enddo ! i  
        enddo ! j

        return
        end

        subroutine bilinear_laps_3d(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                             ,nz_laps,array_3d,result)

cdoc    interpolate 3-d array to find the field values at fractional grid
cdoc    points. hopefully it's more efficient to call bilinear_laps_2d.

        real array_3d(imax,jmax,nz_laps)                 ! i
        real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! i
        real result(nx_laps,ny_laps,nz_laps)             ! o

        do k = 1,nz_laps
           call bilinear_laps_2d(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                          ,array_3d(1,1,k),result(1,1,k))
        enddo ! k

        return
        end

        subroutine bilinear_laps_2d(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                             ,array_2d,result)

cdoc    interpolate 2-d array to find the field values at fractional grid
cdoc    points.

        real array_2d(imax,jmax)                         ! i
        real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! i
        real result(nx_laps,ny_laps)                     ! o

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error in bilinear_laps_2d: stop'
            stop
        endif

        do il = 1,nx_laps
        do jl = 1,ny_laps

          ri = ri_a(il,jl)

          i = int(ri)
          if(i .eq. imax)i=i-1

          rj = rj_a(il,jl)

          j = int(rj)
          if(j .eq. jmax)j=j-1

          if(i .ge. 1 .and. i .le. imax .and.
     1       j .ge. 1 .and. j .le. jmax) then

            fraci = ri - i
            fracj = rj - j

            z1=array_2d(i  , j  )
            z2=array_2d(i+1, j  )
            z3=array_2d(i+1, j+1)
            z4=array_2d(i  , j+1)

            if(  z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                result(il,jl) =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                          - (z2+z4-z3-z1)*fraci*fracj

            else
                result(il,jl) = r_missing_data

            endif

          else
            result(il,jl) = r_missing_data

          endif

        enddo ! j
        enddo ! i  

        return
        end
