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
c
c
        subroutine stats_1d(num_sfc,x,y,title,a_t,b_t,xbar,ybar
     1                     ,bias,std,badflag,istatus)
c
c*******************************************************************************
c
c       routine to calculate regression coefficients from y/x data
c       (formerly t and elev data).

c       values returned are the coefficients for the regression equations
c       as well as mean of each array, bias, rms   
c
c       changes:
c               p.a. stamus     12-01-88        original (from j. mcginley)
c
c       inputs/outputs:
c
c          variable     var type    i/o   description
c         ----------   ----------  ----- -------------
c          num_sfc         i         i    number of surface stations.
c          x               ra        i    
c          y               ra        i    
c          b_t             r         o    intercept (y = ax + b)
c          a_t             r         o    slope
c          xbar            r         o    mean value of the stations.
c          ybar            r         o    mean value of the gridded field
c
c       user notes:
c
c       1. units are not changed in this routine.
c
c*******************************************************************************
c
        real x(num_sfc), y(num_sfc)           

        character title*(*)                

c
c
c.....  set up storage variables.
c
        cnt = 0.
        sumxy = 0.
        sumx = 0.
        sumy = 0.
        sumx2 = 0.
        sumy2 = 0.
c
c.....  gather sums and then calculate the 'a' and 'b' for the regression
c.....  equation y = az + b, for both the temperature and dew point.  the
c.....  'b' is the intercept with sea level, and the 'a' is the lapse rate.
c.....  'y' represents the y value and 'z' is x/analyzed
c.....  also calculate the mean elevation of the stations.
c
        istatus = 0

!       write(6,*)' i  y   x  '
        do 10 i=1,num_sfc
          if(x(i).eq.badflag .or. y(i).eq.badflag) go to 10
!         write(6,5)i,y(i),x(i)
5         format(i6,2f8.3)
          sumxy = (x(i) * y(i)) + sumxy
          sumx = x(i) + sumx
          sumx2 = (x(i) * x(i)) + sumx2
          sumy = y(i) + sumy
          cnt = cnt + 1.
          istatus = 1
10      continue

        if(cnt .eq. 0.)then
            write(6,*)' no data for stats'
            xbar = 0.
            ybar = 0.
            bias = 0.
            std = 0.
            istatus = 0
            return
        endif

!       slope    
        denominator = (cnt*sumx2 - sumx*sumx)
        if(denominator .ne. 0.)then
            a_t = (cnt*sumxy - sumx*sumy) / denominator                
        else
            a_t = 0.
            istatus = 0
        endif

!       intercept
        b_t = (sumy - a_t * sumx) / cnt
c
        xbar = sumx / cnt
        ybar = sumy / cnt

        write(6,*)' cnt/a_t/b_t (slope / intercept) = '
     1           ,int(cnt),a_t,b_t

        write(6,*)' xbar,ybar = ',xbar,ybar

!       calculate rms (stdev) of the ob-background differences
        cnt = 0
        sumsq = 0.
        do i = 1,num_sfc
            if(x(i).ne.badflag .and. y(i).ne.badflag) then     
                cnt = cnt + 1.
                sumsq = sumsq + (y(i)-x(i))**2
            endif
        enddo ! i

        if(cnt .gt. 0.)then
            std = sqrt(sumsq / cnt)
        else
            std = 0.
        endif

        bias = ybar - xbar

!       if(a_t .lt. 0.1 .or. a_t .gt. 10.)then
!          write(6,*)' warning, slope is ill conditioned'
!          istatus = 0
!       endif
c
c.....  end of routine
c
!       compute correlation coefficient
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        do i = 1,num_sfc
            if(x(i).ne.badflag .and. y(i).ne.badflag) then     
                sum1 = sum1 + ( (x(i) - xbar) * (y(i) - ybar) )
                sum2 = sum2 + (x(i) - xbar)**2
                sum3 = sum3 + (y(i) - ybar)**2
            endif
        enddo ! i
        denom = sqrt(sum2) * sqrt(sum3)

        if(denom .ne. 0.)then
            r = sum1 / denom
        else
            r = 0.
        endif

        write(6,*)' regression sums = ',sum1,sum2,sum3

        write(6,900)title,int(cnt),bias,std,r
900     format(/,2x,a,' n/bias/rms/r = ',i7,2f9.2,f9.3)

        return
        end

c
c
        subroutine regression(num_pts,x,y,title,a_t,b_t,xbar,ybar
     1                       ,bias,std,badflag,istatus)
c
c*******************************************************************************
c
c       routine to calculate regression coefficients from y/x data
c       (similar to above stats_1d routine so far)

c       values returned are the coefficients for the regression equations
c       as well as mean of each array, bias, rms   
c
c       changes:
c               p.a. stamus     12-01-88        original (from j. mcginley)
c
c       inputs/outputs:
c
c          variable     var type    i/o   description
c         ----------   ----------  ----- -------------
c          num_pts         i         i    number of surface stations.
c          x               ra        i    
c          y               ra        i    
c          b_t             r         o    intercept (y = ax + b)
c          a_t             r         o    slope
c          xbar            r         o    mean elevation of the stations.
c
c       user notes:
c
c       1. units are not changed in this routine.
c
c*******************************************************************************
c
        real x(num_pts), y(num_pts)           

        character title*(*)                

c
c
c.....  set up storage variables.
c
        cnt = 0.
        sumxy = 0.
        sumx = 0.
        sumy = 0.
        sumx2 = 0.
        sumy2 = 0.
c
c.....  gather sums and then calculate the 'a' and 'b' for the regression
c.....  equation y = az + b, for both the temperature and dew point.  the
c.....  'b' is the intercept with sea level, and the 'a' is the lapse rate.
c.....  'y' represents the y value and 'z' is x/analyzed
c.....  also calculate the mean elevation of the stations.
c
        istatus = 0

!       write(6,*)' i  y   x  '
        do 10 i=1,num_pts
          if(x(i).eq.badflag .or. y(i).eq.badflag) go to 10
!         write(6,5)i,y(i),x(i)
5         format(i6,2f8.3)
          sumxy = (x(i) * y(i)) + sumxy
          sumx = x(i) + sumx
          sumx2 = (x(i) * x(i)) + sumx2
          sumy = y(i) + sumy
          cnt = cnt + 1.
          istatus = 1
10      continue

        if(cnt .eq. 0.)then
            write(6,*)' no data for stats'
            xbar = 0.
            ybar = 0.
            bias = 0.
            std = 0.
            return
        endif

!       slope    
        denominator = (cnt*sumx2 - sumx*sumx)
        if(denominator .ne. 0.)then
            a_t = (cnt*sumxy - sumx*sumy) / denominator                
        else
            a_t = 0.
            istatus = 0
        endif

!       intercept
        b_t = (sumy - a_t * sumx) / cnt
c
        xbar = sumx / cnt
        ybar = sumy / cnt

        write(6,*)' cnt/a_t/b_t (slope / intercept) = '
     1           ,int(cnt),a_t,b_t

        write(6,*)' xbar,ybar = ',xbar,ybar

!       calculate rms (stdev) of the ob-background differences
        cnt = 0
        sumsq = 0.
        do i = 1,num_pts
            if(x(i).ne.badflag .and. y(i).ne.badflag) then     
                cnt = cnt + 1.
                sumsq = sumsq + (y(i)-x(i))**2
            endif
        enddo ! i

        if(cnt .gt. 0.)then
            std = sqrt(sumsq / cnt)
        else
            std = 0.
        endif

        bias = ybar - xbar

!       if(a_t .lt. 0.1 .or. a_t .gt. 10.)then
!          write(6,*)' warning, slope is ill conditioned'
!          istatus = 0
!       endif
c
c.....  end of routine
c
!       compute correlation coefficient
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        do i = 1,num_pts
            if(x(i).ne.badflag .and. y(i).ne.badflag) then     
                sum1 = sum1 + ( (x(i) - xbar) * (y(i) - ybar) )
                sum2 = sum2 + (x(i) - xbar)**2
                sum3 = sum3 + (y(i) - ybar)**2
            endif
        enddo ! i
        denom = sqrt(sum2) * sqrt(sum3)

        if(denom .ne. 0.)then
            r = sum1 / denom
        else
            r = 0.
        endif

        write(6,*)' regression sums = ',sum1,sum2,sum3

        write(6,900)title,int(cnt),bias,std,r
900     format(/,2x,a,' n/bias/rms/r = ',i5,2f9.2,f9.3)

        return
        end
