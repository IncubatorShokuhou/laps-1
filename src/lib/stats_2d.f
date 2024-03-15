cdis
c
c
        subroutine stats_2d(ni,nj,x,y,wt_a,a_t,b_t,xbar,ybar,r
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
c               steve albers    2018            added weights
c
c       inputs/outputs:
c
c          variable     var type    i/o   description
c         ----------   ----------  ----- -------------
c          num_sfc         i         i    number of surface stations.
c          x               ra        i    
c          y               ra        i    
c          wt_a            ra        i    weights
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
        real x(ni,nj), y(ni,nj), wt_a(ni,nj)           

c
c
c.....set up storage variables.
c
        write(6,*)' start stats_2d'

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
        do i=1,ni
        do j=1,nj
          if(x(i,j).eq.badflag .or. y(i,j).eq.badflag) go to 10
!         write(6,5)i,y(i,j),x(i,j)
5         format(i6,2f8.3)
          wt = wt_a(i,j)
          sumxy = (x(i,j) * y(i,j)) * wt + sumxy
          sumx = x(i,j) * wt + sumx
          sumx2 = (x(i,j) * x(i,j)) * wt + sumx2
          sumy = y(i,j) * wt + sumy
          cnt = cnt + wt
          istatus = 1
 10       continue
        enddo
        enddo

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
        do i=1,ni
        do j=1,nj
            if(x(i,j).ne.badflag .and. y(i,j).ne.badflag) then     
                wt = wt_a(i,j)
                cnt = cnt + wt
                sumsq = sumsq + wt * (y(i,j)-x(i,j))**2
            endif
        enddo ! i
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
        do i=1,ni
        do j=1,nj
            if(x(i,j).ne.badflag .and. y(i,j).ne.badflag) then     
                wt = wt_a(i,j)
                sum1 = sum1 + wt * ( (x(i,j) - xbar) * (y(i,j) - ybar) )
                sum2 = sum2 + wt * (x(i,j) - xbar)**2
                sum3 = sum3 + wt * (y(i,j) - ybar)**2
            endif
        enddo ! i
        enddo ! i
        denom = sqrt(sum2) * sqrt(sum3)

        if(denom .ne. 0.)then
            r = sum1 / denom
        else
            r = 0.
        endif

        write(6,*)' regression sums = ',sum1,sum2,sum3
        write(6,*)' ratio of points = ',cnt/(float(ni*nj))

        write(6,900)int(cnt),bias,std,r
900     format(1x,' n/bias/rms/r = ',i8,2f9.2,f9.3)

        if(istatus .eq. 1)then
          write(6,*)' success in stats_2d'
        else
          write(6,*)' warning: istatus in stats_2d =',istatus 
        endif
        write(6,*)

        return
        end

