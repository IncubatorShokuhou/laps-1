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
       subroutine set_missing_sndr(image_in,
     &               ndimx,ndimy,
     &               ismsng,
     &               i2_missing_data,
     &               istatus)
c
c
c
c j. smart    11-13-95   include technique to determine bad data in
c                        addition to missing sat pixels
c
        implicit none
c
        integer i,j,n
        integer ii,jj
        integer ndimx
        integer ndimy
        integer istatus
        integer mstatus
        integer istat_status
        integer imiss_status
        integer image_in(ndimx,ndimy)
        integer image_temp(ndimx,ndimy)
        real      data(125)
        real      ave,adev,sdev,var,skew,curt
        integer ismsng
        integer i2_missing_data
c
c note that this quality control step is performed on sounder counts
c
        istat_status=0
        imiss_status=0

        do j=1,ndimy
        do i=1,ndimx
           image_temp(i,j)=0
        enddo
        enddo

        do j=2,ndimy-1
        do i=2,ndimx-1

           if(image_in(i,j).le.ismsng)then
              n=0
              do jj=j-1,j+1
              do ii=i-1,i+1
                 if(image_in(ii,jj).ne.ismsng)then
                    n=n+1
                    data(n)=float(image_in(ii,jj))
                 endif
              enddo
              enddo
              if(n.ge.2)then
                 call moment(data,n,
     &              ave,adev,sdev,var,skew,curt,
     &              mstatus)

                 if(abs(image_in(i,j)-ave).gt.3.0*sdev.or.
     &              image_in(i,j).eq.ismsng)then
                    image_temp(i,j) = ave
                    istat_status=istat_status-1 
                 else
                    image_temp(i,j) = image_in(i,j)
                 endif
              else
                 image_temp(i,j)=i2_missing_data
                 imiss_status=imiss_status-1
              endif

           else
              image_temp(i,j)=image_in(i,j)
           endif

        end do
        end do

        do j=2,ndimy-1
        do i=2,ndimx-1
           image_in(i,j)=image_temp(i,j)
        enddo
        enddo

        do j=1,ndimy
           if(image_in(1,j).eq.0)then
              image_in(1,j)=i2_missing_data
              imiss_status=imiss_status-1
           endif
           if(image_in(ndimx,j).eq.0)then
              image_in(ndimx,j)=i2_missing_data
              imiss_status=imiss_status-1
           endif
        enddo
        do i=1,ndimx
           if(image_in(i,1).eq.0)then
              image_in(i,1)=i2_missing_data
              imiss_status=imiss_status-1
           endif
           if(image_in(i,ndimy).eq.0)then
              image_in(i,ndimy)=i2_missing_data
              imiss_status=imiss_status-1
           endif
        enddo
        write(6,*)'   # reset due to i2_missing: ',imiss_status
        write(6,*)'   # reset due to statistic: ',istat_status
        istatus=istatus+imiss_status+istat_status
        return
        end

