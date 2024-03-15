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
c   subroutine to read in the soil type infomation
c   created 5/2/92
c   chandran subramaniam
c
c
      subroutine soil_in5(imax,jmax,soiltype,istatus) 

      include 'soilm.inc'
      integer imax,jmax
      integer nf
      integer soiltype(imax,jmax)
      real    r_missing_data
      real,allocatable::static_stl(:,:)
      character*150 c_dir
      character*256 filename
      character*3   var
      character*150 directory
      character*31  ext
      character*10  units
      character*125 comment
c
c current categories for top layer soil texture types
c in static file
c ------------------------------
c  1          sand
c  2          loamy sand
c  3          sandy loam
c  4          silt loam
c  5          silt
c  6          loam
c  7          sandy clay loam
c  8          silty clay loam
c  9          clay loam
c 10          sandy clay
c 11          silty clay
c 12          clay
c 13          organic materials
c 14          water
c 15          bedrock
c 16          other (land-ice)
c
c expected texture types for this lsm
c 1. loamy sand:
c 2. sandy loam:
c 3. loam: 
c 4. sandy clay loam:
c 5. silty clay loam:
c 6. silty clay
c
c mapping between 16 category and 6 category.
c
c if type = 0 then type = 5 ! by jrs. cannot have type = 0;
c                             default to original setting
c if type16 = 1 then type6 = 1
c if type16 = 2 then type6 = 1
c if type16 = 3 then type6 = 2
c if type16 = 4 then type6 = 5
c if type16 = 5 then type6 = 4
c if type16 = 6 then type6 = 3
c if type16 = 7 then type6 = 4
c if type16 = 8 then type6 = 5
c if type16 = 9 then type6 = 5
c if type16 =10 then type6 = 6
c if type16 =11 then type6 = 6
c if type16 =12 then type6 = 6
c if type16 =13 then type6 = 1
c if type16 =14 then type6 = 5 !absurd but was = 5 prior to this
c if type16 =15 then type6 = 5 !absurd but was = 5 prior to this
c if type16 =16 then type6 = 5 !absurd but was = 5 prior to this
c

      istatus = -1

      allocate(static_stl(imax,jmax))
      ext='static'
      call get_directory(ext,c_dir,lend)
      call get_r_missing_data(r_missing_data,istatus)

c     filename=c_dir(1:lend)//'soil/soils.dat'
c     open(unit = 2, file = filename, status = 'old',
c    1  access = 'sequential', iostat = ierr, err = 664)
c     do j = 1 , jmax
c        read(2,*) (soiltype(i,j), i = 1, imax)
c     enddo
c     close(2)
c     nf = index(filename,' ')-1
c     write(6,*) 'got soils data from ',filename(1:nf)
c     istatus = 0
c     return

      var='stl'
      ext='nest7grid'
      call rd_laps_static(c_dir,ext,imax,jmax,1,var,units
     .,comment,static_stl,gridspace,istatus)
      if(istatus.ne.1)then
         print*,'error reading static file for soil type'
         return
      endif

c664   write(6,*)'using default soil types'
      print*,' using static stl soil texture '
      do j = 1 , jmax
         do i = 1, imax
            if(static_stl(i,j).eq.0..or.
     +         static_stl(i,j).eq.r_missing_data)then
             soiltype(i,j) = 5. !default to 5 if no type for this grid point
            elseif(static_stl(i,j).eq.1.)then
             soiltype(i,j) = 1.
            elseif(static_stl(i,j).eq.2.)then
             soiltype(i,j) = 1.
            elseif(static_stl(i,j).eq.3.)then
             soiltype(i,j) = 2.
            elseif(static_stl(i,j).eq.4.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.5.)then
             soiltype(i,j) = 4.
            elseif(static_stl(i,j).eq.6.)then
             soiltype(i,j) = 3.
            elseif(static_stl(i,j).eq.7.)then
             soiltype(i,j) = 4.
            elseif(static_stl(i,j).eq.8.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.9.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.10.)then
             soiltype(i,j) = 6.
            elseif(static_stl(i,j).eq.11.)then
             soiltype(i,j) = 6.
            elseif(static_stl(i,j).eq.12.)then
             soiltype(i,j) = 6.
            elseif(static_stl(i,j).eq.13.)then
             soiltype(i,j) = 1.
            elseif(static_stl(i,j).eq.14.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.15.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.16.)then
             soiltype(i,j) = 5.
            endif
         enddo
      enddo

      do j=1,jmax
      do i=1,imax
         if(soiltype(i,j).eq.0)then
           print*,'i/j= ',i,j,'soil type = 0 at i/j'
           print*,'soil type value = 0 in soil_in5.f'
           print*,'at i/j = ',i,j
           print*,'!! terminating !!'
           return
         endif
      enddo
      enddo

      deallocate (static_stl)

      istatus = 0

      return
      end
