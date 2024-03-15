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
       subroutine calc_evap(imax,jmax,
     &                       laps_u,
     &                       laps_v,
     &                       laps_t,
     &                       laps_td,
     &                       laps_evap,
     &                       istatus)

c       subroutine to calculate the pan evaporation rates for the 
c       laps grid. 10 km every hour
c       uses the assumption of non radiation limiting as in pg 166
c       hydrology for engineers  by linsley, kohler and paulhus
c       2/5/93
c


      integer*4 imax,jmax
      include 'soilm.inc'

      real    laps_u(imax,jmax)
      real    laps_v(imax,jmax)
      real    laps_t(imax,jmax)
      real    laps_td(imax,jmax)
      real    laps_evap(imax,jmax)
      real    val1, val2, val3
      real    windspeed, tempdegf, dewpoint_degf
      integer*4 istatus

      do j = 1, jmax
        do i = 1, imax
           windspeed = (laps_u(i,j)**2 + laps_v(i,j)**2)**0.5 !m/s
           windspeed = windspeed * 53.6979     ! miles per day
           tempdegf = 1.8 * (laps_t(i,j) - 273.15) + 32.0  !deg f
           dewpoint_degf = 1.8 * (laps_td(i,j) - 273.15) + 32.0
           val1 = ( 0.0041 * tempdegf + 0.676) ** 8
           val2 = ( 0.0041 * dewpoint_degf + 0.676) ** 8
           val3 = (0.37 + 0.0041 * windspeed)

           laps_evap(i,j) = (val1 -val2)**0.88 * val3
        enddo
      enddo
      istatus = 1
c     write(6,*) 'completed evaporation calc'

      return
      end
