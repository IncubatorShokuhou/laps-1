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

      real function psatoz(psain)

c*  this routine converts a pressure in a standard atmosphere in millibars
c*  into a height in meters

      psa = psain
      scale_high = 8600.
      scale_low  = 6000.
      scale_height = scale_low + (scale_high - scale_low) * psa/1000.

!     initial guess
      scale_2 = scale_high * 0.35 + scale_height * 0.65
      scale_heights = log(1012./psa)
      psatoz = scale_heights * scale_2

      icount = 0

10    icount = icount + 1

      error = ztopsa(psatoz)/psa - 1.0

      correction = error * scale_height

c     write(6,1)icount,psa,psatoz,error,correction
1     format(i4,f8.1,f10.2,e15.4,f10.2)

      psatoz = psatoz + correction

      if(icount .ge. 20)then
          write(6,*)' too many iterations in psatoz, input= ',psain
          write(6,1)icount,psa,psatoz,error,correction
          stop
      endif

      if(abs(error) .gt. 1e-6)goto10

c     write(6,1)icount,psa,psatoz,error,correction

      return
      end

