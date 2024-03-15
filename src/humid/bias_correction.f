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
      function bias_correction (raw_data, satellite, sounder, 
     1     channel_index)

      implicit none


c     routine to fit goes bias corrections into variational processing
c     note that optran is too warm.  to make data appear as optran values,
c     add the bias.
c     author birkenheuer     9/29/99

c     optran interface code
c     chief application gimpap/afwa

      real bias_correction
      integer channel_index
      integer satellite         ! 8= goes 8,  10 = goes 10 etc
      integer sounder           ! 1= sounder, 0 = imager
      real bias_8 (22), bias_10(22),bias_9(22), bias_12(22)
      namelist /bias_coefficients_nl/ bias_8,bias_10,bias_9,bias_12
      save bias_8, bias_10,bias_9, bias_12, first_time
      real raw_data
      character*200 fname
      integer len
      integer first_time
      data first_time /0/




c     code***********************************************************

      bias_correction = raw_data

      if (first_time .eq. 0) then !get coef

c     read in namelist
         call get_directory('static',fname,len)
         open (23, file=fname(1:len)//'bias_coefficients.nl',
     1        status = 'old', err = 24)

         read(23,bias_coefficients_nl,end=24)
         close (23)

         first_time = 1         ! set to not first time
      endif


c     end read namelist



      if(satellite .eq. 8) then
         if(sounder .eq. 1) then ! correct raw_data to match optran
            bias_correction = raw_data + bias_8(channel_index)
         endif
      endif

      
      if(satellite .eq.10) then
         if(sounder.eq.1) then
            bias_correction = raw_data + bias_10(channel_index)
         endif
      endif

      return

 24   write(6,*) 'warning no bias coefficients found'

      return

      end

