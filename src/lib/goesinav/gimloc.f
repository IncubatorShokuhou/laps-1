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

c***********************************************************************
c**   noaa/nesdis/socc/software branch and isi                                  
c***********************************************************************        
c**                                                                             
c**   project     operations ground equipment for goes-next                     
c**   system      earth location users guide                                    
c**   routine     gimloc                                                        
c**   source      f.gimloc                                                      
c**   load name   any                                                           
c**   programmer  thomas i. babicki                                             
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  --   ---------------------------------------------        
c**    a     5/01/89  tb   initial creation(socc/isi)                           
c**
c**    b     2/19/93  jh   add access to both imager and sounder sets
c**
c***********************************************************************        

c
c aws -- i eliminated the program code here...
c

c***********************************************************************        
c**                                                                             
c**   noaa/nesdis/socc/software branch                                          
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : time50                                                        
c**   source    : f.time50                                                      
c**   load name : any                                                           
c**   programmer: thomas i. babicki                                             
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a      2/17/89  tb   initial creation                                     
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   time50 accepts two words containing date and time                         
c**   and returns time expressed as double precision minutes from               
c**   1950 jan. 1.0                                                             
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: none                                                    
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
      function time50(i)                                                        
c                                                                               
c     calling parameters                                                        
c                                                                               
      integer*4 i(2)                                                            
c                                                                               
c     local variables                                                           
c                                                                               
      integer*4 ny                                                              
c                    year                                                       
      integer*4 nd                                                              
c                    day of year                                                
      integer*4 nh                                                              
c                    hour                                                       
      integer*4 nm                                                              
c                    minute                                                     
      real*8    s                                                               
c                    seconds                                                    
      integer j                                                                 
c                                                                               
c     include files - rec                                                       
c                                                                               
c                                                                               
c***********************************************************************        
c                                                                               
c     convert input year, day of year, hour and minute                          
c     to integer values.                                                        
c                                                                               
      ny=i(1)/10000                                                          
      iaa=i(1)-(ny*10000)                                                    
      nd=(i(1)-(ny*10000))*.1                                                
      iab=(iaa-(nd*10))*10                                                      
      nbc=i(2)/10000000.                                                     
      iac=i(2)-(nbc*10000000)                                                
      nh=iab+nbc                                                                
      def=i(2)-iac                                                           
      nm=iac*.00001                                                             
      s=(i(2)-(def+(nm*100000)))*.001d0
      print 1000,ny,nd,nh,nm,s                                                  
 1000 format (1h1,'year =',i4,/,1x,'jday =',i3,/,1x,                            
     *            'hour =',i2,/,1x,'min  =',i2,/,1x,                            
     *            'sec  =',f6.3)                                                
c                                                                               
c***********************************************************************        
c                                                                               
c     here we convert integer year and day of year to number of                 
c     days from 0 hour ut, 1950 jan. 1.0                                        
c     this convertion is based on an algorithm by fliegel and van               
c     flandern, comm. of acm, vol.11, no. 10, oct. 1968 (p.657)                 
c                                                                               
      j=nd+1461*(ny+4799)/4-3*((ny+4899)/100)/4-2465022                         
c                                                                               
c     compute time in minutes from january 1.0, 1950                            
c                                                                               
      time50=j*1440.d0+nh*60.d0+nm+s/60.d0                                      
      return                                                                    
      end                                                                       
c***********************************************************************        
c**   noaa/nesdis/socc/software branch                                          
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : timex                                                         
c**   source    : f.timex                                                       
c**   load name : any                                                           
c**   programmer: thomas i. babicki                                             
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a      4/25/89  tb   initial creation                                     
c***********************************************************************        
      function timex(ny,nd,nh,nm,s)                                             
c     
      real*8 timex
c                                                                          
c     calling parameters                                                        
c                                                                               
      integer*4 ny                                                              
c                    year                                                       
      integer*4 nd                                                              
c                    day of year                                                
      integer*4 nh                                                              
c                    hour                                                       
      integer*4 nm                                                              
c                    minute                                                     
      real*8    s
c                    seconds                                                    
      integer j                                                                 
c                                                                               
c***********************************************************************        
c                                                                               
      print 1000,ny,nd,nh,nm,s                                                  
 1000 format (/,1x,'year =',i4,/,1x,'jday =',i3,/,1x,                           
     *            'hour =',i2,/,1x,'min  =',i2,/,1x,                            
     *            'sec  =',f6.3)                                                
c                                                                               
c***********************************************************************        
c                                                                               
      j=nd+1461*(ny+4799)/4-3*((ny+4899)/100)/4-2465022                         
c                                                                               
c     compute actual time in minutes from january 1.0, 1950                     
c                                                                               
      timex=j*1440.d0+nh*60.d0+nm+s/60.d0                                       
      return                                                                    
      end                                                                       
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : setcons                                                       
c**   source    : f.setcons                                                     
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a     02/16/89  il   initial creation                                     
c**
c**   b     05/19/94  np   added calculation of instrument elevation and
c**                        scan angle biases based on user input
c***********************************************************************        
c**                                                                             
c**   this subroutine generates constants in common  instcomm                   
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: instco                                                  
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      subroutine setcon(instr,ns_nad_cy,ns_nad_inc,ew_nad_cy,ew_nad_inc)
c                                                                               
c     calling parameters                                                        
c                                                                               
c                                                                               
c     local variables

      integer ns_nad_cy,ns_nad_inc,ew_nad_cy,ew_nad_inc
c                                                                               
c                                                                               
c     include files                                                             
c                                                                               
      include 'elcons.inc'                                                     
        include 'instco.inc'                                                          
c***********************************************************************        
      incmax(1)=6136                                                            
      incmax(2)=2805                                                            
      elvinc(1)=8.0d-6
      elvinc(2)=17.5d-6                                                         
      scninc(1)=16.d-6                                                          
      scninc(2)=35.d-6                                                          
      elvln(1)=28.d-6                                                           
      elvln(2)=280.d-6                                                          
      scnpx(1)=16.d-6                                                           
      scnpx(2)=280.d-6
c     ************************************************************
c     commented out elevation and scan bias constants since instrument
c     earth nadir position is available in gvar data and periodically
c     updated
c
c       
c        elvmax(1)=0.220896d0                         
c        elvmax(2)=0.22089375d0
c        scnmax(1)=0.24544d0
c        scnmax(2)=0.2454375d0

c     recompute elevation and scan biases based on user inputs of                 
c     cycles & increments obtained from gvar

c     elvmax(instr) = (ns_nad_cy*incmax(instr)+ns_nad_inc)*elvinc(instr)

      if(instr.eq.1)then                                                                 
        elvmax(instr)=(ns_nad_cy*incmax(instr)+ns_nad_inc)*elvinc(instr)
      else
        elvmax(instr)=((9-ns_nad_cy)*incmax(instr)-ns_nad_inc)
     +                *elvinc(instr)
      endif
      
      scnmax(instr) = (ew_nad_cy*incmax(instr)+ew_nad_inc)*scninc(instr)
c     ************************************************************
      return                                                                    
      end                                                                       
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : lmodel
c**   source    : f.lmodel                                                      
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   1     01/09/89  il   initial creation                                     
c**   2     06/02/89  il   coordinate axes changed according to                 
c**                        ford's definition in sdaip, drl 504-01               
c**   3     08/21/89  il   corrected orbit angle computations                   
c**   4     03/08/94  sc   s/c compensation applied unconditionally;
c**                        reference radial distance, latitude and
c**                        orbit yaw set to zero if imc disabled.
c**   5     03/08/94  sc   added trap for slat=syaw=0; corrected
c**                        expression for lam.
c***********************************************************************        
c**                                                                             
c**   this subroutine computes the position of the satellite and the            
c**   attitude of the imager or sounder.  the calculations are based            
c**   on the oats orbit and attitude model represented by the o&a               
c**   parameter set in gvar block 0.                                            
c**        inputs:                                                              
c**          time, epoch time, o&a parameter set, imc status.                   
c**                                                                             
c**        outputs:                                                             
c**          the spacecraft position vector in earth fixed coordinates;         
c**          the geometric roll, pitch, yaw angles and the roll,                
c**          pitch misalignments for either the imager or the sounder;          
c**          the earth fixed to instrument frame transformation matrix;         
c**          geographic latitude and longitude at subsatellite point.           
c**                                                                             
c**   description                                                               
c**   lmodel accepts an input double precision time in minutes from             
c**   1950, jan.1.0 and an input set of o&a parameters and computes             
c**   position of the satellite, the attitude angles and attitude               
c**   misalignments and the instrument to earth fixed coordinates               
c**   transformation matrix.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: /elcomm/ xs,q3,pitch,roll,yaw,pma,rma,bt                
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : inst2er,gatt                                            
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      subroutine lmodel(t,tu,rec,imc,rlat,rlon)                                 
c                                                                               
c     calling arguments                                                         
c                                                                               
      real*8 t                                                                  
c                                   input time from 1950, jan 1.0 (min)         
      real*8 tu                                                                 
c                                   epoch time from 1950, jan 1.0 (min)         
      real*8 rec(336)                                                           
c                                   input o&a parameter set                     
      integer imc                                                               
c                                   input imc status: 0 - on, 1 - off           
      real*8 rlat                                                               
c                                   subsatellite geodetic latitude (rad)        
      real*8 rlon                                                               
c                                   subsatellite longitude in radians           
c                                                                               
c     local variables
c                                                                   
      real*8  r                                                                 
c                    normalized satellite distance (in units of kmer9)          
      real*8 ts                                                                 
c                    time from epoch in minutes                                 
      real*8 b(3,3)                                                             
c                    spaccraft to earth fixed coordinates transformation        
c                    matrix                                                     
      real*8 te                                                                 
c                    exponenential time delay from epoch (in minutes)           
      real*8 phi                                                                
c                    subsatellite geocentric latitude in radians                
      real*8 dr                                                                 
c                    radial distance from the nominal (km)                      
      real*8 psi                                                                
c                    orbital yaw (in radians)                                   
      real*8  lam                                                               
c                    imc longitude (in radians)                                 
      real*8  u                                                                 
c                    argument of latitude (in radians)                          
      real*8 su,cu                                                              
c                    sin(u), cos(u)                                             
      real*8 si,ci                                                              
c                    sine and cosine of the orbit inclination                   
      real*8 slat                                                               
c                    sine of geocentric latitude                                
      real*8 asc                                                                
c                    longitude of the ascending node in radians                 
      real*8 sa,ca                                                              
c                    sine and cosine of asc                                     
      real*8 syaw                                                               
c                    sine of the orbit yaw                                      
      real*8 wa                                                                 
c                    solar orbit angle in radians                               
      real*8 w                                                                  
c                    orbit angle in radians                                     
      real*8 sw,cw                                                              
c                    sin(w),  cos(w)                                            
      real*8 s2w,c2w                                                            
c                    sin(2*w),  cos(2*w)                                        
      real*8 sw1,cw1                                                            
c                    sin(0.927*w),  cos(0.927*w)                                
      real*8 sw3,cw3                                                            
c                    sine and cosine of 1.9268*w                                
      real*8 dlat                                                               
c                    change in sine of geocentric latitude                      
      real*8 dyaw                                                               
c                    change in sine of orbit yaw                                
      real*8 gatt                                                               
c                    subroutine function                                        
c      real*8 a1,a2                                                              
c                    work areas                                                 
c                                                                               
c     include files                                                             
c                                                                               
        include 'elcons.inc'                                                          
        include 'elcomm.inc'                                                          
c                                                                               
c***********************************************************************        
c                                                                               
c     assign reference values to the subsatellite longitude and                 
c     latitude, the radial distance and the orbit yaw.                          
c                                                                               
      lam=rec(5)                                                                
      dr=rec(6)                                                                 
      phi=rec(7)                                                                
      psi=rec(8)                                                                
c                                                                               
c     assign reference values to the attitudes and misalignments                
c                                                                               
      roll=rec(9)                                                               
      pitch=rec(10)                                                             
      yaw=rec(11)                                                               
      rma=0.0d0
      pma=0.0d0
c                                                                               
c     if imc is off, compute changes in the satellite orbit                     
c                                                                               
      if (imc.ne.0) then                                                        
c
c     set reference radial distance, latitude and orbit yaw to zero
c
      dr=0.0d0
      phi=0.0d0
      psi=0.0d0
c
c     compute time since epoch (in minutes)                                     
c                                                                               
      ts=t-tu                                                                   
c                                                                               
c     computes orbit angle and the related trigonometric functions.             
c     earth rotational rate=.729115e-4 (rad/s)                                  
c                                                                               
      w=0.729115d-4*60.0d0*ts                                                   
      sw=dsin(w)
      cw=dcos(w)
      sw1=dsin(0.927d0*w)
      cw1=dcos(0.927d0*w)
      s2w=dsin(2.0d0*w)
      c2w=dcos(2.0d0*w)
      sw3=dsin(1.9268d0*w)
      cw3=dcos(1.9268d0*w)
c                                                                               
c     computes change in the imc longitude from the reference                   
c                                                                               
      lam=lam+rec(18)+(rec(19)+rec(20)*w)*w                                     
     1        +(rec(27)*sw1+rec(28)*cw1+rec(21)*sw+rec(22)*cw
     2        +rec(23)*s2w+rec(24)*c2w + rec(25)*sw3+rec(26)*cw3                
     3        +w*(rec(29)*sw+rec(30)*cw))*2.0d0
c                                                                               
c     computes change in radial distance from the reference (km)                
c                                                                               
      dr=dr + rec(31) + rec(32)*cw+rec(33)*sw                                   
     1        +rec(34)*c2w+rec(35)*s2w + rec(36)*cw3+rec(37)*sw3                
     2        +rec(38)*cw1+rec(39)*sw1 + w*(rec(40)*cw+rec(41)*sw)              
c                                                                               
c     computes the sine of the change in the geocentric latitude                
c                                                                               
      dlat=rec(42) + rec(43)*cw+rec(44)*sw                                      
     1        +rec(45)*c2w+rec(46)*s2w                                          
     2        +w*(rec(47)*cw+rec(48)*sw)                                        
     3        +rec(49)*cw1+rec(50)*sw1                                          
c                                                                               
c     computes geocentric latitude by using an expansion for arcsine            
c                                                                               
      phi=phi+dlat*(1.0d0+dlat*dlat/6.0d0)
c                                                                               
c     computes sine of the change in the orbit yaw                              
c                                                                               
      dyaw=rec(51) + rec(52)*sw+rec(53)*cw                                      
     1        +rec(54)*s2w+rec(55)*c2w                                          
     2        +w*(rec(56)*sw+rec(57)*cw)                                        
     3        +rec(58)*sw1+rec(59)*cw1                                          
c                                                                               
c     computes the orbit yaw by using an expansion for arcsine.                 
c                                                                               
      psi=psi+dyaw*(1.0d0+dyaw*dyaw/6.0d0)
c                                                                               
c     calculation of changes in the satellite orbit ends here                   
c                                                                               
      end if                                                                    
c                                                                               
c     conversion of the imc longitude and orbit yaw to the subsatellite         
c     longitude and the orbit inclination (ref: goes-pcc-tm-2473, inputs        
c     required for earth location and gridding by sps,  june 6, 1988)           
c                                                                               
      slat=dsin(phi)
      syaw=dsin(psi)
      si=slat**2+syaw**2                                                        
      ci=dsqrt(1.0d0-si)
      si=dsqrt(si)
      if (slat.eq.0.0d0.and.syaw .eq. 0.0d0) then
      u=0.0d0
      else
      u=datan2(slat,syaw)
      endif
      su=dsin(u)
      cu=dcos(u)
c                                                                               
c     computes longitude of the ascending node                                  
c                                                                               
      asc=lam - u
      sa=dsin(asc)
      ca=dcos(asc)
c                                                                               
c     computes the subsatellite geographic latitude                             
c                                                                               
      rlat=datan(aebe2*dtan(phi))
c                                                                               
c     computes the subsatellite longitude                                       
c                                                                               
      rlon=asc+datan2(ci*su,cu)

c                                                                               
c     computes the spacecraft to earth fixed coordinates transformation         
c     matrix:                                                                   
c         (vector in ecef coordinates) = b * (vector in s/c coordinates)        
c                                                                               
      b(1,2)=-sa*si                                                             
      b(2,2)= ca*si                                                             
      b(3,2)=-ci                                                                
      b(1,3)=-ca*cu+sa*su*ci                                                    
      b(2,3)=-sa*cu-ca*su*ci                                                    
      b(3,3)=-slat                                                              
      b(1,1)=-ca*su-sa*cu*ci                                                    
      b(2,1)=-sa*su+ca*cu*ci                                                    
      b(3,1)= cu*si                                                             
c                                                                               
c     computes the normalized spacecraft position vector in earth fixed         
c     coordinates - xs.                                                         
c                                                                               
      r=(nomorb+dr)/ae                                                          
      xs(1)=-b(1,3)*r                                                           
      xs(2)=-b(2,3)*r                                                           
      xs(3)=-b(3,3)*r                                                           
c                                                                               
c     precomputes q3 (used in lpoint)                                           
c                                                                               
      q3=xs(1)**2+xs(2)**2+aebe2*xs(3)**2-1.0d0
c                                                                               
c     computes the attitudes and misalignments if imc is off                    
c                                                                               
      if (imc.ne.0) then                                                        
c                                                                               
c     computes the solar orbit angle                                            
c                                                                               
         wa=rec(60)*ts                                                          
c                                                                               
c     computes the difference between current time, ts, and the                 
c     exponential time, rec(61). note that both times are since epoch.          
c                                                                               
         te=ts-rec(61)                                                          
c                                                                               
c     computes roll + roll misalignment                                         
c                                                                               
         roll=roll+gatt(62,rec,wa,te)                                           
c                                                                               
c     computes pitch + pitch misalignment                                       
c                                                                               
         pitch=pitch+gatt(117,rec,wa,te)                                        
c                                                                               
c     computes yaw                                                              
c                                                                               
         yaw=yaw+gatt(172,rec,wa,te)                                            
c                                                                               
c     computes roll misalignment                                                
c                                                                               
         rma=gatt(227,rec,wa,te)                                                
c                                                                               
c     computes pitch misalignment                                               
c                                                                               
         pma=gatt(282,rec,wa,te)                                                
c                                                                               
c     apply the spcecraft compensation
c                                                                               

         roll=roll+rec(15)
         pitch=pitch+rec(16)
         yaw=yaw+rec(17)

      end if                                                                    
c                                                                               
c     computes the instrument to earth fixed coordinates transformation         
c     matrix - bt                                                               
c                                                                               
      call inst2er(roll,pitch,yaw,b,bt)
      return                                                                    
      end                                                                       
c***********************************************************************        
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : gpoint                                                        
c**   source    : f.gpoint                                                      
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a     12/10/87  il   initial creation                                     
c**   a     06/10/88  il   replaced asin with atan to save time                 
c**   a     06/02/89  il   coordinate axes changed according to                 
c**                        ford's definition in sdaip, drl 504-01               
c**   4     03/08/94  sc   implemented new formulae for scan angle
c**                        correction due to misalignments
c***********************************************************************        
c**                                                                             
c**   this subroutine converts geographic latitude and longitude                
c**   to the related elevation and scan angles.                                 
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: none                                                    
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      subroutine gpoint(rlat,rlon,alf,gam,ierr)                                 
c                                                                               
c     calling parameters                                                        
c                                                                               
      real*8   rlat                                                             
c                             geographic latitude in radians (input)            
      real*8   rlon                                                             
c                             geographic longitude in radians (input)           
      real*8   alf                                                              
c                             elevation angle in radians (output)               
      real*8   gam                                                              
c                             scan angle in radians (output)                    
      integer ierr                                                              
c                             output status; 0 - successful completion,         
c                             1 - point with given lat/lon is invisible         
c                                                                               
c     local variables                                                           
c                                                                               
      real*8 f(3)                                                               
c                         pointing vector in earth centered coordinates         
      real*8 ft(3)                                                              
c                         pointing vector in instrument coordinates             
      real*8 u(3)                                                               
c                         coordinates of the earth point (km)                   
      real*8 sing,slat,w1,w2                                                    
c                                    work space                                 
c                                                                               
c     include files                                                             
c                                                                               
        include 'elcons.inc'                                                          
        include 'elcomm.inc'                                                          
c***********************************************************************        
c                                                                               
c     computes sinus of geographic (geodetic) latitude                          
c                                                                               
      sing=dsin(rlat)
      w1=aebe4*sing*sing                                                        
c                                                                               
c     sinus of the geocentric latitude                                          
c                                                                               
      slat=((0.375d0*w1-0.5d0)*w1+1.0d0)*sing/aebe2
c                                                                               
c     computes local earth radius at specified point                            
c                                                                               
      w2=slat*slat                                                              
      w1=aebe3*w2                                                               
      w1=(0.375d0*w1-0.5d0)*w1+1.d0
c                                                                               
c     computes cartesian coordinates of the point                               
c                                                                               
      u(3)=slat*w1                                                              
      w2=w1*dsqrt(1.0d0-w2)
      u(1)=w2*dcos(rlon)
      u(2)=w2*dsin(rlon)
c                                                                               
c     pointing vector from satellite to the earth point                         
c                                                                               
      f(1)=u(1)-xs(1)                                                           
      f(2)=u(2)-xs(2)                                                           
      f(3)=u(3)-xs(3)                                                           
      w2=u(1)*sngl(f(1))+u(2)*sngl(f(2))+                                       
     1   u(3)*sngl(f(3))*aebe2                                                  
c                                                                               
c     verifies visibility of the point                                          
c                                                                               
      if (w2.gt.0.0d0) then
c                               invisible point on the earth                    
                   ierr=1                                                       
                   alf=99999.0d0
                   gam=99999.0d0
                   return                                                       
       end if                                                                   
c                                                                               
c     converts pointing vector to instrument coordinates                        
c                                                                               
      ft(1)=bt(1,1)*f(1)+bt(2,1)*f(2)+bt(3,1)*f(3)                              
      ft(2)=bt(1,2)*f(1)+bt(2,2)*f(2)+bt(3,2)*f(3)                              
      ft(3)=bt(1,3)*f(1)+bt(2,3)*f(2)+bt(3,3)*f(3)                              
c                                                                               
c     converts pointing vector to scan and elevation angles and                 
c     corrects for the roll and pitch misalignments                             
c                                                                               
      gam=datan(ft(1)/sqrt(ft(2)**2+ft(3)**2))
      alf=-datan(ft(2)/ft(3))
      w1=dsin(alf)
      w2=dcos(gam)
      alf=alf+rma*(1.0d0-dcos(alf)/w2)+pma*w1*(1.0d0/w2+dtan(gam))
      gam=gam-rma*w1                                                            
      ierr=0                                                                    
      return                                                                    
      end                                                                       
c***********************************************************************        
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : inst2er                                                       
c**   source    : f.inst2er                                                     
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   1     08/16/88  il   initial creation                                     
c**   2     11/11/88  il   trigonometric functions replaced with                
c**                        small angle approximations                           
c**   3    06/02/89   il   coordinate axes changed according to                 
c**                        ford's definition in sdaip, drl 504-01               
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   inst2er accepts the single precision roll, pitch and yaw angles           
c**   of an instrument and returns the double precision instrument to           
c**   earth coordinates transformation matrix.                                  
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: none                                                    
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      subroutine inst2er(r,p,y,a,at)
c                                                                               
c     calling parameters                                                        
c                                                                               
      real*8 r                                                                  
c                                   roll angle in radians                       
      real*8 p                                                                  
c                                   pitch angle in radians                      
      real*8 y                                                                  
c                                   yaw angle in radians                        
      real*8 a(3,3)                                                             
c                                   spacecraft to ecef coordinates              
c                                   transformation matrix                       
      real*8 at(3,3)                                                            
c                                   instrument to ecef coordinates              
c                                   transformation matrix                       
c                                                                               
c     local variables                                                           
c                                                                               
      real*8    rpy(3,3)                                                        
c                                   instrument to body coordinates              
c                                   transformation matrix                       
      integer*4 i,j                                                             
c                                   indices                                     
c                                                                               
c     include files                                                             
c                                                                               
c***********************************************************************        
c                                                                               
c     we compute instrument to body coordinates transformation                  
c     matrix by using a small angle approximation of trigonometric              
c     functions of the roll, pitch and yaw.                                     
c                                                                               
      rpy(1,1)=1.0d0-0.50d0*(p*p+y*y)
      rpy(1,2)=-y                                                               
      rpy(1,3)=p                                                                
      rpy(2,1)=y+p*r                                                            
      rpy(2,2)=1.0d0-0.50d0*(y*y+r*r)
      rpy(2,3)=-r                                                               
      rpy(3,1)=-p+r*y                                                           
      rpy(3,2)=r+p*y                                                            
      rpy(3,3)=1.0d0-0.50d0*(p*p+r*r)
c                                                                               
c     multiplication of matrices a and rpy                                      
c                                                                               
      do 20 i=1,3                                                               
      do 10 j=1,3                                                               
      at(i,j)=a(i,1)*rpy(1,j)+a(i,2)*rpy(2,j)+a(i,3)*rpy(3,j)                   
   10 continue                                                                  
   20 continue                                                                  
      return                                                                    
      end                                                                       
c***********************************************************************        
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : lpoint                                                        
c**   source    : f.lpoint                                                      
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a     01/09/89  il   initial creation                                     
c**   a     06/02/89  il   coordinate axes changed according to                 
c**                        ford's definition in sdaip, drl504-01                
c**   3     03/08/94  sc   implemented new formulae for scan angle
c**                        corrections due to misalignments
c***********************************************************************        
c**                                                                             
c**   this subroutine converts the instrument elevation and scan                
c**   angles to the related geographic latitude and longitude.                  
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: none                                                    
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      subroutine lpoint(alpha,zeta,rlat,rlon,ierr)                              
c                                                                               
c     calling parameters                                                        
c                                                                               
      real*8   alpha                                                            
c                             elevation angle (rad)                             
      real*8   zeta                                                             
c                             scan angle (rad)                                  
      real*8   rlat                                                             
c                             latitude in radians (output)                      
      real*8   rlon                                                             
c                             longitude in radians (output)                     
      integer ierr                                                              
c                             output status; 0 - point on the earth             
c                             found, 1 - instrument points off earth            
c                                                                               
c     local variables                                                           
c                                                                               
      real*8 g1(3)                                                              
c                          pointing vector in earth centered coordinates        
      real*8 h                                                                  
c                          slant distance to the earth point (km)               
      real*8 q1,q2,d                                                            
c                          work space                                           
      real*8 g(3)                                                               
c                          pointing vector in instrument coordinates            
      real*8 u(3)                                                               
c                          coordinates of the earth point (km)                  
      real*8 sa,ca,da,dz,d1,cz                                                  
c                                     work space                                
c                                                                               
c     include files                                                             
c                                                                               
        include 'elcons.inc'                                                          
        include 'elcomm.inc'                                                          
c***********************************************************************        
      ierr=1                                                                    
c                                                                               
c     computes trigonometric functions of the scan and elevation                
c     angles corrected for the roll and pitch misalignments                     
c                                                                               
      ca=dcos(alpha)
      sa=dsin(alpha)
      cz=dcos(zeta)
      da=alpha-pma*sa*(1.0d0/cz+dtan(zeta))-rma*(1.0d0-ca/cz)
      dz=zeta+rma*sa                                                            
c                              corrected scan angle                             
      cz=dcos(dz)
c                                                                               
c     computes pointing vector in instrument coordinates                        
c                                                                               
      g(1)=dsin(dz)
      g(2)=-cz*dsin(da)
      g(3)=cz*dcos(da)
c                                                                               
c     transforms the pointing vector to earth fixed coordinates                 
c                                                                               
      g1(1)=bt(1,1)*g(1)+bt(1,2)*g(2)+bt(1,3)*g(3)                              
      g1(2)=bt(2,1)*g(1)+bt(2,2)*g(2)+bt(2,3)*g(3)                              
      g1(3)=bt(3,1)*g(1)+bt(3,2)*g(2)+bt(3,3)*g(3)                              
c                                                                               
c     computes coefficients and solves a quadratic equation to                  
c     find the intersect of the pointing vector with the earth                  
c     surface                                                                   
c                                                                               
      q1=g1(1)**2+g1(2)**2+aebe2*g1(3)**2                                       
      q2=xs(1)*g1(1)+xs(2)*g1(2)+aebe2*xs(3)*g1(3)                              
      d=q2*q2-q1*q3                                                             
      if (dabs(d).lt.1.d-9) d=0.0d0
c                                                                               
c     if the disciminante of the equation, d, is negative, the                  
c     instrument points off the earth                                           
c                                                                               
      if (d.lt.0.0d0) then
         rlat=999999.0d0
         rlon=999999.0d0
         return                                                                 
      end if                                                                    
      d=dsqrt(d)                                                                
c                                                                               
c     slant distance from the satellite to the earth point                      
c                                                                               
      h=-(q2+d)/q1                                                              
c                                                                               
c     cartesian coordinates of the earth point                                  
c                                                                               
      u(1)=xs(1)+h*g1(1)                                                        
      u(2)=xs(2)+h*g1(2)                                                        
      u(3)=xs(3)+h*g1(3)                                                        
c                                                                               
c     sinus of geocentric latitude                                              
c                                                                               
      d1=u(3)/dsqrt(u(1)**2+u(2)**2+u(3)**2)
c                                                                               
c     geographic (geodetic) coordinates of the point                            
c                                                                               
      rlat=datan(aebe2*d1/dsqrt(1.0d0-d1*d1))
      rlon=datan2(u(2),u(1))
      ierr=0                                                                    
      return                                                                    
      end                                                                       
c***********************************************************************        
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : sndeloc                                                       
c**   source    : f.sndeloc                                                     
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a     02/16/89  il   initial creation                                     
c***********************************************************************        
c**                                                                             
c**   sndeloc accepts the mirror position in cycles and increments,             
c**   servo error values, and the positional offsets for four detectors         
c**   of a selected sounder channel and computes the detector earth             
c**   locations in latitude/longitude coordinates.                              
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: none                                                    
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : lpoint                                                  
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      subroutine sndelo(cyew,incew,cyns,incns,svew,svns,doff,geo)               
c                                                                               
c     calling parameters                                                        
c                                                                               
      integer cyew                                                              
c                         e-w cycles                                            
      integer incew                                                             
c                         e-w increments                                        
      integer cyns                                                              
c                         n-s cycles                                            
      integer incns                                                             
c                         n-s increments                                        
      real*8 svew                                                               
c                         e-w servo error in radians                            
      real*8 svns                                                               
c                         n-s servo error in radians                            
      real*8 doff(4,2)                                                          
c                         offsets for 4 detectors (radians)                     
c                            doff(*,1) = e-w offset                             
c                            doff(*,2) = n-s offset                             
      real*8 geo(4,2)                                                           
c                         geographic coordinates related to 4 detectors         
c                            geo(*,1) = latitude in radians                     
c                            geo(*,2) = longitude in radians                    
c                                                                               
c     local variables                                                           
c                                                                               
c     real*8 e,s,h,ev,sc,alpha,beta,sine,cose,de,ds
      real*8 e,s,h,ev,sc,sine,cose,de,ds
      integer i,ier                                                             
c                                                                               
c     include files                                                             
c                                                                               
        include 'instco.inc'                                                          
c***********************************************************************        
c                                                                               
c     convert the mirror position, given in cycles and increments, to           
c     elevation and scan angles                                                 
c                                                                               
c     e=(cyns*incmax(2)+incns)*elvinc(2)-elvmax(2)
      e=((cyns-9)*incmax(2)+incns)*elvinc(2)-elvmax(2)                              
      s=(cyew*incmax(2)+incew)*scninc(2)-scnmax(2)                              
c                                                                               
c     correct elevation and scan angles for servo errors obtaining the          
c     true mirror pointing                                                      
c                                                                               
      e=e+svns                                                                  
      s=s+svew
      sine=dsin(e)
      cose=dcos(e)
      h=-2.0d0*scnpx(2)
c
c     compute detector rotation offsets for each detector
c
c      alpha = 0.643501d0 + e
c      beta = 0.244979d0 - e
c
c      doff(1,1) = -0.064976d0
c      doff(1,2) = 0.00042d0
c      doff(2,1) = 0.00056d0
c      doff(2,2) = 0.00014d0 
c      doff(3,1) = -0.064976d0
c      doff(3,2) = -0.065396d0
c      doff(4,1) = 0.00056d0
c      doff(4,2) = -0.065116d0
c      doff(1,1) = - 700.0d0*dcos(alpha)*1.0d-6
c      doff(1,2) =   700.0d0*dsin(alpha)*1.0d-6
c      doff(2,1) =   577.23479d0*dcos(beta)*1.d-6
c      doff(2,2) =   577.23479d0*dsin(beta)*1.0d-6
c      doff(3,1) = - 577.23479d0*dcos(beta)*1.0d-6
c      doff(3,2) = - 577.23479d0*dsin(beta)*1.0d-6
c      doff(4,1) =   700.0d0*dcos(alpha)*1.0d-6
c      doff(4,2) = - 700.0d0*dsin(alpha)*1.0d-6
c                                                                               
c     compute earth locations for four detectors
c                                                                               
      do 10 i=1,4

c     compute positional offsets of i-th detector

        de=(2.5-i)*elvln(2)+doff(i,2)
        ds=h+doff(i,1)
c                                                                               
c     compute elevation and scan angles related to i-th detector                
c     and correct them for the detector positional offsets                      
c                                                                               
c           ev=e+doff(i,2)
c           sc=s+doff(i,1)
c
c     convert positional offsets to angular offsets and
c     correct elevation and scan angles

            ev = e + de*cose - ds*sine
            sc = s + de*sine + ds*cose

c     transform detector's pointing angles to geographic coordinates            
c     of the corresponding point on the earth surface.                          
c     note:  if a detector looks off the earth, the related latitude            
c            and longitude are set to 999999.                                   
c                                                                               
           call lpoint(ev,sc,geo(i,1),geo(i,2),ier)                             
           h=-h                                                                 
   10 continue                                                                  
      return                                                                    
      end                                                                       
c***********************************************************************        
c***********************************************************************        
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : evsc2lpf                                                      
c**   source    : f.evsc2lpf                                                    
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a     10/27/88  il   initial creation                                     
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   this subroutine converts elevation and scan angles                        
c**   to the fractional line and pixel numbers.                                 
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: none                                                    
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      subroutine evsc2l(instr,elev,scan,rl,rp)                                  
c                                                                               
c     calling parameters                                                        
c                                                                               
      integer instr                                                             
c                       instrument code (1-imager, 2-sounder)                   
      real*8 elev                                                                 
c                       elevation angle in radians                              
      real*8 scan                                                                 
c                       scan angle in radians                                   
      real*8 rl
c                       line number                                             
      real*8 rp
c                       pixel number                                            
c                                                                               
c     local variables - none                                                    
c                                                                               
c                                                                               
c     include files                                                             
c                                                                               
        include 'instco.inc'                                                          
c**************************************************************                 
c                                                                               
c     compute fractional line number                                            
c                                                                               
      rl=(elvmax(instr)-elev)/elvln(instr)                                      
      if (instr.eq.1) then
           rl=rl+4.5d0
      else
           rl=rl+2.5d0
      end if
c                                                                               
c     compute fractional pixel number                                           
c                                                                               
      rp=(scnmax(instr)+scan)/scnpx(instr)+1.0d0
      return                                                                    
      end                                                                       
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : evln                                                          
c**   source    : f.evln                                                        
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a     10/27/88  il   initial creation                                     
c***********************************************************************        
c**                                                                             
c**   this function converts fractional line number to elevation angle          
c**   in radians.                                                               
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: none                                                    
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      function evln(instr,rline)                                                
c                                                                               
c     calling parameters                                                        
c                                                                               
      integer instr                                                             
c                       instrument code (1-imager, 2-sounder)                   
      real*8  evln,rline
c                       fractional line  number                                 
c                                                                               
c     local variables - none                                                    
c                                                                               
c                                                                               
c     include files                                                             
c                                                                               
        include 'instco.inc'                                                          
c***********************************************************************        
      if (instr.eq.1) then                                                      
       evln=elvmax(instr)*1.0d0 - (rline-4.5 )*(elvln(instr)*1.0d0)
      else                                                                      
           evln=elvmax(instr)*1.0d0 - (rline
     +          -2.5d0)*(elvln(instr)*1.0d0)
      end if                                                                    
      return                                                                    
      end                                                                       
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : scpx                                                          
c**   source    : f.scpx                                                        
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a     09/22/87  il   initial creation                                     
c***********************************************************************        
c**                                                                             
c**   this function converts fractional pixel number to scan angle              
c**   in radians.                                                               
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : any                                                     
c**   commons modified: none                                                    
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      function scpx(instr,pix)                                                  
c                                                                               
c     calling parameters                                                        
c                                                                               
      integer instr                                                             
c                       instrument code (1-imager, 2-sounder)                   
      real*8 scpx,pix
c                       fractional pixel number                                 
c                                                                               
c     local variables                                                           
c                                                                               
c                                                                               
c     include files                                                             
c                                                                               
        include 'instco.inc'                                                    
c***********************************************************************        
      scpx=((pix*1.0d0)-1.0d0)*(scnpx(instr)*1.0d0)
     +     -(scnmax(instr)*1.0d0)
      return                                                                    
      end                                                                       
c***********************************************************************        
c**                                                                             
c**   integral systems, inc.                                                    
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   project   : operations ground equipment for goes-next                     
c**   system    : earth location users guide                                    
c**   routine   : gatt                                                          
c**   source    : f.gatt                                                        
c**   load name : any                                                           
c**   programmer: igor levine                                                   
c**                                                                             
c**   ver.    data    by   comment                                              
c**   ----  --------  ---  ---------------------------------------------        
c**   a     12/01/88  il   initial creation                                     
c**                                                                             
c***********************************************************************        
c**                                                                             
c**    this function computes an attitude/misalignment angle from a             
c**    given subset of the o&a parameters in gvar blok 0.                       
c**    argument k0 indicates the first word of the subset.                      
c**                                                                             
c***********************************************************************        
c**                                                                             
c**   called by       : lmodel                                                  
c**   commons modified: none
c**   inputs          : none                                                    
c**   outputs         : none                                                    
c**   routines called : none                                                    
c**                                                                             
c***********************************************************************        
c***********************************************************************        
      function gatt(k0,rec,wa,te)
c                                                                               
c     calling parameters                                                        
c                                                                               
      integer k0                                                                
c                    starting position of a parameter subset in the             
c                    o&a set                                                    
      real*8 rec(336)                                                           
c                    input o&a parameter set                                    
      real*8 wa
c                    input solar orbit angle in radians                         
      real*8 te
c                    input exponential time delay from epoch (minutes)          
c                                                                               
c     local variables                                                           
c                                                                               
      integer*4 i,j,m,l,ll,k                                                    
      real*8 gatt,ir,jr,mr,att
      equivalence (j,jr),(m,mr)
c                                                                               
c     include files                                                             
c                                                                               
c***********************************************************************        
c                                                                               
c     constant component                                                        
c                                                                               
      k=k0                                                                      
      att=rec(k+2)                                                              
c                                                                               
c     computes the exponential term                                             
c                                                                               
      if ((te.ge.0.0d0).and.(rec(k+1).gt. 0))then
         att=att+rec(k)*dexp(-te/rec(k+1))
      endif
c                                                                               
c     extracts the number of sinusoids                                          
c                                                                               
      ir=rec(k+3)
      i = int(ir)
c                                                                               
c     calculation of sinusoids                                                  
c                                                                               
      do 10 l=1,i                                                               
           att=att+rec(k+2*l+2)*dcos(wa*l+rec(k+2*l+3))
   10 continue                                                                  
c                                                                               
c     pointer to the number of monomial sinusoids                               
c                                                                               
      k=k+34                                                                    
c                                                                               
c     extacts number of monomial sinusoids                                      
c                                                                               
      ir=rec(k)
      i=int(ir)
c     kkk=rec(k)
c                                                                               
c     computes monomial sinusoids
c                                                                               
      do 20 l=1,i
           ll=k+5*l
c                                                                               
c          order of sinusoid                                                    
c                                                                               
           jr=rec(ll-4)                                                         
c                                                                               
c          order of monomial sinusoid                                           
c                                                                               
           mr=rec(ll-3)                                                         
c                                                                               
           att=att+rec(ll-2)*((wa-rec(ll))**m)*dcos(j*wa+rec(ll-1))
   20 continue                                                                  
      gatt=att                                                                  
      return                                                                    
      end                                                                       
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *        
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *        
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *        
c dlat,dlon...... 50.00000  -150.00000       (15x,f9.5,2x,f10.5)
c line...(1),(2). 3584   1230                (15x,i5,2x,i5)
c ipixel.(1),(2).10253   1155                (15x,i5,2x,i5)
c key............0                           (15x,i1)
c imc............1                           (15x,i1)
c instr..........2                           (15x,i1)
c epoch time.....1989032062934.567           (15x,i4,i3,2i2,f6.3)
c start time.....1989032064934.567           (15x,i4,i3,2i2,f6.3)
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* *
c      dlat ---------- latitude  (negative is south)(f9.5)
c      dlon ---------- longitude (negative is west)(f10.5)

c      line(2) ------- line no.  (see instr for line(instr)(i5)
c      ipixel(2) ----- pixel no. (see instr for ipixel(instr)(i5)

c      key ----------- 0 = line/pixel to lat/lon conversion.
c                      1 = lat/lon to line/pixel conversion.

c      imc ----------- imc status (0=imc on, 1=imc off)

c      instr --------- instrument (1=imager, 2=sounder)

c      epoch time----- year, j-day, hour, minute, seconds

c      start time----- year, j-day, hour, minute, seconds
c  test for user:input data block(unnum data block)..........
c  fifty(50) lines available for use in input block..........
c***********************************************************************
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

