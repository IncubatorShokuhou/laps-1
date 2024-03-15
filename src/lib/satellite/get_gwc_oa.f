       subroutine get_gwc_oa(filename,imc_id,rec,ir,istatus)
c
c ***************************************************************************** 
c *****************************************************************************
c	subroutine to read in afgwc's sdhs goes "orbit and attitude" (oa) data,
c	convert from gould floating point values to sun/ibm for ultimate use in 
c	noaa's "gimloc" earth-location software.
c
c	original:	bruce h. thomas, march 1997
c			the aerospace corporation and dmsp spo
c			environmental applications center
c
c	purpose:	the code was developed to allow sdhs goes-gvar data to
c			be used within noaa/fsl's laps model.  the project was
c			undertaken as a prototype for the gtwaps program.
c
c	notes:		code is explicit when decoding each "o+a" word for
c			clarity and documentation; could be simplified with
c			loops if shorter code is desired!!!!
c
c	subroutines called: 
c 
c	1) bcd_to_int - decodes the "binary coded decimal" date/time from o+a
c	2) byteswp4   - swaps "little endian" 4-byte integers into "big endian"
c	3) cnvtglfl   - decodes "gould floating point" values from o+a
c	4) cnvtint    - decodes "gould integer" values from o+a
c
c *****************************************************************************
c *****************************************************************************
c
      implicit none
c
      logical print_flag
        
      character*4 imc_id
      integer oadata(336),i,ir
      integer byteswp4
      real*8 rec(ir)
      real*8 rtime
c
      character*255 filename
      integer     istatus
c
c     data print_flag / .true. /
      data print_flag / .false. /
c     data rec / 336*0.0d0 /
c
c     equivalence (oadata(1), imc_id)
      equivalence (oadata(12), rtime)
c
      print *,'***************************'
      print *,'goes_oa: start of execution'
      print *,'***************************'
      istatus = 0
c
c j.smart: 5-13-97. initialize array rec
c
      do i=1,ir
         rec(i)=0.0d0
      enddo
c
c *****************************************************************************
c     open input file through soft link to "fort.10"
c *****************************************************************************
c
      open(unit=10,file=filename,form='unformatted',access='direct'
     &,status='old',recl=1412,err=1010)
c
c *****************************************************************************
c     read in the orbit and attitude data
c *****************************************************************************
c
      read(10,rec=1,err=300) oadata
c
c *****************************************************************************
c     close input file
c *****************************************************************************
c
      close(10)
c
c *****************************************************************************
c *****************************************************************************
c     swap input data from "little endian (vax)" to "big endian (sun/ibm), aka
c     change the 4-byte integer's byte order.
c *****************************************************************************
c *****************************************************************************
c
      do 100 i=2,336
        oadata(i)=byteswp4(oadata(i))
100   continue
c
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
c     convert the "gould floating point" into sun/ibm double precision using
c     the subroutine "cnvtglfl."
c *****************************************************************************
c *****************************************************************************
c *****************************************************************************
c
c *****************************************************************************
c     word 1 is the four byte character "imc identifier"
c     words 2-4 are spares
c *****************************************************************************
      write(imc_id,112)oadata(1)
112   format(a)

      call cnvtglfl(oadata(005),rec(005))	! reference longitude (+ east)
      call cnvtglfl(oadata(006),rec(006))	! reference radial distance from nominal
      call cnvtglfl(oadata(007),rec(007))	! reference latitude (+ north)
      call cnvtglfl(oadata(008),rec(008))	! reference orbit yaw
      call cnvtglfl(oadata(009),rec(009))	! reference attitude: role
      call cnvtglfl(oadata(010),rec(010))	! reference attitude: pitch
      call cnvtglfl(oadata(011),rec(011))	! reference attitude: yaw
c *****************************************************************************
c     input words 12 and 13 represent the "binary-coded-decimal" epoch time
c     - convert the "binary-coded-decimal" (bcd) into time code required
c *****************************************************************************
      call bcdtoint (oadata(12), oadata(13), rtime)
      rec(12) = rtime
c *****************************************************************************
c     spacecraft compensation roll, pitch, and yaw
c     change in longitude from reference (13 values)
c     change in radial distance (11 values)
c     sine in geocentric latitude (9 values)
c     sine orbit yaw (9 values)
c     daily solar rate
c     exponential start time from epoch
c *****************************************************************************
      call cnvtglfl(oadata(014),rec(014))	! imc set enable time from epoch
      call cnvtglfl(oadata(015),rec(015))	! spacecraft compensation roll
      call cnvtglfl(oadata(016),rec(016))	! spacecraft compensation pitch
      call cnvtglfl(oadata(017),rec(017))	! spacecraft compensation yaw
      call cnvtglfl(oadata(018),rec(018))	! change in long from ref #01
      call cnvtglfl(oadata(019),rec(019))	! change in long from ref #02
      call cnvtglfl(oadata(020),rec(020))	! change in long from ref #03
      call cnvtglfl(oadata(021),rec(021))	! change in long from ref #04
      call cnvtglfl(oadata(022),rec(022))	! change in long from ref #05
      call cnvtglfl(oadata(023),rec(023))	! change in long from ref #06
      call cnvtglfl(oadata(024),rec(024))	! change in long from ref #07
      call cnvtglfl(oadata(025),rec(025))	! change in long from ref #08
      call cnvtglfl(oadata(026),rec(026))	! change in long from ref #09
      call cnvtglfl(oadata(027),rec(027))	! change in long from ref #10
      call cnvtglfl(oadata(028),rec(028))	! change in long from ref #11
      call cnvtglfl(oadata(029),rec(029))	! change in long from ref #12
      call cnvtglfl(oadata(030),rec(030))	! change in long from ref #13
      call cnvtglfl(oadata(031),rec(031))	! change in radial dist from ref #01
      call cnvtglfl(oadata(032),rec(032))	! change in radial dist from ref #02
      call cnvtglfl(oadata(033),rec(033))	! change in radial dist from ref #03
      call cnvtglfl(oadata(034),rec(034))	! change in radial dist from ref #04
      call cnvtglfl(oadata(035),rec(035))	! change in radial dist from ref #05
      call cnvtglfl(oadata(036),rec(036))	! change in radial dist from ref #06
      call cnvtglfl(oadata(037),rec(037))	! change in radial dist from ref #07
      call cnvtglfl(oadata(038),rec(038))	! change in radial dist from ref #08
      call cnvtglfl(oadata(039),rec(039))	! change in radial dist from ref #09
      call cnvtglfl(oadata(040),rec(040))	! change in radial dist from ref #10
      call cnvtglfl(oadata(041),rec(041))	! change in radial dist from ref #11
      call cnvtglfl(oadata(042),rec(042))	! sine geocentric lat, total #01
      call cnvtglfl(oadata(043),rec(043))	! sine geocentric lat, total #02
      call cnvtglfl(oadata(044),rec(044))	! sine geocentric lat, total #03
      call cnvtglfl(oadata(045),rec(045))	! sine geocentric lat, total #04
      call cnvtglfl(oadata(046),rec(046))	! sine geocentric lat, total #05
      call cnvtglfl(oadata(047),rec(047))	! sine geocentric lat, total #06
      call cnvtglfl(oadata(048),rec(048))	! sine geocentric lat, total #07
      call cnvtglfl(oadata(049),rec(049))	! sine geocentric lat, total #08
      call cnvtglfl(oadata(050),rec(050))	! sine geocentric lat, total #09
      call cnvtglfl(oadata(051),rec(051))	! sine orbit yaw, total #01
      call cnvtglfl(oadata(052),rec(052))	! sine orbit yaw, total #02
      call cnvtglfl(oadata(053),rec(053))	! sine orbit yaw, total #03
      call cnvtglfl(oadata(054),rec(054))	! sine orbit yaw, total #04
      call cnvtglfl(oadata(055),rec(055))	! sine orbit yaw, total #05
      call cnvtglfl(oadata(056),rec(056))	! sine orbit yaw, total #06
      call cnvtglfl(oadata(057),rec(057))	! sine orbit yaw, total #07
      call cnvtglfl(oadata(058),rec(058))	! sine orbit yaw, total #08
      call cnvtglfl(oadata(059),rec(059))	! sine orbit yaw, total #09
      call cnvtglfl(oadata(060),rec(060))	! daily solar rate
      call cnvtglfl(oadata(061),rec(061))	! exponential start time from epoch
c *****************************************************************************
c *****************************************************************************
c     block of code applys to roll attitude angle:
c *****************************************************************************
c *****************************************************************************
      call cnvtglfl(oadata(062),rec(062))	! exponential magnitude 
      call cnvtglfl(oadata(063),rec(063))	! exponential time constant
      call cnvtglfl(oadata(064),rec(064))	! constant, mean attitude angle
      call cnvtint(oadata(065))
      rec(065) = dfloat(oadata(065))		! number of sinusoids/angles
c *****************************************************************************
c     magnitude/phase angle of first->fifteenth order sinusoids 
c *****************************************************************************
      call cnvtglfl(oadata(066),rec(066))	! magnitude of 1st order sinusoid
      call cnvtglfl(oadata(067),rec(067))	! phase angle of 1st order sinusoid
      call cnvtglfl(oadata(068),rec(068))	! magnitude of 2nd order sinusoid
      call cnvtglfl(oadata(069),rec(069))	! phase angle of 2nd order sinusoid
      call cnvtglfl(oadata(070),rec(070))	! magnitude of 3rd order sinusoid
      call cnvtglfl(oadata(071),rec(071))	! phase angle of 3rd order sinusoid
      call cnvtglfl(oadata(072),rec(072))	! magnitude of 4th order sinusoid
      call cnvtglfl(oadata(073),rec(073))	! phase angle of 4th order sinusoid
      call cnvtglfl(oadata(074),rec(074))	! magnitude of 5th order sinusoid
      call cnvtglfl(oadata(075),rec(075))	! phase angle of 5th order sinusoid
      call cnvtglfl(oadata(076),rec(076))	! magnitude of 6th order sinusoid
      call cnvtglfl(oadata(077),rec(077))	! phase angle of 6th order sinusoid
      call cnvtglfl(oadata(078),rec(078))	! magnitude of 7th order sinusoid
      call cnvtglfl(oadata(079),rec(079))	! phase angle of 7th order sinusoid
      call cnvtglfl(oadata(080),rec(080))	! magnitude of 8th order sinusoid
      call cnvtglfl(oadata(081),rec(081))	! phase angle of 8th order sinusoid
      call cnvtglfl(oadata(082),rec(082))	! magnitude of 9th order sinusoid
      call cnvtglfl(oadata(083),rec(083))	! phase angle of 9th order sinusoid
      call cnvtglfl(oadata(084),rec(084))	! magnitude of 10th order sinusoid
      call cnvtglfl(oadata(085),rec(085))	! phase angle of 10th order sinusoid
      call cnvtglfl(oadata(086),rec(086))	! magnitude of 11th order sinusoid
      call cnvtglfl(oadata(087),rec(087))	! phase angle of 11th order sinusoid
      call cnvtglfl(oadata(088),rec(088))	! magnitude of 12th order sinusoid
      call cnvtglfl(oadata(089),rec(089))	! phase angle of 12th order sinusoid
      call cnvtglfl(oadata(090),rec(090))	! magnitude of 13th order sinusoid
      call cnvtglfl(oadata(091),rec(091))	! phase angle of 13th order sinusoid
      call cnvtglfl(oadata(092),rec(092))	! magnitude of 14th order sinusoid
      call cnvtglfl(oadata(093),rec(093))	! phase angle of 14th order sinusoid
      call cnvtglfl(oadata(094),rec(094))	! magnitude of 15th order sinusoid
      call cnvtglfl(oadata(095),rec(095))	! phase angle of 15th order sinusoid
      call cnvtint(oadata(096))
      rec(096) = dfloat(oadata(096))		! number of monomial sinusoids
      call cnvtint(oadata(097))
      rec(097) = dfloat(oadata(097))		! order of applicable sinusoid
      call cnvtint(oadata(098))
      rec(098) = dfloat(oadata(098))		! order of 1st monomial sinusoid
      call cnvtglfl(oadata(099),rec(099))	! magnitude of 1st monomial sinusoid
      call cnvtglfl(oadata(100),rec(100))	! phase angle of 1st monomial sinusoid
      call cnvtglfl(oadata(101),rec(101))	! angle from epoch where 1st monomial is zero
      call cnvtint(oadata(102))
      rec(102) = dfloat(oadata(102))		! order of applicable sinusoid
      call cnvtint(oadata(103))
      rec(103) = dfloat(oadata(103))		! order of 2nd monomial sinusoid
      call cnvtglfl(oadata(104),rec(104))	! magnitude of 2nd monomial sinusoid
      call cnvtglfl(oadata(105),rec(105))	! phase angle of 2nd monomial sinusoid
      call cnvtglfl(oadata(106),rec(106))	! angle from epoch where 2nd monomial is zero
      call cnvtint(oadata(107))
      rec(107) = dfloat(oadata(107))		! order of applicable sinusoid
      call cnvtint(oadata(108))
      rec(108) = dfloat(oadata(108))		! order of 3rd monomial sinusoid
      call cnvtglfl(oadata(109),rec(109))	! magnitude of 3rd monomial sinusoid
      call cnvtglfl(oadata(110),rec(110))	! phase angle of 3rd monomial sinusoid
      call cnvtglfl(oadata(111),rec(111))	! angle from epoch where 3rd monomial is zero
      call cnvtint(oadata(112))
      rec(112) = dfloat(oadata(112))		! order of applicable sinusoid
      call cnvtint(oadata(113))
      rec(113) = dfloat(oadata(113))		! order of 4th monomial sinusoid
      call cnvtglfl(oadata(114),rec(114))	! magnitude of 4th monomial sinusoid
      call cnvtglfl(oadata(115),rec(115))	! phase angle of 4th monomial sinusoid
      call cnvtglfl(oadata(116),rec(116))	! angle from epoch where 4th monomial is zero
c *****************************************************************************
c *****************************************************************************
c     block of code applys to pitch attitude angle:
c *****************************************************************************
c *****************************************************************************
      call cnvtglfl(oadata(117),rec(117))	! exponential magnitude 
      call cnvtglfl(oadata(118),rec(118))	! exponential time constant
      call cnvtglfl(oadata(119),rec(119))	! constant, mean attitude angle
      call cnvtint(oadata(120))
      rec(120) = dfloat(oadata(120))		! number of sinusoids/angles
c *****************************************************************************
c     magnitude/phase angle of first->fifteenth order sinusoids 
c *****************************************************************************
      call cnvtglfl(oadata(121),rec(121))	! magnitude of 1st order sinusoid
      call cnvtglfl(oadata(122),rec(122))	! phase angle of 1st order sinusoid
      call cnvtglfl(oadata(123),rec(123))	! magnitude of 2nd order sinusoid
      call cnvtglfl(oadata(124),rec(124))	! phase angle of 2nd order sinusoid
      call cnvtglfl(oadata(125),rec(125))	! magnitude of 3rd order sinusoid
      call cnvtglfl(oadata(126),rec(126))	! phase angle of 3rd order sinusoid
      call cnvtglfl(oadata(127),rec(127))	! magnitude of 4th order sinusoid
      call cnvtglfl(oadata(128),rec(128))	! phase angle of 4th order sinusoid
      call cnvtglfl(oadata(129),rec(129))	! magnitude of 5th order sinusoid
      call cnvtglfl(oadata(130),rec(130))	! phase angle of 5th order sinusoid
      call cnvtglfl(oadata(131),rec(131))	! magnitude of 6th order sinusoid
      call cnvtglfl(oadata(132),rec(132))	! phase angle of 6th order sinusoid
      call cnvtglfl(oadata(133),rec(133))	! magnitude of 7th order sinusoid
      call cnvtglfl(oadata(134),rec(134))	! phase angle of 7th order sinusoid
      call cnvtglfl(oadata(135),rec(135))	! magnitude of 8th order sinusoid
      call cnvtglfl(oadata(136),rec(136))	! phase angle of 8th order sinusoid
      call cnvtglfl(oadata(137),rec(137))	! magnitude of 9th order sinusoid
      call cnvtglfl(oadata(138),rec(138))	! phase angle of 9th order sinusoid
      call cnvtglfl(oadata(139),rec(139))	! magnitude of 10th order sinusoid
      call cnvtglfl(oadata(140),rec(140))	! phase angle of 10th order sinusoid
      call cnvtglfl(oadata(141),rec(141))	! magnitude of 11th order sinusoid
      call cnvtglfl(oadata(142),rec(142))	! phase angle of 11th order sinusoid
      call cnvtglfl(oadata(143),rec(143))	! magnitude of 12th order sinusoid
      call cnvtglfl(oadata(144),rec(144))	! phase angle of 12th order sinusoid
      call cnvtglfl(oadata(145),rec(145))	! magnitude of 13th order sinusoid
      call cnvtglfl(oadata(146),rec(146))	! phase angle of 13th order sinusoid
      call cnvtglfl(oadata(147),rec(147))	! magnitude of 14th order sinusoid
      call cnvtglfl(oadata(148),rec(148))	! phase angle of 14th order sinusoid
      call cnvtglfl(oadata(149),rec(149))	! magnitude of 15th order sinusoid
      call cnvtglfl(oadata(150),rec(150))	! phase angle of 15th order sinusoid
      call cnvtint(oadata(151))
      rec(151) = dfloat(oadata(151))		! number of monomial sinusoids
      call cnvtint(oadata(152))
      rec(152) = dfloat(oadata(152))		! order of applicable sinusoid
      call cnvtint(oadata(153))
      rec(153) = dfloat(oadata(153))		! order of 1st monomial sinusoid
      call cnvtglfl(oadata(154),rec(154))	! magnitude of 1st monomial sinusoid
      call cnvtglfl(oadata(155),rec(155))	! phase angle of 1st monomial sinusoid
      call cnvtglfl(oadata(156),rec(156))	! angle from epoch where 1st monomial is zero
      call cnvtint(oadata(157))
      rec(157) = dfloat(oadata(157))		! order of applicable sinusoid
      call cnvtint(oadata(158))
      rec(158) = dfloat(oadata(158))		! order of 2nd monomial sinusoid
      call cnvtglfl(oadata(159),rec(159))	! magnitude of 2nd monomial sinusoid
      call cnvtglfl(oadata(160),rec(160))	! phase angle of 2nd monomial sinusoid
      call cnvtglfl(oadata(161),rec(161))	! angle from epoch where 2nd monomial is zero
      call cnvtint(oadata(162))
      rec(162) = dfloat(oadata(162))		! order of applicable sinusoid
      call cnvtint(oadata(163))
      rec(163) = dfloat(oadata(163))		! order of 3rd monomial sinusoid
      call cnvtglfl(oadata(164),rec(164))	! magnitude of 3rd monomial sinusoid
      call cnvtglfl(oadata(165),rec(165))	! phase angle of 3rd monomial sinusoid
      call cnvtglfl(oadata(166),rec(166))	! angle from epoch where 3rd monomial is zero
      call cnvtint(oadata(167))
      rec(167) = dfloat(oadata(167))		! order of applicable sinusoid
      call cnvtint(oadata(168))
      rec(168) = dfloat(oadata(168))		! order of 4th monomial sinusoid
      call cnvtglfl(oadata(169),rec(169))	! magnitude of 4th monomial sinusoid
      call cnvtglfl(oadata(170),rec(170))	! phase angle of 4th monomial sinusoid
      call cnvtglfl(oadata(171),rec(171))	! angle from epoch where 4th monomial is zero
c *****************************************************************************
c *****************************************************************************
c     block of code applys to yaw attitude angle:
c *****************************************************************************
c *****************************************************************************
      call cnvtglfl(oadata(172),rec(172))	! exponential magnitude 
      call cnvtglfl(oadata(173),rec(173))	! exponential time constant
      call cnvtglfl(oadata(174),rec(174))	! constant, mean attitude angle
      call cnvtint(oadata(175))
      rec(175) = dfloat(oadata(175))		! number of sinusoids/angles
c *****************************************************************************
c     magnitude/phase angle of first->fifteenth order sinusoids 
c *****************************************************************************
      call cnvtglfl(oadata(176),rec(176))	! magnitude of 1st order sinusoid
      call cnvtglfl(oadata(177),rec(177))	! phase angle of 1st order sinusoid
      call cnvtglfl(oadata(178),rec(178))	! magnitude of 2nd order sinusoid
      call cnvtglfl(oadata(179),rec(179))	! phase angle of 2nd order sinusoid
      call cnvtglfl(oadata(180),rec(180))	! magnitude of 3rd order sinusoid
      call cnvtglfl(oadata(181),rec(181))	! phase angle of 3rd order sinusoid
      call cnvtglfl(oadata(182),rec(182))	! magnitude of 4th order sinusoid
      call cnvtglfl(oadata(183),rec(183))	! phase angle of 4th order sinusoid
      call cnvtglfl(oadata(184),rec(184))	! magnitude of 5th order sinusoid
      call cnvtglfl(oadata(185),rec(185))	! phase angle of 5th order sinusoid
      call cnvtglfl(oadata(186),rec(186))	! magnitude of 6th order sinusoid
      call cnvtglfl(oadata(187),rec(187))	! phase angle of 6th order sinusoid
      call cnvtglfl(oadata(188),rec(188))	! magnitude of 7th order sinusoid
      call cnvtglfl(oadata(189),rec(189))	! phase angle of 7th order sinusoid
      call cnvtglfl(oadata(190),rec(190))	! magnitude of 8th order sinusoid
      call cnvtglfl(oadata(191),rec(191))	! phase angle of 8th order sinusoid
      call cnvtglfl(oadata(192),rec(192))	! magnitude of 9th order sinusoid
      call cnvtglfl(oadata(193),rec(193))	! phase angle of 9th order sinusoid
      call cnvtglfl(oadata(194),rec(194))	! magnitude of 10th order sinusoid
      call cnvtglfl(oadata(195),rec(195))	! phase angle of 10th order sinusoid
      call cnvtglfl(oadata(196),rec(196))	! magnitude of 11th order sinusoid
      call cnvtglfl(oadata(197),rec(197))	! phase angle of 11th order sinusoid
      call cnvtglfl(oadata(198),rec(198))	! magnitude of 12th order sinusoid
      call cnvtglfl(oadata(199),rec(199))	! phase angle of 12th order sinusoid
      call cnvtglfl(oadata(200),rec(200))	! magnitude of 13th order sinusoid
      call cnvtglfl(oadata(201),rec(201))	! phase angle of 13th order sinusoid
      call cnvtglfl(oadata(202),rec(202))	! magnitude of 14th order sinusoid
      call cnvtglfl(oadata(203),rec(203))	! phase angle of 14th order sinusoid
      call cnvtglfl(oadata(204),rec(204))	! magnitude of 15th order sinusoid
      call cnvtglfl(oadata(205),rec(205))	! phase angle of 15th order sinusoid
      call cnvtint(oadata(206))
      rec(206) = dfloat(oadata(206))		! number of monomial sinusoids
      call cnvtint(oadata(207))
      rec(207) = dfloat(oadata(207))		! order of applicable sinusoid
      call cnvtint(oadata(208))
      rec(208) = dfloat(oadata(208))		! order of 1st monomial sinusoid
      call cnvtglfl(oadata(209),rec(209))	! magnitude of 1st monomial sinusoid
      call cnvtglfl(oadata(210),rec(210))	! phase angle of 1st monomial sinusoid
      call cnvtglfl(oadata(211),rec(211))	! angle from epoch where 1st monomial is zero
      call cnvtint(oadata(212))
      rec(212) = dfloat(oadata(212))		! order of applicable sinusoid
      call cnvtint(oadata(213))
      rec(213) = dfloat(oadata(213))		! order of 2nd monomial sinusoid
      call cnvtglfl(oadata(214),rec(214))	! magnitude of 2nd monomial sinusoid
      call cnvtglfl(oadata(215),rec(215))	! phase angle of 2nd monomial sinusoid
      call cnvtglfl(oadata(216),rec(216))	! angle from epoch where 2nd monomial is zero
      call cnvtint(oadata(217))
      rec(217) = dfloat(oadata(217))		! order of applicable sinusoid
      call cnvtint(oadata(218))
      rec(218) = dfloat(oadata(218))		! order of 3rd monomial sinusoid
      call cnvtglfl(oadata(219),rec(219))	! magnitude of 3rd monomial sinusoid
      call cnvtglfl(oadata(220),rec(220))	! phase angle of 3rd monomial sinusoid
      call cnvtglfl(oadata(221),rec(221))	! angle from epoch where 3rd monomial is zero
      call cnvtint(oadata(222))
      rec(222) = dfloat(oadata(222))		! order of applicable sinusoid
      call cnvtint(oadata(223))
      rec(223) = dfloat(oadata(223))		! order of 4th monomial sinusoid
      call cnvtglfl(oadata(224),rec(224))	! magnitude of 4th monomial sinusoid
      call cnvtglfl(oadata(225),rec(225))	! phase angle of 4th monomial sinusoid
      call cnvtglfl(oadata(226),rec(226))	! angle from epoch where 4th monomial is zero
c *****************************************************************************
c *****************************************************************************
c     block of code applys to roll misalignment angle:
c *****************************************************************************
c *****************************************************************************
      call cnvtglfl(oadata(227),rec(227))	! exponential magnitude 
      call cnvtglfl(oadata(228),rec(228))	! exponential time constant
      call cnvtglfl(oadata(229),rec(229))	! constant, mean attitude angle
      call cnvtint(oadata(230))
      rec(230) = dfloat(oadata(230))		! number of sinusoids/angles
c *****************************************************************************
c     magnitude/phase angle of first->fifteenth order sinusoids 
c *****************************************************************************
      call cnvtglfl(oadata(231),rec(231))	! magnitude of 1st order sinusoid
      call cnvtglfl(oadata(232),rec(232))	! phase angle of 1st order sinusoid
      call cnvtglfl(oadata(233),rec(233))	! magnitude of 2nd order sinusoid
      call cnvtglfl(oadata(234),rec(234))	! phase angle of 2nd order sinusoid
      call cnvtglfl(oadata(235),rec(235))	! magnitude of 3rd order sinusoid
      call cnvtglfl(oadata(236),rec(236))	! phase angle of 3rd order sinusoid
      call cnvtglfl(oadata(237),rec(237))	! magnitude of 4th order sinusoid
      call cnvtglfl(oadata(238),rec(238))	! phase angle of 4th order sinusoid
      call cnvtglfl(oadata(239),rec(239))	! magnitude of 5th order sinusoid
      call cnvtglfl(oadata(240),rec(240))	! phase angle of 5th order sinusoid
      call cnvtglfl(oadata(241),rec(241))	! magnitude of 6th order sinusoid
      call cnvtglfl(oadata(242),rec(242))	! phase angle of 6th order sinusoid
      call cnvtglfl(oadata(243),rec(243))	! magnitude of 7th order sinusoid
      call cnvtglfl(oadata(244),rec(244))	! phase angle of 7th order sinusoid
      call cnvtglfl(oadata(245),rec(245))	! magnitude of 8th order sinusoid
      call cnvtglfl(oadata(246),rec(246))	! phase angle of 8th order sinusoid
      call cnvtglfl(oadata(247),rec(247))	! magnitude of 9th order sinusoid
      call cnvtglfl(oadata(248),rec(248))	! phase angle of 9th order sinusoid
      call cnvtglfl(oadata(249),rec(249))	! magnitude of 10th order sinusoid
      call cnvtglfl(oadata(250),rec(250))	! phase angle of 10th order sinusoid
      call cnvtglfl(oadata(251),rec(251))	! magnitude of 11th order sinusoid
      call cnvtglfl(oadata(252),rec(252))	! phase angle of 11th order sinusoid
      call cnvtglfl(oadata(253),rec(253))	! magnitude of 12th order sinusoid
      call cnvtglfl(oadata(254),rec(254))	! phase angle of 12th order sinusoid
      call cnvtglfl(oadata(255),rec(255))	! magnitude of 13th order sinusoid
      call cnvtglfl(oadata(256),rec(256))	! phase angle of 13th order sinusoid
      call cnvtglfl(oadata(257),rec(257))	! magnitude of 14th order sinusoid
      call cnvtglfl(oadata(258),rec(258))	! phase angle of 14th order sinusoid
      call cnvtglfl(oadata(259),rec(259))	! magnitude of 15th order sinusoid
      call cnvtglfl(oadata(260),rec(260))	! phase angle of 15th order sinusoid
      call cnvtint(oadata(261))
      rec(261) = dfloat(oadata(261))		! number of monomial sinusoids
      call cnvtint(oadata(262))
      rec(262) = dfloat(oadata(262))		! order of applicable sinusoid
      call cnvtint(oadata(263))
      rec(263) = dfloat(oadata(263))		! order of 1st monomial sinusoid
      call cnvtglfl(oadata(264),rec(264))	! magnitude of 1st monomial sinusoid
      call cnvtglfl(oadata(265),rec(265))	! phase angle of 1st monomial sinusoid
      call cnvtglfl(oadata(266),rec(266))	! angle from epoch where 1st monomial is zero
      call cnvtint(oadata(267))
      rec(267) = dfloat(oadata(267))		! order of applicable sinusoid
      call cnvtint(oadata(268))
      rec(268) = dfloat(oadata(268))		! order of 2nd monomial sinusoid
      call cnvtglfl(oadata(269),rec(269))	! magnitude of 2nd monomial sinusoid
      call cnvtglfl(oadata(270),rec(270))	! phase angle of 2nd monomial sinusoid
      call cnvtglfl(oadata(271),rec(271))	! angle from epoch where 2nd monomial is zero
      call cnvtint(oadata(272))
      rec(272) = dfloat(oadata(272))		! order of applicable sinusoid
      call cnvtint(oadata(273))
      rec(273) = dfloat(oadata(273))		! order of 3rd monomial sinusoid
      call cnvtglfl(oadata(274),rec(274))	! magnitude of 3rd monomial sinusoid
      call cnvtglfl(oadata(275),rec(275))	! phase angle of 3rd monomial sinusoid
      call cnvtglfl(oadata(276),rec(276))	! angle from epoch where 3rd monomial is zero
      call cnvtint(oadata(277))
      rec(277) = dfloat(oadata(277))		! order of applicable sinusoid
      call cnvtint(oadata(278))
      rec(278) = dfloat(oadata(278))		! order of 4th monomial sinusoid
      call cnvtglfl(oadata(279),rec(279))	! magnitude of 4th monomial sinusoid
      call cnvtglfl(oadata(280),rec(280))	! phase angle of 4th monomial sinusoid
      call cnvtglfl(oadata(281),rec(281))	! angle from epoch where 4th monomial is zero
c *****************************************************************************
c *****************************************************************************
c     block of code applys to pitch misalignment angle:
c *****************************************************************************
c *****************************************************************************
      call cnvtglfl(oadata(282),rec(282))	! exponential magnitude 
      call cnvtglfl(oadata(283),rec(283))	! exponential time constant
      call cnvtglfl(oadata(284),rec(284))	! constant, mean attitude angle
      call cnvtint(oadata(285))
      rec(285) = dfloat(oadata(285))		! number of sinusoids/angles
c *****************************************************************************
c     magnitude/phase angle of first->fifteenth order sinusoids 
c *****************************************************************************
      call cnvtglfl(oadata(286),rec(286))	! magnitude of 1st order sinusoid
      call cnvtglfl(oadata(287),rec(287))	! phase angle of 1st order sinusoid
      call cnvtglfl(oadata(288),rec(288))	! magnitude of 2nd order sinusoid
      call cnvtglfl(oadata(289),rec(289))	! phase angle of 2nd order sinusoid
      call cnvtglfl(oadata(290),rec(290))	! magnitude of 3rd order sinusoid
      call cnvtglfl(oadata(291),rec(291))	! phase angle of 3rd order sinusoid
      call cnvtglfl(oadata(292),rec(292))	! magnitude of 4th order sinusoid
      call cnvtglfl(oadata(293),rec(293))	! phase angle of 4th order sinusoid
      call cnvtglfl(oadata(294),rec(294))	! magnitude of 5th order sinusoid
      call cnvtglfl(oadata(295),rec(295))	! phase angle of 5th order sinusoid
      call cnvtglfl(oadata(296),rec(296))	! magnitude of 6th order sinusoid
      call cnvtglfl(oadata(297),rec(297))	! phase angle of 6th order sinusoid
      call cnvtglfl(oadata(298),rec(298))	! magnitude of 7th order sinusoid
      call cnvtglfl(oadata(299),rec(299))	! phase angle of 7th order sinusoid
      call cnvtglfl(oadata(300),rec(300))	! magnitude of 8th order sinusoid
      call cnvtglfl(oadata(301),rec(301))	! phase angle of 8th order sinusoid
      call cnvtglfl(oadata(302),rec(302))	! magnitude of 9th order sinusoid
      call cnvtglfl(oadata(303),rec(303))	! phase angle of 9th order sinusoid
      call cnvtglfl(oadata(304),rec(304))	! magnitude of 10th order sinusoid
      call cnvtglfl(oadata(305),rec(305))	! phase angle of 10th order sinusoid
      call cnvtglfl(oadata(306),rec(306))	! magnitude of 11th order sinusoid
      call cnvtglfl(oadata(307),rec(307))	! phase angle of 11th order sinusoid
      call cnvtglfl(oadata(308),rec(308))	! magnitude of 12th order sinusoid
      call cnvtglfl(oadata(309),rec(309))	! phase angle of 12th order sinusoid
      call cnvtglfl(oadata(310),rec(310))	! magnitude of 13th order sinusoid
      call cnvtglfl(oadata(311),rec(311))	! phase angle of 13th order sinusoid
      call cnvtglfl(oadata(312),rec(312))	! magnitude of 14th order sinusoid
      call cnvtglfl(oadata(313),rec(313))	! phase angle of 14th order sinusoid
      call cnvtglfl(oadata(314),rec(314))	! magnitude of 15th order sinusoid
      call cnvtglfl(oadata(315),rec(315))	! phase angle of 15th order sinusoid
      call cnvtint(oadata(316))
      rec(316) = dfloat(oadata(316))		! number of monomial sinusoids
      call cnvtint(oadata(317))
      rec(317) = dfloat(oadata(317))		! order of applicable sinusoid
      call cnvtint(oadata(318))
      rec(318) = dfloat(oadata(318))		! order of 1st monomial sinusoid
      call cnvtglfl(oadata(319),rec(319))	! magnitude of 1st monomial sinusoid
      call cnvtglfl(oadata(320),rec(320))	! phase angle of 1st monomial sinusoid
      call cnvtglfl(oadata(321),rec(321))	! angle from epoch where 1st monomial is zero
      call cnvtint(oadata(322))
      rec(322) = dfloat(oadata(322))		! order of applicable sinusoid
      call cnvtint(oadata(323))
      rec(323) = dfloat(oadata(323))		! order of 2nd monomial sinusoid
      call cnvtglfl(oadata(324),rec(324))	! magnitude of 2nd monomial sinusoid
      call cnvtglfl(oadata(325),rec(325))	! phase angle of 2nd monomial sinusoid
      call cnvtglfl(oadata(326),rec(326))	! angle from epoch where 2nd monomial is zero
      call cnvtint(oadata(327))
      rec(327) = dfloat(oadata(327))		! order of applicable sinusoid
      call cnvtint(oadata(328))
      rec(328) = dfloat(oadata(328))		! order of 3rd monomial sinusoid
      call cnvtglfl(oadata(329),rec(329))	! magnitude of 3rd monomial sinusoid
      call cnvtglfl(oadata(330),rec(330))	! phase angle of 3rd monomial sinusoid
      call cnvtglfl(oadata(331),rec(331))	! angle from epoch where 3rd monomial is zero
      call cnvtint(oadata(332))
      rec(332) = dfloat(oadata(332))		! order of applicable sinusoid
      call cnvtint(oadata(333))
      rec(333) = dfloat(oadata(333))		! order of 4th monomial sinusoid
      call cnvtglfl(oadata(334),rec(334))	! magnitude of 4th monomial sinusoid
      call cnvtglfl(oadata(335),rec(335))	! phase angle of 4th monomial sinusoid
      call cnvtglfl(oadata(336),rec(336))	! angle from epoch where 4th monomial is zero
c *****************************************************************************
c     print out the results to the user
c *****************************************************************************
      if (print_flag) then
        do 200 i = 2,336
          write(6,1100) i,rec(i)
200     continue
      endif
c *****************************************************************************
c     print termination message, then skip error message
c *****************************************************************************
      print *,'goes_oa obtained. normal termination!'
      goto 400
c *****************************************************************************
c     write i/o error message for input read, then stop!
c *****************************************************************************
300   write(6,1000)
      stop
c
400   continue
c
      print *,'*************************'
c
c *****************************************************************************
c     format statements:
c *****************************************************************************
c
1000  format(t2,'goes_oa: io error on read of input data....')
1100  format(t2,'rec = ',i3.3,', value = ',e22.15)
c
      goto 9999

1010  write(6,*)'error opening gwc o and a'
      istatus = -1

9999  return
      end
