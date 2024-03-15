      module params
!$$$  subprogram documentation block
!                .      .    .                                       .
! module:    params
!   prgmmr: gilbert         org: w/np11    date: 2001-06-05
!
! abstract: this fortran module contains info on all the available 
!           grib parameters.
!
! program history log:
! 2000-05-11  gilbert
! 2003-08-07  gilbert  -  added more parameters
! 2003-09-26  gilbert  -  added more parameters
! 2005-11-17  gordon   -  added more parameters for the wave & smoke models 
!
! usage:    use params
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$

      integer,parameter :: maxparam=239

      type gribparam
          integer :: g1tblver
          integer :: grib1val
          integer :: grib2dsc
          integer :: grib2cat
          integer :: grib2num
          character(len=8) :: abbrev
      end type gribparam

      type(gribparam),dimension(maxparam) :: paramlist

      data paramlist(1)%g1tblver /2/ 
      data paramlist(1)%grib1val /1/ 
      data paramlist(1)%grib2dsc /0/ 
      data paramlist(1)%grib2cat /3/ 
      data paramlist(1)%grib2num /0/ 
      data paramlist(1)%abbrev   /'pres    '/ 
      data paramlist(2)%g1tblver /2/ 
      data paramlist(2)%grib1val /2/ 
      data paramlist(2)%grib2dsc /0/ 
      data paramlist(2)%grib2cat /3/ 
      data paramlist(2)%grib2num /1/ 
      data paramlist(2)%abbrev   /'prmsl   '/ 
      data paramlist(3)%g1tblver /2/ 
      data paramlist(3)%grib1val /3/ 
      data paramlist(3)%grib2dsc /0/ 
      data paramlist(3)%grib2cat /3/ 
      data paramlist(3)%grib2num /2/ 
      data paramlist(3)%abbrev   /'ptend   '/ 
      data paramlist(4)%g1tblver /2/ 
      data paramlist(4)%grib1val /4/ 
      data paramlist(4)%grib2dsc /0/ 
      data paramlist(4)%grib2cat /2/ 
      data paramlist(4)%grib2num /14/ 
      data paramlist(4)%abbrev   /'pvort   '/ 
      data paramlist(5)%g1tblver /2/ 
      data paramlist(5)%grib1val /5/ 
      data paramlist(5)%grib2dsc /0/ 
      data paramlist(5)%grib2cat /3/ 
      data paramlist(5)%grib2num /3/ 
      data paramlist(5)%abbrev   /'icaht   '/ 
      data paramlist(6)%g1tblver /2/ 
      data paramlist(6)%grib1val /6/ 
      data paramlist(6)%grib2dsc /0/ 
      data paramlist(6)%grib2cat /3/ 
      data paramlist(6)%grib2num /4/ 
      data paramlist(6)%abbrev   /'gp      '/ 
      data paramlist(7)%g1tblver /2/ 
      data paramlist(7)%grib1val /7/ 
      data paramlist(7)%grib2dsc /0/ 
      data paramlist(7)%grib2cat /3/ 
      data paramlist(7)%grib2num /5/ 
      data paramlist(7)%abbrev   /'hgt     '/ 
      data paramlist(8)%g1tblver /2/ 
      data paramlist(8)%grib1val /8/ 
      data paramlist(8)%grib2dsc /0/ 
      data paramlist(8)%grib2cat /3/ 
      data paramlist(8)%grib2num /6/ 
      data paramlist(8)%abbrev   /'dist    '/ 
      data paramlist(9)%g1tblver /2/ 
      data paramlist(9)%grib1val /9/ 
      data paramlist(9)%grib2dsc /0/ 
      data paramlist(9)%grib2cat /3/ 
      data paramlist(9)%grib2num /7/ 
      data paramlist(9)%abbrev   /'hstdv   '/ 
      data paramlist(10)%g1tblver /2/ 
      data paramlist(10)%grib1val /10/ 
      data paramlist(10)%grib2dsc /0/ 
      data paramlist(10)%grib2cat /14/ 
      data paramlist(10)%grib2num /0/ 
      data paramlist(10)%abbrev   /'tozne   '/ 
      data paramlist(11)%g1tblver /2/ 
      data paramlist(11)%grib1val /11/ 
      data paramlist(11)%grib2dsc /0/ 
      data paramlist(11)%grib2cat /0/ 
      data paramlist(11)%grib2num /0/ 
      data paramlist(11)%abbrev   /'tmp     '/ 
      data paramlist(12)%g1tblver /2/ 
      data paramlist(12)%grib1val /12/ 
      data paramlist(12)%grib2dsc /0/ 
      data paramlist(12)%grib2cat /0/ 
      data paramlist(12)%grib2num /1/ 
      data paramlist(12)%abbrev   /'vtmp    '/ 
      data paramlist(13)%g1tblver /2/ 
      data paramlist(13)%grib1val /13/ 
      data paramlist(13)%grib2dsc /0/ 
      data paramlist(13)%grib2cat /0/ 
      data paramlist(13)%grib2num /2/ 
      data paramlist(13)%abbrev   /'pot     '/ 
      data paramlist(14)%g1tblver /2/ 
      data paramlist(14)%grib1val /14/ 
      data paramlist(14)%grib2dsc /0/ 
      data paramlist(14)%grib2cat /0/ 
      data paramlist(14)%grib2num /3/ 
      data paramlist(14)%abbrev   /'epot    '/ 
      data paramlist(15)%g1tblver /2/ 
      data paramlist(15)%grib1val /15/ 
      data paramlist(15)%grib2dsc /0/ 
      data paramlist(15)%grib2cat /0/ 
      data paramlist(15)%grib2num /4/ 
      data paramlist(15)%abbrev   /'t max   '/ 
      data paramlist(16)%g1tblver /2/ 
      data paramlist(16)%grib1val /16/ 
      data paramlist(16)%grib2dsc /0/ 
      data paramlist(16)%grib2cat /0/ 
      data paramlist(16)%grib2num /5/ 
      data paramlist(16)%abbrev   /'t min   '/ 
      data paramlist(17)%g1tblver /2/ 
      data paramlist(17)%grib1val /17/ 
      data paramlist(17)%grib2dsc /0/ 
      data paramlist(17)%grib2cat /0/ 
      data paramlist(17)%grib2num /6/ 
      data paramlist(17)%abbrev   /'dpt     '/ 
      data paramlist(18)%g1tblver /2/ 
      data paramlist(18)%grib1val /18/ 
      data paramlist(18)%grib2dsc /0/ 
      data paramlist(18)%grib2cat /0/ 
      data paramlist(18)%grib2num /7/ 
      data paramlist(18)%abbrev   /'depr    '/ 
      data paramlist(19)%g1tblver /2/ 
      data paramlist(19)%grib1val /19/ 
      data paramlist(19)%grib2dsc /0/ 
      data paramlist(19)%grib2cat /0/ 
      data paramlist(19)%grib2num /8/ 
      data paramlist(19)%abbrev   /'lapr    '/ 
      data paramlist(20)%g1tblver /2/ 
      data paramlist(20)%grib1val /20/ 
      data paramlist(20)%grib2dsc /0/ 
      data paramlist(20)%grib2cat /19/ 
      data paramlist(20)%grib2num /0/ 
      data paramlist(20)%abbrev   /'vis     '/ 
      data paramlist(21)%g1tblver /2/ 
      data paramlist(21)%grib1val /21/ 
      data paramlist(21)%grib2dsc /0/ 
      data paramlist(21)%grib2cat /15/ 
      data paramlist(21)%grib2num /6/ 
      data paramlist(21)%abbrev   /'rdsp1   '/ 
      data paramlist(22)%g1tblver /2/ 
      data paramlist(22)%grib1val /22/ 
      data paramlist(22)%grib2dsc /0/ 
      data paramlist(22)%grib2cat /15/ 
      data paramlist(22)%grib2num /7/ 
      data paramlist(22)%abbrev   /'rdsp2   '/ 
      data paramlist(23)%g1tblver /2/ 
      data paramlist(23)%grib1val /23/ 
      data paramlist(23)%grib2dsc /0/ 
      data paramlist(23)%grib2cat /15/ 
      data paramlist(23)%grib2num /8/ 
      data paramlist(23)%abbrev   /'rdsp3   '/ 
      data paramlist(24)%g1tblver /2/ 
      data paramlist(24)%grib1val /24/ 
      data paramlist(24)%grib2dsc /0/ 
      data paramlist(24)%grib2cat /7/ 
      data paramlist(24)%grib2num /0/ 
      data paramlist(24)%abbrev   /'pli     '/ 
      data paramlist(25)%g1tblver /2/ 
      data paramlist(25)%grib1val /25/ 
      data paramlist(25)%grib2dsc /0/ 
      data paramlist(25)%grib2cat /0/ 
      data paramlist(25)%grib2num /9/ 
      data paramlist(25)%abbrev   /'tmp a   '/ 
      data paramlist(26)%g1tblver /2/ 
      data paramlist(26)%grib1val /26/ 
      data paramlist(26)%grib2dsc /0/ 
      data paramlist(26)%grib2cat /3/ 
      data paramlist(26)%grib2num /8/ 
      data paramlist(26)%abbrev   /'presa   '/ 
      data paramlist(27)%g1tblver /2/ 
      data paramlist(27)%grib1val /27/ 
      data paramlist(27)%grib2dsc /0/ 
      data paramlist(27)%grib2cat /3/ 
      data paramlist(27)%grib2num /9/ 
      data paramlist(27)%abbrev   /'gp a    '/ 
      data paramlist(28)%g1tblver /2/ 
      data paramlist(28)%grib1val /28/ 
      data paramlist(28)%grib2dsc /10/ 
      data paramlist(28)%grib2cat /0/ 
      data paramlist(28)%grib2num /0/ 
      data paramlist(28)%abbrev   /'wvsp1   '/ 
      data paramlist(29)%g1tblver /2/ 
      data paramlist(29)%grib1val /29/ 
      data paramlist(29)%grib2dsc /10/ 
      data paramlist(29)%grib2cat /0/ 
      data paramlist(29)%grib2num /1/ 
      data paramlist(29)%abbrev   /'wvsp2   '/ 
      data paramlist(30)%g1tblver /2/ 
      data paramlist(30)%grib1val /30/ 
      data paramlist(30)%grib2dsc /10/ 
      data paramlist(30)%grib2cat /0/ 
      data paramlist(30)%grib2num /2/ 
      data paramlist(30)%abbrev   /'wvsp3   '/ 
      data paramlist(31)%g1tblver /2/ 
      data paramlist(31)%grib1val /31/ 
      data paramlist(31)%grib2dsc /0/ 
      data paramlist(31)%grib2cat /2/ 
      data paramlist(31)%grib2num /0/ 
      data paramlist(31)%abbrev   /'wdir    '/ 
      data paramlist(32)%g1tblver /2/ 
      data paramlist(32)%grib1val /32/ 
      data paramlist(32)%grib2dsc /0/ 
      data paramlist(32)%grib2cat /2/ 
      data paramlist(32)%grib2num /1/ 
      data paramlist(32)%abbrev   /'wind    '/ 
      data paramlist(33)%g1tblver /2/ 
      data paramlist(33)%grib1val /33/ 
      data paramlist(33)%grib2dsc /0/ 
      data paramlist(33)%grib2cat /2/ 
      data paramlist(33)%grib2num /2/ 
      data paramlist(33)%abbrev   /'u grd   '/ 
      data paramlist(34)%g1tblver /2/ 
      data paramlist(34)%grib1val /34/ 
      data paramlist(34)%grib2dsc /0/ 
      data paramlist(34)%grib2cat /2/ 
      data paramlist(34)%grib2num /3/ 
      data paramlist(34)%abbrev   /'v grd   '/ 
      data paramlist(35)%g1tblver /2/ 
      data paramlist(35)%grib1val /35/ 
      data paramlist(35)%grib2dsc /0/ 
      data paramlist(35)%grib2cat /2/ 
      data paramlist(35)%grib2num /4/ 
      data paramlist(35)%abbrev   /'strm    '/ 
      data paramlist(36)%g1tblver /2/ 
      data paramlist(36)%grib1val /36/ 
      data paramlist(36)%grib2dsc /0/ 
      data paramlist(36)%grib2cat /2/ 
      data paramlist(36)%grib2num /5/ 
      data paramlist(36)%abbrev   /'v pot   '/ 
      data paramlist(37)%g1tblver /2/ 
      data paramlist(37)%grib1val /37/ 
      data paramlist(37)%grib2dsc /0/ 
      data paramlist(37)%grib2cat /2/ 
      data paramlist(37)%grib2num /6/ 
      data paramlist(37)%abbrev   /'mntsf   '/ 
      data paramlist(38)%g1tblver /2/ 
      data paramlist(38)%grib1val /38/ 
      data paramlist(38)%grib2dsc /0/ 
      data paramlist(38)%grib2cat /2/ 
      data paramlist(38)%grib2num /7/ 
      data paramlist(38)%abbrev   /'sgcvv   '/ 
      data paramlist(39)%g1tblver /2/ 
      data paramlist(39)%grib1val /39/ 
      data paramlist(39)%grib2dsc /0/ 
      data paramlist(39)%grib2cat /2/ 
      data paramlist(39)%grib2num /8/ 
      data paramlist(39)%abbrev   /'v vel   '/ 
      data paramlist(40)%g1tblver /2/ 
      data paramlist(40)%grib1val /40/ 
      data paramlist(40)%grib2dsc /0/ 
      data paramlist(40)%grib2cat /2/ 
      data paramlist(40)%grib2num /9/ 
      data paramlist(40)%abbrev   /'dzdt    '/ 
      data paramlist(41)%g1tblver /2/ 
      data paramlist(41)%grib1val /41/ 
      data paramlist(41)%grib2dsc /0/ 
      data paramlist(41)%grib2cat /2/ 
      data paramlist(41)%grib2num /10/ 
      data paramlist(41)%abbrev   /'abs v   '/ 
      data paramlist(42)%g1tblver /2/ 
      data paramlist(42)%grib1val /42/ 
      data paramlist(42)%grib2dsc /0/ 
      data paramlist(42)%grib2cat /2/ 
      data paramlist(42)%grib2num /11/ 
      data paramlist(42)%abbrev   /'abs d   '/ 
      data paramlist(43)%g1tblver /2/ 
      data paramlist(43)%grib1val /43/ 
      data paramlist(43)%grib2dsc /0/ 
      data paramlist(43)%grib2cat /2/ 
      data paramlist(43)%grib2num /12/ 
      data paramlist(43)%abbrev   /'rel v   '/ 
      data paramlist(44)%g1tblver /2/ 
      data paramlist(44)%grib1val /44/ 
      data paramlist(44)%grib2dsc /0/ 
      data paramlist(44)%grib2cat /2/ 
      data paramlist(44)%grib2num /13/ 
      data paramlist(44)%abbrev   /'rel d   '/ 
      data paramlist(45)%g1tblver /2/ 
      data paramlist(45)%grib1val /45/ 
      data paramlist(45)%grib2dsc /0/ 
      data paramlist(45)%grib2cat /2/ 
      data paramlist(45)%grib2num /15/ 
      data paramlist(45)%abbrev   /'vucsh   '/ 
      data paramlist(46)%g1tblver /2/ 
      data paramlist(46)%grib1val /46/ 
      data paramlist(46)%grib2dsc /0/ 
      data paramlist(46)%grib2cat /2/ 
      data paramlist(46)%grib2num /16/ 
      data paramlist(46)%abbrev   /'vvcsh   '/ 
      data paramlist(47)%g1tblver /2/ 
      data paramlist(47)%grib1val /47/ 
      data paramlist(47)%grib2dsc /10/ 
      data paramlist(47)%grib2cat /1/ 
      data paramlist(47)%grib2num /0/ 
      data paramlist(47)%abbrev   /'dir c   '/ 
      data paramlist(48)%g1tblver /2/ 
      data paramlist(48)%grib1val /48/ 
      data paramlist(48)%grib2dsc /10/ 
      data paramlist(48)%grib2cat /1/ 
      data paramlist(48)%grib2num /1/ 
      data paramlist(48)%abbrev   /'sp c    '/ 
      data paramlist(49)%g1tblver /2/ 
      data paramlist(49)%grib1val /49/ 
      data paramlist(49)%grib2dsc /10/ 
      data paramlist(49)%grib2cat /1/ 
      data paramlist(49)%grib2num /2/ 
      data paramlist(49)%abbrev   /'uogrd   '/ 
      data paramlist(50)%g1tblver /2/ 
      data paramlist(50)%grib1val /50/ 
      data paramlist(50)%grib2dsc /10/ 
      data paramlist(50)%grib2cat /1/ 
      data paramlist(50)%grib2num /3/ 
      data paramlist(50)%abbrev   /'vogrd   '/ 
      data paramlist(51)%g1tblver /2/ 
      data paramlist(51)%grib1val /51/ 
      data paramlist(51)%grib2dsc /0/ 
      data paramlist(51)%grib2cat /1/ 
      data paramlist(51)%grib2num /0/ 
      data paramlist(51)%abbrev   /'spf h   '/ 
      data paramlist(52)%g1tblver /2/ 
      data paramlist(52)%grib1val /52/ 
      data paramlist(52)%grib2dsc /0/ 
      data paramlist(52)%grib2cat /1/ 
      data paramlist(52)%grib2num /1/ 
      data paramlist(52)%abbrev   /'r h     '/ 
      data paramlist(53)%g1tblver /2/ 
      data paramlist(53)%grib1val /53/ 
      data paramlist(53)%grib2dsc /0/ 
      data paramlist(53)%grib2cat /1/ 
      data paramlist(53)%grib2num /2/ 
      data paramlist(53)%abbrev   /'mixr    '/ 
      data paramlist(54)%g1tblver /2/ 
      data paramlist(54)%grib1val /54/ 
      data paramlist(54)%grib2dsc /0/ 
      data paramlist(54)%grib2cat /1/ 
      data paramlist(54)%grib2num /3/ 
      data paramlist(54)%abbrev   /'p wat   '/ 
      data paramlist(55)%g1tblver /2/ 
      data paramlist(55)%grib1val /55/ 
      data paramlist(55)%grib2dsc /0/ 
      data paramlist(55)%grib2cat /1/ 
      data paramlist(55)%grib2num /4/ 
      data paramlist(55)%abbrev   /'vapp    '/ 
      data paramlist(56)%g1tblver /2/ 
      data paramlist(56)%grib1val /56/ 
      data paramlist(56)%grib2dsc /0/ 
      data paramlist(56)%grib2cat /1/ 
      data paramlist(56)%grib2num /5/ 
      data paramlist(56)%abbrev   /'sat d   '/ 
      data paramlist(57)%g1tblver /2/ 
      data paramlist(57)%grib1val /57/ 
      data paramlist(57)%grib2dsc /0/ 
      data paramlist(57)%grib2cat /1/ 
      data paramlist(57)%grib2num /6/ 
      data paramlist(57)%abbrev   /'evp     '/ 
      data paramlist(58)%g1tblver /2/ 
      data paramlist(58)%grib1val /58/ 
      data paramlist(58)%grib2dsc /0/ 
      data paramlist(58)%grib2cat /6/ 
      data paramlist(58)%grib2num /0/ 
      data paramlist(58)%abbrev   /'c ice   '/ 
      data paramlist(59)%g1tblver /2/ 
      data paramlist(59)%grib1val /59/ 
      data paramlist(59)%grib2dsc /0/ 
      data paramlist(59)%grib2cat /1/ 
      data paramlist(59)%grib2num /7/ 
      data paramlist(59)%abbrev   /'prate   '/ 
      data paramlist(60)%g1tblver /2/ 
      data paramlist(60)%grib1val /60/ 
      data paramlist(60)%grib2dsc /0/ 
      data paramlist(60)%grib2cat /19/ 
      data paramlist(60)%grib2num /2/ 
      data paramlist(60)%abbrev   /'tstm    '/ 
      data paramlist(61)%g1tblver /2/ 
      data paramlist(61)%grib1val /61/ 
      data paramlist(61)%grib2dsc /0/ 
      data paramlist(61)%grib2cat /1/ 
      data paramlist(61)%grib2num /8/ 
      data paramlist(61)%abbrev   /'a pcp   '/ 
      data paramlist(62)%g1tblver /2/ 
      data paramlist(62)%grib1val /62/ 
      data paramlist(62)%grib2dsc /0/ 
      data paramlist(62)%grib2cat /1/ 
      data paramlist(62)%grib2num /9/ 
      data paramlist(62)%abbrev   /'ncpcp   '/ 
      data paramlist(63)%g1tblver /2/ 
      data paramlist(63)%grib1val /63/ 
      data paramlist(63)%grib2dsc /0/ 
      data paramlist(63)%grib2cat /1/ 
      data paramlist(63)%grib2num /10/ 
      data paramlist(63)%abbrev   /'acpcp   '/ 
      data paramlist(64)%g1tblver /2/ 
      data paramlist(64)%grib1val /64/ 
      data paramlist(64)%grib2dsc /0/ 
      data paramlist(64)%grib2cat /1/ 
      data paramlist(64)%grib2num /12/ 
      data paramlist(64)%abbrev   /'srweq   '/ 
      data paramlist(65)%g1tblver /2/ 
      data paramlist(65)%grib1val /65/ 
      data paramlist(65)%grib2dsc /0/ 
      data paramlist(65)%grib2cat /1/ 
      data paramlist(65)%grib2num /13/ 
      data paramlist(65)%abbrev   /'weasd   '/ 
      data paramlist(66)%g1tblver /2/ 
      data paramlist(66)%grib1val /66/ 
      data paramlist(66)%grib2dsc /0/ 
      data paramlist(66)%grib2cat /1/ 
      data paramlist(66)%grib2num /11/ 
      data paramlist(66)%abbrev   /'sno d   '/ 
      data paramlist(67)%g1tblver /2/ 
      data paramlist(67)%grib1val /67/ 
      data paramlist(67)%grib2dsc /0/ 
      data paramlist(67)%grib2cat /19/ 
      data paramlist(67)%grib2num /3/ 
      data paramlist(67)%abbrev   /'mixht   '/ 
      data paramlist(68)%g1tblver /2/ 
      data paramlist(68)%grib1val /68/ 
      data paramlist(68)%grib2dsc /10/ 
      data paramlist(68)%grib2cat /5/ 
      data paramlist(68)%grib2num /2/ 
      data paramlist(68)%abbrev   /'tthdp   '/ 
      data paramlist(69)%g1tblver /2/ 
      data paramlist(69)%grib1val /69/ 
      data paramlist(69)%grib2dsc /10/ 
      data paramlist(69)%grib2cat /5/ 
      data paramlist(69)%grib2num /0/ 
      data paramlist(69)%abbrev   /'mthd    '/ 
      data paramlist(70)%g1tblver /2/ 
      data paramlist(70)%grib1val /70/ 
      data paramlist(70)%grib2dsc /10/ 
      data paramlist(70)%grib2cat /5/ 
      data paramlist(70)%grib2num /1/ 
      data paramlist(70)%abbrev   /'mth a   '/ 
      data paramlist(71)%g1tblver /2/ 
      data paramlist(71)%grib1val /71/ 
      data paramlist(71)%grib2dsc /0/ 
      data paramlist(71)%grib2cat /6/ 
      data paramlist(71)%grib2num /1/ 
      data paramlist(71)%abbrev   /'t cdc   '/ 
      data paramlist(72)%g1tblver /2/ 
      data paramlist(72)%grib1val /72/ 
      data paramlist(72)%grib2dsc /0/ 
      data paramlist(72)%grib2cat /6/ 
      data paramlist(72)%grib2num /2/ 
      data paramlist(72)%abbrev   /'cdcon   '/ 
      data paramlist(73)%g1tblver /2/ 
      data paramlist(73)%grib1val /73/ 
      data paramlist(73)%grib2dsc /0/ 
      data paramlist(73)%grib2cat /6/ 
      data paramlist(73)%grib2num /3/ 
      data paramlist(73)%abbrev   /'l cdc   '/ 
      data paramlist(74)%g1tblver /2/ 
      data paramlist(74)%grib1val /74/ 
      data paramlist(74)%grib2dsc /0/ 
      data paramlist(74)%grib2cat /6/ 
      data paramlist(74)%grib2num /4/ 
      data paramlist(74)%abbrev   /'m cdc   '/ 
      data paramlist(75)%g1tblver /2/ 
      data paramlist(75)%grib1val /75/ 
      data paramlist(75)%grib2dsc /0/ 
      data paramlist(75)%grib2cat /6/ 
      data paramlist(75)%grib2num /5/ 
      data paramlist(75)%abbrev   /'h cdc   '/ 
      data paramlist(76)%g1tblver /2/ 
      data paramlist(76)%grib1val /76/ 
      data paramlist(76)%grib2dsc /0/ 
      data paramlist(76)%grib2cat /6/ 
      data paramlist(76)%grib2num /6/ 
      data paramlist(76)%abbrev   /'c wat   '/ 
      data paramlist(77)%g1tblver /2/ 
      data paramlist(77)%grib1val /77/ 
      data paramlist(77)%grib2dsc /0/ 
      data paramlist(77)%grib2cat /7/ 
      data paramlist(77)%grib2num /1/ 
      data paramlist(77)%abbrev   /'bli     '/ 
      data paramlist(78)%g1tblver /2/ 
      data paramlist(78)%grib1val /78/ 
      data paramlist(78)%grib2dsc /0/ 
      data paramlist(78)%grib2cat /1/ 
      data paramlist(78)%grib2num /14/ 
      data paramlist(78)%abbrev   /'sno c   '/ 
      data paramlist(79)%g1tblver /2/ 
      data paramlist(79)%grib1val /79/ 
      data paramlist(79)%grib2dsc /0/ 
      data paramlist(79)%grib2cat /1/ 
      data paramlist(79)%grib2num /15/ 
      data paramlist(79)%abbrev   /'sno l   '/ 
      data paramlist(80)%g1tblver /2/ 
      data paramlist(80)%grib1val /80/ 
      data paramlist(80)%grib2dsc /10/ 
      data paramlist(80)%grib2cat /4/ 
      data paramlist(80)%grib2num /0/ 
      data paramlist(80)%abbrev   /'wtmp    '/ 
      data paramlist(81)%g1tblver /2/ 
      data paramlist(81)%grib1val /81/ 
      data paramlist(81)%grib2dsc /2/ 
      data paramlist(81)%grib2cat /0/ 
      data paramlist(81)%grib2num /0/ 
      data paramlist(81)%abbrev   /'land    '/ 
      data paramlist(82)%g1tblver /2/ 
      data paramlist(82)%grib1val /82/ 
      data paramlist(82)%grib2dsc /10/ 
      data paramlist(82)%grib2cat /4/ 
      data paramlist(82)%grib2num /1/ 
      data paramlist(82)%abbrev   /'dsl m   '/ 
      data paramlist(83)%g1tblver /2/ 
      data paramlist(83)%grib1val /83/ 
      data paramlist(83)%grib2dsc /2/ 
      data paramlist(83)%grib2cat /0/ 
      data paramlist(83)%grib2num /1/ 
      data paramlist(83)%abbrev   /'sfc r   '/ 
      data paramlist(84)%g1tblver /2/ 
      data paramlist(84)%grib1val /84/ 
      data paramlist(84)%grib2dsc /0/ 
      data paramlist(84)%grib2cat /19/ 
      data paramlist(84)%grib2num /1/ 
      data paramlist(84)%abbrev   /'albdo   '/ 
      data paramlist(85)%g1tblver /2/ 
      data paramlist(85)%grib1val /85/ 
      data paramlist(85)%grib2dsc /2/ 
      data paramlist(85)%grib2cat /0/ 
      data paramlist(85)%grib2num /2/ 
      data paramlist(85)%abbrev   /'tsoil   '/ 
      data paramlist(86)%g1tblver /2/ 
      data paramlist(86)%grib1val /86/ 
      data paramlist(86)%grib2dsc /2/ 
      data paramlist(86)%grib2cat /0/ 
      data paramlist(86)%grib2num /3/ 
      data paramlist(86)%abbrev   /'soil m  '/ 
      data paramlist(87)%g1tblver /2/ 
      data paramlist(87)%grib1val /87/ 
      data paramlist(87)%grib2dsc /2/ 
      data paramlist(87)%grib2cat /0/ 
      data paramlist(87)%grib2num /4/ 
      data paramlist(87)%abbrev   /'veg     '/ 
      data paramlist(88)%g1tblver /2/ 
      data paramlist(88)%grib1val /88/ 
      data paramlist(88)%grib2dsc /10/ 
      data paramlist(88)%grib2cat /5/ 
      data paramlist(88)%grib2num /3/ 
      data paramlist(88)%abbrev   /'salty   '/ 
      data paramlist(89)%g1tblver /2/ 
      data paramlist(89)%grib1val /89/ 
      data paramlist(89)%grib2dsc /0/ 
      data paramlist(89)%grib2cat /3/ 
      data paramlist(89)%grib2num /10/ 
      data paramlist(89)%abbrev   /'den     '/ 
      data paramlist(90)%g1tblver /2/ 
      data paramlist(90)%grib1val /90/ 
      data paramlist(90)%grib2dsc /2/ 
      data paramlist(90)%grib2cat /0/ 
      data paramlist(90)%grib2num /5/ 
      data paramlist(90)%abbrev   /'watr    '/ 
      data paramlist(91)%g1tblver /2/ 
      data paramlist(91)%grib1val /91/ 
      data paramlist(91)%grib2dsc /10/ 
      data paramlist(91)%grib2cat /2/ 
      data paramlist(91)%grib2num /0/ 
      data paramlist(91)%abbrev   /'ice c   '/ 
      data paramlist(92)%g1tblver /2/ 
      data paramlist(92)%grib1val /92/ 
      data paramlist(92)%grib2dsc /10/ 
      data paramlist(92)%grib2cat /2/ 
      data paramlist(92)%grib2num /1/ 
      data paramlist(92)%abbrev   /'icetk   '/ 
      data paramlist(93)%g1tblver /2/ 
      data paramlist(93)%grib1val /93/ 
      data paramlist(93)%grib2dsc /10/ 
      data paramlist(93)%grib2cat /2/ 
      data paramlist(93)%grib2num /2/ 
      data paramlist(93)%abbrev   /'diced   '/ 
      data paramlist(94)%g1tblver /2/ 
      data paramlist(94)%grib1val /94/ 
      data paramlist(94)%grib2dsc /10/ 
      data paramlist(94)%grib2cat /2/ 
      data paramlist(94)%grib2num /3/ 
      data paramlist(94)%abbrev   /'siced   '/ 
      data paramlist(95)%g1tblver /2/ 
      data paramlist(95)%grib1val /95/ 
      data paramlist(95)%grib2dsc /10/ 
      data paramlist(95)%grib2cat /2/ 
      data paramlist(95)%grib2num /4/ 
      data paramlist(95)%abbrev   /'u ice   '/ 
      data paramlist(96)%g1tblver /2/ 
      data paramlist(96)%grib1val /96/ 
      data paramlist(96)%grib2dsc /10/ 
      data paramlist(96)%grib2cat /2/ 
      data paramlist(96)%grib2num /5/ 
      data paramlist(96)%abbrev   /'v ice   '/ 
      data paramlist(97)%g1tblver /2/ 
      data paramlist(97)%grib1val /97/ 
      data paramlist(97)%grib2dsc /10/ 
      data paramlist(97)%grib2cat /2/ 
      data paramlist(97)%grib2num /6/ 
      data paramlist(97)%abbrev   /'ice g   '/ 
      data paramlist(98)%g1tblver /2/ 
      data paramlist(98)%grib1val /98/ 
      data paramlist(98)%grib2dsc /10/ 
      data paramlist(98)%grib2cat /2/ 
      data paramlist(98)%grib2num /7/ 
      data paramlist(98)%abbrev   /'ice d   '/ 
      data paramlist(99)%g1tblver /2/ 
      data paramlist(99)%grib1val /99/ 
      data paramlist(99)%grib2dsc /0/ 
      data paramlist(99)%grib2cat /1/ 
      data paramlist(99)%grib2num /16/ 
      data paramlist(99)%abbrev   /'sno m   '/ 
      data paramlist(100)%g1tblver /2/ 
      data paramlist(100)%grib1val /100/ 
      data paramlist(100)%grib2dsc /10/ 
      data paramlist(100)%grib2cat /0/ 
      data paramlist(100)%grib2num /3/ 
      data paramlist(100)%abbrev   /'htsgw   '/ 
      data paramlist(101)%g1tblver /2/ 
      data paramlist(101)%grib1val /101/ 
      data paramlist(101)%grib2dsc /10/ 
      data paramlist(101)%grib2cat /0/ 
      data paramlist(101)%grib2num /4/ 
      data paramlist(101)%abbrev   /'wvdir   '/ 
      data paramlist(102)%g1tblver /2/ 
      data paramlist(102)%grib1val /102/ 
      data paramlist(102)%grib2dsc /10/ 
      data paramlist(102)%grib2cat /0/ 
      data paramlist(102)%grib2num /5/ 
      data paramlist(102)%abbrev   /'wvhgt   '/ 
      data paramlist(103)%g1tblver /2/ 
      data paramlist(103)%grib1val /103/ 
      data paramlist(103)%grib2dsc /10/ 
      data paramlist(103)%grib2cat /0/ 
      data paramlist(103)%grib2num /6/ 
      data paramlist(103)%abbrev   /'wvper   '/ 
      data paramlist(104)%g1tblver /2/ 
      data paramlist(104)%grib1val /104/ 
      data paramlist(104)%grib2dsc /10/ 
      data paramlist(104)%grib2cat /0/ 
      data paramlist(104)%grib2num /7/ 
      data paramlist(104)%abbrev   /'swdir   '/ 
      data paramlist(105)%g1tblver /2/ 
      data paramlist(105)%grib1val /105/ 
      data paramlist(105)%grib2dsc /10/ 
      data paramlist(105)%grib2cat /0/ 
      data paramlist(105)%grib2num /8/ 
      data paramlist(105)%abbrev   /'swell   '/ 
      data paramlist(106)%g1tblver /2/ 
      data paramlist(106)%grib1val /106/ 
      data paramlist(106)%grib2dsc /10/ 
      data paramlist(106)%grib2cat /0/ 
      data paramlist(106)%grib2num /9/ 
      data paramlist(106)%abbrev   /'swper   '/ 
      data paramlist(107)%g1tblver /2/ 
      data paramlist(107)%grib1val /107/ 
      data paramlist(107)%grib2dsc /10/ 
      data paramlist(107)%grib2cat /0/ 
      data paramlist(107)%grib2num /10/ 
      data paramlist(107)%abbrev   /'dirpw   '/ 
      data paramlist(108)%g1tblver /2/ 
      data paramlist(108)%grib1val /108/ 
      data paramlist(108)%grib2dsc /10/ 
      data paramlist(108)%grib2cat /0/ 
      data paramlist(108)%grib2num /11/ 
      data paramlist(108)%abbrev   /'perpw   '/ 
      data paramlist(109)%g1tblver /2/ 
      data paramlist(109)%grib1val /109/ 
      data paramlist(109)%grib2dsc /10/ 
      data paramlist(109)%grib2cat /0/ 
      data paramlist(109)%grib2num /12/ 
      data paramlist(109)%abbrev   /'dirsw   '/ 
      data paramlist(110)%g1tblver /2/ 
      data paramlist(110)%grib1val /110/ 
      data paramlist(110)%grib2dsc /10/ 
      data paramlist(110)%grib2cat /0/ 
      data paramlist(110)%grib2num /13/ 
      data paramlist(110)%abbrev   /'persw   '/ 
      data paramlist(111)%g1tblver /2/ 
      data paramlist(111)%grib1val /111/ 
      data paramlist(111)%grib2dsc /0/ 
      data paramlist(111)%grib2cat /4/ 
      data paramlist(111)%grib2num /0/ 
      data paramlist(111)%abbrev   /'nswrs   '/ 
      data paramlist(112)%g1tblver /2/ 
      data paramlist(112)%grib1val /112/ 
      data paramlist(112)%grib2dsc /0/ 
      data paramlist(112)%grib2cat /5/ 
      data paramlist(112)%grib2num /0/ 
      data paramlist(112)%abbrev   /'nlwrs   '/ 
      data paramlist(113)%g1tblver /2/ 
      data paramlist(113)%grib1val /113/ 
      data paramlist(113)%grib2dsc /0/ 
      data paramlist(113)%grib2cat /4/ 
      data paramlist(113)%grib2num /1/ 
      data paramlist(113)%abbrev   /'nswrt   '/ 
      data paramlist(114)%g1tblver /2/ 
      data paramlist(114)%grib1val /114/ 
      data paramlist(114)%grib2dsc /0/ 
      data paramlist(114)%grib2cat /5/ 
      data paramlist(114)%grib2num /1/ 
      data paramlist(114)%abbrev   /'nlwrt   '/ 
      data paramlist(115)%g1tblver /2/ 
      data paramlist(115)%grib1val /115/ 
      data paramlist(115)%grib2dsc /0/ 
      data paramlist(115)%grib2cat /5/ 
      data paramlist(115)%grib2num /2/ 
      data paramlist(115)%abbrev   /'lwavr   '/ 
      data paramlist(116)%g1tblver /2/ 
      data paramlist(116)%grib1val /116/ 
      data paramlist(116)%grib2dsc /0/ 
      data paramlist(116)%grib2cat /4/ 
      data paramlist(116)%grib2num /2/ 
      data paramlist(116)%abbrev   /'swavr   '/ 
      data paramlist(117)%g1tblver /2/ 
      data paramlist(117)%grib1val /117/ 
      data paramlist(117)%grib2dsc /0/ 
      data paramlist(117)%grib2cat /4/ 
      data paramlist(117)%grib2num /3/ 
      data paramlist(117)%abbrev   /'g rad   '/ 
      data paramlist(118)%g1tblver /2/ 
      data paramlist(118)%grib1val /118/ 
      data paramlist(118)%grib2dsc /0/ 
      data paramlist(118)%grib2cat /4/ 
      data paramlist(118)%grib2num /4/ 
      data paramlist(118)%abbrev   /'brtmp   '/ 
      data paramlist(119)%g1tblver /2/ 
      data paramlist(119)%grib1val /119/ 
      data paramlist(119)%grib2dsc /0/ 
      data paramlist(119)%grib2cat /4/ 
      data paramlist(119)%grib2num /5/ 
      data paramlist(119)%abbrev   /'lwrad   '/ 
      data paramlist(120)%g1tblver /2/ 
      data paramlist(120)%grib1val /120/ 
      data paramlist(120)%grib2dsc /0/ 
      data paramlist(120)%grib2cat /4/ 
      data paramlist(120)%grib2num /6/ 
      data paramlist(120)%abbrev   /'swrad   '/ 
      data paramlist(121)%g1tblver /2/ 
      data paramlist(121)%grib1val /121/ 
      data paramlist(121)%grib2dsc /0/ 
      data paramlist(121)%grib2cat /0/ 
      data paramlist(121)%grib2num /10/ 
      data paramlist(121)%abbrev   /'lhtfl   '/ 
      data paramlist(122)%g1tblver /2/ 
      data paramlist(122)%grib1val /122/ 
      data paramlist(122)%grib2dsc /0/ 
      data paramlist(122)%grib2cat /0/ 
      data paramlist(122)%grib2num /11/ 
      data paramlist(122)%abbrev   /'shtfl   '/ 
      data paramlist(123)%g1tblver /2/ 
      data paramlist(123)%grib1val /123/ 
      data paramlist(123)%grib2dsc /0/ 
      data paramlist(123)%grib2cat /2/ 
      data paramlist(123)%grib2num /20/ 
      data paramlist(123)%abbrev   /'blydp   '/ 
      data paramlist(124)%g1tblver /2/ 
      data paramlist(124)%grib1val /124/ 
      data paramlist(124)%grib2dsc /0/ 
      data paramlist(124)%grib2cat /2/ 
      data paramlist(124)%grib2num /17/ 
      data paramlist(124)%abbrev   /'u flx   '/ 
      data paramlist(125)%g1tblver /2/ 
      data paramlist(125)%grib1val /125/ 
      data paramlist(125)%grib2dsc /0/ 
      data paramlist(125)%grib2cat /2/ 
      data paramlist(125)%grib2num /18/ 
      data paramlist(125)%abbrev   /'v flx   '/ 
      data paramlist(126)%g1tblver /2/ 
      data paramlist(126)%grib1val /126/ 
      data paramlist(126)%grib2dsc /0/ 
      data paramlist(126)%grib2cat /2/ 
      data paramlist(126)%grib2num /19/ 
      data paramlist(126)%abbrev   /'wmixe   '/ 
      data paramlist(127)%g1tblver /2/ 
      data paramlist(127)%grib1val /127/ 
      data paramlist(127)%grib2dsc /255/ 
      data paramlist(127)%grib2cat /255/ 
      data paramlist(127)%grib2num /255/ 
      data paramlist(127)%abbrev   /'img d   '/ 
!
!  grib1 parameters in ncep local table version 2
!  added 8/07/2003
!
      data paramlist(128)%g1tblver /2/ 
      data paramlist(128)%grib1val /229/ 
      data paramlist(128)%grib2dsc /0/ 
      data paramlist(128)%grib2cat /0/ 
      data paramlist(128)%grib2num /192/ 
      data paramlist(128)%abbrev   /'snohf   '/ 
      data paramlist(129)%g1tblver /2/ 
      data paramlist(129)%grib1val /153/ 
      data paramlist(129)%grib2dsc /0/ 
      data paramlist(129)%grib2cat /1/ 
      data paramlist(129)%grib2num /22/ 
      data paramlist(129)%abbrev   /'clwmr   '/ 
      data paramlist(130)%g1tblver /2/ 
      data paramlist(130)%grib1val /140/ 
      data paramlist(130)%grib2dsc /0/ 
      data paramlist(130)%grib2cat /1/ 
      data paramlist(130)%grib2num /192/ 
      data paramlist(130)%abbrev   /'crain   '/ 
      data paramlist(131)%g1tblver /2/ 
      data paramlist(131)%grib1val /141/ 
      data paramlist(131)%grib2dsc /0/ 
      data paramlist(131)%grib2cat /1/ 
      data paramlist(131)%grib2num /193/ 
      data paramlist(131)%abbrev   /'cfrzr   '/ 
      data paramlist(132)%g1tblver /2/ 
      data paramlist(132)%grib1val /142/ 
      data paramlist(132)%grib2dsc /0/ 
      data paramlist(132)%grib2cat /1/ 
      data paramlist(132)%grib2num /194/ 
      data paramlist(132)%abbrev   /'cicep   '/ 
      data paramlist(133)%g1tblver /2/ 
      data paramlist(133)%grib1val /143/ 
      data paramlist(133)%grib2dsc /0/ 
      data paramlist(133)%grib2cat /1/ 
      data paramlist(133)%grib2num /195/ 
      data paramlist(133)%abbrev   /'csnow   '/ 
      data paramlist(134)%g1tblver /2/ 
      data paramlist(134)%grib1val /214/ 
      data paramlist(134)%grib2dsc /0/ 
      data paramlist(134)%grib2cat /1/ 
      data paramlist(134)%grib2num /196/ 
      data paramlist(134)%abbrev   /'cprat   '/ 
      data paramlist(135)%g1tblver /2/ 
      data paramlist(135)%grib1val /135/ 
      data paramlist(135)%grib2dsc /0/ 
      data paramlist(135)%grib2cat /1/ 
      data paramlist(135)%grib2num /197/ 
      data paramlist(135)%abbrev   /'mconv   '/ 
      data paramlist(136)%g1tblver /2/ 
      data paramlist(136)%grib1val /194/ 
      data paramlist(136)%grib2dsc /1/ 
      data paramlist(136)%grib2cat /1/ 
      data paramlist(136)%grib2num /193/ 
      data paramlist(136)%abbrev   /'cpofp   '/ 
      data paramlist(137)%g1tblver /2/ 
      data paramlist(137)%grib1val /228/ 
      data paramlist(137)%grib2dsc /0/ 
      data paramlist(137)%grib2cat /1/ 
      data paramlist(137)%grib2num /199/ 
      data paramlist(137)%abbrev   /'pevap   '/ 
      data paramlist(138)%g1tblver /2/ 
      data paramlist(138)%grib1val /136/ 
      data paramlist(138)%grib2dsc /0/ 
      data paramlist(138)%grib2cat /2/ 
      data paramlist(138)%grib2num /192/ 
      data paramlist(138)%abbrev   /'vw sh   '/ 
      data paramlist(139)%g1tblver /2/ 
      data paramlist(139)%grib1val /172/ 
      data paramlist(139)%grib2dsc /0/ 
      data paramlist(139)%grib2cat /2/ 
      data paramlist(139)%grib2num /193/ 
      data paramlist(139)%abbrev   /'m flx   '/ 
      data paramlist(140)%g1tblver /2/ 
      data paramlist(140)%grib1val /196/ 
      data paramlist(140)%grib2dsc /0/ 
      data paramlist(140)%grib2cat /2/ 
      data paramlist(140)%grib2num /194/ 
      data paramlist(140)%abbrev   /'ustm    '/ 
      data paramlist(141)%g1tblver /2/ 
      data paramlist(141)%grib1val /197/ 
      data paramlist(141)%grib2dsc /0/ 
      data paramlist(141)%grib2cat /2/ 
      data paramlist(141)%grib2num /195/ 
      data paramlist(141)%abbrev   /'vstm    '/ 
      data paramlist(142)%g1tblver /2/ 
      data paramlist(142)%grib1val /252/ 
      data paramlist(142)%grib2dsc /0/ 
      data paramlist(142)%grib2cat /2/ 
      data paramlist(142)%grib2num /196/ 
      data paramlist(142)%abbrev   /'cd      '/ 
      data paramlist(143)%g1tblver /2/ 
      data paramlist(143)%grib1val /253/ 
      data paramlist(143)%grib2dsc /0/ 
      data paramlist(143)%grib2cat /2/ 
      data paramlist(143)%grib2num /197/ 
      data paramlist(143)%abbrev   /'fricv   '/ 
      data paramlist(144)%g1tblver /2/ 
      data paramlist(144)%grib1val /130/ 
      data paramlist(144)%grib2dsc /0/ 
      data paramlist(144)%grib2cat /3/ 
      data paramlist(144)%grib2num /192/ 
      data paramlist(144)%abbrev   /'mslet   '/ 
      data paramlist(145)%g1tblver /2/ 
      data paramlist(145)%grib1val /204/ 
      data paramlist(145)%grib2dsc /0/ 
      data paramlist(145)%grib2cat /4/ 
      data paramlist(145)%grib2num /192/ 
      data paramlist(145)%abbrev   /'dswrf   '/ 
      data paramlist(146)%g1tblver /2/ 
      data paramlist(146)%grib1val /211/ 
      data paramlist(146)%grib2dsc /0/ 
      data paramlist(146)%grib2cat /4/ 
      data paramlist(146)%grib2num /193/ 
      data paramlist(146)%abbrev   /'uswrf   '/ 
      data paramlist(147)%g1tblver /2/ 
      data paramlist(147)%grib1val /205/ 
      data paramlist(147)%grib2dsc /0/ 
      data paramlist(147)%grib2cat /5/ 
      data paramlist(147)%grib2num /192/ 
      data paramlist(147)%abbrev   /'dlwrf   '/ 
      data paramlist(148)%g1tblver /2/ 
      data paramlist(148)%grib1val /212/ 
      data paramlist(148)%grib2dsc /0/ 
      data paramlist(148)%grib2cat /5/ 
      data paramlist(148)%grib2num /193/ 
      data paramlist(148)%abbrev   /'ulwrf   '/ 
      data paramlist(149)%g1tblver /2/ 
      data paramlist(149)%grib1val /213/ 
      data paramlist(149)%grib2dsc /0/ 
      data paramlist(149)%grib2cat /6/ 
      data paramlist(149)%grib2num /192/ 
      data paramlist(149)%abbrev   /'cdlyr   '/ 
      data paramlist(150)%g1tblver /2/ 
      data paramlist(150)%grib1val /132/ 
      data paramlist(150)%grib2dsc /0/ 
      data paramlist(150)%grib2cat /7/ 
      data paramlist(150)%grib2num /193/ 
      data paramlist(150)%abbrev   /'4lftx   '/ 
      data paramlist(151)%g1tblver /2/ 
      data paramlist(151)%grib1val /157/ 
      data paramlist(151)%grib2dsc /0/ 
      data paramlist(151)%grib2cat /7/ 
      data paramlist(151)%grib2num /6/ 
      data paramlist(151)%abbrev   /'cape    '/ 
      data paramlist(152)%g1tblver /2/ 
      data paramlist(152)%grib1val /156/ 
      data paramlist(152)%grib2dsc /0/ 
      data paramlist(152)%grib2cat /7/ 
      data paramlist(152)%grib2num /7/ 
      data paramlist(152)%abbrev   /'cin     '/ 
      data paramlist(153)%g1tblver /2/ 
      data paramlist(153)%grib1val /190/ 
      data paramlist(153)%grib2dsc /0/ 
      data paramlist(153)%grib2cat /7/ 
      data paramlist(153)%grib2num /8/ 
      data paramlist(153)%abbrev   /'hlcy    '/ 
      data paramlist(154)%g1tblver /2/ 
      data paramlist(154)%grib1val /131/ 
      data paramlist(154)%grib2dsc /0/ 
      data paramlist(154)%grib2cat /7/ 
      data paramlist(154)%grib2num /192/ 
      data paramlist(154)%abbrev   /'lft x   '/ 
      data paramlist(155)%g1tblver /2/ 
      data paramlist(155)%grib1val /158/ 
      data paramlist(155)%grib2dsc /0/ 
      data paramlist(155)%grib2cat /19/ 
      data paramlist(155)%grib2num /11/ 
      data paramlist(155)%abbrev   /'tke     '/ 
      data paramlist(156)%g1tblver /2/ 
      data paramlist(156)%grib1val /176/ 
      data paramlist(156)%grib2dsc /0/ 
      data paramlist(156)%grib2cat /191/ 
      data paramlist(156)%grib2num /192/ 
      data paramlist(156)%abbrev   /'nlat    '/ 
      data paramlist(157)%g1tblver /2/ 
      data paramlist(157)%grib1val /177/ 
      data paramlist(157)%grib2dsc /0/ 
      data paramlist(157)%grib2cat /191/ 
      data paramlist(157)%grib2num /193/ 
      data paramlist(157)%abbrev   /'elon    '/ 
      data paramlist(158)%g1tblver /2/ 
      data paramlist(158)%grib1val /234/ 
      data paramlist(158)%grib2dsc /1/ 
      data paramlist(158)%grib2cat /0/ 
      data paramlist(158)%grib2num /192/ 
      data paramlist(158)%abbrev   /'bgrun   '/ 
      data paramlist(159)%g1tblver /2/ 
      data paramlist(159)%grib1val /235/ 
      data paramlist(159)%grib2dsc /1/ 
      data paramlist(159)%grib2cat /0/ 
      data paramlist(159)%grib2num /193/ 
      data paramlist(159)%abbrev   /'ssrun   '/ 
      data paramlist(160)%g1tblver /2/ 
      data paramlist(160)%grib1val /144/ 
      data paramlist(160)%grib2dsc /2/ 
      data paramlist(160)%grib2cat /0/ 
      data paramlist(160)%grib2num /192/ 
      data paramlist(160)%abbrev   /'soilw   '/ 
      data paramlist(161)%g1tblver /2/ 
      data paramlist(161)%grib1val /155/ 
      data paramlist(161)%grib2dsc /2/ 
      data paramlist(161)%grib2cat /0/ 
      data paramlist(161)%grib2num /193/ 
      data paramlist(161)%abbrev   /'gflux   '/ 
      data paramlist(162)%g1tblver /2/ 
      data paramlist(162)%grib1val /207/ 
      data paramlist(162)%grib2dsc /2/ 
      data paramlist(162)%grib2cat /0/ 
      data paramlist(162)%grib2num /194/ 
      data paramlist(162)%abbrev   /'mstav   '/ 
      data paramlist(163)%g1tblver /2/ 
      data paramlist(163)%grib1val /208/ 
      data paramlist(163)%grib2dsc /2/ 
      data paramlist(163)%grib2cat /0/ 
      data paramlist(163)%grib2num /195/ 
      data paramlist(163)%abbrev   /'sfexc   '/ 
      data paramlist(164)%g1tblver /2/ 
      data paramlist(164)%grib1val /223/ 
      data paramlist(164)%grib2dsc /2/ 
      data paramlist(164)%grib2cat /0/ 
      data paramlist(164)%grib2num /196/ 
      data paramlist(164)%abbrev   /'cnwat   '/ 
      data paramlist(165)%g1tblver /2/ 
      data paramlist(165)%grib1val /226/ 
      data paramlist(165)%grib2dsc /2/ 
      data paramlist(165)%grib2cat /0/ 
      data paramlist(165)%grib2num /197/ 
      data paramlist(165)%abbrev   /'bmixl   '/ 
      data paramlist(166)%g1tblver /2/ 
      data paramlist(166)%grib1val /154/ 
      data paramlist(166)%grib2dsc /0/ 
      data paramlist(166)%grib2cat /14/ 
      data paramlist(166)%grib2num /192/ 
      data paramlist(166)%abbrev   /'o3mr    '/ 
      data paramlist(167)%g1tblver /2/ 
      data paramlist(167)%grib1val /222/ 
      data paramlist(167)%grib2dsc /0/ 
      data paramlist(167)%grib2cat /3/ 
      data paramlist(167)%grib2num /193/ 
      data paramlist(167)%abbrev   /'5wavh   '/ 
      data paramlist(168)%g1tblver /2/ 
      data paramlist(168)%grib1val /145/ 
      data paramlist(168)%grib2dsc /0/ 
      data paramlist(168)%grib2cat /1/ 
      data paramlist(168)%grib2num /200/ 
      data paramlist(168)%abbrev   /'pevpr   '/ 
      data paramlist(169)%g1tblver /2/ 
      data paramlist(169)%grib1val /146/ 
      data paramlist(169)%grib2dsc /0/ 
      data paramlist(169)%grib2cat /6/ 
      data paramlist(169)%grib2num /193/ 
      data paramlist(169)%abbrev   /'cwork   '/ 
      data paramlist(170)%g1tblver /2/ 
      data paramlist(170)%grib1val /147/ 
      data paramlist(170)%grib2dsc /0/ 
      data paramlist(170)%grib2cat /3/ 
      data paramlist(170)%grib2num /194/ 
      data paramlist(170)%abbrev   /'u-gwd   '/ 
      data paramlist(171)%g1tblver /2/ 
      data paramlist(171)%grib1val /148/ 
      data paramlist(171)%grib2dsc /0/ 
      data paramlist(171)%grib2cat /3/ 
      data paramlist(171)%grib2num /195/ 
      data paramlist(171)%abbrev   /'v-gwd   '/ 
      data paramlist(172)%g1tblver /2/ 
      data paramlist(172)%grib1val /221/ 
      data paramlist(172)%grib2dsc /0/ 
      data paramlist(172)%grib2cat /3/ 
      data paramlist(172)%grib2num /196/ 
      data paramlist(172)%abbrev   /'hpbl    '/ 
      data paramlist(173)%g1tblver /2/ 
      data paramlist(173)%grib1val /230/ 
      data paramlist(173)%grib2dsc /0/ 
      data paramlist(173)%grib2cat /3/ 
      data paramlist(173)%grib2num /197/ 
      data paramlist(173)%abbrev   /'5wava   '/ 
! added 9/26/2003
      data paramlist(174)%g1tblver /130/
      data paramlist(174)%grib1val /160/
      data paramlist(174)%grib2dsc /2/
      data paramlist(174)%grib2cat /3/
      data paramlist(174)%grib2num /192/  
      data paramlist(174)%abbrev   /'soill   '/
      data paramlist(175)%g1tblver /130/
      data paramlist(175)%grib1val /171/
      data paramlist(175)%grib2dsc /2/
      data paramlist(175)%grib2cat /3/
      data paramlist(175)%grib2num /193/  
      data paramlist(175)%abbrev   /'unknown '/
      data paramlist(176)%g1tblver /130/
      data paramlist(176)%grib1val /219/
      data paramlist(176)%grib2dsc /2/
      data paramlist(176)%grib2cat /0/
      data paramlist(176)%grib2num /201/  
      data paramlist(176)%abbrev   /'wilt    '/
      data paramlist(177)%g1tblver /130/
      data paramlist(177)%grib1val /222/
      data paramlist(177)%grib2dsc /2/
      data paramlist(177)%grib2cat /3/
      data paramlist(177)%grib2num /194/  
      data paramlist(177)%abbrev   /'sltyp   '/
      data paramlist(178)%g1tblver /2/
      data paramlist(178)%grib1val /224/
      data paramlist(178)%grib2dsc /2/
      data paramlist(178)%grib2cat /3/
      data paramlist(178)%grib2num /0/      
      data paramlist(178)%abbrev   /'sotyp   '/
      data paramlist(179)%g1tblver /2/
      data paramlist(179)%grib1val /225/
      data paramlist(179)%grib2dsc /2/
      data paramlist(179)%grib2cat /0/
      data paramlist(179)%grib2num /198/    
      data paramlist(179)%abbrev   /'vgtyp   '/
      data paramlist(180)%g1tblver /130/
      data paramlist(180)%grib1val /230/
      data paramlist(180)%grib2dsc /2/
      data paramlist(180)%grib2cat /3/
      data paramlist(180)%grib2num /195/  
      data paramlist(180)%abbrev   /'smref   '/
      data paramlist(181)%g1tblver /130/
      data paramlist(181)%grib1val /231/
      data paramlist(181)%grib2dsc /2/
      data paramlist(181)%grib2cat /3/
      data paramlist(181)%grib2num /196/  
      data paramlist(181)%abbrev   /'smdry   '/
      data paramlist(182)%g1tblver /2/
      data paramlist(182)%grib1val /238/
      data paramlist(182)%grib2dsc /0/
      data paramlist(182)%grib2cat /1/
      data paramlist(182)%grib2num /201/    
      data paramlist(182)%abbrev   /'snowc   '/
      data paramlist(183)%g1tblver /130/
      data paramlist(183)%grib1val /240/
      data paramlist(183)%grib2dsc /2/
      data paramlist(183)%grib2cat /3/
      data paramlist(183)%grib2num /197/  
      data paramlist(183)%abbrev   /'poros   '/
      data paramlist(184)%g1tblver /129/
      data paramlist(184)%grib1val /131/
      data paramlist(184)%grib2dsc /0/
      data paramlist(184)%grib2cat /1/
      data paramlist(184)%grib2num /202/  
      data paramlist(184)%abbrev   /'frain   '/
      data paramlist(185)%g1tblver /129/
      data paramlist(185)%grib1val /132/
      data paramlist(185)%grib2dsc /0/
      data paramlist(185)%grib2cat /6/
      data paramlist(185)%grib2num /199/  
      data paramlist(185)%abbrev   /'fice    '/
      data paramlist(186)%g1tblver /129/
      data paramlist(186)%grib1val /133/
      data paramlist(186)%grib2dsc /0/
      data paramlist(186)%grib2cat /1/
      data paramlist(186)%grib2num /203/  
      data paramlist(186)%abbrev   /'frime   '/
      data paramlist(187)%g1tblver /129/
      data paramlist(187)%grib1val /134/
      data paramlist(187)%grib2dsc /0/
      data paramlist(187)%grib2cat /6/
      data paramlist(187)%grib2num /194/  
      data paramlist(187)%abbrev   /'cuefi   '/
      data paramlist(188)%g1tblver /129/
      data paramlist(188)%grib1val /135/
      data paramlist(188)%grib2dsc /0/
      data paramlist(188)%grib2cat /6/
      data paramlist(188)%grib2num /195/  
      data paramlist(188)%abbrev   /'tcond   '/
      data paramlist(189)%g1tblver /129/
      data paramlist(189)%grib1val /136/
      data paramlist(189)%grib2dsc /0/
      data paramlist(189)%grib2cat /6/
      data paramlist(189)%grib2num /196/  
      data paramlist(189)%abbrev   /'tcolw   '/
      data paramlist(190)%g1tblver /129/
      data paramlist(190)%grib1val /137/
      data paramlist(190)%grib2dsc /0/
      data paramlist(190)%grib2cat /6/
      data paramlist(190)%grib2num /197/  
      data paramlist(190)%abbrev   /'tcoli   '/
      data paramlist(191)%g1tblver /129/
      data paramlist(191)%grib1val /138/
      data paramlist(191)%grib2dsc /0/
      data paramlist(191)%grib2cat /1/
      data paramlist(191)%grib2num /204/  
      data paramlist(191)%abbrev   /'tcolr   '/
      data paramlist(192)%g1tblver /129/
      data paramlist(192)%grib1val /139/
      data paramlist(192)%grib2dsc /0/
      data paramlist(192)%grib2cat /1/
      data paramlist(192)%grib2num /205/  
      data paramlist(192)%abbrev   /'tcols   '/
      data paramlist(193)%g1tblver /129/
      data paramlist(193)%grib1val /140/
      data paramlist(193)%grib2dsc /0/
      data paramlist(193)%grib2cat /6/
      data paramlist(193)%grib2num /198/  
      data paramlist(193)%abbrev   /'tcolc   '/
      data paramlist(194)%g1tblver /130/
      data paramlist(194)%grib1val /159/
      data paramlist(194)%grib2dsc /0/
      data paramlist(194)%grib2cat /19/
      data paramlist(194)%grib2num /192/ 
      data paramlist(194)%abbrev   /'mxsalb  '/
      data paramlist(195)%g1tblver /130/
      data paramlist(195)%grib1val /170/
      data paramlist(195)%grib2dsc /0/
      data paramlist(195)%grib2cat /19/
      data paramlist(195)%grib2num /193/ 
      data paramlist(195)%abbrev   /'snfalb  '/
      data paramlist(196)%g1tblver /2/
      data paramlist(196)%grib1val /170/
      data paramlist(196)%grib2dsc /0/
      data paramlist(196)%grib2cat /1/
      data paramlist(196)%grib2num /24/  
      data paramlist(196)%abbrev   /'rwmr    '/
      data paramlist(197)%g1tblver /2/
      data paramlist(197)%grib1val /171/
      data paramlist(197)%grib2dsc /0/
      data paramlist(197)%grib2cat /1/
      data paramlist(197)%grib2num /25/  
      data paramlist(197)%abbrev   /'snmr    '/
      data paramlist(198)%g1tblver /130/
      data paramlist(198)%grib1val /181/
      data paramlist(198)%grib2dsc /2/
      data paramlist(198)%grib2cat /0/
      data paramlist(198)%grib2num /199/ 
      data paramlist(198)%abbrev   /'ccond   '/
      data paramlist(199)%g1tblver /130/
      data paramlist(199)%grib1val /203/
      data paramlist(199)%grib2dsc /2/
      data paramlist(199)%grib2cat /0/
      data paramlist(199)%grib2num /200/ 
      data paramlist(199)%abbrev   /'rsmin   '/
      data paramlist(200)%g1tblver /130/
      data paramlist(200)%grib1val /246/
      data paramlist(200)%grib2dsc /2/
      data paramlist(200)%grib2cat /0/
      data paramlist(200)%grib2num /202/ 
      data paramlist(200)%abbrev   /'rcs     '/
      data paramlist(201)%g1tblver /130/
      data paramlist(201)%grib1val /247/
      data paramlist(201)%grib2dsc /2/
      data paramlist(201)%grib2cat /0/
      data paramlist(201)%grib2num /203/ 
      data paramlist(201)%abbrev   /'rct     '/
      data paramlist(202)%g1tblver /130/
      data paramlist(202)%grib1val /248/
      data paramlist(202)%grib2dsc /2/
      data paramlist(202)%grib2cat /0/
      data paramlist(202)%grib2num /204/ 
      data paramlist(202)%abbrev   /'rcq     '/
      data paramlist(203)%g1tblver /130/
      data paramlist(203)%grib1val /249/
      data paramlist(203)%grib2dsc /2/
      data paramlist(203)%grib2cat /0/
      data paramlist(203)%grib2num /205/ 
      data paramlist(203)%abbrev   /'rcsol   '/
      data paramlist(204)%g1tblver /2/
      data paramlist(204)%grib1val /254/
      data paramlist(204)%grib2dsc /0/
      data paramlist(204)%grib2cat /7/
      data paramlist(204)%grib2num /194/ 
      data paramlist(204)%abbrev   /'ri      '/
      data paramlist(205)%g1tblver /129/
      data paramlist(205)%grib1val /190/
      data paramlist(205)%grib2dsc /3/
      data paramlist(205)%grib2cat /1/
      data paramlist(205)%grib2num /192/ 
      data paramlist(205)%abbrev   /'usct    '/
      data paramlist(206)%g1tblver /129/
      data paramlist(206)%grib1val /191/
      data paramlist(206)%grib2dsc /3/
      data paramlist(206)%grib2cat /1/
      data paramlist(206)%grib2num /193/ 
      data paramlist(206)%abbrev   /'vsct    '/
      data paramlist(207)%g1tblver /129/
      data paramlist(207)%grib1val /171/
      data paramlist(207)%grib2dsc /0/
      data paramlist(207)%grib2cat /191/
      data paramlist(207)%grib2num /194/ 
      data paramlist(207)%abbrev   /'tsec    '/
      data paramlist(208)%g1tblver /129/
      data paramlist(208)%grib1val /180/
      data paramlist(208)%grib2dsc /0/
      data paramlist(208)%grib2cat /14/
      data paramlist(208)%grib2num /193/ 
      data paramlist(208)%abbrev   /'ozcon   '/
      data paramlist(209)%g1tblver /129/
      data paramlist(209)%grib1val /181/
      data paramlist(209)%grib2dsc /0/
      data paramlist(209)%grib2cat /14/
      data paramlist(209)%grib2num /194/ 
      data paramlist(209)%abbrev   /'ozcat   '/
      data paramlist(210)%g1tblver /2/
      data paramlist(210)%grib1val /193/
      data paramlist(210)%grib2dsc /1/
      data paramlist(210)%grib2cat /1/
      data paramlist(210)%grib2num /2/   
      data paramlist(210)%abbrev   /'pop     '/
      data paramlist(211)%g1tblver /2/
      data paramlist(211)%grib1val /195/
      data paramlist(211)%grib2dsc /1/
      data paramlist(211)%grib2cat /1/
      data paramlist(211)%grib2num /192/ 
      data paramlist(211)%abbrev   /'cpozp   '/
      data paramlist(212)%g1tblver /2/
      data paramlist(212)%grib1val /180/
      data paramlist(212)%grib2dsc /0/
      data paramlist(212)%grib2cat /2/
      data paramlist(212)%grib2num /22/  
      data paramlist(212)%abbrev   /'gust    '/
! added 11/17/2005 - for wave models
      data paramlist(213)%g1tblver /0/
      data paramlist(213)%grib1val /31/
      data paramlist(213)%grib2dsc /0/
      data paramlist(213)%grib2cat /2/
      data paramlist(213)%grib2num /0/   
      data paramlist(213)%abbrev   /'wdir    '/
      data paramlist(214)%g1tblver /0/
      data paramlist(214)%grib1val /32/
      data paramlist(214)%grib2dsc /0/
      data paramlist(214)%grib2cat /2/
      data paramlist(214)%grib2num /1/   
      data paramlist(214)%abbrev   /'wind    '/
      data paramlist(215)%g1tblver /0/
      data paramlist(215)%grib1val /33/
      data paramlist(215)%grib2dsc /0/
      data paramlist(215)%grib2cat /2/
      data paramlist(215)%grib2num /2/   
      data paramlist(215)%abbrev   /'u grd   '/
      data paramlist(216)%g1tblver /0/
      data paramlist(216)%grib1val /34/
      data paramlist(216)%grib2dsc /0/
      data paramlist(216)%grib2cat /2/
      data paramlist(216)%grib2num /3/   
      data paramlist(216)%abbrev   /'v grd   '/
      data paramlist(217)%g1tblver /0/
      data paramlist(217)%grib1val /100/
      data paramlist(217)%grib2dsc /10/
      data paramlist(217)%grib2cat /0/
      data paramlist(217)%grib2num /3/   
      data paramlist(217)%abbrev   /'htsgw   '/
      data paramlist(218)%g1tblver /0/
      data paramlist(218)%grib1val /101/
      data paramlist(218)%grib2dsc /10/
      data paramlist(218)%grib2cat /0/
      data paramlist(218)%grib2num /4/   
      data paramlist(218)%abbrev   /'wvdir   '/
      data paramlist(219)%g1tblver /0/
      data paramlist(219)%grib1val /103/
      data paramlist(219)%grib2dsc /10/
      data paramlist(219)%grib2cat /0/
      data paramlist(219)%grib2num /6/   
      data paramlist(219)%abbrev   /'wvper   '/
      data paramlist(220)%g1tblver /0/
      data paramlist(220)%grib1val /107/
      data paramlist(220)%grib2dsc /10/
      data paramlist(220)%grib2cat /0/
      data paramlist(220)%grib2num /10/  
      data paramlist(220)%abbrev   /'dirpw   '/
      data paramlist(221)%g1tblver /0/
      data paramlist(221)%grib1val /108/
      data paramlist(221)%grib2dsc /10/
      data paramlist(221)%grib2cat /0/
      data paramlist(221)%grib2num /11/  
      data paramlist(221)%abbrev   /'perpw   '/
      data paramlist(222)%g1tblver /0/
      data paramlist(222)%grib1val /109/
      data paramlist(222)%grib2dsc /10/
      data paramlist(222)%grib2cat /0/
      data paramlist(222)%grib2num /12/  
      data paramlist(222)%abbrev   /'dirsw   '/
      data paramlist(223)%g1tblver /0/
      data paramlist(223)%grib1val /110/
      data paramlist(223)%grib2dsc /10/
      data paramlist(223)%grib2cat /0/
      data paramlist(223)%grib2num /13/    
      data paramlist(223)%abbrev   /'persw   '/
! added 11/17/2005 - for wave models
      data paramlist(224)%g1tblver /129/
      data paramlist(224)%grib1val /156/
      data paramlist(224)%grib2dsc /0/
      data paramlist(224)%grib2cat /13/
      data paramlist(224)%grib2num /192/ 
      data paramlist(224)%abbrev   /'pmtc    '/
      data paramlist(225)%g1tblver /129/
      data paramlist(225)%grib1val /157/
      data paramlist(225)%grib2dsc /0/
      data paramlist(225)%grib2cat /13/
      data paramlist(225)%grib2num /193/ 
      data paramlist(225)%abbrev   /'pmtf    '/
      data paramlist(226)%g1tblver /3/
      data paramlist(226)%grib1val /11/
      data paramlist(226)%grib2dsc /0/
      data paramlist(226)%grib2cat /0/
      data paramlist(226)%grib2num /0/   
      data paramlist(226)%abbrev   /'tmp     '/
      data paramlist(227)%g1tblver /2/
      data paramlist(227)%grib1val /129/
      data paramlist(227)%grib2dsc /0/
      data paramlist(227)%grib2cat /3/
      data paramlist(227)%grib2num /198/ 
      data paramlist(227)%abbrev   /'mslma   '/
      data paramlist(228)%g1tblver /129/
      data paramlist(228)%grib1val /163/
      data paramlist(228)%grib2dsc /0/
      data paramlist(228)%grib2cat /13/
      data paramlist(228)%grib2num /194/ 
      data paramlist(228)%abbrev   /'lpmtf   '/
      data paramlist(229)%g1tblver /129/
      data paramlist(229)%grib1val /164/
      data paramlist(229)%grib2dsc /0/
      data paramlist(229)%grib2cat /13/
      data paramlist(229)%grib2num /195/ 
      data paramlist(229)%abbrev   /'lipmf   '/
! added 8/25/2006 - jfb ncar/mmm
      data paramlist(230)%g1tblver /2/
      data paramlist(230)%grib1val /189/
      data paramlist(230)%grib2dsc /0/
      data paramlist(230)%grib2cat /0/
      data paramlist(230)%grib2num /15/  
      data paramlist(230)%abbrev   /'vptmp   '/
      data paramlist(231)%g1tblver /2/
      data paramlist(231)%grib1val /178/
      data paramlist(231)%grib2dsc /0/
      data paramlist(231)%grib2cat /1/
      data paramlist(231)%grib2num /23/  
      data paramlist(231)%abbrev   /'icmr    '/
      data paramlist(232)%g1tblver /2/
      data paramlist(232)%grib1val /179/
      data paramlist(232)%grib2dsc /0/
      data paramlist(232)%grib2cat /1/
      data paramlist(232)%grib2num /32/  
      data paramlist(232)%abbrev   /'grmr    '/
      data paramlist(233)%g1tblver /2/
      data paramlist(233)%grib1val /198/
      data paramlist(233)%grib2dsc /0/
      data paramlist(233)%grib2cat /1/
      data paramlist(233)%grib2num /207/ 
      data paramlist(233)%abbrev   /'ncip    '/
      data paramlist(234)%g1tblver /2/
      data paramlist(234)%grib1val /186/
      data paramlist(234)%grib2dsc /0/
      data paramlist(234)%grib2cat /1/
      data paramlist(234)%grib2num /206/ 
      data paramlist(234)%abbrev   /'tipd    '/
      data paramlist(235)%g1tblver /2/
      data paramlist(235)%grib1val /188/
      data paramlist(235)%grib2dsc /2/
      data paramlist(235)%grib2cat /0/
      data paramlist(235)%grib2num /206/ 
      data paramlist(235)%abbrev   /'rdrip   '/
      data paramlist(236)%g1tblver /2/
      data paramlist(236)%grib1val /239/
      data paramlist(236)%grib2dsc /0/
      data paramlist(236)%grib2cat /1/
      data paramlist(236)%grib2num /208/ 
      data paramlist(236)%abbrev   /'sno t   '/
      data paramlist(237)%g1tblver /130/
      data paramlist(237)%grib1val /171/
      data paramlist(237)%grib2dsc /2/
      data paramlist(237)%grib2cat /3/
      data paramlist(237)%grib2num /193/ 
      data paramlist(237)%abbrev   /'rlyrs   '/
      data paramlist(238)%g1tblver /2/
      data paramlist(238)%grib1val /187/
      data paramlist(238)%grib2dsc /0/
      data paramlist(238)%grib2cat /17/
      data paramlist(238)%grib2num /192/ 
      data paramlist(238)%abbrev   /'ltng    '/
      data paramlist(239)%g1tblver /2/
      data paramlist(239)%grib1val /137/
      data paramlist(239)%grib2dsc /0/
      data paramlist(239)%grib2cat /3/
      data paramlist(239)%grib2num /199/ 
      data paramlist(239)%abbrev   /'tslsa   '/

      contains


         subroutine param_g1_to_g2(g1val,g1ver,g2disc,g2cat,g2num)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    param_g1_to_g2 
!   prgmmr: gilbert         org: w/np11    date: 2001-06-05
!
! abstract: this subroutine returns the corresponding grib2 discipline
!   category and number for a given grib1 parameter value and table version.
!
! program history log:
! 2000-05-11  gilbert
!
! usage:    call param_g1_to_g2(g1val,g1ver,g2disc,g2cat,g2num)
!   input argument list:
!     g1val    - grib1 parameter number for which discipline is requested
!     g1ver    - grib1 parameter table version number
!
!   output argument list:      
!     g2disc   - corresponding grib2 discipline number
!     g2cat    - corresponding grib2 category number
!     g2num    - corresponding grib2 parameter number within category g2cat
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: g1val,g1ver
           integer,intent(out) :: g2disc,g2cat,g2num

           g2disc=255
           g2cat=255
           g2num=255
! for testing
!           g2num=g1val
! for testing

           do n=1,maxparam
              if (paramlist(n)%grib1val.eq.g1val .and.
     &            paramlist(n)%g1tblver.eq.g1ver ) then
                 g2disc=paramlist(n)%grib2dsc
                 g2cat=paramlist(n)%grib2cat
                 g2num=paramlist(n)%grib2num
                 return
              endif
           enddo

           print *,'param_g1_to_g2:grib1 param ',g1val,' not found.',
     &             ' for table version ',g1ver
           return
         end subroutine

         character(len=8) function param_get_abbrev(g2disc,g2cat,g2num)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    param_get_abbrev 
!   prgmmr: gilbert         org: w/np11    date: 2002-01-04
!
! abstract: this function returns the parameter abbreviation for
!   a given grib2 discipline, category and parameter number.
!
! program history log:
! 2001-06-05  gilbert
!
! usage:     abrev=param_get_abbrev(g2disc,g2cat,g2num)
!   input argument list:
!     g2disc   - grib2 discipline number (see code table 0.0)
!     g2cat    - corresponding grib2 category number
!     g2num    - corresponding grib2 parameter number within category g2cat
!
! returns:  ascii paramter abbreviation
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: g2disc,g2cat,g2num

           param_get_abbrev='unknown '

           do n=1,maxparam
              if (paramlist(n)%grib2dsc.eq.g2disc.and.
     &             paramlist(n)%grib2cat.eq.g2cat.and.
     &             paramlist(n)%grib2num.eq.g2num) then
                 param_get_abbrev=paramlist(n)%abbrev
                 return
              endif
           enddo

!           print *,'param_get_abbrev:grib2 param ',g2disc,g2cat,
!     &              g2num,' not found.'
           return
         end function


         subroutine param_g2_to_g1(g2disc,g2cat,g2num,g1val,g1ver)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    param_g2_to_g1 
!   prgmmr: gilbert         org: w/np11    date: 2002-01-04
!
! abstract: this function returns the grib 1 parameter number for 
!   a given grib2 discipline, category and parameter number.
!
! program history log:
! 2001-06-05  gilbert
!
! usage:     call param_g2_to_g1(g2disc,g2cat,g2num,g1val,g1ver)
!   input argument list:
!     g2disc   - grib2 discipline number (see code table 0.0)
!     g2cat    - corresponding grib2 category number
!     g2num    - corresponding grib2 parameter number within category g2cat
!
!   output argument list:      
!     g1val    - grib1 parameter number for which discipline is requested
!     g1ver    - grib1 parameter table version number
!
! remarks: none
!
! attributes:
!   language: fortran 90
!   machine:  ibm sp
!
!$$$
           integer,intent(in) :: g2disc,g2cat,g2num
           integer,intent(out) :: g1val,g1ver

           g1val=255
           g1ver=255

! for testing
!           if ( g2disc.eq.255.and.g2cat.eq.255 ) then
!             g1val=g2num
!             g1ver=2
!             return
!           endif
! for testing

           do n=1,maxparam
              if (paramlist(n)%grib2dsc.eq.g2disc.and.
     &             paramlist(n)%grib2cat.eq.g2cat.and.
     &             paramlist(n)%grib2num.eq.g2num) then
                 g1val=paramlist(n)%grib1val
                 g1ver=paramlist(n)%g1tblver
                 return
              endif
           enddo

           print *,'param_g2_to_g1:grib2 param ',g2disc,g2cat,
     &              g2num,' not found.'
           return
         end subroutine


      end module

