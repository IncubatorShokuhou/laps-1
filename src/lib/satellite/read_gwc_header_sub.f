      subroutine read_gwc_header(filename,l_cell,strpix,strscnl,stppix,
     &stpscnl,reqobstm,imgtype,golatsbp,golonsbp,iwidth,idepth,goalpha,
     &istrbdy1,istrbdy2,istpbdy1,istpbdy2,bepixfc,bescnfc,fsci,idecimat,
     &iostatus)
c
c***********************************************************************
c  purpose:  read and decode the header information for gvar data coming
c  from sdhs.
c
c  method:  read in character data, short and long integer data, floating
c  point data and double precesion data.  all but the character data is 
c  in a vax specific format.  the bits and bytes of this data must be 
c  manipulated to obtain the correct representation.  the method used is
c  platform independent.
c  
c  references:
c  1.  final interface specification for the satellite data handling 
c  system communications network (sdhs-comnet), 01 february 1988
c  2.  afgwc/dons pv-wave program auto_convert.pro 
c  3.  vax fortran reference manual from the sdhs programmer library  
c***********************************************************************

c***********************************************************************
c  begin variable declaration section.  the variables are declared in
c  functional groups. the header variables are declared in two groups in 
c  the order that they are stored. temporary/working variables are then
c  declared.
c***********************************************************************

      implicit none    
      character*(*) filename
      character*3 datatyp  !data type
      character*3 datasbtp !data subtype
cc      character*2 spare1     !spare
      integer   hdrsize    !header size
      integer datastr    !data start
      integer shptptr    !shipping table ptr 
      integer lun 
      integer istatus
      integer strpix     !start pixel 
      integer strscnl    !start scanline 
      integer stppix     !stop pixel 
      integer stpscnl    !stop scanline 
      integer irdttop(2)
      integer irdtbot(2)
c     integer n 
cc    character*6 satid  !satellite identification
cc    character*2 orgdattp !data type
      integer reqobstm   !requested observation time

cc      integer decimat    !decimation resolution 
      integer   idecimat   !returned variable

cc      character*3 orgsnres !original sensor resolution

cc      integer strbdy1    !requested start boundary 1
      integer   istrbdy1   !returned variable

cc      integer strbdy2    !requested start boundary 2
      integer   istrbdy2   !returned variable

cc      integer stpbdy1    !requested stop boundary 1 
      integer   istpbdy1   !returned variable

cc      integer stpbdy2    !requested stop boundary 2
      integer   istpbdy2   !returned variable

cc      integer govismd    !goes vis mode
      real   golonsbp      !goes longitude subpoint
      real   golatsbp      !goes latitude subpoint
      real   goalpha       !goes alpha
      integer golonsbpi
      integer golatsbpi
      integer goalphai

cc      real godelta       !goes delta
cc      real gozeta        !goes zeta
cc      real goeta         !goes eta
cc      real gorho         !goes rho
cc      real gosbscn       !scanline of satellite subpoint
cc      real gosbsamp      !pixel of satellite subpoint
cc      real gogeocal      !goes geoc altitude????     
cc      integer spare2(6)  !spare bytes not currently used

      integer bepixfc    !begin pixel first cell
      integer bescnfc    !begin scaneline first cell

cc      integer width      !tracks in width of image
      integer   iwidth     !returned variable

cc      integer depth      !tracks in depth of image
      integer   idepth     !returned variable

cc      integer lstpix     !last pixel in complete image
cc      integer lstscn     !last scanline in complete imagec
cc      byte reqdc           !required data complete
cc      character*3 ortype   !original type
cc      character*3 orstype  !original subtype
cc      real*8 polarin       !polar inclination
cc      byte orproj          !original reference projection
cc      character*1 orfield  !orientation field
cc      real*8 gsublon       !goes subpoint longitude spare for dmsp
cc      character*1 orsflag  !original request size flag
cc      real*8 asclon        !ascending longitude     
cc      real*8 anommm        !anomolous mean motion
cc      real*8 relerr        !relative earth rotation rate
cc      real*8 argper        !arg perigee
cc      real*8 eccen         !eccentricity
cc      real*8 arglat        !arg latitude
cc      real*8 dttop         !dt top
cc      real*8 dtbot         !dt bottom
        integer idtbot
cc      real*8 midtime       !midpoint time

      integer idttop
      integer fsci       !first scanline of complete image
      character*2 imgtype  !image type
 
cc      integer stci       !start time of complete image
cc      byte chan            !channel 1=ir1, 2=ir2, 3=ir3, 4=water vap
cc      real rotang        !rotation angle between aries and greenwich
cc      real absmag        !absolute magnitude of craft in inertial 
                          !reference frame
cc      real saterth(3,3)  !satellite to earth tranformation matrix
cc      real pointvc       !pointing vector component
cc      real satvec(3)     !satellite position vector
cc      real rma           !role misalignment angle
cc      real pma           !pitch misalignment angle

      real*8   r8dttop
      real*8   r8dtbot
      real*8   r8time1
      real*8   r8time2
      real*8   r8timeavg

      real     rdttop
      real     rdtbot

      integer   i4time_cur
      integer   i4time_now_gg
      integer   i4time1,i4time2
      integer   i4timediff

      character ctime1*9,ctime2*9
      character cjjj1*3,cjjj2*3
      character c_hm1*4,c_hm2*4
      character cfname_cur*9
      character cfname1*9
      character cfname2*9

c     real*8    r8time

      integer   byteswp2, byteswp4
      character input(1024)*1 !just read it all and sort it later!

c      byte junk512(512)     !byte array used to 'read over' unneeded
c                            !data
c      byte junk95(95)       !byte array used to 'read over' unneeded
c                            !data
      integer   iostatus    !i/o status flag
      integer   i        !array indices

      logical   l_cell
      logical   lopen
c***********************************************************************
c  open the file that contains the gvar header and pixel data.  the 
c  record length in the direct access read is 1024 bytes (the number of
c  bytes that make up the header).
c***********************************************************************

      write(6,*)'opening ',filename
      inquire(file=filename,opened=lopen,number=lun)
      if(lopen)write(6,*)'file already open',lun
      open(unit=8, file=filename, err=100, iostat=iostatus,
     &access='direct', recl=1024,status='old')

cc     &access='direct', form='unformatted',recl=1024,status='old')

c***********************************************************************
c  read in the first 20 bytes of header information into the variables 
c  listed.  the rest of the first 512 bytes do not contain useful 
c  information for gvar data.
c***********************************************************************
      read (8,rec=1) input

c                      3       3        2       4       4      4
cc    read (8,rec=1) datatyp,datasbtp,spare1,hdrsize,datastr,shptptr

c***********************************************************************
c  because the multi-byte integer data is coming from a vax machine the 
c  bytes are not in the correct order for a unix machine.  the routine 
c  byteswp4 does the funcitonal equivalent of a byte swap.  the resulting
c  integer representation of the bit string is true for positive integers  
c  but may not be true for negative integer (this depends on whether 
c  the negative integers were stored in complement form).  
c***********************************************************************
      do i=1,3
         datatyp(i:i) = input(i)
         datasbtp(i:i) = input(i+3)
      enddo
      hdrsize = byteswp4(input(9))
      datastr = byteswp4(input(13))
      shptptr = byteswp4(input(17))      
        
      print*, datatyp, ' ',datasbtp, ' ',hdrsize, ' ',
     &   datastr, ' ', shptptr

c***********************************************************************
c  read in the bulk of the header data.  read over the first 512 bytes,
c  then read in the values for the header variables.  these values are
c  in vax binary format for character, short (2 byte) integer, long (
c  4 byte) integer, floating point, and double precision.  to obtain the 
c  proper resprentation on a unix machine the bits and bytes of this infor-
c  mation must be manipulated.  bytes must be swapped for short and long
c  integers.  the sign bit, exponent bits and fraction bits in double and
c  floating point values must be properly interpreted.  all floating 
c  point and double precisions data are first read into working integer
c  variables and then the bit strings for the different components are 
c  "picked" out and the floating point number reconstructed.  two 
c  integer variables must be used to read in the bit strings for double 
c  precision data  
c***********************************************************************
c                     1        513    517      521     525     529
c                     512        4      4       4        4      6
c      read (8,rec=1) junk512, strpix, strscnl, stppix,stpscnl,satid,
c       535      537       541       543       546      548
c       c*2       4          2        3         2       2
c     & orgdattp, reqobstm, decimat, orgsnres, strbdy1, strbdy2, 
c      550      552      554       556       560       564
c       2         2       2         4          4        4
c     & stpbdy1,stpbdy2, govismd, golonsbi, golatsbi, goalphai,
c       568         572     576     580      584      588
c         4           4     4        4         4       4
c     & godeltai, gozetai, goetai, gorhoi, gosbscni, gosbsampi, 
c        592        596      620     624      628     630    632
c          4        24          4      4       2      2       4
c     & gogeocali, spare2, bepixfc, bescnfc, width, depth, lstpix,
c        636     640   641     644      647         651        655
c          4      1      3       3          4        4           2
c     & lstscn, reqdc, ortype, orstype, polarini1, polarini2, imgtype,
c       657     658       659         663        667     668
c        1       1         4           4         1        4
c     & orproj, orfield, gsubloni1, gsubloni2, orsflag, ascloni1,
c         672      676       680       684       688
c           4       4         4         4         4
c     & ascloni2, anommmi1, anommmi2, relerri1, relerri2,   
c        692      696        700     704       708       712
c         4         4         4       4         4         4
c     & argperi1,argperi2, ecceni1, ecceni2, arglati1, arglati2,   
c         716      720    724       728    732     736       740
c         4         4       4        4      4       4          4
c     & dttopi1,dttopi2,dtboti1,dtboti2,midtimei1,midtimei2, fsci,
c        745   749   750   845     849     851      887     891
c         4     1    95     4        4      36        4       12
c     & stci, chan,junk95,rotangi,absmagi,saterthi,pointvci,satveci,
c       903    907
c        4      4
c     & rmai, pmai
      close(8)
c***********************************************************************
c  perform a 4 byte swap for long integers and a two byte swap for short
c  integers
c***********************************************************************
      strpix  = byteswp4(input(513))
      strscnl = byteswp4(input(517))
      stppix  = byteswp4(input(521))
      stpscnl = byteswp4(input(525))
c     reqobstm= byteswp4(input(537))

      idecimat  =  byteswp2(input(541))
      istrbdy1  =  byteswp2(input(546))
      istrbdy2  =  byteswp2(input(548))
      istpbdy1  =  byteswp2(input(550))
      istpbdy2  =  byteswp2(input(552))
 
      irdttop(1)   =  byteswp4(input(716))
      irdttop(2)   =  byteswp4(input(720))
      irdtbot(1)   =  byteswp4(input(724))
      irdtbot(2)   =  byteswp4(input(728))

      golonsbpi =  byteswp4(input(556))
      golatsbpi =  byteswp4(input(560))
      goalphai  =  byteswp4(input(564))

      if(l_cell)then
         call cnvtvxfl(input(556),golonsbp)
         call cnvtvxfl(input(560),golatsbp)
         call cnvtvxfl(input(564),goalpha)
         call cnvtvxfl(input(716),rdttop)
         call cnvtvxfl(input(724),rdtbot)
         idttop =int(rdttop)
         idtbot =int(rdtbot)
      else
         call convert_to_real(golonsbpi,golonsbp)
         call convert_to_real(golatsbpi,golatsbp)
         call convert_to_real(goalphai,goalpha)
         call convert_to_double(irdttop(1),irdttop(2),r8dttop)
         call convert_to_double(irdtbot(1),irdtbot(2),r8dtbot)
         idttop = int(r8dttop)
         idtbot = int(r8dtbot)
c        call convert_to_real(irdttop(1),rdttop)
c        call convert_to_real(irdtbot(1),rdtbot)
      endif

      bepixfc  =  byteswp4(input(620))
      bescnfc  =  byteswp4(input(624))
      iwidth   =  byteswp2(input(628))
      idepth   =  byteswp2(input(630))

      do i=1,2
         imgtype(i:i)=input(654+i)
      enddo

      fsci   =  byteswp4(input(740))
c
      write(ctime1,101)idttop
      write(ctime2,101)idtbot
101   format(i9)

      do i=1,9
         if(ctime1(i:i).eq.' ')ctime1(i:i)='0'
         if(ctime2(i:i).eq.' ')ctime2(i:i)='0'
      enddo

      cjjj1=ctime1(1:3)
      c_hm1=ctime1(4:7)
      cjjj2=ctime2(1:3)
      c_hm2=ctime2(4:7)
      i4time_cur=i4time_now_gg()
      call make_fnam_lp(i4time_cur,cfname_cur,istatus)
      cfname1=cfname_cur(1:2)//cjjj1//c_hm1
      cfname2=cfname_cur(1:2)//cjjj2//c_hm2
      call i4time_fname_lp(cfname1,i4time1,istatus)
      call i4time_fname_lp(cfname2,i4time2,istatus)
      r8time1=float(i4time1)/1.e6
      r8time2=float(i4time2)/1.e6 
      r8timeavg=r8time1+r8time2
      reqobstm=int((r8timeavg/2.0)*1.e6)
      
c     i4timediff=int(abs(i4time2-i4time1)/2.)
c     reqobstm=i4time1+i4timediff 
c
c   no further header variables are required - return here.
c
c***********************************************************************
c  if there is an error opening the data file print out a message and
c  end the program
c***********************************************************************

 100  if (iostatus .ne. 0)then
        print *, 'error reading file ', filename , 'io status is', 
     &  iostatus
      end if


      return 
      end
