      subroutine get_tile_list(min_lat,max_lat,min_lon,max_lon
     1,maxtiles,isbego,iwbego,itilesize,ctiletype
     1,num_tiles_needed,ctile_name_list,iwoc1,iwoc2
     1,isoc1,isoc2,istatus)

      implicit  none

      real      max_lat,min_lat
      real      max_lon,min_lon

      integer   num_tiles_needed
      integer   isbego,iwbego
      integer   itilesize
      integer   itile_ns
      integer   iwoc0
      integer   iwoc,iwoc1,iwoc2
      integer   isoc,isoc1,isoc2
      integer   isocpt,isocpo
      integer   iwocpt,iwocpo,iwocph
      integer   iwocdif
      integer   istatus
      integer   maxtiles
      integer   maxtiles_loc

      parameter (maxtiles_loc = 1000)

      character*1 ctiletype
      character*3 nstitle
      character*4 ewtitle
      character*8 ctilenamelist(maxtiles_loc)
      character*(*) ctile_name_list(maxtiles)

      double precision r8term

      istatus = 1

      r8term=(min_lat-float(isbego))/float(itilesize)
      isoc1=(int(r8term+200.)-200)*itilesize+isbego
      r8term=(min_lon-float(iwbego))/float(itilesize)
      iwoc1=(int(r8term+400.)-400)*itilesize+iwbego
      r8term=(max_lat-float(isbego))/float(itilesize)
      isoc2=(int(r8term+200.)-200)*itilesize+isbego
      r8term=(max_lon-float(iwbego))/float(itilesize)
      iwoc2=(int(r8term+400.)-400)*itilesize+iwbego

      num_tiles_needed=0

      if(iwoc1.lt.-180)then
         if(itilesize.eq.180)then
            iwocdif=iwoc2-iwoc1
            iwoc1=iwoc1+iwocdif
            iwoc2=iwoc2+abs(iwocdif)
c        else
c           iwoc1=360+iwoc1
         endif
c     elseif(iwoc1.gt.180)then
c        iwoc1=360-iwoc1
      endif
      print*,'noddy iwoc1, iwoc2 ',iwoc1,iwoc2
      do iwoc = iwoc1,iwoc2,itilesize

         iwoc0 = iwoc
         if(iwoc.lt.-180)iwoc0=360+iwoc0
c        if(iwoc.gt.+180)iwoc0=360-iwoc0

         iwocph=abs(iwoc0)/100
         iwocpt=(abs(iwoc0)-iwocph*100)/10
         iwocpo=abs(iwoc0)-iwocph*100-iwocpt*10
!
! th: 8 aug 2002 we now allow 180e longitudes (and greater). the only 
! time we want to assign w is when the longitude is less than 0.
!
         if(iwoc0.ge.0.and.iwoc0.lt.180) then
            write(ewtitle,'(3i1,a1)')iwocph,iwocpt,iwocpo,'e'
         else
            write(ewtitle,'(3i1,a1)')iwocph,iwocpt,iwocpo,'w'
         endif

         if(ewtitle(1:1).eq.' ')ewtitle(1:1)='0'
         if(ewtitle(2:2).eq.' ')ewtitle(2:2)='0'

c        ewtitle=ewtitle2//ewtitle1

         do isoc = isoc1,isoc2,itilesize

            isocpt=abs(isoc)/10
            isocpo=abs(isoc)-isocpt*10

            if(isoc.ge.0)then
              write(nstitle,'(2i1,a1)')isocpt,isocpo,'n'
            else
              write(nstitle,'(2i1,a1)')isocpt,isocpo,'s'
            endif

            num_tiles_needed=num_tiles_needed+1

            if(num_tiles_needed .gt. maxtiles_loc)then
                print*,'more tiles needed than array allocation a'
     1                ,num_tiles_needed,maxtiles_loc
                istatus = 0
                return
            endif

            ctilenamelist(num_tiles_needed)=nstitle//ewtitle

         enddo
      enddo

      if(num_tiles_needed.eq.3.and.maxtiles.eq.2)then
         num_tiles_needed = 2
      endif

      if(num_tiles_needed.le.maxtiles)then
         do itile_ns=1,num_tiles_needed
            ctile_name_list(itile_ns)=ctilenamelist(itile_ns)
         enddo
      else
         print*,'more tiles than array allocation b'
     1         ,num_tiles_needed,maxtiles
         istatus = 0
      endif

      return
      end
