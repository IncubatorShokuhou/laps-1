      subroutine loadlapsdata(nx_l,ny_l,maxchannels,n_lvd_fields_max,
     &                        ntm,c_type,r_image_status,csatid,
     &                        image_vis,
     &                        image_39,image_67,
     &                        image_ir,image_12,
     &                        var_lvd,c_lvd,units_lvd,
     &                        nlf,laps_data,istatus)
c
c routine loads image satellite data that has already been
c navigated to the laps domain. currently this is only from afwa.
c
      implicit none

      integer   nx_l,ny_l
      integer   maxchannels
      integer   n_lvd_fields_max
      integer   i
      integer   ispec
      integer   ntm
      integer   nlf
      integer   istatus
      integer   lstatus

      real      r_image_status(maxchannels)
      real      laps_data(nx_l,ny_l,n_lvd_fields_max)
      real      image_vis(nx_l,ny_l)
      real      image_39(nx_l,ny_l)
      real      image_67(nx_l,ny_l)
      real      image_ir(nx_l,ny_l)
      real      image_12(nx_l,ny_l)

      character c_lvd(n_lvd_fields_max)*125
      character c_type(maxchannels)*3
      character var_lvd(n_lvd_fields_max)*3
      character units_lvd(n_lvd_fields_max)*10
      character csatid*(*)

      do i=1,ntm

         istatus=0

         call lvd_file_specifier(c_type(i),ispec,lstatus)
         if(lstatus.ne.0)goto 900

         if(ispec.eq.1)then
           print*,'not ready to deal with gms visible data'

c           call some-routine-to-deal-with-gms-visible-data
c           if(r_image_status(i).lt.0.3333)then
c           nlf=nlf+1
c           call move(visraw,laps_data(1,1,nlf),nx_l,ny_l)
c           var_lvd(nlf) = 'svs'       ! satellite, visible
c           c_lvd(nlf)=csatid//' (visible) satellite - raw'
c           units_lvd(nlf) = 'counts'
c and more
c
c only goes has 3.9u
c
         elseif(ispec.eq.2)then

            if(csatid.ne.'gmssat'.or.csatid.ne.'metsat')then
            if(r_image_status(i).le.0.3333)then
               nlf=nlf+1
               call move(image_39,laps_data(1,1,nlf),nx_l,ny_l)
               var_lvd(nlf) = 's3a'       ! satellite, , averaged
               c_lvd(nlf)=csatid//' (3.9u) ir b-temps'
               nlf=nlf+1
               call move(image_39,laps_data(1,1,nlf),nx_l,ny_l)
               var_lvd(nlf)  = 's3c'       ! satellite, , filtered
               c_lvd(nlf)=csatid//' (3.9u) ir b-temps'
            else
               write(6,*)'39u image not processed: missing data'
            endif
            endif

         elseif(ispec.eq.3)then

         if(r_image_status(i).le.0.3333)then
            nlf=nlf+1
            call move(image_67,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf) = 's4a'       ! satellite, averaged
            c_lvd(nlf)=csatid//' (6.7u) ir b-temps'
            nlf=nlf+1
            call move(image_67,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf) = 's4c'       ! satellite, filtered
            c_lvd(nlf)=csatid//' (6.7u) ir b-temps'
         else
            write(6,*)'wv image not processed: missing data'
         endif

         elseif(ispec.eq.4)then

         if(r_image_status(i).lt.0.3333)then
            nlf=nlf+1
            call move(image_ir,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf)  = 's8a'       ! satellite, channel-4, averaged
            c_lvd(nlf)=csatid//' (11.0u) ir b-temps'
            nlf=nlf+1
            call move(image_ir,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf)='s8w'       ! satellite, channel-4, warm pixel
            c_lvd(nlf)=csatid//' (11.0u) ir b-temps'
            nlf=nlf+1
            call move(image_ir,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf)='s8c'       ! satellite, channel-4, warm pixel
            c_lvd(nlf)=csatid//' (11.0u) ir b-temps'
         else
            write(6,*)'ir image not processed: missing ir data'
         endif

         elseif(ispec.eq.5)then

         if(r_image_status(i).lt.0.3333)then
            nlf=nlf+1
            call move(image_12,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf) = 'sca'       ! satellite, averaged
            c_lvd(nlf)=csatid//' (12.0u) ir b-temps'
            nlf=nlf+1
            call move(image_12,laps_data(1,1,nlf),nx_l,ny_l)
            var_lvd(nlf) = 'scc'       ! satellite, averaged
            c_lvd(nlf)=csatid//' (12.0u) ir b-temps'
         else
            write(6,*)'12u image not processed: missing data'
         endif

         endif
         goto 1000

900      print*,'error returned from lvd_file_specifier'
         istatus=1

1000  enddo

      return
      end
