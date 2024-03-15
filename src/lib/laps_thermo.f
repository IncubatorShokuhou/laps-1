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

        subroutine put_stability(
     1           i4time_needed                   ! input
     1          ,nx_l,ny_l,nz_l                  ! input
     1          ,heights_3d                      ! input
     1          ,lat,lon,topo                    ! input
     1          ,laps_cycle_time                 ! input
     1          ,temp_3d                         ! input
     1          ,rh_3d_pct                       ! input
     1          ,temp_sfc_k                      ! input
     1          ,pres_sfc_pa                     ! input
     1          ,twet_snow                       ! input
     1          ,td_3d_k                         ! output
     1          ,istatus)                        ! output

cdoc    calculate and write out set of 2-d stability grids

!       arrays passed in
        real temp_3d(nx_l,ny_l,nz_l)
        real rh_3d_pct(nx_l,ny_l,nz_l)
        real heights_3d(nx_l,ny_l,nz_l)
        real temp_sfc_k(nx_l,ny_l)
        real pres_sfc_pa(nx_l,ny_l)
        real lat(nx_l,ny_l)
        real lon(nx_l,ny_l)
        real topo(nx_l,ny_l)

!       output
        real td_3d_k(nx_l,ny_l,nz_l)

!       local declarations for stability 
        real t_sfc_f(nx_l,ny_l)
        real td_sfc_k(nx_l,ny_l)
        real td_sfc_f(nx_l,ny_l)

        real pbe_2d(nx_l,ny_l)
        real nbe_2d(nx_l,ny_l)
        real si_2d(nx_l,ny_l)
        real tt_2d(nx_l,ny_l)
        real k_2d(nx_l,ny_l)
        real lcl_2d(nx_l,ny_l)
        real wb0_2d(nx_l,ny_l)
        real wb1_2d(nx_l,ny_l)

        real li(nx_l,ny_l)

        real pres_sfc_mb(nx_l,ny_l)
        real pres_3d(nx_l,ny_l,nz_l)

        integer nfields
        parameter (nfields=8)

        character*10 units_2d_a(nfields)
        character*125 comment_2d_a(nfields)
        character*3 var_2d_a(nfields)
        real out_multi_2d(nx_l,ny_l,nfields)

        character*31 ext

        character*3 var_2d
        character*10  units_2d
        character*125 comment_2d

        real k_to_f

!       read in surface dewpoint data
        var_2d = 'td'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,nx_l,ny_l       
     1                      ,td_sfc_k,0,istatus)
        if(istatus .ne. 1)then
            write(6,*)' laps sfc dewpoint not available'
            write(6,*)' abort put_stability routine'
            return
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

!       call get_pres_3d(i4time_needed,nx_l,ny_l,nz_l,pres_3d,istatus)       
!       if(istatus .ne. 1)return

!       convert rh to td
        do k = 1,nz_l
        do j = 1,ny_l
        do i = 1,nx_l
            if(rh_3d_pct(i,j,k) .ge. 0. .and.
     1         rh_3d_pct(i,j,k) .le. 100.       )then    ! rh in valid range

                t_c = temp_3d(i,j,k)-273.15
                td_c = dwpt(t_c,rh_3d_pct(i,j,k))
                td_3d_k(i,j,k) = td_c + 273.15

!               ew = pres_3d(i,j,k)/100. * sh_3d(i,j,k)
!               td_c = dpt(ew)                        ! check valid input range

!               td_c = make_td(pres_3d(i,j,k)/100.
!    1                        ,temp_3d(i,j,k)-273.15
!    1                        ,sh_3d(i,j,k)*1000.
!    1                        ,-100.)

            elseif(rh_3d_pct(i,j,k) .gt. 100. .and.
     1             rh_3d_pct(i,j,k) .le. 101.       )then ! rh slightly high

                td_3d_k(i,j,k) = temp_3d(i,j,k)
                write(6,*)' warning: rh out of bounds',rh_3d_pct(i,j,k)       
     1                   ,' at ',i,j,k,' setting td = t'

            else ! invalid value of rh to pass into dwpt
                td_3d_k(i,j,k) = r_missing_data
                write(6,*)' error: rh out of bounds',rh_3d_pct(i,j,k)       
     1                   ,' at ',i,j,k
                istatus = 0
                return

            endif

        enddo ! i
        enddo ! j
        enddo ! k            

        call laps_be(nx_l,ny_l,nz_l,twet_snow,lat,lon
     1              ,temp_sfc_k,td_sfc_k,pres_sfc_pa
     1              ,temp_3d,td_3d_k,heights_3d,topo,blayr_thk_pa
     1              ,pbe_2d,nbe_2d,si_2d,tt_2d,k_2d,lcl_2d,wb0_2d
     1              ,wb1_2d,r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' warning, bad istatus returned from laps_be'
            return
        endif

!       fill pres_sfc_mb
        call move(pres_sfc_pa,pres_sfc_mb,nx_l,ny_l)
        call multcon(pres_sfc_mb,.01,nx_l,ny_l)

!       convert t, td to f
        do i = 1,nx_l
        do j = 1,ny_l
            t_sfc_f(i,j)  = k_to_f(temp_sfc_k(i,j))
            td_sfc_f(i,j) = k_to_f(td_sfc_k(i,j))
        enddo ! j
        enddo ! i

        flag = 0.0
        call li_laps(t_sfc_f,td_sfc_f,pres_sfc_mb
     1              ,i4time_needed,nx_l,ny_l,li,flag,istatus)

!       call move
        call move(pbe_2d,out_multi_2d(1,1,1),nx_l,ny_l)
        call move(nbe_2d,out_multi_2d(1,1,2),nx_l,ny_l)
        call move(    li,out_multi_2d(1,1,3),nx_l,ny_l)
        call move( si_2d,out_multi_2d(1,1,4),nx_l,ny_l)
        call move( tt_2d,out_multi_2d(1,1,5),nx_l,ny_l)
        call move(  k_2d,out_multi_2d(1,1,6),nx_l,ny_l)
        call move(lcl_2d,out_multi_2d(1,1,7),nx_l,ny_l)
        call move(wb0_2d,out_multi_2d(1,1,8),nx_l,ny_l)

!       add var arrays
        ext = 'lst'

        var_2d_a(1) = 'pbe'
        var_2d_a(2) = 'nbe'
        var_2d_a(3) = 'li'
        var_2d_a(4) = 'si'
        var_2d_a(5) = 'tt'
        var_2d_a(6) = 'k'
        var_2d_a(7) = 'lcl'
        var_2d_a(8) = 'wb0'

        units_2d_a(1) = 'j/kg'
        units_2d_a(2) = 'j/kg'
        units_2d_a(3) = 'k'
        units_2d_a(4) = 'k'
        units_2d_a(5) = 'k'
        units_2d_a(6) = 'k'
        units_2d_a(7) = 'm'
        units_2d_a(8) = 'm'

        comment_2d_a(1) = 'cape'
        comment_2d_a(2) = 'cin'
        comment_2d_a(3) = 'lifted_index'
        comment_2d_a(4) = 'showalter_index'
        comment_2d_a(5) = 'total_totals'
        comment_2d_a(6) = 'k_index'
        comment_2d_a(7) = 'lcl'
        comment_2d_a(8) = 'wet_bulb_zero'

        call put_laps_multi_2d(i4time_needed,ext,var_2d_a,units_2d_a
     1                        ,comment_2d_a,out_multi_2d,nx_l,ny_l,8    
     1                        ,istatus)
        if(istatus .ne. 1)then
            write(6,*)' lst output error'
        else
            write(6,*)' successfully wrote lst'
        endif

        return

        end


        subroutine laps_be(ni,nj,nk,twet_snow,lat,lon
     1        ,t_sfc_k,td_sfc_k,p_sfc_pa,t_3d_k,td_3d_k,ht_3d_m,topo       
     1        ,blayr_thk_pa,pbe_2d,nbe_2d,si_2d,tt_2d,k_2d,lcl_2d,wb0_2d
     1        ,wb1_2d,r_missing_data,istatus)

!       1991    steve albers
cdoc    returns 2-d pbe and nbe in joules, parcel is lifted from lowest level
!                                                                    i.e. sfc

        real lat(ni,nj)
        real lon(ni,nj)
        real t_sfc_k(ni,nj)
        real td_sfc_k(ni,nj)
        real p_sfc_pa(ni,nj)
        real topo(ni,nj)
        real t_3d_k(ni,nj,nk)
        real td_3d_k(ni,nj,nk)
        real ht_3d_m(ni,nj,nk)
        real p_1d_pa(nk)
        real p_1d_mb(nk)                   ! local
        real pres_3d(ni,nj,nk)             ! local
        real pbe_2d(ni,nj)
        real nbe_2d(ni,nj)

        real si_2d(ni,nj)
        real tt_2d(ni,nj)
        real k_2d(ni,nj)
        real lcl_2d(ni,nj)
        real wb0_2d(ni,nj)
        real wb1_2d(ni,nj)
        
        include 'lapsparms.for'
        integer mxl
        parameter (mxl=max_lvls+1) ! number of 3d levels plus the sfc level

        common/indx/ p(mxl),t(mxl),td(mxl),ht(mxl),pbecr(20,4)
     1  ,tdfcr(20,2),vel(20)
     1  ,temdif(mxl),partem(mxl),pbe(mxl)
     #  ,dd85,ff85,dd50,ff50
        real lcl,li,k_index

!       initialize pbe array
        do i = 1,mxl
            pbe(i) = 0.
        enddo ! i

        call get_systime_i4(i4time,istatus)
        if(istatus .ne. 1)stop

        call get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)
        if(istatus .ne. 1)stop

        if(nj .gt. 600)then
            jint = 40
        else
            jint = 20
        endif

        do j = 1,nj
        do i = 1,ni
c       write(6,*)' i = ',i

            if(j .eq. (j/jint)*jint .and. i .eq. ni/2)then
                idebug = 2 ! 1
                write(6,*)
                write(6,*)
     1          ' --- stability index debugging at gridpoint ---',i,j
            else 
                idebug = 0
            endif
        
            do k = 1,nk
                p_1d_pa(k) = pres_3d(i,j,k)
                p_1d_mb(k) = p_1d_pa(k) / 100.
            enddo ! k

            do k = 1,nk
                if(p_1d_pa(k) .lt. p_sfc_pa(i,j))then ! first level above sfc
                    n_first_level = k
                    goto100
                endif
            enddo

100         p(1) = p_sfc_pa(i,j) / 100.       ! pa to mb
            t(1) = t_sfc_k(i,j) - 273.15      ! k to c
            td(1) = td_sfc_k(i,j) - 273.15    ! k to c
            ht(1) = topo(i,j)

            n = 1

            do k = n_first_level,nk
                n = n + 1
                p(n) = p_1d_mb(k)
                t(n) = t_3d_k(i,j,k)  - 273.15 ! k to c
                td(n)= td_3d_k(i,j,k) - 273.15 ! k to c
!               td(n)= t(n)
                ht(n)= ht_3d_m(i,j,k)
            enddo ! k

            nlevel = n

            io = 0

            call sindx(nlevel,li,si,bli,tt,sweat
     1                ,twet_snow                                        ! i
     1                ,hwb0,hwb_snow                                    ! o
     1                ,plcl,lcl,ccl
     1                ,tconv,io
     1                ,icp,ict,k_index,tmax,pbeneg,pbepos,t500,pbli
     1                ,velneg,water,ihour,idebug,istatus)
            if(istatus .ne. 1)then
                write(6,*)' warning: bad istatus returned from sindx'
                return
            endif

            pbe_2d(i,j) = pbepos
            nbe_2d(i,j) = pbeneg
          
            si_2d(i,j) = si
            tt_2d(i,j) = tt
            k_2d(i,j)  = k_index

            if(lcl .ne. r_missing_data)then
                lcl_2d(i,j) = lcl                        ! m msl
            else
                lcl_2d(i,j) = r_missing_data
            endif

            if(hwb0 .ne. r_missing_data)then
                wb0_2d(i,j) = hwb0*304.8006 + topo(i,j)      ! kft agl to m msl
            else
                wb0_2d(i,j) = r_missing_data
            endif

            if(hwb_snow .ne. r_missing_data)then
                wb1_2d(i,j) = hwb_snow*304.8006 + topo(i,j)  ! kft agl to m msl
            else
                wb1_2d(i,j) = r_missing_data
            endif

            iwarn = 0
            if(nanf(wb0_2d(i,j)) .eq. 1)then
                write(6,*)' warning: hwb0 nan ',i,j
                iwarn = 1
            endif

            if(idebug .ge. 1 .or. iwarn .eq. 1)then

                write(6,*)
                write(6,*)' indices at:',i,j,lat(i,j),lon(i,j)

                write(6,*)' n    p         t         td'
     1                   ,'         ht       tdif       be'
                do n = 1,nlevel
                    write(6,301)n,p(n),t(n),td(n),ht(n),temdif(n),pbe(n)
301                 format(1x,i2,f8.1,f10.2,f10.2,f10.0,f10.2,f10.1)
                enddo

                sweat=0.0
!               if(iflag1*iflag2.eq.0)sweat=0.0

                write(6,420)pbeneg,pbepos
 420            format(' pbeneg',f8.1,'               pbepos',f8.1)

                write(6,430)dd85,ff85,dd50,ff50,t500
 430            format(' 850mb wind',2f5.0,'         500mb wind',2f5.0
     #         ,'    500mb temp= ',f6.1)

!               if(io2.eq.1)goto1000

                write(6,60)li,si,bli
 60             format(' li= ',f5.1,20x,'si= ',f5.1,15x,'bli= ',f5.1)

                write(6,61,err=161)tt,sweat,hwb0
 61             format(' total totals=',f6.1,10x,'sweat= ',f7.1
     1                ,10x,'w. bulb zero=',f5.1,' kft agl')

 161            itconv=int(tconv+0.5)

                write(6,62)lcl*.003281,ccl,itconv
 62             format(' lcl= ',f5.1,' kft msl',11x,'ccl=',f5.1
     1                ,' kft agl',7x,'convective temp=',i4,' deg f')

                itmax=nint(tmax)

                write(6,63,err=163)k_index,itmax,water
 63             format(' k index =',f7.1,13x,'tmax =',i4,' f',12x
     1                ,'precip. water=',f5.2,' in.')
 163            write(6,*)
c
                if(ict+icp.gt.0.or.io.ge.2)write(6,71)p(1),plcl
 71             format(' energy analysis - lifting a parcel from'
     1                ,f6.0,'mb','   lcl=',f7.1,' mb')
c
                do 600 ii=1,icp
 600            write(6,64)pbecr(ii,1),pbecr(ii,2),vel(ii),pbecr(ii,4)
 64             format(' energy at',f6.0,'mb=',f10.2
     1                ,' joules/kg       velocity=',f6.2,'m/s   ht='
     1                ,f7.0,' m')
c
                write(6,*)
                do 700 ii=1,ict
 700            write(6,5)tdfcr(ii,1),tdfcr(ii,2)
 5              format(5x,'environmental minus parcel temperature at'
     1                ,f6.0,'mb  =',f6.1)
c
                if(bli.le.0.0)write(6,65)pbeneg,velneg
                if(bli.le.0.0)write(6,66)pbepos
 65             format('     cap strength =',f10.1,' joules/kg.'
     +                ,' velocity needed',f5.1,'m/s')
 66             format('     positive area=',f10.1,' joules/kg.')
                goto2000
c
c  alternative outputting
 1000           continue
                write(6,1001)sta,t500,velneg,pbepos
 1001           format(4x,a3/f5.1,f7.1,f7.1)

2000        endif

        enddo ! i
        enddo ! j
c
 9999   continue

        return
        end

!
        subroutine sindx(nlevel,li,si,bli,tt,sweat
     1   ,twet_snow,hwb0,hwb_snow
     1   ,plcl_pbe,lcl_pbe_msl,ccl  
     1   ,tconv,io,icp,ict,k,tmax,pbeneg,pbepos,tman50
     1   ,pbli,velneg,water,ihour,idebug,istatus)

cdoc    calculate a variety of stability indices from an input sounding

!       1991    steve albers
!       1999    steve albers    adding more indices to active output

!       note that a surface parcel is currently used for the indices

        include 'lapsparms.for'
        integer mxl
        parameter (mxl=max_lvls+1) ! number of 3d levels plus the sfc level

        dimension q(mxl),w(mxl),wb(mxl)
        common/indx/ p(mxl),t(mxl),td(mxl),ht(mxl),pbecr(20,4)
     1              ,tdfcr(20,2),vel(20),temdif(mxl),partem(mxl)
     1              ,pbe(mxl),dd85,ff85,dd50,ff50
        real li,k,lcl_agl,lcl_pbe_msl
!       esl(x)=6.1078+x*(.443652+x*(.014289+x*(2.65065e-4+x*
!    1 (3.03124e-6+x*(2.034081e-8+x*(6.13682e-11))))))
!       tdew(e)=237.7/((7.5/alog10(e/6.11))-1.)

        esl(x)=6.1121*exp(17.67*x/(x+243.5)) 

        logical l_large_domain

        data epsiln/.6220/,g/9.80665/
        data blthck/50.0/
        data rpd/.0174532925063/
 1      format('    mixed parcel is:',f10.1,'mb',f15.3,'c       '
     1          ,'mixing ratio=',f10.7)

c       write(6,15)
 15     format('          pressure     temp.    dew pt.     '
     1          ,'q      wet bulb')

        iout=io

        l_large_domain = .false.

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            return
        endif

!       fill wb array (imported from old code)
 	do 100 n=1,nlevel                                            
            if(td(n).eq.-99.0)then                                         
                q(n)=-99.                                                      
                wb(n)=-99.                                                     
  	        if(io.ge.2)write(6,3)n,p(n),t(n),td(n),q(n),wb(n)    
                ifl=1                                                
                td(n)=t(n)-15.0                                      
                if(t(n).lt.-40.)td(n)=t(n)                           
            endif                                                           
                                                                          
 	    q(n)=(esl(td(n))*epsiln)/p(n)                             
 	    q(n)=max(q(n),0.0001)                                    
 	    w(n)=q(n)/(1.-q(n))                                      
 	    iout=io                                                  
 	    if(n.le.1)iout=io                                        

 	    wb(n)=wtblb(p(n),t(n),w(n),0,istatus)                            
            if(istatus .ne. 1)then
                write(6,*)' bad istatus return from wtblb: '
                write(6,*)'p,t,td,es,epsiln,q,w = '
                write(6,*)p(n),t(n),td(n),esl(td(n)),epsiln,q(n),w(n)
!               write(6,*)'es calc:',
!    1                         6.1121*exp(17.67*td(n)/(td(n)+243.5))
                return
            endif

 	    if(wb(n).lt.td(n))wb(n)=td(n)                            
 	    if(wb(n).gt.t(n)) wb(n)=t(n)                             
 	    if(io.ge.2.and.ifl.eq.0)
     1         write(6,3)n,p(n),t(n),td(n),q(n),wb(n)     
            ifl=0                                                            
 3	    format(' lvl(',i2,')',3f10.1,f11.6,f10.2)                      
 100	continue                                                             

!       here is the new section imported (from old code) for tt,si,k
c                                                                         
c  calculate lifted index                                                 
 	call itplv(p,t,nlevel,500.,tman50,io,istatus)                         
!wni modified to check istatus of itplv.  when running laps
!wni for east asia, some of the mountains extended above the 500mb
!wni level, which cause a complete abort of put stability.  now,
!wni we just fill indices dependent on 500mb to missing if this happens.
!wni brent shaw, wni, dec 2006
!	wreq50=esl(tp500)/500.                                                   
!	if(wreq50.lt.wmean)goto40                                               

c                                                                         
c  calculate showalter index                                              
        si=r_missing_data
        tt=r_missing_data
        sweat=r_missing_data
        k=r_missing_data
        if (istatus .eq. 1) then  !wni .. 500mb below ground
          if(p(1).ge.850.0)then
            call itplv(p,t ,nlevel,850.,tman85,io,istatus)
 	    call itplv(p,td,nlevel,850.,tdmn85,io,istatus)
 	    thetae=thae(tman85,tdmn85,850.)
 	    call msad5(tp500,500.,thetae,25.,20.,slope,i1,i2,ia,0
     1                ,istatus)    
            if(istatus .ne. 1)then
                write(6,*)' error, skipping showlater in sindx'
            else
!	        if(wreq50.lt.wmean)goto50                                 
!	        rh=wmean/wreq50                                           
 50	        si=tman50-tp500                                           
            endif ! valid results from msad5

c                                                                         
c  calculate total totals and sweat and k indicies                        
 	    call itplv(p,t,nlevel,700.,tman70,io,istatus)        
 	    call itplv(p,td,nlevel,700.,tdmn70,io,istatus)       
 	    tt=tman85+tdmn85-2.*tman50                   
 	    a=max(tdmn85,0.)                             
 	    b=max(tt-49.,0.)                             
 	    c=0.                                         
 	    if(ff85.lt.15.or.ff50.lt.15.)   goto500      
 	    if(dd85.gt.dd50)                goto500      
 	    if(dd85.lt.130..or.dd85.gt.250.)goto500      
 	    if(dd50.lt.210..or.dd50.gt.310.)goto500      
 	    c=sin((dd50-dd85)*rpd)+.2                    
 500	    sweat=12.*a+20.*b+2.*ff85+ff50+125.*c        
 	    k=tman85-tman50+tdmn85-tman70+tdmn70         

          endif
        else  !wni 
          print *,"no 500mb level, some indices set to missing" !wni
          istatus = 1 ! wni
        endif ! wni

c                                                                         
c  calculate wet bulb zero level                                          
        twet0 = 0.       ! deg c
        call wtblb_lvl(twet0,p,t,q,wb,mxl,nlevel,pwb0,hwb0)

        call wtblb_lvl(twet_snow,p,t,q,wb,mxl,nlevel,pwb_snow,hwb_snow)       

!       calculate theta(e) based on sfc parcel
        if(l_large_domain)then
            thetae=oe_fast(t(1),td(1),p(1)) + 273.15
        else
            thetae=oe(t(1),td(1),p(1)) + 273.15
        endif

        if(idebug .ge. 2)then
            if(l_large_domain)then
                thetae2=oe(t(1),td(1),p(1)) + 273.15
                write(6,510)p(1),t(1),td(1),thetae,thetae2
 510            format(' p/t/td/thetae',5f8.2)
                if(abs(thetae-thetae2) .gt. 1.0)then
                    write(6,*)' error: large difference in theta values'
                    stop
                endif
            else
                write(6,510)p(1),t(1),td(1),thetae
            endif
        endif

!       calculate "fast" lcl based on sfc parcel
        call lcl_fast(p(1),t(1),td(1),lcl_agl,tlcl_pbe,plcl_pbe)

!       this lcl is for the surface parcel passed in
        lcl_pbe_msl = lcl_agl + ht(1)

!       calculate cape/cin based on sfc parcel
        call potbe(q,nlevel,p(1),t(1),w(1),plcl_pbe
     1   ,tlcl_pbe,lcl_pbe_msl,thetae,icp,ict,io,pbeneg,pbepos
     1   ,velneg,idebug)
c
        do 600 i=1,icp
            vel(i)=sqrt(abs(2.*pbecr(i,2)))
            if(pbecr(i,2).lt.0.)vel(i)=-vel(i)
c           write(6,4)pbecr(i,1),pbecr(i,2),vel(i)
 4          format(' energy at',f6.0,'mb='
     1            ,f15.2,'joules       velocity=',f6.1,'ms')
 600    continue

c       do 700 i=1,ict
c           write(6,5)tdfcr(i,1),tdfcr(i,2)
c700    continue

 5      format(' environmental - parcel temperature at',f7.0
     1          ,'mb=',f8.2,' deg c')

        return

        end
c
c
c
c        block data
c        common/indx/ p(mxl),t(mxl),td(mxl),ht(mxl),pbecr(20,4),tdfcr(20,2)
c     1              ,vel(20),temdif(mxl),partem(mxl),pbe(mxl)
c     #              ,dd85,ff85,dd50,ff50
c        data pbe/mxl*0./
c        end
c
c
c
c
c
c
c
        subroutine potbe(q,nlevel,pmean,tmean,wmean,plcl
     #   ,tlcl,lcl,blthte,icp,ict,io,pbeneg,pos_area_max
     #   ,velneg,idebug)

cdoc    calculate a pbe/lcl related indices from an input sounding

!       steve albers 1991

        include 'lapsparms.for'
        integer mxl
        parameter (mxl=max_lvls+1) ! number of 3d levels plus the sfc level

        common/indx/ p(mxl),t(mxl),td(mxl),ht(mxl),pbecr(20,4)
     1              ,tdfcr(20,2),vel(20)
     1              ,temdif(mxl),partem(mxl),pbe(mxl)
     #              ,dd85,ff85,dd50,ff50
        real lcl,nbe_min
        dimension q(mxl),deltah(mxl)
        real gamma
        parameter (gamma = .009760) ! dry adiabatic lapse rate deg/m
        data g/9.80665/

        if(idebug .ge. 2)then
            write(6,*)
            write(6,*)' subroutine potbe: blthte = ',blthte
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error reading r_missing_data: stop'
            stop
        endif

        nbe_min=0.
        pos_area_max=0.

        partem(1) = t(1)
        temdif(1) = 0.
        pbe(1) = 0.

c       write(6,606)pmean,tmean,wmean,plcl,tlcl,lcl
c 606   format(' pmean,tmean,wmean,plcl,tlcl,lcl',2f10.2,f10.5,2f10.2
c     +   ,f5.1)

        do 800 n=2,nlevel
            deltap=p(n-1)-p(n)
            n_low=n-1

            deltah(n) = ht(n) - ht(n-1)

            if(plcl .lt. p(n-1))then ! lower level is below lcl
                if(plcl .lt. p(n))then ! upper level is below lcl
                    if(idebug .ge. 2)write(6,*)' dry case'
                    partem(n)=partem(n-1)-gamma*deltah(n)
                    temdif(n)=partem(n)-t(n)
                    pbe(n)=pbe(n-1)+g*(temdif(n)+temdif(n-1))/(t(n)
     #                          +t(n-1)+546.30)*deltah(n)

                else ! upper level is above lcl
c                   bracketing case around lcl - dry adiabatic part
                    if(idebug .ge. 2)write(6,*)' dry adiabatic part'

                    delta_ht_dry = lcl - ht(n-1)

                    if(idebug .ge. 2)write(6,307)tlcl
 307                format(' parcel temp at lcl= ',f10.3)

                    call itplv(p,t,nlevel,plcl,sntlcl,io,istatus)

                    t_dif_lcl=tlcl-sntlcl

                    pbe_dry=g*(t_dif_lcl+temdif(n-1))/
     1                        (sntlcl+t(n-1)+546.30) * delta_ht_dry

 951                format(' pos_area_max,psi,neg,nng,pbecr',5f9.1,i3)

!d                  write(6,777)n,p(n),tlcl,sntlcl,t_dif_lcl
!d    #                    ,temdif(n-1),delta_ht_dry,ht(n),pbe(n-1)+pbe_dry

c                   moist adiabatic part
                    if(idebug .ge. 2)write(6,*)' moist adiabatic part'
                    delta_ht_wet=deltah(n)-delta_ht_dry

                    partem(n) = tmlaps_fast(blthte,p(n))
                    temdif(n)=partem(n)-t(n)

                    pbe_wet = g*(temdif(n)+t_dif_lcl)/
     1                          (t(n)+sntlcl+546.30) * delta_ht_wet

                    pbe(n)=pbe(n-1) + pbe_dry + pbe_wet

                endif ! upper level below lcl (dry or bracket)

            else ! lower level is above lcl
                partem(n) = tmlaps_fast(blthte,p(n))
                if(idebug .ge. 2)then
                  write(6,*)' moist case',blthte,p(n),partem(n)
                endif

          !     add layer
                temdif(n)=partem(n)-t(n)
                pbe(n)=pbe(n-1)+g*(temdif(n)+temdif(n-1))/(t(n)
     #                          +t(n-1)+546.30)*deltah(n)

            endif


            if(idebug .ge. 2)then
                write(6,777)n,p(n),partem(n),t(n)
     #                ,temdif(n),temdif(n-1),deltah(n),ht(n),pbe(n)
!    #                ,blthte
            endif
 777        format(' pbe data',i3,f6.1,6f8.2,f12.3,f8.2)

 800    continue
c
c  flag significant levels for energy
        icp=0
        ict=0
        n1 = 2
c       write(6,382)n1,nlevel
 382    format(' n1,nlevel',2i3)
c
c  determine energy extrema - neutral buoyancy
        do 1100 n=n1,nlevel
            if(idebug .ge. 2)write(6,940)n
 940        format(' looking for neutral buoyancy - energy extremum, lev
     1el',i3)

            if((temdif(n)*temdif(n-1)).lt.0.)then
                ict=ict+1
                slope=(temdif(n)-temdif(n-1))/alog((p(n)/p(n-1)))
                frac=temdif(n-1)/(temdif(n-1)-temdif(n))
                tdfcr(ict,1)=exp(alog(p(n-1))-temdif(n-1)/slope)
                tdfcr(ict,2)=0.
                if(idebug .ge. 2)then
                    write(6,948)ict,n,temdif(n-1),temdif(n),slope,frac
     #                         ,tdfcr(ict,1)
                endif
 948            format(' neut buoy',2i3,5f12.5)
                icp=icp+1
                tmid=(1.-frac)*t(n-1)+frac*t(n)
                pbecr(icp,1)=tdfcr(ict,1)
                pbecr(icp,2)=pbe(n-1)+g*temdif(n-1)/(t(n-1)+tmid+546.3)
     &             *deltah(n)*frac
c
                if(p(n) .ge. 500.)then
                    nbe_min=amin1(pbecr(icp,2),nbe_min)
                endif

                pos_area=pbecr(icp,2)-nbe_min
                pos_area_max=max(pos_area,pos_area_max)

                if(idebug .ge. 2)write(6,951)pos_area_max
     1                       ,pos_area,pbeneg,nbe_min,pbecr(icp,2),icp

                pbecr(icp,3)=pbecr(icp,2)-pbecr(max(1,icp-1),2)
                pbecr(icp,4)=ht(n-1)+frac*deltah(n)
c
                if(idebug .ge. 2)then
                  write(6,949)icp,pbecr(icp,1),pbecr(icp,2),pbecr(icp,3)
     &          ,pbecr(icp,4),pos_area_max,nbe_min
 949              format(' pbecr',i3,4f10.3,2f10.2)
                endif

            endif

            goto1100
c
c  find buoyancy maxima
c           if(n.lt.nlevel)then
c               if(temdif(n) .gt. temdif(n-1)
c       1                               .and. temdif(n) .gt. temdif(n+1))then
c                   if(abs(temdif(n)).ge.abs(temdif(n-1)).or.
c     &                 temdif(n)*temdif(n-1).le.0)then
c                       ict=ict+1
c                       tdfcr(ict,1)=p(n)
c                       tdfcr(ict,2)=-temdif(n)
c                       pbecr(1,3)=pbecr(1,2)
c                       write(6,999)ict,temdif(n)
c 999                   format(' buoyancy max, ict temdif(n)',i3,f5.1)
c                    endif
c                endif
c            endif
c
c  look for an lfc (zero energy)
c 1000      continue

c  find temperature difference as a function of height then
c  find pbe as a function of height (above sig level)
c  use quadratic formula to find zero pbe
c           tave=0.5*(t(n)+t(n-1)+546.3)
c           x=temdif(n-1)
c           y=(temdif(n)-temdif(n-1))/deltah(n)
c           twoa=g*y/tave
c           b=g*x/tave
c           c=pbe(n-1)
c           arg=-b/twoa
c           argdsc=b*b-2.*twoa*c
c
c           if(argdsc.ge.0)then
c               disc=sqrt(argdsc)/twoa
c               write(6,352)twoa,b,c,x,y,arg,disc,argdsc
c 352           format(8f15.5)

c               do isign=-1,1,2
c                   hh=arg+abs(disc)*isign
c                   frac=hh/deltah(n)

c                   if(abs(frac-0.5).le.0.5)then ! zero crossing detected
c                       icp=icp+1
c                       pbecr(icp,1)=
c       1                   exp((alog(p(n))-alog(p(n-1)))*frac+alog(p(n-1)))
c                       pbecr(icp,2)=0.
c                       pbecr(icp,3)=0.
c                       pbecr(icp,4)=ht(n-1)+hh

c                       write(6,1089)n,pbe(n),pbe(n-1),hh,pbecr(icp,1)
c     &                          ,pbecr(icp,4),ht(n),hh,deltah(n),frac
c 1089                  format(i3,9f11.3)
c                   endif

c               enddo

c           endif

1100    continue

c       write(6,464)icp,ict,n1,nlevel
 464    format(' icp,ict,n1,nlevel',4i5)
c

!       case when equlibrium level is above top of domain
        pos_area_max = max(pos_area_max,pbe(nlevel) - nbe_min)


!       at least one region of positive area in sounding
        if(pos_area_max .gt. 0.0)then
            pbeneg=nbe_min
            velneg=sqrt(2.0*abs(pbeneg))

        else ! case when no positive area exists anywhere in sounding
            pos_area_max=-1.0
            pbeneg = r_missing_data ! -1e6
            velneg = 0.

        endif

        if(idebug .ge. 2)then
          write(6,485)pos_area_max,pbeneg,velneg
 485      format(' pos_area_max',f10.1,' pbeneg',f10.1,' velneg',f10.1)
        endif

        return
        end
c
!
c       1991    steve albers
c
        function wtblb(p,tc,w,iout,istatus)

cdoc    calculate wet bulb, given p,t,w

        thetae=thaek(p,tc,w)
        call msad5(wtblb_arg,p,thetae,tc,20.,slope,i1,i2,ia,iout
     1            ,istatus)      

        if(istatus .ne. 1)then
            write(6,*)' error in wtblb',p,tc,w,iout
!           stop
        endif

        wtblb = wtblb_arg
        return
        end
c
c
c
        subroutine blayr(p,t,q,pmean,tmean,wmean,thknes,nlevel,hh,
     +                    lowest,io)

cdoc    calculate boundary layer mean values from an input sounding

!       steve albers 1991

        include 'lapsparms.for'
        integer mxl
        parameter (mxl=max_lvls+1) ! number of 3d levels plus the sfc level

        real intlog
        dimension p(mxl),t(mxl),q(mxl)
        tvirt(tt,qq)=tt/(1.-qq*.37803)
        thick(p1,p2,tc,qq)=alog(p1/p2)*tvirt((tc+273.15),qq)*.09604

        if(nlevel.lt.2)goto9000
        zero=0.
        hh=0.
        sumwt=0.
        sumt=0.
        sumq=0.
        if(io.ge.2)write(6,1)lowest,p(lowest),t(lowest),q(lowest),sumwt
     +                        ,sumwt,sumt,sumq
c
        nlow=lowest+1

        do 100 i=nlow,nlevel
        if((p(lowest)-p(i))-thknes)50,150,200
 50     wt=p(i)-p(i-1)
        sumwt=sumwt+wt
        at=t(i-1)
        aq=q(i-1)
        argi=1./alog(p(i)/p(i-1))
        bt=(t(i)-t(i-1))*argi
        bq=(q(i)-q(i-1))*argi
        intlog=p(i)*(alog(p(i))-1.)-p(i-1)*(alog(p(i-1))-1.)
        pdelt=p(i)-p(i-1)
        argt=(at-(bt*alog(p(i-1))))*pdelt+bt*intlog
        argq=(aq-(bq*alog(p(i-1))))*pdelt+bq*intlog
        sumt=sumt+argt
        sumq=sumq+argq
        argtp=argt/pdelt
        argqp=argq/pdelt
        hh=hh+thick(p(i-1),p(i),argtp,argqp)
        if(io.ge.2)write(6,1)i,p(i),t(i),q(i),wt,sumwt,sumt,sumq,bt,bq
     &  ,zero,hh
 100    continue
        goto9000
c
 150    wt=p(i)-p(i-1)
        goto250
c
 200    continue
        wt=-(thknes-p(lowest)+p(i-1))
 250    continue
        pbound=p(lowest)-thknes
        sumwt=sumwt+wt
        at=t(i-1)
        aq=q(i-1)
        argi=1./alog(p(i)/p(i-1))
        bt=(t(i)-t(i-1))*argi
        bq=(q(i)-q(i-1))*argi
        intlog=pbound*(alog(pbound)-1.)-p(i-1)*(alog(p(i-1))-1.)
        pdelt=pbound-p(i-1)
        argt=(at-(bt*alog(p(i-1))))*pdelt+bt*intlog
        argq=(aq-(bq*alog(p(i-1))))*pdelt+bq*intlog
        sumt=sumt+argt
        sumq=sumq+argq

        if(pdelt .ne. 0.)then
            hh=hh+thick(p(i-1),pbound,argt/pdelt,argq/pdelt)
        endif

        if(io.ge.2)write(6,1)i,p(i),t(i),q(i),wt,sumwt,sumt,sumq,bt,bq
     +              ,intlog,hh
! hongli jiang: w>=d+3 from f6.4 to f7.4 10/14/2013
 1      format(1x,i2,f6.0,f6.1,f7.4,2f7.1,2f11.4,f11.6,f9.6,f9.3,f6.2)
c
        pmean=p(lowest)-.5*thknes
        sumwti=1./sumwt
        tmean=sumt*sumwti
        qmean=sumq*sumwti
        wmean=qmean/(1.-qmean)
        goto9999
c
 9000   if(io.ge.2)write(6,9001)
 9001   format(' not enough levels too obtain a boundary layer')
 9999   return
        end
c
c
c
c
c
        subroutine llcl(lcl,tlcl,plcl,p,tc,w,io) ! in meters

cdoc    calculate lcl properties from an input parcel

!       steve albers 1991

        real lcl,kappa
        esl(x)=6.1121*exp(17.67*x/(x+243.5)) 
        data epsiln/.62197/,gammai/102.4596/ ! lapse rate m/deg
        tk=273.15+tc
        tlcl=tk
        kappa=.28613105*(1.-.23*w)
        slope=5.
        e=(w*p)/(epsiln+w)
        iter=1
        goto140
 100    continue
        iter=iter+1
        if(io.ge.1)write(6,25)tlcl,res,slope,delta
 25     format(' tlcl=',5f12.5)
        if(abs(res)-.05)150,150,130
 130    tlcl=tlcl+delta
 140    res=(esl(tlcl-273.15)/e)**kappa*tk-tlcl
        if(iter.ne.1)slope=(res-resold)/delta
        delta=-res/slope
        resold=res
        goto100
 150    tlcl=tlcl-273.15
        plcl=p*esl(tlcl)/e
        lcl=(tc-tlcl)*gammai
        return
        end
c
c
        subroutine lcl_fast(p,tc,td,hlcl,tlcl,plcl)

cdoc    calculate lcl properties from an input sounding (efficiently)

!       steve albers 1991

        esl(x)=6.1121*exp(17.67*x/(x+243.5)) 
        data epsiln/.62197/,gammai/102.4596/        ! lapse rate m/deg

        tlcl = td - (0.212 + 0.001571 * td - 0.000436 * tc) * (tc - td)
        plcl=p*esl(tlcl)/esl(td)
        hlcl=(tc-tlcl)*gammai

        return
        end
c
c
c
c
c
        subroutine itplv(p,param,nlevel,pint,parman,io,istatus)

cdoc    interpolate any parameter from a pressure sounding to a specific pres

!       steve albers 1991

        include 'lapsparms.for'
        integer mxl
        parameter (mxl=max_lvls+1) ! number of 3d levels plus the sfc level

        dimension p(mxl),param(mxl)

        if(p(1) .lt. pint)then
            write(6,*)' error in itplv: p(1) < pint',p(1),pint
            istatus = 0
            print *, p(:)  ! wnidb
            return
        endif

        do 100 n=1,nlevel
        if(p(n)-pint)300,200,100
 100    continue
c
 200    parman=param(n)
        goto999
c
 300    frac=alog(pint/p(n-1))/alog(p(n)/p(n-1))
        parman=param(n-1)+frac*(param(n)-param(n-1))
        if(io.ge.2)write(6,666)n,pint,frac,parman,p(n),p(n-1),param(n)
     #  ,param(n-1)
 666    format(' interpolating to level',i3,f10.3,f8.5,f8.3,2f7.1,2f6.1)
c
 999    istatus = 1
        return
        end
c
        subroutine newtn(x,xold,y,yold,slope,iter,io,fence,istatus)

cdoc    newton iteration

!       steve albers 1991

        if(iter.gt.0)slope=(y-yold)/(x-xold)
        delta=-y/slope
        delta=amax1(-fence,amin1(delta,fence))
        xold=x
        x=x+delta
        iter=iter+1
        if(iter .gt. 100)then
            write(6,*)' error in newtn: too many iterations',iter
            istatus = 0
            return
        endif

        if(io.ge.2)write(6,1)iter,x,xold,y,yold,slope,delta
 1      format(' newton iter',i4,2f11.4,2f11.7,2e11.3)
        yold=y

        istatus = 1
        return
        end
c
c
        subroutine msad5
     ^  (temnew,presnw,thetae,tguess,slopeg,slope,i1,i2,ia,io,istatus)       

cdoc    calculate along a moist adiabat. solve for t, given thetae and p

!       steve albers 1991

        dimension tempnw(4)

        if(thetae .gt. 700.)then
            write(6,*)' error: passed in thetae too large in msad5'
     1                ,thetae,presnw
            istatus = 0
            return
        endif

        i1=0
        i2=0
        ia=0
c
c  estimate parcel temp at 500mb
        if(presnw.ne.500.)goto725
        tempnw(1)=thetae-(307.260+72.122*exp((thetae-382.635)*.0141672))
        tempnw(1)=tempnw(1)+0.65*exp(-(.077*(tempnw(1)+27.))**2)
        tempnw(1)=tempnw(1)-0.50*exp(-(.200*(tempnw(1)+4.0))**2)
        if(thetae.lt.254.90.or.thetae.gt.361.53)goto725
        temnew=tempnw(1)
        goto900
c
 725    tempnw(1)=tguess
        slope=slopeg
c  enter iterative loop if thetae is out of range of approximate formula
c  determine latent heat released above 500mb
 750    iter=1
 760    epsiln=thae(tempnw(iter),tempnw(iter),presnw)-thetae
        i1=i1+1

        if(i1 .gt. 200000)then
            write(6,*)' error, too many i1 loops in msad5'
            istatus = 0
            return
        elseif(i1 .gt. 199980)then
            write(6,*)' warning, too many i1 loops in msad5'
            write(6,*)thetae,tempnw(iter),presnw,epsiln  
        endif

c  test for convergence
        if(abs(epsiln).lt..01)goto890
        if(io.ge.3.and.iter.eq.1)
     ^  write(6,666)iter,tempnw(iter),slope,epsiln,delta,epsold,delold
     ^  ,slopeg
        if(iter.gt.1)slope=(epsiln-epsold)/delold
        delta=-epsiln/slope
        iter=iter+1
        if(io.ge.3)
     ^  write(6,666)iter,tempnw(iter),slope,epsiln,delta,epsold,delold
 666    format(' msad5',i2,2x,7f11.5)
        tempnw(iter)=tempnw(iter-1)+delta
        i2=i2+1
        if(iter.eq.4)goto850
        epsold=epsiln
        delold=delta
        goto760
c
c  use aitken's formula to accelerate convergence
 850    ratio=delta/delold
        ratabs=abs(ratio)
        ratio=amin1(ratabs,.5)*ratio/ratabs
        if(io.ge.3)write(6,666)
     ^  iter,tempnw(iter),delta,delold,ratio
        tempnw(1)=tempnw(3)+delta/(1.-ratio)
        ia=ia+1
        goto750
c
 890    temnew=tempnw(iter)
        iout=0
        if(iout.eq.0.or.io.lt.3)goto900
        write(6,666)iter,temnew,slope,epsiln
        write(6,203)
 203    format('  ')
 900    continue
c       write(6,901)thetae,presnw,temnew,i1,i2,ia
c901    format(' travelled down moist adiabat, thetae=',f7.2
c     #  ,' new pressure=',f7.2,' new temp=',f7.2,3i3)
        istatus = 1
        return
        end
c
!
c      1991     steve albers
c
       function thae(tc,td,p)

cdoc   calculate theta(e), given t, td, p
c
c   computes the equivalent potential tempurature (k).
c    (using the rossby definition)
        esl(x)=6.1121*exp(17.67*x/(x+243.5)) 
        t=tc+273.15
        e=esl(td)
        pmei=1./(p-e)
        w=(e*.62197)*pmei
        t=t*(1.+w)/(1.+w*.62197)
       cp=0.2396
       thae=thd(p,t,w,pmei) * exp(rl(tc)*w/(cp*t))
       return
       end
c
       function thd(p,t,w,pmei)

cdoc   computes the dry air potential temperature

       real ak
       parameter (ak=.28613105)

       aks=ak * (1.0+1.608*w)/(1.0 + 1.941569*w)
       thd=t * ((1000.0*pmei)**aks)
       return
       end

       function rl(tm2)
cdoc   latent heat of evaporation
c      tm2=t-273.15
c      rl=597.31-0.589533*tm2+0.001005333*(tm2*tm2)
       rl=597.31-((0.589533+0.001005333*tm2)*tm2)
       return
       end
c
c
       subroutine dryad(p1,p2,t1,t2,w,io)

!       steve albers 1991

cdoc   computes the dry air potential temperature (theta), given p and t
       ak=.28613105
c       aks=ak * (1.0+1.608*w)/(1.0 + 1.941569*w)
c       e=w*p1/(0.62197 +w)
c       pmei=1./(p1-e)
c       t2=t1 * ((p2*pmei)**aks)
        t2=t1*(p2/p1)**ak
        if(io.ge.2)write(6,1)p1,p2,t1,t2,w
 1      format(' p1,p2,t1,t2,w',4f10.3,f10.5)
       return
       end
c
       function thete(t,td,alt,hgt)

cdoc   compute theta(e), given t, td, altimeter setting, and elevation (hgt)

c  t temp (c), td dew pt (c), alt (altimeter setting in.)
c  hgt height asl (m).
c  convert pa from inches to mb
       pa=alt*33.8639
c  sfc pressure
       pa=pa*(1.0/(1.0+hgt*2.2222e-05))**5.25530
        thete=thae(t,td,pa)
        return
        end
c
c
c
       function xmxrat(pres,dewp)
cdoc   compute mixing ratio (gm/gm) given dew point temp and the pressure (mb)
       ratmix=exp(21.16-5415.0/dewp)
       ratmix=ratmix/pres
       if(ratmix.lt.(5.0e-05)) ratmix=5.0e-05
       xmxrat=ratmix
       return
       end
c
c
c
       function thaek(p,tc,w)
cdoc   computes the equivalent potential tempurature (k).
cdoc   (using the rossby definition)
        q=w/(1.0+w)
        e=(p*q)/.62197
        pmei=1./(p-e)
        t=tc+273.15
        t=t*(1.+w)/(1.+w*.62197)
       cp=0.2396
       thaek=thd(p,t,w,pmei) * exp(rl(tc)*w/(cp*t))
       return
       end
c

        function oe_fast(t_in,td_in,pres)

!       steve albers 1991

!       t  in deg c  (input) valid input temp range is -60c to +60c
!       td in deg c  (input) max dewpoint depression allowed is 30
!       pres in mb   (input) min allowed is 500mb
!       oe_fast in c (output)

cdoc    quick way to get theta(e) using lookup table

!       1991    steve albers

        character*31 ext
        character*150 directory
        integer len_dir

        integer pres2_low,pres2_high,pres2_interval,n_pres2
        parameter (pres2_low = 250)
        parameter (pres2_high = 1100)
        parameter (pres2_interval = 10)

        parameter (n_pres2 = (pres2_high - pres2_low) / pres2_interval)

        integer t_low,t_high,t_interval,n_t
        parameter (t_low  = -60)
        parameter (t_high = +60)
        parameter (t_interval = 2)

        parameter (n_t = (t_high - t_low) / t_interval)

        integer tdprs_low,tdprs_high,tdprs_interval,n_tdprs
        parameter (tdprs_low  = 0)
        parameter (tdprs_high = +30)
        parameter (tdprs_interval = 1)

        parameter (n_tdprs = (tdprs_high - tdprs_low) / tdprs_interval)

        real thetae_lut(0:n_t,0:n_tdprs,0:n_pres2)

        logical l_write_lut

        save init,thetae_lut
        data init/0/

        l_write_lut = .false.

        if(init .eq. 0)then
            init = 1

            if(l_write_lut)then
                ext = 'dat'
                call get_directory(ext,directory,len_dir)

                write(6,*)' reading thetae_lut.dat from file'
                open(11,file=directory(1:len_dir)//'thetae_lut.dat'
     1                    ,form='unformatted',status='old',err=10)
                read(11,err=10)thetae_lut
                close(11)
                goto20
10              continue
                close(11)

                write(6,*)
     1            ' generating thetae_lut.dat - no valid file exists'
                open(12,file=directory(1:len_dir)//'thetae_lut.dat'
     1                         ,form='unformatted',status='unknown')
            else
                write(6,*)' generating thetae lut in memory...'
            endif

            i = 0
            do t = t_low,t_high,t_interval

                j = 0
                do tdprs = tdprs_low,tdprs_high,tdprs_interval

                    k = 0
                    do pres2 = pres2_low,pres2_high,pres2_interval

                        td = t-tdprs
                        thetae_lut(i,j,k) = oe(t,td,pres2)

                        if(i/10*10 .eq. i .and. j/10*10 .eq. j .and.
     1                                        k/10*10 .eq. k)then
                            write(6,201)t,td,pres2,thetae_lut(i,j,k)
201                         format(1x,2f10.0,f10.0,f10.2)
                        endif

                        k = k + 1

                        enddo
                    j = j + 1

                enddo
                i = i + 1

            enddo

            if(l_write_lut)then
                write(12)thetae_lut
                close(12)
            endif

20      endif

        tdprs = t_in - td_in

        if(t_in  .ge. float(t_high))    t_in  = float(t_high)     - 1e-5
        if(tdprs .ge. float(tdprs_high))tdprs = float(tdprs_high) - 1e-5
        if(pres  .ge. float(pres2_high))pres  = float(pres2_high) - 1e-2

        ri = (t_in  - float(t_low)    ) / float(t_interval)
        rj = (tdprs - float(tdprs_low)) / float(tdprs_interval)
        rk = (pres  - float(pres2_low)) / float(pres2_interval)

        i = int(ri)
        j = int(rj)
        k = int(rk)

        fraci = ri - i
        fracj = rj - j
        frack = rk - k

        z1=thetae_lut(i  , j  ,k)
        z2=thetae_lut(i+1, j  ,k)
        z3=thetae_lut(i+1, j+1,k)
        z4=thetae_lut(i  , j+1,k)

        thetae_low =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                    - (z2+z4-z3-z1)*fraci*fracj

        z1=thetae_lut(i  , j  ,k+1)
        z2=thetae_lut(i+1, j  ,k+1)
        z3=thetae_lut(i+1, j+1,k+1)
        z4=thetae_lut(i  , j+1,k+1)

        thetae_high =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                    - (z2+z4-z3-z1)*fraci*fracj

        oe_fast = thetae_low * (1. - frack) + thetae_high * frack

c       residual = abs(oe_fast - oe(t_in,t_in-tdprs,pres))
c       write(6,101)t_in,t_in-tdprs,pres,oe(t_in,t-tdprs,pres),residual
c101    format(1x,3f10.2,f10.4)

        return
        end

        function tmlaps_fast(thetae_in,pres_in)

!       steve albers 1991

cdoc    quick way to get temperature along a moist adiabat given theta(e)
cdoc    this uses a lookup table

!       thetae in k  (input)
!       pres in mb   (input)
!       t moist in k (output)

!       1991    steve albers

        character*31 ext
        character*150 directory
        integer len_dir

        integer thetae_low,thetae_high,thetae_interval,n_thetae
        parameter (thetae_low  = 210)
        parameter (thetae_high = 410)
        parameter (thetae_interval = 1)

        parameter (n_thetae = (thetae_high - thetae_low) / thetae_interv
     1al)

        integer pres_low,pres_high,pres_interval,n_pres
        parameter (pres_low = 50)
        parameter (pres_high = 1100)
        parameter (pres_interval = 5)

        parameter (n_pres = (pres_high - pres_low) / pres_interval)

        real t_moist(0:n_thetae,0:n_pres)

        logical l_write_lut

        save init,t_moist
        data init/0/

        l_write_lut = .false.

        if(init .eq. 0)then
            init = 1

            if(l_write_lut)then
                ext = 'dat'
                call get_directory(ext,directory,len_dir)

                write(6,*)' reading tmlaps_lut.dat from file'
                open(11,file=directory(1:len_dir)//'tmlaps_lut.dat'
     1                    ,form='unformatted',status='old',err=10)
                read(11,err=10)t_moist
                close(11)
                goto20
10              continue
                close(11)

                write(6,*)
     1            ' generating tmlaps_lut.dat - no valid file exists'
                open(11,file=directory(1:len_dir)//'tmlaps_lut.dat'
     1                             ,form='unformatted',status='new')
            else
                write(6,*)' generating tmlaps lut in memory...'
            endif

            i = 0
            do thetae = thetae_low,thetae_high,thetae_interval
                j = 0

                do pres = pres_low,pres_high,pres_interval

                    t_moist(i,j) = tmlaps(thetae-273.15,pres)

                    if(i/5*5 .eq. i .and. j/5*5 .eq. j)then
!                       write(6,201)pres,thetae,t_moist(i,j)
201                     format(1x,f10.0,f10.0,f10.2)
                    endif
                    j = j + 1

                enddo
                i = i + 1

            enddo

            if(l_write_lut)then
                write(11)t_moist
                close(11)
            endif

20      endif

        ri = (thetae_in - float(thetae_low)) / float(thetae_interval)
        rj = (pres_in   - float(pres_low)  ) / float(pres_interval)

        if(ri .ge. float(n_thetae))ri = float(n_thetae) - 1e-2
        if(rj .ge. float(n_pres))  rj = float(n_pres)   - 1e-2

        i = int(ri)
        j = int(rj)

        fraci = ri - i
        fracj = rj - j

        z1=t_moist(i  , j  )
        z2=t_moist(i+1, j  )
        z3=t_moist(i+1, j+1)
        z4=t_moist(i  , j+1)

        tmlaps_fast =  z1+(z2-z1)*fraci+(z4-z1)*fracj
     1                    - (z2+z4-z3-z1)*fraci*fracj

c       residual = abs(tmlaps_fast - tmlaps(thetae-273.15,pres))
c       write(6,101)tmlaps_fast,tmlaps(thetae-273.15,pres),residual
c101    format(1x,2f10.2,f10.4)

        return
        end

        function dwpt_laps(t,rh)
c
cdoc this function returns the dew point (celsius) given the temperature
cdoc (celsius) and relative humidity (%).
c
c       baker,schlatter 17-may-1982     original version
c
c       the formula is used in the
c   processing of u.s. rawinsonde data and is referenced in parry, h.
c   dean, 1969: "the semiautomatic computation of rawinsondes,"
c   technical memorandum wbtm edl 10, u.s. department of commerce,
c   environmental science services administration, weather bureau,
c   office of systems development, equipment development laboratory,
c   silver spring, md (october), page 9 and page ii-4, line 460.

c   this has been updated to test for out of bounds values (2004 steve albers)

        if(rh .lt. 0. .or. rh .gt. 100.)then
            call get_r_missing_data(r_missing_data,istatus)
            dwpt_laps = r_missing_data
            return
        endif

        if(t .lt. -100. .or. t .gt. +100.)then
            call get_r_missing_data(r_missing_data,istatus)
            dwpt_laps = r_missing_data
            return
        endif

        x = 1.-0.01*rh
c   compute dew point depression.
        dpd =(14.55+0.114*t)*x+((2.5+0.007*t)*x)**3+(15.9+0.117*t)*x**14
        dwpt_laps = t-dpd
        return
        end


        function twet_fast(t_c,td_c,pres_mb)

!       steve albers 1991
cdoc    this is a fast approximate wet bulb routine using lookup tables
!       warning: this routine is only active below 500mb due to size restriction
!       of the lookup table. max dewpoint depression allowed is 30
!       further approximation used (twet = t_c) when t_c is outside valid range

!       twet_fast in c (output)

        if(pres_mb .ge. 500. .and. t_c .ge. -60. 
     1                       .and. t_c .le. +60.)then ! valid ranges to input
            thetae_k = oe_fast(t_c,td_c,pres_mb) + 273.15
            twet_fast = tmlaps_fast(thetae_k,pres_mb)

        else
            twet_fast = t_c

        endif

        return
        end

       subroutine li_laps(tc,td,pr,i4time,imax,jmax,li,flag
     !                   ,istatus)

cdoc   compute 2-d grid of li, given a grid of parcels to launch

!      1991     steve albers

       real tc(imax,jmax) ! input t  in deg f
       real td(imax,jmax) ! input td in deg f
       real pr(imax,jmax) ! input pr in mb
       real li(imax,jmax) ! output li in deg k/c

       real t500laps(imax,jmax) ! used locally only

       character*13 filename13

       call get_r_missing_data(r_missing_data,istatus)
       if(istatus .ne. 1)then
           write(6,*)' error reading r_missing_data'
           return
       endif

       call constant(li,r_missing_data,imax,jmax)
       call get_laps_cycle_time(laps_cycle_time,istatus)
!      get 500 mb temp field
       if(flag .eq. 0.0)then ! use lt1 (or equivalent) file
           write(6,*)' reading 500 temp from lt1 (or equivalent) file'
           k_level = 500
           call get_temp_2d(i4time,laps_cycle_time,i4time_nearest
     1                          ,k_level,imax,jmax,t500laps,istatus)

           if(istatus .ne. 1)then
               write(6,*)' aborting li calculation, no 500 temps'
               goto999
           endif

       elseif(abs(flag) .ge. 1e6)then ! use sounding file

           i4time_round = ((i4time+21600)/43200) * 43200
           open(52,file='disk$laps_s90x:'
     1          //filename13(i4time_round,'dmn'),status='old')
           read(52,*)
           read(52,*)
           read(52,*)
           read(52,*)
           read(52,*)

10         read(52,*,end=20)pres,temp

           if(pres .eq. 500.)then
               write(6,*)' using sounding 500mb temp',temp

               do i=1,imax
               do j=1,jmax
                   t500laps(i,j) = temp + 273.15
               enddo
               enddo

               goto50

           else
               goto10

           endif

20         write(6,*)' could not find 500mb temp'
           istatus = 0
           close(52)
           return

50         close(52)

       else ! use input temp
           write(6,*)' using input 500mb temp'
           do i=1,imax
           do j=1,jmax
                t500laps(i,j) = flag
           enddo
           enddo

       endif

!      calculate fields
       do j=1,jmax
          do i=1,imax

              if(abs(tc(i,j)) .le. 1e5
     1   .and. abs(td(i,j)) .le. 1e5
     1   .and. abs(pr(i,j)) .le. 1e5
     1   .and. abs(t500laps(i,j)) .le. 1e5)then
                  tcij_c = (tc(i,j)-32.)/1.8
                  tdij_c = (td(i,j)-32.)/1.8
                  t500laps_c = t500laps(i,j) - 273.15
                  li(i,j) = func_li(tcij_c,tdij_c,pr(i,j),t500laps_c
     1                             ,r_missing_data)
              else
                  li(i,j) = r_missing_data
              endif

          end do ! j
      end do ! i

999   return
      end

        function func_li(t_c,td_c,psta_mb,t500_c,r_missing_data)

cdoc    calculate li given an input parcel

        td_in = min(td_c,t_c)

        thetae = thae(t_c,td_in,psta_mb)

!       t500_parcel = t500(thetae)
        call thetae_to_t500(thetae,t500_parcel,istatus)

        if(istatus .eq. 1)then
            func_li = t500_c - t500_parcel

        else
            write(6,1,err=2)psta_mb,t_c,td_c
1           format(' warning: setting li to missing - p,t,td= ',3f10.0)       
2           func_li = r_missing_data

        endif

        return
        end


        subroutine thetae_to_t500(thetae,t500,istatus)

cdoc    given theta(e), what is t-500mb?

        real diff(10)
c
c calculate li
        iter=1
        iterc1=1
        iterc2=1

c estimate parcel temp at 500mb
        t500nw=thetae-(307.260+72.122*exp((thetae-382.635)*.0141672))
        t500nw=t500nw+0.65*expm(-(.077*(t500nw+27.))**2)
        t500nw=t500nw-0.50*expm(-(.200*(t500nw+4.0))**2)
!       if(output.eq.1)write(31,*)thetae,' t500 ',t500nw
        if(thetae.lt.254.90.or.thetae.gt.361.53)then
            t500ol=t500nw
c           write(6,*)' using iterative method to get 500mb temp'
c           write(31,*)' using iterative method to get 500mb temp'
c
c enter iterative loop if thetae is out of range of approximate formula
c determine latent heat released above 500mb
750         fudge=thae(t500ol,t500ol,500.)-1.21937*(t500ol+273.15)
            thp5=thetae-fudge
            t500nw=(thp5/1.21937)-273.15
!           write(31,666)
!       1   iter,iok,jok,xok,yok,ius,jus,xus,yus,t500nw,t500ol,fudge,ratio
!666        format(i1,2i3,f6.2,f4.1,2i3,2f7.3,2x,4f10.5)
            iterc1=iterc1+1

            if(iterc1 .gt. 1000)then
                write(6,*)
     1               ' warning in thetae_to_t500, too many iterations, '       
     1              ,' thetae = ',thetae
                istatus = 0
                return
            endif

c test for convergence
            diff(iter)=t500nw-t500ol
            if(abs(diff(iter)).lt..10)goto900
            iter=iter+1
            iterc2=iterc2+1
c
c use aitken's formula to accelerate convergence
            if(iter.eq.3)then
                ratio=diff(2)/diff(1)
!               write(31,666)
!       1       iter,iok,jok,xok,yok,ius,jus,xus,yus,t500nw,diff(2)
!       1           ,diff(1),ratio
                t500nw=t500ol+diff(2)/(1.-ratio)
                t500ol=t500nw
                iter=1
                iterai=iterai+1
                goto750
            endif

800         t500ol=t500nw
            goto750

        endif

900     t500 = t500nw

        istatus = 1
        return
        end


        function expm(x)

cdoc    calculate exp function with check to avoid underflow with large inputs.

        if(x .ge. -70.)then
            expm = exp(x)
        else
            expm = 0.
        endif

        end


        function devirt_td(t_k,td_k,p_pa)

cdoc    this function yields the approximate temperature given the virtual temp
cdoc    tv from the mthermo library is called. 

!       please suggest improvements to this routine to steve albers at fsl.

!       t_k       input temp in k
!       td_k      input dew point temp in k
!       p_pa      input pressure as pascals
!       devirt_td output devirtualized temp as k

        t_c  = t_k  - 273.15
        td_c = td_k - 273.15
        p_mb = p_pa / 100.

        tv_c = tv(t_c,td_c,p_mb)

        tv_k = tv_c + 273.15

        devirt_td = t_k - (tv_k-t_k)

        return
        end

        function devirt_sh(t_k,sh,p_pa)

cdoc    this function yields the approximate temperature given the virtual temp
cdoc    tv from the mthermo library is called. 

!       please suggest improvements to this routine to steve albers at fsl.

!       t_k       input temp in k
!       sh        input specific humidity (dimensionless)
!       p_pa      input pressure as pascals
!       devirt_td output devirtualized temp as k

        t_c  = t_k  - 273.15
        p_mb = p_pa / 100.

        tv_c = tv_sh(t_c,sh,p_mb)

        tv_k = tv_c + 273.15

        devirt_sh = t_k - (tv_k-t_k)

        return
        end

        function devirt_rh(t_k,rh,p_pa)

cdoc    this function yields the approximate temperature given the virtual temp
cdoc    devirt_td from the laps library is called. 

!       please suggest improvements to this routine to steve albers at fsl.

!       t_k       input temp in k
!       rh        input rh as fraction
!       p_pa      input pressure as pascals
!       devirt_rh output devirtualized temp as k

        t_c  = t_k  - 273.15

        rh_pct = rh * 100.

        td_c = dwpt(t_c,rh_pct)
        td_k = td_c + 273.15

        devirt_k = devirt_td(t_k,td_k,p_pa)

        devirt_rh = devirt_k

        return
        end

        function tv_sh(t,sh,p)
c
cdoc    this function returns the virtual temperature tv (celsius) of
cdoc    a parcel of air at temperature t (celsius), dew point td
cdoc    (celsius), and pressure p (millibars). the equation appears
cdoc    in most standard meteorological texts.
c
c       baker,schlatter 17-may-1982     original version
c       albers                 1994     modified for sh input
c
        data cta,eps/273.16,0.62197/
c   cta = difference between kelvin and celsius temperatures.
c   eps = ratio of the mean molecular weight of water (18.016 g/mole)
c         to that of dry air (28.966 g/mole)
        tk = t+cta
c   calculate the dimensionless mixing ratio.
        w = sh / (1. - sh)
        tv_sh = tk*(1.+w/eps)/(1.+w)-cta
        return
        end

      function tsa_fast(os,p)
c
cdoc  this function returns the temperature tsa (celsius) on a saturation
cdoc  adiabat at pressure p (millibars). os is the equivalent potential
cdoc  temperature of the parcel (celsius). sign(a,b) replaces the
cdoc  algebraic sign of a with that of b.
c
c       baker,schlatter 17-may-1982     original version
c       modification for better convergence, keith brewster, feb 1994.
c
c   b is an empirical constant approximately equal to the latent heat
c   of vaporization for water divided by the specific heat at constant
c   pressure for dry air.
c
      real b
      parameter (b=2.6518986)
      a= os+273.15
c
c   above 200 mb figure all the moisture is wrung-out, so
c   the temperature is that which has potential temp of theta-e.
c   otherwise iterate to find combo of moisture and temp corresponding
c   to thetae.
c
      if( p.lt.200.) then
        tq=a*((p/1000.)**.286)
      else
c   d is an initial value used in the iteration below.
        d= 120.
c   tq is the first guess for tsa.
        tq= 253.15
        x = 0.
c
c   iterate to obtain sufficient accuracy....see table 1, p.8
c   of stipanuk (1973) for equation used in iteration.
        do 1 i= 1,25
           tqc= tq-273.15
           d= 0.5*d

           x_last = x

           x= a*exp(-b*w_fast(tqc,p)/tq)-tq*((1000./p)**.286)

c          write(6,*)' tsa_fast: t,err= ',i,tqc,x
           if (abs(x).lt.1e-3) go to 2
c
           if (x_last * x .lt. 0.) then
               slope = (x-x_last) / (tq - tq_last)
               delta = - x / slope
               ad = amin1(abs(delta),d)
               tq_last = tq
               tq = tq + sign(ad,delta)
           else
               tq_last = tq
               tq= tq+sign(d,x)
           end if

 1      continue
        end if
 2      tsa_fast = tq-273.15
c       write(6,*)' tsa_fast: i,t,err= ',i,tsa_fast,x
        return
        end
c

        subroutine get_tw_approx_2d(t_k,td_k,p_pa,ni,nj,tw_k)

cdoc    calculate wet bulb, using a fast approximate method

!       steve albers 1991

!       this routine is fast but only accurate near 0 degrees c (273k)

        real t_k(ni,nj)     ! input
        real td_k(ni,nj)    ! input
        real p_pa(ni,nj)    ! input
        real tw_k(ni,nj)    ! output

        const = alog(0.5)

        if(t_k(1,1) .eq. 0. .or. td_k(1,1) .eq. 0.)then
            write(6,*)' bad input data to get_tw_approx_2d'
            return
        endif

        do j = 1,nj
        do i = 1,ni
            t_f =  t_k(i,j)  * 1.8 - 459.67
            td_f = td_k(i,j) * 1.8 - 459.67
            tmid_f = 0.5 * (t_f + td_f)
            depress50 = 13.4 + tmid_f * 0.1
            rh = 0.5 ** ((t_f - td_f)/depress50)
            start = 240. - (t_k(i,j) - 240.) / 6.
            ratio = (t_k(i,j) - start)/100. * 
     1                                   (1.0 + 0.9 * (1.0 - sqrt(rh)))        
            dtw_dtd = ratio  ! (t = td)
            drh_dtd = const / depress50 ! (t = td)
            dtw_drh = dtw_dtd / drh_dtd
            tw_f = t_f - (rh - 1.0) * dtw_drh
            tw_k(i,j) = (tw_f + 459.67) / 1.8
        enddo ! i
        enddo ! j

        return
        end


        subroutine get_tw_2d(t_k,td_k,p_pa,ni,nj,tw_k)

cdoc    calculate 2-d grid of wet bulb, using 'tw_fast'

!       steve albers 1991
!       warning: this routine may not work because it calls tw_fast

        real t_k(ni,nj)     ! input
        real td_k(ni,nj)    ! input
        real p_pa(ni,nj)    ! input
        real tw_k(ni,nj)    ! output

        do j = 1,nj
        do i = 1,ni
            tw_k(i,j) = 
     1      tw_fast(t_k(i,j)-273.15,td_k(i,j)-273.15,p_pa(i,j)*.01)
     1                                                     + 273.15
        enddo ! i
        enddo ! j

        return
        end

        subroutine get_tw_2d_orig(t_k,td_k,p_pa,ni,nj,tw_k)

cdoc    calculate 2-d grid of wet bulb, using 'tw'

!       steve albers 1991

        real t_k(ni,nj)     ! input
        real td_k(ni,nj)    ! input
        real p_pa(ni,nj)    ! input
        real tw_k(ni,nj)    ! output

        do j = 1,nj
        do i = 1,ni
            tw_k(i,j) = 
     1      tw(t_k(i,j)-273.15,td_k(i,j)-273.15,p_pa(i,j)*.01)
     1                                                     + 273.15
        enddo ! i
        enddo ! j

        return
        end

        function tw_fast(t,td,p)
c
c
cdoc    this function returns the wet-bulb temperature tw (celsius)
cdoc    given the temperature t (celsius), dew point td (celsius)
cdoc    and pressure p (mb).  see p.13 in stipanuk (1973), referenced
cdoc    above, for a description of the technique.

cdoc    warning: this routine may not work because it calls tsa_fast
c
c       baker,schlatter 17-may-1982     original version
c
c   determine the mixing ratio line thru td and p.
        aw = w_fast(td,p)
c
c   determine the dry adiabat thru t and p.
        ao = o(t,p)
        pi = p
c
c   iterate to locate pressure pi at the intersection of the two
c   curves .  pi has been set to p for the initial guess.
        do 4 i= 1,10
           x= .02*(tmr(aw,pi)-tda(ao,pi))
           if (abs(x).lt.0.01) go to 5
 4         pi= pi*(2.**(x))
c   find the temperature on the dry adiabat ao at pressure pi.
 5      ti= tda(ao,pi)
c
c   the intersection has been located...now, find a saturation
c   adiabat thru this point. function os returns the equivalent
c   potential temperature (k) of a parcel saturated at temperature
c   ti and pressure pi.
        aos= os_fast(ti+273.15,pi)-273.15
c   function tsa returns the wet-bulb temperature (c) of a parcel at
c   pressure p whose equivalent potential temperature is aos.
        tw_fast = tsa_fast(aos,p)
        return
        end

        function w_fast(t,p) ! saturation mixing ratio wrt water
c
cdoc    this function returns the mixing ratio (grams of water vapor per
cdoc    kilogram of dry air) given the temperature t (celsius) and pressure
cdoc    (millibars). the formula is quoted in most meteorological texts.
cdoc    note this is a faster version done by steve albers.
c
c       baker,schlatter 17-may-1982     original version
c       albers                 1992     modified for laps
c
        x= eslo(t)
        w_fast= 622.*x/(p-x)
        return
        end

        function wice_fast(t,p) ! saturation mixing ratio wrt ice
c
cdoc    this function returns the mixing ratio (grams of water vapor per
cdoc    kilogram of dry air) given the temperature t (celsius) and pressure
cdoc    (millibars). the formula is quoted in most meteorological texts.
cdoc    note this is a faster version done by steve albers
c
c       baker,schlatter 17-may-1982     original version
c       albers                 1993     modified for laps
c
        x= esilo(t)
        wice_fast= 622.*x/(p-x)
        return
        end

        function os_fast(tk,p)
c
cdoc    this function returns the equivalent potential temperature os
cdoc    (k) for a parcel of air saturated at temperature t (k)
cdoc    and pressure p (millibars).
cdoc    note this is a faster version done by steve albers

c       baker,schlatter 17-may-1982     original version
c
        data b/2.6518986/
c   b is an empirical constant approximately equal to the latent heat
c   of vaporization for water divided by the specific heat at constant
c   pressure for dry air.
        tc = tk - 273.15

!       from w routine
        x= eslo(tc)
        w= 622.*x/(p-x)

        os_fast= tk*((1000./p)**.286)*(exp(b*w/tk))

        return
        end


        subroutine wtblb_lvl(twet_c,p,t,q,wb,mxl,nlevel,pwb0,hwb0)

        real twet_c         ! input:  wet bulb temperature in degrees c 
        real p(mxl)         ! input:  pressure sounding (mb)
        real t(mxl)         ! input:  temperature sounding (c)
        real q(mxl)         ! input:  specific humidity sounding
        real wb(mxl)        ! input:  wet bulb temperature sounding
        integer mxl         ! input:  size of p,t,q arrays
        integer nlevel      ! input:  number of levels in vertical arrays
        real pwb0           ! output: pressure of the wet bulb level
        real hwb0           ! output: height of the wet bulb level (kft agl)

        io = 0
 	iout=min(io,1)                                                          
 	if(wb(1).ge.twet_c)goto150                                                  
!	if(io.ge.2)write(6,87)                                                  
 87	format(' surface wetbulb temp is below twet_c')                          
 	hwb0=0.                                                                 
 	pwb0=0.                                                                 
 	goto390                                                                 
c                                                                         
!       test for bracketing of the wet bulb zero
 150	do 200 n=2,nlevel                                                    
! 	    if(wb(n)*wb(n-1))250,250,200                                            
            if( (wb(n) - twet_c) * (wb(n-1) - twet_c) )250,250,200 
 200    continue                                                           

        write(6,*)' warning, no wet bulb detected in loop'
        goto390

 250    frac = (twet_c - wb(n-1)) / (wb(n) - wb(n-1))
        pwb0 = p(n-1) * (1. - frac) + p(n) * frac

!       integrate hydrostatically to get height from pressure
 	call blayr(p,t,q,dum1,dum2,dum3,p(1)-pwb0,nlevel,hwb0,1,io) 
!	if(io.ge.1)write(6,351)pwb0,hwb0                            
 351	format(' wetbulb zero is at',f6.1,'mb     or',f6.2,'kft  agl')        

        goto400                      ! normal condition

 390    hwb0 = r_missing_data        ! indeterminate condition
        pwb0 = r_missing_data

 400    continue

        return
        end
