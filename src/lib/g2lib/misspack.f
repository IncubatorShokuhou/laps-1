      subroutine misspack(fld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    misspack
!   prgmmr: gilbert          org: w/np11    date: 2000-06-21
!
! abstract: this subroutine packs up a data field using a complex
!   packing algorithm as defined in the grib2 documention.  it
!   supports grib2 complex packing templates with or without
!   spatial differences (i.e. drts 5.2 and 5.3).
!   it also fills in grib2 data representation template 5.2 or 5.3 
!   with the appropriate values.
!   this version assumes that missing value management is being used and that
!   1 or 2 missing values appear in the data.
!
! program history log:
! 2000-06-21  gilbert
! 2004-12-29  gilbert  -  corrected bug when encoding secondary missing values.
!
! usage:    call misspack(fld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
!   input argument list:
!     fld()    - contains the data values to pack
!     ndpts    - the number of data values in array fld()
!     idrsnum  - data representation template number 5.n
!                must equal 2 or 3.
!     idrstmpl - contains the array of values for data representation
!                template 5.2 or 5.3
!                (1) = reference value - ignored on input
!                (2) = binary scale factor
!                (3) = decimal scale factor
!                    .
!                    .
!                (7) = missing value management
!                (8) = primary missing value
!                (9) = secondary missing value
!                    .
!                    .
!               (17) = order of spatial differencing  ( 1 or 2 )
!                    .
!                    .
!
!   output argument list: 
!     idrstmpl - contains the array of values for data representation
!                template 5.3
!                (1) = reference value - set by misspack routine.
!                (2) = binary scale factor - unchanged from input
!                (3) = decimal scale factor - unchanged from input
!                    .
!                    .
!     cpack    - the packed data field (character*1 array)
!     lcpack   - length of packed field cpack().
!
! remarks: none
!
! attributes:
!   language: xl fortran 90
!   machine:  ibm sp
!
!$$$

      integer,intent(in) :: ndpts,idrsnum
      real,intent(in) :: fld(ndpts)
      character(len=1),intent(out) :: cpack(*)
      integer,intent(inout) :: idrstmpl(*)
      integer,intent(out) :: lcpack

      real(4) :: ref
      integer(4) :: iref
      integer,allocatable :: ifld(:),ifldmiss(:),jfld(:)
      integer,allocatable :: jmin(:),jmax(:),lbit(:)
      integer,parameter :: zero=0
      integer,allocatable :: gref(:),gwidth(:),glen(:)
      integer :: glength,grpwidth
      logical :: simple_alg = .false.
      
      alog2=alog(2.0)
      bscale=2.0**real(-idrstmpl(2))
      dscale=10.0**real(idrstmpl(3))
      missopt=idrstmpl(7)
      if ( missopt.ne.1 .and. missopt.ne.2 ) then
         print *,'misspack: unrecognized option.'
         lcpack=-1
         return
      else     !  get missing values
         call rdieee(idrstmpl(8),rmissp,1)
         if (missopt.eq.2) call rdieee(idrstmpl(9),rmisss,1)
      endif
!
!  find min value of non-missing values in the data,
!  and set up missing value mapping of the field.
!
      allocate(ifldmiss(ndpts))
      rmin=huge(rmin)
      if ( missopt .eq. 1 ) then        ! primary missing value only
         do j=1,ndpts
           if (fld(j).eq.rmissp) then
              ifldmiss(j)=1
           else
              ifldmiss(j)=0
              if (fld(j).lt.rmin) rmin=fld(j)
           endif
         enddo
      endif
      if ( missopt .eq. 2 ) then        ! primary and secondary missing values
         do j=1,ndpts
           if (fld(j).eq.rmissp) then
              ifldmiss(j)=1
           elseif (fld(j).eq.rmisss) then
              ifldmiss(j)=2
           else
              ifldmiss(j)=0
              if (fld(j).lt.rmin) rmin=fld(j)
           endif
         enddo
      endif
!
!  allocate work arrays:
!  note: -ifldmiss(j),j=1,ndpts is a map of original field indicating 
!         which of the original data values
!         are primary missing (1), sencondary missing (2) or non-missing (0).
!        -jfld(j),j=1,nonmiss is a subarray of just the non-missing values from
!         the original field.
!
      !if (rmin.ne.rmax) then
        iofst=0
        allocate(ifld(ndpts))
        allocate(jfld(ndpts))
        allocate(gref(ndpts))
        allocate(gwidth(ndpts))
        allocate(glen(ndpts))
        !
        !  scale original data
        !
        nonmiss=0
        if (idrstmpl(2).eq.0) then        !  no binary scaling
           imin=nint(rmin*dscale)
           !imax=nint(rmax*dscale)
           rmin=real(imin)
           do j=1,ndpts
              if (ifldmiss(j).eq.0) then
                nonmiss=nonmiss+1
                jfld(nonmiss)=nint(fld(j)*dscale)-imin
              endif
           enddo
        else                              !  use binary scaling factor
           rmin=rmin*dscale
           !rmax=rmax*dscale
           do j=1,ndpts
              if (ifldmiss(j).eq.0) then
                nonmiss=nonmiss+1
                jfld(nonmiss)=nint(((fld(j)*dscale)-rmin)*bscale)
              endif
           enddo
        endif
        !
        !  calculate spatial differences, if using drs template 5.3
        !
        if (idrsnum.eq.3) then        ! spatial differences
           if (idrstmpl(17).ne.1.and.idrstmpl(17).ne.2) idrstmpl(17)=2
           if (idrstmpl(17).eq.1) then      ! first order
              ival1=jfld(1)
              do j=nonmiss,2,-1
                 jfld(j)=jfld(j)-jfld(j-1)
              enddo
              jfld(1)=0
           elseif (idrstmpl(17).eq.2) then      ! second order
              ival1=jfld(1)
              ival2=jfld(2)
              do j=nonmiss,3,-1
                 jfld(j)=jfld(j)-(2*jfld(j-1))+jfld(j-2)
              enddo
              jfld(1)=0
              jfld(2)=0
           endif
           !
           !  subtract min value from spatial diff field
           !
           isd=idrstmpl(17)+1
           minsd=minval(jfld(isd:nonmiss))
           do j=isd,nonmiss
              jfld(j)=jfld(j)-minsd
           enddo
           !
           !   find num of bits need to store minsd and add 1 extra bit
           !   to indicate sign
           !
           temp=alog(real(abs(minsd)+1))/alog2
           nbitsd=ceiling(temp)+1
           !
           !   find num of bits need to store ifld(1) ( and ifld(2)
           !   if using 2nd order differencing )
           !
           maxorig=ival1
           if (idrstmpl(17).eq.2.and.ival2.gt.ival1) maxorig=ival2
           temp=alog(real(maxorig+1))/alog2
           nbitorig=ceiling(temp)+1
           if (nbitorig.gt.nbitsd) nbitsd=nbitorig
           !   increase number of bits to even multiple of 8 ( octet )
           if (mod(nbitsd,8).ne.0) nbitsd=nbitsd+(8-mod(nbitsd,8))
           !
           !  store extra spatial differencing info into the packed
           !  data section.
           !
           if (nbitsd.ne.0) then
              !   pack first original value
              if (ival1.ge.0) then
                 call sbyte(cpack,ival1,iofst,nbitsd)
                 iofst=iofst+nbitsd
              else
                 call sbyte(cpack,1,iofst,1)
                 iofst=iofst+1
                 call sbyte(cpack,iabs(ival1),iofst,nbitsd-1)
                 iofst=iofst+nbitsd-1
              endif
              if (idrstmpl(17).eq.2) then
               !  pack second original value
                 if (ival2.ge.0) then
                    call sbyte(cpack,ival2,iofst,nbitsd)
                    iofst=iofst+nbitsd
                 else
                    call sbyte(cpack,1,iofst,1)
                    iofst=iofst+1
                    call sbyte(cpack,iabs(ival2),iofst,nbitsd-1)
                    iofst=iofst+nbitsd-1
                 endif
              endif
              !  pack overall min of spatial differences
              if (minsd.ge.0) then
                 call sbyte(cpack,minsd,iofst,nbitsd)
                 iofst=iofst+nbitsd
              else
                 call sbyte(cpack,1,iofst,1)
                 iofst=iofst+1
                 call sbyte(cpack,iabs(minsd),iofst,nbitsd-1)
                 iofst=iofst+nbitsd-1
              endif
           endif
         !print *,'sdp ',ival1,ival2,minsd,nbitsd
        endif     !  end of spatial diff section
        !
        !  expand non-missing data values to original grid.
        !
        miss1=minval(jfld(1:nonmiss))-1
        miss2=miss1-1
        n=0
        do j=1,ndpts
           if ( ifldmiss(j).eq.0 ) then
              n=n+1
              ifld(j)=jfld(n)
           elseif ( ifldmiss(j).eq.1 ) then
              ifld(j)=miss1
           elseif ( ifldmiss(j).eq.2 ) then
              ifld(j)=miss2
           endif
        enddo
        !
        !   determine groups to be used.
        !
        if ( simple_alg ) then
           !  set group length to 10 :  calculate number of groups
           !  and length of last group
           ngroups=ndpts/10
           glen(1:ngroups)=10
           itemp=mod(ndpts,10)
           if (itemp.ne.0) then
              ngroups=ngroups+1
              glen(ngroups)=itemp
           endif
        else
           ! use dr. glahn's algorithm for determining grouping.
           !
           kfildo=6
           minpk=10 
           inc=1
           maxgrps=(ndpts/minpk)+1
           allocate(jmin(maxgrps))
           allocate(jmax(maxgrps))
           allocate(lbit(maxgrps))
           call pack_gp(kfildo,ifld,ndpts,missopt,minpk,inc,miss1,miss2,
     &                  jmin,jmax,lbit,glen,maxgrps,ngroups,ibit,jbit,
     &                  kbit,novref,lbitref,ier)
           !print *,'sagier = ',ier,ibit,jbit,kbit,novref,lbitref
           do ng=1,ngroups
              glen(ng)=glen(ng)+novref
           enddo
           deallocate(jmin)
           deallocate(jmax)
           deallocate(lbit)
        endif
        !  
        !  for each group, find the group's reference value (min)
        !  and the number of bits needed to hold the remaining values
        !
        n=1
        do ng=1,ngroups
           !  how many of each type?
           num0=count(ifldmiss(n:n+glen(ng)-1) .eq. 0)
           num1=count(ifldmiss(n:n+glen(ng)-1) .eq. 1)
           num2=count(ifldmiss(n:n+glen(ng)-1) .eq. 2)
           if ( num0.eq.0 ) then      ! all missing values
              if ( num1.eq.0 ) then       ! all secondary missing
                gref(ng)=-2
                gwidth(ng)=0
              elseif ( num2.eq.0 ) then       ! all primary missing
                gref(ng)=-1
                gwidth(ng)=0
              else                           ! both primary and secondary
                gref(ng)=0
                gwidth(ng)=1
              endif
           else                       ! contains some non-missing data
             !    find max and min values of group
             gref(ng)=huge(n)
             imax=-1*huge(n)
             j=n
             do lg=1,glen(ng)
                if ( ifldmiss(j).eq.0 ) then
                  if (ifld(j).lt.gref(ng)) gref(ng)=ifld(j) 
                  if (ifld(j).gt.imax) imax=ifld(j) 
                endif
                j=j+1
             enddo
             if (missopt.eq.1) imax=imax+1
             if (missopt.eq.2) imax=imax+2
             !   calc num of bits needed to hold data
             if ( gref(ng).ne.imax ) then
                temp=alog(real(imax-gref(ng)+1))/alog2
                gwidth(ng)=ceiling(temp)
             else
                gwidth(ng)=0
             endif
           endif
           !   subtract min from data
           j=n
           mtemp=2**gwidth(ng)
           do lg=1,glen(ng)
              if (ifldmiss(j).eq.0) then       ! non-missing
                 ifld(j)=ifld(j)-gref(ng)
              elseif (ifldmiss(j).eq.1) then    ! primary missing
                 ifld(j)=mtemp-1
              elseif (ifldmiss(j).eq.2) then    ! secondary missing
                 ifld(j)=mtemp-2
              endif
              j=j+1
           enddo
           !   increment fld array counter
           n=n+glen(ng)
        enddo
        !  
        !  find max of the group references and calc num of bits needed 
        !  to pack each groups reference value, then
        !  pack up group reference values
        !
        !write(77,*)'grefs: ',(gref(j),j=1,ngroups)
        igmax=maxval(gref(1:ngroups))
        if (missopt.eq.1) igmax=igmax+1
        if (missopt.eq.2) igmax=igmax+2
        if (igmax.ne.0) then
           temp=alog(real(igmax+1))/alog2
           nbitsgref=ceiling(temp)
           ! restet the ref values of any "missing only" groups.
           mtemp=2**nbitsgref
           do j=1,ngroups
               if (gref(j).eq.-1) gref(j)=mtemp-1
               if (gref(j).eq.-2) gref(j)=mtemp-2
           enddo
           call sbytes(cpack,gref,iofst,nbitsgref,0,ngroups)
           itemp=nbitsgref*ngroups
           iofst=iofst+itemp
           !         pad last octet with zeros, if necessary,
           if (mod(itemp,8).ne.0) then
              left=8-mod(itemp,8)
              call sbyte(cpack,zero,iofst,left)
              iofst=iofst+left
           endif
        else
           nbitsgref=0
        endif
        !
        !  find max/min of the group widths and calc num of bits needed
        !  to pack each groups width value, then
        !  pack up group width values
        !
        !write(77,*)'gwidths: ',(gwidth(j),j=1,ngroups)
        iwmax=maxval(gwidth(1:ngroups))
        ngwidthref=minval(gwidth(1:ngroups))
        if (iwmax.ne.ngwidthref) then
           temp=alog(real(iwmax-ngwidthref+1))/alog2
           nbitsgwidth=ceiling(temp)
           do i=1,ngroups
              gwidth(i)=gwidth(i)-ngwidthref
           enddo
           call sbytes(cpack,gwidth,iofst,nbitsgwidth,0,ngroups)
           itemp=nbitsgwidth*ngroups
           iofst=iofst+itemp
           !         pad last octet with zeros, if necessary,
           if (mod(itemp,8).ne.0) then
              left=8-mod(itemp,8)
              call sbyte(cpack,zero,iofst,left)
              iofst=iofst+left
           endif
        else
           nbitsgwidth=0
           gwidth(1:ngroups)=0
        endif
        !
        !  find max/min of the group lengths and calc num of bits needed
        !  to pack each groups length value, then
        !  pack up group length values
        !
        !write(77,*)'glens: ',(glen(j),j=1,ngroups)
        ilmax=maxval(glen(1:ngroups-1))
        nglenref=minval(glen(1:ngroups-1))
        nglenlast=glen(ngroups)
        if (ilmax.ne.nglenref) then
           temp=alog(real(ilmax-nglenref+1))/alog2
           nbitsglen=ceiling(temp)
           do i=1,ngroups-1
              glen(i)=glen(i)-nglenref
           enddo
           call sbytes(cpack,glen,iofst,nbitsglen,0,ngroups)
           itemp=nbitsglen*ngroups
           iofst=iofst+itemp
           !         pad last octet with zeros, if necessary,
           if (mod(itemp,8).ne.0) then
              left=8-mod(itemp,8)
              call sbyte(cpack,zero,iofst,left)
              iofst=iofst+left
           endif
        else
           nbitsglen=0
           glen(1:ngroups)=0
        endif
        !
        !  for each group, pack data values
        !
        !write(77,*)'iflds: ',(ifld(j),j=1,ndpts)
        n=1
        ij=0
        do ng=1,ngroups
           glength=glen(ng)+nglenref
           if (ng.eq.ngroups ) glength=nglenlast
           grpwidth=gwidth(ng)+ngwidthref
       !write(77,*)'ngp ',ng,grpwidth,glength,gref(ng)
           if ( grpwidth.ne.0 ) then
              call sbytes(cpack,ifld(n),iofst,grpwidth,0,glength)
              iofst=iofst+(grpwidth*glength)
           endif
           do kk=1,glength
              ij=ij+1
        !write(77,*)'sag ',ij,fld(ij),ifld(ij),gref(ng),bscale,rmin,dscale
           enddo
           n=n+glength
        enddo
        !         pad last octet with zeros, if necessary,
        if (mod(iofst,8).ne.0) then
           left=8-mod(iofst,8)
           call sbyte(cpack,zero,iofst,left)
           iofst=iofst+left
        endif
        lcpack=iofst/8
        !
        if ( allocated(ifld) ) deallocate(ifld)
        if ( allocated(jfld) ) deallocate(jfld)
        if ( allocated(ifldmiss) ) deallocate(ifldmiss)
        if ( allocated(gref) ) deallocate(gref)
        if ( allocated(gwidth) ) deallocate(gwidth)
        if ( allocated(glen) ) deallocate(glen)
      !else           !   constant field ( max = min )
      !  nbits=0
      !  lcpack=0
      !  nbitsgref=0
      !  ngroups=0
      !endif

!
!  fill in ref value and number of bits in template 5.2
!
      call mkieee(rmin,ref,1)   ! ensure reference value is ieee format
!      call gbyte(ref,idrstmpl(1),0,32)
      iref=transfer(ref,iref)
      idrstmpl(1)=iref
      idrstmpl(4)=nbitsgref
      idrstmpl(5)=0         ! original data were reals
      idrstmpl(6)=1         ! general group splitting
      idrstmpl(10)=ngroups          ! number of groups
      idrstmpl(11)=ngwidthref       ! reference for group widths
      idrstmpl(12)=nbitsgwidth      ! num bits used for group widths
      idrstmpl(13)=nglenref         ! reference for group lengths
      idrstmpl(14)=1                ! length increment for group lengths
      idrstmpl(15)=nglenlast        ! true length of last group
      idrstmpl(16)=nbitsglen        ! num bits used for group lengths
      if (idrsnum.eq.3) then
         idrstmpl(18)=nbitsd/8      ! num bits used for extra spatial
                                    ! differencing values
      endif

      return
      end
