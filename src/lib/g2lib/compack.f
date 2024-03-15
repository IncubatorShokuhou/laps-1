      subroutine compack(fld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    compack
!   prgmmr: gilbert          org: w/np11    date: 2000-06-21
!
! abstract: this subroutine packs up a data field using a complex
!   packing algorithm as defined in the grib2 documention.  it
!   supports grib2 complex packing templates with or without
!   spatial differences (i.e. drts 5.2 and 5.3).
!   it also fills in grib2 data representation template 5.2 or 5.3 
!   with the appropriate values.
!
! program history log:
! 2000-06-21  gilbert
!
! usage:    call compack(fld,ndpts,idrsnum,idrstmpl,cpack,lcpack)
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
!                (1) = reference value - set by compack routine.
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
      integer,allocatable :: ifld(:)
      integer,allocatable :: jmin(:),jmax(:),lbit(:)
      integer,parameter :: zero=0
      integer,allocatable :: gref(:),gwidth(:),glen(:)
      integer :: glength,grpwidth
      logical :: simple_alg = .false.
      
      alog2=alog(2.0)
      bscale=2.0**real(-idrstmpl(2))
      dscale=10.0**real(idrstmpl(3))
!
!  find max and min values in the data
!
      rmax=fld(1)
      rmin=fld(1)
      do j=2,ndpts
        if (fld(j).gt.rmax) rmax=fld(j)
        if (fld(j).lt.rmin) rmin=fld(j)
      enddo
!
!  if max and min values are not equal, pack up field.
!  if they are equal, we have a constant field, and the reference
!  value (rmin) is the value for each point in the field and
!  set nbits to 0.
!
      if (rmin.ne.rmax) then
        iofst=0
        allocate(ifld(ndpts))
        allocate(gref(ndpts))
        allocate(gwidth(ndpts))
        allocate(glen(ndpts))
        !
        !  scale original data
        !
        if (idrstmpl(2).eq.0) then        !  no binary scaling
           imin=nint(rmin*dscale)
           !imax=nint(rmax*dscale)
           rmin=real(imin)
           do j=1,ndpts
             ifld(j)=nint(fld(j)*dscale)-imin
           enddo
        else                              !  use binary scaling factor
           rmin=rmin*dscale
           !rmax=rmax*dscale
           do j=1,ndpts
             ifld(j)=nint(((fld(j)*dscale)-rmin)*bscale)
           enddo
        endif
        !
        !  calculate spatial differences, if using drs template 5.3
        !
        if (idrsnum.eq.3) then        ! spatial differences
           if (idrstmpl(17).ne.1.and.idrstmpl(17).ne.2) idrstmpl(17)=2
           if (idrstmpl(17).eq.1) then      ! first order
              ival1=ifld(1)
              do j=ndpts,2,-1
                 ifld(j)=ifld(j)-ifld(j-1)
              enddo
              ifld(1)=0
           elseif (idrstmpl(17).eq.2) then      ! second order
              ival1=ifld(1)
              ival2=ifld(2)
              do j=ndpts,3,-1
                 ifld(j)=ifld(j)-(2*ifld(j-1))+ifld(j-2)
              enddo
              ifld(1)=0
              ifld(2)=0
           endif
           !
           !  subtract min value from spatial diff field
           !
           isd=idrstmpl(17)+1
           minsd=minval(ifld(isd:ndpts))
           do j=isd,ndpts
              ifld(j)=ifld(j)-minsd
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
           missopt=0
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
        !  for each group, find the group's reference value
        !  and the number of bits needed to hold the remaining values
        !
        n=1
        do ng=1,ngroups
           !    find max and min values of group
           gref(ng)=ifld(n)
           imax=ifld(n)
           j=n+1
           do lg=2,glen(ng)
              if (ifld(j).lt.gref(ng)) gref(ng)=ifld(j) 
              if (ifld(j).gt.imax) imax=ifld(j) 
              j=j+1
           enddo
           !   calc num of bits needed to hold data
           if ( gref(ng).ne.imax ) then
              temp=alog(real(imax-gref(ng)+1))/alog2
              gwidth(ng)=ceiling(temp)
           else
              gwidth(ng)=0
           endif
           !   subtract min from data
           j=n
           do lg=1,glen(ng)
              ifld(j)=ifld(j)-gref(ng)
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
        if (igmax.ne.0) then
           temp=alog(real(igmax+1))/alog2
           nbitsgref=ceiling(temp)
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
        if ( allocated(gref) ) deallocate(gref)
        if ( allocated(gwidth) ) deallocate(gwidth)
        if ( allocated(glen) ) deallocate(glen)
      else           !   constant field ( max = min )
        nbits=0
        lcpack=0
        nbitsgref=0
        ngroups=0
      endif

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
      idrstmpl(7)=0         ! no internal missing values
      idrstmpl(8)=0         ! primary missing value
      idrstmpl(9)=0         ! secondary missing value
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
