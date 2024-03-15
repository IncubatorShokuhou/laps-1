!dis   
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis    
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis   
!dis 

module horiz_interp

  ! contains routines for doing horizontal interpolation
  implicit none
  private

  ! interp method codes for standard interpolation
  integer, public, parameter          :: method_nearest = 0
  integer, public, parameter          :: method_linear  = 1
  integer, public, parameter          :: method_higher = 2
  
  ! data type codes for masked interpolation
  integer, public, parameter          :: flag_field = 0
  integer, public, parameter          :: category_field = 1
  integer, public, parameter          :: value_field = 2
  public interpolate_standard
  public interpolate_masked_val

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine interpolate_standard(nx_in, ny_in, data_in, &
                                  nx_out, ny_out, xloc, yloc, &
                                  method, data_out)
    !
    ! performs standard interpolation (i.e., non-masked fields)
    implicit none

    ! dimensions and grid we are interpolating from
    integer, intent(in)        :: nx_in
    integer, intent(in)        :: ny_in
    real, intent(in)           :: data_in(nx_in,ny_in)

    ! dimensions of output grid
    integer, intent(in)        :: nx_out
    integer, intent(in)        :: ny_out

    ! pre-computed real i/j indices of input grid
    ! for each point in output grid
    real, intent(in)           :: xloc(nx_out,ny_out)
    real, intent(in)           :: yloc(nx_out,ny_out)

    ! code indicating interp method (see module top)
    integer, intent(in)        :: method
   
    ! output interpolated data
    real, intent(out)          :: data_out(nx_out,ny_out)
   
    ! locals
    integer                    :: i
    integer                    :: j
    integer                    :: ilo, jlo
    real                       :: dlon, dlat
    real                       :: deltalat, deltalon

    select case(method)
      case(method_nearest)
        do j = 1, ny_out
          do i = 1, nx_out
            ilo = nint(xloc(i,j))
            jlo = nint(yloc(i,j))
            data_out(i,j) = data_in( ilo, jlo )
          enddo
        enddo

      case(method_linear)
        do j = 1, ny_out
          do i = 1, nx_out
            ilo = min(floor(xloc(i,j)), nx_in-1)
            jlo = min(floor(yloc(i,j)), ny_in-1)
            dlon = xloc(i,j) - ilo
            dlat = yloc(i,j) - jlo
            deltalat = 1.
            deltalon = 1.
            data_out(i,j) = ((deltalon-dlon)*((deltalat-dlat)*data_in(ilo,jlo) &
                          + dlat*data_in(ilo,jlo+1))   &
                          + dlon*((deltalat-dlat)*data_in(ilo+1,jlo)     &
                          + dlat*data_in(ilo+1,jlo+1))) &
                          / ( deltalat * deltalon ) 
          enddo
        enddo
      case(method_higher) 
        do j = 1, ny_out
          do i = 1, nx_out
            data_out(i,j) = bint( xloc(i,j), yloc(i,j), data_in, nx_in, ny_in, 0)
          enddo
        enddo
      case default
        print '(a,i2)', 'uknown interpolation method code: ', method
        stop 'interpolate_standard'
    end select
    return
  end subroutine interpolate_standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine interpolate_masked_val(nx_src, ny_src, lmask_src, data_src, &
                                 nx_out, ny_out, lmask_out, data_out, &
                                 isrcr, jsrcr, make_srcmask, &
                                 min_val, max_val, def_val, val_mask, &
                                 method)

    implicit none
    
    integer, intent(in)           :: nx_src
    integer, intent(in)           :: ny_src
    real, intent(inout)           :: lmask_src(nx_src,ny_src)
    real, intent(in)              :: data_src(nx_src,ny_src)
    integer, intent(in)           :: nx_out, ny_out
    real, intent(in)              :: lmask_out(nx_out,ny_out)
    real, intent(inout)             :: data_out(nx_out,ny_out)
    real, intent(in)              :: isrcr(nx_out,ny_out)
    real, intent(in)              :: jsrcr(nx_out,ny_out)
    logical, intent(in)           :: make_srcmask
    real, intent(in)              :: min_val, max_val, def_val, val_mask
    integer, intent(in)           :: method

    integer                       :: i,j,k,l,ilo,jlo, ii, jj, kk, ll
    integer                       :: search_rad
    integer                       :: search_rad_max
    real, parameter               :: bad_point = -99999.
    logical                       :: bad_points_flag
    logical                       :: fixed_bad
    integer                       :: total_points
    integer                       :: points_val_mask
    integer                       :: points_skipped
    integer                       :: num_to_fix
    integer                       :: fixed_with_out
    integer                       :: fixed_with_src
    integer                       :: fixed_with_def
    integer                       :: close_pts
    integer                       :: ic,jc,iic,jjc
    real                          :: deltagx, deltagy
    real                          :: dix, djy
    real                          :: data_close(4)
    real                          :: distance(4)
    real                          :: distx,disty
    real                          :: sum_dist
    real                          :: a,b,c,d,e,f,g,h
    real                          :: stl(4,4)
    real                          :: valb
    real                          :: srcmask(nx_src,ny_src)

    ! can we use the source data land mask that was input, or do we need to make it
    ! from the min_val/max_val arguments?  

    if (make_srcmask) then
      print '(a)', 'interpolate_masked_val: making source data land mask...'
    
      ! this code has the flexibility to handled either water data points (lmask =0)
      ! or land points (lmask = 1). so if we are making the landmask using valid 
      ! range, but we do not know anything about the variable other than the 
      ! valid mask value, we need to figure out which points are 0 and which should be
      ! 1.  to do this, initialize the array to the invalid value, which we determine
      ! from the valid value.
      if (val_mask .eq. 1) then
        lmask_src(:,:) = 0
      else
        lmask_src(:,:) = 1
      endif
      
      ! now figure out where to set the mask using a where statement
      where((data_src .ge. min_val).and.(data_src .le. max_val)) lmask_src = val_mask
    else
      print '(a)', 'interpolate_masked: using source landmask field.'
    endif

   
    bad_points_flag = .false.
    
    ! initialize counters
    total_points = 0
    points_val_mask = 0
    points_skipped = 0
    num_to_fix = 0
    fixed_with_out = 0
    fixed_with_src = 0
    fixed_with_def = 0

    ! select interpolation method.  putting the case statement out here
    ! increases the amount of replicated code but should be more efficient
    ! than checking this condition at every grid point.

    select case(method)
  
      case(method_nearest)
        ! use nearest neigbor
        print '(a)', 'masked interpolation using nearest neighbor value...'
        out_j_loop_1: do j = 1, ny_out
          out_i_loop_1: do i = 1, nx_out

            total_points = total_points + 1
            ! we only need to process this point if the lmask_out is equal
            ! to val_mask.  for example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  any output point that is water
            ! (lmask = 0) will then be skipped.  during this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.

            if (lmask_out(i,j) .eq. val_mask) then
              ! process this point 
              points_val_mask = points_val_mask + 1
              ilo = nint(isrcr(i,j))
              jlo = nint(jsrcr(i,j))
        
              ! see if this point can be used
              if (lmask_src(ilo,jlo).eq. lmask_out(i,j)) then
                data_out(i,j) = data_src(ilo,jlo)
              else
                data_out(i,j) = bad_point
                num_to_fix = num_to_fix + 1
                bad_points_flag = .true.
              endif 
            else
              ! the output grid does not require a value for this point
              ! but do not zero out in case this is a field begin
              ! done twice (once for water and once for land, e.g.
              ! skintemp
              !data_out(i,j) = 0.
              points_skipped = points_skipped + 1
            endif
          enddo out_i_loop_1
        enddo out_j_loop_1
      
      case (method_linear)
        ! use a 4-point interpolation
        print '(a)', 'masked interpolation using 4-pt linear interpolation...'
        deltagx = 1.
        deltagy = 1.
        out_j_loop_2: do j = 1, ny_out
          out_i_loop_2: do i = 1, nx_out

            total_points = total_points + 1
            ! we only need to process this point if the lmask_out is equal
            ! to val_mask.  for example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  any output point that is water
            ! (lmask = 0) will then be skipped.  during this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.

            if (lmask_out(i,j) .eq. val_mask) then
              ! process this point
              points_val_mask = points_val_mask + 1
              ilo = min(floor(isrcr(i,j)),nx_src-1)
              jlo = min(floor(jsrcr(i,j)),ny_src-1)
              dix = isrcr(i,j) - float(ilo)
              djy = jsrcr(i,j) - float(jlo)
 
              ! loop around the four surrounding points
              ! and count up the number of points we can use based
              ! on common mask value
              close_pts = 0
              sum_dist = 0.
              outer_four_j: do jc = 0,1
                outer_four_i: do ic = 0,1
                  iic = ilo + ic
                  jjc = jlo + jc
                  if (lmask_src(iic,jjc).eq. lmask_out(i,j)) then
                    close_pts = close_pts + 1
                    data_close(close_pts) = data_src(iic,jjc)
                       
                    ! compute distance to this valid point
                    ! in grid units and add to sum of distances (actually,
                    ! we are doing a weight, which is inversely proportional
                    ! to distance)
                    if (ic .eq. 0) then
                      distx = deltagx - dix
                    else
                      distx =  dix
                    endif
                    if (jc .eq. 0) then
                      disty = deltagy - djy
                    else
                      disty = djy
                    endif  
                    distance(close_pts) = sqrt(distx**2+disty**2)
                    sum_dist = sum_dist + distance(close_pts)
                  endif 
                enddo outer_four_i
              enddo outer_four_j
 
              ! did we find at least one point in the surrounding four 
              ! that was usable?

              if (close_pts .gt. 0) then
               
                ! if we have all four points, then do bilinear interpolation
                if (close_pts .eq. 4) then
                   data_out(i,j) = ((deltagx - dix)*((deltagy-djy)*data_close(1) &
                                 + djy*data_close(3)) &
                                 + dix*((deltagy-djy)*data_close(2) &
                                 + djy*data_close(4))) &
                                 / (deltagx * deltagy)
                else if ((close_pts .gt. 1).and.(close_pts .lt. 4)) then

                  ! simple distance-weighted average by computing
                  ! the sum of all distances to each point and using
                  ! each individual distance divided by the total
                  ! distance as the weighting

                  data_out(i,j) = 0.
                  do k = 1, close_pts
                    data_out(i,j) = data_out(i,j) + &
                                    (distance(k)/sum_dist) * data_close(k)
                  enddo
                else
                  ! set output value = to one point we found
                  data_out(i,j) = data_close(1)
                endif
              else
                bad_points_flag = .true.
                data_out(i,j) = bad_point
                num_to_fix = num_to_fix + 1   
              endif 
            else
              ! the output grid does not require a value for this point
              ! but do not zero out in case this is a field begin                             ! done twice (once for water and once for land, e.g.                            ! skintemp
              !data_out(i,j) = 0.
              points_skipped = points_skipped + 1
            endif
          enddo out_i_loop_2
        enddo out_j_loop_2                                        

      case (method_higher)
        ! 16-point interpolation 
        out_j_loop_3: do j = 1, ny_out
          out_i_loop_3: do i = 1, nx_out

            total_points = total_points + 1
            ! we only need to process this point if the lmask_out is equal
            ! to val_mask.  for example, if one is processing soil parameters,
            ! which are only valid at land points (lmask = 1), then the user
            ! passes in 1. for the val_mask.  any output point that is water
            ! (lmask = 0) will then be skipped.  during this loop, if we have
            ! points that cannot be properly assigned, we will mark them as bad
            ! and set the bad_points_lag.

            if (lmask_out(i,j) .eq. val_mask) then
              ! process this point
              points_val_mask = points_val_mask + 1   
 
              ! do a 4x4 loop around in the input data around the output
              ! point to get our 16 points of influence.  only use those
              ! that are appropriately masked
              valb = 0.
              close_pts = 0
              ilo = int(isrcr(i,j)+0.00001)
              jlo = int(jsrcr(i,j)+0.00001)
              dix = isrcr(i,j) - ilo
              djy = jsrcr(i,j) - jlo
              if ( (abs(dix).gt.0.0001).or.(abs(djy).gt.0.0001) ) then
                ! do the interpolation loop
                stl(:,:) = 0.
                loop_16_1: do k = 1,4
                  kk = ilo + k - 2
                  if ((kk .lt. 1).or.(kk .gt. nx_src)) cycle loop_16_1
                  loop_16_2: do l = 1, 4
                    ll = jlo + l - 2
                    if ((ll .lt. 1).or.(ll .gt. ny_src)) cycle loop_16_2

                    ! check land mask at this source point
                    if (lmask_src(kk,ll).ne. val_mask) cycle loop_16_2
                    ! if we are here, then mask tests passed
                    stl(k,l) = data_src(kk,ll) 
                    if ( (stl(k,l) .eq. 0.).and.(min_val.le.0.).and. &
                                                (max_val.ge.0.) ) then
                      stl = 1.e-5
                    endif
                    close_pts = close_pts + 1
                  enddo loop_16_2
                enddo loop_16_1
  
                ! did we find any valid points?

                if ( (close_pts .gt. 0).and. ( &
                  (stl(2,2).gt.0.).and.(stl(2,3).gt.0.).and. &
                  (stl(3,2).gt.0.).and.(stl(3,3).gt.0.)  ) ) then
                  a = oned(dix,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
                  b = oned(dix,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
                  c = oned(dix,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
                  d = oned(dix,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
                  valb = oned(djy,a,b,c,d)
                  if (close_pts .ne. 16) then
                    e = oned(djy,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
                    f = oned(djy,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
                    g = oned(djy,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
                    h = oned(djy,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
                    valb = (valb+oned(dix,e,f,g,h)) * 0.5
                  endif
                  data_out(i,j) = valb
            
                else
                  bad_points_flag = .true.
                  data_out(i,j) = bad_point
                  num_to_fix = num_to_fix + 1
                endif
              else
                ! we are right on a source point, so try to use it
                if (lmask_src(ilo,jlo).eq.val_mask) then
                  data_out(i,j) = data_src(ilo,jlo)
                else
                  bad_points_flag = .true.
                  data_out(i,j) = bad_point
                  num_to_fix = num_to_fix + 1 
                endif
              endif
            else
              ! the output grid does not require a value for this point
              !data_out(i,j) = 0.
              points_skipped = points_skipped + 1      
            endif
          enddo out_i_loop_3
        enddo out_j_loop_3        

    end select

     ! do we need to correct bad points?
  
    if (bad_points_flag) then
     
      search_rad_max = 10
      fix_bad_j: do j = 1, ny_out
        fix_bad_i: do i = 1, nx_out
 
          if (data_out(i,j).ne. bad_point) cycle fix_bad_i
          
          ! first, search for nearest non-bad point in the output domain
          ! which is usually higher resolution. 
          fixed_bad = .false.
          search_out_loop: do search_rad = 1, search_rad_max
            search_out_j: do ll = -(search_rad-1), (search_rad-1),1
              jj = j + ll
              if ((jj .lt. 1).or.(jj .gt. ny_out)) cycle search_out_j
              search_out_i: do kk = -(search_rad), search_rad, 1
                 ii = i + kk
                 if ((ii .lt. 1).or.(ii .gt. nx_out)) cycle search_out_j
                 if ((data_out(ii,jj).ne.bad_point).and. &
                    (lmask_out(ii,jj) .eq. val_mask) ) then
                  data_out(i,j) = data_out(ii,jj)
                  fixed_bad = .true.
                  fixed_with_out = fixed_with_out + 1
                  exit search_out_loop
                endif
              enddo search_out_i
            enddo search_out_j
          enddo search_out_loop

          ! did we fix the point?  if not, then do same search on src data.
          if (.not. fixed_bad) then
            search_rad_max = 10
            search_src_loop: do search_rad = 1, search_rad_max
              search_src_j: do ll = -(search_rad-1), (search_rad-1),1
                jj = nint(jsrcr(i,j)) + ll
                if ((jj .lt. 1).or.(jj .gt. ny_src)) cycle search_src_j
                search_src_i: do kk = -(search_rad), search_rad, 1
                   ii = nint(isrcr(i,j)) + kk
                   if ((ii .lt. 1).or.(ii .gt. nx_src)) cycle search_src_j
                   if (lmask_src(ii,jj).eq.val_mask) then
                     data_out(i,j) = data_src(ii,jj)
                     fixed_bad = .true.
                     fixed_with_src = fixed_with_src + 1
                     exit search_src_loop
                   endif
                enddo search_src_i
              enddo search_src_j
            enddo search_src_loop
          endif
          ! now is the point fixed?  if not, we have to use a default value.
          if (.not.fixed_bad) then
            fixed_with_def = fixed_with_def + 1
            data_out(i,j) = def_val
            print '(a,f10.3,a,2i5)', 'interpolate_masked: bogus value of ', def_val, &
                ' used at point ', i, j
          endif
         
        enddo fix_bad_i
      enddo fix_bad_j
    endif
    print '(a)',     '----------------------------------------'
    print '(a)',     'masked interpolation summary: '
    print '(a,i10)', '  total points in grid:       ', total_points
    print '(a,i10)', '  points needing values:      ', points_val_mask
    print '(a,i10)', '  points not required:        ', points_skipped
    print '(a,i10)', '  points needing fix:         ', num_to_fix
    print '(a,i10)', '  points fixed with out grid: ', fixed_with_out
    print '(a,i10)', '  points fixed with src grid: ', fixed_with_src
    print '(a,i10)', '  points fixed with def val:  ', fixed_with_def
    print '(a)',     '----------------------------------------'
    return
  end subroutine interpolate_masked_val        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function bint(xx,yy,list,iii,jjj,ibint)
   
      !  16 point interpolation
   
      implicit none
   
      real :: xx , yy
      integer :: ibint , iii, jjj
      real list(iii,jjj),stl(4,4)
   
      integer :: ib , jb, n , i , j , k , kk , l , ll
      real :: bint , x , y , a , b , c , d , e , f , g , h
      ib=iii-ibint
      jb=jjj-ibint
      bint = 0.0
      n = 0
      i = int(xx+0.00001)
      j = int(yy+0.00001)
      x = xx - i
      y = yy-j
     
      if ( ( abs(x).gt.0.0001 ) .or. ( abs(y).gt.0.0001 ) ) then
         stl(:,:)=1.e-10      
         loop_1 : do k = 1,4
            kk = i + k - 2
            if ( ( kk .lt. 1) .or. ( kk .gt. ib ) ) then
               cycle loop_1
            end if
            loop_2 : do l = 1,4
               stl(k,l) = 0.
               ll = j + l - 2
               if ( ( ll .gt. jb ) .or. ( ll .lt. 1 ) ) then
                  cycle loop_2
               end if
               stl(k,l) = list(kk,ll)
               n = n + 1
               if ( stl(k,l) .eq. 0. ) then
                  stl(k,l) = 1.e-20
               end if
            end do loop_2
         end do loop_1
         a = oned(x,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
         b = oned(x,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
         c = oned(x,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
         d = oned(x,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
         bint = oned(y,a,b,c,d)
   
         if(n.ne.16) then
            e = oned(y,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
            f = oned(y,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
            g = oned(y,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
            h = oned(y,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
            bint = (bint+oned(x,e,f,g,h)) * 0.5
         end if
   
      else
         bint = list(i,j)
      end if
   end function bint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function oned(x,a,b,c,d) 
   
      implicit none
   
      real :: x,a,b,c,d,oned
   
      oned = 0.                
   
      if      ( x .eq. 0. ) then
         oned = b      
      else if ( x .eq. 1. ) then
         oned = c      
      end if
   
      if(b*c.ne.0.) then
         if ( a*d .eq. 0. ) then
            if      ( ( a .eq. 0 ) .and. ( d .eq. 0 ) ) then
               oned = b*(1.0-x)+c*x                                        
            else if ( a .ne. 0. ) then
               oned = b+x*(0.5*(c-a)+x*(0.5*(c+a)-b))            
            else if ( d .ne. 0. ) then
               oned = c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)) 
            end if
         else
            oned = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))+ &
                   x*(c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)))                                   
         end if
      end if
   
   end function oned 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                     
end module horiz_interp
