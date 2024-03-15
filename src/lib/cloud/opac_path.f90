
      subroutine opac_wt(s,extc,field,ns,opac_thr,field_wt)
!                        i  i     i    i    i        o

!     calculate field weighted by the opacity along the line of sight

      implicit real*8 (a-h,o-z)

      dimension s(ns),extc(ns),field(ns),od(ns),opac(ns)

!     'opac_thr' isn't yet used, though it could be if needed for efficiency to
!     stop looping earlier      

      call opac_path(s,extc,ns,od,opac)
!                    i  i   i   o   o

      sum_field_wt = 0d0

      do i = 2,ns      
          il = i-1
          d_opac = opac(i) - opac(il)
          field_ave = 0.5d0 * (field(i) + field(il))
          sum_fieldwt = sum_fieldwt + d_opac * field_ave
          sum_opac = opac(i)
      enddo ! i

      if(sum_opac .gt. 0.)then
          field_wt = sum_field_wt / sum_opac
      else
          write(6,'(" error in opac_wt")')
          stop
      endif

      return
      end

      subroutine opac_path(s,extc,ns,od,opac)
!                          i  i   i   o   o

!     calculate optical depth and opacity along a path. if a solid surface
!     such as the ground is in the path a high value of 'extc' is used along
!     with two identical values of s.

      implicit real*8 (a-h,o-z)

      parameter (extc_large=1d10)
      parameter (od_large=1d10)

      trans(x) = exp(-min(x,80d0))
      opacity(x) = 1d0 - trans(x)

      dimension s(ns),extc(ns),od(ns),opac(ns)

      od(1) = 0.
      opac(1) = 0.
      odsum = 0d0

      do i = 2,ns      
          il = i-1
          ds = s(i) - s(il)
          if(extc(i) .lt. extc_large)then ! not the ground
              extc_ave = 0.5d0 * (extc(i) + extc(il))
              od(i) = extc_ave * ds
              odsum = odsum + od(il)
          else                            ! hitting the ground
              od(i) = od_large
              opac(i) = opacity(od(i))
          endif
      enddo ! i
         
      return
      end
      
