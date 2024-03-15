
      subroutine bilinear_fill(t,imax,jmax,nskip,r_missing_data)

      integer lowi_lut(imax)
      integer lowj_lut(jmax)

      dimension t(imax,jmax)

!     bilinearly interpolate to fill in rest of domain
!     fills in final analysis value and weights from obs alone
!     we may have to extrapolate at the n and e edges
      do i = 1,imax
          lowi_lut(i) = ((i-1)/nskip)*nskip + 1
          il = lowi_lut(i)
          ih = il + nskip
          if(ih .gt. imax)lowi_lut(i) = lowi_lut(i) - nskip
!         write(6,*)' i,il,ih',i,il,ih
      enddo ! i
      do j = 1,jmax
          lowj_lut(j) = ((j-1)/nskip)*nskip + 1
          jl = lowj_lut(j)
          jh = jl + nskip
          if(jh .gt. jmax)lowj_lut(j) = lowj_lut(j) - nskip
      enddo ! i

      do j=1,jmax
          jl = lowj_lut(j)
          jh = jl + nskip
          fracj = dble(j-jl)/dble(nskip)

          do i=1,imax

              il = lowi_lut(i)
              ih = il + nskip
              fraci = dble(i-il)/dble(nskip)

!             write(6,*)' i,il,ih',i,il,ih
!             write(6,*)i,j,il,ih,jl,jh

!             calculate interpolated cloud cover
              z1=t(il,jl)
              z2=t(ih,jl)
              z3=t(ih,jh)
              z4=t(il,jh)

              t(i,j) = z1+(z2-z1)*fraci+(z4-z1)*fracj &
                     - (z2+z4-z3-z1)*fraci*fracj

          enddo ! i
      enddo ! j

      return
      end

