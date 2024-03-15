  subroutine compute_tke(p3d, t3d, u3d, v3d, z3d, topo, &
                         nx, ny, nz, tke3d)

     implicit none

     integer, intent(in)  :: nx, ny, nz
     real, intent(in)     :: p3d(nx, ny, nz)
     real, intent(in)     :: t3d(nx, ny, nz)
     real, intent(in)     :: u3d(nx, ny, nz)
     real, intent(in)     :: v3d(nx, ny, nz)
     real, intent(in)     :: z3d(nx, ny, nz)
     real, intent(in)     :: topo(nx, ny)
     real, intent(out)    :: tke3d(nx, ny, nz)

     integer :: i, j
     real    :: tke1d(nz)

     print '(a)', 'entering compute_tke'
     ! initialize the tke array

     tke3d(:, :, :) = 0.

     ! loop over horizontal to compute tke on a column by
     ! column basis
     do j = 1, ny
        do i = 1, nx

           ! use compute_dtf3
           call compute_dtf3(p3d(i, j, :), t3d(i, j, :), u3d(i, j, :), &
                             v3d(i, j, :), z3d(i, j, :), topo(i, j), &
                             tke1d, nz)
           tke3d(i, j, :) = tke1d

        end do
     end do

     print '(a,2f8.2)', 'min/max computed tke = ', minval(tke3d), &
        maxval(tke3d)
  end subroutine compute_tke

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_dtf3(p, t, u, v, z, zter, tke_kh, nz)

!
! adrian marroquin fsl
! version modified for itfa
! 11/09/98
! modified for fortran 90 free-form
!
!-------------------------------------------------------------------------
     use constants
     implicit none
     integer, intent(in)    :: nz
     real, intent(in)       :: p(nz)
     real, intent(in)       :: t(nz)
     real, intent(in)       :: u(nz)
     real, intent(in)       :: v(nz)
     real, intent(in)       :: z(nz)
     real, intent(in)       :: zter
     real, intent(out)      :: tke_kh(nz)

     real :: c1, c2, c3, c13, c23, ce, alinf, akarm, cr, &
             cepn, cepp, iepn3, akm, prands
     real, dimension(nz):: brnt, shr, ri, epsilon

     integer :: klev, k
     real    :: pi1, pi2, th1, th2, pi3, th3, brunt, &
                shru, shrv, beta, ztop, zsfc, zlev, &
                rff, br, dz, alb, als, all

     real, external :: vertirreg, rf, rfkondo
     ! constants from stull (1988), page 219

     c1 = 1.44
     c2 = 1.00
     c3 = 1.92
     c13 = c1/c3
     c23 = c2/c3
     ce = 0.19
     alinf = 200.
     akarm = 0.35
     cr = 0.54
     cepn = 2.5
     cepp = 0.76
     iepn3 = 0
     akm = 75.0
!
!--------------------------------------------------------------
! compute ri, brnt, and shr
!
!
     prands = 2.5
     klev = nz
     pi1 = (p(1)/100000.)**kappa
     pi2 = (p(2)/100000.)**kappa
     th1 = t(1)/pi1
     th2 = t(2)/pi2

     do k = 2, nz - 1
        pi3 = (p(k + 1)/100000.)**kappa
        th3 = t(k + 1)/pi3
!
! vertical derivatives using pressure
!
        brunt = vertirreg(th1, th2, th3, &
                          p(k - 1), p(k), p(k + 1))
        shru = vertirreg(u(k - 1), u(k), u(k + 1), &
                         p(k - 1), p(k), p(k + 1))

        shrv = vertirreg(v(k - 1), v(k), v(k + 1), &
                         p(k - 1), p(k), p(k + 1))

        beta = grav*grav*p(k)/(r*pi2*th2*th2)
        brnt(k) = -beta*brunt
        shr(k) = beta*p(k)*(shru*shru + shrv*shrv)/(r*pi2)
        th1 = th2
        th2 = th3
        pi1 = pi2
        pi2 = pi3
        ri(k) = brnt(k)/(shr(k) + 1.e-10)
     end do

     brnt(1) = brnt(2)
     shr(1) = shr(2)
     ri(1) = ri(2)
     brnt(nz) = brnt(nz - 1)
     shr(nz) = shr(nz - 1)

     ri(nz) = ri(nz - 1)

     do k = 1, nz
        if (ri(k) .gt. 120.) ri(k) = 120.
     end do
!
!--------------------------------------------------------------
!
! now compute dissipation
! ztop equivalent to cpbl
!
     ztop = zter + 3000.
     zsfc = zter
!
     do k = 1, klev
!
        zlev = z(k)
        if (ri(k) .gt. 0.01) then
           rff = rfkondo(ri(k))
        else
           rff = rf(ri(k))
        end if
        epsilon(k) = akm*shr(k)*(c13 - c23*rff)
        if (epsilon(k) .lt. 0.) epsilon(k) = 0.
        if (iepn3 .eq. 0) then                 ! if iepn3 = 1, only epn
           if (brnt(k) .le. 0.) then
              tke_kh(k) = 0.
           else
              br = sqrt(brnt(k))
              tke_kh(k) = 0.7*epsilon(k)/(ce*br)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (zsfc .le. zlev .and. zlev .le. ztop) then                     !
                 dz = zlev - zsfc                                           !
                 alb = alinf*akarm*dz/(akarm*dz + alinf)                      !
                 als = cr*sqrt(tke_kh(k))/br                                !
                 all = amin1(alb, als)                                       !
                 tke_kh(k) = (all*epsilon(k))**.666666                      !
              end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           end if
        end if

     end do    ! end of k-loop (klev)

     return
  end subroutine compute_dtf3
!---------------------------------------------------------------------
  function rfkondo(ri)

     implicit none
     real, parameter  :: c0 = 6.873
     real, parameter  :: c1 = 7.
     real             :: d1, ahm, ri, rfkondo
!
! rfc (critical flux ri) = 0.143
!
     if (ri .gt. 1.) then
        rfkondo = 1./c1
     else
        if (0.01 .lt. ri .and. ri .le. 1.) then
           d1 = 1.+c0*ri
           ahm = 1./(c0*ri + 1./d1)
           rfkondo = ri*ahm
        end if
     end if
!
! for ri < 0.01 use rf (yamada form)
!
     return
  end
!-------------------------------------------------------------------------
  function vertirreg(f1, f2, f3, x1, x2, x3)
     implicit none
     real :: f1, f2, f3, x1, x2, x3
     real :: dx1, dx2, rat1, rat2, sdx
     real :: vertirreg

     dx1 = x2 - x1
     dx2 = x3 - x2
     rat1 = dx1/dx2
     rat2 = 1./rat1
     sdx = 1./(dx1 + dx2)
     vertirreg = ((f3 - f2)*rat1 + (f2 - f1)*rat2)*sdx
     return
  end
!--------------------------------------------------------------------------
  function rf(ri)

     implicit none
     real :: ri, rf
     real, parameter :: c1 = 0.056
     real, parameter :: c2 = 0.300
     real, parameter :: c3 = 0.333
     real, parameter :: a1 = 0.780
     real, parameter :: a2 = 0.790
     real, parameter :: b1 = 15.0
     real, parameter :: b2 = 8.0
!
     real :: e1, e2, e3, e4, e5, f1, f2, f3, f4, f42

     e1 = b1 - 6.*a1
     e2 = b1 + 12.*a1*(1.-c2) + 3.*b2*(1.-c3)
     e3 = b1*(1.-3.*c1) - 6.*a1
     e4 = b1*(1.-3.*c1) + 12.*a1*(1.-c2) + 9.*a2*(1.-c2)
     e5 = b1 + 3.*a1*(1.-c2) + 3.*b2*(1.-c3)

     f1 = 0.5*a2*e5/(a1*e4)
     f2 = a1*e3/(a2*e5)
     f3 = 2.*a1*(e3*e5 - 2.*e1*e4)/(a2*e5*e5)
     f4 = a1*e3/(a2*e5)
     f42 = f4*f4

     rf = f1*(ri + f2 - sqrt(ri*ri + f3*ri + f42))

     return
  end
