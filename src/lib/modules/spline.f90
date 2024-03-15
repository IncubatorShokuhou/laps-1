
! downloaded by steve albers from an iowa state math dept. website produced 
! by john burkhardt: 

! http://orion.math.iastate.edu/burkardt/f_src/spline/spline.html

subroutine basis_function_b_val ( tdata, tval, yval )
!
!*******************************************************************************
!
!! basis_function_b_val evaluates the b spline basis function.
!
!
!  discussion:
!
!    the b spline basis function is a piecewise cubic which
!    has the properties that:
!
!    * it equals 2/3 at tdata(3), 1/6 at tdata(2) and tdata(4);
!    * it is 0 for tval <= tdata(1) or tdata(5) <= tval;
!    * it is strictly increasing from tdata(1) to tdata(3),
!      and strictly decreasing from tdata(3) to tdata(5);
!    * the function and its first two derivatives are continuous
!      at each node tdata(i).
!
!  reference:
!
!    alan davies and philip samuels,
!    an introduction to computational geometry for curves and surfaces,
!    clarendon press, 1996.
!
!  modified:
!
!    06 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real tdata(5), the nodes associated with the basis function.
!    the entries of tdata are assumed to be distinct and increasing.
!
!    input, real tval, a point at which the b spline basis function is
!    to be evaluated.
!
!    output, real yval, the value of the function at tval.
!
  implicit none
!
  integer, parameter :: ndata = 5
!
  integer left
  integer right
  real tdata(ndata)
  real tval
  real u
  real yval
!
  if ( tval <= tdata(1) .or. tval >= tdata(ndata) ) then
    yval = 0.0e+00
    return
  end if
!
!  find the interval [ tdata(left), tdata(right) ] containing tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  u is the normalized coordinate of tval in this interval.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
!
!  now evaluate the function.
!
  if ( tval < tdata(2) ) then
    yval = u**3 / 6.0e+00
  else if ( tval < tdata(3) ) then
    yval = ( 1.0e+00 + 3.0e+00 * u + 3.0e+00 * u**2 - 3.0e+00 * u**3 ) / 6.0e+00
  else if ( tval < tdata(4) ) then
    yval = ( 4.0e+00 - 6.0e+00 * u**2 + 3.0e+00 * u**3 ) / 6.0e+00
  else if ( tval < tdata(5) ) then
    yval = ( 1.0e+00 - u )**3 / 6.0e+00
  end if

  return
end
subroutine basis_function_beta_val ( beta1, beta2, tdata, tval, yval )
!
!*******************************************************************************
!
!! basis_function_beta_val evaluates the beta spline basis function.
!
!
!  discussion:
!
!    with beta1 = 1 and beta2 = 0, the beta spline basis function 
!    equals the b spline basis function.
!
!    with beta1 large, and beta2 = 0, the beta spline basis function
!    skews to the right, that is, its maximum increases, and occurs
!    to the right of the center.
!
!    with beta1 = 1 and beta2 large, the beta spline becomes more like
!    a linear basis function; that is, its value in the outer two intervals
!    goes to zero, and its behavior in the inner two intervals approaches
!    a piecewise linear function that is 0 at the second node, 1 at the
!    third, and 0 at the fourth.
!
!  reference:
!
!    alan davies and philip samuels,
!    an introduction to computational geometry for curves and surfaces,
!    clarendon press, 1996, page 129.
!
!  modified:
!
!    09 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real beta1, the skew or bias parameter.
!    beta1 = 1 for no skew or bias.
!
!    input, real beta2, the tension parameter.
!    beta2 = 0 for no tension.
!
!    input, real tdata(5), the nodes associated with the basis function.
!    the entries of tdata are assumed to be distinct and increasing.
!
!    input, real tval, a point at which the b spline basis function is
!    to be evaluated.
!
!    output, real yval, the value of the function at tval.
!
  implicit none
!
  integer, parameter :: ndata = 5
!
  real a
  real b
  real beta1
  real beta2
  real c
  real d
  integer left
  integer right
  real tdata(ndata)
  real tval
  real u
  real yval
!
  if ( tval <= tdata(1) .or. tval >= tdata(ndata) ) then
    yval = 0.0e+00
    return
  end if
!
!  find the interval [ tdata(left), tdata(right) ] containing tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  u is the normalized coordinate of tval in this interval.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
!
!  now evaluate the function.
!
  if ( tval < tdata(2) ) then

    yval = 2.0e+00 * u**3 

  else if ( tval < tdata(3) ) then

    a = beta2 + 4.0e+00 * beta1 + 4.0e+00 * beta1**2 &
      + 6.0e+00 * ( 1.0e+00 - beta1**2 ) &
      - 3.0e+00 * ( 2.0e+00 + beta2 + 2.0e+00 * beta1 ) &
      + 2.0e+00 * ( 1.0e+00 + beta2 + beta1 + beta1**2 )

    b = - 6.0e+00 * ( 1.0e+00 - beta1**2 ) &
        + 6.0e+00 * ( 2.0e+00 + beta2 + 2.0e+00 * beta1 ) &
        - 6.0e+00 * ( 1.0e+00 + beta2 + beta1 + beta1**2 )

    c = - 3.0e+00 * ( 2.0e+00 + beta2 + 2.0e+00 * beta1 ) &
        + 6.0e+00 * ( 1.0e+00 + beta2 + beta1 + beta1**2 )

    d = - 2.0e+00 * ( 1.0e+00 + beta2 + beta1 + beta1**2 )

    yval = a + b * u + c * u**2 + d * u**3

  else if ( tval < tdata(4) ) then

    a = beta2 + 4.0e+00 * beta1 + 4.0e+00 * beta1**2

    b = - 6.0e+00 * ( beta1 - beta1**3 )

    c = - 3.0e+00 * ( beta2 + 2.0e+00 * beta1**2 + 2.0e+00 * beta1**3 )

    d = 2.0e+00 * ( beta2 + beta1 + beta1**2 + beta1**3 )

    yval = a + b * u + c * u**2 + d * u**3

  else if ( tval < tdata(5) ) then

    yval = 2.0e+00 * beta1**3 * ( 1.0e+00 - u )**3

  end if

  yval = yval / ( 2.0e+00 + beta2 + 4.0e+00 * beta1 + 4.0e+00 * beta1**2 &
    + 2.0e+00 * beta1**3 )

  return
end
subroutine basis_matrix_b_uni ( mbasis )
!
!*******************************************************************************
!
!! basis_matrix_b_uni sets up the uniform b spline basis matrix.
!
!
!  reference:
!
!    foley, van dam, feiner, hughes,
!    computer graphics: principles and practice,
!    page 493.
!
!  modified:
!
!    05 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    output, real mbasis(4,4), the basis matrix.
!
  implicit none
!
  real mbasis(4,4)
!
  mbasis(1,1) = - 1.0e+00 / 6.0e+00
  mbasis(1,2) =   3.0e+00 / 6.0e+00
  mbasis(1,3) = - 3.0e+00 / 6.0e+00
  mbasis(1,4) =   1.0e+00 / 6.0e+00

  mbasis(2,1) =   3.0e+00 / 6.0e+00
  mbasis(2,2) = - 6.0e+00 / 6.0e+00
  mbasis(2,3) =   3.0e+00 / 6.0e+00
  mbasis(2,4) =   0.0e+00

  mbasis(3,1) = - 3.0e+00 / 6.0e+00
  mbasis(3,2) =   0.0e+00
  mbasis(3,3) =   3.0e+00 / 6.0e+00
  mbasis(3,4) =   0.0e+00

  mbasis(4,1) =   1.0e+00 / 6.0e+00
  mbasis(4,2) =   4.0e+00 / 6.0e+00
  mbasis(4,3) =   1.0e+00 / 6.0e+00
  mbasis(4,4) =   0.0e+00

  return
end
subroutine basis_matrix_beta_uni ( beta1, beta2, mbasis )
!
!*******************************************************************************
!
!! basis_matrix_beta_uni sets up the uniform beta spline basis matrix.
!
!
!  discussion:
!
!    if beta1 = 1 and beta2 = 0, then the beta spline reduces to 
!    the b spline.
!
!  reference:
!
!    foley, van dam, feiner, hughes,
!    computer graphics: principles and practice,
!    page 505.
!
!  modified:
!
!    03 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real beta1, the skew or bias parameter.
!    beta1 = 1 for no skew or bias.
!
!    input, real beta2, the tension parameter.
!    beta2 = 0 for no tension.
!
!    output, real mbasis(4,4), the basis matrix.
!
  implicit none
!
  real beta1
  real beta2
  real delta
  integer i
  integer j
  real mbasis(4,4)
!
  mbasis(1,1) = - 2.0e+00 * beta1**3
  mbasis(1,2) =   2.0e+00 * ( beta2 + beta1**3 + beta1**2 + beta1 )
  mbasis(1,3) = - 2.0e+00 * ( beta2 + beta1**2 + beta1 + 1.0e+00 )
  mbasis(1,4) =   2.0e+00

  mbasis(2,1) =   6.0e+00 * beta1**3
  mbasis(2,2) = - 3.0e+00 * ( beta2 + 2.0e+00 * beta1**3 + 2.0e+00 * beta1**2 )
  mbasis(2,3) =   3.0e+00 * ( beta2 + 2.0e+00 * beta1**2 )
  mbasis(2,4) =   0.0e+00

  mbasis(3,1) = - 6.0e+00 * beta1**3
  mbasis(3,2) =   6.0e+00 * beta1 * ( beta1 - 1.0e+00 ) * ( beta1 + 1.0e+00 )
  mbasis(3,3) =   6.0e+00 * beta1
  mbasis(3,4) =   0.0e+00

  mbasis(4,1) =   2.0e+00 * beta1**3
  mbasis(4,2) =   4.0e+00 * beta1 * ( beta1 + 1.0e+00 ) + beta2
  mbasis(4,3) =   2.0e+00
  mbasis(4,4) =   0.0e+00

  delta = beta2 + 2.0e+00 * beta1**3 + 4.0e+00 * beta1**2 &
    + 4.0e+00 * beta1 + 2.0e+00

  mbasis(1:4,1:4) = mbasis(1:4,1:4) / delta

  return
end
subroutine basis_matrix_bezier ( mbasis )
!
!*******************************************************************************
!
!! basis_matrix_bezier_uni sets up the cubic bezier spline basis matrix.
!
!
!  discussion:
!
!    this basis matrix assumes that the data points are stored as
!    ( p1, p2, p3, p4 ).  p1 is the function value at t = 0, while
!    p2 is used to approximate the derivative at t = 0 by
!    dp/dt = 3 * ( p2 - p1 ).  similarly, p4 is the function value
!    at t = 1, and p3 is used to approximate the derivative at t = 1
!    by dp/dt = 3 * ( p4 - p3 ).
!
!  reference:
!
!    foley, van dam, feiner, hughes,
!    computer graphics: principles and practice,
!    page 489.
!
!  modified:
!
!    05 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    output, real mbasis(4,4), the basis matrix.
!
  implicit none
!
  real mbasis(4,4)
!
  mbasis(1,1) = - 1.0e+00
  mbasis(1,2) =   3.0e+00
  mbasis(1,3) = - 3.0e+00
  mbasis(1,4) =   1.0e+00

  mbasis(2,1) =   3.0e+00
  mbasis(2,2) = - 6.0e+00
  mbasis(2,3) =   3.0e+00
  mbasis(2,4) =   0.0e+00

  mbasis(3,1) = - 3.0e+00
  mbasis(3,2) =   3.0e+00
  mbasis(3,3) =   0.0e+00
  mbasis(3,4) =   0.0e+00

  mbasis(4,1) =   1.0e+00
  mbasis(4,2) =   0.0e+00
  mbasis(4,3) =   0.0e+00
  mbasis(4,4) =   0.0e+00

  return
end
subroutine basis_matrix_hermite ( mbasis )
!
!*******************************************************************************
!
!! basis_matrix_hermite sets up the hermite spline basis matrix.
!
!
!  discussion:
!
!    this basis matrix assumes that the data points are stored as
!    ( p1, p2, p1', p2' ), with p1 and p1' being the data value and 
!    the derivative dp/dt at t = 0, while p2 and p2' apply at t = 1.
!
!  reference:
!
!    foley, van dam, feiner, hughes,
!    computer graphics: principles and practice,
!    page 484.
!
!  modified:
!
!    06 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    output, real mbasis(4,4), the basis matrix.
!
  implicit none
!
  real mbasis(4,4)
!
  mbasis(1,1) =   2.0e+00
  mbasis(1,2) = - 2.0e+00
  mbasis(1,3) =   1.0e+00
  mbasis(1,4) =   1.0e+00

  mbasis(2,1) = - 3.0e+00
  mbasis(2,2) =   3.0e+00
  mbasis(2,3) = - 2.0e+00
  mbasis(2,4) = - 1.0e+00

  mbasis(3,1) =   0.0e+00
  mbasis(3,2) =   0.0e+00
  mbasis(3,3) =   1.0e+00
  mbasis(3,4) =   0.0e+00

  mbasis(4,1) =   1.0e+00
  mbasis(4,2) =   0.0e+00
  mbasis(4,3) =   0.0e+00
  mbasis(4,4) =   0.0e+00

  return
end
subroutine basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )
!
!*******************************************************************************
!
!! basis_matrix_overhauser_nonuni sets up the nonuniform overhauser spline basis matrix.
!
!
!  discussion:
!
!    this basis matrix assumes that the data points p1, p2, p3 and
!    p4 are not uniformly spaced in t, and that p2 corresponds to t = 0,
!    and p3 to t = 1.
!
!  modified:
!
!    05 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real alpha, beta.
!    alpha = | p2 - p1 | / ( | p3 - p2 | + | p2 - p1 | )
!    beta  = | p3 - p2 | / ( | p4 - p3 | + | p3 - p2 | ).
!
!    output, real mbasis(4,4), the basis matrix.
!
  implicit none
!
  real alpha
  real beta
  real mbasis(4,4)
!
  mbasis(1,1) = - ( 1.0e+00 - alpha )**2 / alpha
  mbasis(1,2) =   beta + ( 1.0e+00 - alpha ) / alpha
  mbasis(1,3) =   alpha - 1.0e+00 / ( 1.0e+00 - beta )
  mbasis(1,4) =   beta**2 / ( 1.0e+00 - beta )

  mbasis(2,1) =   2.0e+00 * ( 1.0e+00 - alpha )**2 / alpha
  mbasis(2,2) = ( - 2.0e+00 * ( 1.0e+00 - alpha ) - alpha * beta ) / alpha
  mbasis(2,3) = ( 2.0e+00 * ( 1.0e+00 - alpha ) &
    - beta * ( 1.0e+00 - 2.0e+00 * alpha ) ) / ( 1.0e+00 - beta )
  mbasis(2,4) = - beta**2 / ( 1.0e+00 - beta )

  mbasis(3,1) = - ( 1.0e+00 - alpha )**2 / alpha
  mbasis(3,2) =   ( 1.0e+00 - 2.0e+00 * alpha ) / alpha
  mbasis(3,3) =   alpha
  mbasis(3,4) =   0.0e+00

  mbasis(4,1) =   0.0e+00
  mbasis(4,2) =   1.0e+00
  mbasis(4,3) =   0.0e+00
  mbasis(4,4) =   0.0e+00

  return
end
subroutine basis_matrix_overhauser_nul ( alpha, mbasis )
!
!*******************************************************************************
!
!! basis_matrix_overhauser_nul sets up the left nonuniform overhauser spline basis matrix.
!
!
!  discussion:
!
!    this basis matrix assumes that the data points p1, p2, and
!    p3 are not uniformly spaced in t, and that p1 corresponds to t = 0,
!    and p2 to t = 1. (???)
!
!  modified:
!
!    27 august 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real alpha.
!    alpha = | p2 - p1 | / ( | p3 - p2 | + | p2 - p1 | )
!
!    output, real mbasis(3,3), the basis matrix.
!
  implicit none
!
  real alpha
  real mbasis(3,3)
!
  mbasis(1,1) =   1.0e+00 / alpha
  mbasis(1,2) = - 1.0e+00 / ( alpha * ( 1.0e+00 - alpha ) )
  mbasis(1,3) =   1.0e+00 / ( 1.0e+00 - alpha )

  mbasis(2,1) = - ( 1.0e+00 + alpha ) / alpha
  mbasis(2,2) =   1.0e+00 / ( alpha * ( 1.0e+00 - alpha ) )
  mbasis(2,3) = - alpha / ( 1.0e+00 - alpha )

  mbasis(3,1) =   1.0e+00
  mbasis(3,2) =   0.0e+00
  mbasis(3,3) =   0.0e+00

  return
end
subroutine basis_matrix_overhauser_nur ( beta, mbasis )
!
!*******************************************************************************
!
!! basis_matrix_overhauser_nur sets up the right nonuniform overhauser spline basis matrix.
!
!
!  discussion:
!
!    this basis matrix assumes that the data points pn-2, pn-1, and
!    pn are not uniformly spaced in t, and that pn-1 corresponds to t = 0,
!    and pn to t = 1. (???)
!
!  modified:
!
!    27 august 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real beta.
!    beta = | p(n) - p(n-1) | / ( | p(n) - p(n-1) | + | p(n-1) - p(n-2) | )
!
!    output, real mbasis(3,3), the basis matrix.
!
  implicit none
!
  real beta
  real mbasis(3,3)
!
  mbasis(1,1) =   1.0e+00 / beta
  mbasis(1,2) = - 1.0e+00 / ( beta * ( 1.0e+00 - beta ) )
  mbasis(1,3) =   1.0e+00 / ( 1.0e+00 - beta )

  mbasis(2,1) = - ( 1.0e+00 + beta ) / beta
  mbasis(2,2) =   1.0e+00 / ( beta * ( 1.0e+00 - beta ) )
  mbasis(2,3) = - beta / ( 1.0e+00 - beta )

  mbasis(3,1) =   1.0e+00
  mbasis(3,2) =   0.0e+00
  mbasis(3,3) =   0.0e+00

  return
end
subroutine basis_matrix_overhauser_uni ( mbasis )
!
!*******************************************************************************
!
!! basis_matrix_overhauser_uni sets up the uniform overhauser spline basis matrix.
!
!
!  discussion:
!
!    this basis matrix assumes that the data points p1, p2, p3 and
!    p4 are uniformly spaced in t, and that p2 corresponds to t = 0,
!    and p3 to t = 1.
!
!  reference:
!
!    foley, van dam, feiner, hughes,
!    computer graphics: principles and practice,
!    page 505.
!
!  modified:
!
!    05 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    output, real mbasis(4,4), the basis matrix.
!
  implicit none
!
  real mbasis(4,4)
!
  mbasis(1,1) = - 1.0e+00 / 2.0e+00
  mbasis(1,2) =   3.0e+00 / 2.0e+00
  mbasis(1,3) = - 3.0e+00 / 2.0e+00
  mbasis(1,4) =   1.0e+00 / 2.0e+00

  mbasis(2,1) =   2.0e+00 / 2.0e+00
  mbasis(2,2) = - 5.0e+00 / 2.0e+00
  mbasis(2,3) =   4.0e+00 / 2.0e+00
  mbasis(2,4) = - 1.0e+00 / 2.0e+00

  mbasis(3,1) = - 1.0e+00 / 2.0e+00
  mbasis(3,2) =   0.0e+00
  mbasis(3,3) =   1.0e+00 / 2.0e+00
  mbasis(3,4) =   0.0e+00

  mbasis(4,1) =   0.0e+00
  mbasis(4,2) =   2.0e+00 / 2.0e+00
  mbasis(4,3) =   0.0e+00
  mbasis(4,4) =   0.0e+00

  return
end
subroutine basis_matrix_overhauser_uni_l ( mbasis )
!
!*******************************************************************************
!
!! basis_matrix_overhauser_uni_l sets up the left uniform overhauser spline basis matrix.
!
!
!  discussion:
!
!    this basis matrix assumes that the data points p1, p2, and p3
!    are not uniformly spaced in t, and that p1 corresponds to t = 0,
!    and p2 to t = 1.
!
!  modified:
!
!    05 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    output, real mbasis(3,3), the basis matrix.
!
  implicit none
!
  real mbasis(3,3)
!
  mbasis(1,1) =   2.0e+00
  mbasis(1,2) = - 4.0e+00
  mbasis(1,3) =   2.0e+00

  mbasis(2,1) = - 3.0e+00
  mbasis(2,2) =   4.0e+00
  mbasis(2,3) = - 1.0e+00

  mbasis(3,1) =   1.0e+00
  mbasis(3,2) =   0.0e+00
  mbasis(3,3) =   0.0e+00

  return
end
subroutine basis_matrix_overhauser_uni_r ( mbasis )
!
!*******************************************************************************
!
!! basis_matrix_overhauser_uni_r sets up the right uniform overhauser spline basis matrix.
!
!
!  discussion:
!
!    this basis matrix assumes that the data points p(n-2), p(n-1), 
!    and p(n) are uniformly spaced in t, and that p(n-1) corresponds to 
!    t = 0, and p(n) to t = 1.
!
!  modified:
!
!    05 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    output, real mbasis(3,3), the basis matrix.
!
  implicit none
!
  real mbasis(3,3)
!
  mbasis(1,1) =   2.0e+00
  mbasis(1,2) = - 4.0e+00
  mbasis(1,3) =   2.0e+00

  mbasis(2,1) = - 3.0e+00
  mbasis(2,2) =   4.0e+00
  mbasis(2,3) = - 1.0e+00

  mbasis(3,1) =   1.0e+00
  mbasis(3,2) =   0.0e+00
  mbasis(3,3) =   0.0e+00

  return
end
subroutine basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! basis_matrix_tmp computes q = t * mbasis * p
!
!
!  discussion:
!
!    ydata is a vector of data values, most frequently the values of some
!    function sampled at uniformly spaced points.  mbasis is the basis
!    matrix for a particular kind of spline.  t is a vector of the
!    powers of the normalized difference between tval and the left
!    endpoint of the interval.
!
!  modified:
!
!    06 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer left, indicats that tval is in the interval
!    [ tdata(left), tdata(left+1) ], or that this is the "nearest"
!    interval to tval.
!    for tval < tdata(1), use left = 1.
!    for tval > tdata(ndata), use left = ndata - 1.
!
!    input, integer n, the order of the basis matrix.
!
!    input, real mbasis(n,n), the basis matrix.
!
!    input, integer ndata, the dimension of the vectors tdata and ydata.
!
!    input, real tdata(ndata), the abscissa values.  this routine
!    assumes that the tdata values are uniformly spaced, with an
!    increment of 1.0.
!
!    input, real ydata(ndata), the data values to be interpolated or 
!    approximated.
!
!    input, real tval, the value of t at which the spline is to be
!    evaluated.
!
!    output, real yval, the value of the spline at tval.
!
  implicit none
!
  integer, parameter :: maxn = 4
  integer n
  integer ndata
!
  real arg
  integer first
  integer i
  integer j
  integer left
  real mbasis(n,n)
  real tdata(ndata)
  real temp
  real tval
  real tvec(maxn)
  real ydata(ndata)
  real yval
!
  if ( left == 1 ) then
    arg = 0.5 * ( tval - tdata(left) )
    first = left
  else if ( left < ndata - 1 ) then
    arg = tval - tdata(left)
    first = left - 1
  else if ( left == ndata - 1 ) then
    arg = 0.5e+00 * ( 1.0e+00 + tval - tdata(left) )
    first = left - 1
  end if

  do i = 1, n - 1
    tvec(i) = arg**( n - i )
  end do
  tvec(n) = 1.0e+00

  yval = 0.0e+00
  do j = 1, n
    yval = yval + dot_product ( tvec(1:n), mbasis(1:n,j) ) &
      * ydata(first - 1 + j)
  end do

  return
end
subroutine bc_val ( n, t, xcon, ycon, xval, yval )
!
!*******************************************************************************
!
!! bc_val evaluates a parameterized bezier curve.
!
!
!  discussion:
!
!    bc_val(t) is the value of a vector function of the form
!
!      bc_val(t) = ( x(t), y(t) )
!
!    where
!
!      x(t) = sum (i = 0 to n) xcon(i) * bern(i,n)(t)
!      y(t) = sum (i = 0 to n) ycon(i) * bern(i,n)(t)
!
!    where bern(i,n)(t) is the i-th bernstein polynomial of order n
!    defined on the interval [0,1], and where xcon(*) and ycon(*)
!    are the coordinates of n+1 "control points".
!
!  modified:
!
!    02 march 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the order of the bezier curve, which
!    must be at least 0.
!
!    input, real t, the point at which the bezier curve should
!    be evaluated.  the best results are obtained within the interval
!    [0,1] but t may be anywhere.
!
!    input, real xcon(0:n), ycon(0:n), the x and y coordinates
!    of the control points.  the bezier curve will pass through
!    the points ( xcon(0), ycon(0) ) and ( xcon(n), ycon(n) ), but
!    generally not through the other control points.
!
!    output, real xval, yval, the x and y coordinates of the point
!    on the bezier curve corresponding to the given t value.
!
  implicit none
!
  integer n
!
  real bval(0:n)
  integer i
  real t
  real xcon(0:n)
  real xval
  real ycon(0:n)
  real yval
!
  call bp01 ( n, bval, t )

  xval = dot_product ( xcon(0:n), bval(0:n) )
  yval = dot_product ( ycon(0:n), bval(0:n) )

  return
end
function bez_val ( n, x, a, b, y )
!
!*******************************************************************************
!
!! bez_val evaluates a bezier function at a point.
!
!
!  discussion:
!
!    the bezier function has the form:
!
!      bez(x) = sum (i = 0 to n) y(i) * bern(n,i)( (x-a)/(b-a) )
!
!    where bern(n,i)(x) is the i-th bernstein polynomial of order n
!    defined on the interval [0,1], and y(*) is a set of
!    coefficients.
!
!    if we define the n+1 equally spaced
!
!      x(i) = ( (n-i)*a + i*b)/ n,
!
!    for i = 0 to n, then the pairs ( x(i), y(i) ) can be regarded as
!    "control points".
!
!  modified:
!
!    02 march 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the order of the bezier function, which
!    must be at least 0.
!
!    input, real x, the point at which the bezier function should
!    be evaluated.  the best results are obtained within the interval
!    [a,b] but x may be anywhere.
!
!    input, real a, b, the interval over which the bezier function
!    has been defined.  this is the interval in which the control
!    points have been set up.  note bez(a) = y(0) and bez(b) = y(n),
!    although bez will not, in general pass through the other
!    control points.  a and b must not be equal.
!
!    input, real y(0:n), a set of data defining the y coordinates
!    of the control points.
!
!    output, real bez_val, the value of the bezier function at x.
!
  implicit none
!
  integer n
!
  real a
  real b
  real bez_val
  real bval(0:n)
  integer i
  integer ncopy
  real x
  real x01
  real y(0:n)
!
  if ( b - a == 0.0e+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'bez_val - fatal error!'
    write ( *, '(a,g14.6)' ) '  null interval, a = b = ', a
    stop
  end if

  x01 = ( x - a ) / ( b - a )

  ncopy = n
  call bp01 ( ncopy, bval, x01 )

  bez_val = dot_product ( y(0:n), bval(0:n) )

  return
end
subroutine bp01 ( n, bern, x )
!
!*******************************************************************************
!
!! bp01 evaluates the bernstein basis polynomials for [01,1] at a point x.
!
!
!  formula:
!
!    bern(n,i,x) = [n!/(i!*(n-i)!)] * (1-x)**(n-i) * x**i
!
!  first values:
!
!    b(0,0,x) = 1
!
!    b(1,0,x) =      1-x
!    b(1,1,x) =                x
!
!    b(2,0,x) =     (1-x)**2
!    b(2,1,x) = 2 * (1-x)    * x
!    b(2,2,x) =                x**2
!
!    b(3,0,x) =     (1-x)**3
!    b(3,1,x) = 3 * (1-x)**2 * x
!    b(3,2,x) = 3 * (1-x)    * x**2
!    b(3,3,x) =                x**3
!
!    b(4,0,x) =     (1-x)**4
!    b(4,1,x) = 4 * (1-x)**3 * x
!    b(4,2,x) = 6 * (1-x)**2 * x**2
!    b(4,3,x) = 4 * (1-x)    * x**3
!    b(4,4,x) =                x**4
!
!  special values:
!
!    b(n,i,1/2) = c(n,k) / 2**n
!
!  modified:
!
!    12 january 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the degree of the bernstein basis polynomials.
!    for any n greater than or equal to 0, there is a set of n+1 bernstein
!    basis polynomials, each of degree n, which form a basis for
!    all polynomials on [0,1].
!
!    output, real bern(0:n), the values of the n+1 bernstein basis
!    polynomials at x.
!
!    input, real x, the point at which the polynomials are to be
!    evaluated.
!
  implicit none
!
  integer n
!
  real bern(0:n)
  integer i
  integer j
  real x
!
  if ( n == 0 ) then

    bern(0) = 1.0e+00

  else if ( n > 0 ) then

    bern(0) = 1.0e+00 - x
    bern(1) = x

    do i = 2, n
      bern(i) = x * bern(i-1)
      do j = i-1, 1, -1
        bern(j) = x * bern(j-1) + ( 1.0e+00 - x ) * bern(j)
      end do
      bern(0) = ( 1.0e+00 - x ) * bern(0)
    end do

  end if

  return
end
subroutine bp_approx ( n, a, b, ydata, xval, yval )
!
!*******************************************************************************
!
!! bp_approx evaluates the bernstein polynomial for f(x) on [a,b].
!
!
!  formula:
!
!    bern(f)(x) = sum ( 0 <= i <= n ) f(x(i)) * b_base(i,x)
!
!    where
!
!      x(i) = ( ( n - i ) * a + i * b ) / n
!      b_base(i,x) is the value of the i-th bernstein basis polynomial at x.
!
!  discussion:
!
!    the bernstein polynomial bern(f) for f(x) is an approximant, not an
!    interpolant; in other words, its value is not guaranteed to equal
!    that of f at any particular point.  however, for a fixed interval
!    [a,b], if we let n increase, the bernstein polynomial converges
!    uniformly to f everywhere in [a,b], provided only that f is continuous.
!    even if f is not continuous, but is bounded, the polynomial converges
!    pointwise to f(x) at all points of continuity.  on the other hand,
!    the convergence is quite slow compared to other interpolation
!    and approximation schemes.
!
!  modified:
!
!    10 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the degree of the bernstein polynomial to be used.
!
!    input, real a, b, the endpoints of the interval on which the
!    approximant is based.  a and b should not be equal.
!
!    input, real ydata(0:n), the data values at n+1 equally spaced points 
!    in [a,b].  if n = 0, then the evaluation point should be 0.5 * ( a + b).
!    otherwise, evaluation point i should be ( (n-i)*a + i*b ) / n ).
!
!    input, real xval, the point at which the bernstein polynomial 
!    approximant is to be evaluated.  xval does not have to lie in the 
!    interval [a,b].
!
!    output, real yval, the value of the bernstein polynomial approximant
!    for f, based in [a,b], evaluated at xval.
!
  implicit none
!
  integer n
!
  real a
  real b
  real bvec(0:n)
  integer i
  real xval
  real ydata(0:n)
  real yval
!
!  evaluate the bernstein basis polynomials at xval.
!
  call bpab ( n, bvec, xval, a, b )
!
!  now compute the sum of ydata(i) * bvec(i).
!
  yval = dot_product ( ydata(0:n), bvec(0:n) )

  return
end
subroutine bpab ( n, bern, x, a, b )
!
!*******************************************************************************
!
!! bpab evaluates the bernstein basis polynomials for [a,b] at a point x.
!
!
!  formula:
!
!    bern(n,i,x) = [n!/(i!*(n-i)!)] * (b-x)**(n-i) * (x-a)**i / (b-a)**n
!
!  first values:
!
!    b(0,0,x) =   1
!
!    b(1,0,x) = (      b-x                ) / (b-a)
!    b(1,1,x) = (                 x-a     ) / (b-a)
!
!    b(2,0,x) = (     (b-x)**2            ) / (b-a)**2
!    b(2,1,x) = ( 2 * (b-x)    * (x-a)    ) / (b-a)**2
!    b(2,2,x) = (                (x-a)**2 ) / (b-a)**2
!
!    b(3,0,x) = (     (b-x)**3            ) / (b-a)**3
!    b(3,1,x) = ( 3 * (b-x)**2 * (x-a)    ) / (b-a)**3
!    b(3,2,x) = ( 3 * (b-x)    * (x-a)**2 ) / (b-a)**3
!    b(3,3,x) = (                (x-a)**3 ) / (b-a)**3
!
!    b(4,0,x) = (     (b-x)**4            ) / (b-a)**4
!    b(4,1,x) = ( 4 * (b-x)**3 * (x-a)    ) / (b-a)**4
!    b(4,2,x) = ( 6 * (b-x)**2 * (x-a)**2 ) / (b-a)**4
!    b(4,3,x) = ( 4 * (b-x)    * (x-a)**3 ) / (b-a)**4
!    b(4,4,x) = (                (x-a)**4 ) / (b-a)**4
!
!  modified:
!
!    12 january 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the degree of the bernstein basis polynomials.
!    for any n greater than or equal to 0, there is a set of n+1 
!    bernstein basis polynomials, each of degree n, which form a basis 
!    for polynomials on [a,b].
!
!    output, real bern(0:n), the values of the n+1 bernstein basis
!    polynomials at x.
!
!    input, real x, the point at which the polynomials are to be
!    evaluated.  x need not lie in the interval [a,b].
!
!    input, real a, b, the endpoints of the interval on which the
!    polynomials are to be based.  a and b should not be equal.
!
  implicit none
!
  integer n
!
  real a
  real b
  real bern(0:n)
  integer i
  integer j
  real x
!
  if ( b == a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'bpab - fatal error!'
    write ( *, '(a,g14.6)' ) '  a = b = ', a
    stop
  end if

  if ( n == 0 ) then

    bern(0) = 1.0e+00

  else if ( n > 0 ) then

    bern(0) = ( b - x ) / ( b - a )
    bern(1) = ( x - a ) / ( b - a )

    do i = 2, n
      bern(i) = ( x - a ) * bern(i-1) / ( b - a )
      do j = i-1, 1, -1
        bern(j) = ( ( b - x ) * bern(j) + ( x - a ) * bern(j-1) ) / ( b - a )
      end do
      bern(0) = ( b - x ) * bern(0) / ( b - a )
    end do

  end if

  return
end
subroutine data_to_dif ( diftab, ntab, xtab, ytab )
!
!*******************************************************************************
!
!! data_to_dif sets up a divided difference table from raw data.
!
!
!  discussion:
!
!    space can be saved by using a single array for both the diftab and
!    ytab dummy parameters.  in that case, the divided difference table will
!    overwrite the y data without interfering with the computation.
!
!  modified:
!
!    11 april 2000
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    output, real diftab(ntab), the divided difference coefficients 
!    corresponding to the input (xtab,ytab).
!
!    input, integer ntab, the number of pairs of points
!    (xtab(i),ytab(i)) which are to be used as data.  the
!    number of entries to be used in diftab, xtab and ytab.
!
!    input, real xtab(ntab), the x values at which data was taken.
!    these values must be distinct.
!
!    input, real ytab(ntab), the corresponding y values.
!
  implicit none
!
  integer ntab
!
  real diftab(ntab)
  integer i
  integer j
  logical rvec_distinct
  real xtab(ntab)
  real ytab(ntab)
!
  if ( .not. rvec_distinct ( ntab, xtab ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'data_to_dif - fatal error!'
    write ( *, '(a)' ) '  two entries of xtab are equal!'
    return
  end if
!
!  copy the data values into diftab.
!
  diftab(1:ntab) = ytab(1:ntab)
!
!  compute the divided differences.
!
  do i = 2, ntab
    do j = ntab, i, -1

      diftab(j) = ( diftab(j) - diftab(j-1) ) / ( xtab(j) - xtab(j+1-i) )

    end do
  end do
 
  return
end
subroutine dif_val ( diftab, ntab, xtab, xval, yval )
!
!*******************************************************************************
!
!! dif_val evaluates a divided difference polynomial at a point.
!
!
!  modified:
!
!    11 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real diftab(ntab), the divided difference polynomial coefficients.
!
!    input, integer ntab, the number of divided difference
!    coefficients in diftab, and the number of points xtab.
!
!    input, real xtab(ntab), the x values upon which the
!    divided difference polynomial is based.
!
!    input, real xval, a value of x at which the polynomial
!    is to be evaluated.
!
!    output, real yval, the value of the polynomial at xval.
!
  implicit none
!
  integer ntab
!
  real diftab(ntab)
  integer i
  real xtab(ntab)
  real xval
  real yval
!
  yval = diftab(ntab)
  do i = 1, ntab-1
    yval = diftab(ntab-i) + ( xval - xtab(ntab-i) ) * yval
  end do
 
  return
end
subroutine least_set ( ntab, xtab, ytab, ndeg, ptab, array, eps, ierror )
!
!*******************************************************************************
!
!! least_set constructs the least squares polynomial approximation to data.
!
!
!  discussion:
!
!    the routine least_eval must be used to evaluate the approximation at a
!    point.
!
!  modified:
!
!    20 november 2000
!
!  parameters:
!
!    input, integer ntab, the number of data points.
!
!    input, real xtab(ntab), the x data.  the values in xtab
!    should be distinct, and in increasing order.
!
!    input, real ytab(ntab), the y data values corresponding
!    to the x data in xtab.
!
!    input, integer ndeg, the degree of the polynomial which the
!    program is to use.  ndeg must be at least 1, and less than or 
!    equal to ntab-1.
!
!    output, real ptab(ntab).  ptab(i) is the value of the
!    least squares polynomial at the point xtab(i).
!
!    output, real array(2*ntab+3*ndeg), an array containing data about 
!    the polynomial.
!
!    output, real eps, the root-mean-square discrepancy of the
!    polynomial fit.
!
!    output, integer ierror, error flag.
!    zero, no error occurred;
!    nonzero, an error occurred, and the polynomial could not be computed.
!
  implicit none
!
  integer ndeg
  integer ntab
!
  real array(2*ntab+3*ndeg)
  real eps
  real error
  integer i
  integer i0l1
  integer i1l1
  integer ierror
  integer it
  integer k
  integer mdeg
  real ptab(ntab)
  real rn0
  real rn1
  real s
  real sum2
  real xtab(ntab)
  real y_sum
  real ytab(ntab)
!
  ierror = 0
!
!  check ndeg.
!
  if ( ndeg < 1 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'least_set - fatal error!'
    write ( *, '(a)' ) '  ndeg < 1.'
    return
  else if ( ndeg >= ntab ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'least_set - fatal error!'
    write ( *, '(a)' ) '  ndeg >= ntab.'
    return
  end if
!
!  check that the abscissas are strictly increasing.
!
  do i = 1, ntab-1
    if ( xtab(i) >= xtab(i+1) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'least_set - fatal error!'
      write ( *, '(a)' ) '  xtab must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  xtab(', i, ') = ', xtab(i)
      write ( *, '(a,i6,a,g14.6)' ) '  xtab(', i+1, ') = ', xtab(i+1)
      return
    end if
  end do

  i0l1 = 3 * ndeg
  i1l1 = 3 * ndeg + ntab

  y_sum = sum ( ytab )
  rn0 = ntab
  array(2*ndeg) = y_sum / real ( ntab )

  ptab(1:ntab) = y_sum / real ( ntab )

  error = 0.0e+00
  do i = 1, ntab
    error = error + ( y_sum / real ( ntab ) - ytab(i) )**2
  end do

  if ( ndeg == 0 ) then
    eps = sqrt ( error / real ( ntab ) )
    return
  end if

  array(1) = sum ( xtab ) / real ( ntab )

  s = 0.0e+00
  sum2 = 0.0e+00
  do i = 1, ntab
    array(i1l1+i) = xtab(i) - array(1)
    s = s + array(i1l1+i)**2
    sum2 = sum2 + array(i1l1+i) * ( ytab(i) - ptab(i) )
  end do

  rn1 = s
  array(2*ndeg+1) = sum2 / s

  do i = 1, ntab
    ptab(i) = ptab(i) + sum2 * array(i1l1+i) / s
  end do

  error = sum ( ( ptab(1:ntab) - ytab(1:ntab) )**2 )

  if ( ndeg == 1 ) then
    eps = sqrt ( error / real ( ntab ) )
    return
  end if

  do i = 1, ntab
    array(3*ndeg+i) = 1.0e+00
  end do

  mdeg = 2
  k = 2

  do

    array(ndeg-1+k) = rn1 / rn0

    sum2 = 0.0e+00
    do i = 1, ntab
      sum2 = sum2 + xtab(i) * array(i1l1+i)**2
    end do

    array(k) = sum2 / rn1

    s = 0.0e+00
    sum2 = 0.0e+00
    do i = 1, ntab
      array(i0l1+i) = ( xtab(i) - array(k) ) * array(i1l1+i) &
        - array(ndeg-1+k) * array(i0l1+i)
      s = s + array(i0l1+i)**2
      sum2 = sum2 + array(i0l1+i) * ( ytab(i) - ptab(i) )
    end do

    rn0 = rn1
    rn1 = s
    it = i0l1
    i0l1 = i1l1
    i1l1 = it
    array(2*ndeg+k) = sum2 / rn1

    do i = 1, ntab
      ptab(i) = ptab(i) + sum2 * array(i1l1+i) / rn1
    end do

    error = sum ( ( ptab(1:ntab) - ytab(1:ntab) )**2 )

    if ( mdeg >= ndeg ) then
      exit
    end if

    mdeg = mdeg + 1
    k = k + 1

  end do

  eps = sqrt ( error / real ( ntab ) )

  return
end
subroutine least_val ( x, ndeg, array, value )
!
!*******************************************************************************
!
!! least_val evaluates a least squares polynomial defined by least_set.
!
!
!  modified:
!
!    01 march 1999
!
!  parameters:
!
!    input, real x, the point at which the polynomial is to be evaluated.
!
!    input, integer ndeg, the degree of the polynomial fit used.
!    this is the value of ndeg as returned from least_set.
!
!    input, real array(*), an array of a certain dimension.
!    see least_set for details on the size of array.
!    array contains information about the polynomial, as set up by least_set.
!
!    output, real value, the value of the polynomial at x.
!
  implicit none
!
  real array(*)
  real dk
  real dkp1
  real dkp2
  integer k
  integer l
  integer ndeg
  real value
  real x
!
  if ( ndeg <= 0 ) then

    value = array(2*ndeg)

  else if ( ndeg == 1 ) then

    value = array(2*ndeg) + array(2*ndeg+1) * ( x - array(1) )

  else

    dkp2 = array(3*ndeg)
    dkp1 = array(3*ndeg-1) + ( x - array(ndeg) ) * dkp2

    do l = 1, ndeg-2

      k = ndeg - 1 - l

      dk = array(2*ndeg+k) + ( x - array(k+1) ) * dkp1 - array(ndeg+1+k) * dkp2

      dkp2 = dkp1

      dkp1 = dk

    end do

    value = array(2*ndeg) + ( x - array(1) ) * dkp1 - array(ndeg+1) * dkp2

  end if

  return
end
subroutine parabola_val2 ( ndim, ndata, tdata, ydata, left, tval, yval )
!
!*******************************************************************************
!
!! parabola_val2 evaluates a parabolic interpolant through tabular data.
!
!
!  discussion:
!
!    this routine is a utility routine used by overhauser_spline_val.
!    it constructs the parabolic interpolant through the data in
!    3 consecutive entries of a table and evaluates this interpolant
!    at a given abscissa value.
!
!  modified:
!
!    02 march 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndim, the dimension of a single data point.
!    ndim must be at least 1.
!
!    input, integer ndata, the number of data points.
!    ndata must be at least 3.
!
!    input, real tdata(ndata), the abscissas of the data points.  the
!    values in tdata must be in strictly ascending order.
!
!    input, real ydata(ndim,ndata), the data points corresponding to
!    the abscissas.
!
!    input, integer left, the location of the first of the three
!    consecutive data points through which the parabolic interpolant
!    must pass.  1 <= left <= ndata - 2.
!
!    input, real tval, the value of t at which the parabolic interpolant
!    is to be evaluated.  normally, tdata(1) <= tval <= t(ndata), and
!    the data will be interpolated.  for tval outside this range,
!    extrapolation will be used.
!
!    output, real yval(ndim), the value of the parabolic interpolant at tval.
!
  implicit none
!
  integer ndata
  integer ndim
!
  real dif1
  real dif2
  integer i
  integer left
  real t1
  real t2
  real t3
  real tval
  real tdata(ndata)
  real ydata(ndim,ndata)
  real y1
  real y2
  real y3
  real yval(ndim)
!
!  check.
!
  if ( left < 1 .or. left > ndata-2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'parabola_val2 - fatal error!'
    write ( *, '(a)' ) '  left < 1 or left > ndata-2.'
    stop
  end if

  if ( ndim < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'parabola_val2 - fatal error!'
    write ( *, '(a)' ) '  ndim < 1.'
    stop
  end if
!
!  copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t1 >= t2 .or. t2 >= t3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'parabola_val2 - fatal error!'
    write ( *, '(a)' ) '  t1 >= t2 or t2 >= t3.'
    stop
  end if
!
!  construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  do i = 1, ndim

    y1 = ydata(i,left)
    y2 = ydata(i,left+1)
    y3 = ydata(i,left+2)

    dif1 = ( y2 - y1 ) / ( t2 - t1 )
    dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
         - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

    yval(i) = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )

  end do

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! r_random returns a random real in a given range.
!
!
!  modified:
!
!    06 april 2001
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real rlo, rhi, the minimum and maximum values.
!
!    output, real r, the randomly chosen value.
!
  implicit none
!
  real r
  real rhi
  real rlo
  real t
!
!  pick t, a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  set r in ( rlo, rhi ).
!
  r = ( 1.0e+00 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! r_swap swaps two real values.
!
!
!  modified:
!
!    01 may 2000
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input/output, real x, y.  on output, the values of x and
!    y have been interchanged.
!
  implicit none
!
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine rvec_bracket ( n, x, xval, left, right )
!
!*******************************************************************************
!
!! rvec_bracket searches a sorted array for successive brackets of a value.
!
!
!  discussion:
!
!    if the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  modified:
!
!    06 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, length of input array.
!
!    input, real x(n), an array sorted into ascending order.
!
!    input, real xval, a value to be bracketed.
!
!    output, integer left, right, the results of the search.
!    either:
!      xval < x(1), when left = 1, right = 2;
!      xval > x(n), when left = n-1, right = n;
!    or
!      x(left) <= xval <= x(right).
!
  implicit none
!
  integer n
!
  integer i
  integer left
  integer right
  real x(n)
  real xval
!
  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
subroutine rvec_bracket3 ( n, t, tval, left )
!
!*******************************************************************************
!
!! rvec_bracket3 finds the interval containing or nearest a given value.
!
!
!  discussion:
!
!    the routine always returns the index left of the sorted array
!    t with the property that either
!    *  t is contained in the interval [ t(left), t(left+1) ], or
!    *  t < t(left) = t(1), or
!    *  t > t(left+1) = t(n).
!
!    the routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  modified:
!
!    05 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, length of the input array.
!
!    input, real t(n), an array sorted into ascending order.
!
!    input, real tval, a value to be bracketed by entries of t.
!
!    input/output, integer left.
!
!    on input, if 1 <= left <= n-1, left is taken as a suggestion for the
!    interval [ t(left), t(left+1) ] in which tval lies.  this interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  after that, a binary search is used.
!
!    on output, left is set so that the interval [ t(left), t(left+1) ]
!    is the closest to tval; it either contains tval, or else tval
!    lies outside the interval [ t(1), t(n) ].
!
  implicit none
!
  integer n
!
  integer high
  integer left
  integer low
  integer mid
  real t(n)
  real tval
!
!  check the input data.
!
  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'rvec_bracket3 - fatal error!'
    write ( *, '(a)' ) '  n must be at least 2.'
    stop
  end if
!
!  if left is not between 1 and n-1, set it to the middle value.
!
  if ( left < 1 .or. left > n - 1 ) then
    left = ( n + 1 ) / 2
  end if
!
!  case 1: tval < t(left):
!  search for tval in [t(i), t(i+1)] for intervals i = 1 to left-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( tval >= t(left-1) ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...binary search for tval in [t(i), t(i+1)] for intervals i = 2 to left-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( tval >= t(mid) ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  case2: t(left+1) < tval:
!  search for tval in {t(i),t(i+1)] for intervals i = left+1 to n-1.
!
  else if ( tval > t(left+1) ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( tval >= t(n-1) ) then
      left = n - 1
      return
    end if
!
!  ...binary search for tval in [t(i), t(i+1)] for intervals i = left+2 to n-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( tval >= t(mid) ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  case3: t(left) <= tval <= t(left+1):
!  t is in [t(left), t(left+1)], as the user said it might be.
!
  else

  end if

  return
end
function rvec_distinct ( n, x )
!
!*******************************************************************************
!
!! rvec_distinct is true if the entries in a real vector are distinct.
!
!
!  modified:
!
!    16 september 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the number of entries in the vector.
!
!    input, real x(n), the vector to be checked.
!
!    output, logical rvec_distinct is .true. if all n elements of x 
!    are distinct.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  logical rvec_distinct
  real x(n)
!
  rvec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1 
      if ( x(i) == x(j) ) then
        return
      end if
    end do
  end do

  rvec_distinct = .true.

  return
end
subroutine rvec_even ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! rvec_even returns n real values, evenly spaced between alo and ahi.
!
!
!  modified:
!
!    31 october 2000
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real alo, ahi, the low and high values.
!
!    input, integer n, the number of values.
!
!    output, real a(n), n evenly spaced values.
!    normally, a(1) = alo and a(n) = ahi.
!    however, if n = 1, then a(1) = 0.5*(alo+ahi).
!
  implicit none
!
  integer n
!
  real a(n)
  real ahi
  real alo
  integer i
!
  if ( n == 1 ) then

    a(1) = 0.5e+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i ) * alo + real ( i - 1 ) * ahi ) / real ( n - 1 )
    end do

  end if

  return
end
subroutine rvec_order_type ( n, a, order )
!
!*******************************************************************************
!
!! rvec_order_type determines if a real array is (non)strictly ascending/descending.
!
!
!  modified:
!
!    20 july 2000
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the number of entries of the array.
!
!    input, real a(n), the array to be checked.
!
!    output, integer order, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  integer order
!
!  search for the first value not equal to a(1).
!
  i = 1

  do

    i = i + 1

    if ( i > n ) then
      order = 0
      return
    end if

    if ( a(i) > a(1) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  now we have a "direction".  examine subsequent entries.
!
  do

    i = i + 1
    if ( i > n ) then
      exit
    end if

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i) > a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i) > a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do
 
  return
end
subroutine rvec_print ( n, a, title )
!
!*******************************************************************************
!
!! rvec_print prints a real vector.
!
!
!  modified:
!
!    16 december 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the number of components of the vector.
!
!    input, real a(n), the vector to be printed.
!
!    input, character ( len = * ) title, a title to be printed first.
!    title may be blank.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine rvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! rvec_random returns a random real vector in a given range.
!
!
!  modified:
!
!    04 february 2001
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real alo, ahi, the range allowed for the entries.
!
!    input, integer n, the number of entries in the vector.
!
!    output, real a(n), the vector of randomly chosen values.
!
  implicit none
!
  integer n
!
  real a(n)
  real ahi
  real alo
  integer i
!
  do i = 1, n
    call r_random ( alo, ahi, a(i) )
  end do

  return
end
subroutine rvec_sort_bubble_a ( n, a )
!
!*******************************************************************************
!
!! rvec_sort_bubble_a ascending sorts a real array using bubble sort.
!
!
!  discussion:
!
!    bubble sort is simple to program, but inefficient.  it should not
!    be used for large arrays.
!
!  modified:
!
!    01 february 2001
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the number of entries in the array.
!
!    input/output, real a(n).
!    on input, an unsorted array.
!    on output, the array has been sorted.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  integer j
!
  do i = 1, n-1
    do j = i+1, n
      if ( a(i) > a(j) ) then
        call r_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine s3_fs ( a1, a2, a3, n, b, x )
!
!*******************************************************************************
!
!! s3_fs factors and solves a tridiagonal linear system.
!
!
!  note:
!
!    this algorithm requires that each diagonal entry be nonzero.
!
!  modified:
!
!    05 december 1998
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input/output, real a1(2:n), a2(1:n), a3(1:n-1).
!    on input, the nonzero diagonals of the linear system.
!    on output, the data in these vectors has been overwritten
!    by factorization information.
!
!    input, integer n, the order of the linear system.
!
!    input/output, real b(n).
!    on input, b contains the right hand side of the linear system.
!    on output, b has been overwritten by factorization information.
!
!    output, real x(n), the solution of the linear system.
!
  implicit none
!
  integer n
!
  real a1(2:n)
  real a2(1:n)
  real a3(1:n-1)
  real b(n)
  integer i
  real x(n)
  real xmult
!
!  the diagonal entries can't be zero.
!
  do i = 1, n
    if ( a2(i) == 0.0e+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 's3_fs - fatal error!'
      write ( *, '(a,i6,a)' ) '  a2(', i, ') = 0.'
      return
    end if
  end do

  do i = 2, n-1

    xmult = a1(i) / a2(i-1)
    a2(i) = a2(i) - xmult * a3(i-1)

    b(i) = b(i) - xmult * b(i-1)

  end do

  xmult = a1(n) / a2(n-1)
  a2(n) = a2(n) - xmult * a3(n-1)

  x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
  do i = n-1, 1, -1
    x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
  end do

  return
end
subroutine sgtsl ( n, c, d, e, b, info )
!
!*******************************************************************************
!
!! sgtsl solves a general tridiagonal linear system.
!
!
!  reference:
!
!    dongarra, moler, bunch and stewart,
!    linpack user's guide,
!    siam, (society for industrial and applied mathematics),
!    3600 university city science center,
!    philadelphia, pa, 19104-2688.
!    isbn 0-89871-172-x
!
!  modified:
!
!    31 october 2001
!
!  parameters:
!
!    input, integer n, the order of the tridiagonal matrix.
!
!    input/output, real c(n), contains the subdiagonal of the tridiagonal
!    matrix in entries c(2:n).  on output, c is destroyed.
!
!    input/output, real d(n).  on input, the diagonal of the matrix.
!    on output, d is destroyed.
!
!    input/output, real e(n), contains the superdiagonal of the tridiagonal
!    matrix in entries e(1:n-1).  on output e is destroyed.
!
!    input/output, real b(n).  on input, the right hand side.  on output,
!    the solution.
!
!    output, integer info, error flag.
!    0, normal value.
!    k, the k-th element of the diagonal becomes exactly zero.  the
!       subroutine returns if this error condition is detected.
!
  implicit none
!
  integer n
!
  real b(n)
  real c(n)
  real d(n)
  real e(n)
  integer info
  integer k
  real t
!
  info = 0
  c(1) = d(1)

  if ( n >= 2 ) then

    d(1) = e(1)
    e(1) = 0.0e+00
    e(n) = 0.0e+00

    do k = 1, n - 1
!
!  find the larger of the two rows, and interchange if necessary.
!
      if ( abs ( c(k+1) ) >= abs ( c(k) ) ) then
        call r_swap ( c(k), c(k+1) )
        call r_swap ( d(k), d(k+1) )
        call r_swap ( e(k), e(k+1) )
        call r_swap ( b(k), b(k+1) )
      end if
!
!  fail if no nonzero pivot could be found.
!
      if ( c(k) == 0.0e+00 ) then
        info = k
        return
      end if
!
!  zero elements.
!
      t = -c(k+1) / c(k)
      c(k+1) = d(k+1) + t * d(k)
      d(k+1) = e(k+1) + t * e(k)
      e(k+1) = 0.0e+00
      b(k+1) = b(k+1) + t * b(k)

    end do

  end if

  if ( c(n) == 0.0e+00 ) then
    info = n
    return
  end if
!
!  back solve.
!
  b(n) = b(n) / c(n)

  if ( n > 1 ) then

    b(n-1) = ( b(n-1) - d(n-1) * b(n) ) / c(n-1)

    do k = n-2, 1, -1
      b(k) = ( b(k) - d(k) * b(k+1) - e(k) * b(k+2) ) / c(k)
    end do

  end if

  return
end
subroutine spline_b_val ( ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! spline_b_val evaluates a cubic b spline approximant.
!
!
!  discussion:
!
!    the cubic b spline will approximate the data, but is not 
!    designed to interpolate it.
!
!    in effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  modified:
!
!    07 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndata, the number of data values.
!
!    input, real tdata(ndata), the abscissas of the data.
!
!    input, real ydata(ndata), the data values.
!
!    input, real tval, a point at which the spline is to be evaluated.
!
!    output, real yval, the value of the function at tval.
!
  implicit none
!
  integer ndata
!
  real bval
  integer left
  integer right
  real tdata(ndata)
  real tval
  real u
  real ydata(ndata)
  real yval
!
!  find the nearest interval [ tdata(left), tdata(right) ] to tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  evaluate the 5 nonzero b spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
  yval = 0.0e+00
!
!  b function associated with node left - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = ( 1.0e+00 - 3.0e+00 * u + 3.0e+00 * u**2 - u**3 ) / 6.0e+00
  if ( left-1 > 0 ) then
    yval = yval + ydata(left-1) * bval
  else
    yval = yval + ( 2.0e+00 * ydata(1) - ydata(2) ) * bval
  end if
!
!  b function associated with node left,
!  evaluated in its third interval.
!
  bval = ( 4.0e+00 - 6.0e+00 * u**2 + 3.0e+00 * u**3 ) / 6.0e+00
  yval = yval + ydata(left) * bval
!
!  b function associated with node right,
!  evaluated in its second interval.
!
  bval = ( 1.0e+00 + 3.0e+00 * u + 3.0e+00 * u**2 - 3.0e+00 * u**3 ) / 6.0e+00
  yval = yval + ydata(right) * bval
!
!  b function associated with node right+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = u**3 / 6.0e+00
  if ( right+1 <= ndata ) then
    yval = yval + ydata(right+1) * bval
  else
    yval = yval + ( 2.0e+00 * ydata(ndata) - ydata(ndata-1) ) * bval
  end if

  return
end
subroutine spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! spline_beta_val evaluates a cubic beta spline approximant.
!
!
!  discussion:
!
!    the cubic beta spline will approximate the data, but is not 
!    designed to interpolate it.
!
!    if beta1 = 1 and beta2 = 0, the cubic beta spline will be the
!    same as the cubic b spline approximant.
!
!    with beta1 = 1 and beta2 large, the beta spline becomes more like
!    a linear spline.
!
!    in effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  modified:
!
!    12 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, real beta1, the skew or bias parameter.
!    beta1 = 1 for no skew or bias.
!
!    input, real beta2, the tension parameter.
!    beta2 = 0 for no tension.
!
!    input, integer ndata, the number of data values.
!
!    input, real tdata(ndata), the abscissas of the data.
!
!    input, real ydata(ndata), the data values.
!
!    input, real tval, a point at which the spline is to be evaluated.
!
!    output, real yval, the value of the function at tval.
!
  implicit none
!
  integer ndata
!
  real a
  real b
  real beta1
  real beta2
  real bval
  real c
  real d
  real delta
  integer left
  integer right
  real tdata(ndata)
  real tval
  real u
  real ydata(ndata)
  real yval
!
!  find the nearest interval [ tdata(left), tdata(right) ] to tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  evaluate the 5 nonzero beta spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )

  delta = 2.0e+00 + beta2 + 4.0e+00 * beta1 + 4.0e+00 * beta1**2 &
    + 2.0e+00 * beta1**3

  yval = 0.0e+00
!
!  beta function associated with node left - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = ( 2.0e+00 * beta1**3 * ( 1.0e+00 - u )**3 ) / delta

  if ( left-1 > 0 ) then
    yval = yval + ydata(left-1) * bval
  else
    yval = yval + ( 2.0e+00 * ydata(1) - ydata(2) ) * bval
  end if
!
!  beta function associated with node left,
!  evaluated in its third interval.
!
  a = beta2 + 4.0e+00 * beta1 + 4.0e+00 * beta1**2

  b = - 6.0e+00 * beta1 * ( 1.0e+00 - beta1 ) * ( 1.0e+00 + beta1 )

  c = - 3.0e+00 * ( beta2 + 2.0e+00 * beta1**2 + 2.0e+00 * beta1**3 )

  d = 2.0e+00 * ( beta2 + beta1 + beta1**2 + beta1**3 )

  bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta

  yval = yval + ydata(left) * bval
!
!  beta function associated with node right,
!  evaluated in its second interval.
!
  a = 2.0e+00

  b = 6.0e+00 * beta1

  c = 3.0e+00 * beta2 + 6.0e+00 * beta1**2

  d = - 2.0e+00 * ( 1.0e+00 + beta2 + beta1 + beta1**2 )

  bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta

  yval = yval + ydata(right) * bval
!
!  beta function associated with node right+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = 2.0e+00 * u**3 / delta

  if ( right+1 <= ndata ) then
    yval = yval + ydata(right+1) * bval
  else
    yval = yval + ( 2.0e+00 * ydata(ndata) - ydata(ndata-1) ) * bval
  end if

  return
end
subroutine spline_constant_val ( ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! spline_constant_val evaluates a piecewise constant spline at a point.
!
!
!  discussion:
!
!    ndata-1 points tdata define ndata intervals, with the first
!    and last being semi-infinite.
!
!    the value of the spline anywhere in interval i is ydata(i).
!
!  modified:
!
!    16 november 2001
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndata, the number of data points defining the spline.
!
!    input, real tdata(ndata-1), the breakpoints.  the values of tdata should
!    be distinct and increasing.
!
!    input, ydata(ndata), the values of the spline in the intervals
!    defined by the breakpoints.
!
!    input, real tval, the point at which the spline is to be evaluated.
!
!    output, real yval, the value of the spline at tval.  
!
  implicit none
!
  integer ndata
!
  integer i
  real tdata(ndata-1)
  real tval
  real ydata(ndata)
  real yval
!
  do i = 1, ndata-1
    if ( tval <= tdata(i) ) then
      yval = ydata(i)
      return
    end if
  end do

  yval = ydata(ndata)

  return
end
subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )
!
!*******************************************************************************
!
!! spline_cubic_set computes the second derivatives of a cubic spline.
!
!
!  discussion:
!
!    for data interpolation, the user must call spline_cubic_set to 
!    determine the second derivative data, passing in the data to be 
!    interpolated, and the desired boundary conditions.
!
!    the data to be interpolated, plus the spline_cubic_set output, 
!    defines the spline.  the user may then call spline_cubic_val to 
!    evaluate the spline at any point.
!
!    the cubic spline is a piecewise cubic polynomial.  the intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  the cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.  
!
!    for any point t in the interval t(ival), t(ival+1), the form of
!    the spline is
!
!      spl(t) = a(ival)
!             + b(ival) * ( t - t(ival) ) 
!             + c(ival) * ( t - t(ival) )**2
!             + d(ival) * ( t - t(ival) )**3
!
!    if we assume that we know the values y(*) and ypp(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      a(ival) = y(ival)
!      b(ival) = ( y(ival+1) - y(ival) ) / ( t(ival+1) - t(ival) )
!        - ( ypp(ival+1) + 2 * ypp(ival) ) * ( t(ival+1) - t(ival) ) / 6
!      c(ival) = ypp(ival) / 2
!      d(ival) = ( ypp(ival+1) - ypp(ival) ) / ( 6 * ( t(ival+1) - t(ival) ) )
!
!    since the first derivative of the spline is
!
!      spl'(t) =     b(ival)
!              + 2 * c(ival) * ( t - t(ival) )
!              + 3 * d(ival) * ( t - t(ival) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot i results in a total of n-2 equations, of the form:
!
!      b(ival-1) + 2 c(ival-1) * (t(ival)-t(ival-1)) 
!      + 3 * d(ival-1) * (t(ival) - t(ival-1))**2 = b(ival)
!
!    or, setting h(ival) = t(ival+1) - t(ival)
!
!      ( y(ival) - y(ival-1) ) / h(ival-1)
!      - ( ypp(ival) + 2 * ypp(ival-1) ) * h(ival-1) / 6
!      + ypp(ival-1) * h(ival-1)
!      + ( ypp(ival) - ypp(ival-1) ) * h(ival-1) / 2
!      = 
!      ( y(ival+1) - y(ival) ) / h(ival)
!      - ( ypp(ival+1) + 2 * ypp(ival) ) * h(ival) / 6
!
!    or
!
!      ypp(ival-1) * h(ival-1) + 2 * ypp(ival) * ( h(ival-1) + h(ival) )
!      + ypp(ival) * h(ival) 
!      =
!      6 * ( y(ival+1) - y(ival) ) / h(ival)
!      - 6 * ( y(ival) - y(ival-1) ) / h(ival-1)    
!
!    boundary conditions must be applied at the first and last knots.  
!    the resulting tridiagonal system can be solved for the ypp values.
!
!  modified:
!
!    20 november 2000
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the number of data points; n must be at least 2. 
!
!    input, real t(n), the points where data is specified.  
!    the values should be distinct, and increasing.
!
!    input, real y(n), the data values to be interpolated.
!
!    input, integer ibcbeg, the left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be ybcbeg;
!      2: the second derivative at the left endpoint should be ybcbeg.
!
!    input, real ybcbeg, the left boundary value, if needed.
!
!    input, integer ibcend, the right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be ybcend;
!      2: the second derivative at the right endpoint should be ybcend.
!
!    input, real ybcend, the right boundary value, if needed.
!
!    output, real ypp(n), the second derivatives of the cubic spline.
!
  implicit none
!
  integer n
!
  real diag(n)
  integer i
  integer ibcbeg
  integer ibcend
  real sub(2:n)
  real sup(1:n-1)
  real t(n)
  real y(n)
  real ybcbeg
  real ybcend
  real ypp(n)
!
!  check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_cubic_set - fatal error!'
    write ( *, '(a)' ) '  the number of knots must be at least 2.'
    write ( *, '(a,i6)' ) '  the input value of n = ', n
    stop
  end if

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'spline_cubic_set - fatal error!'
      write ( *, '(a)' ) '  the knots must be strictly increasing, but'
      write ( *, '(a,i6,a,g14.6)' ) '  t(',  i,') = ', t(i)
      write ( *, '(a,i6,a,g14.6)' ) '  t(',i+1,') = ', t(i+1)
      stop
    end if
  end do
!
!  set the first equation.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.0e+00
    diag(1) = 1.0e+00
    sup(1) = -1.0e+00
  else if ( ibcbeg == 1 ) then
    ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    diag(1) = ( t(2) - t(1) ) / 3.0e+00 
    sup(1) = ( t(2) - t(1) ) / 6.0e+00
  else if ( ibcbeg == 2 ) then
    ypp(1) = ybcbeg
    diag(1) = 1.0e+00
    sup(1) = 0.0e+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_cubic_set - fatal error!'
    write ( *, '(a)' ) '  the boundary flag ibcbeg must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  the input value is ibcbeg = ', ibcbeg
    stop
  end if
!
!  set the intermediate equations.
!
  do i = 2, n-1
    ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
           - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    sub(i) = ( t(i) - t(i-1) ) / 6.0e+00
    diag(i) = ( t(i+1) - t(i-1) ) / 3.0e+00
    sup(i) = ( t(i+1) - t(i) ) / 6.0e+00
  end do
!
!  set the last equation.
!
  if ( ibcend == 0 ) then
    ypp(n) = 0.0e+00
    sub(n) = -1.0e+00
    diag(n) = 1.0e+00
  else if ( ibcend == 1 ) then
    ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    sub(n) = ( t(n) - t(n-1) ) / 6.0e+00
    diag(n) = ( t(n) - t(n-1) ) / 3.0e+00
  else if ( ibcend == 2 ) then
    ypp(n) = ybcend
    sub(n) = 0.0e+00
    diag(n) = 1.0e+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_cubic_set - fatal error!'
    write ( *, '(a)' ) '  the boundary flag ibcend must be 0, 1 or 2.'
    write ( *, '(a,i6)' ) '  the input value is ibcend = ', ibcend
    stop
  end if
!
!  special case:
!    n = 2, ibcbeg = ibcend = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0e+00
    ypp(2) = 0.0e+00
!
!  solve the linear system.
!
  else

    call s3_fs ( sub, diag, sup, n, ypp, ypp )

  end if

  return
end
subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )
!
!*******************************************************************************
!
!! spline_cubic_val evaluates a cubic spline at a specific point.
!
!
!  discussion:
!
!    spline_cubic_set must have already been called to define the 
!    values of ypp.
!
!    for any point t in the interval t(ival), t(ival+1), the form of
!    the spline is
!
!      spl(t) = a 
!             + b * ( t - t(ival) ) 
!             + c * ( t - t(ival) )**2
!             + d * ( t - t(ival) )**3
!
!    here:
!      a = y(ival)
!      b = ( y(ival+1) - y(ival) ) / ( t(ival+1) - t(ival) )
!        - ( ypp(ival+1) + 2 * ypp(ival) ) * ( t(ival+1) - t(ival) ) / 6
!      c = ypp(ival) / 2
!      d = ( ypp(ival+1) - ypp(ival) ) / ( 6 * ( t(ival+1) - t(ival) ) )
!
!  modified:
!
!    20 november 2000
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the number of data values.
!
!    input, real t(n), the knot values.
!
!    input, real y(n), the data values at the knots.
!
!    input, real ypp(n), the second derivatives of the spline at the knots.
!
!    input, real tval, a point, typically between t(1) and t(n), at 
!    which the spline is to be evalulated.  if tval lies outside 
!    this range, extrapolation is used.
!
!    output, real yval, ypval, yppval, the value of the spline, and
!    its first two derivatives at tval.
!
  implicit none
!
  integer n
!
  real dt
  real h
  integer left
  integer right
  real t(n)
  real tval
  real y(n)
  real ypp(n)
  real yppval
  real ypval
  real yval
!
!  determine the interval [t(left), t(right)] that contains tval.
!  values below t(1) or above t(n) use extrapolation.
!
  call rvec_bracket ( n, t, tval, left, right )
!
!  evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) &
       + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0e+00 + ypp(left) / 3.0e+00 ) * h &
       + dt * ( 0.5e+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0e+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
       - ( ypp(right) / 6.0e+00 + ypp(left) / 3.0e+00 ) * h &
       + dt * ( ypp(left) &
       + dt * ( 0.5e+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h 

  return
end
subroutine spline_cubic_val2 ( n, t, y, ypp, left, tval, yval, ypval, yppval )
!
!*******************************************************************************
!
!! spline_cubic_val2 evaluates a cubic spline at a specific point.
!
!
!  discussion:
!
!    this routine is a modification of spline_cubic_val; it allows the
!    user to speed up the code by suggesting the appropriate t interval
!    to search first.
!
!    spline_cubic_set must have already been called to define the
!    values of ypp.
!
!    in the left interval, let right = left+1.  the form of the spline is
!
!      spl(t) = 
!          a
!        + b * ( t - t(left) )
!        + c * ( t - t(left) )**2
!        + d * ( t - t(left) )**3
!
!    here:
!      a = y(left)
!      b = ( y(right) - y(left) ) / ( t(right) - t(left) )
!        - ( ypp(right) + 2 * ypp(left) ) * ( t(right) - t(left) ) / 6
!      c = ypp(left) / 2
!      d = ( ypp(right) - ypp(left) ) / ( 6 * ( t(right) - t(left) ) )
!
!  modified:
!
!    20 november 2000
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer n, the number of knots.
!
!    input, real t(n), the knot values.
!
!    input, real y(n), the data values at the knots.
!
!    input, real ypp(n), the second derivatives of the spline at
!    the knots.
!
!    input/output, integer left, the suggested t interval to search.
!    left should be between 1 and n-1.  if left is not in this range,
!    then its value will be ignored.  on output, left is set to the
!    actual interval in which tval lies.
!
!    input, real tval, a point, typically between t(1) and t(n), at
!    which the spline is to be evalulated.  if tval lies outside
!    this range, extrapolation is used.
!
!    output, real yval, ypval, yppval, the value of the spline, and
!    its first two derivatives at tval.
!
  implicit none
!
  integer n
!
  real dt
  real h
  integer left
  integer right
  real t(n)
  real tval
  real y(n)
  real ypp(n)
  real yppval
  real ypval
  real yval
!
!  determine the interval [t(left), t(right)] that contains tval.
!  
!  what you want from rvec_bracket3 is that tval is to be computed
!  by the data in interval {t(left), t(right)].
!
  left = 0
  call rvec_bracket3 ( n, t, tval, left )
  right = left + 1
!
!  in the interval left, the polynomial is in terms of a normalized
!  coordinate  ( dt / h ) between 0 and 1.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0e+00 + ypp(left) / 3.0e+00 ) * h &
       + dt * ( 0.5e+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0e+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
      - ( ypp(right) / 6.0e+00 + ypp(left) / 3.0e+00 ) * h &
      + dt * ( ypp(left) &
      + dt * ( 0.5e+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

  return
end
subroutine spline_hermite_set ( ndata, tdata, c )
!
!*************************************************************************
!
!! spline_hermite_set sets up a piecewise cubic hermite interpolant.
!
!
!  reference:
!
!    conte and de boor,
!    algorithm calccf,
!    elementary numerical analysis,
!    1973, page 235.
!
!  modified:
!
!    06 april 1999
!
!  parameters:
!
!    input, integer ndata, the number of data points.
!    ndata must be at least 2.
!
!    input, real tdata(ndata), the abscissas of the data points.
!    the entries of tdata are assumed to be strictly increasing.
!
!    input/output, real c(4,ndata).
!
!    on input, c(1,i) and c(2,i) should contain the value of the
!    function and its derivative at tdata(i), for i = 1 to ndata.
!    these values will not be changed by this routine.
!
!    on output, c(3,i) and c(4,i) contain the quadratic
!    and cubic coefficients of the hermite polynomial
!    in the interval (tdata(i), tdata(i+1)), for i=1 to ndata-1.
!    c(3,ndata) and c(4,ndata) are set to 0.
!
!    in the interval (tdata(i), tdata(i+1)), the interpolating hermite
!    polynomial is given by
!
!    sval(tval) =                 c(1,i)
!       + ( tval - tdata(i) ) * ( c(2,i)
!       + ( tval - tdata(i) ) * ( c(3,i)
!       + ( tval - tdata(i) ) *   c(4,i) ) )
!
  implicit none
!
  integer ndata
!
  real c(4,ndata)
  real divdif1
  real divdif3
  real dt
  integer i
  real tdata(ndata)
!
  do i = 1, ndata-1
    dt = tdata(i+1) - tdata(i)
    divdif1 = ( c(1,i+1) - c(1,i) ) / dt
    divdif3 = c(2,i) + c(2,i+1) - 2.0e+00 * divdif1
    c(3,i) = ( divdif1 - c(2,i) - divdif3 ) / dt
    c(4,i) = divdif3 / dt**2
  end do

  c(3,ndata) = 0.0e+00
  c(4,ndata) = 0.0e+00

  return
end
subroutine spline_hermite_val ( ndata, tdata, c, tval, sval )
!
!*************************************************************************
!
!! spline_hermite_val evaluates a piecewise cubic hermite interpolant.
!
!
!  discussion:
!
!    spline_hermite_set must be called first, to set up the
!    spline data from the raw function and derivative data.
!
!  reference:
!
!    conte and de boor,
!    algorithm pcubic,
!    elementary numerical analysis,
!    1973, page 234.
!
!  modified:
!
!    06 april 1999
!
!  parameters:
!
!    input, integer ndata, the number of data points.
!    ndata is assumed to be at least 2.
!
!    input, real tdata(ndata), the abscissas of the data points.
!    the entries of tdata are assumed to be strictly increasing.
!
!    input, real c(4,ndata), contains the data computed by
!    spline_hermite_set.
!
!    input, real tval, the point where the interpolant is to
!    be evaluated.
!
!    output, real sval, the value of the interpolant at tval.
!
  implicit none
!
  integer ndata
!
  real c(4,ndata)
  real dt
  integer left
  integer right
  real sval
  real tdata(ndata)
  real tval
!
!  find the interval [ tdata(left), tdata(right) ] that contains
!  or is nearest to tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  evaluate the cubic polynomial.
!
  dt = tval - tdata(left)

  sval = c(1,left) + dt * ( c(2,left) + dt * ( c(3,left) + dt * c(4,left) ) )

  return
end
subroutine spline_linear_int ( ndata, tdata, ydata, a, b, int_val )
!
!*******************************************************************************
!
!! spline_linear_int evaluates the integral of a linear spline.
!
!
!  modified:
!
!    01 november 2001
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndata, the number of data points defining the spline.
!
!    input, real tdata(ndata), ydata(ndata), the values of the independent
!    and dependent variables at the data points.  the values of tdata should
!    be distinct and increasing.
!
!    input, real a, b, the interval over which the integral is desired.
!
!    output, real int_val, the value of the integral.
!
  implicit none
!
  integer ndata
!
  real a
  real a_copy
  integer a_left
  integer a_right
  real b
  real b_copy
  integer b_left
  integer b_right
  integer i_left
  real int_val
  real tdata(ndata)
  real tval
  real ydata(ndata)
  real yp
  real yval
!
  int_val = 0.0e+00 

  if ( a == b ) then
    return
  end if

  a_copy = min ( a, b )
  b_copy = max ( a, b )
!
!  find the interval [ tdata(a_left), tdata(a_right) ] that contains, or is
!  nearest to, a.
!
  call rvec_bracket ( ndata, tdata, a_copy, a_left, a_right )
!
!  find the interval [ tdata(b_left), tdata(b_right) ] that contains, or is
!  nearest to, b.
!
  call rvec_bracket ( ndata, tdata, b_copy, b_left, b_right )
!
!  if a and b are in the same interval...
!
  if ( a_left == b_left ) then

    tval = ( a_copy + b_copy ) / 2.0e+00

    yp = ( ydata(a_right) - ydata(a_left) ) / &
         ( tdata(a_right) - tdata(a_left) )

    yval = ydata(a_left) + ( tval - tdata(a_left) ) * yp

    int_val = yval * ( b_copy - a_copy )

    return
  end if
!
!  otherwise, integrate from:
!
!  a               to tdata(a_right),
!  tdata(a_right)  to tdata(a_right+1),...
!  tdata(b_left-1) to tdata(b_left),
!  tdata(b_left)   to b.
!
!  use the fact that the integral of a linear function is the
!  value of the function at the midpoint times the width of the interval.
!
  tval = ( a_copy + tdata(a_right) ) / 2.0e+00

  yp = ( ydata(a_right) - ydata(a_left) ) / &
       ( tdata(a_right) - tdata(a_left) )

  yval = ydata(a_left) + ( tval - tdata(a_left) ) * yp

  int_val = int_val + yval * ( tdata(a_right) - a_copy )

  do i_left = a_right, b_left - 1

    tval = ( tdata(i_left+1) + tdata(i_left) ) / 2.0e+00

    yp = ( ydata(i_left+1) - ydata(i_left) ) / &
         ( tdata(i_left+1) - tdata(i_left) )

    yval = ydata(i_left) + ( tval - tdata(i_left) ) * yp

    int_val = int_val + yval * ( tdata(i_left + 1) - tdata(i_left) )

  end do

  tval = ( tdata(b_left) + b_copy ) / 2.0e+00

  yp = ( ydata(b_right) - ydata(b_left) ) / &
       ( tdata(b_right) - tdata(b_left) )

  yval = ydata(b_left) + ( tval - tdata(b_left) ) * yp

  int_val = int_val + yval * ( b_copy - tdata(b_left) )

  if ( b < a ) then
    int_val = - int_val
  end if

  return
end
subroutine spline_linear_intset ( int_n, int_x, int_v, data_n, data_x, data_y )
!
!*******************************************************************************
!
!! spline_linear_intset sets a linear spline with given integral properties.
!
!
!  discussion:
!
!    the user has in mind an interval, divided by int_n+1 points into
!    int_n intervals.  a linear spline is to be constructed,
!    with breakpoints at the centers of each interval, and extending
!    continuously to the left of the first and right of the last
!    breakpoints.  the constraint on the linear spline is that it is
!    required that it have a given integral value over each interval.
!
!    a tridiagonal linear system of equations is solved for the
!    values of the spline at the breakpoints.
!
!  modified:
!
!    02 november 2001
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer int_n, the number of intervals.
!
!    input, real int_x(int_n+1), the points that define the intervals.
!    interval i lies between int_x(i) and int_x(i+1).
!
!    input, real int_v(int_n), the desired value of the integral of the
!    linear spline over each interval.
!
!    output, integer data_n, the number of data points defining the spline.
!    (this is the same as int_n).
!
!    output, real data_x(data_n), data_y(data_n), the values of the independent
!    and dependent variables at the data points.  the values of data_x are
!    the interval midpoints.  the values of data_y are determined in such
!    a way that the exact integral of the linear spline over interval i
!    is equal to int_v(i).
!
  implicit none
!
  integer int_n
!
  real c(int_n)
  real d(int_n)
  integer data_n
  real data_x(int_n)
  real data_y(int_n)
  real e(int_n)
  integer info
  real int_v(int_n)
  real int_x(int_n+1)
!
!  set up the easy stuff.
!
  data_n = int_n
  data_x(1:data_n) = 0.5e+00 * ( int_x(1:data_n) + int_x(2:data_n+1) )
!
!  set up c, d, e, the coefficients of the linear system.
!
  c(1) = 0.0e+00
  c(2:data_n-1) = 1.0 &
    - ( 0.5 * ( data_x(2:data_n-1) + int_x(2:data_n-1) ) &
    - data_x(1:data_n-2) ) &
    / ( data_x(2:data_n-1) - data_x(1:data_n-2) )
  c(data_n) = 0.0e+00

  d(1) = int_x(2) - int_x(1)

  d(2:data_n-1) = 1.0 &
    + ( 0.5 * ( data_x(2:data_n-1) + int_x(2:data_n-1) ) &
    - data_x(1:data_n-2) ) &
    / ( data_x(2:data_n-1) - data_x(1:data_n-2) ) &
    - ( 0.5 * ( data_x(2:data_n-1) + int_x(3:data_n) ) - data_x(2:data_n-1) ) &
    / ( data_x(3:data_n) - data_x(2:data_n-1) )

  d(data_n) = int_x(data_n+1) - int_x(data_n)

  e(1) = 0.0e+00

  e(2:data_n-1) = ( 0.5 * ( data_x(2:data_n-1) + int_x(3:data_n) ) &
    - data_x(2:data_n-1) ) / ( data_x(3:data_n) - data_x(2:data_n-1) )

  e(data_n) = 0.0e+00
!
!  set up data_y, which begins as the right hand side of the linear system.
!
  data_y(1) = int_v(1)
  data_y(2:data_n-1) = 2.0e+00 * int_v(2:data_n-1) &
    / ( int_x(3:int_n) - int_x(2:int_n-1) )
  data_y(data_n) = int_v(data_n)
!
!  solve the linear system.
!
  call sgtsl ( data_n, c, d, e, data_y, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_linear_intset - fatal error!'
    write ( *, '(a)' ) '  the linear system is singular.'
    stop
  end if

  return
end
subroutine spline_linear_val ( ndata, tdata, ydata, tval, yval, ypval )
!
!*******************************************************************************
!
!! spline_linear_val evaluates a linear spline at a specific point.
!
!
!  discussion:
!
!    because of the extremely simple form of the linear spline,
!    the raw data points ( tdata(i), ydata(i)) can be used directly to
!    evaluate the spline at any point.  no processing of the data
!    is required.
!
!  modified:
!
!    06 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndata, the number of data points defining the spline.
!
!    input, real tdata(ndata), ydata(ndata), the values of the independent
!    and dependent variables at the data points.  the values of tdata should
!    be distinct and increasing.
!
!    input, real tval, the point at which the spline is to be evaluated.
!
!    output, real yval, ypval, the value of the spline and its first
!    derivative dydt at tval.  ypval is not reliable if tval is exactly
!    equal to tdata(i) for some i.
!
  implicit none
!
  integer ndata
!
  integer left
  integer right
  real tdata(ndata)
  real tval
  real ydata(ndata)
  real ypval
  real yval
!
!  find the interval [ tdata(left), tdata(right) ] that contains, or is
!  nearest to, tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  now evaluate the piecewise linear function.
!
  ypval = ( ydata(right) - ydata(left) ) / ( tdata(right) - tdata(left) )

  yval = ydata(left) +  ( tval - tdata(left) ) * ypval

  return
end
subroutine spline_overhauser_nonuni_val ( ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! spline_overhauser_nonuni_val evaluates the nonuniform overhauser spline.
!
!
!  discussion:
!
!    the nonuniformity refers to the fact that the abscissas values
!    need not be uniformly spaced.
!
!  diagnostics:
!
!    the values of alpha and beta have to be properly assigned.
!    the basis matrices for the first and last interval have to
!    be computed.
!
!  modified:
!
!    06 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndata, the number of data points.
!
!    input, real tdata(ndata), the abscissas of the data points.
!    the values of tdata are assumed to be distinct and increasing.
!
!    input, real ydata(ndata), the data values.
!
!    input, real tval, the value where the spline is to
!    be evaluated.
!
!    output, real yval, the value of the spline at tval.
!
  implicit none
!
  integer ndata
!
  real alpha
  real beta
  integer left
  real mbasis(4,4)
  real mbasis_l(3,3)
  real mbasis_r(3,3)
  integer right
  real tdata(ndata)
  real tval
  real ydata(ndata)
  real yval
!
!  find the nearest interval [ tdata(left), tdata(right) ] to tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  evaluate the spline in the given interval.
!
  if ( left == 1 ) then

    alpha = 1.0e+00
    call basis_matrix_overhauser_nul ( alpha, mbasis_l )

    call basis_matrix_tmp ( 1, 3, mbasis_l, ndata, tdata, ydata, tval, yval )

  else if ( left < ndata-1 ) then

    alpha = 1.0e+00
    beta = 1.0e+00
    call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

    call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval, yval )

  else if ( left == ndata-1 ) then

    beta = 1.0e+00
    call basis_matrix_overhauser_nur ( beta, mbasis_r )

    call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, ydata, tval, yval )

  end if

  return
end
subroutine spline_overhauser_uni_val ( ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! spline_overhauser_uni_val evaluates the uniform overhauser spline.
!
!
!  modified:
!
!    06 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndata, the number of data points.
!
!    input, real tdata(ndata), the abscissas of the data points.
!    the values of tdata are assumed to be distinct and increasing.
!    this routine also assumes that the values of tdata are uniformly
!    spaced; for instance, tdata(1) = 10, tdata(2) = 11, tdata(3) = 12...
!
!    input, real ydata(ndata), the data values.
!
!    input, real tval, the value where the spline is to
!    be evaluated.
!
!    output, real yval, the value of the spline at tval.
!
  implicit none
!
  integer ndata
!
  integer left
  real mbasis(4,4)
  real mbasis_l(3,3)
  real mbasis_r(3,3)
  integer right
  real tdata(ndata)
  real tval
  real ydata(ndata)
  real yval
!
!  find the nearest interval [ tdata(left), tdata(right) ] to tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  evaluate the spline in the given interval.
!
  if ( left == 1 ) then

    call basis_matrix_overhauser_uni_l ( mbasis_l )

    call basis_matrix_tmp ( 1, 3, mbasis_l, ndata, tdata, ydata, tval, yval )

  else if ( left < ndata-1 ) then

    call basis_matrix_overhauser_uni ( mbasis )

    call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval, yval )

  else if ( left == ndata-1 ) then

    call basis_matrix_overhauser_uni_r ( mbasis_r )

    call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, ydata, tval, yval )

  end if

  return
end
subroutine spline_overhauser_val ( ndim, ndata, tdata, ydata, tval, yval )
!
!*******************************************************************************
!
!! spline_overhauser_val evaluates an overhauser spline.
!
!
!  discussion:
!
!    over the first and last intervals, the overhauser spline is a 
!    quadratic.  in the intermediate intervals, it is a piecewise cubic.
!    the overhauser spline is also known as the catmull-rom spline.
!
!  reference:
!
!    h brewer and d anderson,
!    visual interaction with overhauser curves and surfaces,
!    siggraph 77, pages 132-137.
!
!    e catmull and r rom,
!    a class of local interpolating splines,
!    in computer aided geometric design,
!    edited by r barnhill and r reisenfeld,
!    academic press, 1974, pages 317-326.
!
!    david rogers and alan adams,
!    mathematical elements of computer graphics,
!    mcgraw hill, 1990, second edition, pages 278-289.
!
!  modified:
!
!   08 april 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndim, the dimension of a single data point.
!    ndim must be at least 1.  there is an internal limit on ndim,
!    called maxdim, which is presently set to 5.
!
!    input, integer ndata, the number of data points.
!    ndata must be at least 3.
!
!    input, real tdata(ndata), the abscissas of the data points.  the
!    values in tdata must be in strictly ascending order.
!
!    input, real ydata(ndim,ndata), the data points corresponding to
!    the abscissas.
!
!    input, real tval, the abscissa value at which the spline
!    is to be evaluated.  normally, tdata(1) <= tval <= t(ndata), and 
!    the data will be interpolated.  for tval outside this range, 
!    extrapolation will be used.
!
!    output, real yval(ndim), the value of the spline at tval.
!
  implicit none
!
  integer, parameter :: maxdim = 5
  integer ndata
  integer ndim
!
  integer i
  integer left
  integer order
  integer right
  real tdata(ndata)
  real tval
  real ydata(ndim,ndata)
  real yl(maxdim)
  real yr(maxdim)
  real yval(ndim)
!
!  check.
!
  call rvec_order_type ( ndata, tdata, order )

  if ( order /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_overhauser_val - fatal error!'
    write ( *, '(a)' ) '  the data abscissas are not strictly ascending.'
    stop
  end if

  if ( ndata < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_overhauser_val - fatal error!'
    write ( *, '(a)' ) '  ndata < 3.'
    stop
  end if

  if ( ndim < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_overhauser_val - fatal error!'
    write ( *, '(a)' ) '  ndim < 1.'
    stop
  end if

  if ( ndim > maxdim ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_overhauser_val - fatal error!'
    write ( *, '(a)' ) '  ndim > maxdim.'
    stop
  end if
!
!  locate the abscissa interval t(left), t(left+1) nearest to or 
!  containing tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  evaluate the "left hand" quadratic defined at t(left-1), t(left), t(right).
!
  if ( left-1 > 0 ) then
    call parabola_val2 ( ndim, ndata, tdata, ydata, left-1, tval, yl )
  end if
!
!  evaluate the "right hand" quadratic defined at t(left), t(right), t(right+1).
!
  if ( right+1 <= ndata ) then
    call parabola_val2 ( ndim, ndata, tdata, ydata, left, tval, yr )
  end if
!
!  average the quadratics.
!
  if ( left == 1 ) then

    yval(1:ndim) = yr(1:ndim)

  else if ( right < ndata ) then

    yval(1:ndim) =  ( ( tdata(right) - tval ) * yl(1:ndim) &
      + ( tval - tdata(left) ) * yr(1:ndim) ) / ( tdata(right) - tdata(left) )

  else

    yval(1:ndim) = yl(1:ndim)

  end if

  return
end
subroutine spline_quadratic_val ( ndata, tdata, ydata, tval, yval, ypval )
!
!*******************************************************************************
!
!! spline_quadratic_val evaluates a quadratic spline at a specific point.
!
!
!  discussion:
!
!    because of the simple form of a piecewise quadratic spline,
!    the raw data points ( tdata(i), ydata(i)) can be used directly to
!    evaluate the spline at any point.  no processing of the data
!    is required.
!
!  modified:
!
!    24 october 1999
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, integer ndata, the number of data points defining the spline.
!    ndata should be odd.
!
!    input, real tdata(ndata), ydata(ndata), the values of the independent
!    and dependent variables at the data points.  the values of tdata should
!    be distinct and increasing.
!
!    input, real tval, the point at which the spline is to be evaluated.
!
!    output, real yval, ypval, the value of the spline and its first
!    derivative dydt at tval.  ypval is not reliable if tval is exactly
!    equal to tdata(i) for some i.
!
  implicit none
!
  integer ndata
!
  real dif1
  real dif2
  integer left
  integer right
  real t1
  real t2
  real t3
  real tdata(ndata)
  real tval
  real y1
  real y2
  real y3
  real ydata(ndata)
  real ypval
  real yval
!
  if ( mod ( ndata, 3 ) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_quadratic_val - fatal error!'
    write ( *, '(a)' ) '  ndata must be odd.'
    stop
  end if
!
!  find the interval [ tdata(left), tdata(right) ] that contains, or is
!  nearest to, tval.
!
  call rvec_bracket ( ndata, tdata, tval, left, right )
!
!  force left to be odd.
!
  if ( mod ( left, 2 ) == 0 ) then
    left = left - 1
  end if
!
!  copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t1 >= t2 .or. t2 >= t3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'spline_quadratic_val - fatal error!'
    write ( *, '(a)' ) '  t1 >= t2 or t2 >= t3.'
    stop
  end if
!
!  construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  y1 = ydata(left)
  y2 = ydata(left+1)
  y3 = ydata(left+2)

  dif1 = ( y2 - y1 ) / ( t2 - t1 )

  dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
       - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

  yval = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )
  ypval = dif1 + dif2 * ( 2.0e+00 * tval - t1 - t2 )

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! timestamp prints the current ymdhms date as a time stamp.
!
!
!  example:
!
!    may 31 2001   9:45:54.872 am
!
!  modified:
!
!    31 may 2001
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    none
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'january  ', 'february ', 'march    ', 'april    ', &
    'may      ', 'june     ', 'july     ', 'august   ', &
    'september', 'october  ', 'november ', 'december ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'am'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'noon'
    else
      ampm = 'pm'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'pm'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'midnight'
      else
        ampm = 'am'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
