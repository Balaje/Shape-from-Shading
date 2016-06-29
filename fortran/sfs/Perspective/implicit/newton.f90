  ! Module to run Newton Raphson Method

module newton

  implicit none
  real(kind = 8), parameter :: tol = 1e-16
  integer, parameter :: maxiter = 1000

contains

  real(kind = 8) function solve(L,delt)
    real(kind = 8), intent(in) :: L, delt
    real(kind = 8) :: x0, error, xn, f, df
    integer :: iterations = 0
    x0 = 0.
    error = 100.
    do while (error > tol .or. iterations == maxiter)
       f = (x0 - delt*exp(-2*x0) - L)
       df = 1 + 2*delt*exp(-2*x0)

       xn = x0 - f/df

       error = abs(xn-x0)

       x0 = xn

       iterations = iterations + 1
    enddo

    solve = xn

  end function solve

end module newton

