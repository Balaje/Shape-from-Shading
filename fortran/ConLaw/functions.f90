! Fortran module that contains functions
!   > Flux f(u)
!   > Derivative f'(u)


module functions

contains

  real (kind = 8) function f(x)
    implicit none
    real(kind = 8), intent(in) :: x

    f = x**2/2
  end function f

end module functions
