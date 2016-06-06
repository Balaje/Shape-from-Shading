  ! Function that contains the flux function and Numerical flux

module func

contains

  !--Flux function--convex
  real (kind = 8) function f(u)
    real(kind = 8), intent(in) :: u
    f = abs(u)
  end function f

  !--Flux function--concave
  real (kind=8) function f1(u)
    real(kind = 8), intent(in) :: u
    f1 = -abs(u)
  end function f1

  !--Godunov Flux Convex
  real(kind = 8) function flux(f,u,v)
    real(kind = 8), external :: f
    real(kind = 8), intent(in) :: u,v
    !Convex flux
    flux = max(f(max(u,0.)),f(min(v,0.)))
  end function flux

  !--Godunov Flux Concave
  real(kind=8) function flux1(f1,u,v)
    real(kind=8), external :: f1
    real(kind=8), intent(in) :: u,v
    !Concave Flux
    flux1 = min(f1(min(u,0.)),f1(max(v,0.)))
  end function flux1

end module func
