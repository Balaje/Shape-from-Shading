! Program to solve the Burger's equation
!  u_t + f(u)_x = 0
!   subject to the initial data
!          u0(x) = { 0, x<0
!                    1, x>0
!          }
! Module file : Contains the declaration of the flux function and the numerical fluxes



module fun
  implicit none

contains

  !flux funtion f
  real (kind = 8) function f(u)
    real(kind = 8), intent(in) :: u

    f = u**2/2
  end function f



  !--------Numerical flux functions
  !--- Godunov Flux
  real (kind=8) function fluxgod(f,u,v)
    real (kind = 8), external :: f
    real (kind = 8), intent(in) :: u,v

    fluxgod = max(f(max(u,0.)), f(min(v,0.)))
  end function fluxgod

  !--- Lax Friedrich flux
  real (kind = 8) function fluxlax(f,u,v,eps)
    real (kind = 8), external :: f
    real (kind = 8), intent(in) :: u,v,eps

    fluxlax = 0.5*(f(u)+f(v) - 1/eps*(v-u))
  end function fluxlax
end module fun
