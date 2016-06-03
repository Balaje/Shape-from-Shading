! Module containing all the flux definied
! > Godunov
! > Lax - Friedrich
! > Murman - Roe
! > Einquist-Osher

module flux

contains
! Godunov Flux
  real (kind = 8) function fgod(f,u,v)

    implicit none
    real (kind = 8), external :: f
    real(kind = 8), intent(in) :: u,v
    real(kind = 8) :: left, right

    left = f(max(u,0.))
    right = f(min(v,0.))

    fgod = max(left, right)

  end function fgod

! Lax-Friedrich Flux

  real (kind = 8) function flax(f,u,v,r)

    implicit none
    real(kind = 8), external :: f
    real(kind = 8), intent(in) :: u,v,r

    flax = 0.5*(f(u) + f(v) - (1/r)*(v-u))
  end function flax

end module flux
