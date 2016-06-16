  ! Fortran code to implement shape from shading
  !     Orthographic Model by Rouy and Tourin


program ortho
  implicit none

  integer, parameter :: Nx = 500, Ny = 500
  real (kind = 8), parameter :: a = 0, b = 1, c = 0, d = 1, tol = 1e-16
  real (kind = 8), dimension(Nx+1,Ny+1) :: intensity, unew, u
  real (kind = 8), dimension(Nx+1) :: x
  real (kind = 8), dimension(Ny+1) :: y
  real (kind = 8) :: error, H, hx, hy, delt, Dx, Dy, Dxp, Dyp, Dxm, Dym, so

  integer :: i,j, iterations = 0

  hx = (b-a)/Nx;
  hy = (d-c)/Ny;

  ! Find out the x and y vectors
  do i = 1,Nx+1
     x(i) = a+(i-1)*hx
  enddo

  do j = 1,Ny+1
     y(j) = c+(j-1)*hy
  enddo

  open(unit = 1, file="vase.txt")
  read(1,*) intensity

  do j = 1,Ny+1
     do i = 1,Nx+1
        u(i,j) = -0.5*log(intensity(j,i));
        unew(i,j) = 0.
     enddo
  enddo

  delt = hx*hy/sqrt(hx**2+hy**2)
  error = 100.

  do while (error > tol)
     !-- i = 1, j = 1
     Dxp = (u(2,1)-u(1,1))/hx
     Dxm = 0.
     Dyp = (u(1,2)-u(1,1))/hx
     Dym = 0.

     Dx = min(0.,Dxp,Dxm)
     Dy = min(0.,Dyp,Dym)

     so = sqrt(1/intensity(1,1)**2 - 1)
     H = sqrt(Dx**2+Dy**2)

     unew(1,1) = u(1,1) - delt*(H-so)

     !-- i = 1, j = Ny+1
     Dxp = (u(2,Ny+1)-u(1,Ny+1))/hx
     Dxm = 0.
     Dyp = 0.
     Dym = (u(1,Ny)-u(1,Ny+1))/hy

     Dx = min(0.,Dxp,Dxm)
     Dy = min(0.,Dyp,Dym)

     so = sqrt(1/intensity(1,Ny+1)**2 - 1)
     H = sqrt(Dx**2+Dy**2)

     unew(1,Ny+1) = u(1,Ny+1) - delt*(H-so)

     !-- Most interior loop
     do j = 2,Ny
        do i = 2,Nx
           Dxp = (u(i+1,j)-u(i,j))/hx
           Dxm = (u(i-1,j)-u(i,j))/hy
           Dyp = (u(i,j+1)-u(i,j))/hx
           Dym = (u(i,j-1)-u(i,j))/hy

           Dx = min(0.,Dxp,Dxm)
           Dy = min(0.,Dyp,Dym)

           so = sqrt(1/intensity(i,j)**2 - 1)
           H = sqrt(Dx**2+Dy**2)

           unew(i,j) = u(i,j) - delt*(H-so)
        enddo
     enddo

     !-- i = Nx+1, j = 1
     Dxp = 0.
     Dxm = (u(Nx,1)-u(Nx+1,1))/hy
     Dyp = (u(Nx+1,2)-u(Nx+1,1))/hx
     Dym = 0.

     Dx = min(0.,Dxp,Dxm)
     Dy = min(0.,Dyp,Dym)

     so = sqrt(1/intensity(Nx+1,1)**2 - 1)
     H = sqrt(Dx**2+Dy**2)

     unew(Nx+1,1) = u(Nx+1,1) - delt*(H-so)

     !-- i = Nx+1, j = Ny+1
     Dxp = 0.
     Dxm = (u(Nx,Ny+1)-u(Nx+1,Ny+1))/hy
     Dyp = 0.
     Dym = (u(Nx+1,Ny)-u(Nx+1,Ny+1))/hy

     Dx = min(0.,Dxp,Dxm)
     Dy = min(0.,Dyp,Dym)

     so = sqrt(1/intensity(Nx+1,Ny+1)**2 - 1)
     H = sqrt(Dx**2+Dy**2)

     unew(Nx+1,Ny+1) = u(Nx+1,Ny+1) - delt*(H-so)

     !-- j = 2,Nx i = 1
     do j = 2,Ny
        Dxp = (u(2,j)-u(1,j))/hx
        Dxm = 0.
        Dyp = (u(1,j+1)-u(1,j))/hy
        Dym = (u(1,j-1)-u(1,j))/hy

        Dx = min(0.,Dxp,Dxm)
        Dy = min(0.,Dyp,Dym)

        so = sqrt(1/intensity(1,j)**2 - 1)
        H = sqrt(Dx**2+Dy**2)

        unew(1,j) = u(1,j) - delt*(H-so)
     enddo

     !-- j  = 2,Ny i = Nx+1
     do j = 2,Ny
        Dxp = 0.
        Dxm = (u(Nx,j)-u(Nx+1,j))/hx
        Dyp = (u(Nx+1,j+1)-u(Nx+1,j))/hy
        Dym = (u(Nx+1,j-1)-u(Nx+1,j))/hy

        Dx = min(0.,Dxp,Dxm)
        Dy = min(0.,Dyp,Dym)

        so = sqrt(1/intensity(Nx+1,j)**2 - 1)
        H = sqrt(Dx**2+Dy**2)

        unew(Nx+1,j) = u(Nx+1,j) - delt*(H-so)
     enddo

     error = maxval(abs(u-unew))
!     print 11, error
!11   format('Error', se22.15)
     u = unew

     iterations = iterations + 1
  enddo

  open(unit = 2, file="op.txt", status="unknown", action="write")

  do j = 1,Ny+1
     do i = 1,Nx+1
        write(unit = 2, fmt = *) x(i), y(j), unew(i,j)
     enddo
  enddo

end program ortho
