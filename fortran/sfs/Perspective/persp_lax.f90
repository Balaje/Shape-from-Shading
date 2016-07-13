! Program to implement the upwinded Godunov Type scheme to
!    solve the perspective shape from shading model
!      introduced by Prados and Faugeras
!
!
! Authors : Balaje and Sanath Keshav
program persp1
implicit none

integer, parameter :: Nx = 201, Ny = 201
real(kind = 8), parameter :: a = 0., b = 1., c = 0., d = 2., tol = 1e-4, f = 1., pi = acos(-1.)
real(kind = 8), dimension(Nx-1,Ny-1) :: u, unew, I, exact, M, Q
real(kind = 8), dimension(Nx-1) :: x
real(kind = 8), dimension(Ny-1) :: y
real(kind = 8), dimension(3) :: Ax, Ay
real(kind = 8) :: hx, hy, delt, error, Dxp, Dxm, Dyp, Dym, Dx, Dy, L, error1, minm, uavg, xavg, yavg
integer :: j,k,Ix,Iy,p, iterations = 0

hx = (b-a)/Nx
hy = (d-c)/Ny

do j = 1,Nx-1
   x(j) = a + (j)*hx
enddo

do j = 1,Ny-1
   y(j) = c + (j)*hy
enddo

open(unit = 1, file = "peaks.txt")
read(1, *) I

! Defining Parameters
minm = 0.01
do k = 1,Ny-1
  do j = 1,Nx-1
    if(I(j,k) > 0 .and. I(j,k) < minm) then
      minm = I(j,k)
    endif
  enddo
enddo

print 13, minm
13  format('Minm = ', se22.15)

do k = 1,Ny-1
   do j = 1,Nx-1
      u(j,k) = -0.5*log(minm*f**2)
      Q(j,k) = sqrt(f**2/(f**2+x(j)**2+y(k)**2))
      M(j,k) = I(j,k)*f**2/Q(j,k)
      unew(j,k) = 0
   enddo
enddo

delt = 0.1*min(hx,hy)
error = 100.

!  do p = 1,00
do while (error> tol)
   !j = 1, k = 1
   Dx = (u(2,1) - u(1,1))/(2*hx)
   Dy = (u(1,2) - u(1,1))/(2*hy)

   uavg = (u(2,1)+u(1,1)+u(1,2)+u(1,1))/4

   unew(1,1) = uavg - delt*(-exp(-2*uavg) + M(1,1)*sqrt(f**2*(Dx**2+Dy**2) + (x(1)*Dx+y(1)*Dy)**2 + Q(1,1)**2))

   ! j = 1, k = N+1
   Dx = (u(2,Ny-1)-u(1,Ny-1))/(2*hx)
   Dy = (u(1,Ny-1) - u(1,Ny-2))/(2*hy)

   uavg = (u(2,Ny-1)+u(1,Ny-1)+u(1,Ny-1)+u(1,Ny-2))/4

   unew(1,Ny-1) = uavg - delt*(-exp(-2*uavg) + M(1,Ny-1)*sqrt(f**2*(Dx**2+Dy**2) + (x(1)*Dx+y(Ny-1)*Dy)**2 + Q(1,Ny-1)**2))


   do k = 2,Ny-2
      do j = 2,Nx-2
        !  if(I(j,k) == 0.) then
        !     unew(j,k) = (u(j-1,k)+u(j+1,k)+u(j,k-1)+u(j,k+1))/4
        !  elseif(I(j,k) == 1.)then
        !     unew(j,k) = 0.
        !  else
         Dx = (u(j+1,k) - u(j-1,k))/(2*hx)
         Dy = (u(j,k+1) - u(j,k-1))/(2*hy)

         uavg = (u(j+1,k)+u(j-1,k)+u(j,k+1)+u(j,k-1))/4

         unew(j,k) = uavg - delt*(-exp(-2*uavg) + M(j,k)*sqrt(f**2*(Dx**2+Dy**2) + (x(j)*Dx+y(k)*Dy)**2 + Q(j,k)**2))

        ! endif
      enddo
   enddo

   ! j = N+1, k = 1
   Dx = (u(nx-1,1) - u(Ny-2,1))/(2*hx)
   Dy = (u(Nx-1,2) - u(nx-1,1))/(2*hy)

   uavg = (u(nx-1,1)+u(Nx-2,1)+u(Nx-1,2)+u(nx-1,1))/4

   unew(nx-1,1) = uavg - delt*(-exp(-2*uavg) + M(nx-1,1)*sqrt(f**2*(Dx**2+Dy**2) + (x(nx-1)*Dx+y(1)*Dy)**2 + Q(nx-1,1)**2))

   ! j = N+1, k = N+1
   Dx = (u(nx-1,nx-1) - u(Nx-2,Ny-1))/(2*hx)
   Dy = (u(nx-1,nx-1) - u(Nx-1,Ny-2))/(2*hy)

   uavg = (u(nx-1,nx-1)+u(Nx-2,Nx-1)+u(nx-1,nx-1)+u(Nx-1,Nx-2))/4

   unew(nx-1,ny-1) = uavg - delt*(-exp(-2*uavg) + M(nx-1,ny-1)*sqrt(f**2*(Dx**2+Dy**2) + (x(nx-1)*Dx+y(ny-1)*Dy)**2 + Q(nx-1,ny-1)**2))

   !j = 2:N, k = 1
   do j=2,Nx-2
      Dx = (u(j+1,1) - u(j-1,1))/(2*hx)
      Dy = (u(j,2) - u(j,1))/(2*hy)

      uavg = (u(j+1,1)+u(j-1,1)+u(j,2)+u(j,1))/4

      unew(j,1) = uavg - delt*(-exp(-2*uavg) + M(j,1)*sqrt(f**2*(Dx**2+Dy**2) + (x(j)*Dx+y(1)*Dy)**2 + Q(j,1)**2))
   enddo

   !j = 2:N,k=N+1
   do j=2,Nx-2
      Dx = (u(j+1,Ny-1) - u(j-1,Ny-1))/(2*hx)
      Dy = (u(j,Ny-1) - u(j,Ny-2))/(2*hy)

      uavg = (u(j+1,Ny-1)+u(j-1,Ny-1)+u(j,Ny-1)+u(j,Ny-2))/4

      unew(j,ny-1) = uavg - delt*(-exp(-2*uavg) + M(j,ny-1)*sqrt(f**2*(Dx**2+Dy**2) + (x(j)*Dx+y(ny-1)*Dy)**2 + Q(j,ny-1)**2))
   enddo

   !k = 2:N, j = 1
   do k = 2,Ny-2
      Dx = (u(2,k) - u(1,k))/(2*hx)
      Dy = (u(1,k+1) - u(1,k-1))/(2*hy)

      uavg = (u(2,k)+u(1,k)+u(1,k+1)+u(1,k-1))/4

      unew(1,k) = uavg - delt*(-exp(-2*uavg) + M(1,k)*sqrt(f**2*(Dx**2+Dy**2) + (x(1)*Dx+y(k)*Dy)**2 + Q(1,k)**2))
   enddo

   !k=2:N, j = N+1
   do k=2,Nx-2
      Dx = (u(Nx-1,k) - u(Nx-2,k))/(2*hx)
      Dy = (u(Nx-1,k+1) - u(Nx-1,k-1))/(2*hy)

      uavg = (u(Nx-1,k)+u(Nx-2,k)+u(Nx-1,k+1)+u(Nx-1,k-1))/4

      unew(nx-1,k) = uavg - delt*(-exp(-2*uavg) + M(nx-1,k)*sqrt(f**2*(Dx**2+Dy**2) + (x(nx-1)*Dx+y(k)*Dy)**2 + Q(nx-1,k)**2))
   enddo

   error = maxval(abs(u-unew))
!   print 11, error
!11 format("error=", se22.15)

 u = unew

 iterations = iterations + 1
 if(mod(iterations,1) == 0) then
   print 12, iterations
   12 format("Total Iterations = ", i10)
   print 11, error
   11 format("error=", se22.15)
 endif
enddo

open(unit = 1,file = "op_neumannf.txt")
do j=1,Nx-1
   do k=1,Ny-1
      write(1,*) x(j), y(k),-exp(unew(j,k))
   enddo
enddo



end program persp1
