%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to implement the upwind scheme to solve        %   
%     2D Eikonal Equation with 0 DBC on a square         %   
%       |\Grad u| = 1 in (0,1) X (0,1)                   %   
%       u = 0 on Boundary                                %
%                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
format long

Nx = 10;
Ny = 10;

a =-1; b = 1;
c = -1; d = 1;

hx = (b-a)/Nx;
hy = (b-a)/Ny;

x = a:hx:b;
y = c:hy:d;

delt = min(hx,hy); % Monotonicity condition

% Tolerance for steady state
eps = 10^-5;
error = 100;

U = zeros(Nx+1,Ny+1);
% for i=1:Nx+1
%     for j=1:Ny+1
%         U(i,j) = 
%     end
% end
Unew = zeros(Nx+1,Ny+1);

iterations = 0;
tic
while error > eps
    for i=2:Nx
        for j=2:Ny
            Dxp = (U(i,j)-U(i+1,j))/hx;
            Dxm = (U(i,j)-U(i-1,j))/hx;
            Dyp = (U(i,j)-U(i,j+1))/hy;
            Dym = (U(i,j)-U(i,j-1))/hy;
            
            Dx = max([0,Dxp,Dxm]);
            Dy = max([0,Dyp,Dym]);
            
            H = sqrt(Dx^2+Dy^2)-1;
            
            Unew(i,j) = U(i,j) - delt*H;
        end
    end
    error = max(max(abs(Unew-U)));
    %error = 0;
    U = Unew;
    surf(x,y,Unew);
    axis([a,b,c,d,0,1]);
    pause(0.01);
    iterations = iterations+1;
end
toc
iterations
%surf(y,x,Unew);
