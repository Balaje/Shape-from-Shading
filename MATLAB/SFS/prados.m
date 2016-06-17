%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Program to implement the Perspective SfS model            
%       proposed by Prados and Faugeras.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

format long

Nx = 187;
Ny = 187;

A = imread('abe.png');
gray = mat2gray(imresize(A,[Nx+1,Ny+1]));

a = 0; b = 1;
c = 0; d = 1;

hx = (b-a)/Nx;
hy = (b-a)/Ny;

x = a:hx:b;
y = c:hx:d;

u = zeros(Nx+1,Ny+1);

unew = zeros(Nx+1,Ny+1);

delt = hx*hy/sqrt(hx^2+hy^2);

error = 100;
tol = 1e-12;

f = 1;

iterations = 0;
while error > tol
    % i = 1, j = 1
    Dxp = (u(1,1)-u(2,1))/hx;
    Dxm = 0;
    Dyp = (u(1,1)-u(1,2))/hy;
    Dym = 0;
    
    Dx = max([0,Dxp,Dxm]);
    Dy = max([0,Dyp,Dym]);
    
    H = gray(1,1,1)*sqrt(f^2*(Dx^2+Dy^2)+(Dx*x(1)+Dy*y(1))^2 + ...
        f^2/(f^2+x(1)^2+y(1)^2)) - f/sqrt(f^2+x(1)^2+y(1)^2);
    
    unew(1,1) = u(1,1) - delt*(H);    
    
    % i = 1, j = Ny+1
    Dxp = (u(1,Ny+1)-u(2,Ny+1))/hx;
    Dxm = 0;
    Dyp = 0;
    Dym = (u(1,Ny+1)-u(1,Ny))/hy;
    
    Dx = max([0,Dxp,Dxm]);
    Dy = max([0,Dyp,Dym]);
            
    H = gray(1,Ny+1,1)*sqrt(f^2*(Dx^2+Dy^2)+(Dx*x(1)+Dy*y(Ny+1))^2 + ...
        f^2/(f^2+x(1)^2+y(Ny+1)^2)) - f/sqrt(f^2+x(1)^2+y(Ny+1)^2);
    
    unew(1,Ny+1) = u(1,Ny+1) - delt*(H);
    for i = 2:Nx
        for j = 2:Ny
            if(abs(gray(i,j,1)) < 0.1)
                unew(i,j) = 0;
            else
                Dxp = (u(i,j)-u(i+1,j))/hx;
                Dxm = (u(i,j)-u(i-1,j))/hx;
                Dyp = (u(i,j)-u(i,j+1))/hy;
                Dym = (u(i,j)-u(i,j-1))/hy;
            
                Dx = max([0,Dxp,Dxm]);
                Dy = max([0,Dyp,Dym]);
            
                H = gray(i,j,1)*sqrt(f^2*(Dx^2+Dy^2)+(Dx*x(i)+Dy*y(j))^2 + ...
                    f^2/(f^2+x(i)^2+y(j)^2)) - f/sqrt(f^2+x(i)^2+y(j)^2);
            
                unew(i,j) = u(i,j) - delt*(H);
            end
        end
    end
    
     % i = Nx+1, j = 1
    Dxp = 0;
    Dxm = (u(Nx+1,1)-u(Nx,1))/hx;
    Dyp = (u(Nx+1,1)-u(Nx+1,2))/hy;
    Dym = 0;
    
    Dx = max([0,Dxp,Dxm]);
    Dy = max([0,Dyp,Dym]);           
    
    H = gray(Nx+1,1,1)*sqrt(f^2*(Dx^2+Dy^2)+(Dx*x(Nx+1)+Dy*y(1))^2 ...
        + f^2/(f^2+x(Nx+1)^2+y(1)^2)) - f/sqrt(f^2+x(Nx+1)^2+y(1)^2);
    
    unew(Nx+1,1) = u(Nx+1,1) - delt*(H);    
                
    % j = Ny+1, i = 2:Nx
    
    % i = Nx+1, j = Ny+1
     Dxp = 0;
     Dxm = (u(Nx+1,Ny+1)-u(Nx,Ny+1))/hx;
     Dyp = 0;
     Dym = (u(Nx+1,Ny+1)-u(Nx+1,Ny))/hy;
     
     Dx = max([0,Dxp,Dxm]);
     Dy = max([0,Dyp,Dym]);
     
    H = gray(Nx+1,Ny+1,1)*sqrt(f^2*(Dx^2+Dy^2)+(Dx*x(Nx+1)+Dy*y(Ny+1))^2 +...
        f^2/(f^2+x(Nx+1)^2+y(Ny+1)^2)) - f/sqrt(f^2+x(Nx+1)^2+y(Ny+1)^2);
            
     unew(Nx+1,Ny+1) = u(Nx+1,Ny+1) - delt*(H);
    
    % j = 2:Nx, i = 1
    for j = 2:Ny
        Dxp = (u(1,j)-u(2,j))/hx;
        Dxm = 0;
        Dyp = (u(1,j)-u(1,j+1))/hy;
        Dym = (u(1,j)-u(1,j-1))/hy;
            
        Dx = max([0,Dxp,Dxm]);
        Dy = max([0,Dyp,Dym]);
            
        H = gray(1,j,1)*sqrt(f^2*(Dx^2+Dy^2)+(Dx*x(1)+Dy*y(j))^2 + ...
            f^2/(f^2+x(1)^2+y(j)^2)) - f/sqrt(f^2+x(1)^2+y(j)^2);
            
        unew(1,j) = u(1,j) - delt*(H);         
    end
    % j = 2:Ny, i = Nx+1
    for j=2:Ny
        Dxp = 0;
        Dxm = (u(Nx+1,j)-u(Nx,j))/hx;
        Dyp = (u(Nx+1,j)-u(Nx+1,j+1))/hy;
        Dym = (u(Nx+1,j)-u(Nx+1,j-1))/hy;
            
        Dx = max([0,Dxp,Dxm]);
        Dy = max([0,Dyp,Dym]);
           
        H = gray(Nx+1,j,1)*sqrt(f^2*(Dx^2+Dy^2)+(Dx*x(Nx+1)+Dy*y(j))^2 + ...
            f^2/(f^2+x(Nx+1)^2+y(j)^2)) - f/sqrt(f^2+x(Nx+1)^2+y(j)^2);
            
        unew(Nx+1,j) = u(Nx+1,j) - delt*(H);  
    end    
    
    error = max(max(abs(unew-u)));
    %mesh(x,y,exp(unew));
    %axis([a,b,c,d,0,1]);
    %pause(0.01);
    u = unew;
    iterations = iterations + 1;
end

iterations
mesh(x,y,exp(unew))