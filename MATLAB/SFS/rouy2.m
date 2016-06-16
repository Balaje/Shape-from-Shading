%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Program to implement the orthographic SfS model            
%       proposed by Rouy and Tourin.
%   Neumann Boundary Condition - Full Neumann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

format long

Nx = 187;
Ny = 187;

A = imread('mozart1.png');
gray = mat2gray(imresize(A,[Nx+1,Ny+1]));

a = 0; b = 1;
c = 0; d = 1;

hx = (b-a)/Nx;
hy = (b-a)/Ny;

x = a:hx:b;
y = c:hx:d;

% Defining the index set
Q = zeros(Nx+1,Ny+1);
u = zeros(Nx+1,Ny+1);
unew = zeros(Nx+1,Ny+1);

for i=1:Nx+1
    for j=1:Ny+1
        u(i,j) = -0.5*log(gray(i,j,1));
    end
end

delt = hx*hy/sqrt(hx^2+hy^2);

error = 100;
tol = 1e-4;

iterations = 0;
while error > tol
    % i = 1, j = 1
    Dxp = (u(1,1)-u(2,1))/hx;
    Dxm = 0;
    Dyp = (u(1,1)-u(1,2))/hy;
    Dym = 0;
    
    Dx = max([0,Dxp,Dxm]);
    Dy = max([0,Dyp,Dym]);
    
    so = sqrt(1/(gray(1,1,1))^2 - 1 );
    H = sqrt(Dx^2+Dy^2);
    
    unew(1,1) = u(1,1) - delt*(H-so);
    
    % j = 1, i = 2:Ny
    for i=2:Ny
        Dxp = (u(i,1)-u(i+1,1))/hx;
        Dxm = (u(i,1)-u(i-1,1))/hx;
        Dyp = (u(i,1)-u(i,2))/hy;
        Dym = 0;
        
        Dx = max([0,Dxp,Dxm]);
        Dy = max([0,Dyp,Dym]);
        
        so = sqrt(1/(gray(i,1,1))^2 - 1 );
        H = sqrt(Dx^2+Dy^2);
        
        unew(i,1) = u(i,1) - delt*(H-so);
       
    end
    
    % i = 1, j = Ny+1
    Dxp = (u(1,Ny+1)-u(2,Ny+1))/hx;
    Dxm = 0;
    Dyp = 0;
    Dym = (u(1,Ny+1)-u(1,Ny))/hy;
    
    Dx = max([0,Dxp,Dxm]);
    Dy = max([0,Dyp,Dym]);
            
    so = sqrt(1/(gray(1,Ny+1,1))^2 - 1 );
    H = sqrt(Dx^2+Dy^2);
    
    unew(1,Ny+1) = u(1,Ny+1) - delt*(H-so);    
    
    
    for i = 2:Nx
        for j = 2:Ny
                Dxp = (u(i,j)-u(i+1,j))/hx;
                Dxm = (u(i,j)-u(i-1,j))/hx;
                Dyp = (u(i,j)-u(i,j+1))/hy;
                Dym = (u(i,j)-u(i,j-1))/hy;
            
                Dx = max([0,Dxp,Dxm]);
                Dy = max([0,Dyp,Dym]);
            
                so = sqrt(1/(gray(i,j,1))^2 - 1 );
                H = sqrt(Dx^2+Dy^2);
            
                unew(i,j) = u(i,j) - delt*(H-so);
        end
    end
    
    % i = Nx+1, j = 1
    Dxp = 0;
    Dxm = (u(Nx+1,1)-u(Nx,1))/hx;
    Dyp = (u(Nx+1,1)-u(Nx+1,2))/hy;
    Dym = 0;
    
    Dx = max([0,Dxp,Dxm]);
    Dy = max([0,Dyp,Dym]);           
    
    so = sqrt(1/(gray(Nx+1,1,1))^2 - 1 );
    H = sqrt(Dx^2+Dy^2);            
    
    unew(Nx+1,1) = u(Nx+1,1) - delt*(H-so);    
                
    % j = Ny+1, i = 2:Nx
    for i=2:Nx
        Dxp = (u(i,Ny+1)-u(i+1,Ny+1))/hx;
        Dxm = (u(i,Ny+1)-u(i-1,Ny+1))/hx;
        Dyp = 0;
        Dym = (u(i+1,Ny+1)-u(i,Ny))/hy;
        
        Dx = max([0,Dxp,Dxm]);
        Dy = max([0,Dyp,Dym]);
        
        so = sqrt(1/(gray(i,Ny+1,1))^2 - 1 );
        H = sqrt(Dx^2+Dy^2);          
                
        unew(i,Ny+1) = u(i,Ny+1) - delt*(H-so);  
    end
    % i = Nx+1, j = Ny+1
     Dxp = 0;
     Dxm = (u(Nx+1,Ny+1)-u(Nx,Ny+1))/hx;
     Dyp = 0;
     Dym = (u(Nx+1,Ny+1)-u(Nx+1,Ny))/hy;
     
     Dx = max([0,Dxp,Dxm]);
     Dy = max([0,Dyp,Dym]);
     
     so = sqrt(1/(gray(Nx+1,Ny+1,1))^2 - 1 );
     H = sqrt(Dx^2+Dy^2);
            
     unew(Nx+1,Ny+1) = u(Nx+1,Ny+1) - delt*(H-so);
    
    % j = 2:Nx, i = 1
    for j = 2:Ny
        Dxp = (u(1,j)-u(2,j))/hx;
        Dxm = 0;
        Dyp = (u(1,j)-u(1,j+1))/hy;
        Dym = (u(1,j)-u(1,j-1))/hy;
            
        Dx = max([0,Dxp,Dxm]);
        Dy = max([0,Dyp,Dym]);
            
        so = sqrt(1/(gray(1,j,1))^2 - 1 );
        H = sqrt(Dx^2+Dy^2);
            
        unew(1,j) = u(1,j) - delt*(H-so);         
    end
    % j = 2:Ny, i = Nx+1
    for j=2:Ny
        Dxp = 0;
        Dxm = (u(Nx+1,j)-u(Nx,j))/hx;
        Dyp = (u(Nx+1,j)-u(Nx+1,j+1))/hy;
        Dym = (u(Nx+1,j)-u(Nx+1,j-1))/hy;
            
        Dx = max([0,Dxp,Dxm]);
        Dy = max([0,Dyp,Dym]);
            
        so = sqrt(1/(gray(Nx+1,j,1))^2 - 1 );
        H = sqrt(Dx^2+Dy^2);
            
        unew(Nx+1,j) = u(Nx+1,j) - delt*(H-so);  
    end
    
    error = max(max(abs(unew-u)));
    
    u = unew;
    iterations = iterations + 1;
end

iterations
mesh(x,y,unew)
