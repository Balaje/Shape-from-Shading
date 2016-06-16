%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Program to implement the orthographic SfS model            
%       proposed by Rouy and Tourin.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

format long

Nx = 187;
Ny = 187;

A = imread('vase1.png');
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
tol = 1e-5;

iterations = 0;
while error > tol

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
    error = max(max(abs(unew-u)));
    %mesh(x,y,unew);
    %axis([a,b,c,d,0,1]);
    %pause(0.01);
    u = unew;
    iterations = iterations + 1;
end

iterations
mesh(x,y,unew)
