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

A = imread('mozart.png');
gray = mat2gray(imresize(A,[Nx+1,Ny+1]));

a = 0; b = 1;
c = 0; d = 1;

hx = (b-a)/Nx;
hy = (b-a)/Ny;

x = a:hx:b;
y = c:hx:d;

% Defining the index set
Q = zeros(Nx+1,Ny+1);
u = ones(Nx+1,Ny+1);
unew = ones(Nx+1,Ny+1);

for i=1:Nx+1
    for j=1:Ny+1
        if(abs(gray(i,j,1)-1) < 1e-12 )
            Q(i,j) = 0;
            u(i,j) = 1;
        else
            Q(i,j) = 1;
        end
    end
end

delt = hx*hy/sqrt(hx^2+hy^2);

error = 100;
tol = 1e-3;

iterations = 0;
while error > tol
    for i = 2:Nx
        for j = 2:Ny
            if(Q(i,j) == 0)
%                 if(x(i) == 0.5 && y(j) == 0.5)
%                     unew(i,j) = 2;
%                 else
%                     unew(i,j) = 1;
%                 end
                unew(i,j) = 1;
            else
                Dxp = (u(i,j)-u(i+1,j))/hx;
                Dxm = (u(i,j)-u(i-1,j))/hx;
                Dyp = (u(i,j)-u(i,j+1))/hy;
                Dym = (u(i,j)-u(i,j-1))/hy;
            
                Dx = max([0,Dxp,Dxm]);
                Dy = max([0,Dyp,Dym]);
            
                H = gray(i,j,1)*sqrt((Dx^2+Dy^2)+(Dx*x(i)+Dy*y(j))^2 + 1/(1+x(i)^2+y(j)^2)) - 1/sqrt(1+x(i)^2+y(j)^2);
            
                unew(i,j) = u(i,j) - delt*(H);
             end
        end
    end
    error = max(max(abs(unew-u)));
    mesh(x,y,exp(unew));
    %axis([a,b,c,d,0,1]);
    pause(0.01);
    u = unew;
    iterations = iterations + 1;
end

iterations
mesh(x,y,exp(unew))
axis([a,b,c,d,0,1]);
