%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Program to implement the orthographic SfS model            
%       proposed by Rouy and Tourin.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

format long

Nx = 50;
Ny = 50;

%A = imread('vase1.png');
%gray = mat2gray(imresize(A,[Nx+1,Ny+1]));

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
        if(I(x(i),y(j))==1)
            u(i,j) = 1;
        end
        %if(x(i)==0.5 && y(j) ==0.5)
        %    u(i,j) = 2;
        %end
    end
end

delt = hx*hy/sqrt(hx^2+hy^2);

error = 100;
tol = 1e-9;

iterations = 0;
while error > tol

    for i = 2:Nx
        for j = 2:Ny
            %if(x(i)==0.5 && y(j) ==0.5)
            %    unew(i,j) = 2;
            %end
            if(I(x(i),y(j))==1)
                unew(i,j) = 1;
            else
                Dxp = (u(i,j)-u(i+1,j))/hx;
                Dxm = (u(i,j)-u(i-1,j))/hx;
                Dyp = (u(i,j)-u(i,j+1))/hy;
                Dym = (u(i,j)-u(i,j-1))/hy;
            
                [Dx,p] = max([Dxm,Dxp,0]);
                [Dy,q] = max([Dym,Dyp,0]);
                %[Dx,p] = max([max(0,Dxm),min(0,Dxp)]);
                %[Dy,q] = max([max(0,Dym),min(0,Dyp)]);
            
                if(p==1 && q == 1)
                    so = sqrt(1/(I(x(i),y(j)))^2 - 1 );
                elseif(p==1 && q==2)
                    so = sqrt(1/(I(x(i),y(j+1)))^2 - 1 );
                elseif(p==2 && q==1)
                    so = sqrt(1/(I(x(i+1),y(j)))^2 - 1 );
                elseif(p==2 && q==2)
                    so = sqrt(1/(I(x(i+1),y(j+1)))^2 - 1 );
                else
                    so = sqrt(1/(I(x(i),y(j)))^2 - 1 );
                end
                
                H = sqrt(Dx^2+Dy^2);
                unew(i,j) = u(i,j) - delt*(H-so);
%                 surf(x,y,unew);
%                 pause(0.01);
                
            end
        end
    end
    error = max(max(abs(unew-u)))
    %mesh(x,y,unew);
    %axis([a,b,c,d,0,1]);
    %pause(0.01);
    u = unew;
    iterations = iterations + 1;
end

iterations
surf(x,y,unew)
