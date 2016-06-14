%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program to solve the 2d Eikonal Equation using conservation laws                                                    
%         u_t + |\Grad u| = 1                                           
%                                                   
%                                                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

format long

N = 100;

a = -1; b = 1; 
c = -1; d = 1;

h = (b-a)/N;

x = a:h:b;
y = c:h:b;

u = zeros(N+1);
unew = zeros(N+1);

delt = h;

error = 100;
eps = 1e-5;

while error > eps
    for i = 2:N
        for j = 2:N
            Dpux = (u(i+1,j)-u(i,j))/h;
            Dmux = (u(i,j)-u(i-1,j))/h;
            Dpuy = (u(i,j+1)-u(i,j))/h;
            Dmuy = (u(i,j)-u(i,j-1))/h;
            
            Dx = fluxfunction(Dmux,Dpux,'m',delt/h);
            Dy = fluxfunction(Dmuy,Dpuy,'m',delt/h);
            
            D = sqrt(Dx^2+Dy^2);
            
            unew(i,j) = u(i,j) - delt*(D-1);
        end
    end
    error = max(max(abs(u-unew)));
    
    u = unew;
    
    surf(x,y,unew);
    axis([a,b,c,d,0,1])
    pause(0.01);
end


