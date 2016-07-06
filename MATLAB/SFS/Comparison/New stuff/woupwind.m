function V = woupwind()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to check the effect of upwinding on   %
%           Discontinous Intensity function     %
%   Run : Without upwinding                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
close all
clear

global unew exactsol a b h

a = 0; b = 1;
N = 50;
h = (b-a)/N;

x = a+h:h:b-h;

u = zeros(N-1,1);
unew = zeros(N-1,1);
exactsol = zeros(N-1,1);

error = 100;
tol = 1e-10;

delt = h;
iterations = 0;

eps = 0;

while error > tol
    % i = 1
    Dxp = (u(2)-u(1))/h;
    Dxm = (u(1)-0)/h;
    D = godunov(Dxm,Dxp,(a+x(1))/2,(a+x(1))/2);
    
    unew(1) = u(1) - delt*(-exp(-2*u(1))+D-(eps*(u(2)-2*u(1))/h^2));
    
    for j=2:N-2
        Dxp = (u(j+1)-u(j))/h;
        Dxm = (u(j)-u(j-1))/h;
        D = godunov(Dxm,Dxp,(x(j-1)+x(j))/2,(x(j-1)+x(j))/2);
        unew(j) = u(j) - delt*(-exp(-2*u(j))+D - eps*(u(j-1)-2*u(j)+u(j+1))/h^2);
    end
    
    % i = N-1
    Dxp = (0-u(N-1))/h;
    Dxm = (u(N-1)-u(N-2))/h;
    
    D = godunov(Dxm,Dxp,(x(N-2)+x(N-1))/2,(x(N-2)+x(N-1))/2);
    unew(N-1) = u(N-1) - delt*(-exp(-2*u(N-1))+D - eps*(u(N-2)-2*u(N-1)+0.15)/h^2);
    
    error = max(abs(u-unew));
    
%         figure(2)
%         plot(x,unew,x,exactsol);
%         axis([a,b,0,0.3])
%         pause(0.001);
    
    iterations = iterations + 1;
    u = unew;
end

V = unew;