%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Implicit Solver to solve the Eikonal Equation
%
%       u_t + |u_x| = 1
%           u(x,0) = 0
%               u(0,t) = u(1,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

format long

a = -1; b = 1;
N = 10;

h = (b-a)/N;
x = a+h:h:b-h;

unew = zeros(N-1,1);
u = zeros(N-1,1);
J = zeros(N-1,N-1);
F = zeros(N-1,1);

error = 100;
tol = 1e-15;

delt = 1;

tic
while error > tol
    % i = 1
    F(1) = unew(1) + delt/h*max([0,(unew(1)-unew(2)),(unew(1))])-u(1)-delt;
    J(1,1) = 1 + delt/h;
    
    %i=2:N-2
    for i=2:N-2
        F(i) = unew(i) + delt/h*max([0,(unew(i)-unew(i+1)),(unew(i)-unew(i-1))])-u(i)-delt;
        J(i,i) = 1+ delt/h;
    end
    
    F(N-1) = unew(N-1) + delt/h*max([0,(unew(N-1)),(unew(N-1)-unew(N-2))])-u(N-1)-delt;
    J(N-1,N-1) = 1+ delt/h;
    
    unew = u - J\F;
    
    error = max(abs(u-unew))
    
    u = unew;
    %plot(x,unew)
    %pause(0.001);
end
toc
plot(x,unew)