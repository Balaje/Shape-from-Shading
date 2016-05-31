%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Program to implement the upwind like scheme          % 
%   for 1D Eikonal Equation                            %
%       |u'| = 1 in \Omega                             % 
%        u(0) = u(1) = 0                               % 
%                                                      % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

format long

N = 20;
a = -1;
b = 1;

h = (b-a)/N;

x = (a:h:b)';

delt = h;

% Tolerance for convergence
eps = 10^-4;
error = 100;

u = (x.^2-1);
unew = u;
iter = 0;
while error > eps
    for i = 2:N
        Dp = (u(i+1)-u(i))/h;
        Dm = (u(i-1)-u(i))/h;
        
        D = min([0,Dp,Dm]);
        
        H = sqrt(D^2)-1;
        
        unew(i) = u(i) - delt*H;
    end
    error = max(abs(unew-u));
    u = unew;
    iter = iter+1;
    plot(x,unew);
    pause(0.1);
end
plot(x,unew)
iter