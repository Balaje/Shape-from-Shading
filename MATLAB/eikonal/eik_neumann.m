%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Program to implement the upwind like scheme          % 
%   for 1D Eikonal Equation (Neumann Problem)          %
%       |u'| = 1 in \Omega                             % 
%        u(0) = u'(1) = 0                              % 
%                                                      % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear
close all

format long

N = 50;
a = -1;
b = 1;

h = (b-a)/N;

x = (a:h:b)';

delt = h;

% Tolerance for convergence
eps = 10^-4;
error = 100;

u = zeros(N+1,1);

unew = u;
iter = 0;
while error > eps
    for i = 2:N
        H = max(abs(max([(u(i)-u(i-1))/h,0])),abs(min([(u(i+1)-u(i))/h,0])));
        if(x(i) < 0)
            unew(i) = u(i) - delt*(H-1);
        else
            unew(i) = u(i) - delt*(H-3);
        end
    end
    error = max(abs(unew-u));
    u = unew;
    iter = iter+1;
    plot(x,unew);
    pause(0.1);
end
plot(x,unew)
iter