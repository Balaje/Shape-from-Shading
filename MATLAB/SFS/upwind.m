%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Solving the 1D SfS orthographic model
%       |u_x| = sqrt(1/I^2-1)
%
%   using an upwind like scheme.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

format long

a = 0; b = 1;
N = 50;

global upx

h = (b-a)/N;

x = (a+h:h:b-h)';

u = zeros(N-1,1);
unew = zeros(N-1,1);

delt = h;

error = 100;
tol = 1e-5;

while error > tol
    % i = 1
    Dm = (u(1))/h;
    Dp = (u(2)-u(1))/h;
    unew(1) = u(1) - delt*(godunov(Dm,Dp,x(1),x(2)) - sqrt(1-J(x(1))^2));
    
    % i = 2:N-2
    for i = 2:N-2
        if(J(x(i)) == 1)
             unew(i) = 1/4;
         else
            Dm = (u(i)-u(i-1))/h;
            Dp = (u(i+1)-u(i))/h;
            unew(i) = u(i) - delt*(godunov(Dm,Dp,x(i),x(i+1)) - sqrt(1-J(x(i))^2));
        end
    end
    
    % i = N-1
    Dm = (u(N-1)-u(N-2))/h;
    Dp = (-u(N-1))/h;
    unew(N-1) = u(N-1) - delt*(godunov(Dm,Dp,x(N-1),b) - sqrt(1-J(x(N-1))^2));
    
    error = max(abs(u-unew));
    u = unew;
    
end

plot(x,unew);


