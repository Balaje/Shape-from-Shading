function V = testdisc(N,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program to check the effect of upwinding on
%
%           Discontinous Intensity function
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
global f

a = 0; b = 1;
h = (b-a)/N;

x = a+h:h:b-h;
f = 1;

J = zeros(N-1,1);
Q = zeros(N-1,1);


for i=1:N-1
    Q(i) = sqrt(f^2/(x(i)^2+f^2));
    J(i) = I(i)*f^2/Q(i);
end
Qb = sqrt(f^2/(b^2+f^2));
Ib = exp(-2*b*(1-b))/(1*sqrt(1+b^2)*sqrt((1*(1-2*b))^2+(b*(1-2*b))^2+(1/(b^2+1))));
Jb = Ib*f^2/Qb;

u = zeros(N-1,1);
unew = zeros(N-1,1);
error = 100;
tol = 1e-5;

delt = 0.5*h;
while error > tol
    % i = 1
    Dpx = (u(2)-u(1))/h;
    Dmx = (u(1)-0)/h;
    
    unew(1) = u(1) - delt*(godunov(Dmx,Dpx,J(1),J(1),x(1)) - exp(-2*u(1)));
    
    % i = 2:N-2
    for j = 2:N-2
        Dpx = (u(j+1)-u(j))/h;
        Dmx = (u(j)-u(j-1))/h;
        
        unew(j) = u(j) - delt*(godunov(Dmx,Dpx,J(j),J(j),x(j)) - exp(-2*u(j)));
    end
    
    % i = N-1
    Dpx = (0-u(N-1))/h;
    Dmx = (u(N-1)-u(N-2))/h;
    
    unew(N-1) = u(N-1) - delt*(godunov(Dmx,Dpx,J(N-1),J(N-1),x(N-1)) - exp(-2*u(N-1)));
    
    
    error = max(abs(u-unew));
    
    u = unew;
end

V = unew;
