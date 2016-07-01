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

a = 0; b = 1;

h = (b-a)/N;

x = a+h:h:b-h;

J = zeros(N-1,1);
Q = zeros(N-1,1);
f = 1;

for i=1:N-1
    Q(i) = sqrt(f^2/(x(i)^2+f^2));
    J(i) = I(i)*f^2/Q(i);
end
u = zeros(N-1,1);
unew = zeros(N-1,1);
error = 100;
tol = 1e-5;

delt = 0.5*h;
while error > tol
    % i = 1
    Dxp = (u(2)-u(1))/h;
    Dxm = (u(1)-0)/h;
    [D,ind] = max([0,-Dxp,Dxm]);
    A = [0,Dxp,Dxm];
    
    H = -exp(-2*u(1)) + J(1)*sqrt( f^2*D^2 + (A(ind)*x(1))^2 + Q(1)^2);
    unew(1) = u(1) - delt*H;
    
    for j=2:N-2
        Dxp = (u(j+1)-u(j))/h;
        Dxm = (u(j)-u(j-1))/h;
        [D,ind] = max([0,-Dxp,Dxm]);
        A = [0,Dxp,Dxm];
        
        H = -exp(-2*u(j)) + J(j)*sqrt( f^2*D^2 + (A(ind)*x(j))^2 + Q(j)^2);
        unew(j) = u(j) - delt*H;
    end
    
    % i = N-1
    Dxp = (-u(N-1))/h;
    Dxm = (u(N-1)-u(N-2))/h;
    [D,ind] = max([0,-Dxp,Dxm]);
    A = [0,Dxp,Dxm];
    
    H = -exp(-2*u(N-1)) + J(N-1)*sqrt( f^2*D^2 + (A(ind)*x(N-1))^2 + Q(N-1)^2);
    unew(N-1) = u(N-1) - delt*H;
    
    
    error = max(abs(u-unew));
    
%     plot(x,unew);
%     pause(0.01);
    
    u = unew;
end

V = unew;

