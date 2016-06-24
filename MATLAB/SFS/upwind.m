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
N = 10;

global upx

for p=1:4
    h(p) = (b-a)/N;
    
    x = (a+h(p):h(p):b-h(p))';
    
    u = zeros(N-1,1);
    for i=1:N-1
        u(i) = 0;
    end
    unew = zeros(N-1,1);
    
    delt = h(p);
    
    error = 100;
    tol = 1e-5;
    
    while error > tol
        % i = 1
        Dm = (u(1))/h(p);
        Dp = (u(2)-u(1))/h(p);
        unew(1) = u(1) - delt*(godunov(Dm,Dp,a,x(2)) - sqrt(1-J(x(1))^2));
        
        % i = 2:N-2
        for i = 2:N-2
            if(J(x(i)) == 1)
                %             if(x(i)==0.25)
                %                 unew(i) = 1;
                %             elseif(x(i) == 0.75)
                %                 unew(i) = -1;
                %             end
                unew(i) = 1;
            else
                Dm = (u(i)-u(i-1))/h(p);
                Dp = (u(i+1)-u(i))/h(p);
                unew(i) = u(i) - delt*(godunov(Dm,Dp,x(i-1),x(i+1)) - sqrt(1-J(x(i))^2));
            end
        end
        
        % i = N-1
        Dm = (u(N-1)-u(N-2))/h(p);
        Dp = (-u(N-1))/h(p);
        unew(N-1) = u(N-1) - delt*(godunov(Dm,Dp,x(N-2),b) - sqrt(1-J(x(N-1))^2));
        
        error = max(abs(u-unew));
        u = unew;
        
    end
    
    exact = zeros(N-1,1);
    for i=1:N-1
        exact(i) = sin(pi*x(i));
        %exact(i) = x(i)*(1-x(i));
    end
    
    figure(p);
    plot(x,unew,x,exact,'*');
    grid on
    
    error1(p) = max(abs(unew(i)-exact(i)));
    
    N = N*2;
end

for j=1:p-1
    order(j) = log(error1(j)/error1(j+1))/log(2.);
end

%order'
error1'


