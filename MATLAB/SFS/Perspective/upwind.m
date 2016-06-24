%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Upwind scheme to solve the perspective sfs
%       model.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

a = 0; b = 1;

N = 10;
global upx

for p=1:5
    h(p) = (b-a)/N;
    
    x = (a+h(p):h(p):b-h(p))';
    
    u = zeros(N-1,1);
    unew = zeros(N-1,1);
    
    error = 100;
    tol = 1e-12;
    f = 1.;
    
    delt = h(p);
    
    while error > tol
        % i = 1
        Dm = (u(1))/h(p);
        Dp = (u(2)-u(1))/h(p);
        
        D = [Dm,Dp,0];
        [Dx,idx] = godunov(Dm,Dp,a,x(2));
        
        unew(1) = u(1) - delt*(sqrt(f^2*Dx + (J(x(1))*D(idx)*x(1))^2 + J(x(1))^2*f^2/(f^2+x(1)^2))...
            - f/sqrt(x(1)^2+f^2));
        
        % i = 2:N-2
        for i = 2:N-2
            if(J(x(i)) == 1)
                unew(i) = 1;
            else
                Dm = (u(i)-u(i-1))/h(p);
                Dp = (u(i+1)-u(i))/h(p);
                
                D = [Dm,Dp,0];
                
                [Dx,idx] = godunov(Dm,Dp,x(i-1),x(i+1));
                
                unew(i) = u(i) - delt*(sqrt(f^2*Dx + (J(x(i))*D(idx)*x(i))^2 + J(x(i))^2*f^2/(f^2+x(i)^2))...
                    - f/sqrt(x(i)^2+f^2));
            end
        end
        
        % i = N-1
        Dm = (u(N-1)-u(N-2))/h(p);
        Dp = (-u(N-1))/h(p);
        
        D = [Dm,Dp,0];
        
        [Dx,idx] = godunov(Dm,Dp,x(N-2),b);
        
        unew(N-1) = u(N-1) - delt*(sqrt(f^2*Dx + (J(x(N-1))*D(idx)*x(N-1))^2 + J(x(N-1))^2*f^2/(f^2+x(N-1)^2))...
            - f/sqrt(x(N-1)^2+f^2));
        
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

order'