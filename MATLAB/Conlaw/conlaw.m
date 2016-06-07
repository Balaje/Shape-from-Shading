%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program to solve the Burger's equation using Godunov Scheme  
%       u_t + f(u)_x = 0, in (-1,2)x(0,T]
%           u(x,0) = {0, x<0                                        
%                     1, x>0  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
format long

N = 100;
a = -5; b = 5;

solver = 'g';

eps = 1;

h = (b-a)/N;

%CFL Condition
delt = eps*h;

tf = 2;
t0 = 0;

M = fix((tf-t0)/delt);

x = zeros(N-1,1);
u = zeros(N-1,1);
unew = zeros(N-1,1);

for i=1:N-1
    x(i) = a + i*h;
    if(x(i) < 0)
        u(i) = -1;
    elseif(x(i) >= 0)
        u(i) = 1;
    end
end

for k = 1:M
    unew(1) = u(1) - eps*(fluxfunction(u(1),u(2),solver,eps) - ...
                            fluxfunction(-1,u(1),solver,eps));
    for j = 2:N-2
        unew(j) = u(j) - eps*(fluxfunction(u(j),u(j+1),solver,eps) - ...
                            fluxfunction(u(j-1),u(j),solver,eps));
    end
    unew(N-1) = u(N-1) - eps*(fluxfunction(u(N-1),1,solver,eps) - ...
                            fluxfunction(u(N-2),u(N-1),solver,eps));
                        
    
    plot(x,unew);
    hold on
    pause(0.01);
    
    u = unew;           
end

figure(2)
plot(x,unew,'o-');
legend('Approximate solution');
title(strcat(solver,' Scheme'));
xlabel('x');
ylabel(strcat('u(x,t = ',num2str(tf),')'));
grid on



