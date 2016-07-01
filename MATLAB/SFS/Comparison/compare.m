%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script to compare the solution plots of 
%       Upwinding vs No-upwinding
%
%   Functions : testdisc.m and testdisc_upw.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

N = 100;

a = 0; b = 1;
h = (b-a)/N;

x = a+h:h:b-h;

I = zeros(N-1,1);
exact = zeros(N-1,1);

for i=1:N-1
     if(x(i) <= 0.5)
        I(i) = exp(-2*x(i))/(1*sqrt(1+x(i)^2)*sqrt((1)^2+(x(i))^2+(1/(x(i)^2+1))));
        exact(i) = x(i);
    else
        I(i) = exp(-2*x(i)*(1-x(i)))/(1*sqrt(1+x(i)^2)*sqrt((1*(1-2*x(i)))^2+(x(i)*(1-2*x(i)))^2+(1/(x(i)^2+1))));
        exact(i) = x(i)*(1-x(i));
    end
end

sol_upw = testdisc_upw(N,I);
sol = testdisc(N,I);

figure(1)
plot(x,I,'*');
title('Intensity Map');

figure(2)
plot(x,sol_upw,x,sol,x,exact);
legend('Upwinded Solution','Not upwinded','Exact Solution')
title('Comparison of solution');