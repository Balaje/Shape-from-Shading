%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script to compare the solution plots of        %
%       Upwinding vs No-upwinding                  % 
%   Functions : testdisc.m and testdisc_upw.m      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

N = 100;
h = 1/N;

x = h:h:1-h;

I = zeros(N-1,1);
exact = zeros(N-1,1);

for i=1:N-1
    if(x(i) < 0.4 || x(i) > 0.6)
        I(i) = exp(-2*x(i)*(1-x(i)))/sqrt(1+((1-2*x(i))*(x(i)^2+1))^2);
        exact(i) = x(i)*(1-x(i));
    else
        I(i) = 1;
        exact(i) = 0;
        
    end
end

sol_upw = testdisc_upw(N,I);
sol = testdisc(N,I);

figure(1)
plot(x,I);
title('Intensity Map');

figure(2)
plot(x,sol_upw,x,sol,x,exact);
legend('Upwinded Solution','Not upwinded')
title('Comparison of solution');