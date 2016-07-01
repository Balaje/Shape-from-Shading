%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script to compare the solution plots of 
%       Upwinding vs No-upwinding
%
%   Functions : testdisc.m and testdisc_upw.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

N = 100;
h = 1/N;

x = h:h:1-h;

I = zeros(N-1,1);

for i=1:N-1
    if(x(i) < 0.3)
        I(i) = x(i)^2;
    elseif(x(i) >= 0.3 && x(i) < 0.6)
        I(i) = 0.2;
    elseif(x(i) >= 0.6 && x(i) < 0.8)
        I(i) = x(i);
    elseif(x(i) == 0.8 )
        I(i) = 0;
    elseif(x(i) > 0.8 && x(i) < 1)
        I(i) = 1-x(i)^2;
    end
end

sol_upw = testdisc_upw(N,I);
sol = testdisc(N,I);

figure(1)
plot(x,I);
title('Intensity Map');

figure(2)
plot(x,sol_upw,x,sol);
legend('Upwinded Solution','Not upwinded')
title('Comparison of solution');