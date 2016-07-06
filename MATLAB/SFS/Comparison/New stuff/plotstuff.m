function plotstuff()

global unew a b exactsol

N = length(unew)+1;

h = (b-a)/N;
x = h:h:1-h;
int = zeros(N-1,1);

figure(1)
plot(x,unew,x,exactsol);
legend('Approximate solution','Exact solution');

for i=1:N-1
    int(i) = I(x(i));
end

figure(2)
plot(x,int,'*');
title('Intensity function')