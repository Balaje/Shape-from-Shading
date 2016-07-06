function plotstuff()

global a b

unew1 = upwind();
unew = woupwind();

N = length(unew)+1;


h = (b-a)/N;
x = h:h:1-h;

exactsol = zeros(N-1,1);

for i=1:N-1
    exactsol(i) = exact(x(i));
end

int = zeros(N-1,1);

figure(1)
plot(x,unew1,x,unew,x,exactsol);
legend('Solution with upwinding','Solution w/o upwinding','Exact solution');

for i=1:N-1
    int(i) = I(x(i));
end

figure(2)
plot(x,int,'*');
title('Intensity function')