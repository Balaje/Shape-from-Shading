clear all;
close all;
clc;
format long;
a = 0;
L = 1;
b = L;
N = 10;
h = (b-a)/N;
u = zeros(N+1,1);
t = 1;
x = a:h:b;
for i = 2:N
    if ((x(i)-t>0) && (x(i)+t<L))
        del = 2*t/10;
        y = x(i)-t:del:x(i)+t;
        u(i) = minv(y) + t;
        
    elseif ((x(i)-t<0) && (x(i)+t<L))
        del = (t+x(i))/10;
        y = 0:del:t+x(i);
        u(i) = min([0,minv(y)] + t);
        
    elseif ((x(i)-t>0) && (x(i)+t>L))
        del = (L+t-x(i))/10;
        y = x(i)-t:del:0;
        u(i) = min([0,minv(y)] + t);  
        
    elseif ((x(i)-t<0) && (x(i)+t>L))
        del = (L)/10;
        y = 0:del:L;
        u(i) = min([0,minv(y)]) + t;
        
    end
end
plot(x,u);