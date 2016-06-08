%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program to solve the open table problem    %                                       
%                                              % 
%                                              % 
%                                              % 
%                                              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear 
close all
clc

format long

a = 0;
b = 1;

N = 100;

h = (b-a)/N;

x = (a+h:h:b-h)';

u = zeros(N-1,1);
v = zeros(N-1,1);
unew = u;
vnew = v;
f = u;
delt = 1.1*h;

%M = 34.965/delt;
M = 1000;

for j = 1:N-1
    if((x(j) >= 0.18 && x(j) <= 0.22) || (x(j) >= 0.78 && x(j) <= 0.82))
        f(j) = 1;
    end
end
for k = 1:M
    % i = 1
    Dpu = (u(2) - u(1))/h;
    Dmu = (u(1))/h;
    arr = [Dpu,Dmu];
    [m1,I] = max(abs(arr));
    Du = arr(I);
    up1 = H(u(1),u(2),v(1),v(2));
    up2 = v(1);
    h1 = -Dpu*up1;
    h2 = -Dmu*up2;
    vnew(1) = v(1) - delt/h*(h1-h2) + delt*(f(1) - (1-abs(Du))*v(1));
    unew(1) = u(1) + delt*(1-abs(Du))*vnew(1);
    
    %2-N-2
    for i=2:N-2
        Dpu = (u(i+1) - u(i))/h;
        Dmu = (u(i) - u(i-1))/h;
        arr = [Dpu,Dmu];
        [m2,I] = max(abs(arr));
        Du = arr(I);
        up1 = H(u(i),u(i+1),v(i),v(i+1));
        up2 = H(u(i-1),u(i),v(i-1),v(i));
        h1 = -Dpu*up1;
        h2 = -Dmu*up2;
        vnew(i) = v(i) - delt/h*(h1-h2) + delt*(f(i) - (1-abs(Du))*v(i));
        unew(i) = u(i) + delt*(1-abs(Du))*vnew(i);
    end
    
    %N-1
    Dpu = ( - u(N-1))/h;
    Dmu = (u(N-1) - u(N-2))/h;
    arr = [Dpu,Dmu];
    [m,I] = max(abs(arr));
    Du = arr(I);
    up1 = v(N-1);
    up2 = H(u(N-2),u(N-1),v(N-2),v(N-1));
    h1 = -Dpu*up1;
    h2 = -Dmu*up2;
    vnew(N-1) = v(N-1) - delt/h*(h1-h2) + delt*(f(N-1) - (1-abs(Du))*v(N-1));
    unew(N-1) = u(N-1) + delt*(1-abs(Du))*vnew(N-1);
    
    v = vnew;
    u = unew;
    
    plot(x,u,'b');
    %plot(x,v,'r');
    axis([a,b,0,0.5]);
    pause(0.001);
    %hold on;
end

