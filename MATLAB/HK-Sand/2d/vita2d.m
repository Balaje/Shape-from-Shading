%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Finzi Vita - 2d finite difference scheme   %
%                                              % 
%                                              % 
%                                              %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

clear 
close all
clc

format long

N = 100;

a = 0; b = 1;
c = 0; d = 1;

h = (b-a)/N;

x = (a+h:h:b-h);
y = (c+h:h:d-h);

u = zeros(N-1);
v = zeros(N-1);
unew = zeros(N-1);
vnew = zeros(N-1);
f = zeros(N-1);

for i=1:N-1
    for j=1:N-1
        if((x(i) >= 0.1 && x(i) <= 0.3) && (y(j) >= 0.5 && y(j) <= 0.7))
            f(i,j) = 1;
        end
    end    
end
delt = 0.5*h;

M = 1000;

for k = 1:M
    % i = 1, j = 1
    Dpx = (u(2,1)-u(1,1))/h;
    Dmx = (u(1,1))/h;
    Dpy = (u(1,2)-u(1,1))/h;
    Dmy = (u(1,1))/h;
    
    [~,Dux] = maxmod(Dpx,Dmx);
    [~,Duy] = maxmod(Dpy,Dmy);
    
    du = Du(Dux,Duy);
    
    up1x = H(u(1,1),u(2,1),v(1,1),v(2,1));
    up2x = v(1,1);
    up1y = H(u(1,1),u(1,2),v(1,1),v(1,2));
    up2y = v(1,1);
    
    h1x = -Dpx*up1x;
    h2x = -Dmx*up2x;
    h1y = -Dpy*up1y;
    h2y = -Dmy*up2y;
    
    vnew(1,1) = v(1,1) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(1,1) - (1-du)*v(1,1));
    unew(1,1) = u(1,1) + delt*(1-du)*vnew(1,1);
    
    % i = 1, j = 2:N-2
    for j = 2:N-2
        Dpx = (u(2,j)-u(1,j))/h;
        Dmx = (u(1,j))/h;
        Dpy = (u(1,j+1)-u(1,j))/h;
        Dmy = (u(1,j)-u(1,j-1))/h;
    
        [~,Dux] = maxmod(Dpx,Dmx);
        [~,Duy] = maxmod(Dpy,Dmy);
    
        du = Du(Dux,Duy);
    
        up1x = H(u(1,j),u(2,j),v(1,j),v(2,j));
        up2x = v(1,j);
        up1y = H(u(1,j),u(1,j+1),v(1,j),v(1,j+1));
        up2y = H(u(1,j-1),u(1,j),v(1,j-1),v(1,j));
    
        h1x = -Dpx*up1x;
        h2x = -Dmx*up2x;
        h1y = -Dpy*up1y;
        h2y = -Dmy*up2y;
    
        vnew(1,j) = v(1,j) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(1,j) - (1-du)*v(1,j));
        unew(1,j) = u(1,j) + delt*(1-du)*vnew(1,j);
    end
    
    % i = 1, j = N-1
    Dpx = (u(2,N-1)-u(1,N-1))/h;
    Dmx = (u(1,N-1))/h;
    Dpy = (-u(1,N-1))/h;
    Dmy = (u(1,N-1)-u(1,N-2))/h;
    
    [~,Dux] = maxmod(Dpx,Dmx);
    [~,Duy] = maxmod(Dpy,Dmy);
    
    du = Du(Dux,Duy);
    
    up1x = H(u(1,N-1),u(2,N-1),v(1,N-1),v(2,N-1));
    up2x = v(1,N-1);
    up1y = v(1,N-1);
    up2y = H(u(1,N-2),u(1,N-1),v(1,N-2),v(1,N-1));
    
    h1x = -Dpx*up1x;
    h2x = -Dmx*up2x;
    h1y = -Dpy*up1y;
    h2y = -Dmy*up2y;
    
    vnew(1,N-1) = v(1,N-1) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(1,N-1) - (1-du)*v(1,N-1));
    unew(1,N-1) = u(1,N-1) + delt*(1-du)*vnew(1,N-1);
    
    %i = 2:N-2, j = 2:N-2
    for i = 2:N-2
        for j = 2:N-2
            Dpx = (u(i+1,j)-u(i,j))/h;
            Dmx = (u(i,j)-u(i-1,j))/h;
            Dpy = (u(i,j+1)-u(i,j))/h;
            Dmy = (u(i,j)-u(i,j-1))/h;
    
            [~,Dux] = maxmod(Dpx,Dmx);
            [~,Duy] = maxmod(Dpy,Dmy);
    
            du = Du(Dux,Duy);
    
            up1x = H(u(i,j),u(i+1,j),v(i,j),v(i+1,j));
            up2x = H(u(i-1,j),u(i,j),v(i-1,j),v(i,j));
            up1y = H(u(i,j),u(i,j+1),v(i,j),v(i,j+1));
            up2y = H(u(i,j-1),u(i,j),v(i,j-1),v(i,j));
    
            h1x = -Dpx*up1x;
            h2x = -Dmx*up2x;
            h1y = -Dpy*up1y;
            h2y = -Dmy*up2y;
    
            vnew(i,j) = v(i,j) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(i,j) - (1-du)*v(i,j));
            unew(i,j) = u(i,j) + delt*(1-du)*vnew(i,j);
        end
    end
    
    % i = N-1, j = 1
    Dpx = (-u(N-1,1))/h;
    Dmx = (u(N-1,1)-u(N-1,1))/h;
    Dpy = (u(N-1,2)-u(N-1,1))/h;
    Dmy = (u(N-1,1))/h;
    
    [~,Dux] = maxmod(Dpx,Dmx);
    [~,Duy] = maxmod(Dpy,Dmy);
    
    du = Du(Dux,Duy);
    
    up1x = v(N-1,1);
    up2x = H(u(N-2,1),u(N-1,1),v(N-2,1),v(N-1,1));
    up1y = H(u(N-1,1),u(N-1,2),v(N-1,2),v(N-1,2));
    up2y = v(N-1,1);
    
    h1x = -Dpx*up1x;
    h2x = -Dmx*up2x;
    h1y = -Dpy*up1y;
    h2y = -Dmy*up2y;
    
    vnew(N-1,1) = v(N-1,1) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(N-1,1) - (1-du)*v(N-1,1));
    unew(N-1,1) = u(N-1,1) + delt*(1-du)*vnew(N-1,1);
    
    % i = N-1, j = 2:N-2
    for j = 2:N-2
        Dpx = (-u(N-1,j))/h;
        Dmx = (u(N-1,j)-u(N-2,j))/h;
        Dpy = (u(N-1,j+1)-u(N-1,j))/h;
        Dmy = (u(N-1,j)-u(N-1,j-1))/h;
    
        [~,Dux] = maxmod(Dpx,Dmx);
        [~,Duy] = maxmod(Dpy,Dmy);
    
        du = Du(Dux,Duy);
        
        up1x = v(N-1,j);
        up2x = H(u(N-2,j),u(N-1,j),v(N-2,j),v(N-1,j));
        up1y = H(u(N-1,j),u(N-1,j+1),v(N-1,j),v(N-1,j+1));
        up2y = H(u(N-1,j-1),u(N-1,j),v(N-1,j-1),v(N-1,j));
    
        h1x = -Dpx*up1x;
        h2x = -Dmx*up2x;
        h1y = -Dpy*up1y;
        h2y = -Dmy*up2y;
        
        vnew(N-1,j) = v(N-1,j) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(N-1,j) - (1-du)*v(N-1,j));
        unew(N-1,j) = u(N-1,j) + delt*(1-du)*vnew(N-1,j);
    end
    % i = N-1, j = N-1
    Dpx = (-u(N-1,N-1))/h;
    Dmx = (u(N-1,N-1)-u(N-2,N-1))/h;
    Dpy = (-u(N-1,N-1))/h;
    Dmy = (u(N-1,N-1)-u(N-1,N-2))/h;
    
    [~,Dux] = maxmod(Dpx,Dmx);
    [~,Duy] = maxmod(Dpy,Dmy);
    
    du = Du(Dux,Duy);
    
    up1x = v(N-1,N-1);
    up2x = H(u(N-2,N-1),u(N-1,N-1),v(N-2,N-1),v(N-1,N-1));
    up1y = v(N-1,N-1);
    up2y = H(u(N-1,N-2),u(N-1,N-1),v(N-1,N-2),v(N-1,N-1));
    
    h1x = -Dpx*up1x;
    h2x = -Dmx*up2x;
    h1y = -Dpy*up1y;
    h2y = -Dmy*up2y;
    
    vnew(N-1,N-1) = v(N-1,N-1) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(N-1,N-1) - (1-du)*v(N-1,N-1));
    unew(N-1,N-1) = u(N-1,N-1) + delt*(1-du)*vnew(N-1,N-1);
    
    % j = 1 ,i = 2:N-2
    for i = 2:N-2
        Dpx = (u(i+1,1)-u(i,1))/h;
        Dmx = (u(i,1)-u(i-1,1))/h;
        Dpy = (u(i,2)-u(i,1))/h;
        Dmy = (u(i,1))/h;
        
        [~,Dux] = maxmod(Dpx,Dmx);
        [~,Duy] = maxmod(Dpy,Dmy);
        
        du = Du(Dux,Duy);
        
        up1x = H(u(i,1),u(i+1,1),v(i,1),v(i+1,1));
        up2x = H(u(i-1,1),u(i,1),v(i-1,1),v(i,1));
        up1y = H(u(i,1),u(i,2),u(i,1),v(i,2));
        up2y = v(i,1);
        
        h1x = -Dpx*up1x;
        h2x = -Dmx*up2x;
        h1y = -Dpy*up1y;
        h2y = -Dmy*up2y;
        
        vnew(i,1) = v(i,1) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(i,1) - (1-du)*v(i,1));
        unew(i,1) = u(i,1) + delt*(1-du)*vnew(i,1);
    end
    
    % j = N-1, i = 2:N-2
    for i = 2:N-2
        Dpx = (u(i+1,N-1)-u(i,N-1))/h;
        Dmx = (u(i,N-1)-u(i-1,N-1))/h;
        Dpy = (-u(i,N-1))/h;
        Dmy = (u(i,N-1)-u(i,N-2))/h;
        
        [~,Dux] = maxmod(Dpx,Dmx);
        [~,Duy] = maxmod(Dpy,Dmy);
        
        du = Du(Dux,Duy);
        
        up1x = H(u(i,N-1),u(i+1,N-1),v(i,N-1),v(i+1,N-1));
        up2x = H(u(i-1,N-1),u(i,N-1),v(i-1,N-1),v(i,N-1));
        up1y = v(i,N-1);
        up2y = H(u(i,N-2),u(i,N-1),v(i,N-2),v(i,N-1));
        
        h1x = -Dpx*up1x;
        h2x = -Dmx*up2x;
        h1y = -Dpy*up1y;
        h2y = -Dmy*up2y;
        
        vnew(i,N-1) = v(i,N-1) - delt/h*(h1x-h2x+h1y-h2y) + delt*(f(i,N-1) - (1-du)*v(i,N-1));
        unew(i,N-1) = u(i,N-1) + delt*(1-du)*vnew(i,N-1);        
    end
    u = unew;
    v = vnew;
    
    surf(x,y,u);
    %plot(x,v,'r');
    axis([a,b,c,d,0,0.4]);
    pause(0.001);
end

