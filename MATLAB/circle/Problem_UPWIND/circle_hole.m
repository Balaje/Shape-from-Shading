%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Upwind scheme for the 2D Eikonal equation    %                                                
%          on a circle with a hole - \Omega       %  
%           |\Grad u| = 1 in \Omega               %
%               u = 0 on \Gamma_1 and \Gamma_2    %  
%                                                 %  
%                                                 %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
format longe;
a = -1;
b = 1;
c = -1;
d = 1;
N = 200;
M = N;
h = (b-a)/N;
k = 0.05*h;
x = a:h:b;
y = c:h:d;
inout = zeros(N+1);
for i = 1:N+1
    for j = 1:N+1
        if (x(i)^2 + y(j)^2 >= 1)
            inout(i,j) = 1;
        end
        if(x(i)^2+y(j)^2 <= 0.01)
            inout(i,j) = 2;
        end
    end
end
u1 = zeros(N+1);
u = zeros(N+1);
error = 50;
tol = 1e-6;
count = 0;

while error>tol
%for l=1:3
    for j = 2:N
        for i = 2:N
            
            if (inout(i,j) == 0)
                % ALL YES
                if (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                %TN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 1 && inout(i,j-1) == 0)
                    eps = sqrt(1- x(i)^2) - y(j);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (- u(i,j))/eps;
                    
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %BN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 1)                    
                    eps = y(j) + sqrt(1 - x(i)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = -(u(i,j))/eps;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %LN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    eps = x(i) + sqrt(1 - y(j)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = -(u(i,j))/eps;
                    
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %RN    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    eps = - x(i) + sqrt(1 - y(j)^2);
                    
                    Dpx = (- u(i,j))/eps;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %TRN    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 1 && inout(i,j-1) == 0)                    
                    eps1 = sqrt(1- x(i)^2) - y(j);
                    eps2 = - x(i) + sqrt(1 - y(j)^2);                   
                    Dpx = (- u(i,j))/eps2;                    
                    Dpy = (- u(i,j))/eps1;                   
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %TLN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 1 && inout(i,j-1) == 0)
                    eps1 = sqrt(1- x(i)^2) - y(j);
                    eps2 = x(i) + sqrt(1 - y(j)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;                    
                    Dpy = (- u(i,j))/eps1;      
                    Dmx = (- u(i,j))/eps2;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                %BRN    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 1)
                    
                    eps1 = sqrt(1- x(i)^2) + y(j);
                    eps2 = -x(i) + sqrt(1 - y(j)^2);  
                    
                    Dpx = (-u(i,j))/eps2;                    
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;                    
                    Dmy = (-u(i,j))/eps1;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                 %BLN
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 0 && inout(i,j-1) == 1)
                    
                    eps1 = sqrt(1- x(i)^2) + y(j);
                    eps2 = x(i) + sqrt(1 - y(j)^2);
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;                    
                    Dmx = (-u(i,j))/eps2;
                    Dmy = (-u(i,j))/eps1;                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %%%--------
                % Inside cirle sweeps
                %TN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 2 && inout(i,j-1) == 0)
                    eps =  sqrt(0.01- x(i)^2) + y(j);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (- u(i,j))/abs(eps);
                    
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %BN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 2)                    
                    eps = y(j) - sqrt(0.01 - x(i)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = -(u(i,j))/abs(eps);
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %LN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 2 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    eps = - x(i) + sqrt(0.01 - y(j)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = -(u(i,j))/abs(eps);
                    
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %RN    
                elseif (inout(i+1,j) == 2 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    eps =  x(i) + sqrt(0.01 - y(j)^2);
                    
                    Dpx = (- u(i,j))/abs(eps);
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %TRN    
                elseif (inout(i+1,j) == 2 && inout(i-1,j) == 0 && inout(i,j+1) == 2 && inout(i,j-1) == 0)                    
                    eps1 = sqrt(0.01- x(i)^2) + y(j);
                    eps2 =  x(i) + sqrt(0.01 - y(j)^2);
                    
                    Dpx = (- u(i,j))/abs(eps2);                    
                    Dpy = (- u(i,j))/abs(eps1);                   
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %TLN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 2 && inout(i,j+1) == 2 && inout(i,j-1) == 0)
                    eps1 = sqrt(0.01- x(i)^2) + y(j);
                    eps2 = - x(i) + sqrt(0.01 - y(j)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;                    
                    Dpy = (- u(i,j))/abs(eps1);      
                    Dmx = (- u(i,j))/abs(eps2);
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                %BRN    
                elseif (inout(i+1,j) == 2 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 2)
                    
                    eps1 = - sqrt(0.01- x(i)^2) + y(j);
                    eps2 = x(i) + sqrt(0.01 - y(j)^2);  
                    
                    Dpx = (-u(i,j))/abs(eps2);                    
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;                    
                    Dmy = (-u(i,j))/abs(eps1);
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                 %BLN
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 2 && inout(i,j+1) == 0 && inout(i,j-1) == 2)
                    
                    eps1 = - sqrt(0.01- x(i)^2) + y(j);
                    eps2 = - x(i) + sqrt(0.01 - y(j)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;                    
                    Dmx = (-u(i,j))/abs(eps2);
                    Dmy = (-u(i,j))/abs(eps1);                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                end
                
                H = sqrt(Dx^2+Dy^2)-1;
                u1(i,j) = u(i,j) - k*H;
            end
            
        end
            
    end
    error = max(max(abs(u1-u)));
    count = count +1;
    mesh(x,y,u1);
    axis([a,b,c,d,0,1]);
    pause(0.01);
    u = u1;
end
%end
%surf(x,y,u1);
