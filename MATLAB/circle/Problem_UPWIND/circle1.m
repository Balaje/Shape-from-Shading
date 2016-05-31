%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Upwind scheme for the 2D Eikonal equation    %                                                
%          on a circle - \Omega                   %  
%           |\Grad u| = 1 in \Omega               %
%               u = 0 on \Gamma                   %  
%                                                 %  
%                                                 %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
format long;

% Domain
a = -1;
b = 1;
c = -1;
d = 1;

% No of Divisions 
N =50;

% Space discretization
h = (b-a)/N;

% CFL Condition
k = 0.1*h;

x = a:h:b;
y = c:h:d;

% Inout matrix describing the domain - circle.
inout = zeros(N+1);
for i = 1:N+1
    for j = 1:N+1
        if (x(i)^2 + y(j)^2 >= 1)
            inout(i,j) = 1;
        end
    end
end

u1 = zeros(N+1);
u = u1;

% Parameters
error = 50;
tol = 1e-9;
count = 0;

tic
while error>tol
    for j = 2:N
        for i = 2:N
            % If inside the domain
            if (inout(i,j) == 0)
                % ALL YES
                if (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                %TOP NO    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 1 && inout(i,j-1) == 0)
                    eps = sqrt(1- x(i)^2) - y(j);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (- u(i,j))/eps;                    
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %BOTTOM NO    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 1)
                    eps = y(j) + sqrt(1 - x(i)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = -(u(i,j))/eps;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %LEFT NO    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    eps = x(i) + sqrt(1 - y(j)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = -(u(i,j))/eps;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %RIGHT NO    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    eps = - x(i) + sqrt(1 - y(j)^2);
                    
                    Dpx = (- u(i,j))/eps;
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %TOP-RIGHT NO    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 1 && inout(i,j-1) == 0)
                    eps1 = sqrt(1- x(i)^2) - y(j);
                    eps2 = - x(i) + sqrt(1 - y(j)^2); 
                    
                    Dpx = (- u(i,j))/eps2;                    
                    Dpy = (- u(i,j))/eps1;                   
                    Dmx = (u(i-1,j) - u(i,j))/h;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %TOP-LEFT NO    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 1 && inout(i,j-1) == 0)
                    eps1 = sqrt(1- x(i)^2) - y(j);
                    eps2 = x(i) + sqrt(1 - y(j)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;                    
                    Dpy = (- u(i,j))/eps1;      
                    Dmx = (- u(i,j))/eps2;
                    Dmy = (u(i,j-1) - u(i,j))/h;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                %BOTTOM-RIGHT NO    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 1)
                    eps1 = sqrt(1- x(i)^2) + y(j);
                    eps2 = -x(i) + sqrt(1 - y(j)^2);  
                    
                    Dpx = (-u(i,j))/eps2;                    
                    Dpy = (u(i,j+1) - u(i,j))/h;
                    Dmx = (u(i-1,j) - u(i,j))/h;                    
                    Dmy = (-u(i,j))/eps1;
                    
                    Dx = min([Dpx,Dmx,0]);
                    Dy = min([Dpy,Dmy,0]);
                    
                 %BOTTOM-LEFT NO
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 0 && inout(i,j-1) == 1)
                    eps1 = sqrt(1- x(i)^2) + y(j);
                    eps2 = x(i) + sqrt(1 - y(j)^2);
                    
                    Dpx = (u(i+1,j) - u(i,j))/h;
                    Dpy = (u(i,j+1) - u(i,j))/h;                    
                    Dmx = (-u(i,j))/eps2;
                    Dmy = (-u(i,j))/eps1;                    
                    
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
    surf(x,y,u1);
    axis([a,b,c,d,0,1]);
    pause(0.01);
   
    u = u1;
end

% for i=1:N+1
%     for j=1:N+1
%         z(i,j) = 1-sqrt(x(i)^2+y(j)^2);
%         if(z(i,j) < 0)
%             z(i,j) = 0;
%         end
%     end
% end
% 
% error1 = max(max(abs(u1-z)));

toc

contour(x,y,u1);