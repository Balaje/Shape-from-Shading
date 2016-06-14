%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Upwind scheme for the 2D Eikonal equation    %                                                
%          on a circle - \Omega                   %  
%           |\Grad u| = 1 in \Omega               %
%               u = 0 on \Gamma                   %  
%                                                 %  
%                                                 %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile on;
clear;
close all;
format long;

% Domain
a = -1;
b = 1;
c = -1;
d = 1;

% No of Divisions 
N = 100;

% Space discretization
h = (b-a)/N;

% CFL Condition
%k = 0.01*h;

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
error = 100;
tol = 1e-16;
count = 0;
tic
while error>tol
    for i = 2:N
        for j = 2:N
            k = h;
            % If inside the domain
            if (inout(i,j) == 0)
                k = h;
                Dpx = (u(i+1,j) - u(i,j))/h;
                Dpy = (u(i,j+1) - u(i,j))/h;
                Dmx = (u(i-1,j) - u(i,j))/h;
                Dmy = (u(i,j-1) - u(i,j))/h;
                % TOP NO
                if (inout(i,j+1) == 1)
                    eps = sqrt(1- x(i)^2) - y(j);
                    if (eps<k)
                        k = eps;
                    end
                    Dpy = (- u(i,j))/eps;                                     
                    
                %BOTTOM NO    
                elseif (inout(i,j-1) == 1)
                    eps = y(j) + sqrt(1 - x(i)^2);
                    if (eps<k)
                        k = eps;
                    end
                    Dmy = -(u(i,j))/eps;                    
                                  
                %LEFT NO    
                elseif (inout(i-1,j) == 1)
                    eps = x(i) + sqrt(1 - y(j)^2); 
                    if (eps<k)
                        k = eps;
                    end
                    Dmx = -(u(i,j))/eps;
                    
                %RIGHT NO    
                elseif (inout(i+1,j) == 1)
                    eps = - x(i) + sqrt(1 - y(j)^2);
                    if (eps<k)
                        k = eps;
                    end
                    Dpx = (- u(i,j))/eps;                                        
                end
%                 Dx = min([Dpx,Dmx,0]);
%                 Dy = min([Dpy,Dmy,0]);
%                 Dx = Dpx;
%                 Dy = Dpy;
%                 if(Dmx < Dpx)
%                     Dx = Dmx;
%                 elseif(0 < Dpx)
%                     Dx = 0;
%                 end
%                 
%                 if(Dmy < Dpy)
%                     Dy = Dmy;
%                 elseif(0 < Dpy)
%                     Dy = 0;
%                 end
                Dx = fluxfunction(-Dmx,Dpx,'m',h);
                Dy = fluxfunction(-Dmy,Dpy,'m',h);
                H = sqrt(Dx^2+Dy^2)-1;
                
                k = 0.7*k;
                u1(i,j) = u(i,j) - k*H;
            end   
        end
    end
    
    error = max(max(abs(u1-u)));    
    u = u1;
    count = count+1;
    surf(x,y,u1);
    pause(0.01);
end
toc
profile viewer;