%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MATLAB Code to implement the level set method to solve             % 
%       u_t + F|\Grad u| = 0 in \Omega = ]-1,1[ x ]-1,1[ and t = (0,T] % 
%           u(x,0) = x^2+y^2-1                                         % 
%               u = 0 on boundary                                      % 
%                                                                      % 
%   Track the surface at t = 4s                                        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% Uniform grid along x and y
N = 100;

a = -1;
b = 1;
c = -1;
d = 1;

h = (b-a)/N;

x = a+h:h:b-h;
y = c+h:h:d-h;

u0 = zeros(N-1);
for i=1:N-1
    for j=1:N-1
        u0(i,j) = sqrt(x(i)^2+y(j)^2)-1^2;
    end
end

t0 = 0;
tf = 2;

% CFL Condition
delt = 0.3*h;

M = fix((tf-t0)/delt);
u = zeros(N-1);

for k = 1:M
    for i=2:N-2
        for j=2:N-2
            Dpx = (u0(i+1,j)-u0(i,j))/h;
            Dmx = (u0(i,j)-u0(i-1,j))/h;
            Dpy = (u0(i,j+1)-u0(i,j))/h;
            Dmy = (u0(i,j)-u0(i,j-1))/h;
            
            Dp = sqrt(max([Dmx,0])^2 + min([Dpx,0])^2 + max([Dmy,0])^2 + min([Dpy,0])^2);
            Dm = sqrt(max([Dpx,0])^2 + min([Dmx,0])^2 + max([Dpy,0])^2 + min([Dmy,0])^2);
            
            u(i,j) = u0(i,j) - delt*( max([F(x(i),y(j)),0])*Dp + min([F(x(i),y(j)),0])*Dm );
        end
    end
    figure(2)
    contour(y,x,u,[0,0],'bl');
    title('Level set front propogation')
    xlabel('x');
    ylabel('y');
    pause(0.01);
    u0 = u;
end