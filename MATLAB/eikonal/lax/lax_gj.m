%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Lax Friedrich like scheme             %  
%       u_t + |u_x| = 0                   %
%           u(x,0) = u_0(x) = |x|         %  
%                                         %  
%   Exact u(x,t) = min{x-t<y<x+t} (u0(y)) %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 6000;
a = -4*pi;
b = 4*pi;

t0 = 0;
tf = 3*pi;

h = (b-a)/N;

delt = h;

u0 = zeros(N-1,1);
u = u0;
x = (a+h:h:b-h)';

for j=1:N-1
    %u0(j) = abs(x(j));
    u0(j) = sin(x(j));
end

M = (tf-t0)/delt;

for k=1:M
    u(1) = min([sin(a),u0(2)]);
    for j=2:N-2
        u(j) = min([u0(j-1),u0(j+1)]);
    end
    u(N-1) = min([u0(N-2),sin(b)]);
%     plot(x,u);
%     pause(0.01);
    u0 = u;
end