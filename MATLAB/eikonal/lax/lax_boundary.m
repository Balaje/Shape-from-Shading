%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lax Friedrich like scheme for a                 %  
%   Initial Boundary Value Problem                %  
%       u_t+|u_x|=1 in (-1,1)                     %  
%           u(x,0) = 0                            %
%           u(-1,t) = u(1,t) = 0                  %  
%                                                 %  
%           At steady state u(x) = 1-|x|          %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10;
a = 0;
b = 1;

h = (b-a)/N;
x = (a+h:h:b-h)';

delt = h;

t0 = 0;
tf = 0.3;

M = ((tf-t0)/delt);
u0 = x.^2;
u = u0;

error = 100;
tol = 10^-3;

iter = 1;
while error > tol
    %for k=1:M+1
        u(1) = min([u0(2),0])+delt;
        for j=2:N-2
            u(j) = min([u0(j+1),u0(j-1)])+delt;
        end
        u(N-1) = min([u0(N-2),0])+delt;
        
        plot(x,u);
        grid on
        axis([0,1,0,0.5])
        title(num2str(iter*delt));
        pause(1);
        
        error = max(abs(u-u0));
        iter = iter+1;
        u0 = u;
    %end
end