%function V = upwind()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to check the effect of upwinding on   %
%           Discontinous Intensity function     %
%   Run : With upwinding                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
close all
clear

global unew exactsol a b h

a = 0; b = 1;
N = 100;
h = (b-a)/N;

x = a+h:h:b-h;

Int = zeros(N-1,1);
for i=1:N-1
    Int(i) = I(x(i));
end
discs = nonzeros(finddisc(Int));
dc = reshape(discs,[2,length(discs)/2])';
if(~isempty(dc))
    discl = dc(:,1);
    discr = dc(:,2);
else
    discl = N;
    discr = N;
end

memberr = ismember(1:N-1,discr);
memberl = ismember(1:N-1,discl);

u = zeros(N-1,1);
unew = zeros(N-1,1);
exactsol = zeros(N-1,1);

error = 100;
tol = 1e-10;

delt = 0.7*h;
iterations = 0;

eps = 0;

while error > tol
    % i = 1
    unew(1) = (0 + u(1))/2 - delt*H((x(1)+a)/2,(u(2)+0)/2,(u(2)-0)/(2*h));
    
    for i=2:N-2
        unew(i) = (u(i-1)+u(i+1))/2 - delt*H((x(i+1)+x(i-1))/2,(u(i+1)+u(i-1))/2,(u(i+1)-u(i-1))/(2*h));
    end
    %i = N-1
    unew(N-1) = (u(N-2)+0)/2 - delt*H((b+x(N-2))/2,(0+u(N-2))/2,(-u(N-2))/(2*h));
       
    
    error = max(abs(u-unew));
    for i=1:N-1
        exactsol(i) = exact(x(i));
    end
    
    figure(2)
    plot(x,unew,'*',x,exactsol);
    %axis([a,b,0,0.3])
    pause(0.001);
    
    iterations = iterations + 1;
    u = unew;
end

plot(x,unew,x,exactsol);
max(abs(unew-exactsol))
V = unew;
