%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program to check the effect of upwinding on
%
%           Discontinous Intensity function
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
close all
clear
a = 0; b = 1;
N = 50;
h = (b-a)/N;

x = a+h:h:b-h;

I = zeros(N-1,1);
exact = zeros(N-1,1);

for i=1:N-1
        if(x(i) < 0.4 || x(i) > 0.6)
            I(i) = exp(-2*x(i)*(1-x(i)))/sqrt(1+((1-2*x(i))*(x(i)^2+1))^2);
            exact(i) = x(i)*(1-x(i));
        else
            I(i) = exp(-4*x(i)*(1-x(i)))/sqrt(1+((2-4*x(i))*(x(i)^2+1))^2);
            exact(i) = 2*x(i)*(1-x(i));
    
        end
%         if(x(i) <= 0.4 || x(i) >= 0.6)
%             I(i) = exp(-2*x(i)*(1-x(i)))/sqrt(1+((1-2*x(i))*(x(i)^2+1))^2);
%             exact(i) = (1-x(i))*x(i);
%         else
%             I(i) = 1;
%             exact(i) = 0;
%     
%         end
%     
end

discs = nonzeros(finddisc(I));
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

J = zeros(N-1,1);
Q = zeros(N-1,1);
f = 1;

for i=1:N-1
    Q(i) = sqrt(f^2/(x(i)^2+f^2));
    J(i) = I(i)*f^2/Q(i);
end
Qb = sqrt(f^2/(b^2+f^2));
Ib = exp(-2*(1-b))/(1*sqrt(1+b^2)*sqrt((1)^2+(b)^2+(1/(b^2+1))));
Jb = Ib*f^2/Qb;

u = 0.5*log(I*f^2);
unew = zeros(N-1,1);
error = 100;
tol = 1e-8;

delt = 0.5*h;
iterations = 0;
while error > tol
    % i = 1
    Dxp = (u(2)-u(1))/h;
    Dxm = (u(1)-0)/h;
    [D,ind] = godunov(J(1)*Dxm,J(2)*Dxp);
    A = [J(1)*Dxm,J(2)*Dxp,0];
    
    H = -exp(-2*u(1)) + sqrt( f^2*D^2 + (A(ind)*x(1))^2 + J(1)^2*Q(1)^2);
    unew(1) = u(1) - delt*H;
    
    for j=2:N-2
        if(memberr(j) == 1)
            D = (u(j+1)-u(j))/h;
            H = -exp(-2*u(j)) + sqrt( J(j)^2*f^2*D^2 + (J(j)*D*x(j))^2 + J(j)^2*Q(j)^2);
            unew(j) = u(j) - delt*H;
        elseif(memberl(j) == 1)
            D = (u(j)-u(j-1))/h;
            H = -exp(-2*u(j)) + sqrt( J(j)^2*f^2*D^2 + (J(j)*D*x(j))^2 + J(j)^2*Q(j)^2);
            unew(j) = u(j)- delt*H;
        else
            Dxp = (u(j+1)-u(j))/h;
            Dxm = (u(j)-u(j-1))/h;
            [D,ind] = godunov(J(j)*Dxm,J(j+1)*Dxp);
            
            A = [J(j)*Dxm,J(j+1)*Dxp,0];
            
            
            H = -exp(-2*u(j)) + sqrt( f^2*D^2 + (A(ind)*x(j))^2 +  J(j)^2*Q(j)^2);
            unew(j) = u(j) - delt*H;
        end
    end
    
    % i = N-1
    Dxp = (-u(N-1))/h;
    Dxm = (u(N-1)-u(N-2))/h;
    [D,ind] = godunov(J(N-1)*Dxm,Jb*Dxp);
    A = [J(N-1)*Dxm,Jb*Dxp,0];
    
    H = -exp(-2*u(N-1)) + sqrt( f^2*D^2 + (A(ind)*x(N-1))^2 +  J(N-1)^2*Q(N-1)^2);
    unew(N-1) = u(N-1) - delt*H;
    
    
    error = max(abs(u-unew));
    
    figure(2)
    plot(x,unew,'-',x,exact);
    pause(0.001);
    iterations = iterations + 1;
    u = unew;
end

figure(1)
plot(x,J,'*')
figure(2)
plot(x,unew,'-',x,exact);
V = unew;
%}