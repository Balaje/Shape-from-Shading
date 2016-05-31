%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Eikonal Equation on a circle    %
%       |\Grad u| = 1 in \Omega     %
%           u = 0 on \Gamma         %
%                                   %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
format long

a = -1;
b = 1;
c = -1;
d = 1;

N = 50;

h = (b-a)/N;

%Condition
k = 0.1*h;

x = a:h:b;
y = c:h:d;

inout = zeros(N+1);

for i=1:N+1
    for j=1:N+1
        if(x(i)^2+y(j)^2 >= 1)
            inout(i,j) = 1;
        end
    end
end

error = 100;
tol = 10^-8;

u0 = zeros(N+1);
u = u0;

count = 0;
while error > tol
    for i=2:N
        for j=2:N
            if(inout(i,j)==0)
                %ALL YES
                if(inout(i-1,j)==0 && inout(i+1,j)==0 && inout(i,j+1)==0 && inout(i,j-1)==0)
                    u(i,j) = f(u0(i+1,j),u0(i-1,j),u0(i,j+1),u0(i,j-1),k,h,h);
                    
                %TN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 1 && inout(i,j-1) == 0)
                    eps = sqrt(1- x(i)^2) - y(j);
                    u(i,j) = f(u0(i+1,j),u0(i-1,j),0,u0(i,j-1),k,h,eps);
                    
                %BN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 1)
                    eps = y(j) + sqrt(1 - x(i)^2);
                    u(i,j) = f(u0(i+1,j),u0(i-1,j),u0(i,j+1),0,k,h,eps);
                                        
                %LN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    eps = x(i) + sqrt(1 - y(j)^2);
                    u(i,j) = f(u0(i+1,j),0,u0(i,j+1),u0(i,j-1),k,eps,h);        
                                        
                %RN    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 0)
                    eps = - x(i) + sqrt(1 - y(j)^2);
                    u(i,j) = f(0,u0(i-1,j),u0(i,j+1),u0(i,j-1),k,eps,h);
                    
                %TRN    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 1 && inout(i,j-1) == 0)
                    eps1 = sqrt(1- x(i)^2) - y(j); %hy
                    eps2 = - x(i) + sqrt(1 - y(j)^2); %hx                  
                    u(i,j) = f(0,u0(i-1,j),0,u0(i,j-1),k,eps2,eps1);
                    
                %TLN    
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 1 && inout(i,j-1) == 0)
                    eps1 = sqrt(1- x(i)^2) - y(j); %hy
                    eps2 = x(i) + sqrt(1 - y(j)^2); %hx
                    u(i,j) = f(u0(i+1,j),0,0,u0(i,j-1),k,eps2,eps1);
                    
                   
                %BRN    
                elseif (inout(i+1,j) == 1 && inout(i-1,j) == 0 && inout(i,j+1) == 0 && inout(i,j-1) == 1)
                    eps1 = sqrt(1- x(i)^2) + y(j); %hy
                    eps2 = -x(i) + sqrt(1 - y(j)^2); %hx                  
                    u(i,j) = f(0,u0(i-1,j),u0(i,j+1),0,k,eps2,eps1);
                    
                 %BLN
                elseif (inout(i+1,j) == 0 && inout(i-1,j) == 1 && inout(i,j+1) == 0 && inout(i,j-1) == 1)
                    eps1 = sqrt(1- x(i)^2) + y(j); %hy
                    eps2 = x(i) + sqrt(1 - y(j)^2); %hx
                    u(i,j) = f(u0(i+1,j),0,u0(i,j+1),0,k,eps2,eps1);                    
               
                end                
            end
        end
    end
    mesh(x,y,u);
    pause(0.01);
    error = max(max(abs(u-u0)));
    u0 = u;
    count = count+1;
end