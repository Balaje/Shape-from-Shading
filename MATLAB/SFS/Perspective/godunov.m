function [V,idx] = godunov(u,v,x,y)
% Function to compute the Godunov Flux
global upx
[V,idx] = max([J(x)^2*abs(max(u,0))^2,J(y)^2*abs(min(v,0))^2]);
if(V==0)
    idx = 3;
end
if(idx==1 || idx==3)
    upx = x;
elseif(idx==2)
    upx = y;
end