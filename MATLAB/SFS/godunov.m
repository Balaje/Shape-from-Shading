function V = godunov(u,v,x,y)
% Function to compute the Godunov Flux
global upx
[V,idx] = max([J(x)*abs(max(u,0)),J(y)*abs(min(v,0))]);
if(idx == 1)
    upx = x;
else
    upx = y;
end
