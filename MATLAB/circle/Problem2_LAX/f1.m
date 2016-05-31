function V = f1(a,b,c,d,k,hx,hy,x,y)
% FUNCTION V = f(a,b,c,d,k,hx,hy)
%   a       : i+1,j
%   b       : i-1,j
%   c       : i,j+1
%   d       : i,j-1
%   k       : Time step
%   hx,hy   : Space step - x and y
  
    cx = (a-b)/(2*hx);
    cy = (c-d)/(2*hy);
    V = (a+b+c+d)/4 - k*(sqrt(cx^2+cy^2)) + k*F(x,y);
end