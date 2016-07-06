function V = godunov(u,v,xi,xi1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the Godunov Flux -  %
%   godunov(Du-,Du+,x(i),x(i+1),u(i),u(i+1))%
%                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   V = max(H(xi,max(u,0)),H(xi1,min(v,0)));