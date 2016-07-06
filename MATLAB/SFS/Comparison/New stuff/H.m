function V = H(x,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the Hamiltonian %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V = J(x)*sqrt(abs(p)^2 + x^2*p^2 + Q(x)^2);
end