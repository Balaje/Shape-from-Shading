function V = I(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the Intensity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    global h
%    x1 = 0.4 - 4*h;
    if(x < 0.4)
        V = exp(-2*x*(1-x))/sqrt(1 + ((1-2*x)*(x^2+1))^2);
%    elseif (x >= x1 && x <= 0.4)
%        V = ((x - 0.4 + 4*h)/(4*h))*(exp(-2*0.15) - exp(-2*x1*(1-x1))/sqrt(1 + ((1-2*x1)*(x1^2+1))^2))...
%            + exp(-2*x1*(1-x1))/sqrt(1 + ((1-2*x1)*(x1^2+1))^2);
    else
        V = exp(-2*0.15);
    end
end