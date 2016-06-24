function V = J(x)
    %V = 1/sqrt(1+(pi*cos(pi*x))^2);
    %V = 1/sqrt(1+((1-2*x))^2);
    V = 1/(sqrt(x^2+1)*sqrt((x^2+1)*pi^2*(cos(pi*x))^2 + 1/(x^2+1)));
end