function V = exact(x)
    if(x <= 0.4)
        V = x*(1-x);
    else
        V = 0.15;
    end
end