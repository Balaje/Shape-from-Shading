function V = H(u,v,x,y)
    if(u <= v)
        V = y;
    elseif(u > v)
        V = x;
    end
end