function [V,ind] = godunov(u,v)
    V = max(abs(max(u,0)),abs(min(v,0)));
    if(V == abs(u))
        ind = 1;
    elseif(V == abs(v));
        ind = 2;
    elseif(V == 0);
        ind = 3;
    end
end


