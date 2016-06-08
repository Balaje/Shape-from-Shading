function [V,val] = maxmod(u,v)
    arr = [u,v];
    [m,I] = max(abs(arr));
    V = arr(I);
    val = m;
end