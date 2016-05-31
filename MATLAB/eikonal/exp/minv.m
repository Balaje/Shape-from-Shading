function V = minv(y)
    V0 = @(x) x.^2;
    V = min(V0(y));
end