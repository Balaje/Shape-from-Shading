function V = H(J,p,x)
    global f
    V = sqrt(J^2*f^2*p^2 + J^2*p^2*x^2 + J^2*f^2/(f^2+x^2));

