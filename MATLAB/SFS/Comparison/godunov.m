function V = godunov(u,v,I1,I2,X)
V = max(H(I1,max(u,0),X),H(I2,min(v,0),X));
    %V = max(H(I1,max(u,0),X1),H(I2,min(v,0),X1));
