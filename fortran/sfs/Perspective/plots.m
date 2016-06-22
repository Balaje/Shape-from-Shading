clear
close all

load op_neumann.txt
load op_dirich.txt

Z1 = op_neumann(:,3);
Z2 = op_dirich(:,3);

z1 = reshape(Z1,[201,201]);
z2 = reshape(Z2,[201,201]);

x = 0:1/200:1;
y = 0:1/200:1;

figure(1)
surf(x,y,z1)
title('Mozart with neumann condition');

figure(2)
surf(x,y,z2)
title('Mozart with dirichlet condition');