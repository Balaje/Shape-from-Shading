clear
close all

load op_neumannf.txt
load op_neumann.txt
load op_dirich.txt

Z1 = op_neumannf(:,3);
Z2 = op_dirich(:,3);
Z3 = op_neumann(:,3);

z1 = reshape(Z1,[201,201]);
z2 = reshape(Z2,[201,201]);
z3 = reshape(Z3,[201,201]);

x = 0:1/200:1;
y = 0:1/200:1;

figure(1)
surf(x,y,z1)
title('Mozart with full neumann condition');

figure(2)
surf(x,y,z2)
title('Mozart with dirichlet condition');

figure(3)
surf(x,y,z3)
title('Mozart with neumann condition');