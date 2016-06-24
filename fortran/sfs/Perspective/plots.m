clear
close all

load op_neumannf.txt
load op_neumann.txt
load op_dirich.txt

load op_neumannf1.txt
load op_neumann1.txt
load op_dirich1.txt

Z1 = op_neumannf(:,3);
Z2 = op_dirich(:,3);
Z3 = op_neumann(:,3);

Z4 = op_neumannf1(:,3);
Z5 = op_dirich1(:,3);
Z6 = op_neumann1(:,3);

z1 = reshape(Z1,[201,201]);
z2 = reshape(Z2,[201,201]);
z3 = reshape(Z3,[201,201]);

z4 = reshape(Z4,[201,201]);
z5 = reshape(Z5,[201,201]);
z6 = reshape(Z6,[201,201]);

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

figure(4)
surf(x,y,z4)
title('(Upwinded)Mozart with full neumann condition');

figure(5)
surf(x,y,z5)
title('(Upwinded)Mozart with dirichlet condition');

figure(6)
surf(x,y,z6)
title('(Upwinded)Mozart with neumann condition');