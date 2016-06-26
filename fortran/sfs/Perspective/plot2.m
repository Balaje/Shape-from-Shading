load op_neumannf.txt
Z1 = op_neumannf(:,3);
z2 = reshape(Z1,[201,201]);
x = 0:1/200:1;
y = 0:1/200:1;
surf(x,y,z2,'EdgeColor','black')