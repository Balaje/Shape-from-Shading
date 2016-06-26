load op_neumannf.txt
Z1 = op_neumannf(:,3);
z2 = reshape(Z1,[201,201]);
surf(x,y,z2,'EdgeColor','black')