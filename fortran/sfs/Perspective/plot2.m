load op_neumannf.txt
Z1 = op_neumannf(:,3);
z2 = reshape(Z1,[201,201]);
x = linspace(0,1,201);
y = linspace(0,1,201);
surf(x,y,z2','EdgeColor','black')
axis([0,1,0,1,-1.2,-0.4])
