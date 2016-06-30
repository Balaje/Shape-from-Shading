load op_neumannf.txt
Z1 = op_neumannf(:,3);
z2 = reshape(Z1,[200,200]);

x = linspace(0,1,200);
y = linspace(0,1,200);

figure(2)
surf(x,y,z2,'EdgeColor','none')
colormap gray
light('Position',[1 1 1],'Style','local')
lighting gouraud
axis([0.1,0.9,0.1,0.9,-1.3,-1])