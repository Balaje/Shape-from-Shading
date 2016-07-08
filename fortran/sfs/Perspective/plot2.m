load op_neumannf.txt
Z1 = op_neumannf(:,3);
z2 = reshape(Z1,[200,200]);

x = linspace(0,2,200);
y = linspace(0,1,200);

figure(1)
surf(x,y,z2','EdgeColor','none')
set(gca,'color','black')
set(gcf,'color','white')

set(gca,'visible','on')
colormap([0.6,0.6,0.6])
material dull
light('Position',[0 0 1],'Style','local')
lighting gouraud
axis equal tight
axis([0,2,0,1,-1.3,-0.4])