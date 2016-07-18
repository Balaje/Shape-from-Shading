load op_neumannf.txt
N = 400;
Z1 = op_neumannf(:,3);
z2 = reshape(Z1,[N,N]);

x = linspace(0,2,N);
y = linspace(0,2,N);

figure(1)
surf(x,y,z2','EdgeColor','none')
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Reconstruction of the figure');

set(gca,'color','black')
set(gcf,'color','white')

set(gca,'visible','on')
colormap([0.6,0.6,0.6])
material dull
light('Position',[0 0 1],'Style','local')
lighting gouraud
%axis equal tight
%axis([0,1,0,1,-0.7,-0.4])