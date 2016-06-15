function V = I(i,j,gray)

%     V = 1/sqrt(1+(2*pi*sin(2*pi*y)*cos(2*pi*x))^2+(2*pi*sin(2*pi*x)*cos(2*pi*y))^2);
%     V = 1/sqrt(5);
%     V = 1/sqrt(1+(16*y*(1-y)*(1-2*x))^2+(16*x*(1-x)*(1-2*y))^2);
        
%     if(gray(i,j,1) == 0)
%         V = 0.01;
%     else
%         V = gray(i,j,1);
%     end
V = gray(i,j,1);

end