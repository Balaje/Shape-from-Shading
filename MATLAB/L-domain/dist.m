x=0:.01:1;
y=x;
m=length(x);
for i=1:m
    for k=1:m
d(i,k)=min([x(i),sqrt(x(i).^2+y(k).^2),y(k),sqrt((1-x(i)).^2+(1-y(k)).^2)]);
    end
end
contour(x,y,d);