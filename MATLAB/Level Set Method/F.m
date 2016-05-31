function V = F(x,y)
    if(y<-1/3 || y > 1/3)
        V = 1*x^0;
    elseif(y>=-1/3 && y<1/3)
        V = -1;
    end
end