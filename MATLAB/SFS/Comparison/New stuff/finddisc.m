function V = finddisc(I)
    ddI = diff(diff(I));
    eps = 0.1;
    V = 0;
    for i=1:length(ddI)
        if(abs(ddI(i)) > eps)
            V = [V i+1];
        end
    end 
end