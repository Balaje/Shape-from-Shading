function V = fluxfunction(u,v,solver,eps)

    if(solver == 'g') % Godunov Flux
        V = max([f1(max([u,0])),f1(min([v,0]))]);
    end
    if(solver == 'm') % Murman Roe Flux
        if(u==v)
            V = f1(u);
        else
            a = abs( (f1(u)-f1(v))/(u-v) );
            V = 0.5*(f1(u)+f1(v)-a*(v-u));
        end
    end
    if(solver == 'l') % Lax Friedrich Flux
        V = 0.5*(f1(u)+f1(v)-1/eps*(v-u));
    end
    if(solver == 'e') % Einguist Osher Flux
        V = f1(max(u,0))+f1(min(v,0)) - f1(0);
    end
    
end