function V = fluxfunction(u,v,solver,eps)

    if(solver == 'g') % Godunov Flux
        V = max([f(max([u,0])),f(min([v,0]))]);
    end
    if(solver == 'm') % Murman Roe Flux
        if(u==v)
            V = f(u);
        else
            a = abs( (f(u)-f(v))/(u-v) );
            V = 0.5*(f(u)+f(v)-a*(v-u));
        end
    end
    if(solver == 'l') % Lax Friedrich Flux
        V = 0.5*(f(u)+f(v)-1/eps*(v-u));
    end
    if(solver == 'e') % Einguist Osher Flux
        V = f(max(u,0))+f(min(v,0)) - f(0);
    end
    
end