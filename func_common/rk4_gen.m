function [y, dx] = rk4_gen(state_eq, modArgs, x, x0, zf, dt)
%p_ta, u_ta, p_tv, u_tv, x0)

if modArgs == 0 % more general case
    x = x';
    k1 = dt * feval(state_eq, x, zf, dt);
    k2 = dt * feval(state_eq, x+k1/2, zf, dt);
    k3 = dt * feval(state_eq, x+k2/2, zf, dt);
    k4 = dt * feval(state_eq, x+k3, zf, dt);
    y = x + k1/6 + k2/3 + k3/3 + k4/6;
    y = y';
    dx =  feval(state_eq,x,zf,dt);

else
    k1 = dt * feval(state_eq, modArgs, x, x0, zf);
    k2 = dt * feval(state_eq, modArgs, x+k1/2, x0, zf);
    k3 = dt * feval(state_eq, modArgs, x+k2/2, x0, zf);
    k4 = dt * feval(state_eq, modArgs, x+k3, x0, zf);

    y = x + k1/6 + k2/3 + k3/3 + k4/6;
    y = y';
    %dx = feval(state_eq, p, u, x, x0)';
    dx =  feval(state_eq,modArgs,x,x0,zf);
end

return