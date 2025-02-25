function x=sys_solver(Ts,x0,par,tols)

option=odeset('RelTol', tols(1), 'AbsTol', tols(2),'Events',@largeangle);

[~,x]=ode45 (@(t,xt) sistema(t,xt,par), [0 Ts], x0, option);
