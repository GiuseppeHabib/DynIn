% shimmy wheel
function xprime = sistema(t,x,par)
% par=[a cpsi kpsi m It b c d e bt ct dt et sigma0 cf ee V ca ka Ia knl];
    xprime=zeros(5,1);
    % ============================= State Space ========================
    sigma=par.sigma0*exp(-x(3)^2);
    Fy=par.d*sin(par.c*atan((1-par.e)*par.b*x(3)+par.e*atan(par.b*x(3)))) ;
    tpy_finale=par.dt*cos(par.ct*atan((1-par.et).*par.bt*x(3)+par.et*atan(par.bt*x(3))));
    Mz=-tpy_finale.*Fy;

    xprime(1) = x(2);
    xprime(2) = -(par.cpsi/par.It)*x(2)-(par.kpsi/par.It)*x(1,1)+((-Fy*par.ee)+Mz)/par.It-par.ca*(x(2)-x(5))/par.It-par.ka*(x(1)-x(4))/par.It-par.knl*(x(1)-x(4))^3/par.It;
    xprime(3) = ((par.ee-par.a)/sigma)*x(2)+(par.V/sigma)*sin(x(1))-(par.V/sigma)*x(3)*cos(x(1))+(x(3)^2*x(2));
    xprime(4) = x(5);
    xprime(5) = par.ca*(x(2)-x(5))/par.Ia+par.ka*(x(1)-x(4))/par.Ia+par.knl*(x(1)-x(4))^3/par.Ia;

    % ===================================================================
end

% par=[a cpsi kpsi m It b c d e bt ct dt et sigma0 cf ee V  ca ka Ia knl];
%      1 2    3    4 5  6 7 8 9 10 11 12 13 14     15 16 17 18 19 20 21