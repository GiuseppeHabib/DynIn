function [parout,weight,x1_q,x2_q] = dyn_machining2DoF(par)

n = 2;  %DoF if modal or dimension if it is odd
x = sym('x',[n,1]);             % state vector (symbolic)
dx = sym('dx',[n,1]);           % first derivative of the state vector (symbolic)
ddx = sym('ddx',[n,1]);         % second derivative of the state vector (symbolic)
xdel = sym('xdel',[n,1]);       % delayed state vector (symbolic)
dxdel = sym('dxdel',[n,1]);     % delayed derivative of the state vector (symbolic)
xcell = {x;dx;ddx;xdel;dxdel};

%% Equations of motion - To customize the system you should modify this part
M0 = [1,0;
      0,par.mu];
C0 = [2*par.zeta1+2*par.zeta2*par.gamma*par.mu,-2*par.zeta2*par.gamma*par.mu;
      -2*par.zeta2*par.gamma*par.mu,2*par.zeta2*par.gamma*par.mu];
K0 = [1+par.gamma^2*par.mu,-par.gamma^2*par.mu;
    -par.gamma^2*par.mu,par.gamma^2*par.mu];
RHS = [-par.alpha3*(x(1)-x(2))^3+par.p*((xdel(1)-x(1))+par.eta2*(xdel(1)-x(1))^2+par.eta3*(xdel(1)-x(1))^3);-par.alpha3*(x(2)-x(1))^3];
% Equations of motion
eq = M0*ddx+C0*dx+K0*x == RHS;

% mass and stiffness matrices
M = M0;
K = K0;
%% modes of the undamped system
[U,D] = eig(inv(M)*K);
weight = [diag(D).',ones(1,n)];   % weights which are used later for the distance defitition
[mapcell,x1_q,x2_q] = semidiscmatr(eq,par.dt,par.r,xcell,'U',U,'var1',1,'var2',2);
parout = par;
parout.mapcell = mapcell;
end