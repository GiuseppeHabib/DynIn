function [parout,weight,x1_q,x2_q] = dyn_invpend(par)
% Inverted pendulum with NLTVA and delayed PD control
% For details about the mechanical system look at: 
% Szaksz, B., Stepan, G., Habib, G. (2024). Dynamical integrity 
% estimation in time delayed systems: a rapid iterative algorithm. 
% Journal of Sound and Vibration, 571, 118045.

n = 2;  %DoF
x = sym('x',[n,1]);             % state vector (symbolic)
dx = sym('dx',[n,1]);           % first derivative of the state vector (symbolic)
ddx = sym('ddx',[n,1]);         % second derivative of the state vector (symbolic)
xdel = sym('xdel',[n,1]);       % delayed state vector (symbolic)
dxdel = sym('dxdel',[n,1]);     % delayed derivative of the state vector (symbolic)
xcell = {x;dx;ddx;xdel;dxdel};

%% Equations of motion - To customize the system you should modify this part
eq1 = ddx(1)-sin(x(1))+2*par.zeta2*par.mu*par.gamma*(dx(1)-dx(2))+par.mu*par.gamma^2*(x(1)-x(2))+par.p*xdel(1)+par.d*dxdel(1) == 0;
eq2 = par.mu*(ddx(2)+2*par.zeta2*par.gamma*(dx(2)-dx(1))+par.gamma^2*(x(2)-x(1)))== 0;
eq = [eq1,eq2];

% mass and stiffness matrices, here we consider the proportional gain p
% as well in the stiffness matrix, since otherwise the natural angular
% frequencies would become complex (otherwise the upright position of the
% inverted pendulum is unstable).
M = [1,0
    0,par.mu];
K = [-1+par.gamma^2*par.mu+par.p, -par.gamma^2*par.mu
    -par.gamma^2*par.mu, par.gamma^2*par.mu];

%% modes of the undamped system
[U,D] = eig(inv(M)*K);
weight = [D(1,1),D(2,2),1,1];   % weights which are used later for the distance defitition

[mapcell,x1_q,x2_q] = semidiscmatr(eq,par.dt,par.r,xcell,'U',U,'var1',1,'var2',2); % x1_q and x2_q are the first and second variables of the state vector
parout = par;
parout.mapcell = mapcell;