function [parout,weight,x1_q,x2_q] = dyn_machining1DoF(par)
% 1DoF nonlinear model of a machining (turning) operation
% For details about the mechanical system look at: 
% Szaksz, B., Stepan, G., Habib, G. (2024). Dynamical integrity 
% estimation in time delayed systems: a rapid iterative algorithm. 
% Journal of Sound and Vibration, 571, 118045.

n = 1;  %DoF
x = sym('x',[n,1]);             % state vector (symbolic)
dx = sym('dx',[n,1]);           % first derivative of the state vector (symbolic)
ddx = sym('ddx',[n,1]);         % second derivative of the state vector (symbolic)
xdel = sym('xdel',[n,1]);       % delayed state vector (symbolic)
dxdel = sym('dxdel',[n,1]);     % delayed derivative of the state vector (symbolic)
xcell = {x;dx;ddx;xdel;dxdel};

%% Equations of motion - To customize the system you should modify this part
eq = ddx(1)+2*par.zeta1*dx(1)+x(1)== par.p*((xdel(1)-x(1))+par.eta2*(xdel(1)-x(1))^2+par.eta3*(xdel(1)-x(1))^3);
        
% mass and stiffness matrices
M = [1];
K = [1];

%% modes of the undamped system
[U,D] = eig(inv(M)*K);
weight = [diag(D).',ones(size(diag(D).'))];   % weights which are used later for the distance defitition

[mapcell,x1_q,x2_q] = semidiscmatr(eq,par.dt,par.r,xcell,'U',U,'var1',1,'var2',2);
parout = par;
parout.mapcell = mapcell;