function [parout,weight,x1_q,x2_q] = dyn_Duffing(par)
% delayed Duffing oscillator
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
eq = ddx== -2*par.zeta*dx+xdel-par.a*x.^3;
 
% mass and stiffness matrices are now simply 1
M = [1];
K = [1];
%%
% modes of the undamped system
[U,D] = eig(inv(M)*K);
weight = [diag(D).',ones(size(diag(D).'))];   % weights which are used later for the distance defitition

% call the semidiscmatr function to get the mapping matrices.
[mapcell,x1_q,x2_q] = semidiscmatr(eq,par.dt,par.r,xcell,'U',U,'var1',1,'var2',2);
parout = par;
parout.mapcell = mapcell;
