function [parout] = dyn_Simplified_Duffing(par)
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
 
% call the semidiscmatr function to get the mapping matrices.
mapcell = semidiscmatr(eq,par.dt,par.r,xcell);
parout = par;
parout.mapcell = mapcell;
