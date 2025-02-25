function [par,weight,x0]=coefficients_chain_of_4_masses(lr,c)

% syms x1dd x2dd x3dd x4dd x1d x2d x3d x4d x1 x2 x3 x4
% 
% eq01a=2*x1+2*c*x1d+x1dd-x2-c*x2d;
% eq01b=(-(1/6))*sqrt(-289+144*lr^2)-x1-c*x1d+4*x2+(2*lr*sqrt(-289+144*lr^2))/(17*sqrt(1+((-(1/17))*sqrt(-289+144*lr^2)+x2)^2))-(2*lr*x2)/sqrt(1+((-(1/17))*sqrt(-289+144*lr^2)+x2)^2)+2*c*x2d+x2dd-x3-c*x3d;
% eq01c=-x2-c*x2d+2*x3+2*c*x3d+x3dd-x4-c*x4d;
% eq01d=-x3-c*x3d+2*x4+2*c*x4d+x4dd;
if lr>=17/12
    M=eye(4);
    C=c*[2,-1,0,0; -1,2,-1,0; 0,-1,2,-1; 0,0,-1,2];
    % K=[2,-1,0,0; -1,double(subs(diff(eq01b,x2),x2,0)),-1,0; 0,-1,2,-1; 0,0,-1,2];
    K=[2,-1,0,0; -1,(-4913+3456*lr^2)/(864*lr^2),-1,0; 0,-1,2,-1; 0,0,-1,2];
    
    [par.vect,par.val]=eig(M\K);
    
    Ms=transpose(par.vect)*M*par.vect;
    par.Cs=transpose(par.vect)*C*par.vect;
    par.Ks=transpose(par.vect)*K*par.vect;
    par.Ut=transpose(par.vect);
    
    weight=[par.val(1,1), 1, par.val(2,2), 1, par.val(3,3), 1, par.val(4,4), 1];
    
    par.lr=lr;
    par.c=c;
    
    x0=[(-(1/34))*sqrt(-289 + 144*lr^2),0 , (-(1/17))*sqrt(-289 + 144*lr^2),0 , (-(2/51))*sqrt(-289 + 144*lr^2),0 , (-(1/51))*sqrt(-289 + 144*lr^2),0 ];
    % par=[Cs(1,1), Cs(1,2), Cs(1,3), Cs(1,4), Cs(2,1), Cs(2,2), Cs(2,3), Cs(2,4), Cs(3,1), Cs(3,2), Cs(3,3), Cs(3,4), Cs(4,1), Cs(4,2), Cs(4,3), Cs(4,4), ...
    %     vact(2,1), vect(2,2), vect(2,3), vect(2,4),
else
    warning('the systm has only one stable point')
end


